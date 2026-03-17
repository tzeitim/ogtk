# Masking Optimization Plan

## Overview

Currently, the masking system uses deterministic but random character-level shuffling to generate scrambled sequences. This can result in:
- Short kmer overlaps with original sequences (e.g., 11bp substring "TGGCAGTCATC" appearing in both)
- Suboptimal assembly performance
- Luck-dependent effectiveness (good scrambling vs bad scrambling)

**Goal**: Create a masking optimization tool that pre-generates optimal scrambled sequences for each META-TARGET combination, evaluated against specific assembly criteria.

---

## Problem Statement

### Current Approach
```python
def scramble_dna_sequence(seed_string: str, sequence: str):
    seed = int(hashlib.md5(seed_string.encode()).hexdigest()[:8], 16)
    random.seed(seed)
    seq_list = list(sequence.upper())
    random.shuffle(seq_list)  # Character-level shuffle
    return ''.join(seq_list)
```

**Issues:**
1. Short kmers can randomly reappear (statistical inevitability)
2. No guarantee of minimal kmer overlap for a given k
3. Cannot optimize for specific assembly parameters
4. GC content distribution may change dramatically

### Desired Approach

Generate **optimized** scrambled sequences that:
- Minimize kmer overlap with original for given k (e.g., k=25, k=50)
- Preserve similar GC content distribution
- Avoid self-complementarity issues
- Are deterministic and reproducible
- Are pre-computed and stored in features CSV

---

## Implementation Plan

### Phase 1: Kmer Analysis Module

**Location**: `ogtk/ltr/fracture/pipeline/kmer_utils.py`

```python
def extract_kmers(sequence: str, k: int) -> set[str]:
    """Extract all kmers of length k from sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def kmer_overlap_score(seq1: str, seq2: str, k: int) -> float:
    """
    Calculate kmer overlap between two sequences.

    Returns:
        Fraction of kmers from seq1 that appear in seq2 (0.0 = no overlap, 1.0 = complete overlap)
    """
    kmers1 = extract_kmers(seq1, k)
    kmers2 = extract_kmers(seq2, k)
    if not kmers1:
        return 0.0
    return len(kmers1 & kmers2) / len(kmers1)

def gc_content(sequence: str) -> float:
    """Calculate GC content of sequence."""
    gc = sequence.count('G') + sequence.count('C')
    return gc / len(sequence) if len(sequence) > 0 else 0.0

def gc_distribution_distance(seq1: str, seq2: str, window: int = 10) -> float:
    """
    Calculate difference in GC content distribution between sequences.
    Uses sliding window to compare local GC content.

    Returns:
        Mean absolute difference in GC content across windows
    """
    gc1 = [gc_content(seq1[i:i+window]) for i in range(0, len(seq1) - window + 1, window)]
    gc2 = [gc_content(seq2[i:i+window]) for i in range(0, len(seq2) - window + 1, window)]

    # Pad shorter list
    max_len = max(len(gc1), len(gc2))
    gc1.extend([0] * (max_len - len(gc1)))
    gc2.extend([0] * (max_len - len(gc2)))

    return sum(abs(g1 - g2) for g1, g2 in zip(gc1, gc2)) / len(gc1)
```

---

### Phase 2: Sequence Shuffling Strategies

**Location**: `ogtk/ltr/fracture/pipeline/shuffle_strategies.py`

```python
import random
from typing import Callable

def character_shuffle(sequence: str, seed: int) -> str:
    """Current approach - shuffle individual characters."""
    random.seed(seed)
    seq_list = list(sequence.upper())
    random.shuffle(seq_list)
    return ''.join(seq_list)

def codon_shuffle(sequence: str, seed: int) -> str:
    """
    Shuffle in triplets (codons) to maintain local structure.
    Better for preserving local GC distribution.
    """
    random.seed(seed)
    seq = sequence.upper()
    # Pad to multiple of 3
    padding = (3 - len(seq) % 3) % 3
    seq += 'N' * padding

    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    random.shuffle(codons)
    result = ''.join(codons)

    # Remove padding
    return result[:len(sequence)]

def block_shuffle(sequence: str, seed: int, block_size: int = 10) -> str:
    """
    Shuffle in blocks to preserve local sequence features.
    Larger blocks = more local structure preserved.
    """
    random.seed(seed)
    seq = sequence.upper()

    blocks = [seq[i:i+block_size] for i in range(0, len(seq), block_size)]
    random.shuffle(blocks)
    result = ''.join(blocks)

    return result[:len(sequence)]

def reverse_complement(sequence: str) -> str:
    """
    Reverse complement as an alternative to shuffling.
    Maintains composition but inverts sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(sequence.upper()))

# Strategy registry
STRATEGIES: dict[str, Callable[[str, int], str]] = {
    'character': character_shuffle,
    'codon': codon_shuffle,
    'block_10': lambda seq, seed: block_shuffle(seq, seed, block_size=10),
    'block_20': lambda seq, seed: block_shuffle(seq, seed, block_size=20),
    'reverse_complement': lambda seq, seed: reverse_complement(seq),
}
```

---

### Phase 3: Optimization Engine

**Location**: `ogtk/ltr/fracture/pipeline/optimize_masking.py`

```python
import polars as pl
from typing import Optional
from .kmer_utils import kmer_overlap_score, gc_distribution_distance
from .shuffle_strategies import STRATEGIES

def score_candidate(
    original: str,
    candidate: str,
    k_values: list[int],
    gc_weight: float = 0.3,
    kmer_weight: float = 0.7
) -> float:
    """
    Score a candidate scrambled sequence.

    Lower score = better (less overlap, similar GC distribution)

    Args:
        original: Original TARGET sequence
        candidate: Candidate scrambled sequence
        k_values: List of kmer sizes to test against (e.g., [25, 50])
        gc_weight: Weight for GC distribution similarity (0-1)
        kmer_weight: Weight for kmer overlap minimization (0-1)

    Returns:
        Combined score (lower is better)
    """
    # Kmer overlap score (average across all k values)
    kmer_scores = [kmer_overlap_score(original, candidate, k) for k in k_values]
    avg_kmer_overlap = sum(kmer_scores) / len(kmer_scores)

    # GC distribution score
    gc_score = gc_distribution_distance(original, candidate)

    # Combined score (weighted)
    return kmer_weight * avg_kmer_overlap + gc_weight * gc_score

def generate_candidates(
    sequence: str,
    n_candidates: int = 1000,
    strategies: Optional[list[str]] = None
) -> list[tuple[str, str, str]]:
    """
    Generate candidate scrambled sequences using various strategies and seeds.

    Args:
        sequence: Original sequence to scramble
        n_candidates: Number of candidates to generate
        strategies: List of strategy names to use (default: all)

    Returns:
        List of (strategy_name, seed_string, scrambled_sequence) tuples
    """
    if strategies is None:
        strategies = list(STRATEGIES.keys())

    candidates = []
    per_strategy = n_candidates // len(strategies)

    for strategy_name in strategies:
        strategy_func = STRATEGIES[strategy_name]
        for i in range(per_strategy):
            seed = i
            scrambled = strategy_func(sequence, seed)
            candidates.append((strategy_name, f"{strategy_name}_{i}", scrambled))

    return candidates

def optimize_masking(
    features_csv: str,
    k_values: list[int] = [25, 50],
    n_candidates: int = 1000,
    output_csv: Optional[str] = None
) -> pl.DataFrame:
    """
    Optimize masking sequences for minimal kmer overlap.

    Args:
        features_csv: Input features CSV with META and TARGET sequences
        k_values: Kmer sizes to optimize for (e.g., [25, 50])
        n_candidates: Number of candidate sequences to test per TARGET
        output_csv: Optional output path for optimized features CSV

    Returns:
        DataFrame with optimized scrambled sequences added
    """
    meta = pl.read_csv(features_csv)

    # Get all META-TARGET combinations
    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'TARGET'),
            suffix="_target",
            how='cross'
        )
    )

    optimized_results = []

    for row in patterns_df.iter_rows(named=True):
        meta_name = row['feature']
        target_name = row['feature_target']
        target_seq = row['seq_target']

        print(f"Optimizing {meta_name} + {target_name}...")

        # Generate candidates
        candidates = generate_candidates(target_seq, n_candidates)

        # Score each candidate
        scored = [
            (strategy, seed, seq, score_candidate(target_seq, seq, k_values))
            for strategy, seed, seq in candidates
        ]

        # Select best (lowest score)
        best = min(scored, key=lambda x: x[3])
        strategy, seed, scrambled_seq, score = best

        print(f"  Best: {strategy} (score={score:.4f})")

        optimized_results.append({
            'meta_feature': meta_name,
            'target_feature': target_name,
            'target_seq': target_seq,
            'optimized_seq': scrambled_seq,
            'strategy': strategy,
            'seed': seed,
            'score': score
        })

    result_df = pl.DataFrame(optimized_results)

    if output_csv:
        result_df.write_csv(output_csv)
        print(f"\nSaved optimized sequences to {output_csv}")

    return result_df
```

---

### Phase 4: CLI Tool

**Location**: `ogtk/ltr/fracture/cli_optimize_masking.py`

```python
#!/usr/bin/env python
"""
CLI tool to optimize masking sequences.

Usage:
    python -m ogtk.ltr.fracture.cli_optimize_masking \\
        --input features.csv \\
        --output optimized_features.csv \\
        --k 25 50 \\
        --candidates 10000
"""
import argparse
from .pipeline.optimize_masking import optimize_masking

def main():
    parser = argparse.ArgumentParser(
        description='Optimize masking sequences for minimal kmer overlap'
    )
    parser.add_argument('--input', required=True,
                       help='Input features CSV')
    parser.add_argument('--output', required=True,
                       help='Output optimized features CSV')
    parser.add_argument('--k', nargs='+', type=int, default=[25, 50],
                       help='Kmer sizes to optimize for (default: 25 50)')
    parser.add_argument('--candidates', type=int, default=1000,
                       help='Number of candidates to test (default: 1000)')

    args = parser.parse_args()

    print(f"Optimizing masking for k={args.k}")
    print(f"Testing {args.candidates} candidates per TARGET sequence...")

    result_df = optimize_masking(
        features_csv=args.input,
        k_values=args.k,
        n_candidates=args.candidates,
        output_csv=args.output
    )

    print("\n=== Optimization Summary ===")
    print(result_df.select(['meta_feature', 'target_feature', 'strategy', 'score']))
    print(f"\nBest average score: {result_df['score'].mean():.4f}")

if __name__ == '__main__':
    main()
```

---

### Phase 5: Update Masking Module to Use Optimized Sequences

**Modify**: `ogtk/ltr/fracture/pipeline/masking.py`

```python
def generate_mask_PEtracer_expression(
    features_csv: str,
    column_name: str = "seq",
    fuzzy_pattern: bool = True,
    fuzzy_kwargs: Optional[dict] = None,
    use_optimized: bool = True,  # NEW parameter
    logger: Optional[CustomLogger] = None
) -> pl.Expr:
    """
    Generate masking expression.

    Args:
        use_optimized: If True, use 'optimized_seq' column from CSV.
                      If False, generate scrambled sequences on-the-fly.
    """
    if fuzzy_kwargs is None:
        fuzzy_kwargs = {}

    meta = pl.read_csv(features_csv)

    # Check if optimized sequences are available
    has_optimized = 'optimized_seq' in meta.columns and use_optimized

    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'TARGET'),
            suffix="_target",
            how='cross'
        )
        .with_columns(
            # Use optimized sequence if available, otherwise scramble
            scrambled_seq = (
                pl.col('optimized_seq') if has_optimized
                else pl.struct(['feature', 'feature_target']).map_elements(
                    lambda x: scramble_dna_sequence(
                        seed_string=f"{x['feature']}_{x['feature_target']}",
                        sequence=meta.filter(pl.col('feature') == x['feature_target']).get_column('seq')[0],
                        fuzzy_pattern=False,
                        fuzzy_kwargs=fuzzy_kwargs
                    ),
                    return_dtype=pl.Utf8
                )
            )
        )
        # ... rest of pattern generation
    )
```

---

## Usage Workflow

### Step 1: Optimize Masking Sequences

```bash
python -m ogtk.ltr.fracture.cli_optimize_masking \
    --input ogtk/ltr/fracture/extensions/PEtracer_metas.csv \
    --output ogtk/ltr/fracture/extensions/PEtracer_metas_optimized.csv \
    --k 25 50 \
    --candidates 10000
```

**Output CSV format**:
```csv
meta_feature,target_feature,target_seq,optimized_seq,strategy,seed,score
META01,RNF2,TGGCAGTCATC...,GATCGTAC...,block_10,block_10_142,0.0234
META01,HEK3,CTTGGGGCCC...,CGGGCTT...,codon,codon_87,0.0189
...
```

### Step 2: Use Optimized Sequences in Pipeline

```python
# Option 1: Load optimized features directly
df = (
    pl.scan_parquet('reads.parquet')
    .pp.mask_repeats('PEtracer_metas_optimized.csv', use_optimized=True)
    .pp.assemble_umis(k=50, ...)
)

# Option 2: Fallback to on-the-fly scrambling if optimized not available
df = (
    pl.scan_parquet('reads.parquet')
    .pp.mask_repeats('PEtracer_metas.csv', use_optimized=False)
    .pp.assemble_umis(k=50, ...)
)
```

---

## Testing Plan

### Test 1: Kmer Overlap Validation
```python
def test_kmer_overlap():
    original = "TGGCAGTCATCTTAGTCATTACGACAGGTGTTCGTTGTAACTCATATA"

    # Current scrambling
    current = scramble_dna_sequence("test", original, False, {})
    current_overlap = kmer_overlap_score(original, current, k=25)

    # Optimized scrambling
    candidates = generate_candidates(original, n_candidates=1000)
    optimized = min(candidates, key=lambda x: score_candidate(original, x[2], [25, 50]))[2]
    optimized_overlap = kmer_overlap_score(original, optimized, k=25)

    assert optimized_overlap < current_overlap, "Optimized should have less overlap"
```

### Test 2: Assembly Performance
```python
def test_assembly_performance():
    # Compare assembly success rates with:
    # 1. No masking
    # 2. Current random masking
    # 3. Optimized masking

    # Measure: success rate, contig length, assembly time
```

---

## Expected Improvements

1. **Kmer overlap**: Reduce from ~15-20% to <5% for k=25
2. **Assembly success**: Increase from 17.91% to >25%
3. **Reproducibility**: 100% deterministic across runs
4. **Performance**: Pre-computed sequences = faster masking

---

## Future Enhancements

1. **Multi-objective optimization**: Balance kmer overlap, GC content, self-complementarity
2. **Evolutionary algorithms**: Use genetic algorithms for better optimization
3. **Target-aware optimization**: Consider specific target sequences in cassette design
4. **Batch optimization**: Optimize multiple features CSVs at once
5. **Validation metrics**: Automated testing of assembly performance improvements
