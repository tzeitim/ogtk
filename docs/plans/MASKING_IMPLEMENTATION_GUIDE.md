# Masking Implementation Guide

## Overview
This guide documents the implementation of sequence masking/unmasking functionality to handle very long cassettes with multiple direct repeats that exceed the 64-base kmer limit in rogtk.

## Background
The masking system consists of two complementary operations:

1. **Masking (before assembly)**: `generate_mask_PEtracer_expression()` replaces repetitive TARGET sequences with deterministic scrambled versions based on their preceding META context. This allows the kmer-based assembler to work on unique/variable regions while preserving structural information.

2. **Unmasking (after assembly)**: `generate_unmask_PEtracer_expression()` restores the original TARGET sequences in assembled contigs by replacing the scrambled versions back to their original form.

---

## Implementation Checklist

### ✅ Step 1: Extract the masking function (COMPLETED)
- [x] Create new module: `ogtk/ltr/fracture/pipeline/masking.py`
- [x] Move `generate_mask_PEtracer_expression()` from `cassiopeia_lineage.py`
- [x] Move helper function `scramble_dna_sequence()`
- [x] Import `fuzzy_match_str` from existing `ogtk.utils.general`
- [x] Update `cassiopeia_lineage.py` to import from new location
- [x] Remove old function definition from `cassiopeia_lineage.py`
- [x] Clean up `cassiopeia_petracer.py` (removed stray `def mask` line)
- [x] Verify imports work correctly

**Files Modified:**
- Created: `ogtk/ltr/fracture/pipeline/masking.py`
- Modified: `ogtk/ltr/fracture/extensions/cassiopeia_lineage.py`
- Modified: `ogtk/ltr/fracture/cassiopeia_petracer.py`

---

### ✅ Step 2: Add `.pp.mask_repeats()` method (COMPLETED)
**Location:** `ogtk/ltr/fracture/pipeline/api_ext.py`

**Tasks:**
- [x] Import `generate_mask_PEtracer_expression` from masking module
- [x] Add `mask_repeats()` method to `PllPipeline` class (LazyFrame namespace)
- [x] Pass through configuration parameters (fuzzy_kwargs, column_name, etc.)
- [x] Ensure logger is passed to the masking function
- [x] Add comprehensive docstring with examples
- [x] Test method registration and chaining

**Implementation:**
```python
# In imports section
from .masking import generate_mask_PEtracer_expression

# In PllPipeline class (after other methods like parse_reads)
def mask_repeats(self,
                 features_csv: str,
                 column_name: str = 'r2_seq',
                 fuzzy_pattern: bool = True,
                 fuzzy_kwargs: dict = None) -> pl.LazyFrame:
    """
    Mask repetitive sequences based on META-TARGET pairs.

    This is useful for very long cassettes with direct repeats that exceed
    the kmer length limit (typically 64 bases in rogtk).

    Args:
        features_csv: Path to CSV file with META and TARGET sequences
        column_name: Column containing sequences to mask (default: 'r2_seq')
        fuzzy_pattern: Whether to use fuzzy matching for sequencing errors
        fuzzy_kwargs: Optional dict with fuzzy matching parameters
                     (wildcard, include_original, sep, max_length)

    Returns:
        LazyFrame with masked sequences

    Example:
        df = (
            pl.scan_parquet('reads.parquet')
            .pp.mask_repeats('features.csv', column_name='r2_seq')
            .pp.assemble_umis(...)
        )
    """
    mask_expr = generate_mask_PEtracer_expression(
        features_csv=features_csv,
        column_name=column_name,
        fuzzy_pattern=fuzzy_pattern,
        fuzzy_kwargs=fuzzy_kwargs,
        logger=self.logger
    )
    return self._ldf.with_columns(mask_expr.alias(column_name))
```

**Testing:** ✅ PASSED
```python
# Test 1: Method registration
import polars as pl
from ogtk.ltr.fracture.pipeline.api_ext import PllPipeline
df = pl.DataFrame({'r2_seq': ['ATCG']}).lazy()
assert hasattr(df.pp, 'mask_repeats')  # ✓ PASSED

# Test 2: Basic functionality
df = pl.DataFrame({'r2_seq': ['AAAATTTTGGGGDDDD']}).lazy()
result = df.pp.mask_repeats('/tmp/test_features.csv', column_name='r2_seq')
collected = result.collect()  # ✓ PASSED

# Test 3: Method chaining
df = pl.DataFrame({'r2_seq': ['AAAATTTTGGGGDDDD'], 'umi': ['UMI001'], 'reads': [100]}).lazy()
result = df.pp.mask_repeats('/tmp/test_features.csv').collect()  # ✓ PASSED
```

**Files Modified:**
- Modified: `ogtk/ltr/fracture/pipeline/api_ext.py` (added import and mask_repeats method)

---

### ✅ Step 3: Integrate into FRACTURE step (COMPLETED)
**Location:** `ogtk/ltr/fracture/pipeline/core.py`

**Tasks:**
- [x] Check for `xp.features_csv` at start of `fracture()` method
- [x] Apply masking to `ldf` if configured
- [x] Log basic statistics (number of reads)
- [x] Add intermediate file saving with `save_intermediate_files` flag
- [x] Add logging for masking operations

**Implementation:**
The masking logic has been integrated into `core.py:779-805`. Key changes:
- Load parquet into a LazyFrame variable (`ldf`)
- Check for `features_csv` configuration
- Apply masking using `.pp.mask_repeats()` method
- Log statistics (number of reads masked)
- Optionally save intermediate masked files if `save_intermediate_files: true`
- Continue with existing assembly logic using the masked `ldf`

The implementation respects all configuration parameters:
- `features_csv`: Path to features CSV (required to enable masking)
- `mask_fuzzy_pattern`: Enable fuzzy matching (default: True)
- `mask_fuzzy_kwargs`: Custom fuzzy matching parameters (optional)
- `save_intermediate_files`: Save masked reads to intermediate/ directory (default: False)

---

### ✅ Step 3b: Implement Unmasking (COMPLETED)

**Why unmasking is needed:**
After assembly, contigs contain scrambled TARGET sequences. These must be restored to their original form for downstream analysis, otherwise the assembled sequences won't match expected targets.

**Implementation:**

1. **Added `generate_unmask_PEtracer_expression()` in `masking.py`** (lines 181-282)
   - Reverses the masking operation
   - Uses the same deterministic scrambling to identify scrambled sequences
   - Replaces scrambled TARGET sequences with original TARGET sequences
   - Pattern: `(META)(.*?)(SCRAMBLED_TARGET)` → `$1${2}<ORIGINAL_TARGET>$4`

2. **Added `unmask_repeats()` method in `api_ext.py`** (lines 673-717)
   - User-facing method in the `.pp` namespace
   - Works on LazyFrame (typically assembled contigs)
   - Parameters must match those used during masking

3. **Integrated into FRACTURE step in `core.py`** (lines 838-850)
   - Applied automatically after assembly if `features_csv` is configured
   - Logs "Restoring original sequences (unmasking contigs)"
   - Recalculates contig length after unmasking

**Critical Note:**
The `fuzzy_pattern` and `fuzzy_kwargs` parameters used during unmasking **must exactly match** those used during masking, otherwise the scrambled sequences won't be recognized.

**Example Usage:**
```python
# Programmatic usage
df_contigs = (
    assembled_contigs
    .pp.unmask_repeats('features.csv', column_name='contig')
)

# Automatic usage in pipeline
# Just set features_csv in config - unmasking happens automatically
features_csv: "${prefix}/conf/features.csv"
```

---

### Step 4: Config integration

**Tasks:**
- [ ] Document `features_csv` parameter
- [ ] Document `save_intermediate_files` flag
- [ ] Document optional masking parameters
- [ ] Add example features.csv format to docs

**Config Parameters:**

```yaml
# Basic masking (required if you want masking enabled)
features_csv: "${prefix}/path/to/features.csv"

# Optional: save intermediate files for debugging
save_intermediate_files: true  # default: false

# Optional: reverse complement sequences before masking (for opposite strand orientation)
mask_reverse_complement: true  # default: false

# Optional: customize fuzzy matching
mask_fuzzy_pattern: true  # default: true
mask_fuzzy_kwargs:
  wildcard: '.{0,2}'        # Allow up to 2 character variations
  include_original: true    # Include exact match in pattern
  sep: '|'                  # Separator for alternatives
  max_length: 150           # Only fuzzify sequences up to this length
```

**Features CSV Format:**
```csv
feature,seq,kind
META01,AGAAGCCGTGTGCCGGTCTA,META
META02,ATCGTGCGGACGAGACAGCA,META
RNF2,TGGCAGTCATCTTAGTCATTACGACAGGTGTTCGTTGTAACTCATATA,TARGET
HEK3,CTTGGGGCCCAGACTGAGCACGACTTGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCA,TARGET
EMX1,GGCCTGAGTCCGAGCAGAAGAACTTGGGCTCCCATCACATCAACCGGTGG,TARGET
```

**Example Features CSV:**
A working example file is available at:
`ogtk/ltr/fracture/extensions/PEtracer_metas.csv`

This file contains 19 META sequences and 3 TARGET sequences (RNF2, HEK3, EMX1) for PEtracer cassettes.

**Package Data Options (for future reference):**

To provide internal example files with the package, consider:

1. **Helper function approach:**
```python
from pathlib import Path
from importlib.resources import files

def get_example_features_csv() -> str:
    return str(files('ogtk.ltr.fracture.extensions').joinpath('PEtracer_metas.csv'))
```

2. **Simple path reference in masking.py:**
```python
EXAMPLE_FEATURES_CSV = Path(__file__).parent.parent / 'extensions' / 'PEtracer_metas.csv'
```

3. **Package configuration (setup.py/pyproject.toml):**
```toml
[tool.setuptools.package-data]
ogtk = ["ltr/fracture/extensions/*.csv"]
```

---

### Step 5: Testing

**Tasks:**
- [ ] Test with `features_csv` enabled
- [ ] Test with `save_intermediate_files: true`
- [ ] Verify masking happens before assembly
- [ ] Check that cassiopeia extension still works after refactoring
- [ ] Test that masked sequences assemble correctly
- [ ] Verify intermediate files are saved in correct location
- [ ] Check logging output is informative

**Test Cases:**

1. **Basic masking test:**
```bash
# Config with features_csv enabled
python ~/src/ogtk/ogtk/ltr/fracture/cli.py \
  --config config_with_masking.yml \
  --target-sample barcode13
```

2. **Debug mode test:**
```yaml
features_csv: "${prefix}/conf/features.csv"
save_intermediate_files: true
```
Check that `workdir/intermediate/masked_reads_valid.parquet` is created.

3. **No masking test:**
```yaml
# features_csv: not set
```
Should run normally without masking.

4. **Verify no regression:**
Test that existing pipelines without masking still work correctly.

---

## Design Rationale

### Why Option 4 (LazyFrame method)?
- **Maximum flexibility**: Can be called from any step or used programmatically
- **Follows existing patterns**: Uses the `.pp.*` namespace already established
- **Chainable**: Works naturally with other Polars operations
- **Optional**: Doesn't force masking on all users

### Why mask in FRACTURE step?
- **Just-in-time**: Masking happens right before it's needed
- **Lazy evaluation**: Efficient with Polars lazy API
- **Configuration-driven**: Automatically enabled based on `features_csv` presence
- **Debuggable**: Optional intermediate file saving with `save_intermediate_files`

### Trade-offs
**Pros:**
- Clean separation of concerns
- No extra I/O in normal operation
- Easy to debug with flag
- Works with existing pipeline architecture

**Cons:**
- Masking logic coupled with FRACTURE step
- Can't reuse masked reads without `save_intermediate_files: true`
- Need to re-mask if running FRACTURE multiple times

---

## Files Modified

1. ✅ `ogtk/ltr/fracture/pipeline/masking.py`
   - Created new module with masking utilities
   - `generate_mask_PEtracer_expression()` - mask repetitive sequences before assembly
   - `generate_unmask_PEtracer_expression()` - restore original sequences after assembly
   - `scramble_dna_sequence()` - deterministic sequence scrambling helper

2. ✅ `ogtk/ltr/fracture/extensions/cassiopeia_lineage.py`
   - Updated to import from new masking module

3. ✅ `ogtk/ltr/fracture/cassiopeia_petracer.py`
   - Cleaned up stray code

4. ✅ `ogtk/ltr/fracture/pipeline/api_ext.py`
   - Added `mask_repeats()` method to PllPipeline class
   - Added `unmask_repeats()` method to PllPipeline class
   - Updated imports

5. ✅ `ogtk/ltr/fracture/pipeline/core.py`
   - Integrated masking into FRACTURE step (before assembly)
   - Integrated unmasking into FRACTURE step (after assembly)
   - Added logging and intermediate file saving

6. ✅ `ogtk/ltr/fracture/extensions/PEtracer_metas.csv`
   - Example features file with 19 META and 3 TARGET sequences

7. [ ] Config files (TO DO - document parameters)

---

## Notes

- The masking function was originally in `cassiopeia_lineage.py` as extension-specific code
- It has been extracted to `masking.py` to make it available to the core pipeline
- The cassiopeia extension can still use it via import
- `cassiopeia_petracer.py` appears to be legacy code and may need deprecation later

---

## Questions/TODOs

- Should we add masking statistics to the pipeline metrics?
- Do we need a separate MASK pipeline step, or is integration in FRACTURE sufficient?
- Should `save_intermediate_files` save ALL intermediate steps or just masking?
- What happens if features_csv file doesn't exist? (add validation)
- ✅ ~~Need to implement unmasking after assembly~~ - COMPLETED

## Summary of Complete Implementation

The masking/unmasking system is now fully integrated:

1. **Pre-assembly**: Reads are masked using `.pp.mask_repeats()`
   - Repetitive TARGET sequences → scrambled versions
   - Enables kmer-based assembly despite direct repeats

2. **Assembly**: rogtk assembles using unique/scrambled sequences
   - Works correctly even with long cassettes (>64bp kmers)

3. **Post-assembly**: Contigs are unmasked using `.pp.unmask_repeats()`
   - Scrambled sequences → original TARGET sequences
   - Final contigs contain correct biological sequences

4. **Automatic**: Both operations happen automatically when `features_csv` is configured
   - No manual intervention required
   - Parameters are preserved between masking and unmasking

## Bug Fixes

### ✅ Fixed Fuzzy Pattern Application (Nov 6, 2025)

**Issue**: Fuzzy patterns were being applied to replacement strings, creating regex patterns instead of plain DNA sequences in masked reads.

**Root cause**: `scramble_dna_sequence()` was applying `fuzzy_match_str()` to create scrambled sequences, resulting in replacements like `A.{0,2}T.{0,2}G...` instead of `ATCG...`

**Fix**:
- Scrambled sequences are now always generated as plain DNA (`fuzzy_pattern=False`)
- Fuzzy matching is applied only to **search patterns** (META and TARGET) for better matching with sequencing errors
- Replacements are always clean DNA sequences

**Files modified**:
- `masking.py` lines 135, 143-156 (mask function)
- `masking.py` lines 250, 259-272 (unmask function)

### ✅ Added Strand Orientation Support (Nov 10, 2025)

**Issue**: Depending on how data was generated, reads relative to the UMI can occur in sense or antisense orientation. The masking pattern `(META)(.*?)(TARGET)` only matches when sequences are in the expected orientation. This resulted in ~78.5% of sequences not being masked when cassettes were in reverse-complement orientation.

**Solution**: Added `mask_reverse_complement` configuration parameter to apply reverse complement before masking and restore original orientation after assembly.

**Implementation**:
- Reverse complement applied using rogtk's Rust-based `.dna.reverse_complement()` Polars extension
- Applied before masking in core.py:790-795
- Sequences remain in reverse-complemented orientation after assembly (not flipped back)
- Default: `false` (maintains backward compatibility with original orientation assumption)

**Usage**:
```yaml
# For data where cassette is in opposite orientation relative to UMI
mask_reverse_complement: true
```

**Design rationale**:
- Step-specific flag avoids breaking assumption that UMI is at read start
- Cleaner than bidirectional masking (which would mask PCR artifacts)
- Uses Rust implementation for performance (not Python-based rev_comp)
- Sequences are NOT flipped back after unmasking - they remain in the orientation that makes anchor sequences work correctly

**Files modified**:
- `core.py` lines 790-795
- `MASKING_IMPLEMENTATION_GUIDE.md` (documentation)

### ✅ Fixed Non-Deterministic Scrambling (Nov 13, 2025)

**Issue**: Pipeline results varied significantly between runs (e.g., 17.91% vs 2.37% vs 2.22% assembly success). The root cause was Python's built-in `hash()` function being non-deterministic across invocations for security reasons (PYTHONHASHSEED randomization).

**Root cause**:
```python
seed = hash(seed_string) % (2**32)  # Non-deterministic across runs!
```
Each pipeline run generated different scrambled sequences, causing assembly performance to vary randomly based on which scrambled sequences happened to work better with the kmer assembler.

**Fix**: Use deterministic MD5 hash instead:
```python
import hashlib
seed = int(hashlib.md5(seed_string.encode()).hexdigest()[:8], 16)
```

**Result**: Scrambled sequences are now consistent across runs, making assembly results reproducible.

**Files modified**:
- `masking.py` line 11 (added hashlib import)
- `masking.py` line 32 (changed to deterministic hash)
