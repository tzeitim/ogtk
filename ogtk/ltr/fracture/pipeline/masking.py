"""
Sequence masking utilities for handling repetitive sequences in long cassettes.

This module provides functions to mask repetitive TARGET sequences based on
preceding META sequences, which helps the kmer-based assembler (rogtk) handle
long cassettes with direct repeats that exceed the 64-base kmer limit.
"""

import polars as pl
import random
from typing import Optional
from ogtk.utils.log import CustomLogger
from ogtk.utils.general import fuzzy_match_str


def scramble_dna_sequence(seed_string: str, sequence: str, fuzzy_pattern: bool,
                          fuzzy_kwargs: dict) -> str:
    """
    Scramble a DNA sequence using a deterministic seed.

    Args:
        seed_string: String to use as seed for random number generator
        sequence: DNA sequence to scramble
        fuzzy_pattern: Whether to apply fuzzy matching to scrambled sequence
        fuzzy_kwargs: Keyword arguments to pass to fuzzy_match_str

    Returns:
        Scrambled DNA sequence (optionally with fuzzy pattern)
    """
    seed = hash(seed_string) % (2**32)  # Ensure positive 32-bit integer
    random.seed(seed)
    seq_list = list(sequence.upper())
    random.shuffle(seq_list)
    scrambled_seq = ''.join(seq_list)

    if fuzzy_pattern:
        return fuzzy_match_str(scrambled_seq, **fuzzy_kwargs)
    return scrambled_seq


def generate_mask_PEtracer_expression(features_csv: str,
                                      column_name: str = "seq",
                                      fuzzy_pattern: bool = True,
                                      fuzzy_kwargs: Optional[dict] = None,
                                      logger: Optional[CustomLogger] = None
                                      ) -> pl.Expr:
    """
    Generate a Polars expression to replace TARGET sequences with scrambled versions
    based on preceding META sequences.

    This function is designed to handle very long cassettes with multiple direct repeats
    that cannot be decomposed by kmer-based assemblers (due to kmer length limits,
    typically 64 bases in rogtk). By replacing repetitive TARGET sequences with
    deterministic scrambled versions based on their preceding META context, the
    assembler can work on unique/variable regions while preserving the structure.

    Expected CSV format:
        feature,seq,kind
        META01,AGAAGCCGTGTGCCGGTCTA,META
        META02,ATCGTGCGGACGAGACAGCA,META
        RNF2,TGGCAGTCATCTTAGTCATTACGACAGGTGTTCGTTGTAACTCATATA,TARGET
        HEK3,CTTGGGGCCCAGACTGAGCACGACTTGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCA,TARGET
        EMX1,GGCCTGAGTCCGAGCAGAAGAACTTGGGCTCCCATCACATCAACCGGTGG,TARGET

    Args:
        features_csv: Path to CSV file containing feature definitions
        column_name: Name of the column to apply replacements to (default: "seq")
        fuzzy_pattern: Whether to perform fuzzy string matching to account for sequencing errors
        fuzzy_kwargs: Dictionary of keyword arguments to pass to fuzzy_match_str function
                     Supported keys: wildcard, include_original, sep, max_length

    Returns:
        pl.Expr: Polars expression with chained replacements

    Usage:
        # Basic usage
        expr = generate_mask_PEtracer_expression('features.csv')

        # With custom fuzzy matching parameters
        expr = generate_mask_PEtracer_expression(
            'features.csv',
            fuzzy_pattern=True,
            fuzzy_kwargs={
                'wildcard': '.{0,2}',      # Allow up to 2 character variations
                'include_original': False,  # Don't include exact match
                'sep': '|',                # Use pipe separator
                'max_length': 150          # Only fuzzify sequences up to 150 chars
            }
        )

        # Apply to a DataFrame
        result = (
            df.lazy()
            .with_columns(expr.alias("masked_seq"))
            .collect()
        )

    Example:
        If your cassette has: META01 -> RNF2 -> META02 -> RNF2
        The second RNF2 will be replaced with a scrambled version (deterministic
        based on META02+RNF2), making it distinguishable for the assembler while
        preserving the information about what was originally there.
    """
    if fuzzy_kwargs is None:
        fuzzy_kwargs = {}

    # Read features CSV
    meta = pl.read_csv(features_csv)

    # Create all META-TARGET combinations and generate scrambled sequences
    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'TARGET'),
            suffix="_target",
            how='cross'
        )
        .with_columns(
            # Generate scrambled sequence for each META-TARGET combination
            pl.struct(['feature', 'feature_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}",
                    sequence=meta.filter(pl.col('feature') == x['feature_target']).get_column('seq')[0],
                    fuzzy_pattern=fuzzy_pattern,
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_seq")
        )
        .with_columns(
            # Build regex pattern: (META)(.*?)(TARGET)(.*$)
            pattern=pl.concat_str([
                pl.lit("("),
                pl.col('seq'),  # META sequence
                pl.lit(")(.*?)("),
                pl.col('seq_target'),  # TARGET sequence
                pl.lit(")(.*$)")
            ]),
            # Build replacement: $1${2}<scrambled>$4
            # This preserves META, the content between META and TARGET,
            # replaces TARGET with scrambled version, and preserves the rest
            replacement=pl.concat_str([
                pl.lit("$1${2}"),
                pl.col('scrambled_seq'),
                pl.lit("$4")
            ])
        )
    )

    # Log pattern generation if logger provided
    if logger is not None:
        logger.debug(f"Generated {len(patterns_df)} replacement patterns")
        if len(patterns_df) > 0:
            logger.debug("Sample patterns:")
            sample = patterns_df.head(3).select(['feature', 'feature_target', 'pattern', 'replacement'])
            for row in sample.iter_rows(named=True):
                logger.debug(f"  {row['feature']} + {row['feature_target']}: {row['pattern'][:60]}...")
                logger.debug(f"    -> {row['replacement'][:60]}...")

    # Build the chained expression by applying all pattern replacements
    expr = pl.col(column_name)

    for row in patterns_df.iter_rows(named=True):
        pattern = row['pattern']
        replacement = row['replacement']
        expr = expr.str.replace(pattern, replacement)

    print(f"Built expression with {len(patterns_df)} chained replacements")
    return expr
