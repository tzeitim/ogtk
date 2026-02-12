"""
Sequence masking utilities for handling repetitive sequences in long cassettes.

This module provides functions to mask repetitive TARGET sequences based on
preceding META sequences, which helps the kmer-based assembler (rogtk) handle
long cassettes with direct repeats that exceed the 64-base kmer limit.
"""

import polars as pl
import random
import hashlib
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
    # Use deterministic hash (Python's hash() is randomized for security)
    seed = int(hashlib.md5(seed_string.encode()).hexdigest()[:8], 16)
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
    .. deprecated::
        Use `generate_mask_flanks_expression()` instead, which handles both
        edited and unedited cassettes universally.

    Generate a Polars expression to replace TARGET sequences with scrambled versions
    based on preceding META sequences.

    This function is designed to handle very long cassettes with multiple direct repeats
    that cannot be decomposed by kmer-based assemblers (due to kmer length limits,
    typically 64 bases in rogtk). By replacing repetitive TARGET sequences with
    deterministic scrambled versions based on their preceding META context, the
    assembler can work on unique/variable regions while preserving the structure.

    Expected CSV format (3 columns: feature, seq, kind):
        feature,seq,kind
        META01,AGAAGCCGTGTGCCGGTCTA,META
        META02,ATCGTGCGGACGAGACAGCA,META
        RNF2,TGGCAGTCATCTTAGTCATTACGACAGGTGTTCGTTGTAACTCATATA,TARGET
        HEK3,CTTGGGGCCCAGACTGAGCACGACTTGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCA,TARGET
        EMX1,GGCCTGAGTCCGAGCAGAAGAACTTGGGCTCCCATCACATCAACCGGTGG,TARGET

    How it works:
        - META sequences are anchors that identify context
        - TARGET sequences are the repetitive elements to be masked
        - When pattern (META)(.*?)(TARGET) is found, TARGET is replaced with a
          scrambled version unique to that META+TARGET combination
        - Example: If sequence contains META01 followed by RNF2, the RNF2 portion
          will be replaced with a scrambled sequence deterministically generated
          from hash("META01_RNF2")

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
            # IMPORTANT: Don't apply fuzzy pattern to replacement - it should be plain DNA
            pl.struct(['feature', 'feature_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}",
                    sequence=meta.filter(pl.col('feature') == x['feature_target']).get_column('seq')[0],
                    fuzzy_pattern=False,  # Replacement must be plain DNA, not regex
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_seq")
        )
        .with_columns(
            # Build regex pattern: (META)(.*?)(TARGET)(.*$)
            # Apply fuzzy matching to META and TARGET for better matching
            pattern=pl.concat_str([
                pl.lit("("),
                pl.when(fuzzy_pattern)
                .then(pl.col('seq').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('seq')),  # META sequence (optionally fuzzified)
                pl.lit(")(.*?)("),
                pl.when(fuzzy_pattern)
                .then(pl.col('seq_target').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('seq_target')),  # TARGET sequence (optionally fuzzified)
                pl.lit(")(.*$)")
            ]),
            # Build replacement: $1${2}<scrambled>$4
            # This preserves META, the content between META and TARGET,
            # replaces TARGET with scrambled version, and preserves the rest
            replacement=pl.concat_str([
                pl.lit("${1}${2}"),
                pl.col('scrambled_seq'),  # Plain scrambled DNA sequence
                pl.lit("${4}")
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

    if logger is not None:
        logger.debug(f"Built expression with {len(patterns_df)} chained replacements")
    return expr


def generate_unmask_PEtracer_expression(features_csv: str,
                                        column_name: str = "contig",
                                        fuzzy_pattern: bool = True,
                                        fuzzy_kwargs: Optional[dict] = None,
                                        logger: Optional[CustomLogger] = None
                                        ) -> pl.Expr:
    """
    .. deprecated::
        Use `generate_unmask_flanks_expression()` instead, which handles both
        edited and unedited cassettes universally.

    Generate a Polars expression to restore original TARGET sequences from scrambled versions.

    This reverses the masking operation performed by generate_mask_PEtracer_expression().
    After assembly, contigs contain scrambled TARGET sequences which need to be restored
    to their original form for downstream analysis.

    Args:
        features_csv: Path to CSV file containing feature definitions (same as used for masking)
        column_name: Name of the column to apply replacements to (default: "contig")
        fuzzy_pattern: Whether fuzzy matching was used during masking
        fuzzy_kwargs: Dictionary of keyword arguments used during masking
                     Must match the parameters used in generate_mask_PEtracer_expression

    Returns:
        pl.Expr: Polars expression with chained replacements to restore original sequences

    Usage:
        # Unmask assembled contigs
        df_contigs = (
            df_contigs
            .with_columns(
                generate_unmask_PEtracer_expression('features.csv').alias('contig')
            )
        )

    Note:
        The fuzzy_pattern and fuzzy_kwargs parameters must match what was used during
        masking, otherwise the scrambled sequences won't be recognized.
    """
    if fuzzy_kwargs is None:
        fuzzy_kwargs = {}

    # Read features CSV
    meta = pl.read_csv(features_csv)

    # Create all META-TARGET combinations and generate scrambled sequences
    # This must match exactly what was done during masking
    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'TARGET'),
            suffix="_target",
            how='cross'
        )
        .with_columns(
            # Generate the same scrambled sequence used during masking
            # MUST use fuzzy_pattern=False to get plain DNA (matching mask behavior)
            pl.struct(['feature', 'feature_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}",
                    sequence=meta.filter(pl.col('feature') == x['feature_target']).get_column('seq')[0],
                    fuzzy_pattern=False,  # Must match masking - plain DNA only
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_seq")
        )
        .with_columns(
            # Build regex pattern: (META)(.*?)(SCRAMBLED_TARGET)(.*$)
            # This finds the scrambled sequences in assembled contigs
            # Apply fuzzy matching to META and scrambled sequence for better matching
            pattern=pl.concat_str([
                pl.lit("("),
                pl.when(fuzzy_pattern)
                .then(pl.col('seq').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('seq')),  # META sequence (optionally fuzzified)
                pl.lit(")(.*?)("),
                pl.when(fuzzy_pattern)
                .then(pl.col('scrambled_seq').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('scrambled_seq')),  # Scrambled TARGET (optionally fuzzified)
                pl.lit(")(.*$)")
            ]),
            # Build replacement: $1${2}<ORIGINAL_TARGET>$4
            # This restores META, content between META and TARGET, original TARGET, and the rest
            replacement=pl.concat_str([
                pl.lit("${1}${2}"),
                pl.col('seq_target'),  # Original TARGET sequence
                pl.lit("${4}")
            ])
        )
    )

    # Log pattern generation if logger provided
    if logger is not None:
        logger.debug(f"Generated {len(patterns_df)} unmasking patterns")
        if len(patterns_df) > 0:
            logger.debug("Sample unmasking patterns:")
            sample = patterns_df.head(3).select(['feature', 'feature_target'])
            for row in sample.iter_rows(named=True):
                logger.debug(f"  {row['feature']} + {row['feature_target']}: restoring original TARGET")

    # Build the chained expression by applying all pattern replacements
    expr = pl.col(column_name)

    for row in patterns_df.iter_rows(named=True):
        pattern = row['pattern']
        replacement = row['replacement']
        expr = expr.str.replace(pattern, replacement)

    if logger is not None:
        logger.debug(f"Built unmasking expression with {len(patterns_df)} chained replacements")
    return expr


def generate_mask_flanks_expression(features_csv: str,
                                     column_name: str = "seq",
                                     fuzzy_pattern: bool = True,
                                     fuzzy_kwargs: Optional[dict] = None,
                                     logger: Optional[CustomLogger] = None
                                     ) -> pl.Expr:
    """
    Generate a Polars expression to mask TARGET flanks based on preceding META sequences.

    This function handles both edited and unedited cassettes by masking the flanking
    sequences around potential lineage mark insertion points. Unlike the original
    generate_mask_PEtracer_expression() which masks whole TARGETs, this function:
    - Splits TARGETs into LEFT_FLANK and RIGHT_FLANK
    - Scrambles each flank independently based on META context
    - Preserves any inserted lineage marks (0-5bp between flanks)

    Expected CSV format (6 columns: feature, seq, left_flank, right_flank, marks, kind):
        feature,seq,left_flank,right_flank,marks,kind
        META01,AGAAGCCGTGTGCCGGTCTA,,,META
        META02,ATCGTGCGGACGAGACAGCA,,,META
        RNF2,CATCTTAGTCATTACGACAGGTGTTCGTTG,CATCTTAGTCATTAC,GACAGGTGTTCGTTG,ACTGT|TAAGT|...,EDITED_TARGET
        HEK3,CAGACTGAGCACGACTTGGCAGAGG...,CAGACTGAGCACG,ACTTGGCAGAGG...,CTATC|CGATT|...,EDITED_TARGET

    How it works:
        - META sequences are anchors that identify context
        - LEFT_FLANK and RIGHT_FLANK define the target boundaries
        - Pattern: (META)(.*?)(LEFT_FLANK)(.{0,5}?)(RIGHT_FLANK)(.*$)
        - Group 4 captures: empty (unedited) or 5bp mark (edited)
        - Flanks are scrambled based on META+TARGET combination
        - Marks are preserved (group 4 passes through unchanged)

    Args:
        features_csv: Path to CSV file containing feature definitions
        column_name: Name of the column to apply replacements to (default: "seq")
        fuzzy_pattern: Whether to perform fuzzy string matching to account for sequencing errors
        fuzzy_kwargs: Dictionary of keyword arguments to pass to fuzzy_match_str function

    Returns:
        pl.Expr: Polars expression with chained replacements

    Usage:
        expr = generate_mask_flanks_expression('features_edited.csv')
        result = df.lazy().with_columns(expr.alias("masked_seq")).collect()
    """
    if fuzzy_kwargs is None:
        fuzzy_kwargs = {}

    # Read features CSV
    meta = pl.read_csv(features_csv)

    # Create all META-EDITED_TARGET combinations
    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'EDITED_TARGET'),
            suffix="_target",
            how='cross'
        )
        .with_columns([
            # Generate scrambled sequence for LEFT_FLANK
            pl.struct(['feature', 'feature_target', 'left_flank_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}_L",
                    sequence=x['left_flank_target'],
                    fuzzy_pattern=False,  # Replacement must be plain DNA
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_left"),
            # Generate scrambled sequence for RIGHT_FLANK
            pl.struct(['feature', 'feature_target', 'right_flank_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}_R",
                    sequence=x['right_flank_target'],
                    fuzzy_pattern=False,  # Replacement must be plain DNA
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_right")
        ])
        .with_columns([
            # Build regex pattern: (META)(.*?)(LEFT_FLANK)(.{0,5}?)(RIGHT_FLANK)(.*$)
            # Groups: 1=META, 2=between, 3=LEFT, 4=MARK(0-5bp), 5=RIGHT, 6=rest
            pl.concat_str([
                pl.lit("("),
                pl.when(fuzzy_pattern)
                .then(pl.col('seq').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('seq')),  # META sequence
                pl.lit(")(.*?)("),
                pl.when(fuzzy_pattern)
                .then(pl.col('left_flank_target').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('left_flank_target')),  # LEFT_FLANK
                pl.lit(")(.{0,5}?)("),  # MARK capture (0-5 bases)
                pl.when(fuzzy_pattern)
                .then(pl.col('right_flank_target').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('right_flank_target')),  # RIGHT_FLANK
                pl.lit(")(.*$)")
            ]).alias("pattern"),
            # Build replacement: ${1}${2}<scrambled_left>${4}<scrambled_right>${6}
            # Preserves META, between content, MARK, and rest
            pl.concat_str([
                pl.lit("${1}${2}"),
                pl.col('scrambled_left'),
                pl.lit("${4}"),
                pl.col('scrambled_right'),
                pl.lit("${6}")
            ]).alias("replacement")
        ])
    )

    # Log pattern generation if logger provided
    if logger is not None:
        logger.debug(f"Generated {len(patterns_df)} flank replacement patterns")
        if len(patterns_df) > 0:
            logger.debug("Sample flank patterns:")
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

    if logger is not None:
        logger.debug(f"Built expression with {len(patterns_df)} chained flank replacements")
    return expr


def generate_unmask_flanks_expression(features_csv: str,
                                       column_name: str = "contig",
                                       fuzzy_pattern: bool = True,
                                       fuzzy_kwargs: Optional[dict] = None,
                                       logger: Optional[CustomLogger] = None
                                       ) -> pl.Expr:
    """
    Generate a Polars expression to restore original TARGET flanks from scrambled versions.

    This reverses the masking operation performed by generate_mask_flanks_expression().
    After assembly, contigs contain scrambled flank sequences which need to be restored
    to their original form for downstream analysis. Lineage marks are preserved.

    Args:
        features_csv: Path to CSV file containing feature definitions (same as used for masking)
        column_name: Name of the column to apply replacements to (default: "contig")
        fuzzy_pattern: Whether fuzzy matching was used during masking
        fuzzy_kwargs: Dictionary of keyword arguments used during masking

    Returns:
        pl.Expr: Polars expression with chained replacements to restore original sequences

    Usage:
        df_contigs = df_contigs.with_columns(
            generate_unmask_flanks_expression('features_edited.csv').alias('contig')
        )
    """
    if fuzzy_kwargs is None:
        fuzzy_kwargs = {}

    # Read features CSV
    meta = pl.read_csv(features_csv)

    # Create all META-EDITED_TARGET combinations (must match masking exactly)
    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'EDITED_TARGET'),
            suffix="_target",
            how='cross'
        )
        .with_columns([
            # Generate the same scrambled sequences used during masking
            pl.struct(['feature', 'feature_target', 'left_flank_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}_L",
                    sequence=x['left_flank_target'],
                    fuzzy_pattern=False,
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_left"),
            pl.struct(['feature', 'feature_target', 'right_flank_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string=f"{x['feature']}_{x['feature_target']}_R",
                    sequence=x['right_flank_target'],
                    fuzzy_pattern=False,
                    fuzzy_kwargs=fuzzy_kwargs),
                return_dtype=pl.Utf8
            ).alias("scrambled_right")
        ])
        .with_columns([
            # Build regex pattern: (META)(.*?)(SCRAMBLED_LEFT)(.{0,5}?)(SCRAMBLED_RIGHT)(.*$)
            pl.concat_str([
                pl.lit("("),
                pl.when(fuzzy_pattern)
                .then(pl.col('seq').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('seq')),  # META sequence
                pl.lit(")(.*?)("),
                pl.when(fuzzy_pattern)
                .then(pl.col('scrambled_left').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('scrambled_left')),  # SCRAMBLED_LEFT
                pl.lit(")(.{0,5}?)("),  # MARK capture (preserved)
                pl.when(fuzzy_pattern)
                .then(pl.col('scrambled_right').map_elements(
                    lambda x: fuzzy_match_str(x, **fuzzy_kwargs),
                    return_dtype=pl.Utf8))
                .otherwise(pl.col('scrambled_right')),  # SCRAMBLED_RIGHT
                pl.lit(")(.*$)")
            ]).alias("pattern"),
            # Build replacement: ${1}${2}<original_left>${4}<original_right>${6}
            pl.concat_str([
                pl.lit("${1}${2}"),
                pl.col('left_flank_target'),  # Original LEFT_FLANK
                pl.lit("${4}"),  # Preserve MARK
                pl.col('right_flank_target'),  # Original RIGHT_FLANK
                pl.lit("${6}")
            ]).alias("replacement")
        ])
    )

    # Log pattern generation if logger provided
    if logger is not None:
        logger.debug(f"Generated {len(patterns_df)} flank unmasking patterns")
        if len(patterns_df) > 0:
            logger.debug("Sample unmasking patterns:")
            sample = patterns_df.head(3).select(['feature', 'feature_target'])
            for row in sample.iter_rows(named=True):
                logger.debug(f"  {row['feature']} + {row['feature_target']}: restoring original flanks")

    # Build the chained expression by applying all pattern replacements
    expr = pl.col(column_name)

    for row in patterns_df.iter_rows(named=True):
        pattern = row['pattern']
        replacement = row['replacement']
        expr = expr.str.replace(pattern, replacement)

    if logger is not None:
        logger.debug(f"Built unmasking expression with {len(patterns_df)} chained flank replacements")
    return expr
