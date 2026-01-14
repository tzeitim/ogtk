"""
Segmentation utilities for splitting long cassette reads at META boundaries.

Core functions:
- segment_by_metas: Sanitize + segment reads at META boundaries + extract edges
- stitch_segments: Rejoin assembled segments into full sequences

All functions operate on LazyFrames for optimal performance with large datasets.
"""

import polars as pl
from typing import Optional
from ogtk.utils.general import fuzzy_match_str
from ogtk.utils.log import CustomLogger


# Default META length (all PEtracer METAs are 20bp)
META_LEN = 20

# Special markers for cassette edge segments
CASSETTE_START_MARKER = "_CASSETTE_START_"
CASSETTE_END_MARKER = "_CASSETTE_END_"


def segment_by_metas(
    ldf: pl.LazyFrame,
    metas: pl.DataFrame,
    cassette_start_anchor: str,
    cassette_end_anchor: str,
    seq_col: str = 'r2_seq',
    keep_cols: Optional[list] = None,
    wildcard: str | None = ".{0,1}",
    sanitize: bool = True,
    logger: Optional[CustomLogger] = None
) -> pl.LazyFrame:
    """
    Split reads at META boundaries into overlapping segments.

    This function performs the complete segmentation pipeline:
    1. Sanitize METAs and cassette anchors (normalize fuzzy matches to canonical sequences)
    2. Extract META-to-META segments
    3. Extract edge segments (cassette start/end)

    Each segment includes both boundary METAs for overlap-based stitching.

    Args:
        ldf: LazyFrame with sequence data
        metas: DataFrame with 'feature', 'seq', and optionally 'kind' columns
        cassette_start_anchor: Sequence marking cassette start (required)
        cassette_end_anchor: Sequence marking cassette end (required)
        seq_col: Name of the column containing sequences
        keep_cols: Additional columns to preserve (default: ['umi'])
        wildcard: Regex pattern for fuzzy matching (default: ".{0,1}" allows 1 error).
                  Set to None for exact matching only.
        sanitize: Whether to sanitize META sequences first (default: True)
        logger: Optional logger for debugging

    Returns:
        LazyFrame with columns: [keep_cols], start_meta, end_meta, segment_seq
    """
    if keep_cols is None:
        keep_cols = ['umi']

    # Filter to only META sequences
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_seqs = meta_only['seq'].to_list()
    meta_names = meta_only['feature'].to_list()

    # ==================== SANITIZE ====================
    # Replace fuzzy matches with canonical sequences for METAs and cassette anchors
    if sanitize and wildcard:
        if logger:
            logger.debug(f"Sanitizing {len(meta_seqs)} META sequences + 2 cassette anchors with wildcard: {wildcard}")

        # Sanitize METAs
        for seq, name in zip(meta_seqs, meta_names):
            fuzzy_pattern = fuzzy_match_str(seq, wildcard=wildcard, include_original=True)
            ldf = ldf.with_columns(pl.col(seq_col).str.replace_all(fuzzy_pattern, seq))

        # Sanitize cassette anchors
        fuzzy_start = fuzzy_match_str(cassette_start_anchor, wildcard=wildcard, include_original=True)
        ldf = ldf.with_columns(pl.col(seq_col).str.replace_all(fuzzy_start, cassette_start_anchor))
        fuzzy_end = fuzzy_match_str(cassette_end_anchor, wildcard=wildcard, include_original=True)
        ldf = ldf.with_columns(pl.col(seq_col).str.replace_all(fuzzy_end, cassette_end_anchor))

    # Keep reference for edge extraction (after sanitization)
    sanitized_ldf = ldf

    if logger:
        logger.info(f"Segmenting by {len(meta_seqs)} META sequences")

    # ==================== STEP 2: META-TO-META SEGMENTS ====================
    # Pre-compute positions for all METAs
    pos_exprs = [
        pl.col(seq_col).str.find(seq).alias(f'_pos_{name}')
        for name, seq in zip(meta_names, meta_seqs)
    ]

    ldf = (
        ldf
        .with_columns(pl.col(seq_col).str.extract_many(meta_seqs).alias('meta_matches'))
        .with_columns(pos_exprs)
    )

    # Create mapping DataFrame for seq -> name conversion
    seq_name_map = pl.DataFrame({
        'meta_seq': meta_seqs,
        'meta_name': meta_names,
    }).lazy()

    # Add row index to preserve order
    ldf = ldf.with_row_index('_row_idx')

    # Explode, join, and create consecutive pairs
    ldf = (
        ldf
        .explode('meta_matches')
        .join(seq_name_map, left_on='meta_matches', right_on='meta_seq', how='left')
        .with_columns(
            pl.col('meta_name').cum_count().over('_row_idx').alias('_meta_idx')
        )
    )

    # Create consecutive pairs by self-join
    ldf_next = (
        ldf
        .select(
            pl.col('_row_idx'),
            (pl.col('_meta_idx') - 1).alias('_meta_idx'),
            pl.col('meta_name').alias('end_meta'),
        )
    )

    ldf = (
        ldf
        .join(ldf_next, on=['_row_idx', '_meta_idx'], how='left')
        .filter(pl.col('end_meta').is_not_null())
        .rename({'meta_name': 'start_meta'})
    )

    # Build position expressions
    start_pos_expr = pl.lit(None).cast(pl.UInt32)
    end_pos_expr = pl.lit(None).cast(pl.UInt32)
    for name in meta_names:
        start_pos_expr = (
            pl.when(pl.col('start_meta') == name)
            .then(pl.col(f'_pos_{name}'))
            .otherwise(start_pos_expr)
        )
        end_pos_expr = (
            pl.when(pl.col('end_meta') == name)
            .then(pl.col(f'_pos_{name}'))
            .otherwise(end_pos_expr)
        )

    # Extract segments
    ldf = (
        ldf
        .with_columns(
            start_pos_expr.alias('start_pos'),
            end_pos_expr.alias('end_pos'),
        )
        .with_columns(
            pl.col(seq_col).str.slice(
                pl.col('start_pos'),
                pl.col('end_pos') - pl.col('start_pos') + META_LEN
            ).alias('segment_seq')
        )
    )

    meta_segments = ldf.select(keep_cols + ['start_meta', 'end_meta', 'segment_seq'])

    # ==================== EDGE SEGMENTS ====================
    # Empirically determine first/last META for each read (handles variable cassette structures)
    if logger:
        logger.debug(f"Extracting edge segments (empirical META detection)")

    # Build position expressions for all METAs
    meta_pos_exprs = [
        pl.col(seq_col).str.find(seq).alias(f'_pos_{name}')
        for name, seq in zip(meta_names, meta_seqs)
    ]

    edge_ldf = sanitized_ldf.with_columns([
        pl.col(seq_col).str.find(cassette_start_anchor).alias('_start_anchor_pos'),
        pl.col(seq_col).str.find(cassette_end_anchor).alias('_end_anchor_pos'),
    ] + meta_pos_exprs)

    # Start edge: cassette_start_anchor -> first META after it
    # Find the META with minimum position that is > start_anchor_pos
    min_meta_after_start = pl.min_horizontal(*[
        pl.when(pl.col(f'_pos_{name}').is_not_null() & (pl.col(f'_pos_{name}') > pl.col('_start_anchor_pos')))
        .then(pl.col(f'_pos_{name}'))
        .otherwise(pl.lit(None).cast(pl.UInt32))
        for name in meta_names
    ])

    # Find which META corresponds to that min position
    first_meta_name_expr = pl.lit(None).cast(pl.Utf8)
    for name in meta_names:
        first_meta_name_expr = (
            pl.when(pl.col(f'_pos_{name}') == pl.col('_first_meta_pos'))
            .then(pl.lit(name))
            .otherwise(first_meta_name_expr)
        )

    start_edge = (
        edge_ldf
        .with_columns(min_meta_after_start.alias('_first_meta_pos'))
        .with_columns(first_meta_name_expr.alias('_first_meta_name'))
        .filter(
            pl.col('_start_anchor_pos').is_not_null() &
            pl.col('_first_meta_pos').is_not_null()
        )
        .with_columns([
            pl.lit(CASSETTE_START_MARKER).alias('start_meta'),
            pl.col('_first_meta_name').alias('end_meta'),
            pl.col(seq_col).str.slice(
                pl.col('_start_anchor_pos'),
                pl.col('_first_meta_pos') - pl.col('_start_anchor_pos') + META_LEN
            ).alias('segment_seq'),
        ])
        .select(keep_cols + ['start_meta', 'end_meta', 'segment_seq'])
    )

    # End edge: last META before cassette_end_anchor -> end_anchor
    # Find the META with maximum position that is < end_anchor_pos
    max_meta_before_end = pl.max_horizontal(*[
        pl.when(pl.col(f'_pos_{name}').is_not_null() & (pl.col(f'_pos_{name}') < pl.col('_end_anchor_pos')))
        .then(pl.col(f'_pos_{name}'))
        .otherwise(pl.lit(None).cast(pl.UInt32))
        for name in meta_names
    ])

    # Find which META corresponds to that max position
    last_meta_name_expr = pl.lit(None).cast(pl.Utf8)
    for name in meta_names:
        last_meta_name_expr = (
            pl.when(pl.col(f'_pos_{name}') == pl.col('_last_meta_pos'))
            .then(pl.lit(name))
            .otherwise(last_meta_name_expr)
        )

    end_anchor_len = len(cassette_end_anchor)
    end_edge = (
        edge_ldf
        .with_columns(max_meta_before_end.alias('_last_meta_pos'))
        .with_columns(last_meta_name_expr.alias('_last_meta_name'))
        .filter(
            pl.col('_end_anchor_pos').is_not_null() &
            pl.col('_last_meta_pos').is_not_null()
        )
        .with_columns([
            pl.col('_last_meta_name').alias('start_meta'),
            pl.lit(CASSETTE_END_MARKER).alias('end_meta'),
            pl.col(seq_col).str.slice(
                pl.col('_last_meta_pos'),
                pl.col('_end_anchor_pos') - pl.col('_last_meta_pos') + end_anchor_len
            ).alias('segment_seq'),
        ])
        .select(keep_cols + ['start_meta', 'end_meta', 'segment_seq'])
    )

    return pl.concat([meta_segments, start_edge, end_edge])


def stitch_segments(
    ldf: pl.LazyFrame,
    metas: pl.DataFrame,
    seq_col: str = 'consensus_seq',
    cassette_start_anchor_len: int | None = None,
    group_cols: list[str] | None = None,
    logger: Optional[CustomLogger] = None
) -> pl.LazyFrame:
    """
    Rejoin assembled segments into full sequences.

    Args:
        ldf: LazyFrame with assembled segments (umi, start_meta, end_meta, consensus_seq)
        metas: DataFrame with 'feature' for ordering
        seq_col: Name of the column containing consensus sequences
        cassette_start_anchor_len: Length of cassette start anchor (for edge segment trimming)
        group_cols: Columns to group by for stitching (default: ['umi'], can include 'sbc')
        logger: Optional logger for debugging

    Returns:
        LazyFrame with columns: [group_cols], stitched_seq
    """
    if group_cols is None:
        group_cols = ['umi']

    # Filter to only META sequences for ordering
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas

    # Build META order map (CASSETTE_START_MARKER gets order -1)
    meta_names = meta_only['feature'].to_list()
    meta_order = {CASSETTE_START_MARKER: -1}
    meta_order.update({m: i for i, m in enumerate(meta_names)})

    if logger:
        logger.info(f"Stitching segments using {len(meta_order)} META order positions")
        logger.info(f"Grouping by: {group_cols}")

    # Build sort order expression
    sort_order_expr = pl.lit(None).cast(pl.Int64)
    for name, order in meta_order.items():
        sort_order_expr = (
            pl.when(pl.col('start_meta') == name)
            .then(pl.lit(order))
            .otherwise(sort_order_expr)
        )

    # Trimming: first segment keeps full, others trim META_LEN from start
    trim_expr = (
        pl.when(pl.col('sort_order') == pl.col('sort_order').min().over(group_cols))
        .then(pl.col(seq_col))
        .otherwise(pl.col(seq_col).str.slice(META_LEN))
    )

    return (
        ldf
        .with_columns(sort_order_expr.alias('sort_order'))
        .sort(*group_cols, 'sort_order')
        .with_columns(trim_expr.alias('contribution'))
        .group_by(group_cols).agg(
            pl.col('contribution').str.concat('').alias('stitched_seq')
        )
    )


def flag_aberrant_molecules(
    ldf: pl.LazyFrame,
    metas: pl.DataFrame,
    seq_col: str = 'r2_seq',
) -> pl.LazyFrame:
    """
    Flag molecules with tandem duplications (any META appearing > 1 time).

    Args:
        ldf: LazyFrame with sequence data
        metas: DataFrame with 'feature' and 'seq' columns
        seq_col: Name of the column containing sequences

    Returns:
        LazyFrame with added 'is_aberrant' boolean column
    """
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas

    count_exprs = [
        pl.col(seq_col).str.count_matches(pl.lit(seq)).alias(f'{name}_count')
        for name, seq in zip(meta_only['feature'], meta_only['seq'])
    ]

    ldf = ldf.with_columns(count_exprs)

    count_cols = [f'{name}_count' for name in meta_only['feature']]
    ldf = ldf.with_columns(
        pl.any_horizontal(pl.col(c) > 1 for c in count_cols).alias('is_aberrant')
    ).drop(count_cols)

    return ldf


def generate_segmentation_report(
    segments_df: pl.DataFrame,
    assembled_df: pl.DataFrame,
    contigs_df: pl.DataFrame,
    metas: pl.DataFrame,
    cassette_start_anchor: str | None = None,
    cassette_end_anchor: str | None = None,
    heterogeneity_threshold: float = 0.20,
) -> dict:
    """Generate a comprehensive diagnostic report for segmentation-based assembly.

    Args:
        heterogeneity_threshold: For smart mode heterogeneity detection, a transition type
            is considered "real" if its frequency is >= this fraction of the dominant type's
            frequency. Default 0.20 (20%). Set to 0 to disable filtering (sensitive mode only).
    """
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_names = meta_only['feature'].to_list()
    meta_seqs = dict(zip(meta_only['feature'], meta_only['seq']))

    report = {
        'segments': {},
        'assembly': {},
        'stitching': {},
        'meta_presence': {},
        'summary': {},
    }

    # Segment stats
    total_segments = segments_df.height
    unique_umis_with_segments = segments_df.select('umi').n_unique()

    # Compute segment_types at MOLECULE level (deduplicate by UMI first)
    # Each UMI contributes at most once per transition type
    umi_transitions = (
        segments_df
        .select(['umi', 'start_meta', 'end_meta', 'segment_seq'])
        .unique(subset=['umi', 'start_meta', 'end_meta'])  # One count per UMI per transition
    )

    segment_types = (
        umi_transitions
        .group_by(['start_meta', 'end_meta'])
        .agg([
            pl.len().alias('count'),  # Now counts molecules, not reads
            pl.col('segment_seq').str.len_chars().mean().alias('mean_length'),
            pl.col('segment_seq').str.len_chars().median().alias('median_length'),
        ])
        .sort('count', descending=True)
    )

    # UMI heterogeneity detection: SENSITIVE vs SMART mode
    # Sensitive: counts ALL unique transition types per UMI (current behavior)
    # Smart: filters out low-frequency transitions (noise) based on adaptive threshold

    # 1. Get counts per transition per UMI (how many reads support each transition)
    umi_transition_counts = (
        segments_df
        .group_by(['umi', 'start_meta', 'end_meta'])
        .len()
        .rename({'len': 'count'})
    )

    # 2. Add total per UMI, frequency, and dominant type frequency
    umi_transition_counts = (
        umi_transition_counts
        .with_columns(
            pl.col('count').sum().over('umi').alias('umi_total'),
            pl.col('count').max().over('umi').alias('umi_max_count'),
        )
        .with_columns(
            (pl.col('count') / pl.col('umi_total')).alias('freq'),
            (pl.col('umi_max_count') / pl.col('umi_total')).alias('dominant_freq'),
        )
    )

    # 3. SENSITIVE mode: all transitions (current behavior)
    transitions_per_umi_sensitive = (
        umi_transition_counts
        .group_by('umi')
        .agg(pl.len().alias('n_transition_types'))
    )

    umi_transition_dist_sensitive = (
        transitions_per_umi_sensitive
        .group_by('n_transition_types')
        .len()
        .sort('n_transition_types')
        .to_dicts()
    )

    umis_single_type_sensitive = transitions_per_umi_sensitive.filter(pl.col('n_transition_types') == 1).height
    umis_multi_type_sensitive = transitions_per_umi_sensitive.filter(pl.col('n_transition_types') > 1).height
    max_types_per_umi_sensitive = transitions_per_umi_sensitive['n_transition_types'].max()

    # 4. SMART mode: adaptive threshold based on dominant type
    # A type is "real" if freq >= heterogeneity_threshold * dominant_freq
    umi_transition_counts_filtered = (
        umi_transition_counts
        .with_columns(
            (pl.col('dominant_freq') * heterogeneity_threshold).alias('adaptive_threshold')
        )
        .filter(pl.col('freq') >= pl.col('adaptive_threshold'))
    )

    transitions_per_umi_smart = (
        umi_transition_counts_filtered
        .group_by('umi')
        .agg(pl.len().alias('n_transition_types'))
    )

    umi_transition_dist_smart = (
        transitions_per_umi_smart
        .group_by('n_transition_types')
        .len()
        .sort('n_transition_types')
        .to_dicts()
    )

    umis_single_type_smart = transitions_per_umi_smart.filter(pl.col('n_transition_types') == 1).height
    umis_multi_type_smart = transitions_per_umi_smart.filter(pl.col('n_transition_types') > 1).height
    max_types_per_umi_smart = transitions_per_umi_smart['n_transition_types'].max()

    expected_pairs = [(meta_names[i], meta_names[i+1]) for i in range(len(meta_names)-1)]
    consecutive_segment_counts = {}
    for start, end in expected_pairs:
        count = segments_df.filter(
            (pl.col('start_meta') == start) & (pl.col('end_meta') == end)
        ).height
        consecutive_segment_counts[f'{start}->{end}'] = count

    edge_start_count = segments_df.filter(pl.col('start_meta') == CASSETTE_START_MARKER).height
    edge_end_count = segments_df.filter(pl.col('end_meta') == CASSETTE_END_MARKER).height

    report['segments'] = {
        'total_segments': total_segments,
        'unique_umis': unique_umis_with_segments,
        'segment_types': segment_types.to_dicts(),
        'consecutive_pairs': consecutive_segment_counts,
        'edge_start_segments': edge_start_count,
        'edge_end_segments': edge_end_count,
        # UMI collision/chimera detection - SMART mode (filtered, default)
        'heterogeneity_threshold': heterogeneity_threshold,
        'umi_transition_distribution': umi_transition_dist_smart,
        'umis_single_transition_type': umis_single_type_smart,
        'umis_multi_transition_type': umis_multi_type_smart,
        'max_transition_types_per_umi': max_types_per_umi_smart,
        # UMI collision/chimera detection - SENSITIVE mode (unfiltered, for debugging)
        'umi_transition_distribution_sensitive': umi_transition_dist_sensitive,
        'umis_single_transition_type_sensitive': umis_single_type_sensitive,
        'umis_multi_transition_type_sensitive': umis_multi_type_sensitive,
        'max_transition_types_per_umi_sensitive': max_types_per_umi_sensitive,
    }

    # Assembly stats
    total_assembled = assembled_df.height
    unique_umis_assembled = assembled_df.select('umi').n_unique()

    assembly_success = {}
    for row in segment_types.iter_rows(named=True):
        key = f"{row['start_meta']}->{row['end_meta']}"
        input_count = row['count']
        input_groups = segments_df.filter(
            (pl.col('start_meta') == row['start_meta']) &
            (pl.col('end_meta') == row['end_meta'])
        ).select('umi').n_unique()

        output_count = assembled_df.filter(
            (pl.col('start_meta') == row['start_meta']) &
            (pl.col('end_meta') == row['end_meta'])
        ).height

        success_rate = output_count / input_groups * 100 if input_groups > 0 else 0
        assembly_success[key] = {
            'input_segments': input_count,
            'input_umi_groups': input_groups,
            'assembled': output_count,
            'success_rate_pct': round(success_rate, 1),
        }

    report['assembly'] = {
        'total_assembled': total_assembled,
        'unique_umis': unique_umis_assembled,
        'by_type': assembly_success,
    }

    # Stitching stats
    total_contigs = contigs_df.height
    contig_lengths = contigs_df.select(pl.col('contig').str.len_chars().alias('length'))

    length_stats = {
        'count': total_contigs,
        'mean': round(contig_lengths['length'].mean() or 0, 1),
        'median': contig_lengths['length'].median() or 0,
        'min': contig_lengths['length'].min() or 0,
        'max': contig_lengths['length'].max() or 0,
        'std': round(contig_lengths['length'].std() or 0, 1),
    }

    segments_per_umi = assembled_df.group_by('umi').len().select('len')
    segments_per_umi_stats = {
        'mean': round(segments_per_umi['len'].mean() or 0, 1),
        'median': segments_per_umi['len'].median() or 0,
        'min': segments_per_umi['len'].min() or 0,
        'max': segments_per_umi['len'].max() or 0,
    }

    report['stitching'] = {
        'total_contigs': total_contigs,
        'contig_length': length_stats,
        'segments_per_umi': segments_per_umi_stats,
    }

    # META presence
    meta_presence = {}
    for meta_name in meta_names:
        meta_seq = meta_seqs[meta_name]
        count = contigs_df.filter(pl.col('contig').str.contains(meta_seq)).height
        pct = count / total_contigs * 100 if total_contigs > 0 else 0
        meta_presence[meta_name] = {'count': count, 'percent': round(pct, 1)}

    if cassette_start_anchor:
        count = contigs_df.filter(pl.col('contig').str.contains(cassette_start_anchor)).height
        pct = count / total_contigs * 100 if total_contigs > 0 else 0
        meta_presence['CASSETTE_START'] = {'count': count, 'percent': round(pct, 1)}

    if cassette_end_anchor:
        count = contigs_df.filter(pl.col('contig').str.contains(cassette_end_anchor)).height
        pct = count / total_contigs * 100 if total_contigs > 0 else 0
        meta_presence['CASSETTE_END'] = {'count': count, 'percent': round(pct, 1)}

    report['meta_presence'] = meta_presence

    # Summary
    avg_meta_presence = sum(v['percent'] for v in meta_presence.values()) / len(meta_presence) if meta_presence else 0

    report['summary'] = {
        'input_segments': total_segments,
        'assembled_segments': total_assembled,
        'output_contigs': total_contigs,
        'mean_contig_length': length_stats['mean'],
        'median_contig_length': length_stats['median'],
        'avg_meta_presence_pct': round(avg_meta_presence, 1),
        'uses_edge_segments': cassette_start_anchor is not None or cassette_end_anchor is not None,
    }

    return report


def format_segmentation_report(report: dict) -> str:
    """Format the segmentation report as a human-readable string."""
    lines = []
    lines.append("=" * 60)
    lines.append("SEGMENTATION ASSEMBLY REPORT")
    lines.append("=" * 60)

    s = report['summary']
    lines.append("\n## Summary")
    lines.append(f"  Input segments:        {s['input_segments']:,}")
    lines.append(f"  Assembled segments:    {s['assembled_segments']:,}")
    lines.append(f"  Output contigs:        {s['output_contigs']:,}")
    lines.append(f"  Mean contig length:    {s['mean_contig_length']:.0f} bp")
    lines.append(f"  Median contig length:  {s['median_contig_length']:.0f} bp")
    lines.append(f"  Avg META presence:     {s['avg_meta_presence_pct']:.1f}%")
    if s['uses_edge_segments']:
        lines.append(f"  Edge segments:         Enabled")

    seg = report['segments']
    lines.append("\n## Segment Extraction")
    lines.append(f"  Total segments:        {seg['total_segments']:,}")
    lines.append(f"  Unique UMIs:           {seg['unique_umis']:,}")
    if seg['edge_start_segments'] > 0:
        lines.append(f"  Start edge segments:   {seg['edge_start_segments']:,}")
    if seg['edge_end_segments'] > 0:
        lines.append(f"  End edge segments:     {seg['edge_end_segments']:,}")

    lines.append("\n  Top segment types (molecule-level):")
    for item in seg['segment_types'][:5]:
        lines.append(f"    {item['start_meta']:>18} -> {item['end_meta']:<10}: {item['count']:>8,} molecules (median {item['median_length']:.0f}bp)")

    # UMI collision/chimera detection - SMART mode (filtered)
    threshold_pct = seg.get('heterogeneity_threshold', 0.20) * 100
    lines.append(f"\n  UMI transition heterogeneity [Smart mode - {threshold_pct:.0f}% of dominant]:")
    lines.append(f"    UMIs with single transition type:   {seg['umis_single_transition_type']:>8,}")
    lines.append(f"    UMIs with multiple transition types:{seg['umis_multi_transition_type']:>8,}")
    total_umis_smart = seg['umis_single_transition_type'] + seg['umis_multi_transition_type']
    if seg['umis_multi_transition_type'] > 0 and total_umis_smart > 0:
        multi_pct = seg['umis_multi_transition_type'] / total_umis_smart * 100
        lines.append(f"    Multi-type rate:                    {multi_pct:>7.1f}%")
        lines.append(f"    Max transition types per UMI:       {seg['max_transition_types_per_umi']:>8}")

    # UMI collision/chimera detection - SENSITIVE mode (unfiltered, for comparison)
    if 'umis_multi_transition_type_sensitive' in seg:
        lines.append(f"  [Sensitive mode - all transitions]:")
        lines.append(f"    UMIs with multiple (unfiltered):    {seg['umis_multi_transition_type_sensitive']:>8,}")
        total_umis_sensitive = seg['umis_single_transition_type_sensitive'] + seg['umis_multi_transition_type_sensitive']
        if total_umis_sensitive > 0:
            multi_pct_sensitive = seg['umis_multi_transition_type_sensitive'] / total_umis_sensitive * 100
            lines.append(f"    Multi-type rate (unfiltered):       {multi_pct_sensitive:>7.1f}%")

    asm = report['assembly']
    lines.append("\n## Assembly")
    lines.append(f"  Total assembled:       {asm['total_assembled']:,}")
    lines.append(f"  Unique UMIs:           {asm['unique_umis']:,}")

    lines.append("\n  Assembly success by type:")
    for seg_type, stats in sorted(asm['by_type'].items(), key=lambda x: -x[1]['assembled']):
        if stats['assembled'] > 0:
            lines.append(f"    {seg_type:>22}: {stats['assembled']:>6,} / {stats['input_umi_groups']:>6,} ({stats['success_rate_pct']:>5.1f}%)")

    meta = report['meta_presence']
    lines.append("\n## META Presence in Contigs")
    for meta_name, stats in meta.items():
        bar_len = int(stats['percent'] / 5)
        bar = "█" * bar_len + "░" * (20 - bar_len)
        lines.append(f"  {meta_name:>15}: {bar} {stats['percent']:>5.1f}% ({stats['count']:,})")

    lines.append("\n" + "=" * 60)

    return "\n".join(lines)


# ==================== QC PLOTTING FUNCTIONS ====================

def plot_meta_transition_heatmap(
    segment_types: pl.DataFrame,
    metas: pl.DataFrame,
    ax=None,
    title: str = "Meta Transition Matrix",
    logger: Optional[CustomLogger] = None,
):
    """
    Plot NxN heatmap of meta->meta transition percentages.

    Args:
        segment_types: DataFrame with start_meta, end_meta, count columns
        metas: DataFrame with feature column for meta ordering
        ax: matplotlib axes (optional)
        title: Plot title
        logger: Optional logger

    Returns:
        matplotlib axes
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LogNorm
    import numpy as np

    # Build ordered meta list including edge markers
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_order = [CASSETTE_START_MARKER] + meta_only['feature'].to_list() + [CASSETTE_END_MARKER]

    # Ensure we have start_meta, end_meta, count columns
    if 'count' not in segment_types.columns:
        # segment_types might be from report dict format
        segment_types = pl.DataFrame(segment_types)

    # Calculate total for percentage conversion
    total_count = segment_types['count'].sum()

    # Create pivot table
    pivot = (
        segment_types
        .select(['start_meta', 'end_meta', 'count'])
        .pivot(on='end_meta', index='start_meta', values='count')
        .fill_null(0)
    )

    # Reorder rows and columns to match meta order
    existing_starts = [m for m in meta_order if m in pivot['start_meta'].to_list()]
    existing_ends = [m for m in meta_order if m in pivot.columns and m != 'start_meta']

    # Build matrix in correct order
    matrix_data = []
    row_labels = []
    for start in existing_starts:
        row = pivot.filter(pl.col('start_meta') == start)
        if row.height > 0:
            row_values = [row[end].item() if end in row.columns else 0 for end in existing_ends]
            matrix_data.append(row_values)
            row_labels.append(start)

    if not matrix_data:
        if logger:
            logger.warning("No transition data to plot")
        return ax

    matrix = np.array(matrix_data, dtype=float)

    # Convert to percentages
    if total_count > 0:
        matrix_pct = (matrix / total_count) * 100
    else:
        matrix_pct = matrix

    # Create axes if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))

    # Use log scale for better visualization
    # Add small value to avoid log(0)
    matrix_plot = matrix_pct + 0.001

    sns.heatmap(
        matrix_plot,
        annot=False,
        cmap='YlOrRd',
        norm=LogNorm(vmin=0.001, vmax=max(matrix_plot.max(), 1)),
        xticklabels=existing_ends,
        yticklabels=row_labels,
        ax=ax,
        cbar_kws={'label': '% of transitions (log scale)'}
    )

    ax.set_xlabel('End Meta')
    ax.set_ylabel('Start Meta')
    ax.set_title(title)

    # Rotate labels for readability
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    return ax


def plot_meta_transition_bars(
    segment_types: pl.DataFrame,
    metas: pl.DataFrame,
    ax=None,
    title: str = "Transition Frequencies by Pair",
    top_n: int = 30,
    logger: Optional[CustomLogger] = None,
):
    """
    Bar chart showing transition percentages for each meta pair.
    Color-coded by skip distance: green=consecutive, yellow=1-skip, orange=2-skip, red=3+ skip

    Args:
        segment_types: DataFrame with start_meta, end_meta, count columns
        metas: DataFrame with feature column for meta ordering
        ax: matplotlib axes (optional)
        title: Plot title
        top_n: Number of top transitions to show
        logger: Optional logger

    Returns:
        matplotlib axes
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Build meta index for calculating skip distance
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_names = meta_only['feature'].to_list()
    meta_idx = {m: i for i, m in enumerate(meta_names)}

    # Add special markers
    meta_idx[CASSETTE_START_MARKER] = -1
    meta_idx[CASSETTE_END_MARKER] = len(meta_names)

    def get_skip_count(start: str, end: str) -> int:
        """Calculate how many metas were skipped."""
        if start not in meta_idx or end not in meta_idx:
            return -1  # Unknown
        start_i = meta_idx[start]
        end_i = meta_idx[end]
        skip = end_i - start_i - 1
        return max(0, skip)

    # Ensure we have the right columns
    if 'count' not in segment_types.columns:
        segment_types = pl.DataFrame(segment_types)

    # Calculate total for percentage conversion
    total_count = segment_types['count'].sum()

    # Add skip count and sort by count
    df = segment_types.select(['start_meta', 'end_meta', 'count'])

    # Calculate skip count for each row
    skip_counts = []
    for row in df.iter_rows(named=True):
        skip_counts.append(get_skip_count(row['start_meta'], row['end_meta']))

    df = df.with_columns(pl.Series('skip_count', skip_counts))
    df = df.sort('count', descending=True).head(top_n)

    if df.height == 0:
        if logger:
            logger.warning("No transition data to plot")
        return ax

    # Create labels and convert to percentages
    labels = [f"{r['start_meta']}→{r['end_meta']}" for r in df.iter_rows(named=True)]
    counts = df['count'].to_list()
    pcts = [(c / total_count * 100) if total_count > 0 else 0 for c in counts]
    skips = df['skip_count'].to_list()

    # Color mapping based on skip count
    colors = []
    for s in skips:
        if s == 0:
            colors.append('#2ecc71')  # Green - consecutive
        elif s == 1:
            colors.append('#f1c40f')  # Yellow - 1 skip
        elif s == 2:
            colors.append('#e67e22')  # Orange - 2 skip
        elif s >= 3:
            colors.append('#e74c3c')  # Red - 3+ skip
        else:
            colors.append('#95a5a6')  # Gray - unknown/edge

    # Create axes if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 8))

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, pcts, color=colors)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()  # Top to bottom
    ax.set_xlabel('% of transitions')
    ax.set_title(title)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2ecc71', label='Consecutive (0 skip)'),
        Patch(facecolor='#f1c40f', label='1 meta skipped'),
        Patch(facecolor='#e67e22', label='2 metas skipped'),
        Patch(facecolor='#e74c3c', label='3+ metas skipped'),
        Patch(facecolor='#95a5a6', label='Edge/Unknown'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    ax.grid(axis='x', alpha=0.3)

    return ax


def plot_meta_coverage_line(
    contigs_df: pl.DataFrame,
    metas: pl.DataFrame,
    ax=None,
    title: str = "Meta Coverage Across Cassette",
    contig_col: str = 'contig',
    logger: Optional[CustomLogger] = None,
):
    """
    Line plot showing presence/coverage of each meta across the cassette.
    X-axis: meta position (meta01, meta02, ..., meta15)
    Y-axis: % of contigs containing that meta

    Args:
        contigs_df: DataFrame with contig sequences
        metas: DataFrame with feature and seq columns
        ax: matplotlib axes (optional)
        title: Plot title
        contig_col: Name of column containing contig sequences
        logger: Optional logger

    Returns:
        matplotlib axes
    """
    import matplotlib.pyplot as plt

    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_names = meta_only['feature'].to_list()
    meta_seqs = dict(zip(meta_only['feature'], meta_only['seq']))

    total_contigs = contigs_df.height
    if total_contigs == 0:
        if logger:
            logger.warning("No contigs to analyze")
        return ax

    coverage = []
    for meta_name in meta_names:
        meta_seq = meta_seqs[meta_name]
        count = contigs_df.filter(pl.col(contig_col).str.contains(meta_seq)).height
        coverage.append(count / total_contigs * 100)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))

    x = range(len(meta_names))
    ax.plot(x, coverage, marker='o', linewidth=2, markersize=8, color='#3498db')
    ax.fill_between(x, coverage, alpha=0.3, color='#3498db')

    ax.set_xticks(x)
    ax.set_xticklabels(meta_names, rotation=45, ha='right')
    ax.set_ylabel('% Contigs with Meta')
    ax.set_xlabel('Meta Position')
    ax.set_title(title)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)

    # Add value labels on points
    for i, (xi, yi) in enumerate(zip(x, coverage)):
        ax.annotate(f'{yi:.0f}%', (xi, yi), textcoords="offset points",
                    xytext=(0, 10), ha='center', fontsize=8)

    return ax


def plot_meta_transition_sankey(
    segment_types: pl.DataFrame,
    metas: pl.DataFrame,
    min_count: int = 10,
    title: str = "Meta Transition Flow",
    logger: Optional[CustomLogger] = None,
):
    """
    Sankey diagram showing flow between meta sequences.

    Args:
        segment_types: DataFrame with start_meta, end_meta, count columns
        metas: DataFrame with feature column for meta ordering
        min_count: Minimum count to include in diagram
        title: Plot title
        logger: Optional logger

    Returns:
        plotly figure object (or None if plotly not available)
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        if logger:
            logger.warning("plotly not installed, skipping Sankey diagram")
        return None

    # Ensure we have the right columns
    if 'count' not in segment_types.columns:
        segment_types = pl.DataFrame(segment_types)

    # Filter to significant transitions
    sig_transitions = segment_types.filter(pl.col('count') >= min_count)

    if sig_transitions.height == 0:
        if logger:
            logger.warning(f"No transitions with count >= {min_count}")
        return None

    # Build node list (unique metas)
    all_metas = set(sig_transitions['start_meta'].to_list() + sig_transitions['end_meta'].to_list())

    # Order nodes by meta sequence order
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_order = [CASSETTE_START_MARKER] + meta_only['feature'].to_list() + [CASSETTE_END_MARKER]
    nodes = [m for m in meta_order if m in all_metas]

    # Add any remaining metas not in order
    for m in all_metas:
        if m not in nodes:
            nodes.append(m)

    node_idx = {n: i for i, n in enumerate(nodes)}

    # Build links
    sources = []
    targets = []
    values = []

    for row in sig_transitions.iter_rows(named=True):
        if row['start_meta'] in node_idx and row['end_meta'] in node_idx:
            sources.append(node_idx[row['start_meta']])
            targets.append(node_idx[row['end_meta']])
            values.append(row['count'])

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color="blue"
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
        )
    )])

    fig.update_layout(title_text=title, font_size=10)

    return fig


def plot_meta_transition_lines(
    segment_types: pl.DataFrame,
    metas: pl.DataFrame,
    fig=None,
    title: str = "Transition Prevalence",
    logger: Optional[CustomLogger] = None,
):
    """
    Plot horizontal lines connecting meta pairs at heights corresponding to prevalence.

    Creates two panels:
    - Top: Adjacent transitions only (no skips)
    - Bottom: Non-adjacent transitions (resections/skips)

    Each transition is shown as a horizontal line spanning from start_meta to end_meta,
    with the y-position (height) indicating the count/prevalence of that transition.

    Args:
        segment_types: DataFrame with start_meta, end_meta, count columns
        metas: DataFrame with feature column for meta ordering
        fig: matplotlib figure (optional, will create if None)
        title: Plot title
        logger: Optional logger

    Returns:
        matplotlib figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np

    # Get metas that actually appear in the segment_types data
    observed_metas = set(segment_types['start_meta'].unique().to_list() +
                        segment_types['end_meta'].unique().to_list())

    # Build ordered meta list - only include metas that appear in the data
    meta_only = metas.filter(pl.col('kind') == 'META')
    meta_names = meta_only['feature'].to_list()

    # Filter to only include metas that appear in segment_types
    meta_order = []
    if CASSETTE_START_MARKER in observed_metas:
        meta_order.append(CASSETTE_START_MARKER)
    for m in meta_names:
        if m in observed_metas:
            meta_order.append(m)
    if CASSETTE_END_MARKER in observed_metas:
        meta_order.append(CASSETTE_END_MARKER)

    meta_to_idx = {m: i for i, m in enumerate(meta_order)}

    # Filter to valid transitions and add position info
    transitions = []
    total_count = 0
    for row in segment_types.iter_rows(named=True):
        start = row['start_meta']
        end = row['end_meta']
        count = row['count']

        if start not in meta_to_idx or end not in meta_to_idx:
            continue

        start_idx = meta_to_idx[start]
        end_idx = meta_to_idx[end]

        # Calculate skip distance (how many metas are skipped)
        skip = abs(end_idx - start_idx) - 1

        transitions.append({
            'start': start,
            'end': end,
            'start_idx': start_idx,
            'end_idx': end_idx,
            'count': count,
            'skip': skip,
        })
        total_count += count

    # Convert counts to percentages
    for t in transitions:
        t['pct'] = (t['count'] / total_count * 100) if total_count > 0 else 0

    # Split into adjacent and non-adjacent
    adjacent = [t for t in transitions if t['skip'] == 0]
    non_adjacent = [t for t in transitions if t['skip'] > 0]

    # Create figure with two panels
    if fig is None:
        fig, (ax_adj, ax_skip) = plt.subplots(2, 1, figsize=(14, 10), height_ratios=[1, 1.2])
    else:
        ax_adj, ax_skip = fig.subplots(2, 1, height_ratios=[1, 1.2])

    fig.suptitle(title, fontsize=12)

    # Color scheme for non-adjacent based on skip distance
    def get_skip_color(skip):
        if skip == 1:
            return '#f39c12'  # orange - 1 skip
        elif skip == 2:
            return '#e74c3c'  # red - 2 skips
        else:
            return '#8e44ad'  # purple - 3+ skips

    def draw_transitions(ax, trans_list, color_func, panel_title):
        if not trans_list:
            ax.text(0.5, 0.5, "No transitions", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(panel_title)
            return

        # Sort by percentage (draw lower values first)
        trans_list = sorted(trans_list, key=lambda x: x['pct'])

        max_pct = max(t['pct'] for t in trans_list)

        for t in trans_list:
            x_start = t['start_idx']
            x_end = t['end_idx']
            y = t['pct']
            color = color_func(t['skip'])

            # Line width proportional to percentage (log scale)
            if max_pct > 0.01:
                lw = 0.8 + 2.0 * (np.log10(t['pct'] + 0.01) / np.log10(max_pct + 0.01))
            else:
                lw = 1.5

            # Draw horizontal line
            ax.hlines(y=y, xmin=x_start, xmax=x_end, colors=color, linewidth=lw, alpha=0.7)

            # Add dots at endpoints
            ax.scatter([x_start, x_end], [y, y], color=color, s=lw*8, alpha=0.7, zorder=5)

        # Set axis properties
        ax.set_yscale('log')
        ax.set_xlim(-0.5, len(meta_order) - 0.5)
        ax.set_xticks(range(len(meta_order)))
        ax.set_xticklabels(meta_order, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('% of transitions')
        ax.set_title(panel_title)
        ax.grid(True, alpha=0.3, axis='y')

    # Panel 1: Adjacent transitions (green)
    draw_transitions(ax_adj, adjacent, lambda s: '#2ecc71', "Adjacent Transitions (no resection)")

    # Panel 2: Non-adjacent transitions (colored by skip distance)
    draw_transitions(ax_skip, non_adjacent, get_skip_color, "Resections (skipped metas)")
    ax_skip.set_xlabel('Meta Position')

    # Legend for bottom panel
    if non_adjacent:
        legend_elements = [
            mpatches.Patch(color='#f39c12', label='1 meta skipped'),
            mpatches.Patch(color='#e74c3c', label='2 metas skipped'),
            mpatches.Patch(color='#8e44ad', label='3+ metas skipped'),
        ]
        ax_skip.legend(handles=legend_elements, loc='upper right', fontsize=8)

    return fig


def plot_umi_transition_heterogeneity(
    report: dict,
    ax=None,
    title: str = "Transition Types per UMI",
    logger: Optional[CustomLogger] = None,
):
    """
    Plot distribution of transition types per UMI (Smart vs Sensitive mode).

    UMIs with multiple transition types may indicate:
    - UMI collisions (different molecules with same UMI)
    - Chimeric reads
    - Assembly errors

    Shows both Smart mode (filtered by adaptive threshold) and Sensitive mode
    (all transitions) side-by-side for comparison.

    Args:
        report: Report dict from generate_segmentation_report()
        ax: matplotlib axes (optional)
        title: Plot title
        logger: Optional logger

    Returns:
        matplotlib axes
    """
    import matplotlib.pyplot as plt
    import numpy as np

    seg = report['segments']
    dist_smart = seg.get('umi_transition_distribution', [])
    dist_sensitive = seg.get('umi_transition_distribution_sensitive', [])

    if not dist_smart:
        if logger:
            logger.warning("No UMI transition distribution data")
        return ax

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    # Get all unique n_types from both distributions
    all_n_types = set()
    for d in dist_smart:
        all_n_types.add(d['n_transition_types'])
    for d in dist_sensitive:
        all_n_types.add(d['n_transition_types'])
    all_n_types = sorted(all_n_types)

    # Build count dicts for easy lookup
    smart_counts = {d['n_transition_types']: d['len'] for d in dist_smart}
    sensitive_counts = {d['n_transition_types']: d['len'] for d in dist_sensitive}

    # Extract aligned data
    counts_smart = [smart_counts.get(n, 0) for n in all_n_types]
    counts_sensitive = [sensitive_counts.get(n, 0) for n in all_n_types]
    total_smart = sum(counts_smart) if counts_smart else 1
    total_sensitive = sum(counts_sensitive) if counts_sensitive else 1

    # Bar positions
    x = np.arange(len(all_n_types))
    width = 0.35

    # Create grouped bar plot
    bars_smart = ax.bar(x - width/2, counts_smart, width, label='Smart (filtered)',
                        color='#3498db', edgecolor='black', linewidth=0.5, alpha=0.8)
    bars_sensitive = ax.bar(x + width/2, counts_sensitive, width, label='Sensitive (all)',
                            color='#e74c3c', edgecolor='black', linewidth=0.5, alpha=0.6)

    # Add percentage labels on bars (only for significant bars)
    for bar, count in zip(bars_smart, counts_smart):
        if count > 0:
            pct = count / total_smart * 100
            height = bar.get_height()
            ax.annotate(f'{pct:.0f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, color='#2980b9')

    for bar, count in zip(bars_sensitive, counts_sensitive):
        if count > 0:
            pct = count / total_sensitive * 100
            height = bar.get_height()
            ax.annotate(f'{pct:.0f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, color='#c0392b')

    ax.set_xlabel('Number of different transition types per UMI')
    ax.set_ylabel('Number of UMIs')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(all_n_types)
    ax.legend(loc='upper right')

    # Add summary stats as text
    single_smart = seg.get('umis_single_transition_type', 0)
    multi_smart = seg.get('umis_multi_transition_type', 0)
    single_sens = seg.get('umis_single_transition_type_sensitive', 0)
    multi_sens = seg.get('umis_multi_transition_type_sensitive', 0)
    threshold = seg.get('heterogeneity_threshold', 0.20)

    stats_text = f"Threshold: {threshold:.0%} of dominant\n"
    if single_smart + multi_smart > 0:
        multi_pct_smart = multi_smart / (single_smart + multi_smart) * 100
        stats_text += f"Smart: {multi_smart:,} multi ({multi_pct_smart:.1f}%)\n"
    if single_sens + multi_sens > 0:
        multi_pct_sens = multi_sens / (single_sens + multi_sens) * 100
        stats_text += f"Sensitive: {multi_sens:,} multi ({multi_pct_sens:.1f}%)"

    ax.text(0.98, 0.50, stats_text,
            transform=ax.transAxes, ha='right', va='center',
            fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Use log scale if range is large
    all_counts = counts_smart + counts_sensitive
    nonzero_counts = [c for c in all_counts if c > 0]
    if nonzero_counts and max(nonzero_counts) / min(nonzero_counts) > 100:
        ax.set_yscale('log')
        ax.set_ylabel('Number of UMIs (log scale)')

    ax.grid(True, alpha=0.3, axis='y')

    return ax


def plot_segmentation_qc(
    report: dict,
    contigs_df: pl.DataFrame,
    metas: pl.DataFrame,
    output_dir,
    sample_name: str,
    contig_col: str = 'contig',
    logger: Optional[CustomLogger] = None,
):
    """
    Generate all segmentation QC plots and save to output directory.

    Generates:
    - {sample}_meta_transition_heatmap.png
    - {sample}_meta_transition_bars.png
    - {sample}_meta_coverage_line.png
    - {sample}_meta_transition_lines.png (horizontal lines at prevalence heights)
    - {sample}_meta_transition_sankey.html (interactive, if plotly available)
    - {sample}_umi_heterogeneity.png (UMI collision/chimera indicator)

    Args:
        report: Report dict from generate_segmentation_report()
        contigs_df: DataFrame with contig sequences
        metas: DataFrame with feature and seq columns
        output_dir: Path to output directory
        sample_name: Sample name for file naming
        contig_col: Name of column containing contig sequences
        logger: Optional logger
    """
    import matplotlib.pyplot as plt
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get segment_types from report
    segment_types = pl.DataFrame(report['segments']['segment_types'])

    if logger:
        logger.info(f"Generating segmentation QC plots for {sample_name}")

    # 1. Heatmap
    try:
        fig, ax = plt.subplots(figsize=(12, 10))
        plot_meta_transition_heatmap(segment_types, metas, ax=ax,
                                      title=f"{sample_name} - Meta Transitions", logger=logger)
        fig.tight_layout()
        heatmap_path = output_dir / f"{sample_name}_meta_transition_heatmap.png"
        fig.savefig(heatmap_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        if logger:
            logger.info(f"Saved {heatmap_path}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to create heatmap: {e}")

    # 2. Bar chart
    try:
        fig, ax = plt.subplots(figsize=(14, 10))
        plot_meta_transition_bars(segment_types, metas, ax=ax,
                                   title=f"{sample_name} - Transition Frequencies", logger=logger)
        fig.tight_layout()
        bars_path = output_dir / f"{sample_name}_meta_transition_bars.png"
        fig.savefig(bars_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        if logger:
            logger.info(f"Saved {bars_path}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to create bar chart: {e}")

    # 3. Coverage line
    try:
        fig, ax = plt.subplots(figsize=(12, 5))
        plot_meta_coverage_line(contigs_df, metas, ax=ax,
                                title=f"{sample_name} - Meta Coverage", contig_col=contig_col, logger=logger)
        fig.tight_layout()
        line_path = output_dir / f"{sample_name}_meta_coverage_line.png"
        fig.savefig(line_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        if logger:
            logger.info(f"Saved {line_path}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to create coverage line: {e}")

    # 4. Transition lines (horizontal lines at prevalence heights, two panels)
    try:
        fig = plot_meta_transition_lines(segment_types, metas,
                                         title=f"{sample_name} - Transition Prevalence", logger=logger)
        fig.tight_layout()
        lines_path = output_dir / f"{sample_name}_meta_transition_lines.png"
        fig.savefig(lines_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        if logger:
            logger.info(f"Saved {lines_path}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to create transition lines: {e}")

    # 5. Sankey (HTML/interactive)
    try:
        sankey_fig = plot_meta_transition_sankey(segment_types, metas,
                                                  title=f"{sample_name} - Transition Flow", logger=logger)
        if sankey_fig is not None:
            sankey_path = output_dir / f"{sample_name}_meta_transition_sankey.html"
            sankey_fig.write_html(str(sankey_path))
            if logger:
                logger.info(f"Saved {sankey_path}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to create Sankey diagram: {e}")

    # 6. UMI transition heterogeneity (collision/chimera detection)
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        plot_umi_transition_heterogeneity(report, ax=ax,
                                          title=f"{sample_name} - UMI Transition Heterogeneity", logger=logger)
        fig.tight_layout()
        hetero_path = output_dir / f"{sample_name}_umi_heterogeneity.png"
        fig.savefig(hetero_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        if logger:
            logger.info(f"Saved {hetero_path}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to create UMI heterogeneity plot: {e}")
