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
) -> dict:
    """Generate a comprehensive diagnostic report for segmentation-based assembly."""
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

    segment_types = (
        segments_df
        .group_by(['start_meta', 'end_meta'])
        .agg([
            pl.len().alias('count'),
            pl.col('segment_seq').str.len_chars().mean().alias('mean_length'),
            pl.col('segment_seq').str.len_chars().median().alias('median_length'),
        ])
        .sort('count', descending=True)
    )

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

    lines.append("\n  Top segment types:")
    for item in seg['segment_types'][:5]:
        lines.append(f"    {item['start_meta']:>18} -> {item['end_meta']:<10}: {item['count']:>8,} (median {item['median_length']:.0f}bp)")

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
