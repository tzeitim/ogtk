"""
Segmentation utilities for splitting long cassette reads at META boundaries.

This module provides functions to:
1. Sanitize METAs (normalize fuzzy matches to canonical sequences)
2. Segment reads at META boundaries into overlapping segments
3. Stitch assembled segments back into full sequences

All functions operate on LazyFrames for optimal performance with large datasets.
"""

import polars as pl
from typing import Optional
from ogtk.utils.general import fuzzy_match_str
from ogtk.utils.log import CustomLogger


# Default META length (all PEtracer METAs are 20bp)
META_LEN = 20


def sanitize_metas(
    ldf: pl.LazyFrame,
    metas: pl.DataFrame,
    seq_col: str = 'r2_seq',
    wildcard: str = ".{0,1}",
    logger: Optional[CustomLogger] = None
) -> pl.LazyFrame:
    """
    Replace fuzzy META matches with canonical sequences to handle sequencing errors.

    Args:
        ldf: LazyFrame with sequence data
        metas: DataFrame with 'feature' (META names) and 'seq' (META sequences)
        seq_col: Name of the column containing sequences
        wildcard: Regex pattern for fuzzy matching (default allows 1 error)
        logger: Optional logger for debugging

    Returns:
        LazyFrame with sanitized sequences
    """
    # Filter to only META sequences (not TARGETs)
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas

    for seq, name in zip(meta_only['seq'], meta_only['feature']):
        fuzzy_pattern = fuzzy_match_str(seq, wildcard=wildcard, include_original=True)
        if logger:
            logger.debug(f"Sanitizing {name}")
        ldf = ldf.with_columns(pl.col(seq_col).str.replace_all(fuzzy_pattern, seq))

    return ldf


# Special markers for cassette edge segments
CASSETTE_START_MARKER = "_CASSETTE_START_"
CASSETTE_END_MARKER = "_CASSETTE_END_"


def segment_by_metas(
    ldf: pl.LazyFrame,
    metas: pl.DataFrame,
    seq_col: str = 'r2_seq',
    keep_cols: Optional[list] = None,
    cassette_start_anchor: Optional[str] = None,
    cassette_end_anchor: Optional[str] = None,
    logger: Optional[CustomLogger] = None
) -> pl.LazyFrame:
    """
    Split reads at META boundaries into overlapping segments.

    Each segment includes both boundary METAs for overlap-based stitching.
    Optionally extracts edge segments (cassette start/end).

    Note: TARGET insertion extraction should be done AFTER assembly (on consensus_seq)
    to benefit from noise reduction. See api_ext.py Step 4b.

    Args:
        ldf: LazyFrame with sequence data (should be sanitized first)
        metas: DataFrame with 'feature' and 'seq' columns
        seq_col: Name of the column containing sequences
        keep_cols: Additional columns to preserve (default: ['umi'])
        cassette_start_anchor: Sequence marking cassette start (for edge segment extraction)
        cassette_end_anchor: Sequence marking cassette end (for edge segment extraction)
        logger: Optional logger for debugging

    Returns:
        LazyFrame with columns: [keep_cols], start_meta, end_meta, segment_seq
    """
    if keep_cols is None:
        keep_cols = ['umi']

    # Keep reference to original ldf for edge segment extraction
    original_ldf = ldf

    # Filter to only META sequences
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas

    meta_seqs = meta_only['seq'].to_list()
    meta_names = meta_only['feature'].to_list()
    seq_to_name = dict(zip(meta_seqs, meta_names))

    if logger:
        logger.info(f"Segmenting by {len(meta_seqs)} META sequences")

    # Pre-compute positions for all METAs (str.find only accepts string literals, not expressions)
    pos_exprs = [
        pl.col(seq_col).str.find(seq).alias(f'_pos_{name}')
        for name, seq in zip(meta_names, meta_seqs)
    ]

    ldf = (
        ldf
        # Extract all METAs in order of appearance (returns list of matched sequences)
        # NOTE: find_many returns positions (u32), extract_many returns the actual strings
        .with_columns(pl.col(seq_col).str.extract_many(meta_seqs).alias('meta_matches'))
        # Add position columns for each META
        .with_columns(pos_exprs)
    )

    # Create mapping DataFrame for seq -> name conversion
    seq_name_map = pl.DataFrame({
        'meta_seq': meta_seqs,
        'meta_name': meta_names,
    }).lazy()

    # Add row index to preserve order after explode/join/group
    ldf = ldf.with_row_index('_row_idx')

    # Explode meta_matches, join to get names, then create consecutive pairs
    # This avoids map_elements and uses native Polars operations
    ldf = (
        ldf
        # Explode to one row per META match
        .explode('meta_matches')
        # Join to convert sequence to name
        .join(seq_name_map, left_on='meta_matches', right_on='meta_seq', how='left')
        # Add index within each row to track position
        .with_columns(
            pl.col('meta_name').cum_count().over('_row_idx').alias('_meta_idx')
        )
    )

    # Create consecutive pairs by self-joining on adjacent indices
    ldf_next = (
        ldf
        .select(
            pl.col('_row_idx'),
            (pl.col('_meta_idx') - 1).alias('_meta_idx'),  # Shift to join with previous
            pl.col('meta_name').alias('end_meta'),
        )
    )

    ldf = (
        ldf
        .join(
            ldf_next,
            on=['_row_idx', '_meta_idx'],
            how='left'
        )
        .filter(pl.col('end_meta').is_not_null())
        .rename({'meta_name': 'start_meta'})
    )

    # Build position selection expressions using when/then chains
    # (since str.find doesn't accept column expressions)
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

    # Get positions and extract segments
    pos_cols = [f'_pos_{name}' for name in meta_names]

    ldf = (
        ldf
        .with_columns(
            start_pos_expr.alias('start_pos'),
            end_pos_expr.alias('end_pos'),
        )
        # Extract segment (includes both boundary METAs)
        .with_columns(
            pl.col(seq_col).str.slice(
                pl.col('start_pos'),
                pl.col('end_pos') - pl.col('start_pos') + META_LEN
            ).alias('segment_seq')
        )
    )

    meta_segments = ldf.select(keep_cols + ['start_meta', 'end_meta', 'segment_seq'])

    # Extract edge segments if cassette anchors are provided
    if not cassette_start_anchor and not cassette_end_anchor:
        return meta_segments

    # Use original_ldf for edge extraction (before META processing modified it)
    edge_segments = extract_edge_segments(
        ldf=original_ldf,
        metas=metas,
        cassette_start_anchor=cassette_start_anchor,
        cassette_end_anchor=cassette_end_anchor,
        seq_col=seq_col,
        keep_cols=keep_cols,
        logger=logger,
    )

    return pl.concat([meta_segments, edge_segments])


def extract_edge_segments(
    ldf: pl.LazyFrame,
    metas: pl.DataFrame,
    cassette_start_anchor: str | None = None,
    cassette_end_anchor: str | None = None,
    seq_col: str = 'r2_seq',
    keep_cols: Optional[list] = None,
    logger: Optional[CustomLogger] = None
) -> pl.LazyFrame:
    """
    Extract edge segments from cassette start/end anchors to first/last METAs.

    This complements segment_by_metas by capturing the regions:
    - [cassette_start_anchor]...[first_META]
    - [last_META]...[cassette_end_anchor]

    Note: TARGET insertion extraction should be done AFTER assembly (on consensus_seq)
    to benefit from noise reduction. See api_ext.py Step 4b.

    Args:
        ldf: LazyFrame with sequence data (should be sanitized first)
        metas: DataFrame with 'feature' and 'seq' columns
        cassette_start_anchor: Sequence marking the start of the cassette
        cassette_end_anchor: Sequence marking the end of the cassette
        seq_col: Name of the column containing sequences
        keep_cols: Additional columns to preserve (default: ['umi'])
        logger: Optional logger for debugging

    Returns:
        LazyFrame with columns: [keep_cols], start_meta, end_meta, segment_seq
        Uses special markers CASSETTE_START_MARKER and CASSETTE_END_MARKER
    """
    if keep_cols is None:
        keep_cols = ['umi']

    if not cassette_start_anchor and not cassette_end_anchor:
        # Return empty LazyFrame with correct schema
        return pl.LazyFrame(schema={
            **{col: pl.String for col in keep_cols},
            'start_meta': pl.String,
            'end_meta': pl.String,
            'segment_seq': pl.String,
        })

    # Filter to only META sequences
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    meta_seqs = meta_only['seq'].to_list()
    meta_names = meta_only['feature'].to_list()

    first_meta_name = meta_names[0]
    last_meta_name = meta_names[-1]
    first_meta_seq = meta_seqs[0]
    last_meta_seq = meta_seqs[-1]

    if logger:
        if cassette_start_anchor:
            logger.debug(f"Extracting start edge: {cassette_start_anchor[:20]}... -> {first_meta_name}")
        if cassette_end_anchor:
            logger.debug(f"Extracting end edge: {last_meta_name} -> ...{cassette_end_anchor[-20:]}")

    edge_segments = []

    # Extract start edge segment: cassette_start_anchor -> first_META
    if cassette_start_anchor:
        start_edge = (
            ldf
            .with_columns([
                pl.col(seq_col).str.find(cassette_start_anchor).alias('_start_anchor_pos'),
                pl.col(seq_col).str.find(first_meta_seq).alias('_first_meta_pos'),
            ])
            # Only include if start anchor is found AND is before first META
            .filter(
                pl.col('_start_anchor_pos').is_not_null() &
                pl.col('_first_meta_pos').is_not_null() &
                (pl.col('_start_anchor_pos') < pl.col('_first_meta_pos'))
            )
            .with_columns([
                pl.lit(CASSETTE_START_MARKER).alias('start_meta'),
                pl.lit(first_meta_name).alias('end_meta'),
                pl.col(seq_col).str.slice(
                    pl.col('_start_anchor_pos'),
                    pl.col('_first_meta_pos') - pl.col('_start_anchor_pos') + META_LEN
                ).alias('segment_seq'),
            ])
            .select(keep_cols + ['start_meta', 'end_meta', 'segment_seq'])
        )
        edge_segments.append(start_edge)

    # Extract end edge segment: last_META -> cassette_end_anchor
    if cassette_end_anchor:
        end_anchor_len = len(cassette_end_anchor)
        end_edge = (
            ldf
            .with_columns([
                pl.col(seq_col).str.find(last_meta_seq).alias('_last_meta_pos'),
                pl.col(seq_col).str.find(cassette_end_anchor).alias('_end_anchor_pos'),
            ])
            # Only include if last META is found AND is before end anchor
            .filter(
                pl.col('_last_meta_pos').is_not_null() &
                pl.col('_end_anchor_pos').is_not_null() &
                (pl.col('_last_meta_pos') < pl.col('_end_anchor_pos'))
            )
            .with_columns([
                pl.lit(last_meta_name).alias('start_meta'),
                pl.lit(CASSETTE_END_MARKER).alias('end_meta'),
                pl.col(seq_col).str.slice(
                    pl.col('_last_meta_pos'),
                    pl.col('_end_anchor_pos') - pl.col('_last_meta_pos') + end_anchor_len
                ).alias('segment_seq'),
            ])
            .select(keep_cols + ['start_meta', 'end_meta', 'segment_seq'])
        )
        edge_segments.append(end_edge)

    if len(edge_segments) == 1:
        return edge_segments[0]
    else:
        return pl.concat(edge_segments)


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

    # Build META order map
    # CASSETTE_START_MARKER gets order -1 (before all METAs)
    meta_names = meta_only['feature'].to_list()
    meta_order = {CASSETTE_START_MARKER: -1}
    meta_order.update({m: i for i, m in enumerate(meta_names)})

    if logger:
        logger.info(f"Stitching segments using {len(meta_order)} META order positions (including edge markers)")
        logger.info(f"Grouping by: {group_cols}")

    # Build sort order expression using when/then chain
    sort_order_expr = pl.lit(None).cast(pl.Int64)
    for name, order in meta_order.items():
        sort_order_expr = (
            pl.when(pl.col('start_meta') == name)
            .then(pl.lit(order))
            .otherwise(sort_order_expr)
        )

    # Build trimming expression:
    # - First segment (min sort_order): keep full
    # - Start edge followed by META segment: trim META_LEN from META segment (normal case)
    # - Other segments: trim META_LEN from start (duplicate META from previous)
    # Use group_cols for the over() clause to handle sbc correctly
    trim_expr = (
        pl.when(pl.col('sort_order') == pl.col('sort_order').min().over(group_cols))
        .then(pl.col(seq_col))  # First segment: keep full
        .otherwise(pl.col(seq_col).str.slice(META_LEN))  # Others: trim start META
    )

    return (
        ldf
        # Add sort order based on META position
        .with_columns(sort_order_expr.alias('sort_order'))
        .sort(*group_cols, 'sort_order')
        # Apply trimming logic
        .with_columns(trim_expr.alias('contribution'))
        # Concatenate within each group (umi or sbc+umi)
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
    # Filter to only META sequences
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas

    # Add count column for each META
    count_exprs = [
        pl.col(seq_col).str.count_matches(pl.lit(seq)).alias(f'{name}_count')
        for name, seq in zip(meta_only['feature'], meta_only['seq'])
    ]

    ldf = ldf.with_columns(count_exprs)

    # Flag if any META appears more than once
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
    """
    Generate a comprehensive diagnostic report for segmentation-based assembly.

    Args:
        segments_df: DataFrame with extracted segments (umi, start_meta, end_meta, segment_seq)
        assembled_df: DataFrame with assembled segments (umi, start_meta, end_meta, consensus_seq)
        contigs_df: DataFrame with stitched contigs (umi, contig)
        metas: DataFrame with META sequences
        cassette_start_anchor: Cassette start anchor (if used)
        cassette_end_anchor: Cassette end anchor (if used)

    Returns:
        Dictionary with diagnostic metrics
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

    # ==================== SEGMENT STATS ====================
    total_segments = segments_df.height
    unique_umis_with_segments = segments_df.select('umi').n_unique()

    # Segments by type
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

    # Expected consecutive pairs
    expected_pairs = [(meta_names[i], meta_names[i+1]) for i in range(len(meta_names)-1)]
    consecutive_segment_counts = {}
    for start, end in expected_pairs:
        count = segments_df.filter(
            (pl.col('start_meta') == start) & (pl.col('end_meta') == end)
        ).height
        consecutive_segment_counts[f'{start}->{end}'] = count

    # Edge segment counts
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

    # ==================== ASSEMBLY STATS ====================
    total_assembled = assembled_df.height
    unique_umis_assembled = assembled_df.select('umi').n_unique()

    # Assembly success by segment type
    assembly_by_type = (
        assembled_df
        .group_by(['start_meta', 'end_meta'])
        .agg([
            pl.len().alias('assembled_count'),
            pl.col('consensus_seq').str.len_chars().mean().alias('mean_consensus_len'),
        ])
        .sort('assembled_count', descending=True)
    )

    # Calculate assembly success rate per type
    assembly_success = {}
    for row in segment_types.iter_rows(named=True):
        key = f"{row['start_meta']}->{row['end_meta']}"
        input_count = row['count']
        # Count unique (umi, start_meta, end_meta) groups in input
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

    # ==================== STITCHING STATS ====================
    total_contigs = contigs_df.height

    # Contig length distribution
    contig_lengths = contigs_df.select(
        pl.col('contig').str.len_chars().alias('length')
    )

    length_stats = {
        'count': total_contigs,
        'mean': round(contig_lengths['length'].mean() or 0, 1),
        'median': contig_lengths['length'].median() or 0,
        'min': contig_lengths['length'].min() or 0,
        'max': contig_lengths['length'].max() or 0,
        'std': round(contig_lengths['length'].std() or 0, 1),
    }

    # How many segments per UMI in assembled data
    segments_per_umi = (
        assembled_df
        .group_by('umi')
        .len()
        .select('len')
    )

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

    # ==================== META PRESENCE ====================
    meta_presence = {}
    for meta_name in meta_names:
        meta_seq = meta_seqs[meta_name]
        count = contigs_df.filter(pl.col('contig').str.contains(meta_seq)).height
        pct = count / total_contigs * 100 if total_contigs > 0 else 0
        meta_presence[meta_name] = {
            'count': count,
            'percent': round(pct, 1),
        }

    # Check cassette anchor presence if provided
    if cassette_start_anchor:
        count = contigs_df.filter(pl.col('contig').str.contains(cassette_start_anchor)).height
        pct = count / total_contigs * 100 if total_contigs > 0 else 0
        meta_presence['CASSETTE_START'] = {'count': count, 'percent': round(pct, 1)}

    if cassette_end_anchor:
        count = contigs_df.filter(pl.col('contig').str.contains(cassette_end_anchor)).height
        pct = count / total_contigs * 100 if total_contigs > 0 else 0
        meta_presence['CASSETTE_END'] = {'count': count, 'percent': round(pct, 1)}

    report['meta_presence'] = meta_presence

    # ==================== SUMMARY ====================
    # Calculate overall metrics
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
    """
    Format the segmentation report as a human-readable string.

    Args:
        report: Dictionary from generate_segmentation_report()

    Returns:
        Formatted string report
    """
    lines = []
    lines.append("=" * 60)
    lines.append("SEGMENTATION ASSEMBLY REPORT")
    lines.append("=" * 60)

    # Summary
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

    # Segment stats
    seg = report['segments']
    lines.append("\n## Segment Extraction")
    lines.append(f"  Total segments:        {seg['total_segments']:,}")
    lines.append(f"  Unique UMIs:           {seg['unique_umis']:,}")
    if seg['edge_start_segments'] > 0:
        lines.append(f"  Start edge segments:   {seg['edge_start_segments']:,}")
    if seg['edge_end_segments'] > 0:
        lines.append(f"  End edge segments:     {seg['edge_end_segments']:,}")

    # Top segment types
    lines.append("\n  Top segment types:")
    for item in seg['segment_types'][:5]:
        lines.append(f"    {item['start_meta']:>18} -> {item['end_meta']:<10}: {item['count']:>8,} (median {item['median_length']:.0f}bp)")

    # Assembly success rates
    asm = report['assembly']
    lines.append("\n## Assembly")
    lines.append(f"  Total assembled:       {asm['total_assembled']:,}")
    lines.append(f"  Unique UMIs:           {asm['unique_umis']:,}")

    lines.append("\n  Assembly success by type:")
    for seg_type, stats in sorted(asm['by_type'].items(), key=lambda x: -x[1]['assembled']):
        if stats['assembled'] > 0:
            lines.append(f"    {seg_type:>22}: {stats['assembled']:>6,} / {stats['input_umi_groups']:>6,} ({stats['success_rate_pct']:>5.1f}%)")

    # META presence
    meta = report['meta_presence']
    lines.append("\n## META Presence in Contigs")
    for meta_name, stats in meta.items():
        bar_len = int(stats['percent'] / 5)  # 20 chars = 100%
        bar = "█" * bar_len + "░" * (20 - bar_len)
        lines.append(f"  {meta_name:>15}: {bar} {stats['percent']:>5.1f}% ({stats['count']:,})")

    lines.append("\n" + "=" * 60)

    return "\n".join(lines)
