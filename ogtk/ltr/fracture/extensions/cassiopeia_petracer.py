import polars as pl
from pathlib import Path
from typing import Set,Dict, Optional, Type, List, Tuple
from enum import Enum
from .base import PostProcessorExtension
from .registry import extension_registry
from .config import ExtensionConfig
from ..pipeline.types import StepResults,FractureXp
from ..pipeline.formats import scan_file, read_file
from ..pipeline import api_ext  # noqa: F401  registers the `dna` namespace
from dataclasses import dataclass, field
from ogtk.utils.log import CustomLogger
from ogtk.utils.general import fuzzy_match_str


def collapse_umis_to_cells(
    allele_df: pl.DataFrame,
    cell_col: str = 'cellBC',
    intbc_col: str = 'intBC',
    min_umis_per_cell: int = 4,
    min_umi_agreement: Optional[float] = 0.5,
    method: str = 'mode_site',
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Collapse multiple UMIs per cell into one consensus allele per (cell, intBC).

    Dispatches to one of two strategies:
      - 'mode_site': per-cutsite mode voting across UMIs (independent per r column)
      - 'mode_umi': pick the allele pattern from the most frequent whole-UMI,
                    using readCount as tiebreaker

    Args:
        allele_df: DataFrame with cellBC, intBC, r001..r015, readCount columns
        cell_col: Column name for cell barcode
        intbc_col: Column name for integration barcode
        min_umis_per_cell: Minimum UMIs required per (cell, intBC) group.
                          Groups with fewer UMIs are filtered out.
        min_umi_agreement: (mode_site only) Minimum fraction for consensus.
                          If the mode has support below this threshold, mark as None.
        method: Collapse strategy — 'mode_site' or 'mode_umi'
        logger: Optional logger for info messages

    Returns:
        DataFrame with one row per (cell, intBC), adds n_umis column.
    """
    if allele_df.height == 0:
        return allele_df.with_columns(pl.lit(0).alias('n_umis'))

    # Identify r columns (allele columns)
    rcols = allele_df.select(pl.col('^r\\d+$')).columns
    if not rcols:
        if logger:
            logger.warning("No r columns found in allele_df, returning as-is with n_umis=1")
        return allele_df.with_columns(pl.lit(1).alias('n_umis'))

    # Step 1: Get group sizes and filter by min_umis_per_cell
    group_sizes = (
        allele_df.group_by([cell_col, intbc_col])
        .agg(pl.len().alias('n_umis'))
    )

    n_groups_before = group_sizes.height
    valid_groups = group_sizes.filter(pl.col('n_umis') >= min_umis_per_cell)
    n_groups_after = valid_groups.height

    if logger and n_groups_before != n_groups_after:
        logger.info(f"UMI collapse: filtered {n_groups_before - n_groups_after} (cell, intBC) groups "
                   f"with < {min_umis_per_cell} UMIs ({n_groups_after} remaining)")

    if valid_groups.height == 0:
        preserve_cols = [c for c in allele_df.columns if c not in rcols + [cell_col, intbc_col]]
        empty_cols = [cell_col, intbc_col, 'n_umis'] + rcols + preserve_cols
        return pl.DataFrame(schema={c: allele_df.schema.get(c, pl.Utf8) for c in empty_cols if c in allele_df.columns or c == 'n_umis'})

    # Filter to valid groups
    df = allele_df.join(valid_groups.select([cell_col, intbc_col]), on=[cell_col, intbc_col], how='semi')

    # Step 2: Dispatch to strategy
    if method == 'mode_site':
        result = _collapse_mode_site(df, rcols, cell_col, intbc_col, min_umi_agreement)
    elif method == 'mode_umi':
        result = _collapse_mode_umi(df, rcols, cell_col, intbc_col)
    else:
        raise ValueError(f"Unknown collapse method: {method!r}. Use 'mode_site' or 'mode_umi'.")

    # Add n_umis column
    result = result.join(valid_groups, on=[cell_col, intbc_col], how='left')

    # Reorder columns: group keys + n_umis + rcols + anything else
    output_cols = [cell_col, intbc_col, 'n_umis'] + rcols
    output_cols += [c for c in result.columns if c not in output_cols]
    result = result.select(output_cols)

    if logger:
        logger.info(f"UMI collapse ({method}): {allele_df.height} rows -> {result.height} (cell, intBC) pairs")

    return result


def _collapse_mode_site(
    df: pl.DataFrame,
    rcols: list[str],
    cell_col: str,
    intbc_col: str,
    min_umi_agreement: Optional[float],
) -> pl.DataFrame:
    """Per-cutsite mode voting.

    For each (cell, intBC, cutsite), independently picks the most frequent
    allele value across UMIs. Positions where the mode has support below
    min_umi_agreement are set to None (missing data).
    """
    # Columns to preserve (take first value per group)
    preserve_cols = [c for c in df.columns if c not in rcols + [cell_col, intbc_col]]

    df_idx = df.with_row_index('_row_idx')

    # Unpivot r columns to long format
    id_cols = [cell_col, intbc_col, '_row_idx'] + preserve_cols
    unpivoted = df_idx.unpivot(
        on=rcols,
        index=id_cols,
        variable_name='_rcol',
        value_name='_allele'
    )

    # For each (cell, intBC, rcol), compute mode and support
    consensus = (
        unpivoted
        .group_by([cell_col, intbc_col, '_rcol'])
        .agg([
            pl.len().alias('_n_total'),
            pl.col('_allele').value_counts(sort=True).first().alias('_mode_struct'),
        ])
        .with_columns([
            pl.col('_mode_struct').struct.field('_allele').alias('_mode_value'),
            pl.col('_mode_struct').struct.field('count').alias('_mode_count'),
        ])
        .with_columns(
            (pl.col('_mode_count') / pl.col('_n_total')).alias('_support')
        )
        .with_columns(
            # When min_umi_agreement is None, keep all mode values (no threshold)
            (pl.col('_mode_value') if min_umi_agreement is None else
             pl.when(pl.col('_support') >= min_umi_agreement)
             .then(pl.col('_mode_value'))
             .otherwise(pl.lit(None)))
            .alias('_consensus_value')
        )
        .select([cell_col, intbc_col, '_rcol', '_consensus_value'])
    )

    # Pivot back to wide format
    result = consensus.pivot(
        on='_rcol',
        index=[cell_col, intbc_col],
        values='_consensus_value'
    )

    # Add preserved columns (take first value per group)
    if preserve_cols:
        preserved = (
            df_idx.group_by([cell_col, intbc_col])
            .agg([pl.col(c).first().alias(c) for c in preserve_cols])
        )
        result = result.join(preserved, on=[cell_col, intbc_col], how='left')

    return result


def _collapse_mode_umi(
    df: pl.DataFrame,
    rcols: list[str],
    cell_col: str,
    intbc_col: str,
) -> pl.DataFrame:
    """Whole-UMI consensus: pick the most frequent allele pattern.

    Instead of voting independently per cutsite, this treats the entire
    r-column vector of each UMI as an atomic unit. For each (cell, intBC)
    group, it finds the most common allele pattern across UMIs.

    Ranking (descending): total_molecule_len > n_umis > total_readCount.
    sum(mol_len) naturally encodes both frequency and quality — many long
    molecules dominate over many short chimeras or few long artifacts.

    This avoids the fragmentation problem where ONT errors in flanking
    context create many unique per-site strings that dilute mode() votes.
    """
    # Build a composite key from all r columns to identify unique patterns
    # Use "|" as separator — safe since allele strings don't contain it
    has_seq = 'Seq' in df.columns
    has_readcount = 'readCount' in df.columns

    df_keyed = df.with_columns(
        pl.concat_str(rcols, separator='|').alias('_pattern'),
        *([pl.col('Seq').str.len_chars().alias('_mol_len')] if has_seq else []),
    )

    # Aggregate: count UMIs per pattern, sum molecule length and readCount as tiebreakers
    agg_exprs = [
        pl.len().alias('_n_umis'),
    ]
    if has_seq:
        agg_exprs.append(pl.col('_mol_len').sum().alias('_total_mol_len'))
    if has_readcount:
        agg_exprs.append(pl.col('readCount').sum().alias('_total_reads'))

    pattern_counts = (
        df_keyed
        .group_by([cell_col, intbc_col, '_pattern'])
        .agg(agg_exprs)
    )

    # For each (cell, intBC), pick the pattern with the highest total
    # molecule length, breaking ties by UMI count, then readCount
    sort_cols = []
    if has_seq:
        sort_cols.append('_total_mol_len')
    sort_cols.append('_n_umis')
    if has_readcount:
        sort_cols.append('_total_reads')

    best_patterns = (
        pattern_counts
        .sort(sort_cols, descending=True)
        .group_by([cell_col, intbc_col])
        .first()
        .select([cell_col, intbc_col, '_pattern'])
    )

    # Explode the pattern back into individual r columns
    result = best_patterns.with_columns(
        pl.col('_pattern').str.split('|').alias('_parts')
    )

    # Assign each element back to its r column
    for i, rcol in enumerate(rcols):
        result = result.with_columns(
            pl.col('_parts').list.get(i).alias(rcol)
        )

    # Replace "null" strings that came from concat_str with actual nulls
    result = result.with_columns([
        pl.when(pl.col(c) == 'null').then(pl.lit(None)).otherwise(pl.col(c)).alias(c)
        for c in rcols
    ])

    result = result.drop(['_pattern', '_parts'])

    return result


def flag_spanning_deletions(pl_allele: pl.DataFrame, rcols: list[str], mode: str = 'unedited') -> pl.DataFrame:
    """Handle alleles spanning multiple cutsites.

    Detects spanning deletions by checking for identical non-empty allele values
    at *adjacent* cutsite positions within the same row. The same allele at
    non-adjacent positions is treated as independent edits and left untouched.

    Args:
        pl_allele: DataFrame with r1, r2, ... columns containing allele values
        rcols: List of r column names in order (e.g., ['r001', 'r002', 'r003', ...])
        mode: How to handle spanning deletions:
            - 'unedited': mark as unedited ("None" string, Cassiopeia state 0)
            - 'missing': mark as missing data (null)
            - 'use': leave as-is (keep the deletion allele value)
            Note: 'both' is resolved by the caller (_build_trees) before reaching here.

    Returns:
        DataFrame with spanning deletions handled according to mode
    """
    valid_modes = ('unedited', 'missing', 'use')
    if mode not in valid_modes:
        raise ValueError(
            f"flag_spanning_deletions: unsupported mode {mode!r}. "
            f"Expected one of {valid_modes}. "
            f"If the caller is resolving a higher-level value like 'both', "
            f"do the expansion before calling this helper."
        )

    if not rcols or mode == 'use':
        return pl_allele

    replacement_value = 'None' if mode == 'unedited' else None

    # For each row, mark positions where the allele equals the adjacent
    # neighbor's allele (non-empty, non-"None" values only)
    # A position is spanning if it matches its left OR right neighbor.
    is_spanning_exprs = []
    for i, rcol in enumerate(rcols):
        neighbors = []
        if i > 0:
            neighbors.append(pl.col(rcol) == pl.col(rcols[i - 1]))
        if i < len(rcols) - 1:
            neighbors.append(pl.col(rcol) == pl.col(rcols[i + 1]))

        # Combine neighbor checks with OR
        is_adj_match = neighbors[0]
        for n in neighbors[1:]:
            is_adj_match = is_adj_match | n

        # Only flag non-null, non-empty, non-"None" alleles
        is_real_allele = (
            pl.col(rcol).is_not_null()
            & (pl.col(rcol) != '')
            & (~pl.col(rcol).str.contains('(?i)^None$'))
        )

        is_spanning_exprs.append(
            (is_real_allele & is_adj_match).alias(f'_span_{rcol}')
        )

    df = pl_allele.with_columns(is_spanning_exprs)

    # Replace spanning positions according to mode
    replace_exprs = []
    for rcol in rcols:
        replace_exprs.append(
            pl.when(pl.col(f'_span_{rcol}'))
            .then(pl.lit(replacement_value))
            .otherwise(pl.col(rcol))
            .cast(pl.Utf8)
            .alias(rcol)
        )

    result = (
        df
        .with_columns(replace_exprs)
        .drop([f'_span_{rcol}' for rcol in rcols])
    )

    return result


def pivot_alleles_wide(allele_table: pl.DataFrame, rcols: list[str]) -> pl.DataFrame:
    """Pivot allele table from long (one row per cell x intBC) to wide (one row per cell).

    Creates columns named like 'TCTGAA..._r1', 'TCTGAA..._r2', etc.
    When multiple UMIs exist per (cell, intBC), takes the mode per cutsite.
    """
    if allele_table.height == 0:
        return pl.DataFrame({'cellBC': []})

    intbc_order = allele_table.get_column('intBC').unique(maintain_order=True).sort().to_list()
    col_order = [f'{ibc}_{r}' for ibc in intbc_order for r in rcols]

    df = (
        allele_table
        .select('cellBC', 'intBC', *rcols)
        .unpivot(index=['cellBC', 'intBC'], on=rcols)
        .with_columns((pl.col('intBC') + '_' + pl.col('variable')).alias('col'))
        .group_by('cellBC', 'col')
        .agg(pl.col('value').mode().first())
        .pivot(index='cellBC', on='col', values='value')
    )
    present = [c for c in col_order if c in df.columns]
    return df.select('cellBC', *present)


def add_depth(
    tdata,
    tree_key: str = 'nj',
    depth_key: str = 'depth',
    normalized_key: str = 'normalized_depth',
) -> None:
    """Compute depth and normalized depth (depth / n_leaves) on TreeData.

    Args:
        tdata: TreeData object with obst containing tree topology
        tree_key: Key in tdata.obst for the tree
        depth_key: Key in tdata.obs for the depth
        normalized_key: Key to store depth / n_leaves in tdata.obs
    """
    import pycea as pc
    if depth_key not in tdata.obs.columns:
        pc.pp.add_depth(tdata, tree=tree_key, key_added=depth_key)

    n_leaves = len(tdata.obs)
    if n_leaves > 0:
        tdata.obs[normalized_key] = tdata.obs[depth_key] / n_leaves


def add_weighted_depth(tdata, tree_key: str, depth_key: str = 'depth', weight: str = 'length') -> None:
    """Compute node depths using weighted edge lengths (e.g., from convexml).

    Unlike add_depth/pycea.pp.add_depth which counts edges, this sums branch
    lengths along the path from root to each node using Dijkstra.
    """
    import networkx as nx
    tree = tdata.obst[tree_key]
    root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]
    depths = nx.single_source_dijkstra_path_length(tree, root, weight=weight)
    nx.set_node_attributes(tree, depths, depth_key)
    tdata.obs[depth_key] = tdata.obs.index.map(lambda x: depths.get(x, None))
    tdata.obs = tdata.obs.copy()  # defragment after column insertions


def write_done(outdir: Path, reason: str = None) -> None:
    """Write done marker file. If reason, writes done_{reason} instead."""
    if reason:
        (outdir / f'done_{reason}').write_text(reason)
    else:
        (outdir / 'done').write_text('success')


def generate_intbc_whitelist(
    ldf: pl.LazyFrame,
    min_umis: int = 100,
    min_proportion_of_sample: float = 0.05,
    min_ratio_to_max: float = 0.1,
    modality: str = 'single-molecule',
    cbc_len: int = 16,
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Generate valid intBC whitelist using per-sample adaptive thresholds.

    For single-molecule: groups by sbc (sample barcode)
    For single-cell: groups by cbc (cell barcode, extracted from compound umi)

    Keeps intBCs that meet BOTH:
    - umis >= max(min_umis, sample_total * min_proportion_of_sample)
    - ratio_to_max >= min_ratio_to_max

    Args:
        ldf: LazyFrame with 'umi', 'sbc', 'intBC' columns
        min_umis: Absolute minimum UMI count threshold
        min_proportion_of_sample: Minimum proportion of sample total UMIs (e.g., 0.05 = 5%)
        min_ratio_to_max: Minimum ratio to the largest intBC per sbc (e.g., 0.1 = 10%)
        modality: 'single-cell' or 'single-molecule'
        cbc_len: Cell barcode length for single-cell (default 16 for 10x)
        logger: Optional logger for info messages

    Returns:
        DataFrame with columns: group_id, intBC, umis, reads, ratio_to_max
        (group_id is sbc for single-molecule, cbc for single-cell)
    """
    # Determine grouping column based on modality
    if modality == 'single-cell':
        # Extract cbc from compound umi (first cbc_len characters)
        work_ldf = ldf.with_columns(
            pl.col('umi').str.slice(0, cbc_len).alias('cbc'),
            pl.col('umi').str.slice(cbc_len).alias('umi_only'),
        )
        group_col = 'cbc'
        umi_col = 'umi_only'
        group_label = 'cbcs'
    else:
        work_ldf = ldf
        group_col = 'sbc'
        umi_col = 'umi'
        group_label = 'sbcs'

    result = (
        work_ldf
        .select(umi_col, group_col, 'intBC')
        .group_by(group_col, 'intBC').agg(
            pl.col(umi_col).n_unique().alias('umis'),
            pl.len().alias('reads')
        )
        .with_columns(
            sample_total_umis=pl.col('umis').sum().over(group_col),
            max_umis=pl.col('umis').max().over(group_col),
        )
        .with_columns(
            min_umis_threshold=pl.max_horizontal(
                pl.col('sample_total_umis') * min_proportion_of_sample,
                pl.lit(min_umis)
            ),
            ratio_to_max=pl.col('umis') / pl.col('max_umis'),
        )
        .filter(
            (pl.col('umis') >= pl.col('min_umis_threshold')) &
            (pl.col('ratio_to_max') >= min_ratio_to_max)
        )
        .select(pl.col(group_col).alias('group_id'), 'intBC', 'umis', 'reads', 'ratio_to_max')
        .sort(['group_id', 'umis'], descending=[False, True])
        .collect()
    )

    if logger:
        n_groups = result['group_id'].n_unique()
        n_intbcs = result.height
        logger.info(f"Generated whitelist: {n_intbcs} valid intBCs across {n_groups} {group_label}")

    return result


def filter_by_whitelist(
    ldf: pl.LazyFrame,
    whitelist: pl.DataFrame,
    modality: str = 'single-molecule',
    cbc_len: int = 16,
    logger: Optional[CustomLogger] = None,
) -> pl.LazyFrame:
    """
    Filter LazyFrame to only include (group, intBC) pairs in whitelist.

    For single-molecule: filters by (sbc, intBC) pair
    For single-cell: filters by (cbc, intBC) pair where cbc is extracted from umi

    Args:
        ldf: LazyFrame with 'umi', 'sbc', and 'intBC' columns
        whitelist: DataFrame with 'group_id' and 'intBC' columns (or legacy 'sbc')
        modality: 'single-cell' or 'single-molecule'
        cbc_len: Cell barcode length for single-cell (default 16 for 10x)
        logger: Optional logger for info messages

    Returns:
        Filtered LazyFrame
    """
    # Handle both new format (group_id) and legacy format (sbc)
    if 'group_id' in whitelist.columns:
        group_col_whitelist = 'group_id'
    elif 'sbc' in whitelist.columns:
        group_col_whitelist = 'sbc'
    else:
        # Global intBC whitelist (no grouping)
        valid_intbcs = set(whitelist['intBC'].to_list())
        n_before = ldf.select(pl.col('intBC').n_unique()).collect().item()
        filtered = ldf.filter(pl.col('intBC').is_in(valid_intbcs))
        n_after = filtered.select(pl.col('intBC').n_unique()).collect().item()
        if logger:
            logger.info(f"Whitelist filter: {n_before} -> {n_after} intBCs (global)")
        return filtered

    # Determine grouping column in ldf based on modality
    if modality == 'single-cell':
        # Extract cbc from compound umi for joining
        work_ldf = ldf.with_columns(
            pl.col('umi').str.slice(0, cbc_len).alias('_filter_group')
        )
        group_label = 'cbc'
    else:
        work_ldf = ldf.with_columns(
            pl.col('sbc').alias('_filter_group')
        )
        group_label = 'sbc'

    # Prepare whitelist for join
    whitelist_for_join = whitelist.select(
        pl.col(group_col_whitelist).alias('_filter_group'),
        'intBC'
    ).lazy()

    n_before = work_ldf.select(pl.struct('_filter_group', 'intBC').n_unique()).collect().item()
    filtered = work_ldf.join(
        whitelist_for_join,
        on=['_filter_group', 'intBC'],
        how='semi'
    ).drop('_filter_group')
    n_after = filtered.select(pl.struct(group_label if modality == 'single-cell' else 'sbc', 'intBC').n_unique()).collect().item() if modality != 'single-cell' else filtered.with_columns(pl.col('umi').str.slice(0, cbc_len).alias('cbc')).select(pl.struct('cbc', 'intBC').n_unique()).collect().item()

    if logger:
        logger.info(f"Whitelist filter: {n_before} -> {n_after} ({group_label}, intBC) pairs")

    return filtered


def plug_cassiopeia(
        ldf: pl.LazyFrame,
        ann_intbc_mod: pl.DataFrame,
        workdir: Path|str ='.',
        logger: None|CustomLogger=  None,
        barcode_interval: List|Tuple = (0, 7),
        cutsite_locations: List =  [40, 67, 94, 121, 148, 175, 202, 229, 256, 283],
        cutsite_width: int = 12,
        context: bool = True,
        context_size: int = 50,
        gap_open_penalty: Optional[int] = None,
        gap_extend_penalty: Optional[int] = None,
        alignment_method: str = 'global',
        # Partition filtering parameters
        intbc_whitelist_path: Optional[str] = None,
        min_molecules_per_group: int = 10,
        min_proportion_of_sample: float = 0.02,
        min_ratio_to_max: float = 0.1,
        top_n_cells: Optional[int] = None,
        # Modality parameters
        modality: str = 'single-molecule',
        cbc_len: int = 16,
        ) -> StepResults:
    """
    readName - A unique identifier for each row/sequence
    cellBC - The cell barcode (cbc for single-cell, sbc for single-molecule)
    UMI - The UMI (Unique Molecular Identifier)
    readCount - The number of reads for this sequence
    seq - The actual sequence to be aligned

    For single-cell data:
        - cellBC is extracted from compound umi (first cbc_len chars)
        - UMI is the remaining part of compound umi (umi_only)
        - Whitelist grouping uses cbc instead of sbc
    """
    import cassiopeia as cas

    allele_params = {
        'barcode_interval': barcode_interval,
        'cutsite_locations': cutsite_locations,
        'cutsite_width': cutsite_width,
        'context': context,
        'context_size': context_size,
    }

    # Build cassiopeia columns based on modality
    if modality == 'single-cell':
        # For single-cell: extract cbc and umi_only from compound umi
        cass_ldf = (
            ldf.with_columns(
                pl.col('umi').str.slice(0, cbc_len).alias('cbc'),
                pl.col('umi').str.slice(cbc_len).alias('umi_only'),
            )
            .with_columns(
                readName=pl.col('umi'),
                cellBC=pl.col('cbc'),           # Use actual cell barcode
                UMI=pl.col('umi_only'),          # Use actual UMI (not compound)
                readCount=pl.col('reads'),
                seq=pl.col('intBC')+pl.col('contig'),
            )
            .select('readName', 'cellBC', 'UMI', 'readCount', 'seq', 'intBC', 'sbc', 'cbc', 'umi')
            .join(ann_intbc_mod.lazy(), left_on='intBC', right_on='intBC', how='inner')
        )
    else:
        # For single-molecule: use compound umi as before
        cass_ldf = (
            ldf.with_columns(
                readName=pl.col('umi'),
                cellBC=pl.col('umi'),
                UMI=pl.col('umi'),
                readCount=pl.col('reads'),
                seq=pl.col('intBC')+pl.col('contig'),
            )
            .select('readName', 'cellBC', 'UMI', 'readCount', 'seq', 'intBC', 'sbc')
            .join(ann_intbc_mod.lazy(), left_on='intBC', right_on='intBC', how='inner')
        )

    # === PRE-FILTERING ===
    filter_metrics = {}

    # Option A: Load pre-generated whitelist
    if intbc_whitelist_path is not None:
        whitelist = read_file(intbc_whitelist_path)
        cass_ldf = filter_by_whitelist(cass_ldf, whitelist, modality=modality, cbc_len=cbc_len, logger=logger)
        filter_metrics['whitelist_source'] = 'file'
        filter_metrics['whitelist_path'] = str(intbc_whitelist_path)

    # Option B: Generate whitelist inline if thresholds set
    elif min_molecules_per_group > 0 or min_proportion_of_sample > 0 or min_ratio_to_max > 0:
        whitelist = generate_intbc_whitelist(
            ldf,  # Use original ldf, not cass_ldf (before cassiopeia columns added)
            min_umis=min_molecules_per_group,
            min_proportion_of_sample=min_proportion_of_sample,
            min_ratio_to_max=min_ratio_to_max,
            modality=modality,
            cbc_len=cbc_len,
            logger=logger,
        )
        cass_ldf = filter_by_whitelist(cass_ldf, whitelist, modality=modality, cbc_len=cbc_len, logger=logger)
        filter_metrics['whitelist_source'] = 'generated'
        filter_metrics['valid_intbc_count'] = whitelist.height

    # Cell-level filter: keep only top N cells by total UMI count (single-cell only)
    if top_n_cells is not None and modality == 'single-cell':
        top_cells = (
            cass_ldf.group_by('cbc')
            .len()
            .sort('len', descending=True)
            .head(top_n_cells)
            .select('cbc')
        )
        n_before = cass_ldf.select(pl.col('cbc').n_unique()).collect().item()
        cass_ldf = cass_ldf.join(top_cells.lazy(), on='cbc', how='semi')
        n_after = cass_ldf.select(pl.col('cbc').n_unique()).collect().item()
        if logger:
            logger.info(f"top_n_cells filter: {n_before} -> {n_after} cells (top {top_n_cells})")
        filter_metrics['top_n_cells'] = top_n_cells
        filter_metrics['cells_before_topn'] = n_before
        filter_metrics['cells_after_topn'] = n_after

    # Filter intBCs with too few UMIs within each mod (across all cells)
    # Always applied, including when using a whitelist
    if min_molecules_per_group > 0:
        n_before = cass_ldf.select(pl.len()).collect().item()
        cass_ldf = cass_ldf.filter(
            pl.len().over(['intBC', 'mod']) >= min_molecules_per_group
        )
        n_after = cass_ldf.select(pl.len()).collect().item()
        if logger:
            logger.info(f"min_molecules_per_group filter ({min_molecules_per_group}): {n_before} -> {n_after} rows")
        filter_metrics['rows_before_partition_filter'] = n_before
        filter_metrics['rows_after_partition_filter'] = n_after

    # === END PRE-FILTERING ===

    res = []
    umi_tables = []
    nones = 0

    # Alignment threading: use half of available CPUs
    import os
    _available_cpus = len(os.sched_getaffinity(0))
    _align_threads = max(1, _available_cpus // 2)
    if logger:
        logger.warning(f"Alignment using {_align_threads} threads (half of {_available_cpus} available CPUs)")

    # Determine partition columns based on modality
    if modality == 'single-cell':
        partition_cols = ['intBC', 'mod']
    else:
        partition_cols = ['intBC', 'mod', 'sbc']

    for partition_key, queries in cass_ldf.collect().partition_by(*partition_cols, as_dict=True).items():
        if modality == 'single-cell':
            intBC, mod = partition_key
            group_id = None
        else:
            intBC, mod, group_id = partition_key
        if mod is None:
            nones+=1
        else:
            if logger:
                logger.debug(f"{mod=} {intBC=} {queries.shape=}")

            # Build alignment kwargs - only include gap penalties if explicitly set
            align_kwargs = {
                'queries': queries.to_pandas(),
                'ref_filepath': f'{workdir}/{mod}.fasta',
                'n_threads': _align_threads,
                'method': alignment_method,
            }
            if gap_open_penalty is not None:
                align_kwargs['gap_open_penalty'] = gap_open_penalty
            if gap_extend_penalty is not None:
                align_kwargs['gap_extend_penalty'] = gap_extend_penalty

            umi_table = cas.pp.align_sequences(**align_kwargs)
    
            allele_table =  cas.pp.call_alleles(
                            umi_table,
                            ref_filepath = f'{workdir}/{mod}.fasta',
                            **allele_params,
                        )


            # Enrich allele columns with actual insertion sequences from CIGAR
            # allele_table already contains Seq and CIGAR columns from umi_table
            pl_allele = pl.DataFrame(allele_table)
            rcols = pl_allele.select(pl.col('^r\\d+$')).columns

            # Zero-pad cutsite column names (r1 -> r001, ... r15 -> r015) so
            # default lexicographic sort gives the biological order. Cassiopeia
            # itself emits r1/r2/... unpadded.
            if rcols:
                rcols_zf = {c: f'r{int(c[1:]):03d}' for c in rcols}
                pl_allele = pl_allele.rename(rcols_zf)
                rcols = list(rcols_zf.values())

            if rcols and 'Seq' in pl_allele.columns and 'CIGAR' in pl_allele.columns:
                pl_allele = pl_allele.with_columns([
                    pl.col(col).cigar.enrich_insertions(pl.col('Seq'), pl.col('CIGAR'))
                    for col in rcols
                ])

            # Fill any remaining null values in r columns with "None" string
            # (nulls come from cas.pp.call_alleles when it can't determine a cutsite allele)
            # Cassiopeia treats "None" as unedited/wildtype (state 0)
            if rcols:
                pl_allele = pl_allele.with_columns([
                    pl.col(col).fill_null("None") for col in rcols
                ])

            allele_table = pl_allele.to_pandas()

            #replace intBC for real intBC since Cassiopeia does something I don't understand yet
            # include colums dropped by cass?
            allele_table['intBC'] = intBC
            allele_table['mod'] = mod
            if modality != 'single-cell':
                allele_table['sbc'] = group_id

            umi_tables.append(umi_table.copy())
            #self.xp.logger.info(f"found {umi_table['intBC'].n_unique()} integrations for {mod}")

            res.append(
                pl.DataFrame(allele_table).with_columns(mod=pl.lit(mod), mols=allele_table.shape[0])
                )

    # Concatenate all partition results. Cassiopeia returns pandas DataFrames
    # with per-partition index naming — when the index has the name 'readName'
    # that becomes a real column; when it is a default RangeIndex polars
    # conversion can materialise an 'index' / 'level_0' column instead. Drop
    # only the truly artificial pandas-index leftovers ('index', 'level_0')
    # and keep 'readName' as legitimate data; diagonal_relaxed fills nulls in
    # partitions that don't have it.
    if res:
        drop_cols = ('index', 'level_0')
        res = [df.drop([c for c in drop_cols if c in df.columns]) for df in res]
        alleles_pl = pl.concat(res, how='diagonal_relaxed')
    else:
        alleles_pl = pl.DataFrame()

    return StepResults(
            results={
                "alleles_pl": alleles_pl,                      # Raw per-UMI alleles
                "alleles_pd": umi_tables,
            },
            metrics={"partitions_processed": len(res), "none_mod_skipped": nones, **filter_metrics}
    )

def save_ref_to_fasta(refs: pl.DataFrame, out_dir: str|Path = '.', field: str = 'mod') -> None:
    """Write reference sequences to individual FASTA files.

    Args:
        refs: DataFrame with 'mod' and 'seq' columns
        out_dir: Output directory for FASTA files
        field: Column to use for file naming and filtering
    """
    if refs.height == 0:
        raise ValueError("refs DataFrame is empty - no references to save")

    unique_values = refs.get_column(field).unique()
    if unique_values.null_count() > 0:
        raise ValueError(f"refs DataFrame has null values in '{field}' column")

    for i in unique_values:
        filtered = refs.filter(pl.col(field) == i)
        if filtered.height == 0:
            raise ValueError(f"No rows found for {field}={i}")

        fasta_content = filtered.dna.to_fasta(read_id_col=field, read_col='seq').get_column('seq_fasta')[0]
        out_path = Path(out_dir) / f"{i}.fasta"
        with open(out_path, 'wt') as out:
            out.write(fasta_content)

def kmer_classify_cassettes(ldf, refs: pl.DataFrame, K: int = 25) -> StepResults:
    """ Annotates contigs based on the top occurring mod based on kmer matches """

    k_ref = refs.kmer.explode_kmers(k=K, seq_field='seq')

    cont_k = (
        ldf
        .kmer.explode_kmers(k=K, seq_field='contig', only_unique=False)
        .filter(pl.col('kmer').is_in(set(k_ref.get_column('kmer').to_list())))
        .join(k_ref.drop('seq').lazy(), left_on='kmer', right_on='kmer', how='inner')
        .group_by('intBC', 'mod', 'sbc')
                .agg(fps=pl.col('umi').n_unique())
        .group_by('intBC', 'mod')
                .agg(fps=pl.col('fps').sum())
        .sort('fps', descending=True)
        .group_by('intBC', maintain_order=True)
        .first()
        .collect(engine="streaming")
        )
    # TODO return how many intBCs didn't get a mod
    # TODO return the number of ties 
    return StepResults(results={'ann_intbc_mod':cont_k},
                       metrics={'n_ann_intbc':cont_k.height})

def parse_contigs(
        ldf: pl.LazyFrame,
        int_anchor1: str,
        int_anchor2 : str,
        sbc_dict: Optional[Dict] = None,
        annotation:str = 'sbc',
        ) ->StepResults:
    """Extracts integration barcodes (intBC) and sample barcodes (sbc)."""
    
    # 1. Extract/use sbcs
    if sbc_dict is not None:
        ldf = ldf.with_columns(pl.col('sbc').replace(sbc_dict).alias(annotation))
    
    # 2. Extract intBC and trim contig to remove the intBC region
    ldf = (
        ldf
        .pp.extract_intbc(int_anchor1, int_anchor2, seq_col='contig')
        .with_columns(
            pl.col('contig').str.replace(f'.*{int_anchor2}', '').alias('contig')
        )
    )
    
    # Get count lazily - will be computed when needed
    n_contigs = ldf.select(pl.len()).collect().item()
    return StepResults(results={'ldf':ldf},
                       metrics={'n_parsed_contigs': n_contigs})


def generate_refs_from_fasta(refs_fasta_path: str|Path, anchor1: str, anchor2: str) -> pl.DataFrame:
    """ Given a FASTA file with references it generates a trimmed version of the cassettes for aligment

    Supports FASTA headers in multiple formats:
    - ">5mer" → mod_5mer
    - ">v5mer_xxx" → mod_5mer
    - ">10mer" → mod_10mer
    """

    refs =  (
            pl.read_csv(refs_fasta_path,
                        has_header=False,
                        )
            .unstack(step=2, how='horizontal')
            .rename({"column_1_0":"mod", "column_1_1":"seq"})
            .with_columns(
                # Handle multiple header formats:
                # ">5mer" -> "5mer" -> "mod_5mer"
                # ">v5mer_xxx" -> "5mer" -> "mod_5mer"
                pl.concat_str([
                    pl.lit("mod_"),
                    pl.col('mod')
                        .str.strip_chars(">")  # Remove leading >
                        .str.extract(r"^v?(\d+mer)", 1)  # Extract Nmer pattern (optional v prefix)
                ]).alias('mod'),
                pl.col('seq').str.to_uppercase())
            # trim sequence
            .with_columns(
                pl.col('seq').str.replace(f".+?({anchor1})", anchor1).str.replace(f"({anchor2}).+?$", anchor2)
                )
            )

    # Validate we got valid mod values
    null_mods = refs.filter(pl.col('mod').is_null()).height
    if null_mods > 0:
        raise ValueError(f"Failed to parse {null_mods} FASTA headers. Expected format: >5mer, >10mer, or >v5mer_xxx")

    return refs


# ============================================================================
# Segment-based allele table generation
# ============================================================================

CASSETTE_START_MARKER = "_CASSETTE_START_"
CASSETTE_END_MARKER = "_CASSETTE_END_"
TARGET_ORDER = ["RNF2", "HEK3", "EMX1"]  # Order within each triplet

CASSETTE_CONFIGS = {
    "5mer": {"n_metas": 4, "n_targets": 12},
    "10mer": {"n_metas": 9, "n_targets": 27},
    "20mer": {"n_metas": 19, "n_targets": 57},
}


def filter_segment_types_for_cassette(
    segment_types: pl.DataFrame,
    metas: pl.DataFrame,
    cassette_type: str,
    min_freq_threshold: float = 0.7,
    total_molecules: int | None = None,
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Filter segment_types to only include expected metas for the cassette type.

    For a 5mer cassette, only META01-META04 are expected. Any other metas
    (META05+) are filtered out UNLESS they appear at high frequency (>min_freq_threshold),
    which may indicate they are real (e.g., a different cassette version).

    Args:
        segment_types: DataFrame with start_meta, end_meta, count columns
        metas: DataFrame with feature column (META01, META02, etc.)
        cassette_type: One of "5mer", "10mer", "20mer"
        min_freq_threshold: Minimum frequency (0-1) for an unexpected meta to be kept.
                           Default 0.7 (70%). Set to 1.0 to strictly filter.
        total_molecules: Total number of molecules for frequency calculation.
                        If None, uses sum of counts in segment_types.
        logger: Optional logger for debugging

    Returns:
        Filtered segment_types DataFrame
    """
    if cassette_type not in CASSETTE_CONFIGS:
        if logger:
            logger.warning(f"Unknown cassette type '{cassette_type}', skipping meta filtering")
        return segment_types

    n_metas = CASSETTE_CONFIGS[cassette_type]["n_metas"]

    # Get ordered meta names from metas DataFrame
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    all_meta_names = meta_only['feature'].to_list()

    # Expected metas for this cassette type (first n_metas)
    expected_metas = set(all_meta_names[:n_metas])
    # Add edge markers as always valid
    expected_metas.add(CASSETTE_START_MARKER)
    expected_metas.add(CASSETTE_END_MARKER)

    if logger:
        logger.debug(f"Cassette {cassette_type}: expecting {n_metas} metas: {sorted(expected_metas)}")

    # Calculate frequency of each meta across all transitions
    if total_molecules is None:
        total_molecules = segment_types['count'].sum()

    if total_molecules == 0:
        return segment_types

    # Find all unique metas in segment_types
    all_starts = set(segment_types['start_meta'].unique().to_list())
    all_ends = set(segment_types['end_meta'].unique().to_list())
    all_observed_metas = all_starts | all_ends

    # Check for unexpected metas with high frequency
    unexpected_metas = all_observed_metas - expected_metas
    high_freq_unexpected = set()

    for meta in unexpected_metas:
        # Count molecules with this meta (in either start or end)
        meta_count = segment_types.filter(
            (pl.col('start_meta') == meta) | (pl.col('end_meta') == meta)
        )['count'].sum()
        freq = meta_count / total_molecules

        if freq >= min_freq_threshold:
            high_freq_unexpected.add(meta)
            if logger:
                logger.info(f"Keeping unexpected meta '{meta}' - appears in {freq:.1%} of molecules (>= {min_freq_threshold:.0%} threshold)")
        elif logger:
            logger.debug(f"Filtering out meta '{meta}' - appears in {freq:.1%} of molecules (< {min_freq_threshold:.0%} threshold)")

    # Final whitelist: expected + high-frequency unexpected
    whitelist = expected_metas | high_freq_unexpected

    # Filter segment_types to only include transitions with whitelisted metas
    filtered = segment_types.filter(
        pl.col('start_meta').is_in(list(whitelist)) &
        pl.col('end_meta').is_in(list(whitelist))
    )

    n_removed = segment_types.height - filtered.height
    if n_removed > 0 and logger:
        logger.info(f"Filtered {n_removed} transitions with unexpected metas (kept {filtered.height})")

    return filtered


def build_target_position_map(
    metas_df: pl.DataFrame,
    cassette_type: str = "5mer",
) -> Dict[Tuple[str, str], List[Tuple[str, int]]]:
    """
    Build mapping from (start_meta, end_meta) segment boundaries to TARGET positions.

    Logic:
    - For 5mer: 4 METAs (META01-META04), creating 4 segments each with 3 TARGETs = 12 positions
    - Segment (_CASSETTE_START_, META01) → positions r1, r2, r3 (RNF2, HEK3, EMX1)
    - Segment (META01, META02) → positions r4, r5, r6
    - etc.

    Returns:
        Dict mapping (start_meta, end_meta) -> [(target_name, position_index), ...]
    """
    raise NotImplementedError("TODO: implement with polars")


def segments_to_allele_table(
    segments_df: pl.DataFrame,
    metas_df: pl.DataFrame,
    int_anchor1: str,
    int_anchor2: str,
    cassette_type: str = "5mer",
    min_consensus_support: float = 0.5,
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Transform segmentation table to cassiopeia allele table format.

    Input columns: [sbc, umi, start_meta, end_meta, segment_seq]

    Output columns:
    - cellBC, UMI, readCount, intBC, sbc, mod
    - r1, r2, r3, ... rN (one per TARGET position)

    Logic:
    1. Extract intBC using: pl.col('segment_seq').str.extract(f'{int_anchor1}(.+?){int_anchor2}', 1)
    2. Map (start_meta, end_meta) to TARGET positions using build_target_position_map()
    3. For each TARGET in segment, extract insertion using flanks from metas_df
       - left_flank and right_flank columns identify TARGET boundaries
       - Insertion = sequence between left_flank and right_flank
       - Empty string = wild-type, None = flanks not found
    4. Group by (sbc, umi) and consensus across multiple observations
    5. Pivot to wide format: r1, r2, r3, ...

    Values in rN columns:
    - "" (empty): Wild-type, no insertion
    - "ACTGT": Insertion sequence (lineage mark)
    - None: Missing data
    """
    raise NotImplementedError("TODO: implement with polars")


@dataclass
class CassiopeiaConfig(ExtensionConfig):
    # Required fields
    int_anchor1: str
    int_anchor2: str
    force: bool = False

    # Optional fields with defaults
    sbc_dict: Optional[Dict] = None
    annotation: str = 'sbc'
    refs_fasta_path: Optional[str] = None
    anchor1: Optional[str] = None
    anchor2: Optional[str] = None

    barcode_interval: Tuple[int, int] = (0, 7)
    cutsite_locations: List[int] = field(default_factory=lambda: [40, 67, 94, 121, 148, 175, 202, 229, 256, 283])
    cutsite_width: int = 12
    context: bool = True
    context_size: int = 50

    # Alignment parameters for cas.pp.align_sequences
    # When None, Cassiopeia's defaults are used
    gap_open_penalty: Optional[int] = None
    gap_extend_penalty: Optional[int] = None
    alignment_method: str = 'global'  # 'local' (Smith-Waterman) or 'global' (Needleman-Wunsch)

    # Spanning deletion handling: how to treat deletions that span multiple cutsites
    # - 'unedited': mark as unedited ("None" string, Cassiopeia state 0)
    # - 'missing': mark as missing data (null)
    # - 'use': leave as-is (keep the deletion allele value)
    # - 'both': run trees for both 'unedited' and 'use' (separate output dirs)
    spanning_deletions: str = 'unedited'

    # Segmented allele extraction fields
    metas_flanks_csv: Optional[str] = None  # Path to PEtracer_metas_flanks.csv
    cassette_type: str = "5mer"  # Options: 5mer, 10mer, 20mer
    min_consensus_support: float = 0.5  # Minimum support for consensus calls

    # Partition filtering parameters - reduce partition explosion from sequencing errors
    # See generate_intbc_whitelist() for threshold logic
    intbc_whitelist_path: Optional[str] = None  # Parquet with 'intBC', 'sbc' columns
    min_molecules_per_group: int = 10           # Absolute minimum UMIs per (sbc, intBC) group
    min_proportion_of_sample: float = 0.02      # % of sample total (e.g., 0.02 = 2%)
    min_ratio_to_max: float = 0.1               # % of largest group per sbc (e.g., 0.1 = 10%)
    top_n_cells: Optional[int] = None            # Keep top N cells by total UMI count (single-cell only)

    # Molecule size filtering
    min_molecule_len: Optional[int] = None      # Min sequence length to include in allele table (None = no filter)

    # Single-cell collapse parameters (SC mode collapses multiple UMIs per cell to one allele per (cell, intBC))
    collapse_method: str = 'mode_site'          # 'mode_site' = per-cutsite mode voting, 'mode_umi' = most frequent whole-UMI pattern
    min_umis_per_cell: int = 4                  # Minimum UMIs for valid consensus (cells with fewer are filtered)
    min_umi_agreement: Optional[float] = 0.5     # Min fraction for consensus (below -> missing/None, null = no threshold)
                                                 # Only used by mode_site

    # Tree generation
    solver: str = 'nj'                          # nj | vanilla | mcgs
    min_cells_for_tree: int = 10
    skip_branch_lengths: bool = False
    min_branch_length: float = 0.01
    tree_mode: str = 'auto'                     # auto | sc | sm
                                                # auto: sc if single-cell modality else sm
    tree_group_by: list = field(default_factory=lambda: ['intBC', 'sbc'])
                                                # SM mode only: columns to partition jobs by
                                                # ['intBC'] = one tree per intBC (pool sbcs)
                                                # ['intBC', 'sbc'] = one tree per (intBC, sbc)
                                                # Ignored in SC mode
    allele_rep_thresh: float = 1.0              # allele representation threshold for character matrix
    tree_test_mode: bool = False                  # Enable subsampling for quick testing
    tree_subsample: list = field(default_factory=list)
    # Generic subsampling factors. Each entry is a dict:
    #   column: str   — column to group by (e.g. 'intBC', 'cellBC', 'UMI')
    #   top_x: int    — keep top X groups ranked by size
    #   sample_n: int — then sample N from that pool (optional, keeps all if omitted)

    # LSF submission for tree building (follows dorado pattern)
    tree_use_lsf: bool = True
    tree_lsf_queue: str = 'gsla-cpu'
    tree_lsf_cores: int = 1
    tree_lsf_mem: str = '4G'
    tree_wait_for_completion: bool = False      # fire-and-forget by default
    tree_conda_env: Optional[str] = None        # conda env for job script

    # Runtime fields (populated during processing, not from config file)
    ann_intbc_mod: Optional[pl.DataFrame] = None  # intBC → mod mapping from classify_cassettes
    

class CassiopeiaStep(Enum):
    """Steps within the Cassiopeia lineage extension"""
    PARSE_CONTIGS = 'parse_contigs'
    CLASSIFY_CASSETTES = "classify_cassettes"
    REGENERATE_FILTERED_QC = "regenerate_filtered_qc"  # Regenerate QC with cassette-type filtering
    PLUG_CASSIOPEIA = "plug_cassiopeia"
    EXTRACT_BARCODES = "extract_barcodes"
    GENERATE_MATRIX = "generate_matrix"
    GENERATE_METADATA = "generate_metadata"
    SEGMENTED_ALLELE = "segmented_allele"  # Direct allele table from segments
    BUILD_TREES = "build_trees"

class CassiopeiaLineageExtension(PostProcessorExtension):
    xp: FractureXp
    config: ExtensionConfig
    def __init__(self, xp: FractureXp):
        super().__init__(xp)
        self.temp_data = {}
    
    def get_config_class(self) -> Type[ExtensionConfig]:
        return CassiopeiaConfig

    @property
    def required_params(self) -> Set[str]:
        return self.config.get_required_fields()
    
    @property
    def name(self) -> str:
        return "cassiopeia_petracer"
    
    def process(self, contigs_path: Path) -> StepResults:
        """Main entry point - orchestrates all sub-steps"""
        
        # Initialize
        #self.temp_data['contigs'] = read_file(contigs_path)
        config = self.config
        force = getattr(config, 'force', False)  
        
        # Run steps based on configuration
        self.workdir = contigs_path.parent
        parsed_path = contigs_path.with_stem(contigs_path.stem + '_parsed')
        cass_mols_path = contigs_path.with_stem(contigs_path.stem + '_cass_mols')
        cass_allele_path = contigs_path.with_stem(contigs_path.stem + '_cass_allele')
        refs_path = contigs_path.parent / "refs.parquet"
        ann_intbc_mod_path = contigs_path.parent / "ann_intbc_mod_path.parquet"

        self.ldf = scan_file(contigs_path).filter(pl.col('contig').str.len_chars() > 0)

        if self.config.refs_fasta_path is not None:
            self.refs = generate_refs_from_fasta(
                    **self.config.get_function_config(generate_refs_from_fasta)
            )
            self.refs.write_parquet(refs_path)
            save_ref_to_fasta(self.refs, out_dir=self.workdir, field='mod')

        #self.temp_data['xp'] = xp
        #self.temp_data['outputs_dir'] = contigs_path.parent / "cassiopeia_outputs"
        #self.temp_data['outputs_dir'].mkdir(exist_ok=True)
        #xp.logger.info(f"{self.temp_data['outputs_dir']=}")

        final_results = {}
        final_metrics = {}
        

        if not parsed_path.exists() or force:
            self.xp.logger.info("Running parse_contigs step")
            result = self._parse_contigs()
            self.ldf = result.results['ldf']
            self.ldf.sink_parquet(parsed_path)
            self.xp.logger.info(f"Writing parsed contigs to {parsed_path}")
            final_metrics.update(result.metrics)
        else:
            self.xp.logger.info(f"Loading parsed contigs from {parsed_path}")
            self.ldf = scan_file(parsed_path)
        
        if self.should_run_step(CassiopeiaStep.CLASSIFY_CASSETTES.value):
            self.xp.logger.info("Running classify_cassettes step")

            #if annotation_path is None:
            if True:
                result = self._kmer_classify_cassettes()
            else:
                pass
                #result = use a given intbc_mod_map
            self.xp.logger.io(f"Saving intBC mod annotations to {ann_intbc_mod_path}")
            # result.results['ann_intbc_mod'] # maps intBC to mod
            self.ann_intbc_mod = result.results['ann_intbc_mod']
            self.ann_intbc_mod.write_parquet(ann_intbc_mod_path)

            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.REGENERATE_FILTERED_QC.value):
            self.xp.logger.info("Running regenerate_filtered_qc step")
            result = self._regenerate_filtered_qc()
            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.PLUG_CASSIOPEIA.value):
            self.xp.logger.info("Running plug_cassiopeia step")
            result = self._plug_cassiopeia()

            # Save raw per-UMI alleles (collapse happens downstream in _build_trees)
            self.alleles_pl = result.results['alleles_pl']
            self.alleles_pl.write_parquet(f"{self.workdir}/alleles_pl.parquet")
            self.xp.logger.io(f"Saved raw per-UMI alleles to {self.workdir}/alleles_pl.parquet")

            self.xp.logger.io(f"Saving cassiopeia allele tables to {cass_allele_path}")
            self.ldf.sink_parquet(cass_allele_path)

        if self.should_run_step(CassiopeiaStep.SEGMENTED_ALLELE.value):
            self.xp.logger.info("Running segmented_allele step")
            result = self._segmented_allele()
            self.alleles_segmented = result.results['allele_table_segmented']
            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.BUILD_TREES.value):
            self.xp.logger.info("Running build_trees step")
            result = self._build_trees()
            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.GENERATE_MATRIX.value):
            self.xp.logger.info("Running generate_matrix step")
            result = self._generate_matrix()
            final_results.update(result.results)
            final_metrics.update(result.metrics)
        
        if self.should_run_step(CassiopeiaStep.GENERATE_METADATA.value):
            self.xp.logger.info("Running generate_metadata step")
            result = self._generate_metadata()
            final_results.update(result.results)
            final_metrics.update(result.metrics)
        
        return StepResults(results=final_results, metrics=final_metrics)

    def _template(self) -> StepResults:
        """ """
        return StepResults(
            results={"":[]},
            metrics={"":[]},
        )

    def _parse_contigs(self) -> StepResults:
        return parse_contigs(ldf=self.ldf, 
                             **self.config.get_function_config(parse_contigs))
    
    def _kmer_classify_cassettes(self) -> StepResults:
        """ Guesses what reference a given intBC corresponds to based on kmer composition analysis"""
        return kmer_classify_cassettes(
                ldf=self.ldf,
                refs=self.refs,
                **self.config.get_function_config(kmer_classify_cassettes)
                )

    def _regenerate_filtered_qc(self) -> StepResults:
        """
        Regenerate segmentation QC plots with cassette-type filtering.

        After cassette classification, we know the cassette type (5mer, 10mer, etc.)
        and can filter out spurious meta matches in the QC plots.

        This step:
        1. Loads segmentation data from intermediate files
        2. Determines cassette_type from refs (e.g., mod_5mer -> "5mer")
        3. Regenerates QC plots with filtered segment_types
        """
        from ..pipeline.segmentation import (
            generate_segmentation_report,
            plot_segmentation_qc,
        )

        # Determine cassette_type from refs
        # refs has 'mod' column like "mod_5mer", "mod_10mer"
        if not hasattr(self, 'refs') or self.refs is None:
            self.xp.logger.warning("No refs available, skipping filtered QC regeneration")
            return StepResults(results={}, metrics={})

        # Get the dominant cassette type from refs (or use config)
        mods = self.refs['mod'].unique().to_list()
        # Extract cassette type: "mod_5mer" -> "5mer"
        cassette_types = [m.replace("mod_", "") for m in mods if m.startswith("mod_")]

        if not cassette_types:
            self.xp.logger.warning("Could not determine cassette type from refs")
            return StepResults(results={}, metrics={})

        # Use the first (or only) cassette type
        cassette_type = cassette_types[0]
        self.xp.logger.info(f"Regenerating QC with cassette_type={cassette_type} filtering")

        # Load intermediate files
        intermediate_dir = self.workdir / "intermediate"
        segments_path = intermediate_dir / "segments.parquet"
        assembled_path = intermediate_dir / "assembled.parquet"

        # Also try legacy path (with _debug suffix in intermediate dir)
        if not segments_path.exists():
            segments_path = intermediate_dir / "segments_debug.parquet"

        if not segments_path.exists():
            self.xp.logger.warning(f"Segments file not found at {segments_path}, skipping filtered QC")
            return StepResults(results={}, metrics={})

        # Load metas
        metas_csv = self.xp.fracture.get('metas_csv') or getattr(self.xp, 'features_csv', None)
        if not metas_csv:
            self.xp.logger.warning("No metas_csv configured, skipping filtered QC")
            return StepResults(results={}, metrics={})

        metas = pl.read_csv(metas_csv)
        segments_df = read_file(segments_path)

        # Load assembled if available. assemble_segmented writes with the
        # suffix "assembled_debug.parquet" (derived from the segments debug
        # path), so check the legacy plain name and the _debug variant.
        if not assembled_path.exists():
            assembled_path = intermediate_dir / "assembled_debug.parquet"

        if assembled_path.exists():
            assembled_df = read_file(assembled_path)
        else:
            # Create minimal assembled_df from segments, keeping whichever
            # molecule-identifying columns segments_df has (sbc + umi, or just
            # umi). This keeps the schema consistent with segments_df so that
            # downstream per-molecule joins in generate_segmentation_report
            # (e.g. the missing-segments rollup) find the columns they expect.
            mol_cols = ['sbc', 'umi'] if 'sbc' in segments_df.columns else ['umi']
            assembled_df = segments_df.select(mol_cols + ['start_meta', 'end_meta']).unique()

        # Load contigs (check for both IPC and Parquet formats)
        # Use specific pattern to avoid matching intermediate files like _parsed or _cass_allele
        contigs_path = list(self.workdir.glob("contigs_segmented_valid.arrow")) or \
                       list(self.workdir.glob("contigs_segmented_valid.parquet")) or \
                       list(self.workdir.glob("contigs_*.arrow")) or \
                       list(self.workdir.glob("contigs_*.parquet"))
        if contigs_path:
            contigs_df = read_file(contigs_path[0])
            if 'stitched_seq' in contigs_df.columns and 'contig' not in contigs_df.columns:
                contigs_df = contigs_df.rename({'stitched_seq': 'contig'})
        else:
            # Minimal contigs from ldf
            contigs_df = self.ldf.select(['umi', 'contig']).collect()

        # Get anchors from config
        cassette_start_anchor = getattr(self.xp, 'start_anchor', None)
        cassette_end_anchor = getattr(self.xp, 'end_anchor', None)

        # Get filtering threshold from config (default 70%)
        min_meta_freq = self.xp.fracture.get('min_meta_freq_threshold', 0.7)

        # Get heterogeneity threshold from config (default 0.20 = 20% of dominant)
        heterogeneity_threshold = self.xp.fracture.get('heterogeneity_threshold', 0.20)

        # Generate unfiltered report first
        report = generate_segmentation_report(
            segments_df=segments_df,
            assembled_df=assembled_df,
            contigs_df=contigs_df,
            metas=metas,
            cassette_start_anchor=cassette_start_anchor,
            cassette_end_anchor=cassette_end_anchor,
            heterogeneity_threshold=heterogeneity_threshold,
        )

        # Now filter the segment_types in the report
        segment_types_df = pl.DataFrame(report['segments']['segment_types'])
        unique_umis = report['segments']['unique_umis']

        filtered_segment_types = filter_segment_types_for_cassette(
            segment_types=segment_types_df,
            metas=metas,
            cassette_type=cassette_type,
            min_freq_threshold=min_meta_freq,
            total_molecules=unique_umis,
            logger=self.xp.logger,
        )

        # Update the report with filtered segment_types
        report['segments']['segment_types'] = filtered_segment_types.to_dicts()

        # Generate filtered QC plots - use central figures directory per Xp
        figures_dir = Path(self.xp.sample_figs)
        figures_dir.mkdir(parents=True, exist_ok=True)

        sample_name = getattr(self.xp, 'target_sample', 'sample')
        plot_segmentation_qc(
            report=report,
            contigs_df=contigs_df,
            metas=metas,
            output_dir=figures_dir,
            sample_name=f"{sample_name}_filtered",
            contig_col='contig',
            logger=self.xp.logger,
        )

        self.xp.logger.info(f"Saved filtered QC plots to {figures_dir}")

        return StepResults(
            results={"filtered_qc_dir": str(figures_dir)},
            metrics={"cassette_type": cassette_type},
        )

    def _convert_to_allele_table(self) -> StepResults:
        """ """
        return StepResults(
            results={"":[]},
            metrics={"":[]},
        )
    
    def _plug_cassiopeia(self) -> StepResults:
        """
        readName - A unique identifier for each row/sequence
        cellBC - The cell barcode
        UMI - The UMI (Unique Molecular Identifier)
        readCount - The number of reads for this sequence
        seq - The actual sequence to be aligned

        """
        # Update config with runtime parameter
        self.config.ann_intbc_mod = self.ann_intbc_mod
        config = self.config.get_function_config(plug_cassiopeia)

        # Get modality from experiment config
        modality = getattr(self.xp, 'modality', 'single-molecule')
        cbc_len = getattr(self.xp, 'cbc_len', 16)

        return plug_cassiopeia(ldf=self.ldf,
                               workdir=self.workdir,
                               logger=self.xp.logger,
                               modality=modality,
                               cbc_len=cbc_len,
                               **config)

    def _segmented_allele(self) -> StepResults:
        """
        Generate allele table directly from segments (no assembly needed).

        Reads segments_debug.parquet from workdir/intermediate and generates allele_table_segmented.parquet.
        """
        intermediate_dir = self.workdir / "intermediate"
        segments_path = intermediate_dir / "segments_debug.parquet"

        # Also try new naming convention
        if not segments_path.exists():
            segments_path = intermediate_dir / "segments.parquet"

        if not segments_path.exists():
            raise FileNotFoundError(
                f"Segments file not found in {intermediate_dir}. "
                "Run fracture with use_segmentation=True first."
            )

        segments_df = read_file(segments_path)

        if self.config.metas_flanks_csv is None:
            raise ValueError("metas_flanks_csv must be set in config for segmented_allele step")

        metas_df = pl.read_csv(self.config.metas_flanks_csv)

        allele_table = segments_to_allele_table(
            segments_df=segments_df,
            metas_df=metas_df,
            int_anchor1=self.config.int_anchor1,
            int_anchor2=self.config.int_anchor2,
            cassette_type=self.config.cassette_type,
            min_consensus_support=self.config.min_consensus_support,
            logger=self.xp.logger,
        )

        output_path = self.workdir / "allele_table_segmented.parquet"
        allele_table.write_parquet(output_path)
        self.xp.logger.io(f"Saved segmented allele table to {output_path}")

        return StepResults(
            results={'allele_table_segmented': allele_table},
            metrics={
                'total_molecules': allele_table.height,
                'unique_intbcs': allele_table['intBC'].n_unique() if allele_table.height > 0 else 0,
            }
        )

    def _get_solver(self, solver_name: str):
        """Instantiate a Cassiopeia solver by name."""
        import cassiopeia as cas
        match solver_name:
            case 'nj':
                return cas.solver.NeighborJoiningSolver(
                    dissimilarity_function=cas.solver.dissimilarity.weighted_hamming_distance,
                    add_root=True,
                )
            case 'vanilla':
                return cas.solver.VanillaGreedySolver()
            case 'mcgs':
                return cas.solver.MaxCutGreedySolver()
            case _:
                raise ValueError(f"Unknown solver: {solver_name}")

    @staticmethod
    def _sanitize_states(states):
        """Flatten tuple/list character states to plain ints for convexml compatibility.

        Cassiopeia can produce tuple states like (1, 2) for ambiguous sites.
        Convexml needs plain integers in a numpy float64 array.
        Strategy: take first element of tuples; leave scalars as-is.
        """
        return [s[0] if isinstance(s, (tuple, list)) else s for s in states]

    def _estimate_branch_lengths(self, cas_tree, config):
        """Estimate branch lengths via convexml."""
        import convexml
        leaf_sequences = {
            l: self._sanitize_states(cas_tree.get_character_states(l))
            for l in cas_tree.leaves
        }
        tree_newick = convexml.convexml(
            tree_newick=cas_tree.get_newick(record_branch_lengths=True),
            leaf_sequences=leaf_sequences,
            minimum_branch_length=config.min_branch_length,
        )['tree_newick']
        import cassiopeia as cas
        return cas.data.CassiopeiaTree(
            character_matrix=cas_tree.character_matrix,
            tree=tree_newick,
            missing_state_indicator=cas_tree.missing_state_indicator,
        )

    @staticmethod
    def _build_allele_palette(obs_df, ann_cols):
        """Build allele -> hex color palette from obs DataFrame annotation columns.

        Normalizes values to strings, assigns ColorHash colors,
        then overrides unedited ('None') to white and missing (NaN) to light gray.
        """
        from colorhash import ColorHash
        import pandas as pd

        palette = {}
        for col in ann_cols:
            for v in obs_df[col].dropna().unique():
                s = str(v)
                if s in palette or s == 'nan':
                    continue
                palette[s] = '#FFFFFF' if 'None' in s else ColorHash(s).hex
        # Ensure missing/nan entries are gray
        palette['nan'] = '#E0E0E0'

        # Stringify obs values so they match palette keys
        for col in ann_cols:
            obs_df[col] = obs_df[col].astype(str)

        return palette

    @staticmethod
    def _save_legend(palette, ann_cols, obs_df, outpath, ncol=1):
        """Save a standalone legend plot mapping allele values to colors.

        Groups legend entries by annotation column so each cutsite/intBC
        has its own section.
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        # Collect per-column unique values (preserving order)
        col_values = {}
        for col in ann_cols:
            vals = [str(v) for v in obs_df[col].dropna().unique() if str(v) != 'nan']
            col_values[col] = sorted(set(vals))

        handles = []
        labels = []
        for col, vals in col_values.items():
            # Section header
            handles.append(Patch(facecolor='none', edgecolor='none'))
            labels.append(f'— {col} —')
            for v in vals:
                color = palette.get(v, '#CCCCCC')
                handles.append(Patch(facecolor=color, edgecolor='#666666', linewidth=0.5))
                labels.append(v)

        n_entries = len(handles)
        fig_height = max(2, n_entries * 0.18)
        fig, ax = plt.subplots(figsize=(3 * ncol, fig_height))
        ax.axis('off')
        ax.legend(handles, labels, loc='center', ncol=ncol,
                  fontsize=5, frameon=False, handlelength=1.2, handleheight=0.8)
        fig.savefig(outpath, bbox_inches='tight', dpi=300)
        plt.close(fig)

    def _build_tree_sc(
        self,
        allele_table: pl.DataFrame,
        rcols: list[str],
        outdir: Path,
    ) -> StepResults:
        """Build a single tree from the full allele table (single-cell mode).

        One tree per sample — all intBCs become characters in the matrix.
        Follows the tested logic from build_tree_sc.py.
        """
        import cassiopeia as cas
        import treedata as td
        import pycea as pc
        import hashlib
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from colorhash import ColorHash
        import re as _re

        config = self.config
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        n_cells = allele_table['cellBC'].n_unique()
        if n_cells < config.min_cells_for_tree:
            self.xp.logger.warning(f"SC tree: only {n_cells} cells, need {config.min_cells_for_tree}. Skipping.")
            write_done(outdir, f'too_few_entries_{n_cells}')
            return StepResults(results={}, metrics={'skipped': True, 'reason': 'too_few_entries'})

        # Fill nulls for cassiopeia (spanning deletions already handled by _build_trees orchestrator)
        allele_table = allele_table.with_columns([pl.col(c).fill_null('None') for c in rcols])
        pd_allele = allele_table.to_pandas()

        # Character matrix
        character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
            pd_allele, allele_rep_thresh=config.allele_rep_thresh,
        )
        self.xp.logger.info(f"SC character matrix: {character_matrix.shape[0]} cells x {character_matrix.shape[1]} characters")
        # Cast to str to avoid pyarrow ArrowInvalid with mixed list/non-list values in object columns
        cm_out = character_matrix.reset_index()
        cm_out = cm_out.astype({c: str for c in cm_out.columns if cm_out[c].dtype == object})
        cm_out.to_parquet(outdir / 'character_matrix.parquet')

        cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)

        if cas_tree.n_cell < config.min_cells_for_tree:
            write_done(outdir, f'too_few_entries_after_filter_{cas_tree.n_cell}')
            return StepResults(results={}, metrics={'skipped': True, 'reason': 'too_few_after_filter'})

        # Solve
        self.xp.logger.info(f"SC: solving tree ({config.solver}, {cas_tree.n_cell} cells)...")
        solver = self._get_solver(config.solver)
        solver.solve(cas_tree, collapse_mutationless_edges=True)

        # Branch length estimation
        if not config.skip_branch_lengths:
            self.xp.logger.info("SC: estimating branch lengths (convexml)...")
            cas_tree = self._estimate_branch_lengths(cas_tree, config)
        else:
            self.xp.logger.info("SC: skipping branch length estimation")

        # Save newick
        newick = cas_tree.get_newick(record_branch_lengths=True)
        (outdir / 'newick.txt').write_text(newick)

        # Build wide allele obs for TreeData
        allele_wide = pivot_alleles_wide(allele_table, rcols)
        allele_wide_pd = allele_wide.to_pandas().set_index('cellBC')
        leaves_in_obs = [l for l in cas_tree.leaves if l in allele_wide_pd.index]
        obs = allele_wide_pd.loc[leaves_in_obs]

        tdata = td.TreeData(
            X=None, allow_overlap=True, obs=obs.copy(),
            obst={config.solver: cas_tree.get_tree_topology()},
        )
        add_weighted_depth(tdata, tree_key=config.solver)
        tdata.write_h5td(outdir / 'tdata.h5td')

        md5 = hashlib.md5(str(tdata).encode()).hexdigest()[:8]
        (outdir / 'tdata_md5.txt').write_text(md5)
        self.xp.logger.info(f"SC: saved TreeData ({cas_tree.n_cell} cells, md5={md5}) -> {outdir}")

        # Annotation columns and colors
        ann_cols = [c for c in tdata.obs.columns if _re.search(r'_r\d+$', c)]
        allele_colors_hex = self._build_allele_palette(tdata.obs, ann_cols)

        sample_name = getattr(self.xp, 'target_sample', 'sample')
        n_cells_tree = cas_tree.n_cell
        n_ann = len(ann_cols)
        dendrogram_ratio = 0.05
        annotation_width = (1 - dendrogram_ratio) / (dendrogram_ratio * max(n_ann, 1))
        base_name = f'sc_{n_cells_tree}_{sample_name}_{config.solver}_span-{config.spanning_deletions}'

        # Circular
        self.xp.logger.info(f"SC: plotting circular tree ({n_cells_tree} cells, {n_ann} annotations)...")
        fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=900, subplot_kw={'projection': 'polar'})
        pc.pl.tree(tdata, tree=config.solver, keys=ann_cols, polar=True,
                   depth_key='depth', extend_branches=False, branch_linewidth=0.15,
                   annotation_width=0.01, palette=allele_colors_hex, ax=ax, legend=False)
        fig.savefig(outdir / f'{base_name}_circ.png', bbox_inches='tight')
        plt.close(fig)

        # Linear
        height_inches = max(5, (n_cells_tree / 1000) * 750 / 100)
        self.xp.logger.info(f"SC: plotting linear tree ({height_inches:.0f}\" tall)...")
        fig, ax = plt.subplots(1, 1, figsize=(5, height_inches), dpi=900)
        pc.pl.tree(tdata, tree=config.solver, keys=ann_cols, polar=False,
                   depth_key='depth', extend_branches=False, branch_linewidth=0.15,
                   annotation_width=annotation_width, palette=allele_colors_hex, ax=ax, legend=False)
        fig.savefig(outdir / f'{base_name}_linear.png', bbox_inches='tight')
        plt.close(fig)

        # Standalone legend
        self._save_legend(allele_colors_hex, ann_cols, tdata.obs, outdir / f'{base_name}_legend.png')

        write_done(outdir)
        self.xp.logger.info(f"SC: done — {n_cells_tree} cells, {cas_tree.n_character} characters, output: {outdir}")
        return StepResults(
            results={'tdata_path': str(outdir / 'tdata.h5td')},
            metrics={'n_cells': n_cells_tree, 'n_characters': cas_tree.n_character},
        )

    def _build_tree_sm(
        self,
        allele_table: pl.DataFrame,
        rcols: list[str],
        outdir: Path,
        filter_dict: dict,
    ) -> StepResults:
        """Build a tree for a single group (single-molecule mode).

        Args:
            filter_dict: Column->value mapping to filter the allele table
                         e.g. {'intBC': 'TCTGAA'} or {'intBC': 'TCTGAA', 'sbc': 'TCAAGT'}
        """
        import cassiopeia as cas
        import treedata as td
        import pycea as pc
        import hashlib
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from colorhash import ColorHash

        config = self.config
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        group_label = '_'.join(str(v) for v in filter_dict.values())

        # Filter to this group
        sub = allele_table
        for col, val in filter_dict.items():
            sub = sub.filter(pl.col(col) == val)

        # Per-group subsampling (SM test mode)
        if config.tree_test_mode and config.tree_subsample:
            sub = self._subsample(sub, config.tree_subsample, self.xp.logger)

        if sub.height < config.min_cells_for_tree:
            self.xp.logger.info(f"SM tree {group_label}: only {sub.height} molecules (need {config.min_cells_for_tree}), skipping")
            write_done(outdir, f'too_few_entries_{sub.height}')
            return StepResults(results={}, metrics={'skipped': True, 'reason': 'too_few_molecules'})

        # Fill nulls for cassiopeia (spanning deletions already handled by _build_trees orchestrator)
        sub = sub.with_columns([pl.col(c).fill_null('None') for c in rcols])
        pd_allele = sub.to_pandas()

        # Character matrix
        self.xp.logger.info(f"SM tree {group_label}: {sub.height} molecules, building character matrix...")
        character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
            pd_allele, allele_rep_thresh=config.allele_rep_thresh,
        )
        self.xp.logger.info(f"SM tree {group_label}: character matrix {character_matrix.shape[0]} x {character_matrix.shape[1]}")
        # Cast to str to avoid pyarrow ArrowInvalid with mixed list/non-list values in object columns
        cm_out = character_matrix.reset_index()
        cm_out = cm_out.astype({c: str for c in cm_out.columns if cm_out[c].dtype == object})
        cm_out.to_parquet(outdir / 'character_matrix.parquet')

        cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix)

        if cas_tree.n_cell < config.min_cells_for_tree:
            self.xp.logger.info(f"SM tree {group_label}: only {cas_tree.n_cell} cells after filtering, skipping")
            write_done(outdir, f'too_few_entries_after_filter_{cas_tree.n_cell}')
            return StepResults(results={}, metrics={'skipped': True, 'reason': 'too_few_after_filter'})

        # Solve
        self.xp.logger.info(f"SM tree {group_label}: solving ({config.solver}, {cas_tree.n_cell} molecules)...")
        solver = self._get_solver(config.solver)
        solver.solve(cas_tree, collapse_mutationless_edges=True)

        # Branch lengths
        if not config.skip_branch_lengths:
            self.xp.logger.info(f"SM tree {group_label}: estimating branch lengths (convexml)...")
            cas_tree = self._estimate_branch_lengths(cas_tree, config)

        # Save newick
        newick = cas_tree.get_newick(record_branch_lengths=True)
        (outdir / 'newick.txt').write_text(newick)

        # Build obs from r columns
        cell_id_col = 'UMI' if 'UMI' in sub.columns else 'cellBC'
        allele_matrix = sub.select(cell_id_col, *rcols).to_pandas().set_index(cell_id_col)

        tdata = td.TreeData(
            X=None, allow_overlap=True,
            obs=allele_matrix.loc[cas_tree.leaves].copy(),
            obst={config.solver: cas_tree.get_tree_topology()},
        )
        add_weighted_depth(tdata, tree_key=config.solver)
        tdata.write_h5td(outdir / 'tdata.h5td')

        md5 = hashlib.md5(str(tdata).encode()).hexdigest()[:8]
        (outdir / 'tdata_md5.txt').write_text(md5)
        self.xp.logger.info(f"SM tree {group_label}: saved TreeData ({cas_tree.n_cell} molecules, md5={md5})")

        # Plot
        allele_colors_hex = self._build_allele_palette(tdata.obs, rcols)

        n_mols = cas_tree.n_cell
        dendrogram_ratio = 0.05
        annotation_width = (1 - dendrogram_ratio) / (dendrogram_ratio * len(rcols))
        base_name = f'sm_{n_mols}_{group_label}_{config.solver}'

        # Circular
        self.xp.logger.info(f"SM tree {group_label}: plotting ({n_mols} molecules)...")
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=900, subplot_kw={'projection': 'polar'})
        pc.pl.tree(tdata, tree=config.solver, keys=rcols, polar=True,
                   depth_key='depth', extend_branches=True, branch_linewidth=0.05,
                   annotation_width=0.02, palette=allele_colors_hex, ax=ax, legend=False)
        fig.savefig(outdir / f'{base_name}_circ.png', bbox_inches='tight')
        plt.close(fig)

        # Linear
        height_inches = max(5, (n_mols / 1000) * 750 / 100)
        fig, ax = plt.subplots(1, 1, figsize=(2, height_inches), dpi=900)
        pc.pl.tree(tdata, tree=config.solver, keys=rcols, polar=False,
                   depth_key='depth', extend_branches=True, branch_linewidth=0.05,
                   annotation_width=annotation_width, palette=allele_colors_hex, ax=ax, legend=False)
        fig.savefig(outdir / f'{base_name}_linear.png', bbox_inches='tight')
        plt.close(fig)

        # Standalone legend
        self._save_legend(allele_colors_hex, rcols, tdata.obs, outdir / f'{base_name}_legend.png')

        write_done(outdir)
        self.xp.logger.info(f"SM tree {group_label}: done — {n_mols} molecules, {cas_tree.n_character} characters, output: {outdir}")
        return StepResults(
            results={'tdata_path': str(outdir / 'tdata.h5td')},
            metrics={'n_molecules': n_mols, 'n_characters': cas_tree.n_character},
        )
                
    def _subsample(self, allele_table: pl.DataFrame, factors: list, logger=None) -> pl.DataFrame:
        """Apply subsampling factors to an allele table.

        Each factor is a dict with keys:
          column: str           — column to filter/subsample on
          values: list          — explicit whitelist (skips top_x/sample_n)
          top_x: int            — keep top X groups ranked by size
          sample_n: int         — then sample N from that pool (optional)
        Factors are applied in order.
        """
        for factor in factors:
            col = factor['column']
            if col not in allele_table.columns:
                if logger:
                    logger.info(f"  {col}: column not found, skipping")
                continue

            n_unique = allele_table.select(col).n_unique()
            values = factor.get('values')

            if values is not None:
                # Explicit whitelist
                pool = pl.DataFrame({col: values})
                allele_table = allele_table.join(pool, on=col, how='semi')
                if logger:
                    matched = allele_table.select(col).n_unique()
                    logger.info(f"  {col}: {n_unique} unique | values: {len(values)} requested, {matched} matched | -> {allele_table.height} rows")
            else:
                # Statistical subsampling
                top_x = factor.get('top_x')
                sample_n = factor.get('sample_n')
                pool = (
                    allele_table.group_by(col).len()
                    .sort('len', descending=True)
                    .head(top_x)
                    .select(col)
                )
                n_after_top = pool.height
                if sample_n is not None:
                    actual_n = min(sample_n, pool.height)
                    pool = pool.sample(n=actual_n, seed=42)
                allele_table = allele_table.join(pool, on=col, how='semi')
                if logger:
                    parts = [f"{col}: {n_unique} unique"]
                    parts.append(f"top {top_x} -> {n_after_top}")
                    if sample_n is not None:
                        parts.append(f"sampled {actual_n}/{sample_n}")
                    parts.append(f"-> {allele_table.height} rows")
                    logger.info(f"  {' | '.join(parts)}")
        return allele_table

    def _build_trees(self) -> StepResults:
        """BUILD_TREES step: prepare tree jobs and optionally submit via LSF."""
        import re

        config = self.config
        workdir = self.workdir

        # Always load raw per-UMI alleles — filters and collapse happen here
        allele_path = workdir / 'alleles_pl.parquet'
        allele_table_raw = pl.read_parquet(allele_path)
        rcols = [c for c in allele_table_raw.columns if re.match(r'^r\d+$', c)]

        # Apply min_molecule_len filter on raw per-UMI data (before collapse)
        if config.min_molecule_len is not None and 'Seq' in allele_table_raw.columns:
            n_before = allele_table_raw.height
            allele_table_raw = allele_table_raw.filter(pl.col('Seq').str.len_chars() >= config.min_molecule_len)
            self.xp.logger.info(
                f"min_molecule_len filter ({config.min_molecule_len}): {n_before} -> {allele_table_raw.height} rows"
            )

        # Determine tree mode
        modality = getattr(self.xp, 'modality', 'single-molecule')
        mode = config.tree_mode
        if mode == 'auto':
            mode = 'sc' if modality == 'single-cell' else 'sm'

        # Resolve spanning deletion modes to iterate over
        span_modes = ['unedited', 'use'] if config.spanning_deletions == 'both' else [config.spanning_deletions]

        tree_base = workdir / 'trees'
        if config.tree_test_mode:
            tree_base = tree_base / 'test'

        all_results = {}
        all_metrics: Dict[str, int] = {'trees_built': 0, 'trees_cached': 0, 'trees_submitted': 0}

        for span_mode in span_modes:
            self.xp.logger.info(
                f"Tree mode: {mode} | solver: {config.solver} | "
                f"spanning: {span_mode} | "
                f"branch_lengths: {'skip' if config.skip_branch_lengths else 'convexml'} | "
                f"execution: {'lsf' if config.tree_use_lsf else 'inline'}"
            )

            # Handle spanning deletions on raw data (before collapse so mode() isn't biased)
            allele_table = flag_spanning_deletions(allele_table_raw, rcols, mode=span_mode)

            # Collapse UMIs to cells for SC mode (on filtered data)
            if mode == 'sc' and allele_table.height > 0:
                n_before = allele_table.height
                allele_table = collapse_umis_to_cells(
                    allele_table,
                    cell_col='cellBC',
                    intbc_col='intBC',
                    min_umis_per_cell=config.min_umis_per_cell,
                    min_umi_agreement=config.min_umi_agreement,
                    method=config.collapse_method,
                    logger=self.xp.logger,
                )
                self.xp.logger.info(
                    f"Collapsed UMIs to cells: {n_before} -> {allele_table.height} rows"
                )
                # Save collapsed table for postmortem (tagged by span mode)
                suffix = f'_span_{span_mode}' if len(span_modes) > 1 else ''
                collapsed_path = workdir / f'alleles_pl_collapsed{suffix}.parquet'
                allele_table.write_parquet(collapsed_path)
                self.xp.logger.io(f"Saved collapsed alleles to {collapsed_path}")

            n_cells = allele_table['cellBC'].n_unique() if 'cellBC' in allele_table.columns else allele_table.height
            n_intbcs = allele_table['intBC'].n_unique() if 'intBC' in allele_table.columns else 0
            self.xp.logger.info(
                f"Loaded allele table: {allele_table.height} rows, {n_cells} cells, "
                f"{n_intbcs} intBCs, {len(rcols)} cutsites — from {allele_path.name}"
            )

            tree_outdir = tree_base / f'span_{span_mode}'
            tree_outdir.mkdir(parents=True, exist_ok=True)
            self.xp.logger.info(f"Output dir: {tree_outdir}")

            # Subsample for testing (SC: global, SM: per-group in _build_tree_sm)
            if config.tree_test_mode and config.tree_subsample:
                self.xp.logger.info("Tree test mode enabled — subsampling data"
                                    + (" (globally for SC)" if mode == 'sc' else " (per-group for SM)"))

            if mode == 'sc':
                allele_for_trees = allele_table
                if config.tree_test_mode and config.tree_subsample:
                    allele_for_trees = self._subsample(allele_table, config.tree_subsample, self.xp.logger)
                jobs = [{'mode': 'sc', 'outdir': tree_outdir}]
            else:
                allele_for_trees = allele_table
                group_cols = [c for c in config.tree_group_by if c in allele_table.columns]
                if not group_cols:
                    group_cols = ['intBC']
                groups = allele_table.select(group_cols).unique()
                self.xp.logger.info(f"SM grouping by {group_cols}: {groups.height} groups")
                jobs = []
                for row in groups.iter_rows(named=True):
                    label = '_'.join(str(row[c]) for c in group_cols)
                    job_outdir = tree_outdir / label
                    job = {'mode': 'sm', 'outdir': job_outdir, 'filter': {c: row[c] for c in group_cols}}
                    jobs.append(job)

            # Filter out already-done jobs
            force = getattr(config, 'force', False)
            pending = [j for j in jobs if not (Path(j['outdir']) / 'done').exists() or force]

            if not pending:
                self.xp.logger.info(f"All {len(jobs)} trees already built (cached) for span_{span_mode}")
                all_metrics['trees_cached'] += len(jobs)
                continue

            self.xp.logger.info(f"Tree building: {len(pending)} pending, {len(jobs) - len(pending)} cached")

            if config.tree_use_lsf:
                job_ids = self._submit_tree_jobs_lsf(pending, allele_path, rcols)
                all_results[f'lsf_job_ids_{span_mode}'] = job_ids
                all_metrics['trees_submitted'] += len(pending)
                all_metrics['trees_cached'] += len(jobs) - len(pending)
            else:
                # Inline execution
                for job in pending:
                    if job['mode'] == 'sc':
                        self._build_tree_sc(allele_for_trees, rcols, job['outdir'])
                    else:
                        self._build_tree_sm(allele_for_trees, rcols, job['outdir'], job['filter'])
                all_metrics['trees_built'] += len(pending)
                all_metrics['trees_cached'] += len(jobs) - len(pending)

            all_results[f'tree_outdir_{span_mode}'] = str(tree_outdir)

        return StepResults(
            results=all_results,
            metrics={k: v for k, v in all_metrics.items() if v > 0},
            )

    def _make_tree_cmd(self, job: dict, allele_path: Path, rcols: list[str]) -> str:
        """Build the CLI command for a single tree-building LSF job."""
        config = self.config
        n_rcols = len(rcols)
        parts = [
            'python', '-m', 'ogtk.ltr.fracture.extensions.cassiopeia_petracer',
            'build-tree',
            '--input', str(allele_path),
            '--outdir', str(job['outdir']),
            '--mode', job['mode'],
            '--solver', config.solver,
            '--n-rcols', str(n_rcols),
            '--spanning-deletions', config.spanning_deletions,
            '--allele-rep-thresh', str(config.allele_rep_thresh),
            '--min-cells', str(config.min_cells_for_tree),
            '--min-branch', str(config.min_branch_length),
        ]
        if config.skip_branch_lengths:
            parts.append('--skip-branch-lengths')
        if job['mode'] == 'sm' and 'filter' in job:
            import json
            parts.extend(['--filter', json.dumps(job['filter'])])
        if hasattr(self.xp, 'target_sample'):
            parts.extend(['--sample-name', self.xp.target_sample])
        if config.tree_test_mode and config.tree_subsample:
            import json
            parts.append('--test-mode')
            parts.extend(['--subsample', json.dumps(config.tree_subsample)])

        return ' '.join(parts)

    def _submit_tree_jobs_lsf(self, jobs: list, allele_path: Path, rcols: list[str]) -> list[str]:
        """Generate and submit LSF jobs for tree building."""
        import re
        import subprocess

        config = self.config

        # Ensure log directory exists
        log_dir = self.workdir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        job_ids = []
        for job in jobs:
            cmd = self._make_tree_cmd(job, allele_path, rcols)

            job_name = Path(job['outdir']).name if job['mode'] == 'sm' else 'sc_tree'

            # Build job script (follows tree_qc pattern: heredoc via stdin)
            conda_activate = ''
            if config.tree_conda_env:
                conda_activate = (
                    f'eval "$(/home/projects/nyosef/pedro/miniforge3/bin/conda shell.bash hook)"\n'
                    f'source ~/miniforge3/etc/profile.d/mamba.sh\n'
                    f'conda activate {config.tree_conda_env}\n'
                )

            job_script = f'#!/bin/bash\n{conda_activate}{cmd}\n'

            bsub_cmd = [
                'bsub',
                '-q', config.tree_lsf_queue,
                '-n', str(config.tree_lsf_cores),
                '-R', 'span[hosts=1]',
                '-R', f'rusage[mem={config.tree_lsf_mem}]',
                '-o', str(log_dir / f'tree_{job_name}.out'),
                '-e', str(log_dir / f'tree_{job_name}.err'),
                '-J', f'tree_{job_name}',
            ]

            result = subprocess.run(bsub_cmd, input=job_script, capture_output=True, text=True)
            if result.returncode == 0:
                match = re.search(r'Job <(\d+)>', result.stdout)
                if match:
                    job_ids.append(match.group(1))
                    self.xp.logger.info(f"Submitted tree job {job_name}: {match.group(1)}")
            else:
                self.xp.logger.error(f"Failed to submit tree job {job_name}: {result.stderr}")

        if config.tree_wait_for_completion and job_ids:
            from ...ltr.fracture.post.tree_qc import _wait_for_lsf_jobs
            self.xp.logger.info(f"Waiting for {len(job_ids)} tree jobs to complete...")
            _wait_for_lsf_jobs(job_ids)

        return job_ids

    def _extract_barcodes(self) -> StepResults:
        """Extract integration and sample barcodes"""
        xp = self.temp_data['xp']
        df_contigs = self.temp_data['contigs']
        
        # Your barcode extraction logic here
        df_annotated = (
            df_contigs
            .with_columns([
                pl.col('contig')
                .str.extract(f'({xp.intbc_5prime}[ATCG]{{10,20}})')
                .alias('integration_barcode'),
                
                pl.col('contig')
                .str.slice(0, xp.sbc_len)
                .alias('sample_barcode')
            ])
            .filter(pl.col('integration_barcode').is_not_null())
        )
        
        output_path = self.temp_data['outputs_dir'] / "barcodes_extracted.parquet"

        df_annotated.write_parquet(output_path)
        self.temp_data['annotated_contigs'] = df_annotated
        
        return StepResults(
            results={"barcodes_extracted": str(output_path)},
            metrics={
                "contigs_with_integration_bc": df_annotated.height,
                "unique_integration_bcs": df_annotated.select('integration_barcode').n_unique()
            }
        )
    
    
    def _generate_matrix(self) -> StepResults:
        # Implementation here  
        pass
        
    def _generate_metadata(self) -> StepResults:
        # Implementation here
        pass

# Register the extension
extension_registry.register(CassiopeiaLineageExtension)


def _build_tree_standalone(args):
    """Entry point for LSF-submitted tree jobs. Mirrors build_tree_sc.py / build_tree_single.py."""
    import re
    import hashlib
    import cassiopeia as cas
    import treedata as td
    import pycea as pc
    import convexml
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from colorhash import ColorHash

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    allele_table = pl.read_parquet(args.input)
    rcols = [f'r{i}' for i in range(1, args.n_rcols + 1)]
    rcols = [c for c in rcols if c in allele_table.columns]

    # Filter for SM mode
    if args.mode == 'sm' and args.filter:
        import json
        filter_dict = json.loads(args.filter) if isinstance(args.filter, str) else args.filter
        for col, val in filter_dict.items():
            allele_table = allele_table.filter(pl.col(col) == val)

    # Subsample for testing
    if args.test_mode and args.subsample:
        import json
        factors = json.loads(args.subsample) if isinstance(args.subsample, str) else args.subsample
        for factor in factors:
            col = factor['column']
            top_x = factor.get('top_x')
            sample_n = factor.get('sample_n')
            if col not in allele_table.columns:
                print(f"  {col}: column not found, skipping")
                continue
            pool = (
                allele_table.group_by(col).len()
                .sort('len', descending=True)
                .head(top_x)
                .select(col)
            )
            if sample_n is not None:
                pool = pool.sample(n=min(sample_n, pool.height), seed=42)
            allele_table = allele_table.join(pool, on=col, how='semi')
            print(f"  {col}: top {top_x}{f' -> sampled {sample_n}' if sample_n else ''} -> {allele_table.height} rows")

    n_cells = allele_table['cellBC'].n_unique() if args.mode == 'sc' else allele_table.height
    if n_cells < args.min_cells:
        print(f"Only {n_cells} rows, need {args.min_cells}. Skipping.")
        write_done(outdir, f'too_few_entries_{n_cells}')
        return

    # Handle spanning deletions
    if args.spanning_deletions != 'use':
        allele_table = flag_spanning_deletions(allele_table, rcols, mode=args.spanning_deletions)

    # Fill nulls for cassiopeia
    allele_table = allele_table.with_columns([pl.col(c).fill_null('None') for c in rcols])
    pd_allele = allele_table.to_pandas()

    # Character matrix
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
        pd_allele, allele_rep_thresh=args.allele_rep_thresh,
    )
    print(f"Character matrix: {character_matrix.shape[0]} x {character_matrix.shape[1]}")
    # Cast to str to avoid pyarrow ArrowInvalid with mixed list/non-list values in object columns
    cm_out = character_matrix.reset_index()
    cm_out = cm_out.astype({c: str for c in cm_out.columns if cm_out[c].dtype == object})
    cm_out.to_parquet(outdir / 'character_matrix.parquet')

    cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)

    if cas_tree.n_cell < args.min_cells:
        print(f"Only {cas_tree.n_cell} cells after filtering. Skipping.")
        write_done(outdir, f'too_few_entries_after_filter_{cas_tree.n_cell}')
        return

    # Solve
    match args.solver:
        case 'nj':
            solver = cas.solver.NeighborJoiningSolver(
                dissimilarity_function=cas.solver.dissimilarity.weighted_hamming_distance,
                add_root=True,
            )
        case 'vanilla':
            solver = cas.solver.VanillaGreedySolver()
        case 'mcgs':
            solver = cas.solver.MaxCutGreedySolver()
        case _:
            raise ValueError(f"Unknown solver: {args.solver}")

    print(f"Solving tree ({args.solver}, {cas_tree.n_cell} cells)...")
    solver.solve(cas_tree, collapse_mutationless_edges=True)

    # Branch length estimation
    if not args.skip_branch_lengths:
        print("Estimating branch lengths...")
        sanitize = lambda states: [s[0] if isinstance(s, (tuple, list)) else s for s in states]
        tree_newick = convexml.convexml(
            tree_newick=cas_tree.get_newick(record_branch_lengths=True),
            leaf_sequences={l: sanitize(cas_tree.get_character_states(l)) for l in cas_tree.leaves},
            minimum_branch_length=args.min_branch,
        )['tree_newick']
        print("Updating branch lengths...")
        cas_tree = cas.data.CassiopeiaTree(
            character_matrix=cas_tree.character_matrix,
            tree=tree_newick,
            missing_state_indicator=cas_tree.missing_state_indicator,
        )

    newick = cas_tree.get_newick(record_branch_lengths=True)
    (outdir / 'newick.txt').write_text(newick)

    # Build obs
    if args.mode == 'sc':
        allele_wide = pivot_alleles_wide(allele_table, rcols)
        obs_df = allele_wide.to_pandas().set_index('cellBC')
        leaves_in_obs = [l for l in cas_tree.leaves if l in obs_df.index]
        obs_df = obs_df.loc[leaves_in_obs]
    else:
        cell_id_col = 'UMI' if 'UMI' in allele_table.columns else 'cellBC'
        obs_df = allele_table.select(cell_id_col, *rcols).to_pandas().set_index(cell_id_col)
        obs_df = obs_df.loc[cas_tree.leaves]

    tdata = td.TreeData(
        X=None, allow_overlap=True, obs=obs_df.copy(),
        obst={args.solver: cas_tree.get_tree_topology()},
    )
    add_weighted_depth(tdata, tree_key=args.solver)
    tdata.write_h5td(outdir / 'tdata.h5td')

    md5 = hashlib.md5(str(tdata).encode()).hexdigest()[:8]
    (outdir / 'tdata_md5.txt').write_text(md5)

    n_cells_tree = cas_tree.n_cell
    sample = getattr(args, 'sample_name', 'sample') or 'sample'

    # Annotation columns and colors
    if args.mode == 'sc':
        ann_cols = [c for c in tdata.obs.columns if re.search(r'_r\d+$', c)]
    else:
        ann_cols = rcols

    allele_colors_hex = CassiopeiaLineageExtension._build_allele_palette(tdata.obs, ann_cols)

    n_ann = len(ann_cols)
    dendrogram_ratio = 0.05
    annotation_width = (1 - dendrogram_ratio) / (dendrogram_ratio * max(n_ann, 1))
    prefix = f'{"sc" if args.mode == "sc" else "sm"}_{n_cells_tree}_{sample}_{args.solver}_span-{args.spanning_deletions}'

    # Circular
    print(f"Plotting circular tree ({n_cells_tree} cells, {n_ann} annotation columns)...")
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=900, subplot_kw={'projection': 'polar'})
    pc.pl.tree(tdata, tree=args.solver, keys=ann_cols, polar=True,
               depth_key='depth', extend_branches=False, branch_linewidth=0.15,
               annotation_width=0.01, palette=allele_colors_hex, ax=ax, legend=False)
    fig.savefig(outdir / f'{prefix}_circ.png', bbox_inches='tight')
    plt.close(fig)

    # Linear
    height_inches = max(5, (n_cells_tree / 1000) * 750 / 100)
    print(f"Plotting linear tree ({height_inches:.0f} inches tall)...")
    fig, ax = plt.subplots(1, 1, figsize=(5, height_inches), dpi=900)
    pc.pl.tree(tdata, tree=args.solver, keys=ann_cols, polar=False,
               depth_key='depth', extend_branches=False, branch_linewidth=0.15,
               annotation_width=annotation_width, palette=allele_colors_hex, ax=ax, legend=False)
    fig.savefig(outdir / f'{prefix}_linear.png', bbox_inches='tight')
    plt.close(fig)

    # Standalone legend
    CassiopeiaLineageExtension._save_legend(allele_colors_hex, ann_cols, tdata.obs, outdir / f'{prefix}_legend.png')

    write_done(outdir)
    print(f"Done: {outdir}")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Cassiopeia PEtracer tree building')
    sub = parser.add_subparsers(dest='command')

    tree_p = sub.add_parser('build-tree', help='Build a phylogenetic tree from allele table')
    tree_p.add_argument('--input', required=True, help='Path to allele table parquet')
    tree_p.add_argument('--outdir', required=True, help='Output directory')
    tree_p.add_argument('--mode', choices=['sc', 'sm'], required=True, help='sc=single-cell, sm=single-molecule')
    tree_p.add_argument('--filter', type=str, default=None,
                        help='JSON dict of column filters for SM mode: {"intBC":"TCTGAA","sbc":"TCAAGT"}')
    tree_p.add_argument('--solver', default='nj', choices=['nj', 'vanilla', 'mcgs'])
    tree_p.add_argument('--n-rcols', type=int, default=15, help='Number of r columns')
    tree_p.add_argument('--spanning-deletions', default='unedited')
    tree_p.add_argument('--allele-rep-thresh', type=float, default=1.0)
    tree_p.add_argument('--min-cells', type=int, default=10, help='Minimum cells/molecules for tree')
    tree_p.add_argument('--min-branch', type=float, default=0.01, help='Minimum branch length')
    tree_p.add_argument('--skip-branch-lengths', action='store_true')
    tree_p.add_argument('--sample-name', default='sample', help='Sample name for plot filenames')
    tree_p.add_argument('--test-mode', action='store_true', help='Enable subsampling for testing')
    tree_p.add_argument('--subsample', type=str, default=None,
                        help='JSON list of subsample factors: [{"column":"intBC","top_x":10,"sample_n":3}, ...]')

    args = parser.parse_args()

    if args.command == 'build-tree':
        _build_tree_standalone(args)
    else:
        parser.print_help()
