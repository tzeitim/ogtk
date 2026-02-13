import polars as pl
from pathlib import Path
from typing import Set,Dict, Optional, Type, List, Tuple
from enum import Enum
from .base import PostProcessorExtension
from .registry import extension_registry
from .config import ExtensionConfig
from ..pipeline.types import StepResults,FractureXp
from ..pipeline.formats import scan_file, read_file
from dataclasses import dataclass, field
from ogtk.utils.log import CustomLogger
from ogtk.utils.general import fuzzy_match_str


def collapse_umis_to_cells(
    allele_df: pl.DataFrame,
    cell_col: str = 'cellBC',
    intbc_col: str = 'intBC',
    min_umis_per_cell: int = 4,
    min_umi_agreement: float = 0.5,
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Collapse multiple UMIs per cell into one consensus allele per (cell, intBC).

    For single-cell lineage tracing, each cell may have multiple UMIs for the same
    integration barcode. This function collapses them to one consensus row per
    (cell, intBC) using mode-based voting.

    Args:
        allele_df: DataFrame with cellBC, intBC, r1, r2, ..., readCount columns
        cell_col: Column name for cell barcode
        intbc_col: Column name for integration barcode
        min_umis_per_cell: Minimum UMIs required for valid consensus (default 4).
                          Groups with fewer UMIs are filtered out.
        min_umi_agreement: Minimum fraction for consensus (default 0.5).
                          If the mode has support below this threshold, mark as None (missing).
        logger: Optional logger for info messages

    Returns:
        DataFrame with one row per (cell, intBC), adds n_umis column.
        Cells with fewer than min_umis_per_cell UMIs are excluded.
    """
    if allele_df.height == 0:
        return allele_df.with_columns(pl.lit(0).alias('n_umis'))

    # Identify r columns (allele columns)
    rcols = allele_df.select(pl.col('^r\\d+$')).columns
    if not rcols:
        if logger:
            logger.warning("No r columns found in allele_df, returning as-is with n_umis=1")
        return allele_df.with_columns(pl.lit(1).alias('n_umis'))

    # Columns to preserve (take first value per group)
    preserve_cols = [c for c in allele_df.columns if c not in rcols + [cell_col, intbc_col]]

    # Add row index for tracking
    df = allele_df.with_row_index('_row_idx')

    # Step 1: Get group sizes and filter by min_umis_per_cell
    group_sizes = (
        df.group_by([cell_col, intbc_col])
        .agg(pl.len().alias('n_umis'))
    )

    n_groups_before = group_sizes.height
    valid_groups = group_sizes.filter(pl.col('n_umis') >= min_umis_per_cell)
    n_groups_after = valid_groups.height

    if logger and n_groups_before != n_groups_after:
        logger.info(f"UMI collapse: filtered {n_groups_before - n_groups_after} (cell, intBC) groups "
                   f"with < {min_umis_per_cell} UMIs ({n_groups_after} remaining)")

    if valid_groups.height == 0:
        # Return empty DataFrame with expected columns
        empty_cols = [cell_col, intbc_col, 'n_umis'] + rcols + preserve_cols
        return pl.DataFrame(schema={c: allele_df.schema.get(c, pl.Utf8) for c in empty_cols if c in allele_df.columns or c == 'n_umis'})

    # Filter to only valid groups
    df = df.join(valid_groups.select([cell_col, intbc_col]), on=[cell_col, intbc_col], how='semi')

    # Step 2: For each r column, compute mode and support using native Polars
    # Use unpivot -> group -> compute -> pivot approach for efficiency

    # Unpivot r columns to long format
    id_cols = [cell_col, intbc_col, '_row_idx'] + preserve_cols
    unpivoted = df.unpivot(
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
            # Get value counts and extract mode info
            pl.col('_allele').value_counts(sort=True).first().alias('_mode_struct'),
        ])
        .with_columns([
            # Extract mode value
            pl.col('_mode_struct').struct.field('_allele').alias('_mode_value'),
            # Extract mode count
            pl.col('_mode_struct').struct.field('count').alias('_mode_count'),
        ])
        .with_columns(
            # Calculate support fraction
            (pl.col('_mode_count') / pl.col('_n_total')).alias('_support')
        )
        .with_columns(
            # Apply threshold: if support < min_umi_agreement, set to None
            pl.when(pl.col('_support') >= min_umi_agreement)
            .then(pl.col('_mode_value'))
            .otherwise(pl.lit(None))
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

    # Add n_umis column
    result = result.join(valid_groups, on=[cell_col, intbc_col], how='left')

    # Add preserved columns (take first value per group)
    if preserve_cols:
        preserved = (
            df.group_by([cell_col, intbc_col])
            .agg([pl.col(c).first().alias(c) for c in preserve_cols])
        )
        result = result.join(preserved, on=[cell_col, intbc_col], how='left')

    # Reorder columns to match expected output
    output_cols = [cell_col, intbc_col, 'n_umis'] + rcols + preserve_cols
    output_cols = [c for c in output_cols if c in result.columns]
    result = result.select(output_cols)

    if logger:
        logger.info(f"UMI collapse: {allele_df.height} rows -> {result.height} (cell, intBC) pairs")

    return result


def flag_spanning_deletions(pl_allele: pl.DataFrame, rcols: list[str], mode: str = 'unedited') -> pl.DataFrame:
    """Handle alleles spanning multiple cutsites.

    Detects when the same non-empty allele value appears in multiple r columns
    for the same row, indicating a deletion that spans multiple cutsites.

    Args:
        pl_allele: DataFrame with r1, r2, ... columns containing allele values
        rcols: List of r column names (e.g., ['r1', 'r2', 'r3', ...])
        mode: How to handle spanning deletions:
            - 'unedited': mark as unedited ("None" string, Cassiopeia state 0)
            - 'missing': mark as missing data (null)
            - 'use': leave as-is (keep the deletion allele value)

    Returns:
        DataFrame with spanning deletions handled according to mode
    """
    if not rcols or mode == 'use':
        return pl_allele

    # Cassiopeia treats "NONE" or strings containing "None" as uncut (state 0)
    # Missing data requires passing missing_data_allele to convert_alleletable_to_character_matrix
    replacement_value = 'None' if mode == 'unedited' else None

    # Add row index for tracking
    df = pl_allele.with_row_index('_row_idx')

    # Unpivot r columns to long format
    unpivoted = df.unpivot(
        on=rcols,
        index='_row_idx',
        variable_name='cutsite',
        value_name='allele'
    )

    # Find alleles that appear in multiple cutsites for the same row
    # (non-empty values only)
    spanning = (
        unpivoted
        .filter(pl.col('allele').is_not_null())
        .filter(pl.col('allele') != '')
        .group_by(['_row_idx', 'allele'])
        .agg(pl.len().alias('count'))
        .filter(pl.col('count') > 1)
        .select(['_row_idx', 'allele'])
        .with_columns(pl.lit(True).alias('_is_spanning'))
    )

    # Mark spanning deletions according to mode
    unpivoted = (
        unpivoted
        .join(spanning, on=['_row_idx', 'allele'], how='left')
        .with_columns(
            pl.when(pl.col('_is_spanning') == True)
            .then(pl.lit(replacement_value))
            .otherwise(pl.col('allele'))
            .alias('allele')
        )
        .with_columns(pl.col('allele').cast(pl.Utf8))  # Ensure string type preserved through pivot
        .drop('_is_spanning')
    )

    # Pivot back to wide format
    pivoted = unpivoted.pivot(
        on='cutsite',
        index='_row_idx',
        values='allele'
    )

    # Restore original column order and merge with non-r columns
    other_cols = [c for c in pl_allele.columns if c not in rcols]
    result = (
        df.select(['_row_idx'] + other_cols)
        .join(pivoted, on='_row_idx')
        .drop('_row_idx')
        .select(pl_allele.columns)  # restore original column order
    )

    return result


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
        spanning_deletions: str = 'unedited',
        # Partition filtering parameters
        intbc_whitelist_path: Optional[str] = None,
        min_molecules_per_group: int = 10,
        min_proportion_of_sample: float = 0.02,
        min_ratio_to_max: float = 0.1,
        # Modality parameters
        modality: str = 'single-molecule',
        cbc_len: int = 16,
        # Single-cell collapse parameters
        collapse_to_cells: bool = True,
        min_umis_per_cell: int = 4,
        min_umi_agreement: float = 0.5,
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

    # === END PRE-FILTERING ===

    res = []
    umi_tables = []
    nones = 0

    # Determine partition columns based on modality
    if modality == 'single-cell':
        partition_cols = ['intBC', 'mod', 'cbc']
    else:
        partition_cols = ['intBC', 'mod', 'sbc']

    for partition_key, queries in cass_ldf.collect().partition_by(*partition_cols, as_dict=True).items():
        intBC, mod, group_id = partition_key  # group_id is cbc for single-cell, sbc for single-molecule
        if mod is None:
            nones+=1
        else:
            if logger:
                logger.debug(f"{mod=} {intBC=} {queries.shape=}")

            # Build alignment kwargs - only include gap penalties if explicitly set
            align_kwargs = {
                'queries': queries.to_pandas(),
                'ref_filepath': f'{workdir}/{mod}.fasta',
                'n_threads': 1,
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
            if rcols and 'Seq' in pl_allele.columns and 'CIGAR' in pl_allele.columns:
                pl_allele = pl_allele.with_columns([
                    pl.col(col).cigar.enrich_insertions(pl.col('Seq'), pl.col('CIGAR'))
                    for col in rcols
                ])

            # Handle deletions spanning multiple cutsites
            pl_allele = flag_spanning_deletions(pl_allele, rcols, mode=spanning_deletions)

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
            # group_id is cbc for single-cell, sbc for single-molecule
            if modality == 'single-cell':
                allele_table['cbc'] = group_id
            else:
                allele_table['sbc'] = group_id

            umi_tables.append(umi_table.copy())
            #self.xp.logger.info(f"found {umi_table['intBC'].n_unique()} integrations for {mod}")

            res.append(
                pl.DataFrame(allele_table).with_columns(mod=pl.lit(mod), mols=allele_table.shape[0])
                )

    # Concatenate all partition results
    alleles_pl = pl.concat(res) if res else pl.DataFrame()

    # For single-cell: collapse UMIs to one allele per (cell, intBC)
    # Keep both raw (per-UMI) and collapsed (per-cell) versions
    collapse_metrics = {}
    alleles_pl_collapsed = None
    if modality == 'single-cell' and collapse_to_cells and alleles_pl.height > 0:
        n_rows_before = alleles_pl.height
        alleles_pl_collapsed = collapse_umis_to_cells(
            alleles_pl,
            cell_col='cellBC',
            intbc_col='intBC',
            min_umis_per_cell=min_umis_per_cell,
            min_umi_agreement=min_umi_agreement,
            logger=logger,
        )
        collapse_metrics = {
            'collapse_enabled': True,
            'rows_before_collapse': n_rows_before,
            'rows_after_collapse': alleles_pl_collapsed.height,
            'min_umis_per_cell': min_umis_per_cell,
            'min_umi_agreement': min_umi_agreement,
        }
    else:
        collapse_metrics = {'collapse_enabled': False}

    return StepResults(
            results={
                "alleles_pl": alleles_pl,                      # Raw per-UMI alleles
                "alleles_pl_collapsed": alleles_pl_collapsed,  # Collapsed per-cell (single-cell only)
                "alleles_pd": umi_tables,
            },
            metrics={"partitions_processed": len(res), "none_mod_skipped": nones, **filter_metrics, **collapse_metrics}
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
        .with_columns(
            pl.col('contig').str
            .extract(f'{int_anchor1}(.+?){int_anchor2}', 1).alias('intBC'),
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

    # Single-cell collapse parameters
    # For single-cell lineage tracing, collapse multiple UMIs per cell to one allele per (cell, intBC)
    collapse_to_cells: bool = True              # Enable collapse for single-cell modality
    min_umis_per_cell: int = 4                  # Minimum UMIs for valid consensus (cells with fewer are filtered)
    min_umi_agreement: float = 0.5              # Min fraction for consensus (below -> missing/None)

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

            # Save raw per-UMI alleles
            self.alleles_pl = result.results['alleles_pl']
            self.alleles_pl.write_parquet(f"{self.workdir}/alleles_pl.parquet")
            self.xp.logger.io(f"Saved raw per-UMI alleles to {self.workdir}/alleles_pl.parquet")

            # Save collapsed per-cell alleles (single-cell only)
            self.alleles_pl_collapsed = result.results.get('alleles_pl_collapsed')
            if self.alleles_pl_collapsed is not None:
                self.alleles_pl_collapsed.write_parquet(f"{self.workdir}/alleles_pl_collapsed.parquet")
                self.xp.logger.io(f"Saved collapsed per-cell alleles to {self.workdir}/alleles_pl_collapsed.parquet")

            self.xp.logger.io(f"Saving cassiopeia allele tables to {cass_allele_path}")
            self.ldf.sink_parquet(cass_allele_path)

        if self.should_run_step(CassiopeiaStep.SEGMENTED_ALLELE.value):
            self.xp.logger.info("Running segmented_allele step")
            result = self._segmented_allele()
            self.alleles_segmented = result.results['allele_table_segmented']
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

        # Load assembled if available
        if assembled_path.exists():
            assembled_df = read_file(assembled_path)
        else:
            # Create minimal assembled_df from segments
            assembled_df = segments_df.select(['umi', 'start_meta', 'end_meta']).unique()

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

    def _pycea_explore(self):
        
       """ """ 
       solvers = ['vanilla', 'mcgs', 'nj']
       solvers = ['vanilla']
       #for index, n_cells in [(i, ii.shape) for i,ii in enumerate(umi_tables) if ii.shape[0]>200]:
       for (index, intbc, mod, n_ori, well), allele_table in alg_df.partition_by('xid', 'intBC', 'mod', 'n_ori', 'well', as_dict=True).items():
           allele_table = allele_table.to_pandas()
           tdata = None
           character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
               allele_table,
               allele_rep_thresh = 1.0)
           
           cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix)
           
           print(f"{index = } {intbc=}\tnumber of cells {cas_tree.n_cell}, number of characters {cas_tree.n_character} ")
           meta = '-'.join(pl.DataFrame(allele_table).select('xid', 'well', 'mod', 'intBC', 'n_ori').unique()[0].to_dummies(separator="_").columns)
           collapse = False
           collapse = True
           
           allele_colors_hex = dict(map(lambda x: (x, ColorHash(x).hex), pl.DataFrame(allele_table).unpivot(on=rcols, value_name="allele").get_column('allele').unique()))
           
           allele_matrix = pl.DataFrame(allele_table).select('UMI', '^r\d+$').to_pandas().set_index('UMI')
           allele_matrix
           for solver in solvers :
               match solver:
                   case "shared":
                       solve = cas.solver.SharedMutationJoiningSolver()
                   case "vanilla":
                       solve = cas.solver.VanillaGreedySolver()
                   case "UPGMA":
                       solve = cas.solver.SharedMutationJoiningSolver()
                   case "mcgs":
                       solve = cas.solver.MaxCutGreedySolver()
                   case "mc": # this one is very slow
                       solve = cas.solver.MaxCutSolver()
                   case "nj":
                       solve = cas.solver.NeighborJoiningSolver(add_root=True)

               print(solver)
               solve.solve(cas_tree, collapse_mutationless_edges=collapse)
               if tdata is None:
                   tdata = td.TreeData(
                       X=None,  
                       allow_overlap=True,
                       obs=allele_matrix.loc[cas_tree.leaves], 
                      # obsm={"alleles": character_matrix.loc[cas_tree.leaves].values}, # can't use this since must be encoded, e.g. character_matrix
                       obst={solver: cas_tree.get_tree_topology()}
                   )
               
               name = f'pctree_{cas_tree.n_cell}_{meta}_{solver}_{collapse}'

               tdata.obst[solver] = cas_tree.get_tree_topology()
               pycea.pp.add_depth(tdata, tree=solver)

               if True:
                   fig, ax = plt.subplots(1,1, figsize=(10, 100))
                   pc.pl.tree(tdata, 
                          tree=solver,
                     keys=rcols,
                     polar=False, 
                     extend_branches=True,
                     palette=allele_colors_hex,
                              branch_linewidth=0.15,
                     ax=ax)
                   fig.savefig(f'{name}.png')
                   
                   if False:
                       fig, ax = plt.subplots(1,1, figsize=(10, 10), dpi=1200, subplot_kw={"projection": "polar"})
                       pc.pl.tree(tdata, 
                              tree=solver,
                         keys=rcols,
                         polar=True, 
                         extend_branches=True,
                                  branch_linewidth=0.15,
                         palette=allele_colors_hex,
                         ax=ax)
                       fig.savefig(f'{name}_circ.png')

                   plt.close('all')
                
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
