import polars as pl
# rogtk includes the "dna" namespace
import rogtk

from ogtk.utils.log import Rlogger, call
from ogtk.utils.general import fuzzy_match_str
from .masking import (
    generate_mask_PEtracer_expression,
    generate_unmask_PEtracer_expression,
    generate_mask_flanks_expression,
    generate_unmask_flanks_expression,
)
from .segmentation import (
    segment_by_metas,
    stitch_segments,
    flag_aberrant_molecules,
    generate_segmentation_report,
    format_segmentation_report,
    plot_segmentation_qc,
    CASSETTE_START_MARKER,
    CASSETTE_END_MARKER,
)

__all__ = [
        'PlDNA',
        'PlPipeline',
        'LazyKmer',
        'EagerKmer',
]

@pl.api.register_lazyframe_namespace("kmer")
class LazyKmer:
    """ """
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf

    def explode_kmers(self, seq_field: str, k = 15, only_unique = True, max_str_len: int = 500):
        """Returns lazy dataframe with the field ``kmer``"""
        return (
            self._ldf
            .with_columns(
                kmer= pl.concat_list([pl.col(seq_field).str.slice(s, k) for s in range(0, max_str_len)]))
            .explode('kmer')
            .filter(pl.col('kmer').str.len_chars()==k)
            .filter(~pl.lit(only_unique) | pl.col('kmer').is_unique())
        )

@pl.api.register_dataframe_namespace("kmer")
class EagerKmer:
    """ """
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df

    def explode_kmers(self, seq_field: str, k = 15, only_unique = True, max_str_len: int= 500):
        """Returns dataframe with the field ``kmer``"""
        return self._df.lazy().kmer.explode_kmers(seq_field, k, only_unique, max_str_len).collect() #pyright: ignore


@pl.api.register_dataframe_namespace("dna")
class PlDNA:
    """Provides DNA sequence manipulation methods for Polars DataFrames.

    Methods are registered under the 'dna' namespace and can be accessed via df.dna.*
    """
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df

    def to_fasta(self,
         read_id_col: str,
         read_col: str) -> pl.DataFrame:
        """Converts reads to FASTA format.

        Generates FASTA formatted strings by combining read IDs and sequences.

        Args:
            read_id_col (str): Column containing read identifiers
            read_col (str): Column containing DNA sequences

        Returns:
            pl.DataFrame: DataFrame with new column '{read_col}_fasta' containing FASTA formatted strings
        """

        return  self._df.with_columns(
                (">"+pl.col(read_id_col)\
                 +"\n"+pl.col(read_col)
                 )
                 .alias(f'{read_col}_fasta')
        )
    def to_fastq(self,
         read_id_col: str,
         read_qual_col: str|None,
         read_col: str)-> pl.DataFrame:
        """Converts reads to FASTQ format.

        Generates FASTQ formatted strings by combining read IDs, sequences, and quality scores.

        Args:
            read_id_col (str): Column containing read identifiers
            read_qual_col (str|None): Column containing quality scores. If None, generates fake scores
            read_col (str): Column containing DNA sequences

        Returns:
            pl.DataFrame: DataFrame with new column '{read_col}_fastq' containing FASTQ formatted strings
        """

        if read_qual_col is None:
            read_qual_col = 'fake_qual'

            return (
                    self._df
                    .with_columns(pl.col(read_col).str.replace_all('.', "I").alias('fake_qual'))
                    .with_columns(
                    ("@"+pl.col(read_id_col)\
                     +"\n"+pl.col(read_col)\
                     +"\n+"\
                     +"\n"+pl.col(read_qual_col)
                    )
                     .alias(f'{read_col}_fastq')
                )
            )

        return self._df.with_columns(
                ("@"+pl.col(read_id_col)\
                 +"\n"+pl.col(read_col)\
                 +"\n+"\
                 +"\n"+pl.col(read_qual_col)
                )
                 .alias(f'{read_col}_fastq')
        )


@pl.api.register_lazyframe_namespace("dna")
class PllDNA:
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf

    def kmers_explode(self, seq_column:str ,k: int = 10, max_str_len=500):
        return(
                self._ldf
                .with_columns(
                    pl.concat_list(
                        [pl.col(seq_column).str.slice(s, k) for s in range(0, max_str_len)]
                        )
                    )
            .explode('seq').filter(pl.col('seq').is_unique())
            )

@pl.api.register_dataframe_namespace("pp")
class PlPipeline:
    """Provides sequence assembly and optimization methods for Polars DataFrames.

    Methods are registered under the 'pp' namespace and can be accessed via df.pp.*
    """
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df
   
    @call
    def assembly_with_opt(self,
                         start_k: int = 25,
                         start_min_coverage: int = 17,
                         method: str = "shortest_path",
                         start_anchor: str | None = None,
                         end_anchor: str | None = None,
                         min_length: int | None = None,
                         export_graphs: bool = False,
                         prioritize_length: bool = False,
                         min_reads: int = 100) -> pl.DataFrame:
        """Optimize assembly parameters for sequences with given anchors.

        Performs De Bruijn graph assembly with parameter optimization.

        Args:
            start_k (int): Initial k-mer size for graph construction
            start_min_coverage (int): Initial minimum k-mer coverage threshold
            method (str): Assembly method ('compression' or 'shortest_path')
            start_anchor (str|None): Start sequence anchor for shortest_path method
            end_anchor (str|None): End sequence anchor for shortest_path method 
            min_length (int|None): Minimum contig length to return
            export_graphs (bool): Whether to export graph visualization files
            prioritize_length (bool): Optimize for length over anchor presence
            min_reads (int): Minimum number of reads to attempt assembly

        Returns:
            pl.DataFrame: DataFrame with assembly results including contig sequences and parameters
        """
        return (
            self._df
             .with_columns(pl.col('r2_seq').str.replace(f'^.*{start_anchor}', start_anchor))
             .with_columns(pl.col('r2_seq').str.replace(f'{end_anchor}.*$', end_anchor))
            .filter(pl.col('reads') > min_reads)
            .group_by(['umi']).agg(
                rogtk.optimize_assembly(
                    expr=pl.col('r2_seq'),
                    start_k=start_k,
                    start_min_coverage=start_min_coverage,
                    method=method,
                    start_anchor=start_anchor,
                    end_anchor=end_anchor,
                    min_length=min_length,
                    export_graphs=export_graphs,
                    prioritize_length=prioritize_length))
            .unnest('r2_seq')
            .with_columns((pl.col('length')==0).alias('failed'))
        )

    def assemble_umi(self,
                    target_umi: str,
                    k: int = 15,
                    min_coverage: int = 20, 
                    method: str = "shortest_path",
                    start_anchor: str | None = None,
                    end_anchor: str | None = None,
                    min_length: int | None = None,
                    auto_k: bool = False,
                    export_graphs: bool = False,
                    only_largest: bool = True,
                    prefix: str | None = None) -> pl.DataFrame:
        """Assemble sequences for a specific UMI using de Bruijn graphs.

        Performs targeted assembly of reads sharing a UMI.

        Args:
            target_umi (str): UMI sequence to assemble
            k (int): K-mer size for graph construction (used if auto_k=False)
            min_coverage (int): Minimum k-mer coverage threshold
            method (str): Assembly method ('compression' or 'shortest_path')
            start_anchor (str|None): Start sequence anchor for shortest_path method
            end_anchor (str|None): End sequence anchor for shortest_path method
            min_length (int|None): Minimum contig length to return
            auto_k (bool): Automatically estimate optimal k-mer size
            export_graphs (bool): Whether to export graph visualization files
            only_largest (bool): Return only the largest contig
            prefix (str|None): Prefix for output files

        Returns:
            pl.DataFrame: DataFrame with assembly results per UMI
        """

        return (
            self._df
            .filter(pl.col('umi')==target_umi)
             .with_columns(pl.col('r2_seq').str.replace(f'^.*{start_anchor}', start_anchor))
             .with_columns(pl.col('r2_seq').str.replace(f'{end_anchor}.*$', end_anchor))
            .group_by(['umi']).agg(
                rogtk.assemble_sequences(
                    expr=pl.col("r2_seq"),
                    k=k,
                    min_coverage=min_coverage,
                    method=method,
                    start_anchor=start_anchor,
                    end_anchor=end_anchor,
                    min_length=min_length,
                    auto_k=auto_k,
                    only_largest=only_largest,
                    export_graphs=export_graphs,
                    prefix=prefix or f"{target_umi}_"
                ).alias('contig')
            )
        )

    def assemble_umis(self,
                    k: int = 15,
                    min_coverage: int = 20,
                    method: str = "shortest_path",
                    start_anchor: str | None = None,
                    end_anchor: str | None = None,
                    min_length: int | None = None,
                    auto_k: bool = False,
                    export_graphs: bool = False,
                    groups = ['umi', 'sbc'],
                    modality: str = 'single-molecule',
                    do_trim: bool = True,
                    only_largest: bool = True) -> pl.DataFrame:
        """Assemble UMI sequences (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.assemble_umis(
                k=k,
                min_coverage=min_coverage,
                method=method,
                start_anchor=start_anchor,
                end_anchor=end_anchor,
                min_length=min_length,
                auto_k=auto_k,
                export_graphs=export_graphs,
                groups=groups,
                modality=modality,
                do_trim=do_trim,
                only_largest=only_largest
            )
            .collect()
        )

    def sweep_assembly_params(self,
                             target_umi: str,
                             k_start: int = 5,
                             k_end: int = 35,
                             k_step: int = 1,
                             cov_start: int = 1,
                             cov_end: int = 250,
                             cov_step: int = 1,
                             method: str = "shortest_path",
                             start_anchor: str | None = None,
                             end_anchor: str | None = None,
                             min_length: int | None = None,
                             export_graphs: bool = False,
                             auto_k: bool = False,
                             prefix: str | None = None) -> pl.DataFrame:

        """Run sequence assembly across ranges of k-mer size and coverage parameters.

        Performs parameter sweep to identify optimal assembly settings.

        Args:
            target_umi (str): UMI sequence to assemble
            k_start (int): Starting k-mer size
            k_end (int): Ending k-mer size
            k_step (int): Step size for k-mer iteration
            cov_start (int): Starting minimum coverage
            cov_end (int): Ending minimum coverage
            cov_step (int): Step size for coverage iteration
            method (str): Assembly method ('compression' or 'shortest_path')
            start_anchor (str|None): Start sequence anchor for shortest_path method
            end_anchor (str|None): End sequence anchor for shortest_path method
            min_length (int|None): Minimum contig length to return
            export_graphs (bool): Whether to export graph visualization files
            auto_k (bool): Automatically estimate optimal k-mer size
            prefix (str|None): Prefix for output files

        Returns:
            pl.DataFrame: DataFrame with assembly results for each parameter combination
        """
        return (
            self._df
            .filter(pl.col('umi')==target_umi)
             .with_columns(pl.col('r2_seq').str.replace(f'^.*{start_anchor}', start_anchor))
             .with_columns(pl.col('r2_seq').str.replace(f'{end_anchor}.*$', end_anchor))
            .group_by(['umi']).agg(
                rogtk.sweep_assembly_params(
                    expr=pl.col("r2_seq"),
                    k_start=k_start,
                    k_end=k_end,
                    k_step=k_step,
                    cov_start=cov_start,
                    cov_end=cov_end,
                    cov_step=cov_step,
                    method=method,
                    start_anchor=start_anchor,
                    end_anchor=end_anchor,
                    min_length=min_length,
                    export_graphs=export_graphs,
                    auto_k=auto_k,
                    prefix=prefix or f"{target_umi}_"
                ).alias("assembly_results")
            )
            .explode("assembly_results")
            .unnest("assembly_results")
        )

    def ont_to_paired_format(self, umi_len: int, sbc_len: int, anchor: str, anchor_ont:str, modality: str = 'single-molecule', cbc_len: int = 16) -> pl.DataFrame:
        """Convert ONT BAM data to FRACTURE's expected R1/R2 format (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.ont_to_paired_format(umi_len=umi_len, sbc_len=sbc_len, anchor_orient=anchor, anchor_ont=anchor_ont, modality=modality, cbc_len=cbc_len)
            .collect()
        )

    # TODO: Consider using **kwargs instead of explicit parameters to avoid duplication.
    # Tradeoff: **kwargs is DRYer but loses IDE autocomplete and type hints.
    # Current approach: Explicit parameters provide better developer experience.
    def mask_repeats(self,
                     features_csv: str,
                     column_name: str = 'r2_seq',
                     fuzzy_pattern: bool = True,
                     fuzzy_kwargs: dict = None,
                     alias: str | None = None) -> pl.DataFrame:
        """Deprecated: Use mask_flanks() instead. Masks whole TARGET sequences."""
        return (
            self._df
            .lazy()
            .pp.mask_repeats(features_csv=features_csv, column_name=column_name, fuzzy_pattern=fuzzy_pattern, fuzzy_kwargs=fuzzy_kwargs, alias=alias)
            .collect()
        )

    def unmask_repeats(self,
                       features_csv: str,
                       column_name: str = 'contig',
                       fuzzy_pattern: bool = True,
                       fuzzy_kwargs: dict = None,
                       alias: str | None = None) -> pl.DataFrame:
        """Deprecated: Use unmask_flanks() instead. Restores whole TARGET sequences."""
        return (
            self._df
            .lazy()
            .pp.unmask_repeats(features_csv=features_csv, column_name=column_name, fuzzy_pattern=fuzzy_pattern, fuzzy_kwargs=fuzzy_kwargs, alias=alias)
            .collect()
        )

    def mask_flanks(self,
                    features_csv: str,
                    column_name: str = 'r2_seq',
                    fuzzy_pattern: bool = True,
                    fuzzy_kwargs: dict = None,
                    alias: str | None = None) -> pl.DataFrame:
        """Mask TARGET flanks based on META context (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.mask_flanks(features_csv=features_csv, column_name=column_name, fuzzy_pattern=fuzzy_pattern, fuzzy_kwargs=fuzzy_kwargs, alias=alias)
            .collect()
        )

    def unmask_flanks(self,
                      features_csv: str,
                      column_name: str = 'contig',
                      fuzzy_pattern: bool = True,
                      fuzzy_kwargs: dict = None,
                      alias: str | None = None) -> pl.DataFrame:
        """Restore original TARGET flanks from scrambled versions (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.unmask_flanks(features_csv=features_csv, column_name=column_name, fuzzy_pattern=fuzzy_pattern, fuzzy_kwargs=fuzzy_kwargs, alias=alias)
            .collect()
        )

    def segment_by_metas(self,
                         metas: pl.DataFrame,
                         cassette_start_anchor: str,
                         cassette_end_anchor: str,
                         seq_col: str = 'r2_seq',
                         keep_cols: list | None = None) -> pl.DataFrame:
        """Split reads at META boundaries into overlapping segments (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.segment_by_metas(
                metas=metas,
                cassette_start_anchor=cassette_start_anchor,
                cassette_end_anchor=cassette_end_anchor,
                seq_col=seq_col,
                keep_cols=keep_cols
            )
            .collect()
        )

    def stitch_segments(self,
                        metas: pl.DataFrame,
                        seq_col: str = 'consensus_seq',
                        cassette_start_anchor_len: int | None = None,
                        group_cols: list[str] | None = None) -> pl.DataFrame:
        """Rejoin assembled segments into full sequences (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.stitch_segments(metas=metas, seq_col=seq_col, cassette_start_anchor_len=cassette_start_anchor_len, group_cols=group_cols)
            .collect()
        )

    def flag_aberrant_molecules(self,
                                metas: pl.DataFrame,
                                seq_col: str = 'r2_seq') -> pl.DataFrame:
        """Flag molecules with tandem duplications (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.flag_aberrant_molecules(metas=metas, seq_col=seq_col)
            .collect()
        )

@pl.api.register_lazyframe_namespace("pp")
class PllPipeline:
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf
        self.logger = Rlogger().get_logger()

    def assemble_umis(self,
                    col_name: str = "r2_seq",
                    k: int = 15,
                    min_coverage: int = 20,
                    method: str = "shortest_path",
                    start_anchor: str | None = None,
                    end_anchor: str | None = None,
                    min_length: int | None = None,
                    auto_k: bool = False,
                    export_graphs: bool = False,
                    groups = ['umi', 'sbc'],
                    modality: str = 'single-molecule',
                    do_trim: bool = True,
                    only_largest: bool = True) -> pl.LazyFrame:
        """
        Assemble UMI sequences with modality-aware grouping.

        For single-cell: groups by 'umi' only (compound cbc+umi)
        For single-molecule: groups by 'umi' and 'sbc'
        """
        # Adjust grouping based on modality
        if modality == 'single-cell':
            # For single-cell, group only by UMI (which contains cbc+umi compound)
            effective_groups = ['umi']
        else:
            # Use provided groups or default for single-molecule
            effective_groups = groups

        # Handle anchor trimming based on method and anchor availability
        ldf = self._ldf

        # Set anchors to None for compression method
        assembly_start_anchor = start_anchor
        assembly_end_anchor = end_anchor

        if method == "compression":
            assembly_start_anchor = None
            assembly_end_anchor = None
        elif method == "shortest_path" and do_trim and start_anchor and end_anchor:
            # Only trim if method requires it AND anchors are provided
            ldf = (
                ldf
                .with_columns(pl.col(col_name).str.replace(f'^.*{start_anchor}', start_anchor))
                .with_columns(pl.col(col_name).str.replace(f'{end_anchor}.*$', end_anchor))
            )

        return (
            ldf
            .group_by(effective_groups).agg(
                rogtk.assemble_sequences(
                    expr=pl.col(col_name),
                    k=k,
                    min_coverage=min_coverage,
                    method=method,
                    start_anchor=assembly_start_anchor,
                    end_anchor=assembly_end_anchor,
                    min_length=min_length,
                    auto_k=auto_k,
                    only_largest=only_largest,
                    export_graphs=export_graphs,
                ).alias('contig')
            )
        )

    # [12N][6N]
    # [  18N  ]
    def parse_reads(self, umi_len, anchor_ont, sbc_len, modality='single-molecule', cbc_len=16):
        '''
        Parse reads from merged R1/R2 format, extracting UMI and barcode information.

        Args:
            umi_len: UMI length
            anchor_ont: ONT anchor sequence
            sbc_len: Sample barcode length (0 for single-cell)
            modality: 'single-cell', 'single-molecule', or 'single-molecule-v2'
                - single-cell: [CBC][UMI][anchor_ont] - compound UMI = CBC + UMI
                - single-molecule: [UMI][SBC][anchor_ont]
                - single-molecule-v2: [SBC][UMI][anchor_ont]
            cbc_len: Cell barcode length (for single-cell data)
        '''
        # metric_field_name represent QC columns
        if modality == 'single-cell':
            # For single-cell: r1_seq structure is cbc + umi
            # Extract compound UMI (cell barcode + umi) for downstream processing
            expected_len = cbc_len + umi_len
            return(
                self._ldf
                    .with_columns(pl.col('r1_seq').str.slice(0, cbc_len).alias('cbc'))
                    .with_columns(pl.col('r1_seq').str.slice(cbc_len, umi_len).alias('umi_only'))
                    .with_columns((pl.col('cbc') + pl.col('umi_only')).alias('umi'))  # Compound UMI
                    .with_columns(pl.lit('').alias('sbc'))  # No sample barcode for single-cell
                    .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
                    # read level metrics
                    .with_columns(metric_reads_ont = pl.col('ont').sum())
                    .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
                    .with_columns(metric_fraction_ont = pl.col('ont').mean())
                    .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
                    # UMI validation - check if prefix length matches expected
                    .with_columns(pl.col('r1_seq')
                                  .str.extract(f'(.*?){anchor_ont}$',1) 
                                  .str.len_chars()
                                  .fill_null(0)
                                  .alias('valid_umi')
                                  )
                    .with_columns((pl.col('valid_umi')==expected_len).alias('valid_umi'))
                    .with_columns(pl.len().over('umi').alias('reads'))
            )
        elif modality == 'single-molecule-v2':
            # For single-molecule-v2: r1_seq structure is sbc + umi (reversed order)
            expected_len = sbc_len + umi_len
            return(
                self._ldf
                    .with_columns(pl.col('r1_seq').str.slice(0, sbc_len).alias('sbc'))
                    .with_columns(pl.col('r1_seq').str.slice(sbc_len, umi_len).alias('umi'))
                    .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
                    # read level metrics
                    .with_columns(metric_reads_ont = pl.col('ont').sum())
                    .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
                    .with_columns(metric_fraction_ont = pl.col('ont').mean())
                    .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
                    # UMI validation
                    .with_columns(pl.col('r1_seq')
                                  .str.extract(f'(.*?){anchor_ont}$',1)
                                  .str.len_chars()
                                  .fill_null(0)
                                  .alias('valid_umi')
                                  )
                    .with_columns((pl.col('valid_umi')==expected_len).alias('valid_umi'))
                    .with_columns(pl.len().over('umi').alias('reads'))
            )
        else:
            # For single-molecule: r1_seq structure is umi + sbc
            expected_len = umi_len + sbc_len
            return(
                self._ldf
                    .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
                    .with_columns(pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc'))
                    .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
                    # read level metrics
                    .with_columns(metric_reads_ont = pl.col('ont').sum())
                    .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
                    .with_columns(metric_fraction_ont = pl.col('ont').mean())
                    .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
                    # UMI validation
                    .with_columns(pl.col('r1_seq')
                                  .str.extract(f'(.*?){anchor_ont}$',1)
                                  .str.len_chars()
                                  .fill_null(0)
                                  .alias('valid_umi')
                                  )
                    .with_columns((pl.col('valid_umi')==expected_len).alias('valid_umi'))
                    .with_columns(pl.len().over('umi').alias('reads'))
            )    

    def parse_read1(self, anchor_ont):
        """ Optional step for parsing also read 1. 
            This is needed to extract additional information when R1 is long
            or when fragmentation didn't ocurr and the pipeline is used for error correction
            TODO: currently it's just a hack that replaces r2_seq keeping the rest of the fields equal
                e.g. the r2_qual won't match the r2_seq field
                The output of this function should be concatenated to the original df

        """
        return (
                self._ldf
                    .with_columns(pl.col('r1_seq')
                                  .str.extract(f'.*?{anchor_ont}(.*)$')
                                  .alias('r2_seq')
                                  )
                )

    def ont_to_paired_format(self, umi_len: int, sbc_len: int, anchor_orient: str, anchor_ont:str, modality: str = 'single-molecule', cbc_len: int = 16, wildcard= '.{0,4}', do_rev_comp: bool = False) -> pl.LazyFrame:
        """Convert ONT BAM data to FRACTURE's expected R1/R2 format (LazyFrame)

        Args:
            umi_len: Length of UMI sequence
            sbc_len: Length of sample barcode (set to 0 for single-cell data)
            anchor_orient: Anchor sequence for read orientation
            anchor_ont: ONT anchor sequence for trimming
            modality: 'single-cell', 'single-molecule', or 'single-molecule-v2'
            cbc_len: Length of cell barcode (only used for single-cell data)
            wildcard: Regex wildcard pattern for fuzzy matching
            do_rev_comp: Apply reverse complement to R2 sequence
        """
        anchor_ont_l = len(anchor_ont)
        
        # For single-cell data, the compound UMI structure is: cell_barcode + umi
        # For single-molecule data, the structure is: umi + sample_barcode
        if modality == 'single-cell':
            # For single-cell: cbc + umi (no sample barcode)
            trim_range = cbc_len + umi_len + anchor_ont_l
            self.logger.debug(f"Single-cell mode: {cbc_len=} {umi_len=} {anchor_ont_l=} {trim_range=}")
        else:
            # For single-molecule: umi + sbc
            trim_range = umi_len + sbc_len + anchor_ont_l
            self.logger.debug(f"Single-molecule mode: {umi_len=} {sbc_len=} {anchor_ont_l=} {trim_range=}")
        
        self.logger.debug(f"{modality=} {anchor_ont=}")

        fuzzy_anchor_pattern = fuzzy_match_str(anchor_orient, wildcard=wildcard)

        return(
        self._ldf
            # Step 1: Orient reads using fuzzy_anchor_pattern 
            .with_columns([
                pl.when(pl.col('sequence').str.contains(fuzzy_anchor_pattern))
                .then(pl.col('sequence'))
                .when(pl.col('sequence').dna.reverse_complement().str.contains(fuzzy_anchor_pattern)) #pyright:ignore
                .then(pl.col('sequence').dna.reverse_complement()) #pyright:ignore
                .otherwise(pl.col('sequence'))
                .alias('oriented_sequence'),
                pl.when(pl.col('sequence').str.contains(fuzzy_anchor_pattern))
                .then(pl.col('quality_scores'))
                .when(pl.col('sequence').dna.reverse_complement().str.contains(fuzzy_anchor_pattern)) #pyright:ignore
                .then(pl.col('quality_scores').str.reverse())
                .otherwise(pl.col('quality_scores'))
                .alias('oriented_quality')
            ])

            # Step 2: Find adapter position - extract everything after the match
            .with_columns([
                pl.col('oriented_sequence').str.extract(
                    f'^(.*)({fuzzy_anchor_pattern})(.*)$',
                    group_index=3
                ).alias('after_adapter')
            ])
            .with_columns([
                # Calculate trim position: total_length - after_adapter_length
                (pl.col('oriented_sequence').str.len_chars() - pl.col('after_adapter').str.len_chars()).alias('trim_pos')
            ])

            # Step 3: Create normalized sequence for visualization (optional)
           # .with_columns([
           #     pl.when(pl.col('after_adapter').is_not_null())
           #     .then(
           #         pl.col('oriented_sequence').str.slice(0, pl.col('trim_pos') - len(anchor_orient)) +
           #         anchor_orient +
           #         pl.col('after_adapter')
           #     )
           #     .otherwise(pl.col('oriented_sequence'))
           #     .alias('normalized_sequence')
           # ])

            # Step 4: Use after_adapter directly as trimmed sequence
            .with_columns([
                pl.when(pl.col('after_adapter').is_not_null())
                .then(pl.col('after_adapter'))
                .otherwise(pl.col('oriented_sequence'))
                .alias('trimmed_sequence'),
                # Trim quality scores from the same position
                pl.when(pl.col('trim_pos').is_not_null())
                .then(pl.col('oriented_quality').str.slice(pl.col('trim_pos')))
                .otherwise(pl.col('oriented_quality'))
                .alias('trimmed_quality')
            ])

            # Step 5: Extract UMI/SBC from clean sequences
            .with_columns([
                pl.col('trimmed_sequence').str.slice(0, trim_range).alias('r1_seq'),
                pl.when(do_rev_comp)
                .then(pl.col('trimmed_sequence').str.slice(trim_range).dna.reverse_complement()) #pyright:ignore
                .otherwise(pl.col('trimmed_sequence').str.slice(trim_range))
                .alias('r2_seq'),
                pl.col('trimmed_quality').str.slice(0, trim_range).alias('r1_qual'),
                pl.when(do_rev_comp)
                .then(pl.col('trimmed_quality').str.slice(trim_range).str.reverse())
                .otherwise(pl.col('trimmed_quality').str.slice(trim_range))
                .alias('r2_qual')
            ])
            .rename({'name':'read_id'})
            #.select(['read_id', 'r1_seq', 'r1_qual', 'r2_seq', 'r2_qual', 'normalized_sequence', 'trim_pos'])
            .select(['read_id', 'r1_seq', 'r1_qual', 'r2_seq', 'r2_qual', 'trim_pos'])
        )

    def mask_repeats(self,
                     features_csv: str,
                     column_name: str = 'r2_seq',
                     fuzzy_pattern: bool = True,
                     fuzzy_kwargs: dict = None,
                     alias: str | None = None) -> pl.LazyFrame:
        """
        .. deprecated::
            Use `mask_flanks()` instead, which handles both edited and unedited
            cassettes universally.

        Mask repetitive sequences based on META-TARGET pairs.

        This is useful for very long cassettes with direct repeats that exceed
        the kmer length limit (typically 64 bases in rogtk). By replacing
        repetitive TARGET sequences with deterministic scrambled versions based
        on their preceding META context, the kmer-based assembler can work on
        unique/variable regions while preserving structural information.

        Args:
            features_csv: Path to CSV file with META and TARGET sequences

                         Expected CSV format (3 columns: feature, seq, kind):

                         feature,seq,kind
                         META01,AGAAGCCGTGTGCCGGTCTA,META
                         META02,ATCGTGCGGACGAGACAGCA,META
                         RNF2,TGGCAGTCATCTTAGTCATTACGACAGGTGTTCGTTGTAACTCATATA,TARGET
                         HEK3,CTTGGGGCCCAGACTGAGCACGACTTGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCA,TARGET

                         - META sequences are anchors that precede TARGET sequences
                         - TARGET sequences will be replaced with scrambled versions when preceded by a META
                         - Scrambling is deterministic: same META+TARGET pair always produces same result

            column_name: Column containing sequences to mask (default: 'r2_seq')
            fuzzy_pattern: Whether to use fuzzy matching for sequencing errors (default: True)
            fuzzy_kwargs: Optional dict with fuzzy matching parameters:
                         - wildcard: Regex for character substitution (default: '.{0,1}')
                         - include_original: Include exact match (default: True)
                         - sep: Separator for alternatives (default: '|')
                         - max_length: Only fuzzify sequences up to this length (default: 200)
            alias: Optional output column name. If None, replaces column_name in place (default: None)

        Returns:
            LazyFrame with masked sequences

        Example:
            # Basic usage - replaces r2_seq in place
            df = (
                pl.scan_parquet('parsed_reads.parquet')
                .pp.mask_repeats('features.csv', column_name='r2_seq')
                .pp.assemble_umis(
                    k=15,
                    min_coverage=20,
                    start_anchor='AAGGTTAAAGAACGACTTCC',
                    end_anchor='GCCTAAAACTGCTCACCTAT'
                )
            )

            # Create new column - preserves original
            df = (
                pl.scan_parquet('parsed_reads.parquet')
                .pp.mask_repeats('features.csv', column_name='r2_seq', alias='masked')
                .with_columns(r2_seq=pl.col('masked'))  # Use masked for assembly
                .pp.assemble_umis(
                    k=15,
                    min_coverage=20,
                    start_anchor='AAGGTTAAAGAACGACTTCC',
                    end_anchor='GCCTAAAACTGCTCACCTAT'
                )
            )

        See Also:
            ogtk.ltr.fracture.pipeline.masking.generate_mask_PEtracer_expression
        """
        mask_expr = generate_mask_PEtracer_expression(
            features_csv=features_csv,
            column_name=column_name,
            fuzzy_pattern=fuzzy_pattern,
            fuzzy_kwargs=fuzzy_kwargs,
            logger=self.logger
        )
        output_name = alias if alias is not None else column_name
        return self._ldf.with_columns(mask_expr.alias(output_name))

    def unmask_repeats(self,
                       features_csv: str,
                       column_name: str = 'contig',
                       fuzzy_pattern: bool = True,
                       fuzzy_kwargs: dict = None,
                       alias: str | None = None) -> pl.LazyFrame:
        """
        .. deprecated::
            Use `unmask_flanks()` instead, which handles both edited and unedited
            cassettes universally.

        Restore original TARGET sequences from scrambled versions after assembly.

        This reverses the masking operation performed by mask_repeats(). After assembly,
        contigs contain scrambled TARGET sequences which need to be restored to their
        original form for downstream analysis.

        Args:
            features_csv: Path to CSV file with META and TARGET sequences
                         (same file used during masking)
            column_name: Column containing sequences to unmask (default: 'contig')
            fuzzy_pattern: Whether fuzzy matching was used during masking (default: True)
            fuzzy_kwargs: Optional dict with fuzzy matching parameters used during masking
                         Must match the parameters used in mask_repeats()
            alias: Optional output column name. If None, replaces column_name in place (default: None)

        Returns:
            LazyFrame with restored original sequences

        Example:
            # After assembly, restore original sequences in place
            df_contigs = (
                assembled_contigs
                .pp.unmask_repeats('features.csv', column_name='contig')
            )

            # Create new column preserving scrambled version
            df_contigs = (
                assembled_contigs
                .pp.unmask_repeats('features.csv', column_name='contig', alias='unmasked')
            )

        Note:
            The fuzzy_pattern and fuzzy_kwargs parameters must match what was used
            during masking, otherwise the scrambled sequences won't be recognized.

        See Also:
            ogtk.ltr.fracture.pipeline.masking.generate_unmask_PEtracer_expression
        """
        unmask_expr = generate_unmask_PEtracer_expression(
            features_csv=features_csv,
            column_name=column_name,
            fuzzy_pattern=fuzzy_pattern,
            fuzzy_kwargs=fuzzy_kwargs,
            logger=self.logger
        )
        output_name = alias if alias is not None else column_name
        return self._ldf.with_columns(unmask_expr.alias(output_name))

    def mask_flanks(self,
                    features_csv: str,
                    column_name: str = 'r2_seq',
                    fuzzy_pattern: bool = True,
                    fuzzy_kwargs: dict = None,
                    alias: str | None = None) -> pl.LazyFrame:
        """
        Mask TARGET flanks based on META context for edited and unedited cassettes.

        This function handles both edited (with lineage marks) and unedited cassettes
        by masking the flanking sequences around potential insertion points. Unlike
        mask_repeats() which masks whole TARGETs, this function:
        - Splits TARGETs into LEFT_FLANK and RIGHT_FLANK
        - Scrambles each flank independently based on META context
        - Preserves any inserted lineage marks (0-5bp between flanks)

        Args:
            features_csv: Path to CSV file with META and EDITED_TARGET sequences

                         Expected CSV format (6 columns):

                         feature,seq,left_flank,right_flank,marks,kind
                         META01,AGAAGCCGTGTGCCGGTCTA,,,META
                         RNF2,CATCTTAGTCATTACGACAGGTGTTCGTTG,CATCTTAGTCATTAC,GACAGGTGTTCGTTG,ACTGT|TAAGT|...,EDITED_TARGET

                         - META sequences are anchors that precede TARGET flanks
                         - left_flank and right_flank define the target boundaries
                         - marks column is informational (pattern uses .{0,5}? to match any mark)

            column_name: Column containing sequences to mask (default: 'r2_seq')
            fuzzy_pattern: Whether to use fuzzy matching for sequencing errors (default: True)
            fuzzy_kwargs: Optional dict with fuzzy matching parameters
            alias: Optional output column name. If None, replaces column_name in place

        Returns:
            LazyFrame with masked flank sequences

        Example:
            df = (
                pl.scan_parquet('parsed_reads.parquet')
                .pp.mask_flanks('features_edited.csv', column_name='r2_seq')
                .pp.assemble_umis(k=15, min_coverage=20)
            )

        See Also:
            ogtk.ltr.fracture.pipeline.masking.generate_mask_flanks_expression
        """
        mask_expr = generate_mask_flanks_expression(
            features_csv=features_csv,
            column_name=column_name,
            fuzzy_pattern=fuzzy_pattern,
            fuzzy_kwargs=fuzzy_kwargs,
            logger=self.logger
        )
        output_name = alias if alias is not None else column_name
        return self._ldf.with_columns(mask_expr.alias(output_name))

    def unmask_flanks(self,
                      features_csv: str,
                      column_name: str = 'contig',
                      fuzzy_pattern: bool = True,
                      fuzzy_kwargs: dict = None,
                      alias: str | None = None) -> pl.LazyFrame:
        """
        Restore original TARGET flanks from scrambled versions after assembly.

        This reverses the masking operation performed by mask_flanks(). After assembly,
        contigs contain scrambled flank sequences which need to be restored to their
        original form for downstream analysis. Lineage marks are preserved.

        Args:
            features_csv: Path to CSV file with META and EDITED_TARGET sequences
                         (same file used during masking)
            column_name: Column containing sequences to unmask (default: 'contig')
            fuzzy_pattern: Whether fuzzy matching was used during masking (default: True)
            fuzzy_kwargs: Optional dict with fuzzy matching parameters used during masking
            alias: Optional output column name. If None, replaces column_name in place

        Returns:
            LazyFrame with restored original flank sequences

        Example:
            df_contigs = (
                assembled_contigs
                .pp.unmask_flanks('features_edited.csv', column_name='contig')
            )

        See Also:
            ogtk.ltr.fracture.pipeline.masking.generate_unmask_flanks_expression
        """
        unmask_expr = generate_unmask_flanks_expression(
            features_csv=features_csv,
            column_name=column_name,
            fuzzy_pattern=fuzzy_pattern,
            fuzzy_kwargs=fuzzy_kwargs,
            logger=self.logger
        )
        output_name = alias if alias is not None else column_name
        return self._ldf.with_columns(unmask_expr.alias(output_name))

    def segment_by_metas(self,
                         metas: pl.DataFrame,
                         cassette_start_anchor: str,
                         cassette_end_anchor: str,
                         seq_col: str = 'r2_seq',
                         keep_cols: list | None = None) -> pl.LazyFrame:
        """
        Split reads at META boundaries into overlapping segments.

        Each segment includes both boundary METAs for overlap-based stitching.
        Segments are suitable for independent assembly.

        Args:
            metas: DataFrame with 'feature', 'seq', and optionally 'kind' columns
            cassette_start_anchor: Sequence marking cassette start (required)
            cassette_end_anchor: Sequence marking cassette end (required)
            seq_col: Name of the column containing sequences (default: 'r2_seq')
            keep_cols: Additional columns to preserve (default: ['umi'])

        Returns:
            LazyFrame with columns: umi, start_meta, end_meta, segment_seq
        """
        return segment_by_metas(
            ldf=self._ldf,
            metas=metas,
            cassette_start_anchor=cassette_start_anchor,
            cassette_end_anchor=cassette_end_anchor,
            seq_col=seq_col,
            keep_cols=keep_cols,
            logger=self.logger
        )

    def stitch_segments(self,
                        metas: pl.DataFrame,
                        seq_col: str = 'consensus_seq',
                        cassette_start_anchor_len: int | None = None,
                        group_cols: list[str] | None = None) -> pl.LazyFrame:
        """
        Rejoin assembled segments into full sequences.

        Uses META order to sort segments and trims overlapping META boundaries.

        Args:
            metas: DataFrame with 'feature' for ordering
            seq_col: Name of the column containing consensus sequences (default: 'consensus_seq')
            cassette_start_anchor_len: Length of cassette start anchor (for edge segment handling)
            group_cols: Columns to group by for stitching (default: ['umi'], can include 'sbc')

        Returns:
            LazyFrame with columns: [group_cols], stitched_seq
        """
        return stitch_segments(
            ldf=self._ldf,
            metas=metas,
            seq_col=seq_col,
            cassette_start_anchor_len=cassette_start_anchor_len,
            group_cols=group_cols,
            logger=self.logger
        )

    def flag_aberrant_molecules(self,
                                metas: pl.DataFrame,
                                seq_col: str = 'r2_seq') -> pl.LazyFrame:
        """
        Flag molecules with tandem duplications (any META appearing > 1 time).

        Useful for pre-filtering before segmentation to exclude aberrant reads.

        Args:
            metas: DataFrame with 'feature' and 'seq' columns
            seq_col: Name of the column containing sequences (default: 'r2_seq')

        Returns:
            LazyFrame with added 'is_aberrant' boolean column
        """
        return flag_aberrant_molecules(
            ldf=self._ldf,
            metas=metas,
            seq_col=seq_col
        )

    def assemble_segmented(self,
                           metas: pl.DataFrame,
                           cassette_start_anchor: str,
                           cassette_end_anchor: str,
                           k: int = 15,
                           min_coverage: int = 20,
                           seq_col: str = 'r2_seq',
                           debug_path: str | None = None,
                           method: str = 'compression',
                           output_dir: str | None = None,
                           sample_name: str | None = None,
                           heterogeneity_threshold: float = 0.20,
                           ) -> pl.LazyFrame:
        """
        Assemble reads using segmentation strategy for long cassettes.

        Segments reads at META boundaries, assembles each segment independently,
        then stitches results.

        Args:
            metas: DataFrame with META sequences ('feature', 'seq', 'kind' columns)
            cassette_start_anchor: Sequence marking cassette start (required)
            cassette_end_anchor: Sequence marking cassette end (required)
            k: K-mer size for assembly (default: 15)
            min_coverage: Minimum coverage for k-mers (default: 20)
            seq_col: Name of the column containing sequences (default: 'r2_seq')
            debug_path: If provided, save segments to this path before assembly
            method: Assembly method - 'compression' (default, recommended for segments)
                    or 'shortest_path'. Note: shortest_path often fails on short
                    segments (<300bp) because the de Bruijn graph may not have a
                    unique path between anchors.
            output_dir: If provided, generate segmentation QC plots in this directory
            sample_name: Sample name for plot titles (required if output_dir is provided)
            heterogeneity_threshold: For UMI heterogeneity detection, a transition type
                is considered "real" if its frequency is >= this fraction of the dominant
                type's frequency. Default 0.20 (20%). Set to 0 for sensitive mode only.

        Returns:
            LazyFrame with assembled contigs (umi, contig columns)
        """
        import rogtk

        self.logger.info("Using segmented assembly strategy")

        # Flag and filter aberrant molecules (tandem duplications)
        self.logger.debug("Flagging aberrant molecules")
        ldf = self._ldf.pp.flag_aberrant_molecules(metas=metas, seq_col=seq_col)
        n_aberrant = ldf.filter(pl.col('is_aberrant')).select(pl.col('umi').n_unique()).collect().item()
        self.logger.info(f"Found {n_aberrant} UMIs with tandem duplications (will be skipped)")
        ldf = ldf.filter(~pl.col('is_aberrant')).drop('is_aberrant')

        # Segment by METAs (sanitize + META→META segments + edge segments)
        # Include 'sbc' in keep_cols to prevent UMI collisions across samples
        self.logger.debug("Segmenting reads at META boundaries")
        keep_cols = ['umi', 'sbc'] if 'sbc' in ldf.collect_schema().names() else ['umi']
        segments_ldf = segment_by_metas(
            ldf=ldf,
            metas=metas,
            cassette_start_anchor=cassette_start_anchor,
            cassette_end_anchor=cassette_end_anchor,
            seq_col=seq_col,
            keep_cols=keep_cols,
            logger=self.logger,
        )

        # TODO remove: optional checkpoint to save segments before assembly
        if debug_path:
            self.logger.info(f"Saving segments to {debug_path}")
            segments_ldf.collect().write_parquet(debug_path)
            segments_ldf = pl.scan_parquet(debug_path)
        # END TODO remove

        # Log segment type count
        segment_types = segments_ldf.select('start_meta', 'end_meta').unique().collect()
        self.logger.info(f"Found {len(segment_types)} unique segment types to assemble")

        # Assemble each segment type
        self.logger.debug(f"Assembling segments with method={method}")
        # Group by sbc+umi to prevent collisions across samples
        group_cols = ['sbc', 'umi', 'start_meta', 'end_meta'] if 'sbc' in keep_cols else ['umi', 'start_meta', 'end_meta']

        if method == 'shortest_path':
            # For shortest_path: assemble each segment type with appropriate anchors
            # Build META sequence lookup (including cassette anchors as pseudo-METAs)
            meta_seq_map = dict(zip(metas['feature'], metas['seq']))
            meta_seq_map[CASSETTE_START_MARKER] = cassette_start_anchor
            meta_seq_map[CASSETTE_END_MARKER] = cassette_end_anchor

            # Get unique segment types
            segment_types_df = segments_ldf.select('start_meta', 'end_meta').unique().collect()
            assembled_parts = []

            for row in segment_types_df.iter_rows(named=True):
                start_meta = row['start_meta']
                end_meta = row['end_meta']

                # Derive anchors from boundary METAs (use k bp from each end)
                start_seq = meta_seq_map.get(start_meta, '')
                end_seq = meta_seq_map.get(end_meta, '')
                start_anchor = start_seq[-k:] if len(start_seq) >= k else start_seq
                end_anchor = end_seq[:k] if len(end_seq) >= k else end_seq

                if not start_anchor or not end_anchor:
                    self.logger.warning(f"Skipping {start_meta}->{end_meta}: missing anchor sequences")
                    continue

                self.logger.debug(f"Assembling {start_meta}->{end_meta} with anchors: {start_anchor[:10]}.../{end_anchor[:10]}...")

                # Filter to this segment type, trim to anchors, and assemble
                segment_assembled = (
                    segments_ldf
                    .filter((pl.col('start_meta') == start_meta) & (pl.col('end_meta') == end_meta))
                    # Trim sequences to anchor boundaries
                    .with_columns(
                        pl.col('segment_seq')
                        .str.replace(f'^.*?({start_anchor})', f'{start_anchor}')
                        .str.replace(f'({end_anchor}).*$', f'{end_anchor}')
                        .alias('segment_seq_trimmed')
                    )
                    .group_by(group_cols)
                    .agg(
                        rogtk.assemble_sequences(
                            expr=pl.col('segment_seq_trimmed'),
                            k=k,
                            min_coverage=min_coverage,
                            method='shortest_path',
                            start_anchor=start_anchor,
                            end_anchor=end_anchor,
                        ).alias('consensus_seq')
                    )
                    .filter(pl.col('consensus_seq').str.len_chars() > 0)
                )
                assembled_parts.append(segment_assembled)

            # Combine all assembled segments
            if assembled_parts:
                assembled = pl.concat(assembled_parts)
            else:
                self.logger.error("No segments could be assembled")
                assembled = pl.LazyFrame(schema={c: pl.Utf8 for c in group_cols} | {'consensus_seq': pl.Utf8})
        else:
            # For compression: single pass assembly (original behavior)
            assembled = (
                segments_ldf
                .group_by(group_cols)
                .agg(
                    rogtk.assemble_sequences(
                        expr=pl.col('segment_seq'),
                        k=k,
                        min_coverage=min_coverage,
                        method=method,
                    ).alias('consensus_seq')
                )
                # Filter out failed assemblies before stitching
                .filter(pl.col('consensus_seq').str.len_chars() > 0)
            )

        # Collect intermediate data for reporting
        self.logger.debug("Collecting data for report generation")
        segments_df = segments_ldf.collect()
        assembled_df = assembled.collect()

        # Validate assemblies: check that consensus contains expected anchors
        # Build anchor lookup for validation
        meta_seq_map = dict(zip(metas['feature'], metas['seq']))
        meta_seq_map[CASSETTE_START_MARKER] = cassette_start_anchor
        meta_seq_map[CASSETTE_END_MARKER] = cassette_end_anchor

        # Add validation columns: does consensus contain start/end anchors?
        validation_exprs = []
        for col_name in ['start_meta', 'end_meta']:
            # Build a when/then chain for each META
            anchor_check = pl.lit(False)
            for meta_name, meta_seq in meta_seq_map.items():
                if meta_seq:
                    # Use first/last k bp as the anchor to check
                    check_seq = meta_seq[-k:] if col_name == 'start_meta' else meta_seq[:k]
                    anchor_check = (
                        pl.when(pl.col(col_name) == meta_name)
                        .then(pl.col('consensus_seq').str.contains(pl.lit(check_seq), literal=True))
                        .otherwise(anchor_check)
                    )
            validation_exprs.append(anchor_check.alias(f'_has_{col_name}_anchor'))

        assembled_df = assembled_df.with_columns(validation_exprs)
        assembled_df = assembled_df.with_columns(
            (pl.col('_has_start_meta_anchor') & pl.col('_has_end_meta_anchor')).alias('_is_valid_assembly')
        )

        # Count invalid assemblies
        n_invalid = assembled_df.filter(~pl.col('_is_valid_assembly')).height
        if n_invalid > 0:
            self.logger.warning(f"Found {n_invalid} truncated assemblies (missing anchors)")

        # Fallback for failed/invalid assemblies: use most common raw segment
        # Identify UMI groups that have segments but no valid assembly
        segment_groups = segments_df.select(group_cols).unique()
        valid_assembled_groups = assembled_df.filter(pl.col('_is_valid_assembly')).select(group_cols).unique()
        missing_groups = segment_groups.join(valid_assembled_groups, on=group_cols, how='anti')

        if missing_groups.height > 0:
            self.logger.info(f"Falling back to raw segments for {missing_groups.height} failed/invalid assemblies")

            # For each missing group, select the most common segment sequence
            # If tied, the first one (effectively random based on data order) is selected
            fallback_df = (
                segments_df
                .join(missing_groups, on=group_cols, how='inner')
                .group_by(group_cols + ['segment_seq'])
                .agg(pl.len().alias('_seq_count'))
                .sort('_seq_count', descending=True)
                .group_by(group_cols)
                .first()  # Take most common (or first if tied)
                .rename({'segment_seq': 'consensus_seq'})
                .with_columns([
                    pl.lit(True).alias('_is_fallback'),
                    pl.lit(True).alias('_has_start_meta_anchor'),  # Raw segments have anchors by construction
                    pl.lit(True).alias('_has_end_meta_anchor'),
                    pl.lit(True).alias('_is_valid_assembly'),
                ])
                .drop('_seq_count')
            )

            # Remove invalid assemblies and merge with fallbacks
            assembled_df = (
                assembled_df
                .filter(pl.col('_is_valid_assembly'))
                .with_columns(pl.lit(False).alias('_is_fallback'))
            )
            assembled_df = pl.concat([assembled_df, fallback_df], how='diagonal')
            self.logger.info(f"Total assembled after fallback: {assembled_df.height}")
        else:
            # No fallbacks needed, just add the column
            assembled_df = assembled_df.with_columns(pl.lit(False).alias('_is_fallback'))

        # Add segment count per molecule (how many segments does this UMI have?)
        molecule_cols = ['sbc', 'umi'] if 'sbc' in group_cols else ['umi']
        assembled_df = assembled_df.with_columns(
            pl.len().over(molecule_cols).alias('_n_segments')
        )

        # Extract TARGET insertions from assembled consensus sequences
        # Done after assembly to use clean consensus (removes sequencing noise)
        #
        # Strategy: Use shorter anchor portions (last/first N chars of flanks) to handle
        # indels that may delete part of the flank.
        if 'left_flank' in metas.columns and 'right_flank' in metas.columns:
            self.logger.debug("Extracting TARGET insertions from assembled segments")
            targets_df = metas.filter(pl.col('kind') == 'EDITED_TARGET')
            anchor_len = 8  # Use last/first 8bp of flanks as anchors
            max_insertion = 6  # Max insertion length to capture (PE marks are 5bp)

            for row in targets_df.iter_rows(named=True):
                target_name = row['feature']
                left_flank = row['left_flank']
                right_flank = row['right_flank']
                if left_flank and right_flank:
                    # Use partial anchors: last N chars of left_flank, first N chars of right_flank
                    # This allows detection even when part of the flank is deleted at cut site
                    left_anchor = left_flank[-anchor_len:] if len(left_flank) > anchor_len else left_flank
                    right_anchor = right_flank[:anchor_len] if len(right_flank) > anchor_len else right_flank

                    # Pattern: left_anchor + (0 to max_insertion bp) + right_anchor
                    # The captured group contains the insertion (empty string = wild-type)
                    pattern = f'{left_anchor}(.{{0,{max_insertion}}}){right_anchor}'
                    col_name = f'{target_name}_ins'
                    assembled_df = assembled_df.with_columns(
                        pl.col('consensus_seq').str.extract(pattern, 1).alias(col_name)
                    )

        # Save assembled segments with insertions
        if debug_path:
            from pathlib import Path
            debug_path_obj = Path(debug_path)
            assembled_path = debug_path_obj.with_stem(debug_path_obj.stem.replace('segments', 'assembled'))
            self.logger.info(f"Saving assembled segments with insertions to {assembled_path}")
            assembled_df.write_parquet(assembled_path)

        # Stitch segments back together
        self.logger.debug("Stitching assembled segments")
        cassette_start_anchor_len = len(cassette_start_anchor) if cassette_start_anchor else None
        # Use same grouping columns (with sbc if present) for stitching
        stitch_group_cols = ['sbc', 'umi'] if 'sbc' in keep_cols else ['umi']
        stitched = assembled_df.lazy().pp.stitch_segments(
            metas=metas,
            seq_col='consensus_seq',
            cassette_start_anchor_len=cassette_start_anchor_len,
            group_cols=stitch_group_cols,
        )

        # Collect contigs for reporting
        contigs_df = stitched.rename({'stitched_seq': 'contig'}).collect()

        # Generate and log diagnostic report
        self.logger.debug("Generating segmentation report")
        report = generate_segmentation_report(
            segments_df=segments_df,
            assembled_df=assembled_df,
            contigs_df=contigs_df,
            metas=metas,
            cassette_start_anchor=cassette_start_anchor,
            cassette_end_anchor=cassette_end_anchor,
            heterogeneity_threshold=heterogeneity_threshold,
        )
        report_str = format_segmentation_report(report)
        self.logger.info(f"\n{report_str}")

        # Generate QC plots if output directory provided
        if output_dir is not None:
            if sample_name is None:
                sample_name = "sample"
            from pathlib import Path
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            plot_segmentation_qc(
                report=report,
                contigs_df=contigs_df,
                metas=metas,
                output_dir=output_path,
                sample_name=sample_name,
                contig_col='contig',
                logger=self.logger,
            )
            self.logger.info(f"Saved segmentation QC plots to {output_path}")

        # Return as LazyFrame for consistency
        return contigs_df.lazy()


# TODO
# ppp is the temporary expansion of the routines in a chunked manner with split logics
# under development
@pl.api.register_dataframe_namespace("ppp")
class Chunked:
    """Provides chunked processing methods for large datasets that don't fit in memory.
    
    Methods are registered under the 'ppp' namespace and can be accessed via df.ppp.*
    This is designed for DataFrames that are already in memory but need chunked operations.
    """
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df
        self.logger = Rlogger().get_logger()
    
    def parse_reads_chunked(self, umi_len: int, anchor_ont: str, sbc_len: int, 
                           chunk_size: int = 100_000_000) -> pl.DataFrame:
        """Parse reads in chunks to avoid memory issues.
        
        Args:
            umi_len (int): Length of UMI sequence
            anchor_ont (str): Anchor sequence to search for
            sbc_len (int): Sample barcode length  
            chunk_size (int): Maximum rows per chunk (default 100M)
            
        Returns:
            pl.DataFrame: Processed dataframe with all parse_reads columns
        """
        total_rows = self._df.height
        
        if total_rows <= chunk_size:
            self.logger.debug(f"Dataset has {total_rows:,} rows, processing without chunking")
            return self._parse_reads_single_chunk(self._df, umi_len, anchor_ont, sbc_len)
        
        self.logger.info(f"Processing {total_rows:,} rows in chunks of {chunk_size:,}")
        
        chunks = []
        for start_idx in range(0, total_rows, chunk_size):
            end_idx = min(start_idx + chunk_size, total_rows)
            chunk_rows = end_idx - start_idx
            
            self.logger.debug(f"Processing chunk {start_idx:,} to {end_idx:,} ({chunk_rows:,} rows)")
            
            chunk_df = self._df.slice(start_idx, chunk_rows)
            processed_chunk = self._parse_reads_without_globals(chunk_df, umi_len, anchor_ont, sbc_len)
            chunks.append(processed_chunk)
        
        # Concatenate all chunks
        self.logger.debug("Concatenating chunks and computing global metrics")
        full_df = pl.concat(chunks)
        
        # Add global metrics and window functions
        return self._add_global_metrics_and_windows(full_df)
    
    def _parse_reads_single_chunk(self, df: pl.DataFrame, umi_len: int, anchor_ont: str, sbc_len: int) -> pl.DataFrame:
        """Process a single chunk that fits in memory using the original logic."""
        return (
            df
            .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
            .with_columns(pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc'))
            .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
            # read level metrics
            .with_columns(metric_reads_ont = pl.col('ont').sum())
            .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
            .with_columns(metric_fraction_ont = pl.col('ont').mean())
            .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
            # UMI validation
            .with_columns(pl.col('r1_seq')
                          .str.extract(f'(.*?){anchor_ont}.*$',1)
                          .str.len_chars()
                          .fill_null(0)
                          .alias('valid_umi')
                          )
            .with_columns((pl.col('valid_umi')==umi_len+sbc_len).alias('valid_umi'))
            .with_columns(pl.len().over('umi').alias('reads'))
        )
    
    def _parse_reads_without_globals(self, df: pl.DataFrame, umi_len: int, anchor_ont: str, sbc_len: int) -> pl.DataFrame:
        """Process reads without global aggregations or window functions."""
        return (
            df
            .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
            .with_columns(pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc'))
            .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
            # UMI validation
            .with_columns(pl.col('r1_seq')
                          .str.extract(f'(.*?){anchor_ont}.*$',1)
                          .str.len_chars()
                          .fill_null(0)
                          .alias('valid_umi')
                          )
            .with_columns((pl.col('valid_umi')==umi_len+sbc_len).alias('valid_umi'))
        )
    
    def _add_global_metrics_and_windows(self, df: pl.DataFrame) -> pl.DataFrame:
        """Add global metrics and window functions to the concatenated dataframe."""
        return (
            df
            # Global read-level metrics
            .with_columns(metric_reads_ont = pl.col('ont').sum())
            .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
            .with_columns(metric_fraction_ont = pl.col('ont').mean())
            .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
            # Window function for reads per UMI
            .with_columns(pl.len().over('umi').alias('reads'))
        )


@pl.api.register_lazyframe_namespace("ppp")
class LazyChunked:
    """Provides chunked processing methods for LazyFrames to handle large parquet files.
    
    Methods are registered under the 'ppp' namespace and can be accessed via ldf.ppp.*
    This processes data in chunks directly from parquet files without loading everything into memory.
    """
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf
        self.logger = Rlogger().get_logger()
    
    def parse_reads_chunked(self, umi_len: int, anchor_ont: str, sbc_len: int, 
                           chunk_size: int = 100_000_000, 
                           total_rows: int = None) -> pl.LazyFrame:
        """Parse reads in chunks from a LazyFrame (typically from parquet scan).
        
        Args:
            umi_len (int): Length of UMI sequence
            anchor_ont (str): Anchor sequence to search for
            sbc_len (int): Sample barcode length
            chunk_size (int): Maximum rows per chunk (default 100M)
            total_rows (int): Total number of rows (if known, avoids counting)
            
        Returns:
            pl.LazyFrame: Processed lazy dataframe ready for sinking
        """
        # Get total row count if not provided
        if total_rows is None:
            self.logger.debug("Counting total rows in dataset")
            total_rows = self._ldf.select(pl.len()).collect().item()
        
        if total_rows <= chunk_size:
            self.logger.debug(f"Dataset has {total_rows:,} rows, processing without chunking")
            return self._parse_reads_single_lazy(self._ldf, umi_len, anchor_ont, sbc_len)
        
        self.logger.info(f"Processing {total_rows:,} rows in chunks of {chunk_size:,}")
        
        # Process chunks and collect them
        chunk_dfs = []
        for start_idx in range(0, total_rows, chunk_size):
            end_idx = min(start_idx + chunk_size, total_rows)
            chunk_rows = end_idx - start_idx
            
            self.logger.debug(f"Processing chunk {start_idx:,} to {end_idx:,} ({chunk_rows:,} rows)")
            
            # Slice the lazy frame and process the chunk
            chunk_ldf = self._ldf.slice(start_idx, chunk_rows)
            processed_chunk = self._parse_reads_without_globals_lazy(chunk_ldf, umi_len, anchor_ont, sbc_len)
            
            # Collect this chunk to avoid memory issues with many lazy frames
            chunk_df = processed_chunk.collect()
            chunk_dfs.append(chunk_df)
        
        # Concatenate all chunks and create a new lazy frame
        self.logger.debug("Concatenating chunks and preparing final lazy frame")
        full_df = pl.concat(chunk_dfs)
        
        # Convert back to lazy frame and add global metrics and window functions
        return self._add_global_metrics_and_windows_lazy(pl.LazyFrame(full_df))
    
    def _parse_reads_single_lazy(self, ldf: pl.LazyFrame, umi_len: int, anchor_ont: str, sbc_len: int) -> pl.LazyFrame:
        """Process a single lazy frame that should fit in memory."""
        return (
            ldf
            .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
            .with_columns(pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc'))
            .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
            # read level metrics
            .with_columns(metric_reads_ont = pl.col('ont').sum())
            .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
            .with_columns(metric_fraction_ont = pl.col('ont').mean())
            .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
            # UMI validation
            .with_columns(pl.col('r1_seq')
                          .str.extract(f'(.*?){anchor_ont}.*$',1)
                          .str.len_chars()
                          .fill_null(0)
                          .alias('valid_umi')
                          )
            .with_columns((pl.col('valid_umi')==umi_len+sbc_len).alias('valid_umi'))
            .with_columns(pl.len().over('umi').alias('reads'))
        )
    
    def _parse_reads_without_globals_lazy(self, ldf: pl.LazyFrame, umi_len: int, anchor_ont: str, sbc_len: int) -> pl.LazyFrame:
        """Process reads without global aggregations or window functions."""
        return (
            ldf
            .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
            .with_columns(pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc'))
            .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
            # UMI validation
            .with_columns(pl.col('r1_seq')
                          .str.extract(f'(.*?){anchor_ont}.*$',1)
                          .str.len_chars()
                          .fill_null(0)
                          .alias('valid_umi')
                          )
            .with_columns((pl.col('valid_umi')==umi_len+sbc_len).alias('valid_umi'))
        )
    
    def _add_global_metrics_and_windows_lazy(self, ldf: pl.LazyFrame) -> pl.LazyFrame:
        """Add global metrics and window functions to the lazy dataframe."""
        return (
            ldf
            # Global read-level metrics  
            .with_columns(metric_reads_ont = pl.col('ont').sum())
            .with_columns(metric_reads_offt = pl.col('ont').not_().sum())
            .with_columns(metric_fraction_ont = pl.col('ont').mean())
            .with_columns(metric_fraction_offt = pl.col('ont').not_().mean())
            # Window function for reads per UMI
            .with_columns(pl.len().over('umi').alias('reads'))
        )
    
    def sink_parquet_chunked(self, path: str, chunk_size: int = 100_000_000) -> None:
        """Sink a large LazyFrame to parquet in chunks to avoid memory issues.
        
        This is useful when even sinking fails due to memory constraints.
        
        Args:
            path (str): Output parquet file path
            chunk_size (int): Maximum rows per chunk when sinking
        """
        self.logger.info(f"Sinking LazyFrame to {path} in chunks of {chunk_size:,}")
        
        # Get total row count
        total_rows = self._ldf.select(pl.len()).collect().item()
        
        if total_rows <= chunk_size:
            self.logger.debug("Dataset fits in single chunk, using regular sink")
            self._ldf.sink_parquet(path)
            return
        
        # Process and write chunks
        first_chunk = True
        for start_idx in range(0, total_rows, chunk_size):
            end_idx = min(start_idx + chunk_size, total_rows)
            chunk_rows = end_idx - start_idx
            
            self.logger.debug(f"Sinking chunk {start_idx:,} to {end_idx:,} ({chunk_rows:,} rows)")
            
            chunk_ldf = self._ldf.slice(start_idx, chunk_rows)
            
            if first_chunk:
                # First chunk creates the file
                chunk_ldf.sink_parquet(path)
                first_chunk = False
            else:
                # Subsequent chunks append (this requires collecting and re-writing)
                # Note: Polars doesn't support true append mode for parquet
                # So we need to collect previous data and append
                existing_df = pl.scan_parquet(path).collect()
                chunk_df = chunk_ldf.collect()
                combined_df = pl.concat([existing_df, chunk_df])
                combined_df.write_parquet(path)
        
        self.logger.info(f"Successfully wrote {total_rows:,} rows to {path}")
