import polars as pl
# rogtk includes the "dna" namespace
import rogtk

from ogtk.utils.log import Rlogger, call

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
                    auto_k: bool = True,
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
                    auto_k: bool = True,
                    export_graphs: bool = False,
                    groups = ['umi', 'sbc'],
                    only_largest: bool = True) -> pl.DataFrame:
        """
        """
        return (
            self._df
             .with_columns(pl.col('r2_seq').str.replace(f'^.*{start_anchor}', start_anchor))
             .with_columns(pl.col('r2_seq').str.replace(f'{end_anchor}.*$', end_anchor))
            .group_by(groups).agg(
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
                ).alias('contig')
            )
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

    def ont_to_paired_format(self, umi_len: int, sbc_len: int, anchor: str) -> pl.DataFrame:
        """Convert ONT BAM data to FRACTURE's expected R1/R2 format (DataFrame wrapper)"""
        return (
            self._df
            .lazy()
            .pp.ont_to_paired_format(umi_len=umi_len, sbc_len=sbc_len, anchor=anchor)
            .collect()
        )

@pl.api.register_lazyframe_namespace("pp")
class PllPipeline:
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf

    # [12N][6N]
    # [  18N  ]
    def parse_reads(self, umi_len, anchor_ont, sbc_len):
        ''' sbc_len : sample barcode length'''
        # metric_field_name represent QC columns
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
                    #.filter(pl.col('ont'))
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


    def ont_to_paired_format(self, umi_len: int, sbc_len: int, anchor: str) -> pl.LazyFrame:
        """Convert ONT BAM data to FRACTURE's expected R1/R2 format (LazyFrame)"""
        return (
            self._ldf
            .with_columns([
                pl.when(pl.col('sequence').str.contains(anchor))
                .then(pl.col('sequence'))
                .when(pl.col('sequence').dna.reverse_complement().str.contains(anchor))
                .then(pl.col('sequence').dna.reverse_complement())
                .otherwise(pl.col('sequence'))
                .alias('oriented_sequence'),
                
                pl.when(pl.col('sequence').str.contains(anchor))
                .then(pl.col('quality_scores'))
                .when(pl.col('sequence').dna.reverse_complement().str.contains(anchor))
                .then(pl.col('quality_scores').str.reverse())
                .otherwise(pl.col('quality_scores'))
                .alias('oriented_quality')
            ])
            .with_columns([
                pl.col('oriented_sequence').str.slice(0, umi_len + sbc_len).alias('r1_seq'),
                pl.col('oriented_sequence').str.slice(umi_len + sbc_len).alias('r2_seq'),
                pl.col('oriented_quality').str.slice(0, umi_len + sbc_len).alias('r1_qual'),
                pl.col('oriented_quality').str.slice(umi_len + sbc_len).alias('r2_qual')
            ])
            .with_columns([
                pl.col('r1_seq').str.slice(0, umi_len).alias('umi'),
                pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc')
            ])
            .select(['name', 'r1_seq', 'r1_qual', 'r2_seq', 'r2_qual', 'umi', 'sbc'])
        )

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
