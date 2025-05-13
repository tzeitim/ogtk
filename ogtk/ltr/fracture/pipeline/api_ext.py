import polars as pl
import rogtk

from ogtk.utils.log import Rlogger, call

__all__ = [
        'PlDNA',
        'PlPipeline'
]

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
@pl.api.register_lazyframe_namespace("pp")
class PllPipeline:
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf

    # [12N][6N]
    # [  18N  ]
    def parse_reads(self, umi_len, anchor_ont, sbc_len):
        ''' sbc_len : sample barcode length'''
        # _field_name represent QC columns
        return(
                self._ldf
                    .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
                    .with_columns(pl.col('r1_seq').str.slice(umi_len, sbc_len).alias('sbc'))
                    .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
                    .with_columns(_reads_ont = pl.col('ont').sum())
                    .with_columns(_reads_offt = pl.col('ont').not_().sum())
                    .with_columns(_fraction_ont = pl.col('ont').mean())
                    .with_columns(_fraction_offt = pl.col('ont').not_().mean())
                    .filter(pl.col('ont'))
                    # UMI validation
                    .with_columns(pl.col('r1_seq')
                                  .str.extract(f'(.*?){anchor_ont}.*$',1)
                                  .str.len_chars()
                                  .alias('valid_umi')
                                  )
                    .with_columns((pl.col('valid_umi')==umi_len+sbc_len))
                    .with_columns(_valid_umis = pl.col('valid_umi').mean())
                    .with_columns(_n_valid_umis = pl.col('valid_umi').sum())
                    .with_columns(_n_invalid_umis = pl.col('valid_umi').not_().sum())
                    .filter(pl.col('valid_umi'))
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
                    .with_columns(pl.col('r1_seq').str.extract(f'.*?{anchor_ont}(.*)$').alias('r2_seq'))
                )
