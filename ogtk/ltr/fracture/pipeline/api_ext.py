import polars as pl
import rogtk


@pl.api.register_dataframe_namespace("dna")
class PlDNA:
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df

    def to_fasta(self,
         read_id_col: str,
         read_col: str) -> pl.DataFrame:
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
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df
   
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
        
        Args:
            start_k: Initial k-mer size for graph construction
            start_min_coverage: Initial minimum k-mer coverage threshold
            method: Assembly method ('compression' or 'shortest_path')
            start_anchor: Start sequence anchor for shortest_path method
            end_anchor: End sequence anchor for shortest_path method
            min_length: Minimum contig length to return
            export_graphs: Whether to export graph visualization files
            prioritize_length: Optimize for length over anchor presence
            min_reads: Minimum number of reads to attempt assembly
            
        Returns:
            DataFrame with assembly results including contig sequences and parameters
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
        
        Args:
            target_umi: UMI sequence to assemble
            k: K-mer size for graph construction (used if auto_k=False)
            min_coverage: Minimum k-mer coverage threshold
            method: Assembly method ('compression' or 'shortest_path')
            start_anchor: Start sequence anchor for shortest_path method
            end_anchor: End sequence anchor for shortest_path method
            min_length: Minimum contig length to return
            auto_k: Automatically estimate optimal k-mer size
            export_graphs: Whether to export graph visualization files
            only_largest: Return only the largest contig
            prefix: Prefix for output files
            
        Returns:
            DataFrame with assembly results per UMI
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
        
        Args:
            target_umi: UMI sequence to assemble
            k_start: Starting k-mer size
            k_end: Ending k-mer size
            k_step: Step size for k-mer iteration
            cov_start: Starting minimum coverage
            cov_end: Ending minimum coverage  
            cov_step: Step size for coverage iteration
            method: Assembly method ('compression' or 'shortest_path')
            start_anchor: Start sequence anchor for shortest_path method
            end_anchor: End sequence anchor for shortest_path method
            min_length: Minimum contig length to return
            export_graphs: Whether to export graph visualization files
            auto_k: Automatically estimate optimal k-mer size
            prefix: Prefix for output files
            
        Returns:
            DataFrame with assembly results for each parameter combination
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

    def parse_reads(self, umi_len, anchor_ont, intbc_5prime):
        return(
                self._ldf
                    .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
                    .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
                    .filter(pl.col('ont'))
                    .with_columns(pl.len().over('umi').alias('reads'))
                    # strip out the UMI from R1
                    .with_columns(pl.col('r1_seq').str.replace(f"^.+?{anchor_ont}", anchor_ont))
                    # strip out the sequences up to intbc_5prime 
        #            .with_columns(pl.col('r2_seq').str.replace(f'^.+?{intbc_5prime}', intbc_5prime))
                )    

