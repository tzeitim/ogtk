import polars as pl
import rogtk


@pl.api.register_dataframe_namespace("dna")
class PlDNA:
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df

    def to_fasta(self, read_id_col: str, read_col: str) -> pl.DataFrame:
        return  self._df.with_columns(
                (">"+pl.col(read_id_col)\
                 +"\n"+pl.col(read_col)
                 )
                 .alias(f'{read_col}_fasta')
        )
    def to_fastq(self, read_id_col: str, read_qual_col: str, read_col: str)-> pl.DataFrame:
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
   
    def assembly_with_opt(self, start_k=25, start_min_coverage=17, start_anchor="CCGATCT", end_anchor="TTTAGTGAGGGT"):
        return (
                self._df
                .filter(pl.col('reads')>100)
                .group_by(['umi']).agg(
                      rogtk.optimize_assembly(
                      expr=pl.col('r2_seq'), 
                      start_k=start_k, 
                      start_min_coverage=start_min_coverage, 
                      start_anchor=start_anchor, 
                      end_anchor=end_anchor,
                      max_iterations=250))
                .unnest('r2_seq')
                .with_columns((pl.col('length')==0).alias('failed'))
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

