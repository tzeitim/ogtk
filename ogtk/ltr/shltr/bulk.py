import polars as pl

def qc_ibars_per_sample(df: pl.DataFrame)-> pl.DataFrame:
    '''
        compute molecule stats at the level of samples and integration
            - total reads per ibar (`ibar_reads`)
            - total umis per ibar (`ibar_umis`)
            - log10 total reads per ibar (`ibar_reads_log10`)
            - log10 total umis per ibar (`ibar_umis_log10`)
        Computes at the whole df level the quantiled counts
            - -log10 quantile conversion total reads per ibar (`ibar_reads_q`)
            - -log10 quantile conversion total umis per ibar (`ibar_umis_q`)
    '''
    groups = ['sample_id', 'raw_ibar']
    df = (df
         .with_columns([
             (pl.col('umi_dom_reads').sum().over(groups)).alias('ibar_reads'),
             (pl.col('umi').n_unique().over(groups)).alias('ibar_umis')])
         .with_columns(
             [pl.col('ibar_reads').log10().alias('ibar_reads_log10'), 
              pl.col('ibar_umis').log10().alias('ibar_umis_log10')])
        .with_columns((-1 * (pl.col('ibar_reads')/(pl.col('ibar_reads').sum())).log10()).alias('ibar_reads_q'))
        .with_columns((-1 * (pl.col('ibar_umis')/(pl.col('ibar_umis').sum())).log10()).alias('ibar_umis_q'))
         )
    return df
