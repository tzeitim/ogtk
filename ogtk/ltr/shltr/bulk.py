from logging import warning
from typing import Sequence,Optional
from sys import prefix
import pysam
import regex
import ogtk 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl
from . import ibars as ib
#from . import plot_bulk as pt

def lalla():
    print("lalal")
    return    ib.load_wl()

def shgrna(sample_id: str,
           valid_ibars: Sequence, 
           parquet_ifn: str|None=None, 
           tbxifn: str|None=None, 
           min_reads: int=int(1), 
           max_reads: int=int(1e6),
           clone: str | None=None,
           downsample: None=None,
           ) -> pl.DataFrame:
    ''' Analytical routine to process UMIfied shRNA from bulk assays.
        If `clone` is not provided it assumes that can be found as the first element of a dash-splitted `sample_id` string. 
        Off-targets are defined by reads without an ibar pattern match.
        UMIs are filtered based on the number of reads
    '''
    assert not ((tbxifn is None) and (parquet_ifn is None) ), "You must provide a tabix or parquet input filename"
    assert not ((tbxifn is not None) and (parquet_ifn is not None) ), "You must provide either a tabix or parquet input filename, not both"

    if parquet_ifn is not None:
        rdf = ib.extract_read_grammar(parquet_ifn = parquet_ifn, batch = sample_id, sample = downsample)

    if tbxifn is not None:
        df = pl.read_csv(tbxifn, separator='\t', has_header=False)
        df.columns=['readid',  'start' ,'end'  , 'cbc' , 'umi' , 'seq' , 'qual']
        rdf = ib.extract_read_grammar(batch = sample_id, df=df.drop('cbc'), sample= downsample)

    tot_reads = rdf.shape[0]
    rdf = ib.ibar_reads_to_molecules(rdf, modality='single-molecule')
    rdf = ib.count_essential_patterns(rdf)
    rdf = ib.noise_properties(rdf)

    # guess from the first field of the sample_id the clone of origin
    # if `clone` is not specified
    if clone is None:
        clone = sample_id.split('_')[0]

    # QC plots
    ib.noise_spectrum(sample_id, rdf, index = ['dss'], columns='tsos')    
    ib.plot_noise_properties(rdf)

    offtarget = rdf.filter(pl.col("raw_ibar").is_null()).shape[0]/rdf.shape[0]
    tot_umis = rdf.select('umi').n_unique()

    print(f'{offtarget=:.2%}')
    print(f'{tot_umis=}')
    print(f'{tot_reads=}')

    print(f'umis per 10k reads {10000*tot_umis/tot_reads}')

    return(rdf
           .filter(pl.col('raw_ibar').is_not_null())
           .filter(pl.col('raw_ibar').is_in(valid_ibars))
           .filter(pl.col('umi_dom_reads')>=min_reads)
           .filter(pl.col('umi_dom_reads')<=max_reads)
           .with_columns(pl.lit(sample_id).alias('sample_id'))
           .with_columns(pl.lit(clone).alias('clone'))
           .with_columns(pl.lit(offtarget).alias('qc_pc_offt'))
           .with_columns(pl.lit(tot_umis).alias('qc_tot_umis'))
           .with_columns(pl.lit(tot_reads).alias('qc_tot_reads'))
           .with_columns(pl.col(pl.Int32).cast(pl.Int64))
           .with_columns((10000.0*pl.col('qc_tot_umis')/pl.col('qc_tot_reads')).alias('qc_umis_per_10kreads'))
        
          )

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
