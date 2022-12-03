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


def lalla():
    print("lalal")
    return    ib.load_wl()

def shgrna(sample_id: str,
           tbxifn: str, 
           valid_ibars: Sequence, 
           min_reads: int=1, 
           max_reads: int=1e6,
           clone: str | None=None,
           ) -> pl.DataFrame:
    ''' Analytical routine to process UMIfied shRNA from bulk assays.
        If `clone` is not provided it assumes that can be found as the first element of a dash-splitted `sample_id` string. 
        Off-targets are defined by reads without an ibar pattern match.
        UMIs are filtered based on the number of reads
    '''
    df = pl.read_csv(tbxifn, sep='\t', has_header=False)
    df.columns=['readid',  'start' ,'end'  , 'cbc' , 'umi' , 'seq' , 'qual']
    
    rdf = ib.extract_read_grammar_new(sample_id, df=df.drop('cbc'))
    tot_reads = rdf.shape[0]
    rdf = ib.ibar_reads_to_molecules(rdf, modality='single-molecule')
    rdf = ib.count_essential_patterns(rdf)

    # guess from the first field of the sample_id the clone of origin
    # if `clone` is not specified
    if clone is None:
        clone = sample_id.split('_')[0]

    ib.noise_spectrum(sample_id, rdf, index = ['dss'], columns='tsos')    

    offtarget = rdf.filter(pl.col("raw_ibar").is_null()).shape[0]/rdf.shape[0]
    tot_umis = rdf.select('umi').n_unique()

    print(f'{offtarget=:.2%}')
    print(f'{tot_umis=}')
    print(f'{tot_reads=}')

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
           .with_columns(10000*(pl.col('qc_tot_umis')/pl.col('qc_tot_reads')).alias('qc_umis_per_10kreads'))

          )
