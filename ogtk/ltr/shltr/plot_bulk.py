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

def umi_read_cov(df: pl.DataFrame,
			 min_reads: int=0,
			 max_reads: int=int(1e6),
			 n_downsample: int=1000,
             over: Sequence|str='sample_id',
                 )-> None:
    ''' Plots a log-log histogram of the reads per umi of a mol-level data frame.

        min_reads: minimum number of reads for a UMI to be considered
        max_reads: top limit of number of reads to be considered
        n_downsample: number of entries (molecules) for the group for which data is aggregate on (`.over`)
    '''
    fg = sns.displot(data=(df
                      .filter(pl.arange(0, pl.count()).shuffle().over(over) < n_downsample)
                      .sort(over)
                      .filter((pl.col('umi_dom_reads')<=max_reads) & (pl.col('umi_dom_reads')>=min_reads))
                      .to_pandas()), 
                x='umi_reads', 
                kind='hist', 
                element='step',
                aspect=2,
                hue='sample_id',
                col='sample_id',
                col_wrap=4,
                log_scale=(10, 10),
                facet_kws={'sharey':True})

        
    if len(fg.axes_dict)>0:
        for ax in fg.axes_dict.values():
            ax.grid()
    else:
        ax = fg.ax
        ax.grid()
        ax.set_title(title)

