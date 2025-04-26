import polars as pl
import numpy as np
import random
from typing import Dict, Any, Optional, Tuple, List, Union
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from scipy.optimize import curve_fit

from ogtk.utils.log import Rlogger, call

logger = Rlogger().get_logger()

@call

def compute_double_anchor(ifn, sample_n = 50, reps=10, start_anchor = "GAGACTGCATGG", end_anchor="TTTAGTGAGGGT"):
    """
    Returns ``reps`` groups of size ``sample_n`` to assess the number molecules that contain both anchor sequences.
    """
    return [
        pl.scan_parquet(ifn)
        .filter(pl.col('ont'))
        .filter(pl.col('umi').is_in(pl.col('umi').sample(sample_n)))
        .group_by('umi')
            .agg(start=pl.col('r2_seq').str.contains(start_anchor).sum(), 
                   end=pl.col('r2_seq').str.contains(end_anchor).sum())
        .with_columns(possible=(pl.col('start')>0) & (pl.col('end')>0))
        .collect()
        .get_column("possible").mean()
    for _ in range(reps)]

    
def compute_anchor_stats(ifn, sample_n = 50, reps=1, start_anchor = "GAGACTGCATGG", end_anchor="TTTAGTGAGGGT"):

    return [
        pl.scan_parquet(ifn)
        .filter(pl.col('ont'))
        .filter(pl.col('umi').is_in(pl.col('umi').sample(sample_n)))
        .group_by('umi')
            .agg(start=pl.col('r2_seq').str.contains(start_anchor).sum()/pl.len(), 
                   end=pl.col('r2_seq').str.contains(end_anchor).sum()/pl.len(),
                len=pl.len()
                )
        .with_columns(pl.lit(_).alias('_'))
        .collect()
    for _ in range(reps)]

def generate_sample_sizes(max_sample_size, reps=1, log_scale=False, samples_rep=20):

    if log_scale:
        early_points = np.array([10, 25, 50, 100, 250, 500])
        base_points = np.logspace(np.log2(1000), np.log2(max_sample_size), num=samples_rep - len(early_points), base=2)
        
        all_points = np.unique(np.concatenate([early_points, base_points]))
        all_points = all_points[all_points <= max_sample_size]
        
    else:
        all_points = np.linspace(1, max_sample_size, num=samples_rep) 
        all_points = all_points[all_points <= max_sample_size]
        
    return np.hstack([ np.array(all_points, dtype=int) for _ in range(reps)])


def compute_saturation_curve(ifn, name=None, max_sample_size=250_000, reps=3, kw_sizes={}, threshold=1, velocity=True):
    """
    ifn = "merged_reads.parquet"

    Output can be plotted using sns

    sns.lineplot(data=curves, x='reads', y='total_umis', hue='name')
    plt.xscale('log')
    plt.yscale('log')

    
    or 
    sns.relplot(data=curves, x='reads', y='total_umis', hue='name',kind='line')

    """
    if name is None:
        name = ifn.split('/')[-1]
        
    lazy_df = pl.scan_parquet(ifn)
    dataset_size = lazy_df.select(pl.len()).collect().item()

    sizes = generate_sample_sizes(max_sample_size, reps, **kw_sizes)  
    
    offsets = [random.randint(0, dataset_size - max(sizes)) for _ in range(len(sizes))]
    
    df_saturation = pl.concat(
        [
            lazy_df.slice(offset, int(size))
             .with_columns(
                 batch=pl.lit(offset), 
                 reads=pl.lit(size),
             )
            .with_columns(umi_cov = pl.len().over('umi'))
            .filter(pl.col('umi_cov')>=threshold)
            .collect()
            .group_by('reads', 'batch').agg(
                total_umis=pl.col('umi').n_unique(),
                mean_reads_umi =pl.col('umi_cov').mean(),
                median_reads_umi =pl.col('umi_cov').median(),
                max_reads_umi =pl.col('umi_cov').max(),
                min_reads_umi =pl.col('umi_cov').min(),
            )
         for _, (offset, size) in enumerate(zip(offsets, sizes))]
        )
        
    df_saturation = (
            df_saturation
            .with_columns(name=pl.lit(name))
            .sort('name', 'reads')
            .group_by('reads', 'name', maintain_order=True)
            .agg(pl.col('total_umis').mean())
            )

    if not velocity:
        return df_saturation

    else:
        return (
                df_saturation
                .with_columns(
                        (pl.col('total_umis').diff() / pl.col('reads').diff()).alias('velocity'))
                .sort('name', 'reads')
                .with_columns(pl.col('velocity').rolling_mean(min_samples=1, window_size=10).over('name'))
                .sort('name', 'reads')
        )

