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

def _generate_sample_sizes(max_sample_size, reps=1, log_scale=False, samples_rep=20):
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

    sizes = _generate_sample_sizes(max_sample_size, reps, **kw_sizes)  
    
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

def find_read_count_threshold(df, min_reads=10, method="kneedle"):
    """
    Find a threshold on the y-axis (read count) to separate real UMIs from noise
    using the knee/elbow detection method.
    
    Args:
        df: Polars DataFrame with 'umi' and 'reads' columns
        min_reads: Minimum read count to consider (default=1)
        plot: Whether to plot the result
        method: Detection method ('kneedle' or 'kmeans')
        
    Returns:
        read_count_threshold: The read count threshold to filter UMIs
    """
    # Get unique UMIs with their read counts and sort
    import os 
    os.environ['OPENBLAS_NUM_THREADS'] = '4'

    ranked_umis = (pl.scan_parquet(df)
                   .filter(pl.col('reads') >= min_reads)
                   .select('umi', 'reads')
                   .unique()
                   .sort('reads', descending=True)
                   .with_columns(
                      rank=pl.col('reads').rank('ordinal', descending=True)
                   )
                   .collect()
                  )
    
    if method == "kneedle":
        try:
            from kneed import KneeLocator
            
            # Extract data for knee detection and convert to float
            x = ranked_umis.get_column('rank').to_numpy().astype(float)
            y = ranked_umis.get_column('reads').to_numpy().astype(float)
            
            # Find the knee point (in this case, it's actually an elbow)
            kn = KneeLocator(
                x, y, 
                curve='convex', 
                direction='decreasing',
                interp_method='polynomial',
                online=True
            )
            
            # Get the threshold value
            if kn.knee is not None:
                knee_index = np.where(x == kn.knee)[0][0]
                read_count_threshold = y[knee_index]
            else:
                # Fallback if knee detection fails
                log_y = np.log10(y)
                log_diff = np.diff(log_y)
                
                # Find the point with the steepest drop
                steepest_point = np.argmin(log_diff) + 1  # +1 because diff reduces length by 1
                read_count_threshold = y[steepest_point]
        
        except ImportError:
            print("Kneed package not found, falling back to manual detection")
            # Manual detection of the elbow
            log_y = np.log10(ranked_umis.get_column('reads').to_numpy().astype(float))
            log_diff = np.diff(log_y)
            
            # Find the point with the steepest drop
            steepest_point = np.argmin(log_diff) + 1  # +1 because diff reduces length by 1
            read_count_threshold = ranked_umis.get_column('reads').to_numpy()[steepest_point]
            
    elif method == "kmeans":
        from sklearn.cluster import KMeans
        
        # Take log of both reads and rank for better clustering
        ranked_umis = ranked_umis.with_columns(
            log_reads=pl.col('reads').log10(),
            log_rank=pl.col('rank').log10()
        )
        
        # Extract data for clustering
        X = np.column_stack([
            ranked_umis.get_column('log_rank').to_numpy(),
            ranked_umis.get_column('log_reads').to_numpy()
        ])
        
        # Fit KMeans with 2 clusters
        kmeans = KMeans(n_clusters=2, random_state=42)
        clusters = kmeans.fit_predict(X)
        
        # Add cluster assignments back to the DataFrame
        ranked_umis = ranked_umis.with_columns(pl.lit(clusters).alias('cluster'))
        
        # Determine which cluster is signal (higher read counts)
        cluster_means = (ranked_umis
                        .group_by('cluster')
                        .agg(pl.col('log_reads').mean())
                        .sort('log_reads', descending=True))
        
        signal_cluster = cluster_means.get_column('cluster')[0]
        noise_cluster = 1 - signal_cluster  # The other cluster
        
        # Find boundary UMIs between clusters
        signal_umis = ranked_umis.filter(pl.col('cluster') == signal_cluster)
        noise_umis = ranked_umis.filter(pl.col('cluster') == noise_cluster)
        
        if len(signal_umis) > 0 and len(noise_umis) > 0:
            # Find the minimum read count in the signal cluster
            min_signal_reads = signal_umis.get_column('reads').min()
            # Find the maximum read count in the noise cluster
            max_noise_reads = noise_umis.get_column('reads').max()
            
            # Threshold is halfway between the two in log space
            log_min_signal = np.log10(min_signal_reads)
            log_max_noise = np.log10(max_noise_reads)
            log_threshold = (log_min_signal + log_max_noise) / 2
            read_count_threshold = 10 ** log_threshold
        else:
            # Fallback if clustering fails
            read_count_threshold = ranked_umis.get_column('reads').median()
    
    else:
        raise ValueError(f"Unknown method: {method}. Choose 'kneedle' or 'kmeans'.")
    
    return read_count_threshold
