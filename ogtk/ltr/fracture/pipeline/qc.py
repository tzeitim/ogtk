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

def plot_enhanced_saturation_curve(saturation_data, output_path=None, sample_name=None):
    """
    Visualize enhanced saturation curves with median coverage tracking.
    
    Parameters:
    -----------
    saturation_data : pl.DataFrame
        DataFrame with enhanced saturation curve data
    output_path : str, optional
        Path to save the plot image to. If None, plot is not saved.
    sample_name : str, optional
        Sample name to include in the plot title
        
    Returns:
    --------
    Tuple[plt.Figure, List[plt.Axes]]
        Figure and axes objects for further customization
    """
    # Convert to pandas if it's a polars DataFrame
    if isinstance(saturation_data, pl.DataFrame):
        df = saturation_data.to_pandas()
    else:
        df = saturation_data
    
    # Define saturation model (Michaelis-Menten kinetics)
    def saturation_model(x, vmax, km):
        return vmax * x / (km + x)
    
    # Create a 3-panel figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 18), 
                                        gridspec_kw={'height_ratios': [2, 1, 1]})
    
    # ---------- PANEL 1: UMI Saturation Curves ----------
    
    # Plot raw UMIs (all windows as semi-transparent points)
    sns.scatterplot(
        data=df, 
        x='size', 
        y='raw_umis',
        hue='window',  # Color by window
        alpha=0.5,
        s=30,
        palette='viridis',
        ax=ax1,
        legend=None  # Hide individual window legend
    )
    
    # Plot filtered UMIs
    if 'filtered_umis' in df.columns:
        sns.scatterplot(
            data=df, 
            x='size', 
            y='filtered_umis',
            color='red',
            alpha=0.7,
            s=30,
            ax=ax1,
            label='Filtered UMIs'
        )
    
    # Fit curve to raw UMI data
    x_data = df['size'].values
    y_data = df['raw_umis'].values
    raw_umis_est = None
    
    try:
        popt, _ = curve_fit(saturation_model, x_data, y_data)
        
        # Generate prediction points for smooth curve
        x_pred = np.logspace(1, np.log10(df['size'].max()), 100)
        y_pred = saturation_model(x_pred, *popt)
        
        # Plot the fitted curve
        ax1.plot(x_pred, y_pred, '--', color='blue', linewidth=2,
                label=f'Fitted curve (Raw UMIs, Vmax={int(popt[0])})')
        
        # Calculate and show estimated total UMIs
        raw_umis_est = int(popt[0])
        ax1.axhline(y=raw_umis_est, color='blue', linestyle=':', alpha=0.5)
        
        # Try to fit filtered UMI data if available
        if 'filtered_umis' in df.columns:
            try:
                # Filter out zeros or NaNs
                mask = (df['filtered_umis'] > 0) & (~df['filtered_umis'].isna())
                x_data_f = df.loc[mask, 'size'].values
                y_data_f = df.loc[mask, 'filtered_umis'].values
                
                if len(x_data_f) > 3:  # Need enough points for curve fitting
                    popt_f, _ = curve_fit(saturation_model, x_data_f, y_data_f)
                    y_pred_f = saturation_model(x_pred, *popt_f)
                    
                    ax1.plot(x_pred, y_pred_f, '--', color='red', linewidth=2,
                            label=f'Fitted curve (Filtered UMIs, Vmax={int(popt_f[0])})')
                    
                    filtered_umis_est = int(popt_f[0])
                    ax1.axhline(y=filtered_umis_est, color='red', linestyle=':', alpha=0.5)
                    
                    # Add saturation differential information
                    ratio = filtered_umis_est / raw_umis_est if raw_umis_est > 0 else 0
                    ax1.text(df['size'].min()*2, raw_umis_est*0.8, 
                            f'Estimated UMIs:\nRaw: {raw_umis_est:,}\nFiltered: {filtered_umis_est:,}\nRatio: {ratio:.2f}', 
                            color='black', fontweight='bold', bbox=dict(facecolor='white', alpha=0.7))
            except Exception as e:
                print(f"Could not fit filtered curve: {str(e)}")
            
    except Exception as e:
        print(f"Could not fit saturation curve: {str(e)}")
    
    # Calculate mean values per size across windows
    mean_df = df.groupby('size').agg({
        'raw_umis': 'mean',
        'filtered_umis': 'mean' if 'filtered_umis' in df.columns else None, 
        'raw_median_coverage': 'mean',
        'filtered_median_coverage': 'mean' if 'filtered_median_coverage' in df.columns else None,
    }).reset_index()
    
    # Remove any None columns that weren't available
    mean_df = mean_df.dropna(axis=1, how='all')
    
    # Plot the mean values
    ax1.plot(mean_df['size'], mean_df['raw_umis'], 'o-', color='black', 
            linewidth=2, markersize=6, label='Mean raw UMIs')
    
    if 'filtered_umis' in mean_df.columns:
        ax1.plot(mean_df['size'], mean_df['filtered_umis'], 'o-', color='darkred', 
                linewidth=2, markersize=6, label='Mean filtered UMIs')
    
    # Format main plot
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of reads', fontsize=14)
    ax1.set_ylabel('Number of UMIs', fontsize=14)
    title = 'UMI Saturation Curve'
    if sample_name:
        title += f' - {sample_name}'
    ax1.set_title(title, fontsize=16)
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=11)
    
    # ---------- PANEL 2: Efficiency curve ----------
    
    # Calculate efficiency (UMIs per read)
    raw_efficiency = mean_df['raw_umis'] / mean_df['size']
    filtered_efficiency = mean_df['filtered_umis'] / mean_df['size'] if 'filtered_umis' in mean_df.columns else None
    
    # Plot efficiency
    ax2.plot(mean_df['size'], raw_efficiency * 100, 'o-', color='blue', 
             linewidth=2, markersize=6, label='Raw UMIs per 100 reads')
    
    if filtered_efficiency is not None:
        ax2.plot(mean_df['size'], filtered_efficiency * 100, 'o-', color='red', 
                linewidth=2, markersize=6, label='Filtered UMIs per 100 reads')
    
    # Format efficiency plot
    ax2.set_xscale('log')
    ax2.set_xlabel('Number of reads', fontsize=14)
    ax2.set_ylabel('UMIs per 100 reads', fontsize=14)
    ax2.set_title('Sequencing Efficiency', fontsize=14)
    ax2.grid(True, which='both', linestyle='--', alpha=0.5)
    ax2.legend(loc='upper right')
    
    # ---------- PANEL 3: Median Coverage Scaling ----------
    
    # Plot median coverage vs sample size
    ax3.plot(mean_df['size'], mean_df['raw_median_coverage'], 'o-', color='blue', 
             linewidth=2, markersize=6, label='Raw median coverage')
    
    if 'filtered_median_coverage' in mean_df.columns:
        ax3.plot(mean_df['size'], mean_df['filtered_median_coverage'], 'o-', color='red', 
                linewidth=2, markersize=6, label='Filtered median coverage')
    
    # Try to fit power law to median coverage
    try:
        # Define power law model: y = a * x^b
        def power_law(x, a, b):
            return a * np.power(x, b)
        
        # Fit power law to median coverage data
        x_data = mean_df['size'].values
        y_data = mean_df['raw_median_coverage'].values
        
        # Log-transform to linearize for better fitting
        log_x = np.log(x_data)
        log_y = np.log(y_data)
        valid_mask = ~np.isnan(log_x) & ~np.isnan(log_y) & ~np.isinf(log_x) & ~np.isinf(log_y)
        
        if np.sum(valid_mask) > 2:  # Need at least 3 valid points
            params = np.polyfit(log_x[valid_mask], log_y[valid_mask], 1)
            a = np.exp(params[1])
            b = params[0]
            
            # Generate and plot fitted curve
            x_pred = np.logspace(np.log10(x_data.min()), np.log10(x_data.max()), 100)
            y_pred = power_law(x_pred, a, b)
            ax3.plot(x_pred, y_pred, '--', color='green', linewidth=2,
                    label=f'Power law fit: {a:.2f}*x^{b:.2f}')
            
            # Add annotation about scaling
            if b > 0.9:
                scaling_note = "Linear scaling (strong)"
            elif b > 0.7:
                scaling_note = "Near-linear scaling"
            elif b > 0.5:
                scaling_note = "Sublinear scaling (moderate)"
            elif b > 0.3:
                scaling_note = "Sublinear scaling (weak)"
            else:
                scaling_note = "Minimal scaling with depth"
                
            ax3.text(0.05, 0.95, f"Coverage scaling: {scaling_note}\nExponent: {b:.2f}", 
                    transform=ax3.transAxes, fontsize=12, fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.7), verticalalignment='top')
    except Exception as e:
        print(f"Could not fit power law to median coverage: {str(e)}")
    
    # If auto threshold was used, show threshold points
    if df['filtering_method'].iloc[0] == 'auto':
        ax3.scatter(df['size'], df['threshold_used'], color='black', marker='x', s=40, 
                   label='Auto threshold values')
    
    # Format median coverage plot
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('Number of reads', fontsize=14)
    ax3.set_ylabel('Median reads per UMI', fontsize=14)
    ax3.set_title('Median Coverage Scaling with Sequencing Depth', fontsize=14)
    ax3.grid(True, which='both', linestyle='--', alpha=0.5)
    ax3.legend(loc='upper left')
    
    # Finalize and save plot
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300)
        print(f"Saved enhanced saturation curve plot to {output_path}")
    
    return fig, [ax1, ax2, ax3]

def compute_enhanced_saturation_curve(ifn, umi_len=12, error_handling='none', error_threshold=0, window_pct=0.1, max_sample_size=500000):
    """
    Compute saturation curve with enhanced thresholding options and median tracking.
    
    Parameters:
    -----------
    ifn : str
        Input parquet file path
    umi_len : int, default=12
        Length of UMI sequence
    error_handling : str, default='none'
        Method to handle sequencing errors:
        - 'none': No error handling
        - 'threshold': Remove UMIs with read count below fixed threshold
        - 'auto': Remove UMIs with coverage below median for that sample size
    error_threshold : int, default=0
        Threshold for UMI filtering when error_handling='threshold'
        Set to 0 to use median-based filtering when error_handling='auto'
    window_pct : float, default=0.1
        Fraction of dataset to use for each sample window (0.1 = 10%)
    max_sample_size : int, default=500000
        Maximum sample size to use for saturation analysis
        
    Returns:
    --------
    pl.DataFrame
        DataFrame with saturation curve data including median coverage metrics
    """
    def generate_saturation_curve_sizes():
        # Generate logarithmically spaced points up to the specified maximum
        max_exponent = np.log10(max_sample_size)
        
        # Create log-spaced points from 10 to max_sample_size
        base_points = np.logspace(1, max_exponent, 15)
        
        # Ensure we include some key early points for the steep part of the curve
        early_points = np.array([10, 25, 50, 100, 250, 500])
        
        # Combine and sort all points
        all_points = np.unique(np.concatenate([early_points, base_points]))
        all_points = all_points[all_points <= max_sample_size]  # Limit to max size
        
        # Round to integers for readability
        return np.array(all_points, dtype=int)
    
    # Load the data and get dataset size
    lazy_df = pl.scan_parquet(ifn)
    dataset_size = lazy_df.select(pl.len()).collect().item()
    
    # Generate sampling sizes
    sizes = generate_saturation_curve_sizes()
    
    # Check if dataset is large enough and filter sizes if needed
    max_size = max(sizes)
    if dataset_size < max_size:
        logger.warning(f"Dataset size {dataset_size} is smaller than maximum sample size {max_size}")
        sizes = sizes[sizes <= dataset_size]
    
    # Calculate window size based on percentage of dataset
    window_size = int(dataset_size * window_pct)
    
    # Ensure window is large enough for largest sample size
    if window_size < max(sizes):
        window_size = max(sizes)
        logger.warning(f"Window size adjusted to {window_size} to accommodate largest sample size")
    
    # Calculate number of windows that can fit in the dataset
    num_windows = int(max(1, (dataset_size - window_size) // (window_size // 2)))
    
    # Use overlapping windows with 50% overlap
    window_starts = [i * (window_size // 2) for i in range(num_windows)]
    
    results = []
    
    # For each window
    for window_idx, window_start in enumerate(window_starts):
        # Load the window of data
        window_df = lazy_df.slice(int(window_start), int(window_size)).collect()
        
        # Skip empty windows
        if window_df.height == 0:
            continue
            
        # Add window identifier
        window_df = window_df.with_columns(
            window=pl.lit(window_idx),
            window_start=pl.lit(window_start),
            window_size=pl.lit(window_size)
        )
        
        # Extract UMI sequence
        window_df = window_df.with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
        
        # For each sample size within this window
        for size in sizes:
            # Skip if size is larger than window data
            if size > window_df.height:
                continue
                
            # Sample data progressively from the window
            sample_df = window_df.head(int(size))
            
            # Calculate raw UMI counts
            umi_counts = sample_df.group_by('umi').agg(count=pl.len())
            
            # Calculate overall median coverage before filtering
            raw_median_coverage = umi_counts.get_column('count').median()
            raw_mean_coverage = umi_counts.get_column('count').mean()
            raw_umis = umi_counts.height
            
            # Store original UMI counts for reference
            filtered_umi_counts = umi_counts.clone()
            
            # Apply error handling if specified
            filtered_umis = raw_umis  # Default to raw count
            actual_threshold = 0  # Track what threshold was actually used
            
            if error_handling == 'threshold' and error_threshold > 0:
                # Filter UMIs with counts below fixed threshold
                filtered_umi_counts = umi_counts.filter(pl.col('count') > error_threshold)
                filtered_umis = filtered_umi_counts.height
                actual_threshold = error_threshold
                
            elif error_handling == 'auto':
                # Use median as threshold
                if raw_median_coverage > 0:
                    filtered_umi_counts = umi_counts.filter(pl.col('count') >= raw_median_coverage)
                    filtered_umis = filtered_umi_counts.height
                    actual_threshold = raw_median_coverage
            
            # Calculate metrics after filtering
            filtered_median_coverage = filtered_umi_counts.get_column('count').median() if filtered_umis > 0 else 0
            filtered_mean_coverage = filtered_umi_counts.get_column('count').mean() if filtered_umis > 0 else 0
            
            # Calculate total reads before and after filtering
            total_raw_reads = umi_counts.get_column('count').sum()
            total_filtered_reads = filtered_umi_counts.get_column('count').sum() if filtered_umis > 0 else 0
            
            # Calculate coverage distribution metrics
            if raw_umis > 0:
                q1_coverage = umi_counts.get_column('count').quantile(0.25)
                q3_coverage = umi_counts.get_column('count').quantile(0.75)
                iqr_coverage = q3_coverage - q1_coverage
            else:
                q1_coverage = q3_coverage = iqr_coverage = 0
            
            # Create result with all metrics
            result = pl.DataFrame({
                'window': [window_idx],
                'window_start': [window_start],
                'size': [size],
                'raw_umis': [raw_umis],
                'filtered_umis': [filtered_umis],
                # Make sure these are explicitly converted to the right types
                'raw_median_coverage': [float(raw_median_coverage)],  # Enforce float type
                'raw_mean_coverage': [float(raw_mean_coverage)],  # Enforce float type
                'filtered_median_coverage': [float(filtered_median_coverage)],  # Enforce float type
                'filtered_mean_coverage': [float(filtered_mean_coverage)],  # Enforce float type
                'threshold_used': [float(actual_threshold)],  # Enforce float type
                'total_raw_reads': [int(total_raw_reads)],  # Enforce integer type
                'total_filtered_reads': [int(total_filtered_reads)],  # Enforce integer type
                'q1_coverage': [float(q1_coverage)],  # Enforce float type
                'q3_coverage': [float(q3_coverage)],  # Enforce float type
                'iqr_coverage': [float(iqr_coverage)],  # Enforce float type
                'filtering_method': [error_handling]
            })

            
            results.append(result)
    
    # Combine results from all windows
    if not results:
        logger.warning("No valid results generated for saturation curve")
        return pl.DataFrame({
            "size": [], 
            "raw_umis": [],
            "filtered_umis": [],
            "raw_median_coverage": [],
            "filtered_median_coverage": [],
            "threshold_used": [],
            "window": []
        })
        
    return pl.concat(results)

def plot_weighted_saturation_curve(saturation_data, output_path=None, sample_name=None):
    """
    Visualize the enhanced UMI saturation curve with weighted metrics.
    
    Parameters:
    -----------
    saturation_data : pl.DataFrame
        DataFrame with weighted saturation curve data
    output_path : str, optional
        Path to save the plot image to. If None, plot is not saved.
    sample_name : str, optional
        Sample name to include in the plot title
        
    Returns:
    --------
    Tuple[plt.Figure, List[plt.Axes]]
        Figure and axes objects for further customization
    """
    # Convert to pandas if it's a polars DataFrame
    if isinstance(saturation_data, pl.DataFrame):
        df = saturation_data.to_pandas()
    else:
        df = saturation_data
    
    # Define saturation model (Michaelis-Menten kinetics)
    def saturation_model(x, vmax, km):
        return vmax * x / (km + x)
    
    # Create a 3-panel figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 18), 
                                        gridspec_kw={'height_ratios': [2, 1, 1]})
    
    # ---------- PANEL 1: Main saturation curve ----------
    
    # Plot all windows as semi-transparent points
    sns.scatterplot(
        data=df, 
        x='size', 
        y='unique_umis',
        hue='window',  # Color by window
        alpha=0.5,
        s=30,
        palette='viridis',
        ax=ax1,
        legend=None  # Hide individual window legend to reduce clutter
    )
    
    # Plot weighted UMI count
    sns.scatterplot(
        data=df, 
        x='size', 
        y='weighted_umi_count',
        color='red',
        alpha=0.7,
        s=30,
        ax=ax1,
        label='Weighted UMI count'
    )
    
    # Plot effective UMI count
    sns.scatterplot(
        data=df, 
        x='size', 
        y='effective_umi_count',
        color='orange',
        alpha=0.7,
        s=30,
        ax=ax1,
        label='Effective UMI count'
    )
    
    # Fit curve to raw UMI data
    x_data = df['size'].values
    y_data = df['unique_umis'].values
    total_umis_est = None
    
    try:
        popt, _ = curve_fit(saturation_model, x_data, y_data)
        
        # Generate prediction points for smooth curve
        x_pred = np.logspace(1, np.log10(df['size'].max()), 100)
        y_pred = saturation_model(x_pred, *popt)
        
        # Plot the fitted curve
        ax1.plot(x_pred, y_pred, '--', color='blue', linewidth=2,
                label=f'Fitted curve (Raw UMIs, Vmax={int(popt[0])})')
        
        # Calculate and show estimated total UMIs
        total_umis_est = int(popt[0])
        ax1.axhline(y=total_umis_est, color='blue', linestyle=':', alpha=0.5)
        
        # Try to fit weighted UMI data
        try:
            # Filter out zeros or NaNs
            mask = (df['weighted_umi_count'] > 0) & (~df['weighted_umi_count'].isna())
            x_data_w = df.loc[mask, 'size'].values
            y_data_w = df.loc[mask, 'weighted_umi_count'].values
            
            if len(x_data_w) > 3:  # Need enough points for curve fitting
                popt_w, _ = curve_fit(saturation_model, x_data_w, y_data_w)
                y_pred_w = saturation_model(x_pred, *popt_w)
                
                ax1.plot(x_pred, y_pred_w, '--', color='red', linewidth=2,
                        label=f'Fitted curve (Weighted UMIs, Vmax={int(popt_w[0])})')
                
                weighted_umis_est = int(popt_w[0])
                ax1.axhline(y=weighted_umis_est, color='red', linestyle=':', alpha=0.5)
                
                # Add saturation differential information
                ratio = weighted_umis_est / total_umis_est if total_umis_est > 0 else 0
                ax1.text(df['size'].min()*2, total_umis_est*0.8, 
                        f'Estimated UMIs:\nRaw: {total_umis_est:,}\nWeighted: {weighted_umis_est:,}\nRatio: {ratio:.2f}', 
                        color='black', fontweight='bold', bbox=dict(facecolor='white', alpha=0.7))
        except Exception as e:
            print(f"Could not fit weighted curve: {str(e)}")
            
    except Exception as e:
        print(f"Could not fit saturation curve: {str(e)}")
    
    # Calculate mean values per size across windows
    mean_df = df.groupby('size').agg({
        'unique_umis': 'mean',
        'weighted_umi_count': 'mean', 
        'effective_umi_count': 'mean'
    }).reset_index()
    
    # Plot the mean values
    ax1.plot(mean_df['size'], mean_df['unique_umis'], 'o-', color='black', 
            linewidth=2, markersize=6, label='Mean raw UMIs')
    
    # Format main plot
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of reads', fontsize=14)
    ax1.set_ylabel('Number of UMIs', fontsize=14)
    title = 'UMI Saturation Curve'
    if sample_name:
        title += f' - {sample_name}'
    ax1.set_title(title, fontsize=16)
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=11)
    
    # ---------- PANEL 2: Efficiency curve ----------
    
    # Calculate efficiency (UMIs per read)
    efficiency = mean_df['unique_umis'] / mean_df['size']
    weighted_efficiency = mean_df['weighted_umi_count'] / mean_df['size']
    
    # Plot efficiency
    ax2.plot(mean_df['size'], efficiency * 100, 'o-', color='blue', 
             linewidth=2, markersize=6, label='Raw UMIs per 100 reads')
    ax2.plot(mean_df['size'], weighted_efficiency * 100, 'o-', color='red', 
             linewidth=2, markersize=6, label='Weighted UMIs per 100 reads')
    
    # Format efficiency plot
    ax2.set_xscale('log')
    ax2.set_xlabel('Number of reads', fontsize=14)
    ax2.set_ylabel('UMIs per 100 reads', fontsize=14)
    ax2.set_title('Sequencing Efficiency', fontsize=14)
    ax2.grid(True, which='both', linestyle='--', alpha=0.5)
    ax2.legend(loc='upper right')
    
    # ---------- PANEL 3: Weighting ratio plot ----------
    
    # Calculate weighted/raw UMI ratio
    df_clean = df.copy()
    df_clean['weighting_ratio'] = df_clean['weighted_umi_count'] / df_clean['unique_umis']
    
    # Plot weighting ratio
    sns.boxplot(
        data=df_clean, 
        x='size', 
        y='weighting_ratio',
        color='lightblue',
        ax=ax3
    )
    
    # Add best fit line for ratio trend
    mean_ratio = df_clean.groupby('size')['weighting_ratio'].mean().reset_index()
    ax3.plot(range(len(mean_ratio)), mean_ratio['weighting_ratio'], 'r-', linewidth=2)
    
    # Format ratio plot
    size_values = sorted(df_clean['size'].unique())
    ax3.set_xticks(range(len(size_values)))
    ax3.set_xticklabels([f"{x:,}" for x in size_values], rotation=45)
    ax3.set_xlabel('Number of reads', fontsize=14)
    ax3.set_ylabel('Weighted/Raw UMI ratio', fontsize=14)
    ax3.set_title('Error Assessment (Decreasing ratio indicates more errors)', fontsize=14)
    ax3.grid(True, which='both', linestyle='--', alpha=0.5)
    
    # Add ratio trend annotation
    first_ratio = mean_ratio['weighting_ratio'].iloc[0] if len(mean_ratio) > 0 else 0
    last_ratio = mean_ratio['weighting_ratio'].iloc[-1] if len(mean_ratio) > 0 else 0
    ratio_change = last_ratio - first_ratio
    
    if ratio_change < -0.05:
        trend_msg = "Significant error accumulation detected"
        color = "red"
    elif ratio_change < 0:
        trend_msg = "Mild error accumulation detected"
        color = "orange"
    else:
        trend_msg = "No significant error accumulation"
        color = "green"
        
    ax3.text(0.5, 0.05, trend_msg, transform=ax3.transAxes, 
            ha='center', fontsize=14, fontweight='bold', color=color,
            bbox=dict(facecolor='white', alpha=0.7))
    
    # Finalize and save plot
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300)
        print(f"Saved weighted saturation curve plot to {output_path}")
    
    return fig, [ax1, ax2, ax3]

def compute_weighted_saturation_curve(ifn, umi_len=12, error_handling='none', error_threshold=1, window_pct=0.1, max_sample_size=500000):
    """
    Compute enhanced saturation curve with weighted UMI metrics to account for sequencing errors.
    
    Parameters:
    -----------
    ifn : str
        Input parquet file path
    umi_len : int, default=12
        Length of UMI sequence
    error_handling : str, default='none'
        Method to handle sequencing errors:
        - 'none': No error handling
        - 'threshold': Remove UMIs with read count below threshold
        - 'adaptive': Use an adaptive threshold based on median coverage
    error_threshold : int, default=1
        Threshold for UMI filtering when error_handling='threshold'
    window_pct : float, default=0.1
        Fraction of dataset to use for each sample window (0.1 = 10%)
    max_sample_size : int, default=500000
        Maximum sample size to use for saturation analysis
        
    Returns:
    --------
    pl.DataFrame
        DataFrame with saturation curve data including weighted UMI metrics
    """
    def generate_saturation_curve_sizes():
        # Generate logarithmically spaced points up to the specified maximum
        max_exponent = np.log10(max_sample_size)
        
        # Create log-spaced points from 10 to max_sample_size
        base_points = np.logspace(1, max_exponent, 15)
        
        # Ensure we include some key early points for the steep part of the curve
        early_points = np.array([10, 25, 50, 100, 250, 500])
        
        # Combine and sort all points
        all_points = np.unique(np.concatenate([early_points, base_points]))
        all_points = all_points[all_points <= max_sample_size]  # Limit to max size
        
        # Round to integers for readability
        return np.array(all_points, dtype=int)
    
    def calculate_shannon_entropy(umi_counts, total_reads):
        """Calculate Shannon entropy of UMI distribution"""
        probabilities = umi_counts.get_column('count') / total_reads
        return -sum(p * np.log2(p) for p in probabilities if p > 0)

    def calculate_effective_umi_count(umi_counts):
        """Calculate effective number of UMIs (exponential of Shannon entropy)"""
        if umi_counts.height == 0:
            return 0
        total = umi_counts.get_column('count').sum()
        probabilities = umi_counts.get_column('count') / total
        shannon = -sum(p * np.log2(p) for p in probabilities if p > 0)
        return 2**shannon  # Effective number of UMIs

    def calculate_weighted_umi_count(umi_counts, median_coverage):
        """Calculate UMI count weighted by coverage relative to median"""
        if umi_counts.height == 0 or median_coverage == 0:
            return 0
        # Formula: sum(min(1, count/median)^alpha) where alpha controls weighting strength
        alpha = 2  # Adjust to control how quickly low-coverage UMIs are downweighted
        weights = np.minimum(1, (umi_counts.get_column('count') / median_coverage)**alpha)
        return weights.sum()
    
    # Load the data and get dataset size
    lazy_df = pl.scan_parquet(ifn)
    dataset_size = lazy_df.select(pl.len()).collect().item()
    
    # Generate sampling sizes
    sizes = generate_saturation_curve_sizes()
    
    # Check if dataset is large enough and filter sizes if needed
    max_size = max(sizes)
    if dataset_size < max_size:
        logger.warning(f"Dataset size {dataset_size} is smaller than maximum sample size {max_size}")
        sizes = sizes[sizes <= dataset_size]
    
    # Calculate window size based on percentage of dataset
    window_size = int(dataset_size * window_pct)
    
    # Ensure window is large enough for largest sample size
    if window_size < max(sizes):
        window_size = max(sizes)
        logger.warning(f"Window size adjusted to {window_size} to accommodate largest sample size")
    
    # Calculate number of windows that can fit in the dataset
    num_windows = int(max(1, (dataset_size - window_size) // (window_size // 2)))
    
    # Use overlapping windows with 50% overlap
    window_starts = [i * (window_size // 2) for i in range(num_windows)]
    
    results = []
    
    # For each window
    for window_idx, window_start in enumerate(window_starts):
        # Load the window of data
        window_df = lazy_df.slice(int(window_start), int(window_size)).collect()
        
        # Skip empty windows
        if window_df.height == 0:
            continue
            
        # Add window identifier
        window_df = window_df.with_columns(
            window=pl.lit(window_idx),
            window_start=pl.lit(window_start),
            window_size=pl.lit(window_size)
        )
        
        # Extract UMI sequence
        window_df = window_df.with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
        
        # For each sample size within this window
        for size in sizes:
            # Skip if size is larger than window data
            if size > window_df.height:
                continue
                
            # Sample data progressively from the window
            sample_df = window_df.head(int(size))
            
            # Calculate UMI counts
            umi_counts = sample_df.group_by('umi').agg(count=pl.len())
            
            # Apply error handling if specified
            if error_handling == 'threshold':
                # Filter UMIs with counts below fixed threshold
                umi_counts = umi_counts.filter(pl.col('count') > error_threshold)
            elif error_handling == 'adaptive':
                # Calculate median coverage
                median_coverage = umi_counts.get_column('count').median()
                # Set adaptive threshold as 10% of median coverage (minimum 1)
                adaptive_threshold = max(1, int(median_coverage * 0.1))
                # Filter UMIs with counts below adaptive threshold
                umi_counts = umi_counts.filter(pl.col('count') > adaptive_threshold)
            
            # Calculate standard metrics
            unique_umis = umi_counts.height
            total_reads = umi_counts.get_column('count').sum()
            median_coverage = umi_counts.get_column('count').median()
            
            # Skip if we have no valid UMIs
            if unique_umis == 0:
                continue
                
            # Calculate distribution metrics
            shannon_entropy = calculate_shannon_entropy(umi_counts, total_reads)
            effective_umi_count = calculate_effective_umi_count(umi_counts)
            weighted_umi_count = calculate_weighted_umi_count(umi_counts, median_coverage)
            
            # Calculate the ratio between weighted and raw UMI counts
            weighting_ratio = weighted_umi_count / unique_umis if unique_umis > 0 else 0
            
            # Create result with all metrics
            result = pl.DataFrame({
                'window': [window_idx],
                'window_start': [window_start],
                'size': [size],
                'unique_umis': [unique_umis],
                'total_reads': [total_reads],
                'median_coverage': [median_coverage],
                'shannon_entropy': [shannon_entropy],
                'effective_umi_count': [effective_umi_count],
                'weighted_umi_count': [weighted_umi_count],
                'weighting_ratio': [weighting_ratio],
                'error_handling': [error_handling]
            })
            
            results.append(result)
    
    # Combine results from all windows
    if not results:
        logger.warning("No valid results generated for saturation curve")
        return pl.DataFrame({
            "size": [], 
            "unique_umis": [], 
            "window": [],
            "shannon_entropy": [],
            "effective_umi_count": [],
            "weighted_umi_count": [],
            "median_coverage": [],
            "weighting_ratio": []
        })
        
    return pl.concat(results)

def compute_saturation_curve(ifn,
		umi_len=12,
		error_handling='none',
		error_threshold=1,
		window_pct=0.1,
		max_sample_size=1_000_000):
    """
    Compute saturation curve showing relationship between sampling depth and UMI diversity.
    Uses a sliding window approach to sample data at different depths.
     The function has been successfully updated with a sliding window approach. Here's a summary of the changes I made:

      1. Added the window_pct parameter to control the size of each sampling window (default 10% of dataset)
      2. Replaced random offset sampling with systematic sliding window approach:
        - Creates overlapping windows with 50% overlap
        - Each window represents a continuous chunk of reads from the dataset
        - Windows capture different regions of the sequencing run
      3. For each window:
        - Samples increasing subsets of reads (from smallest to largest)
        - Counts unique UMIs for each sample size
        - Applies error handling if specified
      4. Improved error handling:
        - Better dataset size checks
        - Window size adjustments based on maximum sample size
        - Empty window skipping

      This approach addresses the issue of UMI counts growing monotonically with read numbers by:

      1. Using progressive sampling within each window to more accurately model saturation
      2. Allowing error threshold filtering to remove low-count UMIs (likely sequencing errors)
      3. Providing multiple independent samples from different parts of the data, enabling statistical analysis

      The sliding window approach provides a more robust measurement of UMI saturation by reducing the impact of random sampling variation and capturing the true
      relationship between sequencing depth and molecular diversity.
    
    Parameters:
    -----------
    ifn : str
        Input parquet file path
    umi_len : int, default=12
        Length of UMI sequence
    error_handling : str, default='none'
        Method to handle sequencing errors:
        - 'none': No error handling
        - 'threshold': Remove UMIs with read count below threshold
        - 'cluster': Cluster similar UMIs together (not implemented yet)
    error_threshold : int, default=1
        Threshold for UMI filtering when error_handling='threshold'
    window_pct : float, default=0.1
        Fraction of dataset to use for each sample window (0.1 = 10%)
        
    Returns:
    --------
    pl.DataFrame
        DataFrame with saturation curve data
    """
    def generate_saturation_curve_sizes():
        max_exponent = np.log10(max_sample_size)
        base_points = np.logspace(1, max_exponent, 15)

        early_points = np.array([10, 25, 50, 100, 250, 500])
        
        all_points = np.unique(np.concatenate([early_points, base_points]))
        all_points = all_points[all_points <= max_sample_size]

        return np.array(all_points, dtype=int)
    
    lazy_df = pl.scan_parquet(ifn)
    dataset_size = lazy_df.select(pl.len()).collect().item()
    
    sizes = generate_saturation_curve_sizes()
    
    # Check if dataset is large enough and filter sizes if needed
    max_size = max(sizes)
    if dataset_size < max_size:
        logger.warning(f"Dataset size {dataset_size} is smaller than maximum sample size {max_size}")
        sizes = sizes[sizes <= dataset_size]
    
    # Calculate window size based on percentage of dataset
    window_size = int(dataset_size * window_pct)
    
    # Ensure window is large enough for largest sample size
    if window_size < max(sizes):
        window_size = max(sizes)
        logger.warning(f"Window size adjusted to {window_size} to accommodate largest sample size")
    
    # Calculate number of windows that can fit in the dataset
    num_windows =int(max(1, (dataset_size - window_size) // (window_size // 2)))
    
    # Use overlapping windows with 50% overlap
    window_starts = [i * (window_size // 2) for i in range(num_windows)]
    
    results = []
    
    # For each window
    for window_idx, window_start in enumerate(window_starts):
        # Load the window of data
        window_df = lazy_df.slice(int(window_start), int(window_size)).collect()
        
        # Skip empty windows
        if window_df.height == 0:
            continue
            
        # Add window identifier
        window_df = window_df.with_columns(
            window=pl.lit(window_idx),
            window_start=pl.lit(window_start),
            window_size=pl.lit(window_size)
        )
        
        # Extract UMI sequence
        window_df = window_df.with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
        
        # Apply error handling if specified
        if error_handling == 'threshold':
            # Count reads per UMI
            umi_counts = window_df.group_by('umi').agg(count=pl.count())
            
            # Filter UMIs with counts below threshold
            valid_umis = umi_counts.filter(pl.col('count') > error_threshold).get_column('umi')
            
            # Filter original dataframe to only include valid UMIs
            window_df = window_df.filter(pl.col('umi').is_in(valid_umis))
        
        # For each sample size within this window
        for size in sizes:
            # Skip if size is larger than window data
            if size > window_df.height:
                continue
                
            # Sample data progressively from the window
            # This ensures we're measuring saturation as we add more data
            sample_df = window_df.head(int(size))
            
            # Calculate unique UMIs at this sample size
            result = sample_df.group_by('window', 'window_start').agg(
                size=pl.lit(size),
                unique_umis=pl.col('umi').n_unique(),
                total_reads=pl.count()
            )
            
            results.append(result)
    
    # Combine results from all windows
    if not results:
        logger.warning("No valid results generated for saturation curve")
        return pl.DataFrame({"size": [], "unique_umis": [], "window": []})
        
    return pl.concat(results)


def plot_saturation_curve(saturation_data: pl.DataFrame, output_path: str = None, sample_name: str = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Visualize the UMI saturation curve data with a cleaner two-panel approach.
    
    Parameters:
    -----------
    saturation_data : pl.DataFrame
        DataFrame with saturation curve data from compute_saturation_curve()
    output_path : str, optional
        Path to save the plot image to. If None, plot is not saved.
    sample_name : str, optional
        Sample name to include in the plot title
        
    Returns:
    --------
    Tuple[plt.Figure, List[plt.Axes]]
        Figure and axes objects for further customization
    """
    # Convert to pandas if it's a polars DataFrame
    if isinstance(saturation_data, pl.DataFrame):
        df = saturation_data.to_pandas()
    else:
        df = saturation_data
    
    # Define saturation model (Michaelis-Menten kinetics)
    def saturation_model(x, vmax, km):
        return vmax * x / (km + x)
    
    # Create a 2-panel figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), gridspec_kw={'height_ratios': [2, 1]})
    
    # ---------- PANEL 1: Main saturation curve ----------
    
    # Plot all windows as semi-transparent points
    sns.scatterplot(
        data=df, 
        x='size', 
        y='unique_umis',
        hue='window',  # Color by window
        alpha=0.5,
        s=30,
        palette='viridis',
        ax=ax1,
        legend=None  # Hide individual window legend to reduce clutter
    )
    
    # Fit curve to all data
    x_data = df['size'].values
    y_data = df['unique_umis'].values
    total_umis_est = None
    
    try:
        popt, _ = curve_fit(saturation_model, x_data, y_data)
        
        # Generate prediction points for smooth curve
        x_pred = np.logspace(1, np.log10(df['size'].max()), 100)
        y_pred = saturation_model(x_pred, *popt)
        
        # Plot the fitted curve
        ax1.plot(x_pred, y_pred, '--', color='red', linewidth=2,
                label=f'Fitted curve (Vmax={int(popt[0])})')
        
        # Calculate and show estimated total UMIs
        total_umis_est = int(popt[0])
        ax1.axhline(y=total_umis_est, color='red', linestyle=':', alpha=0.5)
        ax1.text(df['size'].min()*2, total_umis_est*1.05, 
                f'Estimated total UMIs: {total_umis_est:,}', 
                color='red', fontweight='bold')
        
    except Exception as e:
        logger.warning(f"Could not fit saturation curve: {str(e)}")
    
    # Calculate mean values per size across windows
    mean_df = df.groupby('size').agg({'unique_umis': 'mean'}).reset_index()
    
    # Plot the mean values
    ax1.plot(mean_df['size'], mean_df['unique_umis'], 'o-', color='black', 
            linewidth=2, markersize=6, label='Mean across windows')
    
    # Format main plot
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of reads', fontsize=14)
    ax1.set_ylabel('Number of unique UMIs', fontsize=14)
    title = 'UMI Saturation Curve'
    if sample_name:
        title += f' - {sample_name}'
    ax1.set_title(title, fontsize=16)
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=11)
    
    # Add sequencing depth efficiency markers at 50%, 75%, 90% of estimated total UMIs
    # Position text labels to prevent overlap
    y_positions = []
    
    try:
        if total_umis_est:
            total_est = total_umis_est
            depths = [0.5, 0.75, 0.9]
            colors = ['green', 'orange', 'red']
            
            # Calculate y-position spacing for labels
            y_min = df['unique_umis'].min()
            
            for i, (depth, color) in enumerate(zip(depths, colors)):
                # Calculate reads needed to reach this UMI depth
                reads_needed = (depth * total_est * popt[1]) / (total_est - depth * total_est)
                if reads_needed <= df['size'].max():
                    # Draw marker
                    ax1.axvline(x=reads_needed, color=color, linestyle='--', alpha=0.5)
                    
                    # Calculate non-overlapping y position (staggered vertically)
                    y_pos = y_min * (2 + i*1.5)  # Stagger vertically
                    
                    ax1.text(reads_needed*1.1, y_pos, 
                            f'{int(depth*100)}% UMIs\n({int(reads_needed):,} reads)', 
                            color=color, fontweight='bold', rotation=0)
    except Exception as e:
        pass  # Skip efficiency markers if there's an error
    
    # ---------- PANEL 2: Efficiency curve ----------
    
    # Calculate efficiency (UMIs per read)
    efficiency = mean_df['unique_umis'] / mean_df['size']
    
    # Plot efficiency
    ax2.plot(mean_df['size'], efficiency * 100, 'o-', color='darkred', 
             linewidth=2, markersize=6, label='UMIs per 100 reads')
    
    # Format efficiency plot
    ax2.set_xscale('log')
    ax2.set_xlabel('Number of reads', fontsize=14)
    ax2.set_ylabel('UMIs per 100 reads', fontsize=14)
    ax2.set_title('Sequencing Efficiency', fontsize=14)
    ax2.grid(True, which='both', linestyle='--', alpha=0.5)
    ax2.legend(loc='upper right')
    
    # Add efficiency threshold markers if possible
    if total_umis_est:
        try:
            # Add annotated reference lines showing efficiency at key depths
            for depth, color in zip(depths, colors):
                reads_needed = (depth * total_est * popt[1]) / (total_est - depth * total_est)
                if reads_needed <= mean_df['size'].max():
                    eff_at_depth = saturation_model(reads_needed, *popt) / reads_needed * 100
                    ax2.axvline(x=reads_needed, color=color, linestyle='--', alpha=0.5)
                    ax2.text(reads_needed*1.1, eff_at_depth*0.9, 
                            f'{int(depth*100)}% saturation\n({eff_at_depth:.2f} UMIs/100 reads)', 
                            color=color, fontsize=9)
        except Exception:
            pass
    
    # Finalize and save plot
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300)
        logger.info(f"Saved saturation curve plot to {output_path}")
    
    return fig, [ax1, ax2]

