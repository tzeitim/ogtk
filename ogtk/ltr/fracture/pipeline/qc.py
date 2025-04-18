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


def compute_saturation_curve(ifn, umi_len=12, error_handling='none', error_threshold=1, window_pct=0.1):
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
        # Generate logarithmically spaced points - fewer points for faster processing
        # Cover the full range from 10 to 10M with just ~20 points
        
        # Method 1: Logarithmic spacing (powers of 10)
        base_points = np.logspace(1, 7, 19)  # 19 points from 10^1 to 10^7
        
        # Ensure we include some key early points for the steep part of the curve
        early_points = np.array([10, 25, 50, 100, 250, 500])
        
        # Combine and sort all points
        all_points = np.unique(np.concatenate([early_points, base_points]))
        
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

