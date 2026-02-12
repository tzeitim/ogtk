
import os
os.environ['OPENBLAS_NUM_THREADS'] = '4'

import polars as pl
from pathlib import Path
from ogtk.utils.log import call, Rlogger
from .formats import scan_file, read_file

logger = Rlogger().get_logger()


@call
def plot_segmentation_qc_standalone(
    segments_parquet: str,
    assembled_parquet: str,
    contigs_parquet: str,
    metas: pl.DataFrame,
    output_dir: str,
    sample_name: str,
    cassette_start_anchor: str = None,
    cassette_end_anchor: str = None,
    contig_col: str = 'contig',
    heterogeneity_threshold: float = 0.20,
):
    """
    Generate segmentation QC plots from existing parquet files.

    Can be called without running the full pipeline:

    >>> from ogtk.ltr.fracture.pipeline.plotting import plot_segmentation_qc_standalone
    >>> metas = pl.read_csv("PEtracer_metas.csv")
    >>> plot_segmentation_qc_standalone(
    ...     segments_parquet="workdir/sample/intermediate/segments.parquet",
    ...     assembled_parquet="workdir/sample/intermediate/assembled.parquet",
    ...     contigs_parquet="workdir/sample/contigs_segmented_valid.parquet",
    ...     metas=metas,
    ...     output_dir="workdir/sample/figures",
    ...     sample_name="group_001",
    ... )

    Args:
        segments_parquet: Path to segments parquet file
        assembled_parquet: Path to assembled segments parquet file
        contigs_parquet: Path to contigs parquet file
        metas: DataFrame with feature, seq, and optionally 'kind' columns
        output_dir: Output directory for plots
        sample_name: Sample name for file naming
        cassette_start_anchor: Optional cassette start anchor sequence
        cassette_end_anchor: Optional cassette end anchor sequence
        contig_col: Name of column containing contig sequences (default: 'contig')
        heterogeneity_threshold: For UMI heterogeneity detection, a transition type
            is considered "real" if its frequency is >= this fraction of the dominant
            type's frequency. Default 0.20 (20%). Set to 0 for sensitive mode only.

    Returns:
        dict: The generated segmentation report
    """
    from .segmentation import generate_segmentation_report, plot_segmentation_qc

    # Load data from parquet files
    segments_df = read_file(segments_parquet)
    assembled_df = read_file(assembled_parquet)
    contigs_df = read_file(contigs_parquet)

    # Rename contig column if needed (stitched_seq -> contig)
    if 'stitched_seq' in contigs_df.columns and contig_col not in contigs_df.columns:
        contigs_df = contigs_df.rename({'stitched_seq': contig_col})

    # Generate report from data
    report = generate_segmentation_report(
        segments_df=segments_df,
        assembled_df=assembled_df,
        contigs_df=contigs_df,
        metas=metas,
        cassette_start_anchor=cassette_start_anchor,
        cassette_end_anchor=cassette_end_anchor,
        heterogeneity_threshold=heterogeneity_threshold,
    )

    # Generate plots
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    plot_segmentation_qc(
        report=report,
        contigs_df=contigs_df,
        metas=metas,
        output_dir=output_path,
        sample_name=sample_name,
        contig_col=contig_col,
        logger=logger,
    )

    logger.info(f"Segmentation QC plots saved to {output_path}")

    return report


class PlotDB():
    @call
    # Pipeline instance = ppi
    def plot_preprocess(self, ppi, results):
        ''' ''' 
        import numpy as np

        sns = ppi.sns
        plt = ppi.plt
        xp = ppi.xp

        # TODO add a few sampling steps to include some info about saturation
        # TODO expose these do_ flagss
        do_invalid = False
        do_simple_cov  = False

        keys_to_process = ['parsed_reads']
        if do_invalid:
            keys_to_process.append('parsed_reads_invalid')

        for key in keys_to_process:
            ifn = results.results[key]

            # simple coverage
            # saturation coverage
            fig = sns.displot(data=
                scan_file(ifn)
                    .group_by('umi')
                              .len().rename({'len':'reads'})
                    .collect(engine='streaming'),
                    y='reads', 
                    log_scale=(10, 10), 
                    kind='ecdf', 
                    complementary=True, 
                    stat='count')

            plt.grid()
            plt.title(f"Reads per UMI\n{xp.target_sample}")
            plt.xlim((1,1e7))
            plt.ylim((1,1e4))

            #th_kmeans = qc.find_read_count_threshold(ifn, method='kmeans')
            #th_kneedle = qc.find_read_count_threshold(ifn, method='kneedle')

            #plt.axhline(y=th_kmeans, color='r', linestyle='--', label="kmeans")
            #plt.axhline(y=th_kneedle, color='g', linestyle='--', label="kneedle")
            if hasattr(xp, 'fracture'):
                plt.axhline(y=xp.fracture['min_reads'], color='g', linestyle='--', label="min reads")

            if do_simple_cov:
                out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_coverage.png'
                fig.savefig(out_path)
                xp.logger.info(f"saved {out_path}")
            
            out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_sat-coverage.png'
            
            total_reads = scan_file(ifn).select(pl.len()).collect().item()

            _=[
                sns.ecdfplot(scan_file(ifn)
                                .head(int(i))
                                .group_by('umi')
                                .len().rename({'len':'reads'})
                                .collect(engine="streaming"), 
                             y='reads', 
                             complementary=True, 
                             stat='count',
                             linestyle='dotted',
                             )
                #for i in np.logspace(np.log10(1e3), np.log10(total_reads), num=5, base=10)
                for i in np.linspace(0, total_reads, num=5)
            ]

            #plt.title(f'{scan_file(mfn_pcr).select(pl.len()).collect().item()/1e6:.2f}')
            plt.title(f'{xp.target_sample} {total_reads/1e6:.2f}M reads')

            plt.xscale('log')
            plt.yscale('log')
            #plt.xlim((1,2e3))
            plt.xlim((1,1e7))
            plt.ylim((1,1e4))
            plt.grid(True, which='both', ls='--', alpha=0.3)

            fig = plt.gcf()
            fig.savefig(out_path, bbox_inches='tight')
            xp.logger.info(f"saved {out_path}")
            # anchor analysis

            out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_anchors.png'
            from . import qc
            reps = 10
            sample_n=50
            anchor_stats_df = (
                    pl.concat(qc.compute_anchor_stats(ifn, sample_n=sample_n, reps=reps))
                        .with_columns(ifn=pl.lit(xp.target_sample))
                    )
            fig = sns.displot(data=anchor_stats_df, x='start', y='end',hue='ifn', col='ifn', bins=25)
            plt.ylim(0,1)
            plt.xlim(0,1)

            fig.savefig(out_path)
            xp.logger.info(f"saved {out_path}")
            
            # anchor analysis

            out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_feasible.png'

            from . import qc
            reps = 30
            fig = sns.catplot([
             qc.compute_double_anchor(ifn, 10, reps),
             qc.compute_double_anchor(ifn, 50, reps),
             qc.compute_double_anchor(ifn, 100, reps),
             qc.compute_double_anchor(ifn, 500, reps),
            ], kind='box')

            plt.ylim(0,1)
            fig.savefig(out_path)
            xp.logger.info(f"saved {out_path}")

        # Checkhealth QC - screen raw reads for expected features
        self._plot_checkhealth(ppi, results)

    def _plot_checkhealth(self, ppi, results):
        """Run checkhealth QC on raw reads to verify expected features are present."""
        from . import qc

        xp = ppi.xp
        plt = ppi.plt

        # Get checkhealth config early (needed for extras)
        checkhealth_cfg = getattr(xp, 'checkhealth', {})

        # Get features CSV path - same as metas
        features_csv = getattr(xp, 'features_csv', None)
        if features_csv is None and hasattr(xp, 'fracture'):
            features_csv = xp.fracture.get('metas_csv')

        # Build features from CSV (filtered by kind) and/or extras
        features_dfs = []
        feature_kinds = checkhealth_cfg.get('feature_kinds', ['ANCHOR'])

        if features_csv and Path(features_csv).exists():
            all_features = pl.read_csv(features_csv)
            if 'kind' in all_features.columns:
                filtered_features = all_features.filter(pl.col('kind').is_in(feature_kinds))
                if filtered_features.height > 0:
                    features_dfs.append(filtered_features.select([
                        pl.col('feature').alias('name'),
                        pl.col('seq').alias('sequence')
                    ]))

        # Add extras from config
        extras = checkhealth_cfg.get('extras', [])
        if extras:
            extras_df = pl.DataFrame({
                'name': [f'extra{i+1:02d}' for i in range(len(extras))],
                'sequence': extras
            })
            features_dfs.append(extras_df)

        if not features_dfs:
            xp.logger.info("No ANCHOR features or extras configured, skipping checkhealth")
            return

        features_df = pl.concat(features_dfs)

        # Find raw_bam file - use from results if provided, otherwise search
        raw_bam_path = None
        if hasattr(results, 'results') and results.results.get('raw_bam'):
            raw_bam_path = Path(results.results['raw_bam'])

        if not raw_bam_path or not raw_bam_path.exists():
            raw_bam_path = Path(xp.sample_wd) / 'intermediate' / 'raw_bam.parquet'
            if not raw_bam_path.exists():
                raw_bam_path = Path(xp.sample_wd) / 'intermediate' / 'raw_bam.arrow'

        if not raw_bam_path.exists():
            xp.logger.warning("Raw BAM file not found, skipping checkhealth")
            return

        xp.logger.info(f"Running checkhealth on {raw_bam_path}")

        # Get config options - reuse anchor_orient from main config as default
        sample_n = checkhealth_cfg.get('sample_n', 10000)
        orient_anchor = checkhealth_cfg.get('orient_anchor', getattr(xp, 'anchor_orient', None))
        orient_fuzzy = checkhealth_cfg.get('orient_fuzzy', 0)
        feature_fuzzy = checkhealth_cfg.get('feature_fuzzy', 0)

        # Sample reads (no grouping in auto-run)
        raw_ldf = scan_file(str(raw_bam_path)).head(sample_n)

        # Run checkhealth
        health_df = qc.checkhealth(
            raw_ldf,
            features_df,
            orient_anchor=orient_anchor,
            orient_fuzzy=orient_fuzzy,
            feature_fuzzy=feature_fuzzy,
        ).collect()

        xp.logger.info(f"Checkhealth processed {len(health_df):,} reads")

        # Save tabulated results
        out_dir = Path(xp.sample_figs)
        out_dir.mkdir(parents=True, exist_ok=True)

        health_table = qc.tabulate_health(health_df)
        health_table.write_csv(str(out_dir / f'{xp.target_sample}_checkhealth.csv'))
        xp.logger.info(f"Saved checkhealth table to {out_dir / f'{xp.target_sample}_checkhealth.csv'}")

        # Generate upset plot
        min_subset_size = max(1, int(len(health_df) * 0.01))
        qc.plot_upset(health_df, min_subset_size=min_subset_size, show_percentages=True)

        out_path = out_dir / f'{xp.target_sample}_checkhealth_upset.png'
        plt.savefig(str(out_path), bbox_inches='tight')
        xp.logger.info(f"Saved checkhealth upset plot to {out_path}")

    @call
    def plot_checkhealth(self, ppi, results):
        """
        Standalone checkhealth plot generation for regeneration.

        This is called by PlotRegenerator and uses results['raw_bam'] as input.
        """
        # Delegate to the internal implementation
        self._plot_checkhealth(ppi, results)

    @call
    def plot_segmentation(self, ppi, results):
        """
        Generate segmentation QC plots from pipeline results.

        Args:
            ppi: Pipeline instance with xp (experiment) config
            results: Results object with paths to intermediate files
        """
        xp = ppi.xp

        # Get file paths from results
        segments_parquet = results.results.get('segments')
        assembled_parquet = results.results.get('assembled')
        contigs_parquet = results.results.get('contigs_segmented_valid')

        if not all([segments_parquet, assembled_parquet, contigs_parquet]):
            logger.warning("Missing segmentation results, skipping QC plots")
            return

        # Load metas from config
        metas_csv = xp.fracture.get('metas_csv') if hasattr(xp, 'fracture') else None
        if not metas_csv:
            metas_csv = getattr(xp, 'features_csv', None)
        if not metas_csv:
            logger.warning("No metas_csv or features_csv in config, skipping segmentation plots")
            return

        metas = pl.read_csv(metas_csv)

        # Determine contig column name
        contig_col = 'contig'

        plot_segmentation_qc_standalone(
            segments_parquet=segments_parquet,
            assembled_parquet=assembled_parquet,
            contigs_parquet=contigs_parquet,
            metas=metas,
            output_dir=xp.sample_figs,
            sample_name=xp.target_sample,
            cassette_start_anchor=xp.fracture.get('cassette_start_anchor') if hasattr(xp, 'fracture') else None,
            cassette_end_anchor=xp.fracture.get('cassette_end_anchor') if hasattr(xp, 'fracture') else None,
            contig_col=contig_col,
        )
