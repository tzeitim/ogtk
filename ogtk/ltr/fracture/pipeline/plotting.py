
import os 
os.environ['OPENBLAS_NUM_THREADS'] = '4'

import polars as pl
from ogtk.utils.log import call

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
        for key in ['parsed_reads', 'parsed_reads_invalid']:
            ifn = results.results[key]

            out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_coverage.png'
            fig = sns.displot(data=
                pl.scan_parquet(ifn)
                    .group_by('umi')
                              .len().rename({'len':'reads'})
                    .collect(),
                    y='reads', 
                    log_scale=(10, 10), 
                    kind='ecdf', 
                    complementary=True, 
                    stat='count')

            plt.grid()
            plt.title(f"Reads per UMI\n{xp.target_sample}")
            plt.xlim((1,2e6))
            plt.ylim((1,1e5))

            #th_kmeans = qc.find_read_count_threshold(ifn, method='kmeans')
            #th_kneedle = qc.find_read_count_threshold(ifn, method='kneedle')

            #plt.axhline(y=th_kmeans, color='r', linestyle='--', label="kmeans")
            #plt.axhline(y=th_kneedle, color='g', linestyle='--', label="kneedle")
            plt.axhline(y=10, color='g', linestyle='--', label="hard coded")

            fig.savefig(out_path)
            xp.logger.info(f"saved {out_path}")
            
            # saturation coverage
            out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_sat-coverage.png'
            
            total_reads = pl.scan_parquet(ifn).collect().height # replace for direct extraction from the summary metrics

            _=[
                sns.ecdfplot(pl.scan_parquet(ifn)
                                .head(int(i))
                                .group_by('umi')
                                .len().rename({'len':'reads'})
                                .collect(), 
                             y='reads', 
                             complementary=True, 
                             stat='count',
                             linestyle='dotted',
                             )
                #for i in np.logspace(np.log10(1e3), np.log10(total_reads), num=5, base=10)
                for i in np.linspace(0, total_reads, num=5)
            ]

            #plt.title(f'{pl.scan_parquet(mfn_pcr).collect().height/1e6:.2f}')
            plt.title(f'{xp.target_sample} {total_reads/1e6:.2f}M reads')

            plt.xscale('log')
            plt.yscale('log')
            #plt.xlim((1,2e3))
            plt.xlim((1,2e6))
            plt.ylim((1,1e5))
            plt.grid(True, which='both', ls='--', alpha=0.3)

            fig.savefig(out_path)
            xp.logger.info(f"saved {out_path}")
            # anchor analysis

            out_path = f'{xp.sample_figs}/{xp.target_sample}_{key}_anchors.png'
            from . import qc
            anchor_stats_df = (
                    pl.concat(qc.compute_anchor_stats(ifn, sample_n=100, reps=16))
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
            fig = sns.catplot([
             qc.compute_double_anchor(ifn, 10, 30),
             qc.compute_double_anchor(ifn, 50,  30),
             qc.compute_double_anchor(ifn, 100,  30),
             qc.compute_double_anchor(ifn, 500,  30),
            ], kind='box')

            plt.ylim(0,1)
            fig.savefig(out_path)
            xp.logger.info(f"saved {out_path}")



