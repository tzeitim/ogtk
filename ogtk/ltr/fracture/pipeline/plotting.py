import polars as pl
from ogtk.utils.log import call

class PlotDB():
    @call
    # Pipeline instance = ppi
    def plot_preprocess(self, ppi, results):
        ''' ''' 
        sns = ppi.sns
        plt = ppi.plt
        xp = ppi.xp

        ifn = results.results['parsed_reads']
        out_path = f'{xp.sample_figs}/{xp.target_sample}_coverage.png'
        fig = sns.displot(data=
            pl.scan_parquet(ifn)
                .select('umi', 'reads').unique()
                .collect(),
                y='reads', 
                log_scale=(10, 10), 
                kind='ecdf', 
                complementary=True, 
                stat='count')

        plt.grid()
        plt.title("Reads per UMI")

        xp.logger.io(f"saved {out_path}")
        fig.savefig(out_path)

