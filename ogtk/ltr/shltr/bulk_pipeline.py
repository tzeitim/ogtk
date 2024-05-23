import ogtk.utils as ut
from ogtk.utils import db
from ogtk.utils import log
import ogtk.ltr.shltr.ibars as ibars
import ogtk.ltr.shltr.bulk as bulk
import polars as pl
from pyseq_align import NeedlemanWunsch, SmithWaterman
import matplotlib.pyplot as plt
from functools import wraps
import os
import seaborn as sns

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()

def _list_files(self, pattern):
    return ut.sfind(self.wd_xp, pattern = pattern)

def _list_raw_fastqs(self):
    self.input_files = [i 
                        for i in 
                        ut.sfind(self.wd_fastq, pattern = "*_R1_*fastq.gz") 
                        if not i.split('/')[-1].startswith("Undetermined") 
                        ]


@log.call
def _filter_ibars(df: pl.DataFrame, filters:dict)->pl.DataFrame:
    '''
    Filters according to the keys provided in the filters dictionary.
    
    Supported filters:
     - global_min_mols_ibar
     - min_mols_ibar
     - whitelist

    '''
    assert isinstance(filters, dict), \
            "Please specify a dictionary with: 'min_mols_ibar' or 'whitelist' "

    if 'global_min_mols_ibar' in filters:
        assert "ibar_umisg" in df.columns, \
                "``ibar_umisg`` field missing compute using `.qc_ibars_per_sample`"
        df = df.filter(pl.col('ibar_umisg')>=filters['global_min_mols_ibar'])

    if 'min_mols_ibar' in filters:
        assert "ibar_umis" in df.columns, \
                "``ibar_umis`` field missing compute using `.qc_ibars_per_sample`"
        df = df.filter(pl.col('ibar_umis')>=filters['min_mols_ibar'])

    if 'whitelist' in filters:
        df = df.filter(pl.col('raw_ibar').is_in(filters['whitelist']))

    return df

@log.call
def _generate_molecule_dfs(
        xps,
        cache_dir,
        min_reads=1,
        max_reads=1e4,
        do_sampling=None,
        suffix:str = 'shrna',
        *args,
        **kwargs
    )-> pl.DataFrame:
    ''' xps maps a sample_id to a corresponding parquet path
    '''
    assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"

    df = []
    for path, clone, sample_id  in zip(xps['path'], xps['clone'], xps['sample_id']):
        df.append(
                ibars.reads_to_molecules(
                    sample_id=sample_id,
                    modality='single-molecule',
                    zombie=suffix == 'zshrna',
                    parquet_ifn=path,
                    min_reads=int(min_reads),
                    max_reads=int(max_reads),
                    clone=clone,
                    cache_dir=cache_dir,
                    do_sampling=do_sampling,
                    *args,
                    **kwargs,
                )
        )

    df = pl.concat(df)
    return df

def _qc_ibars_per_sample(df: pl.DataFrame)-> pl.DataFrame:
    '''
        Computes molecule stats at the level of samples and ibar
            - total reads per ibar (`ibar_reads`)
            - total umis per ibar (`ibar_umis`)
            - log10 total reads per ibar (`ibar_reads_log10`)
            - log10 total umis per ibar (`ibar_umis_log10`)
            - fraction of reads in sample
            - fraction of umis in sample

        Computes at the whole df level the quantiled counts
            - total reads per ibar (`ibar_readsg`)
            - total umis per ibar (`ibar_umisg`)
            - -log10 quantile conversion total reads per ibar (`ibar_reads_q`)
            - -log10 quantile conversion total umis per ibar (`ibar_umis_q`)
    '''
    groups = ['sample_id', 'raw_ibar']
    df = (
        df
         .with_columns(
                ibar_reads=pl.col('umi_dom_reads').sum().over(groups),
                ibar_umis=pl.col('umi').n_unique().over(groups),
                ibar_readsg=pl.col('umi_dom_reads').sum().over('raw_ibar'), 
                #ibar_umisg=pl.col('umi').n_unique().over('raw_ibar'),# this is a very good approximation but duplicated UMIs are not resolved 
             )

         .with_columns(
                ibar_reads_log10=pl.col('ibar_reads').log10(), 
                ibar_umis_log10=pl.col('ibar_umis').log10(),

                ibar_readsf=pl.col('ibar_reads')/pl.col('umi_dom_reads').sum().over('sample_id'),
                ibar_umisf=pl.col('ibar_umis')/pl.col('umi').n_unique().over('sample_id'),

                ibar_umisg=pl.col('ibar_umis').unique().sum().over('raw_ibar'),
            )

        .with_columns(
            ibar_reads_q=(-1 * (pl.col('ibar_reads')/(pl.col('umi_dom_reads').sum())).log10()),

            ibar_umisg_log10=pl.col('ibar_umisg').log10(),
            ibar_readsg_log10=pl.col('ibar_readsg').log10(), 
            )

        .with_columns(
            ibar_umis_q=(-1 * (pl.col('ibar_umis')/(pl.col('umi').n_unique())).log10())
            )
         )
    return df

##plotting
def _plot_ibar_qc(
        self,
        df:None|pl.DataFrame=None,
        groups=['sample_id','cluster'],
        top_n = 100,
    ):
    ''' Plots related to ibar coverage 
    '''
    if df is None:
        df = self.mols

    assert all([i in df.columns for i in groups]),\
            f"{','.join([i for i in groups if i not in self.mols.columns])} not found in .mols"
    
    # visualize allocation of reads and molecules across samples as a function of groups
    fig, (preads, pmols) = plt.subplots(1,2, figsize=(5*3.5, 5))

    sns.heatmap(
        df.sort(groups)
        .pivot(index=[groups[0]], 
               columns=groups[1],
               values='umi_reads', 
               aggregate_function='sum')
        .to_pandas().set_index('sample_id'),
        ax=preads
    ).set_title("Reads")

    sns.heatmap(
        df.sort(groups)
        .pivot(index=[groups[0]], 
               columns=groups[1],
               values='umi', 
               aggregate_function='len')
        .to_pandas().set_index('sample_id'),
        ax=pmols,
    ).set_title("Molecules")
    plt.tight_layout()
    plt.show()


    # visualization of top_n integrations 
    if 'cluster' in df.columns:
        #In this example we extract the top n ibars by cluster 
        data = (
            df
            .select(['raw_ibar', 'ibar_umis', 'ibar_umisg', 'ibar_umisg_log10', 'clone', 'cluster'])
            .unique()
            .filter(
               pl.col('ibar_umisg').rank(descending=True, method='dense').over('cluster') <= top_n
            )
        )


        g = sns.displot(
            data = data.to_pandas(),
            x='ibar_umis', hue='clone', 
            bins = 50, 
            col='cluster',
        )
        g.figure.suptitle(f"Top {top_n}")
        plt.show()
    
        g = sns.displot(
                data=data
                    .filter(pl.col('ibar_umisg')>=100) # visualize only integrations with more than 100 molecules
                    .select("^raw_ibar|cluster|ibar_.*sg.*$").unique()    
                    .sort('ibar_umisg', descending=True) 
                .to_pandas(),
                complementary=True,
                y='ibar_umisg_log10',
                stat='count',
                kind='ecdf',
                col='cluster',
        )
        plt.show()

        print(data.select('cluster','raw_ibar').unique().group_by('cluster').len())


def _plot_qc_seq(
        self,
        meta_df,
        xp_labels='sample_id',
    )->None:
    ''' Plots all fields with a "qc" prefix to visualize the different sequencing metrics.
        It additonally plots a 2D histogram of how molecules distributed across ``groups``
    '''
    paths = self.list_files("*qc_stats.parquet")
    df_qc = pl.concat([pl.scan_parquet(i) for i in paths]).collect()

    qc_stats = [i for i in df_qc.columns if i.startswith('qc')]

    for i in qc_stats:
        fg=sns.catplot(
           data=df_qc.join(meta_df, left_on='sample_id', right_on='sample_id', how='left')
                     .sort([xp_labels, 'clone'])
                     .select([xp_labels, i])
                     .unique(maintain_order=True)
               .to_pandas(), 
           y=i, 
           x=xp_labels, 
           kind='bar', 
           aspect=1.5) 

        if 'pc' in i:
            fg.ax.set_ylim((0,1))

        fg.ax.grid()
        fg.ax.set_title(i)
        fg.ax.tick_params(axis='x', rotation=90)

        plt.show()


def _qc_ibar_expression_levels() -> None:
    ''' characterize levels of expression of individual ibars across groups '''
    return None

def _load_ibar_cluster_annotation(path):
    ibar_clusters = (pl
       .scan_parquet(path)
       .drop(['sample_id']).unique()
       .rename({'clone':'clone_ann'})
       .collect()
      )
    return ibar_clusters

# Steps:
# - Load molecules
# - annotate ibars (if meta data provided)
# - filter integrations
# - analysis

class Xp(db.Xp):
    sample_id = ''
    wd_shrna = ''
    mols = None|pl.DataFrame
    emp_ann = None|pl.DataFrame
    path_ibar_clusters = None
    ibar_clusters = None
    input_files = ''
    wd_raw_input = ''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.conf_keys = [i for i in vars(self)]
        # paths
        self.path_mols =  f'{self.wd_xp}/{self.sample_id}.mols.parquet' 

    def init_experiment(self, xps, meta_df, filters, force=False):
        self.set_metadata(meta_df)
        self.init_mols(xps=xps, filters=filters, force=force)
        self.plot_qc_seq()

    def init_mols(self, xps:None|pl.DataFrame=None, force=False, do_annotation=True, filters=None, *args, **kwargs):
        ''' 
        Populates the .mols attribute by:
            - collapsing data from reads to molecules
            - performs QC at the level of samples
            - filters according to the ``filters`` argument (dictionary)
            - annotates ibars (if ``do_annotation`` is passed)

        Arguments:
            - do_annotation: dictionary with any key in ['min_mols_ibar', 'whitelist']
        '''

        for suffix in self.valid_tab_suffixes():
            path_mols = self.path_mols.replace('mols', f'mols.{suffix}')
            if os.path.exists(path_mols) and not force:
                self.logger.io(f"loading {path_mols}") #pyright:ignore
                self.mols = pl.read_parquet(path_mols)
            else:
                assert xps is not None, "Please provide a xps data frame"
                self.mols = _generate_molecule_dfs(
                     xps=xps,
                     cache_dir=getattr(self, f'wd_{suffix}'),
                     force=force,
                     suffix=suffix,
                     *args,
                     **kwargs) 

                # compute sample-specific qc metrics
                self.mols = self.qc_ibars_per_sample()

                if filters is not None:
                    self.mols = self.filter_ibars(filters=filters)

                if do_annotation and self.ibar_clusters is not None:
                    self.metadata_annotate_ibars()

                self.logger.io(f"saving mols to {path_mols}")
                self.mols.write_parquet(path_mols)
        return None
    
    def metadata_annotate_ibars(self, df:None|pl.DataFrame=None, on=['raw_ibar'], how='left'):
        ''' 
        Annotates ibars according to ``df``. If none is passed and
        .ibar_clusters is populated, then .ibar_clusters is used for the annotation 

        If applied recursively, redundant fields are eliminated.
        '''
        if not isinstance(on, list):
            on = [on]

        if df is None and self.ibar_clusters is not None:
            df = self.ibar_clusters

        i = set(self.mols.columns).intersection(df.columns).difference(on) #pyright: ignore
        assert 'raw_ibar' in df.columns, "Missing 'raw_ibar' field from annotation data frame" #pyright:ignore

        self.logger.step('annotating ibar clusters with .ibar_clusters') #pyright: ignore

        df_ann = (                                  #pyright: ignore
                self.mols
                     .join(                            #pyright: ignore
                       other=df.drop(i).unique(),      #pyright: ignore
                       left_on=on,
                       right_on=on,
                       how=how,
                        )
        )

        if all([ i in self.mols.columns for i in ['clone', 'clone_ann']]):
            mismatch_clone = self.mols.filter(pl.col('clone')!=pl.col('clone_ann')).shape[0]
            self.logger.critical(f'a total of {mismatch_clone} didnt match the expected clone ')

        if any(["_right" in i for i in self.mols.columns]):
            self.logger.info("removing redundant fields")
            df_ann = self.mols.select(pl.all().exclude('^.+_right$'))

        return df_ann

    def quick_screen(self, min_umis_ibar:int=100, top_n:int=100)->None:
        '''
        Performs a quick screen of integration identity based on the provided exp.ibar_clusters table
        It assists on determining the expected number of integrations according the provided annotations. 
        In addition it helps identifying meaningful integrations that are not included in the exp.ibar_clusters data frame
        '''

        assert self.ibar_clusters is not None, "A .ibar_clusters attribute must be defined in the experiment"
        
        quick_pass = (
            self.mols
            .filter(pl.col('ibar_umis')>=min_umis_ibar)
            .join(
                other=self.ibar_clusters.select('raw_ibar', 'clone_ann', 'cluster').unique(), 
                left_on=['raw_ibar', 'clone'], 
                right_on=['raw_ibar', 'clone_ann'], how='left'
                ) # polars eliminates the second field F
        )
        
        quick_pass = (
            quick_pass.with_columns(
                cluster=pl
                    .when(pl.col('cluster').is_null())
                    .then(pl.lit('unknown'))
                    .otherwise(pl.col('cluster'))
            )
        )
        
        # let's visualize how many integration of each type there are by cluster
        quick_pass['cluster'].value_counts().sort('count').to_pandas().plot.bar(x='cluster', y='count')
        quick_pass.group_by('cluster', 'raw_ibar').len().group_by('cluster').len().to_pandas().plot.bar(x='cluster', y='len')

        self.plot_ibar_qc(df=quick_pass, top_n=top_n) #pyright:ignore

    def empirical_annotation(self, emp_ann_input_df:pl.DataFrame|None=None, pattern='0p', n_exp_ibars=90):
        '''
        Generates a data frame that determines the spacer-ibar pairs and stores it in the .em_ann slot.

        If several samples are provided, pattern allows to filter for the ones where uncut states are expected.
        
        n_exp_ibars: number of expected integrations at a sample level

        Note: Collision may arise when inputing more data from more than one clone. Special care needs to be 
        taken in order to handle this properly.
        '''
        
        if emp_ann_input_df is None: 
            emp_ann_input_df = (
                self.mols
                .filter(pl.col('sample_id').str.contains(pattern))# This filter selects for the un-induced samples
                .filter(pl.col('ibar_umisg').rank(descending=True, method='dense').over('sample_id') <= n_exp_ibars)
        )

        assert emp_ann_input_df.shape[0]>0, "Empty data frame. Adjust filters or patterns for un-induced samples"
        
        n_samples = emp_ann_input_df['sample_id'].n_unique()

        emp_ann_df = ibars.empirical_metadata_annotation(
            df= emp_ann_input_df, 
            drop_counts=False, 
            top_n=1,
        ).with_columns(pl.col('sample_id').str.split('_').list.get(0).alias('clone'))
                          
        sns.displot(data=emp_ann_df.to_pandas(),
                 x='nlen',
                 hue='sample_id',
                 log_scale=10,
                 bins=60,
                 multiple='stack',
                 aspect=1.5)

        plt.show()

        # There can be collisions across clones. 
        # one option is to filter based on the nlen e.g. >=1e1
        # second option assign the sample with the higher number of molecules 

        # let's use the second option
        emp_ann_df = emp_ann_df.with_columns(pl.len().over('raw_ibar').alias('coll'))
        self.logger.info(f"Total collisions across samples: {emp_ann_df.filter(pl.col('coll')>1).shape[0]}/{n_samples}")

        emp_ann_df = emp_ann_df.filter(pl.col('len').rank(method='dense', descending=True).over('raw_ibar')==1)
        self.logger.info(f"Total integrations kept after filtering: {emp_ann_df.shape[0]}")

        self.emp_ann = emp_ann_df
        self.logger.info("Populated .emp_ann")

    def save_mols(self, suffix):
        ''' ''' 
        assert suffix in self.valid_tab_suffixes(), \
                f"Please provide a valid suffix from {self.valid_tab_suffixes()}"

        path_mols = self.path_mols.replace('mols', f'mols.{suffix}')
        self.mols.write_parquet(path_mols)

        return path_mols 

    @wraps(_qc_ibars_per_sample)
    def qc_ibars_per_sample(self,  *args, **kwargs):
        return _qc_ibars_per_sample(self.mols, *args, **kwargs)

    @wraps(_filter_ibars)
    def filter_ibars(self, *args, **kwargs):
        return _filter_ibars(self.mols, *args, **kwargs)

    def fastq_to_parquet(self, cbc_len=0, umi_len=25, force=False):
        # TODO improve the caching scheme. Currently force only has an effect
        # on the tabulation itself but not on the feature calling
        self.parquets = db.tabulate_xp(self, modality='single-molecule',
                                       cbc_len=cbc_len, umi_len=umi_len,
                                       force=force)

    @wraps(ibars.align_alleles)
    def align_alleles(self, *args, **kwargs):
        ''' 
        '''
        self.mols = ibars.align_alleles(self.mols, *args, **kwargs)

    @wraps(_load_ibar_cluster_annotation)
    def load_ibar_cluster_annotation(self, path: str | None =None):
        ''' '''
        if self.path_ibar_clusters is None and path is None:
            raise ValueError("please provide a path")
       
        path = self.path_ibar_clusters if path is None else path
        self.ibar_clusters = _load_ibar_cluster_annotation(path=path)

    @wraps(_plot_qc_seq)
    def plot_qc_seq(self, df:pl.DataFrame|None=None):
        if df is None:
            return _plot_qc_seq(self, meta_df=self.meta_df)
        return _plot_qc_seq(self, meta_df=self.meta_df)

    @wraps(_plot_ibar_qc) 
    def plot_ibar_qc(self, df:pl.DataFrame|None=None, *args, **kwargs):
        '''
        '''
        _plot_ibar_qc(self, df=df, *args, **kwargs)

    @wraps(_list_files)
    def list_files(self, pattern:str):
        return _list_files(self, pattern)

    @wraps(_list_raw_fastqs)
    def list_raw_fastqs(self):
        _list_raw_fastqs(self)


    def set_metadata(self, meta_df):
        ''' '''
        self.meta_df = meta_df
