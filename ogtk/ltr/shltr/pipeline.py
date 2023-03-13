from ogtk.utils import db
import ogtk.utils as ut
import ogtk.ltr.shltr as shltr
import scanpy as sc
from typing import Sequence,Optional
import polars as pl
import anndata as ad
import os
from functools import wraps


def _plot_cc(exp, 
    max_cell_size = 1.22e4,
    tmax_cell_size = 1000,
    lim= 2.5,
    cell_size_metric = "n_genes_by_counts",
    vmax=100,
    cluster_suffix='',
    pseudo_c=0,
    clone_numer='h1',
    clone_denom='g10',
    ):
    ''' cell size vs clone composition 
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    data = (exp.cc.join(ut.sc.adobs_pd_to_df(exp.ad_sc), 
                       left_on='cbc', 
                       right_on='cbc',
                       how='outer')
            .with_columns(pl.col('total_counts').fill_null(0))
            .with_columns(((pseudo_c+pl.col(clone_numer))/(pseudo_c+pl.col(clone_denom))).log10().alias('log_ratio')
                          ))

    data = data.filter(pl.col('total_counts')<max_cell_size)
    data = data.filter(pl.col('umis_cell')<tmax_cell_size)
    data = data.filter(pl.col('log_ratio')<=lim).filter(pl.col('log_ratio')>=-lim)


    #data = data.filter(pl.col('total_counts')>1e3)
    #data = data.filter(pl.col('umis_cell')<500)


    #ax.set_xlim((0, 5000))
    #ax.set_ylim((0, 500))

    fig, (ax, ax2, ax3) = plt.subplots(1,3, figsize=(3*15,15))
    sns.histplot(data=data.to_pandas(),
                 x=cell_size_metric,
                 y='log_ratio',
                 bins=100, ax= ax, cmap='RdYlBu_r', vmax=vmax)


    #ax.set_xlim((0, max_cell_size))
    ax.set_ylim((-lim, lim))
    ax.grid()
    ax.set_ylabel('H1/G10')

    sns.histplot(data=data.to_pandas(),
                 x='umis_cell',
                 y='log_ratio',
                 bins=100, ax= ax2, cmap='RdYlBu_r', vmax=vmax)

    ax2.set_xlim((0, tmax_cell_size))
    ax2.set_ylim((-lim, lim))
    ax2.grid()

    #ax.set_ylim((0, 500))
    ax2.set_ylabel('H1/G10')

    sns.histplot(data=data.to_pandas(),
                x=cell_size_metric,
                 y='umis_cell', 
                 bins=100,
                 log_scale=(False, False),
                 ax = ax3, cmap='RdYlBu_r', vmax=vmax)
    ax3.grid()
    plt.show()
    plt.figure()

    sns.displot(data=data.to_pandas(),
                x=cell_size_metric,
                y='log_ratio',
                bins=100,
                log_scale=(False, False),
                col=f'leiden{cluster_suffix}', 
                col_wrap=6, 
                cmap='RdYlBu_r', 
                vmax=50)

class Xp(db.Xp):
    def __init__(self, conf_fn):
        super().__init__(conf_fn)
        self.adata = None

    def load_10_mtx(self, cache = True, batch = None):
        ''' Loads the corresponding 10x matrix for a given experiment
        '''
        matrix_dir = f'{self.wd_scrna}/{self.sample_id}/outs/filtered_feature_bc_matrix/'
        self.print(f"importing matrix from:\n {matrix_dir}")

        adata = sc.read_10x_mtx(matrix_dir,
            var_names = 'gene_symbols',
            cache = cache
            )

        if batch is None:
            batch = self.sample_id

        adata.obs['batch'] = batch
        return(adata)

    def list_guide_tables(self,  suffix = 'shrna'):
        ''' Returns the file names of the available guide reads
        '''
        path = f'{self.wd_samplewd}/{suffix}'
        files = ut.sfind(path=path, pattern='*.parquet')
        return(files)

    def load_guide_molecules(
            self,
            clone,
            sample_id = None,
            suffix = 'shrna',
            index = 0,
            valid_ibars = None,
            down_sample = None,
            use_cache = True,
            corr_dir_fn=None):
        ''' Returns data frame at the molecule level
        '''
        if sample_id is None:
            sample_id = self.sample_id

        files = self.list_guide_tables(suffix)
        parquet_ifn = files[index]

        df = shltr.sc.reads_to_molecules(
            sample_id=sample_id,
            parquet_ifn=parquet_ifn,
            valid_ibars = valid_ibars, 
            use_cache=use_cache,
            corr_dict_fn=corr_dir_fn,
            down_sample=down_sample,
            #min_cells_per_ibar=int(len(obs27['cbc']) * 0.2),
            min_reads=2, 
            max_reads=1e6, 
            clone=clone)
        return(df)

    def init_mols(self,
                  valid_ibars: Sequence,
                  clone: str, 
                  ibar_ann: str = '/local/users/polivar/src/artnilet/workdir/scv2/all_parquetino_clusters',
                  min_cell_size: int = 0,
                  ):
        ''' Creates the .mols attribute
            Annotates the ibars using a 'cluster' table
            # TODO : improve the ibar annotation functionality
        '''

        self.mols = self.load_guide_molecules(valid_ibars = valid_ibars, clone = clone)
        
        dfc = (pl
           .scan_parquet(ibar_ann)
           .drop('sample_id').unique()
            #.filter(pl.col('sample_id')==sample_id)
            #.filter(pl.col('cluster')!=0)
            .collect()
            )

        self.mols = (self.mols.join(
                    dfc, 
                    left_on =['raw_ibar'], 
                    right_on=['raw_ibar'], 
                    how='left')
                    )

        self.score_clone(min_cell_size=min_cell_size)

    def init_adata(self):
        ''' Loads the raw cell ranger counts matrix
        '''
        self.raw_adata = self.load_10_mtx()
        #self.adata.obs['bath'] = self.sample_id

    @wraps(shltr.sc.compute_clonal_composition)
    def score_clone(self, min_cell_size=0, *args, **kwargs):
        fn = shltr.sc.compute_clonal_composition
        self.cc = fn(self.mols.filter(pl.col('umis_cell')>min_cell_size), *args, **kwargs)


    def init_object(self, sample_name: str| None= None):
        ''' Main function to load full experiment object
        '''
        if sample_name is None:
            sample_name = self.sample_id

        self.init_workdir()
        db.run_cranger(self, localcores=75, uiport=7777)

        db.tabulate_xp(self)

        self.load_final_ad(sample_name=sample_name)

        self.init_mols(valid_ibars = self.default_ibars(), clone=self.clone)
        #self.init_adata()

    def load_final_ad(self, sample_name):
        ''' Loads final version of mcs and scs adatas
        ''' 
        mcs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells_f.h5ad'
        scs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells_f.h5ad'
        
        if os.path.exists(mcs_ad_path) and os.path.exists(scs_ad_path):
            self.print(f'loading final adatas:\n{mcs_ad_path}\n{scs_ad_path}', 'bold white')

            self.ad_mc = ad.read_h5ad(mcs_ad_path)
            self.ad_sc = ad.read_h5ad(scs_ad_path)
        else:
            self.print(f'No final adatas found, compute them manually using .do_mc(explore=True), select best parameters and re-run .do_mc(explore=False)', 'bold red')

    def save_final_ad(self, sample_name):
        '''
        '''
        mcs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells_f.h5ad'
        scs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells_f.h5ad'

        self.print(f'saving final adatas:\n{mcs_ad_path}\n{scs_ad_path}', 'bold cyan')

        self.ad_mc.write_h5ad(mcs_ad_path)
        self.ad_sc.write_h5ad(scs_ad_path)

    def do_mc(self, sample_name: str | None= None, explore = True, forbidden_mods:Sequence | None=None, force = False):
        '''
        '''
        if sample_name is None:
            sample_name = self.sample_id

        mcs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells.h5ad'
        scs_ad_path = mcs_ad_path.replace('mcells', 'scells')
        outl_ad_path = mcs_ad_path.replace('mcells', 'outliers')

        if os.path.exists(scs_ad_path) and os.path.exists(mcs_ad_path) and not force:
            self.print('loading cached adatas')
            mcs = ad.read_h5ad(mcs_ad_path)
            scs = ad.read_h5ad(scs_ad_path)
            outliers = ad.read_h5ad(outl_ad_path)
            metadata = None
                                    
        else:

            if self.adata is None:
                self.init_adata()


            res = ut.sc.metacellize(set_name = sample_name,
                              adata=self.raw_adata, 
                              adata_workdir=self.wd_scrna, 
                              explore=explore, 
                              return_adatas=True,
                              forbidden_mods=forbidden_mods)
            if not explore:
                scs, mcs, outliers, metadata = res
        
        # cluster metacells
        mcs = ut.sc.scanpyfi(mcs.copy())
        # cluster single-cells
        scs = ut.sc.scanpyfi(scs.copy())

        # propagate metacell results to single-cells
        scs.obs = scs.obs.merge(
                right=mcs.obs.reset_index(drop=True), 
                left_on='metacell', 
                right_index=True, 
                how='outer', 
                suffixes=['', '_mc'])

        self.ad_sc = scs
        self.ad_mc = mcs
        self.outliers = outliers
        self.metadata = metadata
       
        self.save_final_ad(sample_name)

    def default_ibars(self):
            '''
                Loads pre-determined list of valid ibars
            '''
            self.print(':red_square: :red_square: Loading pre-computed valid ibars :red_square: :red_square:', 'bold red')
            valid_ibars = pl.read_csv('/local/users/polivar/src/artnilet/conf/valid_ibars.csv')['valid_ibars'].to_list()      
            return(valid_ibars)

    @wraps(_plot_cc)
    def plot_cc(self, *args, **kwargs):
        _plot_cc(self, *args, **kwargs)

