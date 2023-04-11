from ogtk.utils import db
import matplotlib.pyplot as plt
import seaborn as sns
import ogtk.utils as ut
import ogtk.ltr.shltr as shltr
import scanpy as sc
from typing import Sequence,Optional,Iterable
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
    pseudo_counts=0,
    clone_numerator='h1',
    clone_denominator='g10',
    expr = None,
    ):
    ''' cell size vs clone composition 
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    cn = clone_numerator
    cd = clone_denominator
    a = pseudo_counts
    
    data = (pl.DataFrame(exp.ad_sc.obs)
            .with_columns(pl.col('total_counts').fill_null(0))
            .with_columns(((a+pl.col(cn))/(a+pl.col(cd))).log10().alias('log_ratio'))
            )

    if expr is not None:
        data = data.filter(expr)

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
        self.raw_adata = None
        self.explored = False

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

        adata.obs['sample_id'] = batch
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
            min_reads = 2,
            valid_ibars = None,
            down_sample = None,
            use_cache = True,
            corr_dir_fn=None):
        ''' Returns data frame at the molecule level

            Invokes ``shltr.sc.reads_to_molecules``
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
            min_reads=min_reads, 
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
            Annotates the ibars using a ibar 'cluster' table
            # TODO : improve the ibar annotation functionality
        '''

        self.mols = self.load_guide_molecules(valid_ibars = valid_ibars, clone = clone)
        self.mols = self.annotate_ibars(ibar_ann=ibar_ann)
        self.score_clone(min_cell_size=min_cell_size)
        

    def annotate_ibars(self, 
                  ibar_ann: str = '/local/users/polivar/src/artnilet/workdir/scv2/all_parquetino_clusters',
       ):

        dfc = (pl.scan_parquet(ibar_ann)
               .drop('sample_id').unique()
                #.filter(pl.col('sample_id')==sample_id)
                #.filter(pl.col('cluster')!=0)
                .collect()
            )

        dfc = (self.mols.join(
                    dfc, 
                    left_on =['raw_ibar'], 
                    right_on=['raw_ibar'], 
                    how='left')
                    )
        return(dfc)

    def init_adata(self, force:bool=True):
        ''' Loads the raw cell ranger counts matrix
        '''
        if self.raw_adata is None or force:
            self.raw_adata = self.load_10_mtx()
        #self.adata.obs['bath'] = self.sample_id

    @wraps(db.run_bcl2fq)
    def demux(self, *args, **kwargs):
        ''' demultiplex
        '''
        db.run_bcl2fq(self, *args, **kwargs)
        
    @wraps(db.run_cranger)
    def cellranger(self, *args, **kwargs):
        ''' cell ranger wrapper
        '''
        db.run_cranger(self, *args, **kwargs)

    @wraps(shltr.sc.compute_clonal_composition)
    def score_clone(self, min_cell_size=0, *args, **kwargs):
        ''' Compute clonal composition 

            1. Invoques ``shltr.sc.compute_clonal_composition`` applying the
            ``min_cell_size`` argument which determines the size of a cell in
            total umi counts from the targeted library and not the gene
            expression.

            2. Merges the clonal composition data frame into the single-cell anndata object. 

        '''

        # when computing clonal composition directly merge this information to the ad_sc.obs 
        fn = shltr.sc.compute_clonal_composition

        cc = fn(self.mols.filter(pl.col('umis_cell')>min_cell_size), *args, **kwargs)
        #cc = cc.with_columns(pl.col('total_counts').fill_null(0))
        self.cc = cc.clone()

        cc = cc.with_columns((pl.col('cbc') + '-1')).rename({'cbc':'index'}).to_pandas().set_index('index')

        self.ad_sc.obs = (
                self.ad_sc.obs.merge(
                    right=cc,
                    left_index=True,
                    right_index=True,
                    how='left')
        )

    def init_object(self, sample_name: str| None= None):
        ''' Main function to load full experiment object
        '''
        if sample_name is None:
            sample_name = self.sample_id

        self.init_workdir()
        db.run_cranger(self, localcores=75, uiport=7777)

        res = self.load_final_ad(sample_name=sample_name)
        if res >0:
            return

        db.tabulate_xp(self)

        self.init_mols(valid_ibars = self.default_ibars(), clone=self.clone)
        #self.init_adata()

    def load_final_ad(self, sample_name):
        ''' Loads final version of mcs and scs adatas
        ''' 
        mcs_fad_path = self.return_path('mcs_ad_path')
        scs_fad_path = self.return_path('scs_ad_path')
        
        if os.path.exists(mcs_fad_path) and os.path.exists(scs_fad_path):
            self.print(f'loading final adatas:\n{mcs_fad_path}\n{scs_fad_path}', 'bold white')
            self.ad_mc = ad.read_h5ad(mcs_fad_path)
            self.ad_sc = ad.read_h5ad(scs_fad_path)
            return 0

        else:
            self.print(f'No final adatas found, compute them manually using .do_mc(explore=True), select best parameters and re-run .do_mc(explore=False)', 'bold red', force = True)
            return 1
            raise ValueError 

    def save_final_ad(self, sample_name):
        '''
        '''
        mcs_fad_path = self.return_path('mcs_ad_path')
        scs_fad_path = self.return_path('scs_ad_path')

        self.print(f'saving final adatas:\n{mcs_fad_path}\n{scs_fad_path}', 'bold cyan')

        self.ad_mc.write_h5ad(mcs_fad_path)
        self.ad_sc.write_h5ad(scs_fad_path)

    def return_path(self, kind, sample_name:str|None=None):
        if sample_name is None:
            sample_name = self.sample_id
        #_f for final
        mcs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells.h5ad'
        mcs_fad_path = mcs_ad_path.replace('mcells', 'mcells_f')
        
        scs_ad_path = mcs_ad_path.replace('mcells', 'scells')
        scs_fad_path = mcs_fad_path.replace('mcells', 'scells')
        
        otl_ad_path = mcs_ad_path.replace('mcells', 'outliers')
        
        return locals()[kind]
    
    def clear_mc(self, 
                 kinds:Sequence|None=None,
                ):
        if kinds is None:
            kinds = ['mcs_ad_path', 'scs_ad_path', 'otl_ad_path', 'mcs_fad_path', 'scs_fad_path']
        for i in kinds:
            i = self.return_path(i)
            if os.path.exists(i):
                print(f'removing {i}')
                os.system(f'rm {i}') 
        
        
    @wraps(ut.sc.metacellize)
    def do_mc(
        self,
           sample_name:str|None= None,
           lateral_mods:Sequence|None=None,
           force:bool=False,
           explore:bool=True,
           full_cpus:int=8,
           target_metacell_size:int=100,
           *args,
           **kwargs,
           ):
        if sample_name is None:
            sample_name = self.sample_id

        mcs_ad_path = self.return_path('mcs_ad_path')
        scs_ad_path = self.return_path('scs_ad_path')
        otl_ad_path = self.return_path('otl_ad_path')

        if os.path.exists(scs_ad_path) and os.path.exists(mcs_ad_path) and not force:
            self.print('loading cached adatas')
            mcs = ad.read_h5ad(mcs_ad_path)
            scs = ad.read_h5ad(scs_ad_path)
            outliers = ad.read_h5ad(otl_ad_path)
            metadata = None
                                    
        else:
            if full_cpus is not None:
                kwargs['cpus']={'full':full_cpus, 'moderate':16}

            self.init_adata(force)
            res = ut.sc.metacellize(
                    set_name=sample_name,
                    adata=self.raw_adata, 
                    adata_workdir=self.wd_scrna,
                    explore = explore,
                    *args,
                    **kwargs
                  )

            if not explore:
                scs, mcs, metadata = res
            else:
                self.explored = True
                print(f'{explore=}')
                return res
        # cluster metacells
        #mcs = ut.sc.scanpyfi(mcs.copy())
        # cluster single-cells
        #scs = ut.sc.scanpyfi(scs.copy())

        
        # propagate metacell results to single-cells
        #scsm = scs.obs.merge(
        #        right=mcs.obs.reset_index(drop=True), 
        #        left_on='metacell', 
        #        right_index=True, 
        #        how='outer', 
        #        suffixes=['', '_mc'])

        # rigth way to overwrite .obs ?
        #scs.obs = scsm.loc[scs.obs.index]
        self.ad_sc = scs
        self.ad_mc = mcs
        self.metadata = metadata
       
        self.save_final_ad(sample_name)

    @wraps(shltr.sc.allele_calling)
    def init_alleles(self, 
                     expr: None | pl.Expr=None,
                     valid_ibars: Sequence| None= None,
               *args, **kwargs):
        ''' 
        '''
        if valid_ibars is None:
            valid_ibars = self.default_ibars()

        #ax = sns.histplot(self.mols['db_norm_umis_allele'].log10(),
        #                  element='step', fill=False)
        if expr is not None:
            plt.figure()
         #   ars = sns.histplot(self.mols.filter(expr)['db_norm_umis_allele'].log10(), 
         #                      ax = ax, color='orange', fill=True,
         #                     element='step')
            self.alleles = shltr.sc.allele_calling(self.mols.filter(expr), *args, **kwargs)
        else:
            self.alleles = shltr.sc.allele_calling(self.mols, *args, **kwargs)

        self.alleles = (self.alleles
                        .with_column(pl.col('raw_ibar').is_in(valid_ibars)
                            .alias('valid_ibar'))
                        .sort(['cluster', 'raw_ibar'])
                        )

        self.alleles = self.alleles.join(ut.sc.adobs_pd_to_df(self.ad_sc), 
                                 left_on='cbc', 
                                 right_on='cbc', 
                                 how='left')

        self.alleles = self.alleles.filter(pl.col('umis_cell').is_not_null())


    @wraps(shltr.sc.to_matlin)
    def init_matlin(self, cells= 100, cores=4, *args, **kwargs):

        plt.rcParams['figure.dpi'] = 100
        fn = shltr.sc.to_matlin
        self.matl = fn(self.alleles, cells=cells, *args, **kwargs)
        self.matl.plot_mat(rows=range(0, self.matl.df.shape[0]))
        plt.figure()
        self.matl.allele_distance(cores=cores)    
        self.matl.hclust(optimal_ordering=False)
        self.matl.plot_mat(rows=range(0, self.matl.df.shape[0]))
        plt.show()
        return()

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

    def ingest_xps(self, xps: Iterable):
        '''

        '''
        adatas = [i.ad_sc for i in xps]
        self.ad_sc = ad.AnnData.concatenate(*adatas,
                                 join='outer', 
                                 batch_key='batch_id')
#
        self.ad_sc = self.ad_sc[~self.ad_sc.obs['excluded_cell'], ].copy()
        self.ad_sc.var = self.ad_sc.var.drop(self.ad_sc.var.columns[2:].to_list(), axis='columns')
        #self.ad_sc.obs.groupby(['sample_id', 'batch_id']).head(1).reset_index().loc[:, ['sample_id', 'batch_id']]
        





