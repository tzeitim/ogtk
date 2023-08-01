from ogtk.utils import db
import matplotlib.pyplot as plt
import seaborn as sns
import ogtk.utils as ut
import ogtk.ltr.shltr as shltr
import numpy as np
import scanpy as sc
from typing import Sequence,Optional,Iterable
import polars as pl
import anndata as ad
import os
from functools import wraps


def _cassit(exp, alleles, clusters= [-1, 1, 2],
    allele_rep_thresh = 0.9,
           ):
    '''
    '''
    import cassiopeia as cas
    import pandas as pd
    import tempfile
    import os
    import numpy as np

    cas_dict = {
            'cbc':'cellBC',
            'raw_ibar':'intBC',
            'seq':'r1',
            'cluster':'LineageGroup',
            'sample_id':'sampleID',
            'umi_reads':'readCount',
            'umis_allele':'UMI',
           }
    
    pallele_table = (alleles
         #   .filter(pl.col('valid_ibar'))
            .drop('valid_ibar')
            .rename(cas_dict)
            .filter(pl.col('UMI')>1)
           )

    allele_table = (
                    pallele_table
                    .with_columns(
                        pl.when(pl.col('wt'))
                        .then('NONE')
                        .otherwise(pl.col('r1'))
                        .alias('r1'))
                    .to_pandas()
                    )

    indel_priors = cas.pp.compute_empirical_indel_priors(
            allele_table, 
            grouping_variables=['intBC', 'LineageGroup'])

    clone_allele_table = allele_table[allele_table['LineageGroup'].isin(clusters)]

    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
            clone_allele_table,
            allele_rep_thresh = allele_rep_thresh,
            mutation_priors = indel_priors) 
    
    cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)

    agg_dict ={"intBC": 'nunique', 
               'UMI': 'sum', 
               'sampleID': 'unique', 
               'metacell_name':'unique'} 
    cell_meta = clone_allele_table.groupby('cellBC').agg(agg_dict)
    cell_meta['sampleID'] = [x[0] for x in cell_meta['sampleID']]
    cell_meta['metacell_name'] = [x[0] for x in cell_meta['metacell_name']]


    missing_proportion = (character_matrix == -1).sum(axis=0) / character_matrix.shape[0]
    uncut_proportion = (character_matrix == 0).sum(axis=0) / character_matrix.shape[0]
    n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis=0)

    character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 'uncut_prop', 'n_unique_states']).T

    cas_tree.cell_meta = cell_meta
    cas_tree.character_meta = character_meta
    
    vanilla = True
    if vanilla:
        # create a basic vanilla greedy solver
        vanilla_greedy = cas.solver.VanillaGreedySolver()

        # reconstruct the tree
        vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)

    exp.tree = cas_tree
    exp.tree.clone_allele_table = clone_allele_table
    exp.tree.indel_to_char = state_2_indel

def _cas_return_colors(clone_allele_table, csp='hsv', do_plot=False):
    '''
    '''
    from colorhash import ColorHash
    import colorsys
    import numpy as np

    #clone_allele_table['rgb'] = clone_allele_table['r1'].map(lambda x: ColorHash(x).rgb)
    # simple capture of color hash attributes
    # with a normalization for rgb values
    clone_allele_table['hex'] = clone_allele_table['r1'].map(lambda x: ColorHash(x).hex)
    clone_allele_table['hsl'] = clone_allele_table['r1'].map(lambda x: ColorHash(x).hsl)
    clone_allele_table['rgb'] = clone_allele_table['r1'].map(lambda x: tuple(np.array(ColorHash(x).rgb)/255))
    # further conversion to hsv via colorsys
    clone_allele_table['hsv'] = clone_allele_table['rgb'].map(lambda x: colorsys.rgb_to_hsv(*x))

    #clone_allele_table['rgb'] = clone_allele_table['r1'].map(lambda x: colorsys.hls_to_rgb(*ColorHash(x).hsl))

    selected_columns =['r1', 'hex', 'hsl', 'hsv', 'rgb'] 
    colors = clone_allele_table.loc[:, selected_columns].drop_duplicates().set_index('r1')

    if do_plot:
        # control plot where the first element of the color triplet is used as Y axis
        # coloring by the selected color space (csp)
        x=range(0, colors.shape[0])
        y=colors[csp].map(lambda x: x[0])
        plt.scatter(x=x, y=y, c=colors[csp], s=55)

    colors = colors.rename(columns={csp:'color'}).loc[:, ['color']]
    return colors

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
    
    data = (pl.DataFrame(exp.scs.obs)
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.conf_keys = [i for i in vars(self)]
        self.explored = False
        
        self.raw_adata = None
        self.batch_map = None
        self.attr_d= {'shrna':'mols', 'zshrna':'zmols'}
        #self.is_chimera = 'is_chimera' in vars(self).keys()


    def init_wd(self):
        '''
            Generally useful for pooled Xps, e.g. three 10x lanes from the sample original sample.
            It iterates through defined ``wd_*`` attributes and creates the corresponding directories
        '''
        for i in [getattr(self, i) for i in vars(self).keys() if i.startswith('wd_')]:
            self.print(f"creating {i}")
            os.system(f'mkdir -p {i}')

    def load_10_h5(self, h5_path, batch = None, *args, **kwargs):
        ''' Loads cellranger's h5 output file and annotates ``batch`` as a field in .obs
        '''
        self.print(f"importing h5 from:\n {h5_path}")

        adata = sc.read_10x_h5(filename=h5_path, *args, **kwargs)

        if batch is None:
            batch = self.sample_id

        adata.obs['sample_id'] = batch
        adata.var_names_make_unique()
        return(adata)


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
        if not isinstance(files, list):
            files=[files]
        return(files)

    def load_guide_molecules(
            self,
            clone,
            sample_id = None,
            suffix = 'shrna',
            pattern: str|None=None,
            min_reads = 2,
            valid_ibars = None,
            downsample = None,
            use_cache = True,
            corr_dir_fn=None,
            filter_valid_cells=False):
        ''' Returns data frame at the molecule level

            Invokes ``shltr.sc.reads_to_molecules``
            
            When more than one parquet file is found by the function, all of them are loaded into a single data frame 
        '''
        if sample_id is None:
            sample_id = self.sample_id

        # looks for tabulated parquet files
        parquet_ifn = self.list_guide_tables(suffix)

        if pattern is not None:
            parquet_ifn = [i for i in parquet_ifn if pattern in i]

        self.print('Reading molecule parquets')
        self.print(parquet_ifn)

        df = shltr.sc.reads_to_molecules(
                sample_id=sample_id,
                parquet_ifn=parquet_ifn,
                valid_ibars=valid_ibars, 
                use_cache=use_cache,
                cache_dir=self.return_cache_path('mols', suffix),
                corr_dict_fn=corr_dir_fn,
                downsample=downsample,
                #min_cells_per_ibar=int(len(obs27['cbc']) * 0.2),
                min_reads=min_reads, 
                max_reads=1e6, 
                clone=clone)

        if filter_valid_cells:
            self.load_final_ad()
            df = df.filter(pl.col('cbc').is_in(self.scs.obs.index.to_list()))

        return(df)

    def init_mols(self,
                  valid_ibars: Sequence,
                  clone: str, 
                  suffix: str='shrna',
                  ibar_ann: str|None=None,
                  min_cell_size: int=0,
                  min_reads: int=2,
                  downsample: str|None=None,
                  pattern: str|None=None,
                  sample_id: str|None=None,
                  ):
        ''' Creates the .mols attribute
            Annotates the ibars using a ibar 'cluster' table
            # TODO : improve the ibar annotation functionality
        '''

        if sample_id is None:
            sample_id = self.sample_id

        if downsample is not None:
            self.print(f'Downsampling to {downsample} reads', style='bold #ff0000')
            downsample = int(downsample)

        mols = self.load_guide_molecules(valid_ibars=valid_ibars, 
                                         min_reads=min_reads, 
                                         clone=clone, 
                                         suffix=suffix, 
                                         downsample=downsample, 
                                         pattern=pattern)

        # annotate .mols
        mols = self.annotate_ibars(mols=mols, ibar_ann=ibar_ann)

        # store .mols
        setattr(self, self.attr_d[suffix], mols)

        if self.is_chimera:
            self.score_clone(min_cell_size=min_cell_size)

    def annotate_ibars(
            self, 
            mols,
            ibar_ann: str|None=None,
            ):

        ibar_ann = '/local/users/polivar/src/artnilet/workdir/scv2/ibar_clusters.parquet' if ibar_ann is None else ibar_ann

        dfc = (pl.scan_parquet(ibar_ann)
                .collect()
            )

        dfc = (mols.join(
                    dfc, 
                    left_on=['raw_ibar'], 
                    right_on=['raw_ibar'], 
                    how='left')
                    )
        return(dfc)

    def init_adata(self, force:bool=True, h5_path: str|None= None, cellbender: bool=False, *args, **kwargs):
        ''' Loads the raw cell ranger counts matrix
            When ``h5_path`` is provided it overrides the default routine that uses the matrix dir tree from 10x
        '''
        if h5_path is not None:
            self.raw_adata = self.load_10_h5(h5_path, *args, **kwargs)
            if cellbender:
                print('getting filtered cells from cellbender')
                sc.pp.calculate_qc_metrics(self.raw_adata, inplace=True)
                self.raw_adata = self.raw_adata[self.raw_adata.obs.total_counts>0,: ].copy()
                self.raw_adata.X = self.raw_adata.X.astype(np.float32)
            return None


        if self.raw_adata is None or force:
            self.raw_adata = self.load_10_mtx()
            return None




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

        # when computing clonal composition directly merge this information to the scs.obs 
        fn = shltr.sc.compute_clonal_composition

        cc = fn(self.mols.filter(pl.col('umis_cell')>min_cell_size), *args, **kwargs)
        #cc = cc.with_columns(pl.col('total_counts').fill_null(0))
        self.cc = cc.clone()

        # TODO clean the mess regarding '-1'
        #cc = cc.with_columns((pl.col('cbc') + '-1')).rename({'cbc':'index'}).to_pandas().set_index('index')
        cc = cc.rename({'cbc':'index'}).to_pandas().set_index('index')

        self.scs.obs = (
                self.scs.obs.merge(
                    right=cc,
                    left_index=True,
                    right_index=True,
                    how='left')
        )

    def init_object(self, sample_name: str| None= None, uiport=7777, skip_mols=False, min_reads_per_mol: int=2,):
        ''' Main function to load full experiment object
        '''
        if sample_name is None:
            sample_name = self.sample_id

        self.init_workdir()

        if 'pp' in vars(self):
            db.run_cranger(self, localcores=75, uiport=uiport)

        res = self.load_final_ad()

        if self.tabulate is not None:
            db.tabulate_xp(self)


        if not skip_mols:
            self.init_mols(valid_ibars = self.default_ibars(), clone=self.clone, min_reads=min_reads_per_mol)

    def load_final_ad(self):
        ''' Loads final version of mcs and scs adatas
        ''' 
        mcs_fad_path = self.return_path('mcs_ad_path')
        scs_fad_path = self.return_path('scs_ad_path')
        
        if os.path.exists(mcs_fad_path) and os.path.exists(scs_fad_path):
            self.print(f'loading final adatas:\n{mcs_fad_path}\n{scs_fad_path}', 'bold white')
            self.mcs = ad.read_h5ad(mcs_fad_path)
            self.scs = ad.read_h5ad(scs_fad_path)
            return 0

        else:
            self.print(f'No final adatas found, compute them manually using .do_mc(explore=True), select best parameters and re-run .do_mc(explore=False)', 'bold #ff0000', force = True)
            return 1
            raise ValueError 

    def save_final_ad(self, sample_name):
        '''
        '''
        mcs_fad_path = self.return_path('mcs_ad_path')
        scs_fad_path = self.return_path('scs_ad_path')

        self.print(f'saving final adatas:\n{mcs_fad_path}\n{scs_fad_path}', 'bold cyan')

        self.mcs.write_h5ad(mcs_fad_path)
        self.scs.write_h5ad(scs_fad_path)

    def return_path(self, kind, sample_name:str|None=None, suffix=None):
        if sample_name is None:
            sample_name = self.sample_id
        #_f for final
        mcs_ad_path = f'{self.wd_scrna}/{sample_name}.mcells.h5ad'
        mcs_fad_path = mcs_ad_path.replace('mcells', 'mcells_f')
        
        scs_ad_path = mcs_ad_path.replace('mcells', 'scells')
        scs_fad_path = mcs_fad_path.replace('mcells', 'scells')
        
        otl_ad_path = mcs_ad_path.replace('mcells', 'outliers')
        
        if suffix is not None:
            mols = self.return_cache_path('mols', suffix=suffix)
            mols = f'{mols}/{sample_name}_r2mols.parquet'

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
        use_cache:bool=True,
        explore:bool=True,
        full_cpus:int=8,
        moderate_cpus: int=8,
        target_metacell_size:int=100,
        max_parallel_piles=None,
        adata=None,
        *args,
        **kwargs,
           ):

        if sample_name is None:
            sample_name = self.sample_id

        mcs_ad_path = self.return_path('mcs_ad_path')
        scs_ad_path = self.return_path('scs_ad_path')
        otl_ad_path = self.return_path('otl_ad_path')

        if os.path.exists(scs_ad_path) and os.path.exists(mcs_ad_path) and use_cache:
            self.print('loading cached adatas')
            mcs = ad.read_h5ad(mcs_ad_path)
            scs = ad.read_h5ad(scs_ad_path)
            outliers = ad.read_h5ad(otl_ad_path)
            metadata = None
                                    
        else:
            if full_cpus is not None:
                kwargs['cpus']={'full':full_cpus, 'moderate':moderate_cpus}

            if adata is None and self.raw_adata is None:
                self.init_adata(force)

            res = ut.sc.metacellize(
                    set_name=sample_name,
                    adata=adata if adata is not None else self.raw_adata, 
                    adata_workdir=self.wd_scrna,
                    explore = explore,
                    force = force,
                    max_parallel_piles=max_parallel_piles,
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
        self.scs = scs
        self.mcs = mcs
        self.metadata = metadata
       
        self.save_final_ad(sample_name)

    @wraps(shltr.sc.allele_calling)
    def init_alleles(self, 
                     suffix='shrna',
                     expr:None|pl.Expr=None,
                     valid_ibars:Sequence|None=None,
                     *args,
                     **kwargs):
        ''' 
        '''
        if valid_ibars is None:
            valid_ibars = self.default_ibars()

        mols = getattr(self, self.attr_d[suffix] )

        if expr is not None:
            alleles = shltr.sc.allele_calling(mols.filter(expr), *args, **kwargs)
        else:
            alleles = shltr.sc.allele_calling(mols, *args, **kwargs)
        
        # add mask for valid ibars
        alleles = (alleles
                        .with_columns(pl.col('raw_ibar').is_in(valid_ibars)
                            .alias('valid_ibar'))
                        .sort(['cluster', 'raw_ibar'])
                        )

        # merge with scs.obs
        # TODO a bug lurches here
        strip_pattern = '(.+?)-(.)' if self.batch_map is not None else '(.+)'
        strip_pattern = '(.+)'
        self.print(f'{strip_pattern=}', style='bold #00ff00')

        # the majority of the raw cell barcodes do not belong to any real cell so we constrain the merge to the scs
        # how='outer'
        alleles = alleles.join(
                ut.sc.adobs_pd_to_df(self.scs, strip_pattern), 
                left_on='cbc', 
                right_on='cbc', 
                how='inner') 

        # since polars lacks a 'how=right' merging, it is needed to drop nulls
        # in a field that is meaningful for the alleles, such as 'umis_cell'
        #alleles = alleles.filter(pl.col('umis_cell').is_not_null())

        setattr(self, suffix[0]+'alleles', alleles)
        if self.is_chimera:
            # TODO clean this block, e.g. log2 or 10?
            import numpy as np
            self.salleles = self.salleles.with_columns((pl.col('h1')/pl.col('g10')).log(base=10).alias('cs'))
            self.salleles = self.salleles.with_columns(pl.col(pl.Categorical).cast(pl.Utf8))

            self.scs.obs['cs'] = np.log10(self.scs.obs['g10']/self.scs.obs['h1'])


    @wraps(shltr.sc.to_matlin)
    def init_matlin(self, do_cluster = False, subset=['cluster', 'raw_ibar'], cells= 100, cores=4, *args, **kwargs):

        plt.rcParams['figure.dpi'] = 100
        fn = shltr.sc.to_matlin
        self.matl = fn(self.alleles, cells=cells, subset=subset,  *args, **kwargs)

        self.matl.plot_mat(rows=range(0, self.matl.df.shape[0]))
        plt.figure()

        if do_cluster:
            self.matl.allele_distance(cores=cores)    
            self.matl.hclust(optimal_ordering=False)
            self.matl.plot_mat(rows=range(0, self.matl.df.shape[0]))
            plt.show()
        return()

    def default_ibars(self):
            '''
                Loads pre-determined list of valid ibars
            '''
            self.print(':red_square: :red_square: Loading pre-computed valid ibars :red_square: :red_square:', 'bold #ff0000')
            valid_ibars = pl.read_csv('/local/users/polivar/src/artnilet/conf/valid_ibars.csv')['valid_ibar'].to_list()      
            return(valid_ibars)

    @wraps(_plot_cc)
    def plot_cc(self, *args, **kwargs):
        _plot_cc(self, *args, **kwargs)

    def ingest_xps(self, xps: Iterable, suffix='shrna', force=False, skip_mols=False):
        ''' Merge compatible experiments and integrate them into an invidiual one.
            a ``batch_map`` keeps track of an appended identifier to single-cell barcodes
        '''
        if self.mols is not None and self.raw_adata is not None:
            self.print('This experiment seems to be have ingested others. No need to ingest unless force=True')
            return None



        adatas = [i.scs for i in xps]
        
        adata = ad.concat(adatas, join='outer', label='batch_id', index_unique="_")
        adata = adata[~adata.obs['excluded_cell'], ].copy()
        adata.var = adata.var.drop(adata.var.columns[2:].to_list(), axis='columns')
        adata.obs = adata.obs[['sample_id', 'batch_id']].copy()
        
        batch_map = adata.obs.set_index('sample_id')['batch_id'].to_dict()

        #batch_map = (
        #        adata.obs
        #        .groupby(['sample_id', 'batch_id'])
        #        .head(1)
        #        .reset_index()
        #        .loc[:, ['sample_id', 'batch_id']]
        #        )

        #batch_map = dict(list(zip(batch_map.sample_id, batch_map.batch_id)))

        self.batch_map = batch_map
        assert len(batch_map) == len(adata.obs.sample_id.unique()), "It seems that one or more batches are duplicated"

        self.raw_adata = adata
        self.batch_map = batch_map
        self.conf_keys.append('batch_map')

        if skip_mols:
            self.export_xpconf(xp_conf_keys = set(self.conf_keys))
            return None

        mols = (
           pl.concat([i.load_guide_molecules(clone=i.clone, filter_valid_cells=True) for i in xps])
        )
        #batch_dict = .scs.obs.set_index('sample_id')['batch_id'].to_dict()
        self.mols = (
                mols
                .with_columns(pl.col('sample_id')
                .map_dict(batch_map)
                .alias('batch_id'))        
        )

        ## merge molecules from individual experiments  
        #def pl_exp(xp):
        #    return xp.mols.with_columns(pl.col('cbc')+"_"+batch_map[xp.mols['sample_id'].unique()[0]])
        #
        #self.mols = pl.concat(
        #    [ pl_exp(xp) for xp in xps ]
        #)

        # save parquet file
        self.print(f"saving molecules to {self.return_path('mols', suffix=suffix)}")
        self.mols.write_parquet(self.return_path('mols', suffix=suffix))

        self.export_xpconf(xp_conf_keys = set(self.conf_keys))

    def return_cache_path(self, key, suffix=None):
        '''
        '''
        if key =='mols':
            assert suffix is not None, "A suffix for the type of molecule is needed, e.g., 'shrna', 'zhrna'"
            return f'{self.wd_samplewd}/{suffix}/'

    @wraps(_cassit)
    def cassit(self, *arg, **kwargs):
        _cassit(self, *arg, **kwargs)



@wraps(db.print_template)
def print_template(*args, **kwargs):
    db.print_template(*args, **kwargs)
