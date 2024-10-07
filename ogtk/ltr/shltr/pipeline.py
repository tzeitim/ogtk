import ogtk
from ogtk.utils import db
from ogtk.ltr.shltr import cluster
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import ogtk.utils as ut
import ogtk.ltr.shltr as shltr
import numpy as np
import scanpy as sc
import subprocess
from typing import Sequence,Optional,Iterable,List
import polars as pl
import metacells as mc
import anndata as ad
import os
import time
from functools import wraps
from pyaml import yaml

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()

def _check_required_attributes(self, case):
    """
    Check if the required attributes for a specific case are present in the instance.

    Args:
        case (str): The case for which to check the required attributes. Supported cases
                    include 'scvi_solo', 'cellbender', 'scvi'.

    Raises:
        AttributeError: If any required attribute for the specified case is missing.
    """

    attribute_sets = {
        "scvi_solo": ["path_scvi_solo_sh", "path_scvi_solo_sge", "wd_sge"],
        "cellbender": ["cellbender_attribute_1", "cellbender_attribute_2"],
        "scvi": ["scvi_attribute_1", "scvi_attribute_2"]
    }

    # Select the required attributes based on the case
    required_attrs = attribute_sets.get(case)

    if not required_attrs:
        raise ValueError(f"Unsupported case '{case}'. Supported cases are: {', '.join(attribute_sets.keys())}")

    for attr in required_attrs:
        if not hasattr(self, attr):
            raise AttributeError(f"experiment not configured for '{case}'. Add '{attr}' to the experiment's template")
    return True


def _cassit(exp, 
            alleles, 
            clusters=List, 
            allele_rep_thresh = 1.0,
            solver:str='nj',
            threads=66,
            collapse_mutationless_edges=False,
) -> None:
    '''
    - expects a filtered allele table, e.g. no further decisions are made at the allele level
    '''
    import pandas as pd
    import tempfile

    cas_dict = {
            'cbc':'cellBC',
            'raw_ibar':'intBC',
            'seq':'r1',
            'cluster':'LineageGroup',
            'sample_id':'sampleID',
            'umi_reads':'readCount',
            'umis_allele':'UMI',
           }

    oalleles = alleles.clone()

    if 'cellBC' in alleles:
        alleles = alleles.drop('cellBC')

    pallele_table = (alleles
            .rename(cas_dict)
   )

    allele_table = (
            pallele_table
            .with_columns(
                pl.when(pl.col('wt'))
                .then(pl.lit('NONE'))
                .otherwise(pl.col('r1'))
                .alias('r1'))
            .to_pandas()
    )
    import cassiopeia as cas


    intree = set(allele_table.reset_index()['cellBC'].unique())


    indel_priors = cas.pp.compute_empirical_indel_priors(
            allele_table, 
            grouping_variables=['intBC', 'LineageGroup'])

    clone_allele_table = allele_table[allele_table['LineageGroup'].isin(clusters)]

    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
            clone_allele_table,
            allele_rep_thresh = allele_rep_thresh,
            mutation_priors = indel_priors) 
    
    intree2 = set(character_matrix.reset_index()['index'].unique())

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
    

    if solver == 'vanilla':
        # create a basic vanilla greedy solver
        vanilla_greedy = cas.solver.VanillaGreedySolver()
        # reconstruct the tree
        vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)

    if solver == 'ilp':
        ilp_solver = cas.solver.ILPSolver(convergence_time_limit=500, 
                                        maximum_potential_graph_layer_size=500, 
                                        weighted=True, 
                                        seed=1234)
        ilp_solver.solve(cas_tree)

    if solver == 'hybrid':
        vanilla_greedy = cas.solver.VanillaGreedySolver()

        ilp_solver = cas.solver.ILPSolver(convergence_time_limit=500, 
                                          maximum_potential_graph_layer_size=500, 
                                          weighted=True, 
                                          seed=1234)

        hybrid_solver = cas.solver.HybridSolver(top_solver=vanilla_greedy, 
                                                bottom_solver=ilp_solver, 
                                                cell_cutoff=40, 
                                                threads=threads)
        hybrid_solver.solve(cas_tree, logfile=f'{exp.path_trees}/example_hybrid.log')

    if solver == 'nj':
        nj_solver = cas.solver.NeighborJoiningSolver(
                dissimilarity_function=cas.solver.dissimilarity.weighted_hamming_distance, 
                add_root=True)
        nj_solver.solve(cas_tree, collapse_mutationless_edges=collapse_mutationless_edges)

    exp.tree = cas_tree
    exp.tree.clone_allele_table = clone_allele_table
    exp.tree.indel_to_char = state_2_indel

    # annotate cell metadata with single-cell data from anndata
    # this might include the metacell already

    column_iset = list(set(exp.tree.cell_meta.columns).difference(exp.scs.obs.columns))


    exp.tree.cell_meta = (
            exp.tree.cell_meta.loc[:, column_iset].merge(
                exp.scs.obs,
                left_index=True,
                right_index=True, 
                how='left')
            )

    if exp.mc_ann is not None:
        _pd_annotate_mc(exp.tree.cell_meta, mc_ann=exp.mc_ann)


def _cas_update_mc_ann(exp, mc_ann: pl.DataFrame):
    ''' Annotates a cassiopeia tree cell metadata with a meta cell annotation from MCview
    '''
    exp.load_mc_ann(mc_ann)
    exp.tree.cell_meta = _pd_annotate_mc(exp.tree.cell_meta, mc_ann=exp.mc_ann)
    

def _pd_annotate_mc(df ,
                    mc_ann: pl.DataFrame, 
                    left_on='metacell_name', 
                    right_on='metacell'
                    ):
    ''' Annotates a cassiopeia tree cell metadata with a meta cell annotation from MCview
    '''
    column_iset = list(set(df.columns).intersection(mc_ann.columns))
    index_name = df.index.name
    
    if left_on in column_iset:
        column_iset.remove(left_on)
        
    df = df.loc[:, df.columns.difference(column_iset)]

    df =(
            df.reset_index()
            .merge(
                mc_ann.to_pandas(),
                left_on=left_on,
                right_on=right_on,
                how='left')
            )

    df['cell_type'] = df['cell_type'].astype(str)

    if index_name not in df.columns:
        df = df.set_index('index')
    else:
        df = df.set_index(index_name)

    return df

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
    lim= (-2.5 , 2.5),
    cell_size_metric = "n_genes_by_counts",
    vmax=100,
    cluster_suffix='',
    pseudo_counts=1,
    clone_numerator='ibc_1',
    clone_denominator='ibc_2',
    expr = None,
    ):
    ''' cell size vs clone composition 
    This function is used to explore the relationship between cell size and clone composition.
    the expr argument can be passed to filter the data frame before plotting
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    if 'total_counts' not in exp.scs.obs.columns:
        sc.pp.calculate_qc_metrics(exp.scs, inplace=True)

    cn = clone_numerator
    cd = clone_denominator
    a = pseudo_counts
    
    data = (pl.DataFrame(exp.scs.obs)
            .with_columns(pl.col('total_counts').fill_null(0))
            .with_columns(((a+pl.col(cn))/(a+pl.col(cd))).log(base=2).alias('log_ratio'))
            )

    if expr is not None:
        data = data.filter(expr)

    data = data.filter(pl.col('total_counts')<max_cell_size)
    data = data.filter(pl.col('umis_cell')<tmax_cell_size)
    data = data.filter(pl.col('log_ratio')<=lim[1]).filter(pl.col('log_ratio')>=lim[0])

    fig, (ax, ax2, ax3) = plt.subplots(1,3, figsize=(3*15,15))
    sns.histplot(data=data.to_pandas(),
                 x=cell_size_metric,
                 y='log_ratio',
                 bins=100, ax= ax, cmap='RdYlBu_r', vmax=vmax)


    ax.set_ylim(lim)
    ax.grid()
    ax.set_ylabel(f'{clone_numerator}/{clone_denominator}')

    sns.histplot(data=data.to_pandas(),
                 x='umis_cell',
                 y='log_ratio',
                 bins=100, ax= ax2, cmap='RdYlBu_r', vmax=vmax)

    ax2.set_xlim((0, tmax_cell_size))
    ax2.set_ylim(lim)
    ax2.grid()

    ax.set_ylabel(f'{clone_numerator}/{clone_denominator}')

    sns.histplot(data=data.to_pandas(),
                x=cell_size_metric,
                 y='umis_cell', 
                 bins=100,
                 log_scale=(False, False),
                 ax = ax3, cmap='RdYlBu_r', vmax=vmax)
    ax3.grid()
    plt.show()
    plt.figure()

    if 'leiden' in exp.scs.obs:
        sns.displot(data=data.to_pandas(),
                    x=cell_size_metric,
                    y='log_ratio',
                    bins=100,
                    log_scale=(False, False),
                    col=f'leiden{cluster_suffix}', 
                    col_wrap=6, 
                    cmap='RdYlBu_r', 
                    vmax=50)
def _color_hist(data, vmax=None, vmin=None, bins=25, edge_line=0.77, ax=None, xlim=(-10, 10), cmap='bwr'):
    from matplotlib import colors as mcolors
    from matplotlib import colormaps as cmaps
    import matplotlib.pyplot as plt

    # Set vmax and vmin
    vmax = max(data) if vmax is None else vmax
    vmin = min(data) if vmin is None else vmin

    # Get the 'bwr' colormap
    cm = cmaps.get_cmap(cmap)

    # Create a Normalize object with your chosen vmin and vmax
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Create figure and axis if not provided
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(14,7))
        ax.set_xlim(xlim)


    # Plot histogram and get bins and patches
    n, bins, patches = ax.hist(data, bins, color='green')

    # Apply colormap and set border properties to each patch (bar)
    for bin, patch in zip(bins, patches):
        # Normalize bin value
        bin_norm = norm(bin)
        # Set color for each patch
        patch.set_facecolor(cm(bin_norm))
        # Set border color and thickness of the patch
        patch.set_edgecolor('black')
        patch.set_linewidth(edge_line)

    # Show the plot if the axis was not passed
    if ax is None:
        plt.show()

def _plot_alleles(alleles, *args, **kwargs):
    '''
    '''
    alleles = alleles
    g = sns.displot(
            data=alleles
            .filter(pl.col('umis_top_allele')<=5)    
            .filter(pl.col('ties')==0)
            .filter(pl.col('cluster')>0)
            .filter(pl.col('doublet_prediction')=='singlet')
            .sort('clone')
            .to_pandas(), 
            col='clone',
            row='umis_top_allele',
            x='norm_counts', 
            hue='cluster', 
            #multiple='stack',
            height=1.5,
            aspect=1.5,
            *args, **kwargs)

    return g
    #for ax in g.axes.flat:
    #    ax.set_xlim((4, 8))

    g = sns.displot(
            data=alleles
            .filter(pl.col('umis_top_allele')<=5)    
            .filter(pl.col('ties')==0)
            .filter(pl.col('cluster')>0)
            .filter(pl.col('doublet_prediction')=='singlet')
            .sort('clone')
            .to_pandas(), 
            row='umis_top_allele',
            x='norm_counts', 
            multiple='stack',
            height=1.5,
            aspect=1.5)

class Xp(db.Xp):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.conf_keys = [i for i in vars(self)]
        self.explored = False
        
        self.raw_adata = None
        self.scs: AnnData|None = None
        self.mcs = None
        self.batch_map = None
        self.mols = None
        self.attr_d= {'shrna':'mols', 'zshrna':'zmols'}
        self.ibar_ann=None
        self.cc=None

        if 'is_chimera' not in vars(self).keys():
            self.is_chimera = False
        # paths
        self.path_mcs_h5ad = f'{self.wd_mc}/{self.sample_id}.mcells.h5ad' # pyright: ignore
        self.path_raw_h5ad = f"{self.wd_mc}/{self.sample_id}.full.h5ad" # pyright: ignore
        self.path_cleansc_h5ad = self.return_path("clean_scs_ad_path")#f"{self.wd_mc}/{self.sample_id}.clean.h5ad" # pyright: ignore
        self.path_iterationx_h5ad = f"{self.wd_mc}/{self.sample_id}.iteration-XXX.h5ad" # pyright: ignore
        self.path_final_h5ad = self.path_iterationx_h5ad.replace('iteration-XXX', 'final')
        self.path_trees = f'{self.wd_cas}' # pyright: ignore
        self.path_cr_outs = f'{self.wd_scrna}/{self.sample_id}/outs/'

        # mc ann
        self.mc_colormap = None

        # lineage
        self.tree:None|cas.data.CassiopeiaTree.CassiopeiaTree = None

        # sge TODO remove my credentials
        self.sge_conf = ut.sge.SGE_CONF(user='polivar', host='max-login2')

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
        matrix_dir = f'{self.path_cr_outs}/filtered_feature_bc_matrix/'
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
        assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"

        path = f'{self.wd_xp}/{suffix}'
        files = ut.sfind(path=path, pattern='*raw_reads.parquet')
        if not isinstance(files, list):
            files=[files]
        return(files)

    @ogtk.utils.log.call
    def load_guide_mols(
            self,
            clone:str,
            sample_id:str|None = None,
            suffix:str = 'shrna',
            pattern:str|None=None,
            min_reads:int = 1,
            min_mols_per_ibar:int=1000,
            do_sampling = None,
            normalize=False,
            use_cache:bool = True,
            corr_dir_fn:str|None=None,
            filter_valid_cells=False):
        ''' Returns data frame at the molecule level

            Invokes ``shltr.sc.reads_to_molecules``
            
            When more than one parquet file is found by the function, all of them are loaded into a single data frame 
        '''
        if sample_id is None:
            sample_id = self.sample_id #type: ignore

        # looks for tabulated parquet files
        assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"
        parquet_ifn = self.list_guide_tables(suffix)

        if pattern is not None:
            parquet_ifn = [i for i in parquet_ifn if pattern in i]

        logger.debug('Reading molecule parquets')
        logger.debug(parquet_ifn)

        df = shltr.ibars.reads_to_molecules(
                sample_id=sample_id, #type: ignore
                parquet_ifn=parquet_ifn,#type: ignore
                zombie=suffix == 'zshrna',
                cache_dir=self.return_cache_path('mols', suffix), #type: ignore
                corr_dict_fn=corr_dir_fn,
                do_sampling=do_sampling,
                #min_cells_per_ibar=int(len(obs27['cbc']) * 0.2),
                min_reads=min_reads, 
                max_reads=int(1e6), #type: ignore
                clone=clone)
        logger.info(f'loaded mols {df.shape=}')

        if filter_valid_cells:
            if self.scs is None:
                raise ValueError("Please populate .scs before to filter valid cells")
            if 'excluded_cell' not in self.scs.obs.columns:
                raise ValueError("Expected field `excluded_cell` in .obs")
            #self.load_latest_ad()
            #df = df.filter(pl.col('cbc').is_in(self.scs.obs.index.to_list()))
            logger.info('filtering molecules for only valid cells')
            df = df.filter(pl.col('cbc').is_in(self.scs[~self.scs.obs.excluded_cell].obs.index.to_list()))
    
        logger.info(f'filtering ibars {min_mols_per_ibar=}')

        if normalize:
            df = self.normalize_mols(df)

        return(df.filter(pl.count().over('raw_ibar')>=min_mols_per_ibar))

    def normalize_mols(
            self,
            df:pl.DataFrame,
            min_umis_cell:int=10,
           over_cells:Sequence[str]=['doublet_prediction'],
           over_ibars:Sequence[str]=['cluster'],
           bg_cells_expr:pl.Expr=pl.col('doublet_prediction')=='lowq'
        )->pl.DataFrame:
        ''' Normalizes the molecule counts by the total number of molecules per cell
        In its current state it only works when 'doublet_prediction' is present in the .obs and mols
        '''
        #TODO: still need to smoothen the excluded_cells management

        if 'doublet_prediction' in df.columns:
            df=df.with_columns(pl.col('doublet_prediction').cast(pl.Utf8).fill_null("lowq"))
            logger.critical(f'using Expr:{bg_cells_expr} to normalize ibar counts')

            df=shltr.sc.normalize_mol_counts(
                    df=df.filter(pl.col('umi').count().over('cbc')>=min_umis_cell),
                    over_cells=over_cells,
                    over_ibars=over_ibars,
                    bg_cells_expr=bg_cells_expr)

        return df

    @ogtk.utils.log.call
    def init_mols(self,
                  clone: str, 
                  suffix: str='shrna',
                  ibar_ann:pl.DataFrame|None=None,
                  min_cell_size: int=0,
                  min_reads: int=1,
                  do_sampling: str|None=None,
                  pattern: str|None=None,
                  sample_id: str|None=None,
                  ):
        ''' Creates the .mols attribute
            Annotates the ibars using a ibar 'cluster' table
            # TODO : improve the ibar annotation functionality
        '''

        assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"

        if sample_id is None:
            sample_id = self.sample_id #type: ignore

        if do_sampling is not None:
            self.print(f'Downsampling to {do_sampling} reads', style='bold #ff0000')
            do_sampling = int(do_sampling) #type: ignore

        mols = self.load_guide_mols(
                min_reads=min_reads, 
                clone=clone, 
                suffix=suffix, 
                do_sampling=do_sampling, 
                pattern=pattern)

        # annotate .mols
        if self.ibar_ann is None and ibar_ann is None: 
            ibar_ann = self.cluster_ibars(mols=mols)

        self.ibar_ann = ibar_ann
        mols = self.annotate_ibars(mols=mols, ibar_ann=self.ibar_ann)

        # store .mols
        setattr(self, self.attr_d[suffix], mols)

        #if self.is_chimera:
        #    self.score_clone(min_cell_size=min_cell_size)

    

    @ogtk.utils.log.call
    def cluster_ibars(
            self,
            mols,
            min_cov_log10=2,
            plot = False)->pl.DataFrame:
        ''' Empirical annotation of ibars.
            
            Invokes ogtk.shlter.cluster_ibars
        '''
        clustered_ibars = shltr.cluster.cluster_ibars(
                df=mols, 
                plot=plot, 
                min_cov_log10=min_cov_log10, 
                umap=False, 
                )

        logger.debug(clustered_ibars['cluster'].value_counts())
        return clustered_ibars

    @ogtk.utils.log.call
    def annotate_ibars(
            self, 
            mols:pl.DataFrame,
            ibar_ann:pl.DataFrame,
            )->pl.DataFrame:
        '''
        applies a `.join` to a specified cluster-annotated ibar table.
        TODO: There is an inconsistency that requires unique in the end
        otherwise some columns are duplicated

        '''
        dfc = None

        if isinstance(ibar_ann, str):
            dfc = (pl.scan_parquet(ibar_ann)
                   .drop('clone')
                   .collect()
                )

        if isinstance(ibar_ann, pl.DataFrame):
            dfc = ibar_ann.drop('clone').unique()

        dfc = (mols
               .drop(['cluster'])
               .join(
                    dfc, 
                    left_on=['raw_ibar'], 
                    right_on=['raw_ibar'], 
                    how='left')
              )
        #TODO why unique?
        # the hardwired pre-computed ibar annotation table (now removed) had
        # some ibar collisions across different samples
        return(dfc.unique())

    def init_adata_std(self, kind, attr='raw_adata'):
        ''' Populates `attr` (def: .raw_adata) from a standard `kind`
        '''
        adata = ad.read_h5ad(self.return_path(kind))
        if kind == 'dedoublets_scs_ad_path':
            adata.obs['excluded_cell']  = adata.obs.doublet_prediction != 'singlet'
            adata.obs['sample_id'] = self.sample_id

            print(adata.obs.excluded_cell.mean())
            print(adata.obs.doublet_prediction.value_counts())
            print("excluding doublets\n")

            adata = adata[adata.obs.doublet_prediction == 'singlet'].copy()

        setattr(self, attr, adata)
        adata = None

    @ogtk.utils.log.call
    def init_adata(self,
        force:bool=True,
        h5_path: str|None= None,
        cellbender: bool=False,
        *args,
        **kwargs):
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
            # TODO check if there is a previous instance of the anndatas if 
            if os.path.exists(self.path_raw_h5ad) and not force:
                self.raw_adata = ad.read_h5ad(self.path_raw_h5ad)
                return None
            self.raw_adata = self.load_10_mtx()
            return None



    @ogtk.utils.log.call
    def run_scvi_solo(self,
        input_cells_path=None,
        min_counts=800,
        min_cells_gene=3,
        force=False,
        wait_time=10,
        timeout=1000,
        use_console=False,
        remove_doublets=True):

        self.is_supported('scvi_solo')

        sge_conf = self.sge_conf

        self.job_scvi_solo = ut.sge.SGE_JOB(
             job_template_path = self.path_scvi_solo_sge, 
             id = self.sample_id,
             wd = self.wd_sge, 
             sge_conf=sge_conf,
             console=self.console)
        job = self.job_scvi_solo

        input_cells_path = input_cells_path if input_cells_path is not None else f'{self.path_cr_outs}/raw_feature_bc_matrix.h5'


        dedoublets_scs_ad_path = self.return_path("dedoublets_scs_ad_path") 

        if os.path.exists(dedoublets_scs_ad_path):
            if force:
                print(f'Removed {job.wd}')
                subprocess.run(["rm", "-fr", job.wd])
            else:
                self.print(f"previous computation found. loading {dedoublets_scs_ad_path}")
                self.print(f"Change force=True to re-run")
                self.raw_adata = sc.read_h5ad(dedoublets_scs_ad_path)
                return


        done_file = f'{job.wd}/adata_scvi_solo.h5ad'

        job.fill_job_template({
            'SCRIPT': 'run_scvi_solo.py', 
            'INPUT_H5': input_cells_path,
            'MIN_COUNTS': str(min_counts), 
            'MIN_CELLS': str(min_cells_gene), 
        })
        
        files = [
                 input_cells_path,
                 self.path_scvi_solo_sh
                 ]

        if os.path.exists(done_file):
            os.remove(done_file)

        job.set_use_console(use_console)
        job.submit_job(files=files)
                
        self.print(f"Processing {input_cells_path}")
        start_time = time.time()
        while not os.path.exists(done_file):
            if timeout is not None and time.time() - start_time > timeout:
                print("Timeout reached, file not found.")
                break

            # TODO improve this loop. The while never find the done_file
            print(f"Waiting for job {job.id} to finish")
            job.qstat(times=1, sleep=wait_time)


        cmd = ["rsync", done_file, dedoublets_scs_ad_path]
        subprocess.run(cmd, check=True)
        
        self.print(f"Saved to {dedoublets_scs_ad_path} without removal!")

        adata = sc.read_h5ad(dedoublets_scs_ad_path)
        adata.obs['sample_id'] = self.sample_id
        adata.obs['excluded_cell']  = adata.obs.doublet_prediction != 'singlet'
        adata.obs['sample_id'] = self.sample_id

        self.raw_adata = adata[adata.obs.doublet_prediction == 'singlet'].copy()




    @ogtk.utils.log.call
    @wraps(db.run_cranger)
    def cellranger(self, *args, **kwargs):
        ''' cell ranger wrapper
        '''
        db.run_cranger(self, *args, **kwargs)
    
    @ogtk.utils.log.call
    @wraps(shltr.sc.compute_clonal_composition)
    def score_clone(self, min_cell_size=0, clone_dict:dict|None=None, *args, **kwargs):
        ''' Compute clonal composition 

            1. Invoques ``shltr.sc.compute_clonal_composition`` applying the
            ``min_cell_size`` argument which determines the size of a cell in
            total umi counts from the targeted library and not the gene
            expression.

            2. Merges the clonal composition data frame into the single-cell anndata object. 

            3. Populates the xp.cc (clonal composition) attribute.

        '''

        fn = shltr.sc.compute_clonal_composition

        cc = fn(self.mols.filter(pl.col('umis_cell')>min_cell_size), *args, **kwargs)
        #cc = cc.with_columns(pl.col('total_counts').fill_null(0))
        self.cc = cc.clone()

        # TODO clean the mess regarding '-1'
        #cc = cc.with_columns((pl.col('cbc') + '-1')).rename({'cbc':'index'}).to_pandas().set_index('index')
        cc = cc.rename({'cbc':'index'}).to_pandas().set_index('index')

        # annotate cells taking care of not repeating columns of ibar clusters
        self.scs.obs = (
                self.scs.obs[[i for i in self.scs.obs.columns if i not in cc.columns]]
                    .merge(
                    right=cc,
                    left_index=True,
                    right_index=True,
                    how='left')
        )

    @ogtk.utils.log.call
    def init_object(self,
        sample_name: str| None= None,
        uiport=7777,
        skip_mols=False,
        min_reads_per_mol: int=1,
        kind: str| None = None,
        forced_pattern:str|None=None,
)->None:
        ''' Main function to load full experiment object
        '''
        if sample_name is None:
            sample_name = self.sample_id #type: ignore

        self.init_workdir()

        if 'pp' in vars(self):
            db.run_cranger(self, localcores=75, uiport=uiport)

        if self.raw_adata is None and self.scs is None:
            self.print("No pre-populated cells found. Looking on disk")
            #self.init_adata()
            self.load_latest_ad(forced_pattern=forced_pattern, kind=kind)

        if self.tabulate is not None:  #type: ignore
            db.tabulate_xp(self, modality='singl-cell', cbc_len=16, umi_len=10)

        if not skip_mols:
            self.init_mols(clone=self.clone, min_reads=min_reads_per_mol) #type: ignore

    def load_latest_ad(self, forced_pattern:str|None=None, kind:str|None=None):
        ''' Loads latest version of mcs and scs adatas
        #TODO add option to fetch a specific one
        ''' 
        # check for pre-cleaned cells
        mcs_fad_path = None
        scs_fad_path = None

        if kind is not None:
            self.init_adata_std(kind)
            return

        if os.path.exists(self.path_cleansc_h5ad):
            print(f"found {self.path_cleansc_h5ad}")
            scs_fad_path = self.path_cleansc_h5ad
            self.scs = ad.read_h5ad(scs_fad_path)

        # check for final h5ad, else get latest, else get raw
        # single-cells
        if os.path.exists(self.path_cleansc_h5ad):
            scs_fad_path = self.path_cleansc_h5ad
        else:
            self.init_adata()

        if forced_pattern is None:
            # check for final h5ad, else get latest, else get raw
            if os.path.exists(self.path_final_h5ad):
                mcs_fad_path = self.path_final_h5ad
            else:
                iteration_paths = sorted(ut.sfind(self.wd_mc, "*iteration*.h5ad"))
                if len(iteration_paths)==0:
                    mcs_fad_path = self.path_raw_h5ad
                    self.print(f'No final adatas found for metacells', 'bold #ff0000', force = True)
                    return 1
                else:
                    # get latest iteration
                    mcs_fad_path = iteration_paths[-1]
        else:
            forced_paths = ut.sfind(self.wd_mc, f"*{forced_pattern}*.h5ad")
            if len(forced_paths) ==1:
                mcs_fad_path = forced_paths[0]
            elif len(forced_paths) == 0:
                # no match
                raise ValueError("No files match the provided {forced_pattern=}")
                return 1
            elif len(forced_paths) >1:
                # error there are many options
                raise ValueError("Many files matched the {forced_pattern=}")
                return 1
        
        self.print(f'loading final adatas:\n{mcs_fad_path}\n{scs_fad_path}', 'bold white')
        if mcs_fad_path is not None:
            self.mcs = ad.read_h5ad(mcs_fad_path)

        if scs_fad_path is not None:
            self.scs = ad.read_h5ad(scs_fad_path)
        return 0

        #else:
        #    self.print(f'No final adatas found, compute them manually using .do_mc(explore=True), select best parameters and re-run .do_mc(explore=False)', 'bold #ff0000', force = True)
        #    return 1
        #    raise ValueError 

    def save_final_ad(self, sample_name):
        '''
        '''
        mcs_fad_path = self.return_path('mcs_ad_path')
        scs_fad_path = self.return_path('scs_ad_path')

        self.print(f'saving final adatas:\n{mcs_fad_path}\n{scs_fad_path}', 'bold cyan')

        self.mcs.write_h5ad(mcs_fad_path)
        self.scs.write_h5ad(scs_fad_path)

    def return_path(self, kind:str, sample_name:str|None=None, suffix=None):
        """ returns a standardized path for a given `kind` of file.
            kinds are local variables that define the dataset to return the path to (needs improvement)
            e.g. mols, mcs_ad_path, mcs_fad_path, raw_scs_ad_path, clean_scs_ad_path, dedoublets_scs_ad_path
            scs_ad_path, scs_fad_path, otl_ad_path

        f"{sample_name}.iteration-1.clean"
        """
        if sample_name is None:
            sample_name = self.sample_id

        kinds="mols,mcs_ad_path,mcs_fad_path,raw_scs_ad_path,clean_scs_ad_path,scs_ad_path,scs_fad_path,otl_ad_path".split(',')


        #_f for final
        mcs_ad_path = self.path_mcs_h5ad 
        iteration_mcs_ad_path = mcs_ad_path.replace('mcells', 'iteration-ITER_mcells')
        mcs_fad_path = mcs_ad_path.replace('mcells', 'mcells_f')

        raw_scs_ad_path = mcs_ad_path.replace('mcells', 'raw_scells')
        clean_scs_ad_path = mcs_ad_path.replace('mcells', 'clean_scells')
        dedoublets_scs_ad_path = mcs_ad_path.replace('mcells', 'dedoublet_scells')

        iteration_scs_ad_path = mcs_ad_path.replace('mcells', 'iteration-ITER_scells')
        scs_ad_path = mcs_ad_path.replace('mcells', 'scells')
        scs_fad_path = mcs_fad_path.replace('mcells', 'scells')
        
        otl_ad_path = mcs_ad_path.replace('mcells', 'outliers')
        
        if suffix is not None:
            assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"
            mols = self.return_cache_path('mols', suffix=suffix)
            mols = f'{mols}/{sample_name}_r2mols.parquet'

        valid_kinds = [i for i in locals() if '_path' in i or i == 'mols']
        assert kind in valid_kinds, f"{kind} is invalid. Please provide a `kind` of the following {valid_kinds}"
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
        
    def mc_clean_plots(self,
              raw_adata: ad.AnnData, 
              properly_sampled_min_cell_total: int|float, 
              properly_sampled_max_cell_total: int|float,
              properly_sampled_max_excluded_genes_fraction: int|float,
              sample_name: str | None,
              adata_workdir: None | str = None,
              )-> None:
        """
        """

        if adata_workdir is None:
            adata_workdir=self.wd_mc


        if sample_name is None:
            sample_name =self.wd_mc

        total_umis_of_cells = mc.ut.get_o_numpy(raw_adata, name='__x__', sum=True)
        # fig total umis per cell histogram
        fg = sns.displot(total_umis_of_cells, bins=800, aspect=3, element='step')
        fg.ax.set(xlabel='UMIs', ylabel='Density', yticks=[])
        fg.ax.axvline(x=properly_sampled_min_cell_total, color='darkgreen')
        fg.ax.axvline(x=properly_sampled_max_cell_total, color='crimson')
        fg.ax.set_xlim((100, 1e5))
        fg.ax.set_title(f'{sample_name}')
        fg.savefig(f'{adata_workdir}/{sample_name}_umi_dplot.png')

        # fig total umis per cell histogram (log)
        fg = sns.displot(total_umis_of_cells, bins=800, aspect=3, element='step')
        fg.ax.set(xlabel='UMIs', ylabel='Density', yticks=[])
        fg.ax.axvline(x=properly_sampled_min_cell_total, color='darkgreen')
        fg.ax.axvline(x=properly_sampled_max_cell_total, color='crimson')
        fg.ax.set_xlim((100, 1e5))
        fg.ax.set_xscale('log')
        fg.ax.set_title(f'{sample_name} x-log')
        fg.savefig(f'{adata_workdir}/{sample_name}_umi_dlogplot.png')


        # fig total umis of excluded genes
        too_small_cells_count = sum(total_umis_of_cells < properly_sampled_min_cell_total)
        too_large_cells_count = sum(total_umis_of_cells > properly_sampled_max_cell_total)
        
        too_small_cells_percent = 100.0 * too_small_cells_count / len(total_umis_of_cells)
        too_large_cells_percent = 100.0 * too_large_cells_count / len(total_umis_of_cells)
        
        self.print(f"Will exclude {too_small_cells_count} small ({too_small_cells_percent:.2f}%)\
                cells with less than {properly_sampled_min_cell_total} UMIs")

        self.print(f"Will exclude {too_large_cells_count} large ({too_large_cells_percent:.2f}%)\
                cells with less than {properly_sampled_max_cell_total} UMIs")

        excluded_genes_data = mc.tl.filter_data(raw_adata, var_masks=[f'&excluded_gene'])

        if excluded_genes_data is not None:
            excluded_genes_data = excluded_genes_data[0]
            excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
            excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells

            too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)
            too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)
            
            self.print(f"Will exclude {too_excluded_cells_count} excluded (e.g. mito) ({too_excluded_cells_percent:.2f}%)\
                cells with less than {properly_sampled_max_excluded_genes_fraction * 100.0} UMIs")

            fg = sns.displot(excluded_fraction_of_umis_of_cells + 1e-5,
                    bins=200,
                    aspect=3,
                    element="step",
                    color='orange',
                    log_scale=(10,None))

            fg.ax.set(xlabel="Fraction of excluded gene UMIs", ylabel='Density', yticks=[])
            fg.ax.axvline(x=properly_sampled_max_excluded_genes_fraction, color='crimson')
            fg.ax.set_title(f'{sample_name}')

            print(f'{adata_workdir}/{sample_name}_fr_excluded.png')
            fg.savefig(f'{adata_workdir}/{sample_name}_fr_excluded.png')
        else:
            excluded_fraction_of_umis_of_cells = 0
        plt.show()

    def ccc(self):
        ''' this function populate the .ccc attribute with a table of clonal composition with different scores based on the ratio of integrations and molecules
        '''
        ccc = (
            self.mols
            .rename({'cluster':'ibc'})
            .filter(pl.col('ibc').is_not_null())
            .with_columns(pl.n_unique('raw_ibar').over('cbc', 'ibc', 'doublet_prediction').alias('integrations'))
            .with_columns(pl.count().over('cbc', 'ibc', 'doublet_prediction').alias('molecules'))
            .with_columns(pl.count().over('cbc').alias('mols_per_cell'))
            .select('doublet_prediction', 'cbc', 'integrations', 'ibc', 'molecules', 'mols_per_cell')
            .unique()
            .pivot(index=['doublet_prediction','cbc', 'mols_per_cell'], columns='ibc', values=['integrations', 'molecules'])
            .fill_null(0)
            #.sort("mols_per_cell", 'cbc' )
            #.filter(pl.col("mols_per_cell")<1000)
            #.filter(pl.col("mols_per_cell")>100)
            .filter(pl.col('doublet_prediction')=="singlet")
            .with_columns(
                        score_raw=pl.col('molecules_ibc_2')/pl.col('molecules_ibc_1'),
                        score_ints=pl.col('integrations_ibc_2')/pl.col('integrations_ibc_1'),
                        score_norm=(pl.col('molecules_ibc_2')/pl.col('integrations_ibc_2'))/(pl.col('molecules_ibc_1')/pl.col('integrations_ibc_1')),
            )
            .with_columns(
                        pl.col('score_raw').log10(),
                        pl.col('score_ints').log10(),
                        pl.col('score_norm').log10(),
                        )
            .with_columns(
                        clone=pl.when(pl.col('score_norm')<0).then(pl.lit("left")).otherwise(pl.lit("right"))
            )
            .with_columns(
                        clone=pl.when(
                            (pl.col('score_norm')<0.25) &(pl.col('clone')=='right'))
                        .then(pl.lit("mid")).otherwise(pl.col("clone"))
            )
        )
        self.cc = ccc

    def mc_clean(
        self,
        sample_name:str|None=None, 
        raw_adata:ad.AnnData|None=None, 
        excluded_gene_patterns:Sequence|None=None, 
        excluded_gene_names:Sequence|None=None, 
        suspect_gene_names = Sequence | None,
        suspect_gene_patterns = Sequence | None,
        manual_ban: Sequence | None=[],
        lateral_modules: Sequence | None=None,
        return_adatas=True, 
        log_debug=False, 
        explore=True, 
        cpus = {'full':56, 'moderate':8},
        mc_cpus_key='moderate',
        var_cpus_key='moderate',
        dpi=90,
        properly_sampled_max_excluded_genes_fraction=0.03,
        properly_sampled_min_cell_total=500,
        properly_sampled_max_cell_total=20000,
        force:bool=False,
        random_seed=123456,
        grid:bool=True,
        max_parallel_piles=None,
                    ):
        '''
        This function is the main entry point for the metacells pipeline

        '''
        import metacells.utilities.typing as utt

        if raw_adata is None:
            raw_adata = self.raw_adata

        if sample_name is None:
            sample_name = self.sample_id

        if suspect_gene_names is None:
            suspect_gene_names = self.suspect_gene_names

        if suspect_gene_patterns is None:
            suspect_gene_patterns = self.suspect_gene_patterns

        print(f'mc2 v{mc.__version__}')

        if log_debug:
            import logging
            mc.ut.setup_logger(level=logging.DEBUG)
            print(np.show_config())

        plt.rcParams['figure.dpi'] = dpi

        # sanitize andata
        utt.sum_duplicates(raw_adata.X)
        utt.sort_indices(raw_adata.X)

        mc.ut.set_name(raw_adata, sample_name)

        if lateral_modules is None:
            lateral_modules = []

        if 'excluded_gene' in raw_adata.var.columns and not force:
            print('This anndata does not look raw since some genes have been marked already as excluded, use force=True to start from scratch')
            return raw_adata
            
        if excluded_gene_patterns is None:
            excluded_gene_patterns = self.excluded_gene_patterns

        if excluded_gene_names is None:
            excluded_gene_names = self.excluded_gene_names

        # .var gets masks: bursty_lonely_gene, properly_sampled_gene, excluded_gene
        mc.pl.exclude_genes(
            raw_adata, 
            excluded_gene_patterns=excluded_gene_patterns,
            excluded_gene_names=excluded_gene_names,
            random_seed=random_seed)

        # .obs gets masks properly_sampled_cell, excluded_cell 
        mc.pl.exclude_cells(
            raw_adata, 
            properly_sampled_min_cell_total=properly_sampled_min_cell_total, 
            properly_sampled_max_cell_total=properly_sampled_max_cell_total, 
            properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)
        
        if 'mc_clean' in raw_adata.uns and not force:
            print('This anndata has the "mc_clean" token, use force=True to start from scratch')
            return raw_adata

        ##### 
        self.mc_clean_plots(raw_adata, 
                            properly_sampled_min_cell_total, 
                            properly_sampled_max_cell_total, 
                            properly_sampled_max_excluded_genes_fraction,
                            sample_name)
        #####
        raw_adata.uns['mc_clean'] = True
        clean = mc.pl.extract_clean_data(raw_adata, name=f"{sample_name}.iteration-1.clean")
        ###

        # save anndatas
        raw_adata.write_h5ad(self.path_raw_h5ad)
        clean.write_h5ad(self.path_cleansc_h5ad)

        self.print(f"Saved {self.path_raw_h5ad}")
        self.print(f"Saved {self.path_cleansc_h5ad}")

        self.scs = clean
        raw_adata = None  
        self.print('Clean cells have populated the .scs attribute')

    @wraps(mc.tl.convey_obs_fractions_to_group)
    def mc_convey_cell_annotations_to_metacells(self, property_name="sample_id") -> None:
        mc.tl.convey_obs_fractions_to_group(adata=self.scs, gdata=self.mcs, property_name=property_name)

    def mc_compute_next_iteration(self, iteration: int, cores=10, random_seed=123456) -> None:
        
        max_parallel_piles = mc.pl.guess_max_parallel_piles(self.scs)

        mc.pl.set_max_parallel_piles(max_parallel_piles)
        mc.ut.set_processors_count(cores)

        self.print("# DIVIDE AND CONQUER...")
        global metacells
        self.mcs = None # So can be gc-ed

        mc.pl.divide_and_conquer_pipeline(self.scs, random_seed=random_seed)

        self.print("# COLLECT METACELLS...")
        self.mcs = mc.pl.collect_metacells(
            self.scs, name=f"{self.sample_id}.iteration-{iteration}.metacells", 
            random_seed=random_seed
        )
        self.print(f"Iteration {iteration}: {self.mcs.n_obs} metacells, {self.mcs.n_vars} genes")

        self.print("# CONVEY CELL ANNOTATIONS...")
        self.mc_convey_cell_annotations_to_metacells()

    def mc_finalize_next_iteration(self, iteration: int, *, with_types: bool, cores=10, random_seed=123456) -> None:
        max_parallel_piles = mc.pl.guess_max_parallel_piles(self.scs)
        mc.pl.set_max_parallel_piles(max_parallel_piles)
        mc.ut.set_processors_count(cores)

        self.print("# COMPUTE FOR MCVIEW...")
        mc.pl.compute_for_mcview(adata=self.scs, gdata=self.mcs, random_seed=random_seed)

        self.print("# PLOT UMAP...")
        if with_types:
            type_annotation = f"type.iteration-{iteration}.auto"
        else:
            type_annotation = None

        self.mc_plot_umap(type_annotation=type_annotation)

        self.print("# SAVE CELLS...")
        scs_path =self.return_path('iteration_scs_ad_path').replace('ITER', str(iteration))
        self.scs.write_h5ad(scs_path)

        self.print("# SAVE METACELLS...")
        mcs_path =self.return_path('iteration_mcs_ad_path').replace('ITER', str(iteration))
        self.mcs.write_h5ad(mcs_path)
        
        self.print(mcs_path)
        self.print(scs_path)
        #print("# IMPORT TO MCVIEW...")
        #os.system(
        #    f"Rscript ../scripts/import_dataset.r hca_bm iterative/iteration-{iteration} "
        #    f"'HCABM IT|{iteration}'"
        #)

    def mc_next_iteration_without_types(self, iteration: int, cores = 10) -> None:
        '''
        '''
        self.mc_compute_next_iteration(iteration, cores=cores)
        self.mc_finalize_next_iteration(iteration, with_types=False, cores=cores)

    @wraps(_color_hist)
    def color_hist(self, *args, **kwargs):
        _color_hist(*args, **kwargs)

    @wraps(ut.sc.mc_plot_umap)
    def mc_plot_umap(self, type_annotation: str|None = None) -> None:
        ut.sc.mc_plot_umap(self.mcs, type_annotation)

    @wraps(ut.sc.mc_relate_to_lateral_genes)
    def mc_relate_to_lateral_genes(self, force:bool=False, *args, **kwargs):
        ut.sc.mc_relate_to_lateral_genes(self.scs, force, *args, **kwargs)

    @wraps(ut.sc.mc_update_lateral_genes)
    def mc_update_lateral_genes(self, *args, **kwargs):
        ut.sc.mc_update_lateral_genes(cells=self.scs, *args, **kwargs)

    @wraps(ut.sc.mc_compute_lateral_module_similarity)
    def mc_compute_lateral_module_similarity(self, *args, **kwargs):
        ut.sc.mc_compute_lateral_module_similarity(cells=self.scs, *args, **kwargs) 

    @wraps(ut.sc.mc_plot_associated_lmodules)
    def mc_plot_associated_lmodules(self, *args, **kwargs):
        if 'similarity_of_modules' not in self.scs.uns:
            self.mc_compute_lateral_module_similarity()
        ut.sc.mc_plot_associated_lmodules(adata=self.scs, fig_dir=self.wd_figs)
        
    @wraps(ut.sc.mc_update_lateral_flags)
    def mc_update_lateral_flags(self, attr, *args, **kwargs):
        adata = getattr(self, attr) 
        ut.sc.mc_update_lateral_flags(adata=adata, *args, **kwargs)

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
                    sample_name=sample_name,
                    adata=adata if adata is not None else self.raw_adata, 
                    adata_workdir=self.wd_mc,
                    explore=explore,
                    force=force,
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
       
        self.scs.write_h5ad(self.path_cleansc_h5ad) #type: ignore
        self.mcs.write_h5ad(self.path_iterationx_h5ad.replace('XXX', '1')) #type: ignore

        #self.save_final_ad(sample_name)

    @wraps(shltr.sc.allele_calling)
    def init_alleles(self, 
                     suffix='shrna',
                     expr:None|pl.Expr=None,
                     ann_source:str|None|pl.DataFrame=None,
                     *args,
                     **kwargs):

        #mols = getattr(self, self.attr_d[suffix] )
        assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"


        alleles = shltr.sc.allele_calling(self.mols, *args, **kwargs)

        if self.cc is None:
            self.ccc()
        
        alleles = (
                alleles
                .join(self.cc, left_on='cbc', right_on='cbc', how='left')
                .join(self.mols.select('raw_ibar', 'cluster').unique(), left_on='raw_ibar', right_on='raw_ibar')
                )

        #for ax in g.axes.flat:
        #    ax.set_xlim((4, 8))

        self.alleles = alleles
        return 
        # merge with scs.obs
        # TODO a bug lurks here
        strip_pattern = '(.+?)-(.)' if self.batch_map is not None else '(.+)'
        strip_pattern = '(.+)'


        # the majority of the raw cell barcodes do not belong to any real cell so we constrain the merge to the scs
        # how='outer'
        # adobs_pd_to_df Converts an adata.obs pd data frame into pl data frame...
        alleles = alleles.join(
                ut.sc.adobs_pd_to_df(self.scs, strip_pattern), 
                left_on='cbc', 
                right_on='cbc', 
                how='inner') 

        # since polars lacks a 'how=right' merging, it is needed to drop nulls
        # in a field that is meaningful for the alleles, such as 'umis_cell'
        #alleles = alleles.filter(pl.col('umis_cell').is_not_null())

        setattr(self, suffix[0]+'alleles', alleles)
        #if self.is_chimera:
        #    # TODO clean this block, e.g. log2 or 10?
        #    # the chimera flag should be complemented by the definition of which clones comprise it
        #    # and its corresponding annotation functions, e.g. fiels h1 and g10 in the following lines
        #    self.salleles = self.salleles.with_columns((pl.col('h1')/pl.col('g10')).log(base=10).alias('cs'))
        #    self.salleles = self.salleles.with_columns(pl.col(pl.Categorical).cast(pl.Utf8))

        #    self.scs.obs['cs'] = np.log10(self.scs.obs['g10']/self.scs.obs['h1'])

        if ann_source is not None:
            self.annotate_alleles(ann_source)

    @wraps(_plot_alleles)
    def plot_alleles(self, *args, **kwargs)->sns.FacetGrid:
        logger.step("plotting allele histograms")
        return _plot_alleles(*args, **kwargs)

    @wraps(_cas_update_mc_ann)
    def update_scs_from_mc_ann(self, mc_ann=None):
        ''' Imports MCView output into .scs
        ''' 
        if mc_ann is not None:
            self.load_mc_ann(mc_ann)
        self.scs.obs = _pd_annotate_mc(df=self.scs.obs, mc_ann=self.mc_ann)


    @wraps(shltr.sc.to_matlin)
    def init_matlin(self, do_cluster = False, sort_by=['cluster', 'raw_ibar'], cells= 100, cores=4, *args, **kwargs):

        plt.rcParams['figure.dpi'] = 100
        fn = shltr.sc.to_matlin
        self.matl = fn(self.alleles, cells=cells, sort_by=sort_by,  *args, **kwargs)

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
            logger.critical('Loading pre-computed valid ibars')
            valid_ibars = pl.read_csv('/local/users/polivar/src/artnilet/conf/valid_ibars.csv')['valid_ibar'].to_list()      
            return(valid_ibars)

    @wraps(_plot_cc)
    def plot_cc(self, *args, **kwargs):
        _plot_cc(self, *args, **kwargs)

    def ingest_xps(self,
        xps: Iterable,
        suffix='shrna',
        force=False,
        skip_mols=False,
        filter_valid_cells=False,
        ibar_ann: str|None =None):
        ''' Merge compatible experiments and integrate them into an invidiual one.
            a ``batch_map`` keeps track of an appended identifier to single-cell barcodes

            ibar_ann: path to ibar annotations. None (default) resources to the default path specified on `.annotate_ibar()`
        '''
        assert suffix in ['zshrna', 'shrna'], "Invalid suffix, please use 'zshrna' or 'shrna'"

        if self.mols is not None and self.raw_adata is not None:
            self.print('This experiment seems to be have ingested others. No need to ingest unless force=True')
            return None

        adatas = [i.scs for i in xps]
        
        if not any([True for c in [x.obs.columns for x in adatas] if 'excluded_cell' in c]):
            raise ValueError("The field 'excluded_cell' is not present in some adatas.obs columns. Please include (the metacell scripts do it automatically).")

        adata = ad.concat(adatas, join='outer', label='batch_id', index_unique="_")
        print(adata.obs.excluded_cell.value_counts())
        adata = adata[~adata.obs['excluded_cell'], ].copy()
        adata.var = adata.var.drop(adata.var.columns[2:].to_list(), axis='columns')
        adata.obs = adata.obs[['sample_id', 'batch_id']].copy()
        
        batch_map = adata.obs.set_index('sample_id')['batch_id'].to_dict()

        self.batch_map = batch_map
        assert len(batch_map) == len(adata.obs.sample_id.unique()), "One or more batches are duplicated"

        self.scs = adata
        
        ad_path = 'raw_scs_ad_path'
        self.scs.write(self.return_path(ad_path))

        logger.io(f"saved raw ingested cells as {self.return_path(ad_path)}")
        self.batch_map = batch_map
        self.conf_keys.append('batch_map')

        if skip_mols:
            self.export_xpconf(xp_conf_keys = set(self.conf_keys))
            return None

        #assert not any([xp.mols is None for xp in xps]) , "Molecules need to be initialized manually"
        #self.mols = pl.concat([i.mols for i in xps])

        self.mols = pl.concat(
               [i.load_guide_mols(clone=i.clone, filter_valid_cells=filter_valid_cells) 
                for i in xps]
        )
       
        logger.info(f"{self.mols.shape=}")

        self.mols = (
                self.mols
                .with_columns(
                    pl.col('sample_id')
                    .map_dict(batch_map)
                    .alias('batch_id')
                )        
                .with_columns(pl.col('cbc')+"_"+pl.col('batch_id').alias('cbc'))
        )

        logger.info(f"{self.mols.shape=}")

        # annotate .mols
        if self.ibar_ann is None and ibar_ann is None: 
            ibar_ann = self.cluster_ibars(mols=self.mols)

        self.ibar_ann = ibar_ann
        self.mols = self.annotate_ibars(mols=self.mols, ibar_ann=self.ibar_ann)

        logger.info(f"{self.mols.shape=}")

        # save parquet file
        logger.io(f"saving molecules to {self.return_path('mols', suffix=suffix)}")
        self.mols.write_parquet(self.return_path('mols', suffix=suffix))

        self.export_xpconf(xp_conf_keys = set(self.conf_keys))

    def return_cache_path(self, key, suffix=None):
        '''
        '''
        if key =='mols':
            assert suffix is not None, "A suffix for the type of molecule is needed, e.g., 'shrna', 'zhrna'"
            return f'{self.wd_xp}/{suffix}/'


    @wraps(_cassit)
    def cassit(self, *arg, **kwargs):
        _cassit(self, *arg, **kwargs)

    @wraps(_cas_return_colors)
    def return_allele_colors(self):
        _cas_return_colors(self.tree.clone_allele_table)

    @wraps(_cas_update_mc_ann)
    def tree_update_mc_ann(self, mc_ann: pl.DataFrame):
        _cas_update_mc_ann(self, mc_ann)

    def load_mc_ann(self, source: str|pl.DataFrame):
        ''' Loads MCView-exported metacell annotation and creates a color map based on the MCView annotation
        Populates the experiment's:
         - .mc_ann
         - .mc_color_dict
         - .mc_colormap 

         It additionally registers .mc_colormap into matplotlib's cmap registry
        '''
        match source:
            case pl.DataFrame():
                self.mc_ann = source
            case str():
                self.mc_ann = pl.read_csv(source, null_values=['NA'])

        self.mc_color_dict = (self.mc_ann[:, ['cell_type', 'color']]
                              .unique()
                              .sort('cell_type')
                              .to_pandas()
                              .set_index('color')
                              .to_dict()['cell_type']
                              )

        if self.mc_colormap is not None:
            plt.colormaps.unregister('mc_colormap')

        self.mc_colormap = ListedColormap(self.mc_color_dict.keys())
        plt.colormaps.register(name='mc_colormap', cmap=self.mc_colormap, force=True)

    def annotate_alleles(self, source:str|None=None):
        ''' 
        '''
        if source is not None:
            self.load_mc_ann(source)

        if source is None and self.mc_ann is None:
            raise ValueError("Load an MCView annotation first using `.load_mc_ann`")

        self.salleles = self.salleles.with_columns(pl.col('metacell_name').cast(pl.Utf8()))
        self.salleles = self.salleles.join(self.mc_ann, left_on='metacell_name', right_on='metacell', how='left')

    def lineage_analysis(
            self,
            exp, 
            clone = 'g10',
            cells=1000, 
            clc_dict = {'h1':[-1, 1, 2], 'g10':[0, 3]}):
        import cassiopeia as cas

        exp.is_chimera = True
        exp.init_alleles(by='umis_allele', descending=True)
        expr = (
            (pl.col('umis_allele')>2) & \
            (pl.col('metacell_name')!='Outliers')  &\
            #(pl.col('raw_ibar').is_in(top_cut_ibars['raw_ibar'])) &\
            (pl.col('total_umis')>=10) 
            )
        exp.alleles = exp.salleles
        alleles = exp.alleles.filter(expr)
        
        exp.init_matlin(cells=cells, expr=expr)
        
        exp.x(alleles=alleles, clusters=clc_dict[clone], allele_rep_thresh=1)

        # plot tree
        cas.pl.plot_matplotlib(exp.tree, 
                          orient='right', 
                           allele_table=exp.tree.clone_allele_table, 
                           indel_colors=pp._cas_return_colors(exp.tree.clone_allele_table),
                           figsize=(1*7*1, 1*7*1),
                           colorstrip_spacing = 0, 
                       colorstrip_width = 4, 
                          )
        plt.show()
        return exp


    @wraps(_check_required_attributes)
    def is_supported(self, case):
        return _check_required_attributes(self, case)

@wraps(db.print_template)
def print_template(*args, **kwargs):
    db.print_template(*args, **kwargs)

