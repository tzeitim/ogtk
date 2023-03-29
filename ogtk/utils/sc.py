from anndata import AnnData
import os
from functools import wraps
import matplotlib.pyplot as plt
import metacells as mc
import metacells.utilities.typing as utt
import numpy as np
import os
import pandas as pd
import polars as pl
import regex
import rich
import scipy.sparse as sp
import scanpy as sc
import seaborn as sns
import time
from typing import Sequence

def len_intersect(adata, rs):
    good_cbs = set([i.replace('-1', '') for i in adata.obs.index])
    lib_cbs = set(rs.umis.keys())
    iset = len(good_cbs.intersection(lib_cbs))
    return(np.array([iset, iset/adata.obs.shape[0],iset/len(rs.umis.keys()) ]))

def print_stats(adata, rs):
    iset = len_intersect(adata, rs)
    print(f"A total of {iset[0]} cbs were found on both sets {iset[1]:0.2f} {iset[2]:0.2f} ")

def metacellize(
    set_name:str, 
    adata:AnnData, 
    adata_workdir:str,
    excluded_gene_patterns:Sequence=[], 
    excluded_gene_names: Sequence|None=None, 
    target_metacell_size=100,
    suspect_gene_names=['shRNA', 'rtTA', 'Fun','Cas9', 'Hist1h2ap'], 
    suspect_gene_patterns = ['mt-.*', 'Dnaj.*'],
    usual_suspects=['Pcna', 'Pclaf', 'Jun', 'Top2a', 'Txn', 'Hsp90ab1', 'Fos', 'Dnaj', 'Ube2c'],
    force_usual_suspects=True,
    manual_ban: Sequence | None=[],
    lateral_modules: Sequence | None=None,
    return_adatas=True, 
    log_debug=False, 
    explore=True, 
    cpus = {'full':56, 'moderate':16},
    mc_cpus_key='moderate',
    var_cpus_key='moderate',
    properly_sampled_max_excluded_genes_fraction=0.03,
    properly_sampled_min_cell_total=500,
    properly_sampled_max_cell_total=20000,
    force:bool=False,
    grid:bool=True,
                ):
    '''
    '''

    print(f'mc2 v{mc.__version__}')

    if log_debug:
        import logging
        mc.ut.setup_logger(level=logging.DEBUG)
        print(np.show_config())

    plt.rcParams['figure.dpi'] = 180

    # sanitize andata
    utt.sum_duplicates(adata.X)
    utt.sort_indices(adata.X)

    mc.ut.set_name(adata, set_name)

    if lateral_modules is None:
        lateral_modules = []

    clean_adata(adata, properly_sampled_min_cell_total=properly_sampled_min_cell_total, force=force)
    ##
    adata = qc_masks(adata, 'test', '.', force=force).copy()
    print(adata)

    pl_var = (
            analyze_related_genes(
                adata,
                adata_workdir=adata_workdir,
                suspect_gene_names=suspect_gene_names,
                usual_suspects=usual_suspects,
                suspect_gene_patterns=suspect_gene_patterns, 
                set_name=set_name, 
                lateral_modules=lateral_modules,
                manual_ban=manual_ban,
                grid=grid,
                explore = explore)
            )
    if explore:
        return(pl_var)

    mcs, scs, mdt = invoque_mc(
        adata,
        set_name=set_name,
        adata_workdir=adata_workdir,
        cpus=cpus,
        mc_cpus_key=mc_cpus_key,
        var_cpus_key=var_cpus_key,
    )
    return([mcs, scs, mdt])


def adobs_pd_to_df(adata):
    ''' Converts an adata.obs pd data frame to a pl data frame. It strips out the '-1' cell ranger suffix.
    '''
    import polars as pl
    adata.obs.index.name = None
    df = adata.obs.reset_index()
    
    return(pl.DataFrame(df)
           .rename({'index':'cbc'})
           .with_columns(pl.col('cbc').str.extract('(.+?)-(.)',1))
           )


def scanpyfi(
       adata,
       max_counts = None, 
       min_counts = None, 
       min_n_genes_by_counts = None, 
       max_n_genes_by_counts = None, 
       max_pct_counts_mt = None,
       blacklisted_genes:Sequence|None=None,
       n_pcs=40,
       s=10,
      qc = False):
    '''Run vanilla scanpy clustering 
    '''
    import scanpy as sc
    import polars as pl
    import numpy as np
    import matplotlib.pyplot as plt

    print(adata)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    if qc:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.0, multi_panel=True)

        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    
        print(pl.DataFrame(adata.obs).sort('total_counts'))
    
    #adata = adata[adata.obs.total_counts <= 60000, :].copy()
    
    if max_counts is not None:
        adata = adata[adata.obs.total_counts <= max_counts, :].copy()
        
    if min_counts is not None:
        adata = adata[adata.obs.total_counts >= min_counts, :].copy()
        
    if max_n_genes_by_counts is not None:
        adata = adata[adata.obs.n_genes_by_counts <= max_n_genes_by_counts, :].copy()
        
    if min_n_genes_by_counts is not None:
        adata = adata[adata.obs.n_genes_by_counts >= min_n_genes_by_counts, :].copy()     
        
    if max_pct_counts_mt is not None:
        adata = adata[adata.obs.pct_counts_mt <= max_pct_counts_mt, :].copy()
        
    #sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    #sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    
    if blacklisted_genes is None:
        raw_b = return_raw_gene_str()
        blacklisted_genes = raw_b.split(' ')

    black_mask = adata.var_names.isin(blacklisted_genes)

    np.sum(black_mask)
    
    sc.pp.normalize_total(adata, target_sum=1e4)

    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)
    
    adata.var.highly_variable = adata.var.highly_variable & ~black_mask

    np.sum(adata.var.highly_variable)
    
    sc.tl.pca(adata, svd_solver='arpack', n_comps=2*n_pcs)

    #sc.pl.pca(adata, color='total_counts', s=60)

    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=2*n_pcs)

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

    sc.tl.leiden(adata)

    sc.tl.paga(adata)

    sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph

    sc.tl.umap(adata, init_pos='paga')#, min_dist=0.125)#, n_components=3)


    blood= [ i for i in adata.var_names if i.startswith('Hb')]
    sox= [ i for i in adata.var_names if i.startswith('Sox')]
    hox= [ i for i in adata.var_names if i.startswith('Hox')]

    with plt.rc_context({'figure.dpi':100}):
        sc.pl.umap(adata, color=['total_counts', 'leiden'], s=s, ncols=2, )
        plt.figure()
        sc.pl.umap(adata, color=blood, s=s, ncols=2)
        plt.figure()
        sc.pl.umap(adata, color=['shRNA', 'Cas9'], s=s, ncols=2)
        plt.figure()
        
        #sc.pl.umap(adata, color=blood, s=3, ncols=2, )
        
    #sc.tl.rank_genes_groups(sstenex20c, 'leiden', method='wilcoxon')
    #sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

    #adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X

    #sc.pl.rank_genes_groups_matrixplot(adata, n_genes=10, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled')
    return(adata)

def clean_adata(
    adata,
    force:bool=False,
    lateral_modules = None,
    excluded_gene_patterns: Sequence | None= None,
    random_seed = 12345,
    properly_sampled_max_excluded_genes_fraction=0.03,
    properly_sampled_min_cell_total=500,
    properly_sampled_max_cell_total=20000,
               
               ):
    ''' 
    '''
    if 'excluded_gene' in adata.var.columns and not force:
        print('use the force')
        return adata
    if lateral_modules is None:
        lateral_modules = []
        
    if excluded_gene_patterns is None:
        excluded_gene_patterns = []
    excluded_gene_patterns.append('mt-.*')

    # .var gets masks: bursty_lonely_gene, properly_sampled_gene, excluded_gene
    mc.pl.exclude_genes(
        adata, 
        excluded_gene_patterns=excluded_gene_patterns,
        random_seed=random_seed)
    # .obs gets masks properly_sampled_cell, excluded_cell 
    mc.pl.exclude_cells(
        adata, 
        properly_sampled_min_cell_total=properly_sampled_min_cell_total, 
        properly_sampled_max_cell_total=properly_sampled_max_cell_total, 
        properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)
    
    adata.uns['properly_sampled_min_cell_total'] = properly_sampled_min_cell_total
    adata.uns['properly_sampled_max_cell_total'] = properly_sampled_max_cell_total
    adata.uns['properly_sampled_max_excluded_genes_fraction'] = properly_sampled_max_excluded_genes_fraction
    
    return adata

def qc_masks(adata, 
             set_name, 
             adata_workdir,
             force:bool=False,
            ):
    ''' wraps mc.pl.extract_clean_data(adata) together with histograms for cell, gene and gene module exclusion.
    '''
    total_umis_of_cells = mc.ut.get_o_numpy(adata, name='__x__', sum=True)

    if 'mc_clean' in adata.uns and not force:
        print('use the force')
        return adata
    
    properly_sampled_min_cell_total = adata.uns['properly_sampled_min_cell_total']
    properly_sampled_max_cell_total = adata.uns['properly_sampled_max_cell_total']
    properly_sampled_max_excluded_genes_fraction = adata.uns['properly_sampled_max_excluded_genes_fraction']
   
    # fig total umis per cell histogram (log)
    fg = sns.displot(total_umis_of_cells, bins=800, aspect=3, element='step')
    fg.ax.set(xlabel='UMIs', ylabel='Density', yticks=[])
    fg.ax.axvline(x=properly_sampled_min_cell_total, color='darkgreen')
    fg.ax.axvline(x=properly_sampled_max_cell_total, color='crimson')
    fg.ax.set_xlim((1, 1e5))
    fg.ax.set_xscale('log')
    fg.ax.set_title(f'{set_name}')
    fg.savefig(f'{adata_workdir}/{set_name}_umi_dplot.png')


    # fig total umis of excluded genes
    too_small_cells_count = sum(total_umis_of_cells < properly_sampled_min_cell_total)
    too_large_cells_count = sum(total_umis_of_cells > properly_sampled_max_cell_total)
    
    too_small_cells_percent = 100.0 * too_small_cells_count / len(total_umis_of_cells)
    too_large_cells_percent = 100.0 * too_large_cells_count / len(total_umis_of_cells)
    
    print(f"Will exclude %s (%.2f%%) cells with less than %s UMIs"
            % (too_small_cells_count,
                too_small_cells_percent,
                properly_sampled_min_cell_total))
    print(f"Will exclude %s (%.2f%%) cells with more than %s UMIs"
            % (too_large_cells_count,
                too_large_cells_percent,
                properly_sampled_max_cell_total))

    excluded_genes_data = mc.tl.filter_data(adata, var_masks=['excluded_gene'])[0]
    excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
    excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells

    too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)
    too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)
    
    print(f"Will exclude %s (%.2f%%) cells with more than %.2f%% excluded gene UMIs"
            % (too_excluded_cells_count,
                too_excluded_cells_percent,
                100.0 * properly_sampled_max_excluded_genes_fraction))
    
    fg = sns.displot(excluded_fraction_of_umis_of_cells, bins=200, aspect=3, element="step")
    fg.ax.set(xlabel="Fraction of excluded gene UMIs", ylabel='Density', yticks=[])
    fg.ax.axvline(x=properly_sampled_max_excluded_genes_fraction, color='crimson')
    fg.ax.set_title(f'{set_name}')
    fg.savefig(f'{adata_workdir}/{set_name}_fr_excluded.png')

    adata.uns['mc_clean'] = True
    clean = mc.pl.extract_clean_data(adata)
    return clean
    ###
    
    
def analyze_related_genes(
      adata, 
      adata_workdir,
      set_name,
      usual_suspects,
      suspect_gene_names,
      force_usual_suspects=True,
      explore=True,
      suspect_gene_patterns: None | Sequence = None,
      lateral_modules=[],
      force:bool=False,
      grid=False,
      dpi=900,
      columns=4,
      unit=6,
      aspect_f=0.75,
      cluster = True,
      manual_ban:Sequence=[],
     ):
    
    if 'related_genes_similarity' not in adata.varp or force:
        mc.pl.relate_genes(adata, random_seed=123456)
    else:
        print('use the force')
    
    for i in usual_suspects:
        suspect_gene_names.append(i)

    if suspect_gene_patterns is None:
        suspect_gene_patterns = []

    for i in ['Rpl', 'Rps', 'Mcm', 'Hsp', 'Hist', 'mt-']:
        suspect_gene_patterns.append(i)
    print(f'{suspect_gene_patterns=}')

    suspect_genes_mask = mc.tl.find_named_genes(
                                adata,
                                names=suspect_gene_names,
                                patterns=suspect_gene_patterns)

    suspect_gene_names = sorted(adata.var_names[suspect_genes_mask])
    
    ###
    pl_var = pl.DataFrame({"gene_name":adata.var_names.tolist()})
    pl_var = (pl_var.join(
                    pl.DataFrame(adata.var.reset_index()), 
                    left_on='gene_name', 
                    right_on='index')
                   )
    
    pl_var =  pl_var.with_columns((
                        (pl.col("gene_name").is_in(usual_suspects)) |
                        (pl.col("gene_name").str.contains('|^'.join(suspect_gene_patterns)) )
                        ).alias('suspect_gene')
                    )
    pl_var =(pl_var
             .with_columns(pl.col('suspect_gene')
                           .any().over('related_genes_module')
             .alias('suspect_module'))
            )
    
    suspect_gene_names_pl = (
            pl_var.filter(pl.col('suspect_gene')
                          )['gene_name'].sort().to_list()
            )
    suspect_gene_modules_pl = (
            pl_var.filter(
                        (pl.col('suspect_gene')) & (pl.col('related_genes_module')>=0)
            )['related_genes_module'].unique()                                        
             )

    module_of_genes = adata.var['related_genes_module']
    suspect_gene_modules = np.unique(module_of_genes[suspect_genes_mask])
    suspect_gene_modules = suspect_gene_modules[suspect_gene_modules >= 0]

    #for i in sorted(suspect_gene_modules):
    #    soc = module_of_genes[module_of_genes == i]
        #print(f'm{i}::{len(soc)}::\t{"  ".join(sorted(soc.index.to_list()))}') 
    
    suspect_gene_modules_pl = (
            pl_var
            .filter(pl.col('suspect_module'))
            .groupby(['related_genes_module'])
            .agg(pl.col('gene_name'))
            .with_columns(pl.col('gene_name').arr.join(", "))
            .sort('related_genes_module')
          )
    
    all_modules = [i for i in np.unique(module_of_genes) if int(i) >0]
    all_modules_pl = pl_var.filter(pl.col('related_genes_module')>=0)['related_genes_module']
    
    ###
    similarity_of_genes = mc.ut.get_vv_frame(adata, 'related_genes_similarity')

    row_factor = len(suspect_gene_modules_pl)
    rows =(row_factor//columns) + (row_factor % columns >0 )
 
    # plotting
    if grid:
        plt.rcParams.update(plt.rcParamsDefault)
        fig, axes = plt.subplots(rows, columns, dpi=dpi, figsize=(unit, aspect_f*unit * rows ))
        iaxes = iter(axes.flat)

    def heatmap_mod(x, cmap='RdYlBu'):
        ''' polars-based module plotting
        ''' 
        gene_module = x['related_genes_module'].unique()[0]
        modg = x['gene_name'].unique().to_list()

        similarity_of_module = similarity_of_genes.loc[modg, modg]

        if cluster:
            from scipy.cluster import hierarchy 
            Z = hierarchy.linkage(similarity_of_module, 'ward')
            zl = hierarchy.leaves_list(Z)
            similarity_of_module = similarity_of_module.iloc[zl, zl]

        labels = (x.with_columns(
                    pl.when(pl.col('gene_name').is_in(suspect_gene_names))
                    .then(pl.col('gene_name').str.replace('^', '(*)'))
                    .otherwise(pl.col('gene_name'))
                    ))['gene_name']

        similarity_of_module.index = \
        similarity_of_module.columns = labels

        lateral_txt ="(ignored)" if (gene_module in lateral_modules and not explore)  else " " 

        sns.set(font_scale=0.5)
        if grid:
            ax = next(iaxes)
            ax.tick_params(labelsize=2, width=0, pad=-2, axis='both', which='major')
            ars = sns.heatmap(similarity_of_module,
                              xticklabels=1,
                              vmin=0,
                              vmax=1,
                              cmap=cmap,
                              linewidths=0.01,
                              square=True,
                              ax = ax,
                              cbar=False)

        else:
            ax = sns.heatmap(similarity_of_module,
                              xticklabels=1,
                              vmin=0,
                              vmax=1,
                              cmap=cmap,
                              linewidths=0.01,
                              square=True,
                              cbar=False)

        title_fontsize = 2 if grid else 7.5
        ax.set_title(f'{set_name} Gene Module {gene_module} {lateral_txt}', fontsize=title_fontsize)
        if not grid:
            plt.show()
        #fig.savefig(f'{adata_workdir}/{set_name}_module_{gene_module}.png')
        sns.set(font_scale=1)

        return pl.DataFrame()

    # apply heatmap plotting function
    ret = (pl_var
           .filter(
                (pl.col('suspect_module')) &
                (pl.col('related_genes_module')>=0))
           .sort('related_genes_module')
           .groupby('related_genes_module', maintain_order=True).apply(lambda x: heatmap_mod(x))
     )

    if grid:
        for ax in iaxes:
            ars = sns.heatmap(np.array([[0,0],[0,0]]), cbar=False, cmap='Greys', vmax=1, vmin=0, square = True, ax=ax)
            ax.tick_params(labelsize=0, width=0, pad=-2, axis='both', which='major')

    plt.tight_layout(pad=-0.125)
    
    if explore:
        print(f"exiting {explore=}")
        return
        raise ValueError("Run again with explore=False")

    # define lateral gene list
    # Do we really want to exclude initially all genes that are related to a given module?
    lateral_gene_names = pl_var.filter(pl.col('related_genes_module').is_in(lateral_modules))['gene_name'].sort().to_list()
    for i in manual_ban:
        lateral_gene_names.append(i)
    if force_usual_suspects:
        for i in usual_suspects:
            rich.print(f":vampire:{i}")
            lateral_gene_names.append(i)

    print(f"{len(lateral_gene_names)=}")
    print(' '.join(lateral_gene_names))
    
    ## TODO any other genes to mark?
    mc.pl.mark.mark_lateral_genes(adata, lateral_gene_names=lateral_gene_names)


def invoque_mc(adata, 
               adata_workdir,
               cpus,
               mc_cpus_key,
               var_cpus_key,
               set_name,
               return_adatas=True,
               max_parallel_piles:int|None=50,
               ):
    if max_parallel_piles is None:
        max_parallel_piles = mc.pl.guess_max_parallel_piles(adata)
    
    print(f"{max_parallel_piles=}")

    mc.pl.set_max_parallel_piles(max_parallel_piles)

    mc.utilities.parallel.set_processors_count(cpus[mc_cpus_key])
    print(f"computing metacells with {cpus[mc_cpus_key]} CPUs")

    mc.pl.divide_and_conquer_pipeline(
            adata,
       #    target_metacell_size=target_metacell_size,
            random_seed=123456)
    
    print(f"conquered")
    metacells = mc.pl.collect_metacells(
            adata, 
            name=f'{set_name}.metacells',
            random_seed=12345)

    #mc.pl.compute_umap_by_features(metacells, \
    #                               max_top_feature_genes=1000, \
    #                               min_dist=5, \
    #                               random_seed=123456, \
    #                               umap_k=15)
    ## add layer with normalized log2 expression
    #lfc = metacells[:, :].to_df()
    #lfc = lfc.apply(lambda x: np.log2(1e-5 + (x/np.sum(x))), axis=1)
    #lfc =  lfc.apply(lambda x: x - np.median(x), axis=0)
    #metacells.layers['lfc'] = lfc 

    plt.rcParams.update(plt.rcParamsDefault)
    

    mc.utilities.parallel.set_processors_count(cpus[var_cpus_key])
    print(f"computing outliers with {cpus[var_cpus_key]} CPUs")
    

    outliers = mc.pl.compute_for_mcview(
        adata=adata, 
        gdata=metacells, 
        random_seed=123456, 
        compute_var_var_similarity=dict(top=50, bottom=50)
        )


    mc.utilities.parallel.set_processors_count(cpus['full'])
    print(f"done, back to {cpus['full']} CPUs")

    ofn_single_cells = f'{adata_workdir}/{set_name}.scells.h5ad'
    ofn_meta_cells = f'{adata_workdir}/{set_name}.mcells.h5ad'
    ofn_metadata = f'{adata_workdir}/{set_name}_metadata.csv'
    
    print('writing h5ds')
    print('\n'.join([set_name, ofn_single_cells, ofn_meta_cells, ofn_metadata]))

    adata.write_h5ad(ofn_single_cells)
    metacells.write_h5ad(ofn_meta_cells)

    metadata_df=None

    if 'batch' in adata.obs.columns:
        # TODO switch to pl
        import pandas as pd
        if len(adata.obs['batch'].unique())>1:
            print('storing metadata')

            batch_counts = adata.obs.groupby('metacell').apply(lambda x: x.batch.value_counts(normalize=False)).unstack()

            print(batch_counts.head())
            batch_frac = batch_counts.apply(lambda x: x/sum(x), axis=1)
            batch_frac.columns = [ f'{i}_frac' for i in batch_frac.columns]

            batch_ln = batch_counts.apply(lambda x: np.log10(1+x), axis=1)
            batch_ln.columns = [ f'{i}_ln' for i in batch_ln.columns]

            metadata_df=pd.concat([batch_ln, batch_frac], axis=1, )
            metadata_df.iloc[2:,:].to_csv(ofn_metadata)
    print('DONE') 

    if return_adatas:
        return([adata, metacells, metadata_df])    
    else:
        print("returning output file names")
        return([ofn_single_cells, ofn_meta_cells, ofn_metadata])

def basure():
    print(f'Original {mc.utilities.parallel.get_processors_count()} CPUs')
    mc.utilities.parallel.set_processors_count(cpus['full'])

    print(f'Running on {mc.utilities.parallel.get_processors_count()} CPUs')
    print(f'set shape {adata.shape}')

    excluded_gene_patterns.append('mt-.*')

    # .var gets masks: bursty_lonely_gene, properly_sampled_gene, excluded_gene
    mc.pl.exclude_genes(adata, 
                        excluded_gene_patterns=excluded_gene_patterns,
                        random_seed = 12345)

    mc.pl.exclude_cells(adata, 
                    properly_sampled_min_cell_total=properly_sampled_min_cell_total, 
                    properly_sampled_max_cell_total=properly_sampled_max_cell_total, 
                    properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)

    full = adata

    total_umis_of_cells = mc.ut.get_o_numpy(full, name='__x__', sum=True)

    # fig total umis per cell histogram (log)
    fg = sns.displot(total_umis_of_cells, bins=800, aspect=3, element='step')
    fg.ax.set(xlabel='UMIs', ylabel='Density', yticks=[])
    fg.ax.axvline(x=properly_sampled_min_cell_total, color='darkgreen')
    fg.ax.axvline(x=properly_sampled_max_cell_total, color='crimson')
    fg.ax.set_xscale('log')
    fg.ax.set_title(f'{set_name}')
    fg.savefig(f'{adata_workdir}/{set_name}_umi_dplot.png')


    # fig total umis of excluded genes
    too_small_cells_count = sum(total_umis_of_cells < properly_sampled_min_cell_total)
    too_large_cells_count = sum(total_umis_of_cells > properly_sampled_max_cell_total)
    
    too_small_cells_percent = 100.0 * too_small_cells_count / len(total_umis_of_cells)
    too_large_cells_percent = 100.0 * too_large_cells_count / len(total_umis_of_cells)
    
    print(f"Will exclude %s (%.2f%%) cells with less than %s UMIs"
            % (too_small_cells_count,
                too_small_cells_percent,
                properly_sampled_min_cell_total))
    print(f"Will exclude %s (%.2f%%) cells with more than %s UMIs"
            % (too_large_cells_count,
                too_large_cells_percent,
                properly_sampled_max_cell_total))

    excluded_genes_data = mc.tl.filter_data(full, var_masks=['excluded_gene'])[0]
    excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
    excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells

    too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)
    too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)
    
    print(f"Will exclude %s (%.2f%%) cells with more than %.2f%% excluded gene UMIs"
            % (too_excluded_cells_count,
                too_excluded_cells_percent,
                100.0 * properly_sampled_max_excluded_genes_fraction))
    
    fg = sns.displot(excluded_fraction_of_umis_of_cells, bins=200, aspect=3, element="step")
    fg.ax.set(xlabel="Fraction of excluded gene UMIs", ylabel='Density', yticks=[])
    fg.ax.axvline(x=properly_sampled_max_excluded_genes_fraction, color='crimson')
    fg.ax.set_title(f'{set_name}')
    fg.savefig(f'{adata_workdir}/{set_name}_fr_excluded.png')


    ###
    clean = mc.pl.extract_clean_data(full)
    mc.pl.relate_genes(clean, random_seed=123456)
    
    print(f'added hardwired suspects to suspect_gene_names')

    for i in usual_suspects:
    #for i in ['shRNA', 'rtTA', 'Fun','Cas9', "Sox4", "Oct4", "Nanog"]:
        print(f'{i}', end= ',')
        suspect_gene_names.append(i)

    print('\n', end ='')

    suspect_gene_patterns = ["Rpl", "Rps", "Mcm", 'Hsp', "Hist", "hist", 'hsp']

    suspect_genes_mask = mc.tl.find_named_genes(clean,
                                                names=suspect_gene_names,
                                                patterns=suspect_gene_patterns)
    suspect_gene_names = sorted(clean.var_names[suspect_genes_mask])
    
    print(f'{len(suspect_gene_names)}\n{" ".join(suspect_gene_names)}') 
    
    module_of_genes = clean.var['related_genes_module']
    suspect_gene_modules = np.unique(module_of_genes[suspect_genes_mask])
    suspect_gene_modules = suspect_gene_modules[suspect_gene_modules >= 0]

    print("suspect_gene_modules")
    for i in sorted(suspect_gene_modules):
        soc = module_of_genes[module_of_genes == i]
        print(f'm{i}::{len(soc)}::\t{"  ".join(sorted(soc.index.to_list()))}') 

    all_modules = [i for i in np.unique(module_of_genes) if int(i) >0]
    sns.set(font_scale=0.65)

    with plt.rc_context({'figure.dpi':'150', 'figure.figsize':(7,6)}):
        similarity_of_genes = mc.ut.get_vv_frame(clean, 'related_genes_similarity')
        for gene_module in suspect_gene_modules:
            fig, ax = plt.subplots(1, 1)
            module_genes_mask = module_of_genes == gene_module
            similarity_of_module = similarity_of_genes.loc[module_genes_mask, module_genes_mask]
            
            similarity_of_module.index = \
            similarity_of_module.columns = [
                '(*) ' + name if name in suspect_gene_names else name
                for name in similarity_of_module.index
            ]
            sns.heatmap(similarity_of_module, xticklabels=1, vmin=0, vmax=1, ax=ax, cmap="YlGnBu", linewidths=0.25)
            lateral_txt ="(ignored)" if (gene_module in lateral_mods and not explore)  else " " 
            ax.set_title(f'{set_name} Gene Module {gene_module} {lateral_txt}')
            fig.savefig(f'{adata_workdir}/{set_name}_module_{gene_module}.png')
            plt.show()
    sns.set(font_scale=1)

    # Do we really want to exclude initially all genes that are related to a given module?
    lateral_genes_mask = suspect_genes_mask
    #lateral_genes_mask = [i in usual_suspects for i in clean.var_names]
    
    for gene_module in suspect_gene_modules:
        print(f'{gene_module=}')
        if gene_module in lateral_mods:
            print(f"gene module: {gene_module}")
            # What is |= ?
            # It's a compound operator, when you say: x |= y it's equivalent to x = x | y
            module_genes_mask = module_of_genes == gene_module
            #lateral_genes_mask = forbiden_genes_mask | module_genes_mask
            lateral_genes_mask |= module_genes_mask
    
    lateral_gene_names = sorted(clean.var_names[lateral_genes_mask])

    if force_usual_suspects:
        for i in usual_suspects:
            rich.print(f":vampire:{i}")
            lateral_gene_names.append(i)

    lateral_gene_names = sorted(list(set(lateral_gene_names)))
    

    print(f"{len(lateral_gene_names)=}")
    print(' '.join(lateral_gene_names))
    
    if explore:
        print("exiting")
        raise ValueError("Run again with explore=False")

    ## TODO any other genes to mark?
    mc.pl.mark.mark_lateral_genes(clean, lateral_gene_names=lateral_gene_names)

    max_parallel_piles = mc.pl.guess_max_parallel_piles(clean)
    print(f"{max_parallel_piles=}")

    mc.pl.set_max_parallel_piles(max_parallel_piles)
    
    mc.utilities.parallel.set_processors_count(cpus[mc_cpus_key])
    print(f"computing metacells with {cpus[mc_cpus_key]} CPUs")

    mc.pl.divide_and_conquer_pipeline(clean,
                       #               target_metacell_size=target_metacell_size,
                                      random_seed=123456)
    
    plt.rcParams.update(plt.rcParamsDefault)
    print(f"conquered")


def return_raw_gene_str():
        return "AY036118 Acta1 Actb Actc1 Actg1 Afp Ahdc1 Aldh1a3 Arhgap28 Arl6ip1 Aspm Bhlhe40 Bmp2 Camk1d Car4 Cas9 Ccn2 Ccnd1 Ccnd2 Cdh11 Cdk8 Cdx2 Cdx4 Cenpf Chchd2 Cited2 Clcn3 Clu Cmss1 Cnn1 Col23a1 Cox6c Cox7c Cp Cped1 Crabp1 Cxcl14 Cyp26a1 Cyp51 D10Wsu102e Dbi Dcc Ddit4 Dkk1 Dlc1 Dlx1 Dlx2 Dmrt2 Dynlt1b Egfem1 Eno1 Erh Fam110b Fat4 Fau Fbn2 Fgf14 Fgf8 Flrt2 Fos Frem1 Fst Gapdh Gas1 Gm10076 Gm10260 Gm19951 Gm20628 Gm32061 Gm42418 Gm45889 Gm53 Gng5 Gpc3 Gphn Gpt2 H1f0 Hebp2 Hist1h1a Hist1h1b Hist1h1c Hist1h1d Hist1h1e Hist1h1t Hist1h2aa Hist1h2ab Hist1h2ac Hist1h2ad Hist1h2ae Hist1h2af Hist1h2ag Hist1h2ah Hist1h2ai Hist1h2ak Hist1h2an Hist1h2ao Hist1h2ap Hist1h2bb Hist1h2bc Hist1h2be Hist1h2bf Hist1h2bg Hist1h2bh Hist1h2bj Hist1h2bk Hist1h2bl Hist1h2bm Hist1h2bn Hist1h2bp Hist1h2bq Hist1h2br Hist1h3a Hist1h3b Hist1h3c Hist1h3d Hist1h3e Hist1h3f Hist1h3g Hist1h3h Hist1h3i Hist1h4a Hist1h4b Hist1h4c Hist1h4d Hist1h4f Hist1h4h Hist1h4i Hist1h4j Hist1h4k Hist1h4m Hist1h4n Hist2h2aa1 Hist2h2ab Hist2h2ac Hist2h2bb Hist2h2be Hist2h3b Hist2h3c1 Hist2h3c2 Hist2h4 Hist3h2a Hist3h2ba Hist4h4 Hmgcr Hmgcs1 Hoxa10 Hoxa11os Hoxa9 Hoxc10 Hsp90aa1 Hsp90ab1 Hsp90b1 Hspa12a Hspa12b Hspa13 Hspa14 Hspa1a Hspa1b Hspa1l Hspa2 Hspa4 Hspa4l Hspa5 Hspa8 Hspa9 Hspb1 Hspb11 Hspb2 Hspb3 Hspb6 Hspb7 Hspb8 Hspb9 Hspbap1 Hspbp1 Hspd1 Hspe1 Hspe1-rs1 Hspg2 Hsph1 Id1 Id3 Idi1 Ier3 Igfbp2 Il17rd Insig1 Jun Krt18 Krt7 Krt8 Lix1 Macf1 Mafb Mcm10 Mcm2 Mcm3 Mcm3ap Mcm4 Mcm5 Mcm6 Mcm7 Mcm8 Mcm9 Mcmbp Mcmdc2 Mest Mif Msgn1 Myl3 Myl4 Myl7 Nckap5 Ndnf Ndufa4 Nes Nkx2-5 Nkx3-1 Nlgn1 Nme2 Nnat Nppa Nrg3 P3h2 Pax1 Pcdh9 Pcna Pcsk5 Pdlim3 Perp Pf4 Pfkfb3 Phlda2 Pim1 Plod2 Polr2l Ppia Ptn Ptprn2 Qk Rbfox2 Rbms1 Rbp4 Reln Rgmb Rmst Rnaset2a Rnf128 Rpl10 Rpl10-ps3 Rpl10a Rpl10l Rpl11 Rpl12 Rpl13 Rpl13a Rpl14 Rpl15 Rpl17 Rpl18 Rpl18a Rpl19 Rpl21 Rpl22 Rpl22l1 Rpl23 Rpl23a Rpl24 Rpl26 Rpl27 Rpl27a Rpl28 Rpl29 Rpl3 Rpl30 Rpl31 Rpl32 Rpl34 Rpl35 Rpl35a Rpl36 Rpl36-ps4 Rpl36a Rpl36a-ps1 Rpl36al Rpl37 Rpl37a Rpl38 Rpl39 Rpl39l Rpl3l Rpl4 Rpl41 Rpl5 Rpl6 Rpl7 Rpl7a Rpl7l1 Rpl8 Rpl9 Rpl9-ps1 Rpl9-ps6 Rplp0 Rplp1 Rplp2 Rps10 Rps11 Rps12 Rps13 Rps14 Rps15 Rps15a Rps16 Rps17 Rps18 Rps19 Rps19bp1 Rps2 Rps20 Rps21 Rps23 Rps24 Rps25 Rps26 Rps27 Rps27a Rps27l Rps27rt Rps28 Rps29 Rps3 Rps3a1 Rps4x Rps5 Rps6 Rps6ka1 Rps6ka2 Rps6ka3 Rps6ka4 Rps6ka5 Rps6ka6 Rps6kb1 Rps6kb2 Rps6kc1 Rps6kl1 Rps7 Rps8 Rps9 Rpsa Rspo3 Sec61b Sec61g Serpinh1 Sfrp1 Sfrp5 Sh3bgr Slc2a1 Slc2a3 Smoc1 Snrpg Sox4 Sox9 Sp5 Srrm2 T Tbx6 Tcf15 Tead1 Tinagl1 Tmod1 Tmsb10 Tmsb4x Tnnc1 Tnni1 Tnni3 Tomm6 Top2a Trim30a Tuba1a Tuba1b Tuba1c Uba52 Ube2c Uncx Utrn Wfdc1 Wls Wnt6 Zfp36l1 Zfp36l2 rtTA shRNA"

