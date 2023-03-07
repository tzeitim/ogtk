
def len_intersect(adata, rs):
    good_cbs = set([i.replace('-1', '') for i in adata.obs.index])
    lib_cbs = set(rs.umis.keys())
    iset = len(good_cbs.intersection(lib_cbs))
    return(np.array([iset, iset/adata.obs.shape[0],iset/len(rs.umis.keys()) ]))

def print_stats(adata, rs):
    iset = len_intersect(adata, rs)
    print(f"A total of {iset[0]} cbs were found on both sets {iset[1]:0.2f} {iset[2]:0.2f} ")

def metacellize (set_name, 
                 adata, 
                 adata_workdir, 
                 excluded_gene_patterns = [], 
                 excluded_gene_names=None, 
                 suspect_gene_names=['shRNA', 'rtTA', 'Fun','Cas9', 'hist1h2ap', "Hist1h2ap"], 
                 forbidden_mods = [18, 26, 27, 34], 
                 return_adatas=True, 
                 log_debug=False, 
                 explore=True, 
                 full_cpus=56, 
                 moderate_cpus=8, 
                 properly_sampled_max_excluded_genes_fraction = 0.03):
    import anndata as ad
    import matplotlib.pyplot as plt
    import metacells as mc
    if log_debug:
        import logging
        mc.ut.setup_logger(level=logging.DEBUG)

    import numpy as np
    if log_debug:
        print(np.show_config())
    import os
    import pandas as pd
    import scipy.sparse as sp
    import regex
    
    import metacells.utilities.typing as utt
    
    import scanpy as sc
    import seaborn as sns
    sns.set_theme()
    #sns.set_style("whitegrid")
    #from rich import print 
    sns.set_context("talk")

    utt.sum_duplicates(adata.X)
    utt.sort_indices(adata.X)

    mc.ut.set_name(adata, set_name)

    #full_cpus = 24 #56
    #moderate_cpus = 8


    print(f'Original {mc.utilities.parallel.get_processors_count()} CPUs')
    mc.utilities.parallel.set_processors_count(full_cpus)

    print(f'Running on {mc.utilities.parallel.get_processors_count()} CPUs')
    print(f'set shape {adata.shape}')

    excluded_gene_patterns.append('mt-.*')

    mc.pl.analyze_clean_genes(adata,
                          excluded_gene_names=excluded_gene_names,
                          excluded_gene_patterns=excluded_gene_patterns,
                          random_seed=12345)

    print("picking clean genes")
    mc.pl.pick_clean_genes(adata)

    full = adata

    properly_sampled_min_cell_total = 500
    properly_sampled_max_cell_total = 20000

    total_umis_of_cells = mc.ut.get_o_numpy(full, name='__x__', sum=True)

    with plt.rc_context({"figure.figsize":(15,5)}):
        fig, ax = plt.subplots(1, 1)
        plot = sns.distplot(total_umis_of_cells, bins=200,  )
        plot.set(xlabel='UMIs', ylabel='Density', yticks=[])
        plot.axvline(x=properly_sampled_min_cell_total, color='darkgreen')
        plot.axvline(x=properly_sampled_max_cell_total, color='crimson')
        plot.set_xlim((0, 20000))    

        ax.set_title(f'{set_name}')
        fig.savefig(f'{adata_workdir}/{set_name}_umi_dplot.png')

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

    with plt.rc_context({"figure.figsize":(15,5)}):
        fig, ax = plt.subplots(1, 1)
        
        excluded_genes_data = mc.tl.filter_data(full, var_masks=['~clean_gene'])[0]
        excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
        excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells
        
        plot = sns.distplot(excluded_fraction_of_umis_of_cells, bins =200)
        plot.set(xlabel='Fraction of excluded gene UMIs', ylabel='Density', yticks=[])
        plot.axvline(x=properly_sampled_max_excluded_genes_fraction, color='crimson')

        ax.set_title(f'{set_name}')
        fig.savefig(f'{adata_workdir}/{set_name}_fr_excluded.png')

    too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)
    
    too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)
    
    print(f"Will exclude %s (%.2f%%) cells with more than %.2f%% excluded gene UMIs"
            % (too_excluded_cells_count,
                too_excluded_cells_percent,
                100.0 * properly_sampled_max_excluded_genes_fraction))

    print("Analyzing clean cells")
    mc.pl.analyze_clean_cells(
        full,
        properly_sampled_min_cell_total=properly_sampled_min_cell_total,
        properly_sampled_max_cell_total=properly_sampled_max_cell_total,
        properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)
    
    
    mc.pl.pick_clean_cells(full)
    
    clean = mc.pl.extract_clean_data(full)
    
    print(f'added hardwired suspects to suspect_gene_names')
    usual_suspects =['shRNA', 'rtTA', 'Fun', 'Cas9', 'Pcna', 'Jun', 'Top2a', 'Txn', 'Hsp90ab1', 'Fos'] 
    for i in usual_suspects:
    #for i in ['shRNA', 'rtTA', 'Fun','Cas9', "Sox4", "Oct4", "Nanog"]:
        print(f'{i}', end= ',')
        suspect_gene_names.append(i)

    print('\n', end ='')

    suspect_gene_patterns = ["Rpl", "Rps", "Mcm", 'Hsp', "Hist", "hist", 'hsp']
    suspect_genes_mask = mc.tl.find_named_genes(clean, names=suspect_gene_names,
                                                patterns=suspect_gene_patterns)
    suspect_gene_names = sorted(clean.var_names[suspect_genes_mask])
    
    print(f'{len(suspect_gene_names)}\n{" ".join(suspect_gene_names)}') 
    
    clean.var.reset_index()['index'].value_counts().value_counts()
    
    mc.pl.relate_genes(clean, random_seed=123456)
    
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
            #similarity_of_module.columns = lll
            sns.heatmap(similarity_of_module, xticklabels=1, vmin=0, vmax=1, ax=ax, cmap="YlGnBu", linewidths=0.25)
            forbidden_txt ="(ignored)" if (gene_module in forbidden_mods and not explore)  else " " 
            ax.set_title(f'{set_name} Gene Module {gene_module} {forbidden_txt}')
            fig.savefig(f'{adata_workdir}/{set_name}_module_{gene_module}.png')
            plt.show()
    sns.set(font_scale=1)

    # Do we really want to exclude initially all genes that are related to a given module?
    forbidden_genes_mask = suspect_genes_mask
    #forbidden_genes_mask = [i in usual_suspects for i in clean.var_names]
    
    for gene_module in suspect_gene_modules:
        if gene_module in forbidden_mods:
            print(f"gene module: {gene_module}")
            # What is |= ?
            # It's a compound operator, when you say: x |= y it's equivalent to x = x | y
            module_genes_mask = module_of_genes == gene_module
            #forbidden_genes_mask = forbiden_genes_mask | module_genes_mask
            forbidden_genes_mask |= module_genes_mask

    forbidden_gene_names = sorted(clean.var_names[forbidden_genes_mask])

    print(f"forbidden genes {len(forbidden_gene_names)}")
    print(' '.join(forbidden_gene_names))
    
    if explore:
        return()
        
    max_parallel_piles = mc.pl.guess_max_parallel_piles(clean)
    print(f"max parallel piles {max_parallel_piles}")

    mc.pl.set_max_parallel_piles(max_parallel_piles)
    mc.pl.divide_and_conquer_pipeline(clean,
                                      forbidden_gene_names=forbidden_gene_names,
                                      #target_metacell_size=...,
                                      random_seed=123456)
    
    print(f"conquered")
    metacells = mc.pl.collect_metacells(clean, name=f'{set_name}.metacells')

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

    mc.utilities.parallel.set_processors_count(moderate_cpus)
    print(f"computing outliers with {moderate_cpus} CPUs")

    outliers = mc.pl.compute_for_mcview(
        adata=clean, 
        gdata=metacells, 
        random_seed=123456, 
        compute_var_var_similarity=dict(top=50, bottom=50)
        )

    mc.utilities.parallel.set_processors_count(full_cpus)
    print(f"done, back to {full_cpus} CPUs")

    ofn_single_cells = f'{adata_workdir}/{set_name}.scells.h5ad'
    ofn_meta_cells = f'{adata_workdir}/{set_name}.mcells.h5ad'
    ofn_outliers = f'{adata_workdir}/{set_name}.outliers.h5ad'
    ofn_metadata = f'{adata_workdir}/{set_name}_metadata.csv'
    
    print('writing h5ds')
    print('\n'.join([set_name, ofn_single_cells, ofn_meta_cells, ofn_outliers, ofn_metadata]))

    clean.write_h5ad(ofn_single_cells)
    metacells.write_h5ad(ofn_meta_cells)
    outliers.write_h5ad(ofn_outliers)

    metadata_df=None
    if 'batch' in clean.obs.columns:
        print('storing metadata')

        batch_counts = clean.obs.groupby('metacell').apply(lambda x: x.batch.value_counts(normalize=False)).unstack()

        print(batch_counts.head())
        batch_frac = batch_counts.apply(lambda x: x/sum(x), axis=1)
        batch_frac.columns = [ f'{i}_frac' for i in batch_frac.columns]

        batch_ln = batch_counts.apply(lambda x: np.log10(1+x), axis=1)
        batch_ln.columns = [ f'{i}_ln' for i in batch_ln.columns]

        metadata_df=pd.concat([batch_ln, batch_frac], axis=1, )
        metadata_df.iloc[2:,:].to_csv(ofn_metadata)
    print('DONE') 

    if return_adatas:
        return([clean, metacells, outliers, metadata_df])    
    else:
        return([ofn_single_cells, ofn_meta_cells, ofn_outliers, ofn_metadata])
    col = metacells[:, :].to_df()
    plt.figure()

    umap_x = mc.ut.get_o_numpy(metacells, 'umap_x')
    umap_y = mc.ut.get_o_numpy(metacells, 'umap_y')
    plot = sns.scatterplot(x=umap_x, y=umap_y)
    plt.figure()
    umap_x = mc.ut.get_o_numpy(metacells, 'umap_x')
    umap_y = mc.ut.get_o_numpy(metacells, 'umap_y')
    plot = plt.scatter(x=umap_x, y=umap_y, c =col['Cas9'], cmap ='bwr' )
    return([clean, metacells])
