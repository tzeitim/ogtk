
def len_intersect(adata, rs):
    good_cbs = set([i.replace('-1', '') for i in adata.obs.index])
    lib_cbs = set(rs.umis.keys())
    iset = len(good_cbs.intersection(lib_cbs))
    return(np.array([iset, iset/adata.obs.shape[0],iset/len(rs.umis.keys()) ]))

def print_stats(adata, rs):
    iset = len_intersect(adata, rs)
    print(f"A total of {iset[0]} cbs were found on both sets {iset[1]:0.2f} {iset[2]:0.2f} ")

def metacellize (adata, excluded_gene_patterns = [], excluded_gene_names= None, suspect_gene_names=[], forbidden_mods = [18, 26, 27, 34]):
    import anndata as ad
    import matplotlib.pyplot as plt
    import metacells as mc
    import numpy as np
    import os
    import pandas as pd
    import scipy.sparse as sp
    import regex
    
    from math import hypot
    from matplotlib.collections import LineCollection
    from IPython.display import set_matplotlib_formats
    
    import metacells.utilities.typing as utt
    
    from importlib import reload
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    import scanpy as sc
    import pandas as pd
    import seaborn as sns
    import seaborn as sb
    sns.set_theme()
    #sns.set_style("whitegrid")
    
    sns.set_context("talk")
    #set_matplotlib_formats('svg')
    #sb.set_style("white")

    utt.sum_duplicates(adata.X)
    utt.sort_indices(adata.X)

    mc.ut.set_name(adata, 'scv2')

    print(adata.shape)

    excluded_gene_patterns.append('mt-.*')

    mc.pl.analyze_clean_genes(adata,
                          excluded_gene_names=excluded_gene_names,
                          excluded_gene_patterns=excluded_gene_patterns,
                          random_seed=12345)

    mc.pl.pick_clean_genes(adata)

    adata_workdir = '/local/users/polivar/artlinet/workdir/scv2_c3_unind'


    adata.write(f'{adata_workdir}/mc2.full')

    full = adata

    properly_sampled_min_cell_total = 500
    properly_sampled_max_cell_total = 20000

    total_umis_of_cells = mc.ut.get_o_numpy(full, name='__x__', sum=True)

    with plt.rc_context({"figure.figsize":(15,5)}):
        plt.figure()
        plot = sb.distplot(total_umis_of_cells, bins=200)
        plot.set(xlabel='UMIs', ylabel='Density', yticks=[])
        plot.axvline(x=properly_sampled_min_cell_total, color='darkgreen')
        plot.axvline(x=properly_sampled_max_cell_total, color='crimson')
        plot.set_xlim((0, 20000))    
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
        plt.figure()
        properly_sampled_max_excluded_genes_fraction = 0.03
        
        excluded_genes_data = mc.tl.filter_data(full, var_masks=['~clean_gene'])[0]
        excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
        excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells
        
        plot = sb.distplot(excluded_fraction_of_umis_of_cells, bins =200)
        plot.set(xlabel='Fraction of excluded gene UMIs', ylabel='Density', yticks=[])
        plot.axvline(x=properly_sampled_max_excluded_genes_fraction, color='crimson')
        #plot.savefig('ex1.pdf')
        too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)
        
        too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)
        
        print(f"Will exclude %s (%.2f%%) cells with more than %.2f%% excluded gene UMIs"
              % (too_excluded_cells_count,
                 too_excluded_cells_percent,
                 100.0 * properly_sampled_max_excluded_genes_fraction))
    
    mc.pl.analyze_clean_cells(
        full,
        properly_sampled_min_cell_total=properly_sampled_min_cell_total,
        properly_sampled_max_cell_total=properly_sampled_max_cell_total,
        properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)
    
    
    mc.pl.pick_clean_cells(full)
    
    clean = mc.pl.extract_clean_data(full)
    clean
    
    clean.var_names
    
    for i in ['shRNA', 'rtTA', 'Fun','Cas9', "Sox4", "Oct4", "Nanog"]:
        suspect_gene_names.append(i)
    suspect_gene_patterns = ["Rpl", "Rps", "Mcm", 'Hsp', "Hist"]
    suspect_genes_mask = mc.tl.find_named_genes(clean, names=suspect_gene_names,
                                                patterns=suspect_gene_patterns)
    suspect_gene_names = sorted(clean.var_names[suspect_genes_mask])
    
    print(f'{len(suspect_gene_names)}\n{" ".join(suspect_gene_names)}') 
    
    clean.var.reset_index()['index'].value_counts().value_counts()
    
    mc.pl.relate_genes(clean, random_seed=123456)
    
    

    module_of_genes = clean.var['related_genes_module']
    suspect_gene_modules = np.unique(module_of_genes[suspect_genes_mask])
    suspect_gene_modules = suspect_gene_modules[suspect_gene_modules >= 0]
    print(suspect_gene_modules)
    print(module_of_genes)
    
    all_modules = [i for i in np.unique(module_of_genes) if int(i) >0]
    
    with plt.rc_context({'figure.dpi':'100', 'figure.figsize':(7,6)}):
        plt.figure()

        similarity_of_genes = mc.ut.get_vv_frame(clean, 'related_genes_similarity')
        for gene_module in suspect_gene_modules:
            module_genes_mask = module_of_genes == gene_module
            similarity_of_module = similarity_of_genes.loc[module_genes_mask, module_genes_mask]
            
            similarity_of_module.index = \
            similarity_of_module.columns = [
                '(*) ' + name if name in suspect_gene_names else name
                for name in similarity_of_module.index
            ]
            #similarity_of_module.columns = lll
            ax = plt.axes()
            sb.heatmap(similarity_of_module, vmin=0, vmax=1, ax=ax, cmap="YlGnBu")
            ax.set_title(f'Gene Module {gene_module}')
            plt.show()
    
    forbidden_genes_mask = suspect_genes_mask
    forbidden_genes_mask = forbidden_genes_mask == 'reset the vector'
    for gene_module in suspect_gene_modules:
        #if gene_module not in [1, 7, 13, 16, 29, 38, 41, 36]:
        if gene_module not in forbidden_mods:
            print(gene_module)
            module_genes_mask = module_of_genes == gene_module
            forbidden_genes_mask |= module_genes_mask
    forbidden_gene_names = sorted(clean.var_names[forbidden_genes_mask])
    print(len(forbidden_gene_names))
    print(' '.join(forbidden_gene_names))
    
    max_parallel_piles = mc.pl.guess_max_parallel_piles(clean)
    print(max_parallel_piles)
    mc.pl.set_max_parallel_piles(max_parallel_piles)
    
    len(forbidden_gene_names)
    
    
    mc.pl.divide_and_conquer_pipeline(clean,
                                      forbidden_gene_names=forbidden_gene_names,
                                      #target_metacell_size=...,
                                      random_seed=123456)
    
    metacells = mc.pl.collect_metacells(clean, name='c3.metacells')
    
    mc.pl.compute_umap_by_features(metacells, \
                                   max_top_feature_genes=1000, \
                                   min_dist=5, \
                                   random_seed=123456, \
                                   umap_k=15)
    # add layer with normalized log2 expression
    lfc = metacells[:, :].to_df()
    lfc = lfc.apply(lambda x: np.log2(1e-5 + (x/np.sum(x))), axis=1)
    lfc =  lfc.apply(lambda x: x - np.median(x), axis=0)
    metacells.layers['lfc'] = lfc 
    
    return([clean, metacells])    
    col = metacells[:, :].to_df()
    plt.figure()

    umap_x = mc.ut.get_o_numpy(metacells, 'umap_x')
    umap_y = mc.ut.get_o_numpy(metacells, 'umap_y')
    plot = sb.scatterplot(x=umap_x, y=umap_y)
    plt.figure()
    umap_x = mc.ut.get_o_numpy(metacells, 'umap_x')
    umap_y = mc.ut.get_o_numpy(metacells, 'umap_y')
    plot = plt.scatter(x=umap_x, y=umap_y, c =col['Cas9'], cmap ='bwr' )
    return([clean, metacells])