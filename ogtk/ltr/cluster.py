import matplotlib.pyplot as plt
import polars as pl
import hdbscan
import numpy as np
import seaborn as sns
from .. import UM
from .. import utils as ut
import anndata as ad 
import os


class Exp():
    def __init__(self, name, input_dir, adh5_path, anchor3 ="AATCCAGCTAGCTGTGCAG", anchor5="TCACTGTCTTCCGCAAT"):
        print(f"{name}")
        self.name = name
        self.input_dir = input_dir
        self.adh5_path = adh5_path
        self.anchor3 = anchor3 
        self.anchor5 = anchor5 
        if os.path.exists(self.adh5_path):
            print(f"loading valid cells from {self.adh5_path}")
            self.valid_cells = ad.read_h5ad(self.adh5_path).obs_names
            print(f"{len(self.valid_cells)=}")
            
        self.fq_files = ut.sfind(self.input_dir, pattern = "*_R1_*fastq.gz")
        
    def _pp(self):
        # took around 50 min and 283 GB of RAM. one core
        self.parquets = []
        for file in self.fq_files:
            out_fn = f"{file}".replace('.fastq.gz', '.mols.parquet')
            self.parquets.append(
                ut.tabulate_paired_10x_fastqs_rs(
                file_path=file,
                out_fn=out_fn,
                force=True,
                cbc_len=16,
                umi_len=12,
                modality='single-cell',
                do_rev_comp=False,
            ))

    def _ls_pq(self):
        self.parquets = ut.sfind(self.input_dir, pattern = "*_R1_*.mols.parquet")
        return self.parquets
        ifn = "/home/projects/nyosef/pedro/projects/wlt/workdir/20240515/caspp/05_mol_table.parquet"
    
    def load_cas(self, ifn, min_clone_size=0.01, intbc_filter_expr=pl.col('intBC').str.len_chars()==14):
        print(f"Keeping only clones as little as {min_clone_size:.2%}")
        print(f"Filtering intBCs by {intbc_filter_expr}")
        self.cdf = (
            pl.scan_parquet(ifn)
            .with_columns(ot=pl.col('Seq').str.contains(self.anchor3) & pl.col('Seq').str.contains(self.anchor5))
            .with_columns(csize=pl.col('UMI').n_unique().over('cellBC'))
            .with_columns(pl.col('Seq').str.extract(f'CTAGCTGTGCAGC(.+?)ATTCAACTGCAGT', 1).alias('intbc'))
            .with_columns(clone_size=pl.col('cellBC').n_unique().over('intBC')/pl.col('cellBC').n_unique())
            .with_columns(cl_size=pl.col('cellBC').n_unique().over('intbc')/pl.col('cellBC').n_unique())
            .filter(pl.col('clone_size')>=min_clone_size)
            .filter(intbc_filter_expr)
            .with_columns(pl.col('UMI').n_unique().over('allele', 'cellBC', 'intBC').alias('allele_umis'))
            .collect()
        )
        
    def cluster_cas(self, min_mols_allele, jobs=1, clusterer=None):
        cl = cluster(self.cdf, 
                          name = f"{self.name}_cass", 
                          min_mols_allele=min_mols_allele, 
                          jobs=jobs, 
                          clusterer=clusterer)
        
        self.cl = self.cdf.join(cl, left_on="intbc", right_on="int_bc")

        sns.catplot(
                data=self.cl.group_by('cellBC')
                .agg(pl.col('cluster').n_unique())
                ['cluster']
                .value_counts()
                .sort('cluster')
                .to_pandas(), 
            kind='bar',
            x='cluster', 
            y='count')

        plt.title("clusters per cell")

        self.annotate_cells()

        sns.catplot(data=self.ad.obs.groupby('top_cluster', observed=True).size().reset_index(name="cells"), 
                    x='top_cluster', 
                    y="cells", 
                    kind='bar')
        plt.title("cells per cluster")
        plt.show()
        
        sns.catplot(data=self.cl.group_by('cellBC', 'cluster')
                    .len()
                    .sort(pl.col('len'))
                    .group_by("cellBC", maintain_order=True)
                    .head(2)
                    .with_columns(pl.col('len').rank(descending=True).over('cellBC').alias('rank'))
                    .to_pandas(),
                    x='rank', 
                    y='len', 
                    #log_scale=10, 
                    kind='boxen'

            )       
        plt.show()

        vs_1_2 =(
            self.cl
                .group_by('cellBC', 'cluster').len().sort(pl.col('len'))
                .group_by("cellBC", maintain_order=True).head(2)
                .with_columns(pl.col('len').rank(descending=True).over('cellBC').alias('rank'))
                .pivot(index='cellBC', on='rank', values='len', aggregate_function=pl.sum('len').log10(),)
                .fill_null(0)
        )

        sns.displot(data = vs_1_2.to_pandas(),
            x='1.0', y='2.0', bins=(50,50), cmap='magma',
        )
        plt.show()

        self.vs_1_2 = vs_1_2
    
    def load_ad(self):
        import anndata as ad
        self.ad = ad.read_h5ad(self.adh5_path)

    def annotate_cells(self):
        self.load_ad()
        self.ad.obs = (
        self.ad.obs.merge(
            right=self.cl.group_by('cellBC', 'cluster').len()
                    .group_by("cellBC").agg(pl.col('cluster').top_k_by(by='len', k=1))
                        .explode('cluster')
                        .with_columns(pl.col('cluster').cast(pl.Utf8))
                        .rename({'cluster':'top_cluster'})
                        .with_columns(pl.col('cellBC')+'-1')
                        
                    .to_pandas()
                    .set_index('cellBC'),
            left_index=True, 
            right_index=True, 
            how='left')
        )

        #allclusters 
        self.ad.obs = (
            self.ad.obs.merge(
                right=self.cl.group_by('cellBC')
                            .agg(pl.col('cluster').unique().sort().cast(pl.Utf8))
                            .with_columns(pl.col('cluster').list.join("-"),
                                        pl.col('cluster').list.len().cast(pl.Utf8).alias('n_clusters'))
                            .with_columns(pl.col('cellBC')+'-1')
                        .to_pandas()
                        .set_index('cellBC'),
                left_index=True,
                right_index=True,
                how='left')
        )

        import scanpy as sc

        genes = [i for i in 
                 ['Cas9', 'Mcherry-lt', 'Hygr', 'Egfp-dox', 'Tet-on-bsd', 'mCherry', 'Crispr'] 
                 if i in self.ad.var_names
                 ]
        print(genes)
        sc.pl.umap(self.ad, 
                   color=np.hstack([genes, ['pct_counts_mt', 'top_cluster', 'n_clusters', ]]), 
                   color_map='magma', 
                   #s=20, 
                   ncols=3)

        sc.pl.violin(self.ad, 
                     groupby='top_cluster', 
                     keys=genes,
                     legend=False)


def cluster(df, name, min_mols_allele=0, jobs=10, clusterer=None):
    '''
    Expects intBC cellBC allele 
    '''
    if clusterer is None:
        clusterer = hdbscan.HDBSCAN(core_dist_n_jobs=50)

    data = (
        #df.group_by('cellBC', 'intBC', 'allele').len().sort('')
        df
        .group_by('cellBC', 'intBC', 'allele')
        .len(name='allele_umis')
        .sort('allele_umis')
        .filter(pl.col('allele_umis')>min_mols_allele)
        .filter(pl.col('intBC').str.len_chars() == 14)
    )


    intbc_cov = (
        df
        .filter(pl.col('allele_umis')>min_mols_allele)
        .sort('allele_umis', descending=True)
        .group_by('intBC', 'cellBC', maintain_order=True).head(1)
        .select('intBC', 'cellBC', 'allele_umis').unique().sort('allele_umis')
        .group_by('intBC').agg(mols_intbc=pl.col('allele_umis').sum())
        .sort('mols_intbc')
    )

    print(intbc_cov.sort('mols_intbc'))

    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(17.5, 5))
    sns.heatmap(UM.hdist_all(intbc_cov.sort('intBC')['intBC'].sort(), jobs = jobs), vmax=14, ax=ax1)
    ax1.set_title("hamming distance (sorted alphabetically)")

    sns.heatmap(UM.hdist_all(intbc_cov.sort('mols_intbc')['intBC'], jobs = jobs), vmax =14, ax=ax2)
    ax2.set_title("hamming distance (sorted frequency)")

    sns.heatmap(UM.hdist_all(intbc_cov.sort('mols_intbc')['intBC'][0:50][0:50], jobs = jobs), vmax =14, ax=ax3)
    ax3.set_title("hamming distance (sorted frequency; top 50)")
    fig.savefig(f'{name}_024_hdist_x3.png', bbox_inches='tight')
    plt.show()

    print(data)
    sns.displot(
        data=data.select('intBC', 'allele_umis').unique().to_pandas(),
        x='allele_umis'
        , log_scale=10
    )

    plt.title("molecules per allele")
    plt.grid()
    plt.savefig(f'{name}_07_dist_mols_per_int.png', bbox_inches='tight')
    plt.tight_layout()
    plt.show()
    
    zz = (data      
              .pivot(
                  values='allele_umis', 
                  aggregate_function='sum',
                  index=['cellBC'], 
                  on='intBC')
              .fill_null(0)
         )
    
    z = zz.drop('cellBC')
    z = z.select(np.array(z.columns)[np.sum(z.to_numpy(), axis=0)>0])
    zn = z.transpose().fill_null(0).to_numpy()/1.0
    zn = zn[np.nansum(zn, axis=1)>0]
    
    print("cell x intbc")
    print(zz)
    
    candidate = zz.drop(['cellBC']).columns
    
    cor_mat = np.corrcoef(zn)
    
    print(f"{cor_mat.shape=}")
    print(f"{zn.shape=}")

    from matplotlib.colors import SymLogNorm

    print(zn.min())
    print(zn.max())
    pcm=plt.pcolormesh(zn, norm=SymLogNorm(linthresh=0.01), shading='auto')
    plt.colorbar(pcm, label='molecules')
    plt.xlabel('cells')
    plt.ylabel('intBCs')

    plt.title("raw counts")
    #plt.grid()
    plt.savefig(f'{name}_08_intbc_vs_cell_raw_counts.png', bbox_inches='tight')
    plt.show()
    
    
    clusterer.fit(cor_mat)
    
    labs = clusterer.labels_
    pl.Series(labs).value_counts(sort=True)
    
    cluster_labs =pl.DataFrame(dict(int_bc=candidate, cluster=labs)) 
    yyy = df.join(pl.DataFrame(dict(int_bc=candidate, cluster=labs)), left_on='intBC', right_on='int_bc')
    final_df = (
            yyy
            .group_by(['intBC', 'cluster'])
            .agg(pl.len())
            #.drop('allele_umis')
            #.with_columns(pl.col('sample_id').str.replace('e.?$', '')).rename({'sample_id':'clone'})
            # name clusters based on size. -1 remains as noise and the largest cluster starts with 1
            .with_columns(pl.len().over('cluster').rank('dense', descending=True).alias('rank'))
            .with_columns(cluster=pl.when(pl.col('cluster')==-1).then('cluster').otherwise(pl.col('rank')))
            .drop('rank')
        )
    
    slabs = np.argsort(labs)
    
    fig, (heat, lines) = plt.subplots(1,2, figsize=(10, 5))
    cmap= heat.pcolormesh(cor_mat, cmap='RdYlBu_r')
    lines.scatter(labs, range(len(labs)), s=1)
    lines.set_ylim(0, len(labs))
    lines.set_xlim(min(labs)-0.5, max(labs)+0.5)

    plt.colorbar(cmap, label='PCC')
    plt.title("unsorted cor mat")
    #plt.grid()

    plt.savefig(f'{name}_09_cormat_unsorted.png', bbox_inches='tight')
    plt.show()

    
    slabs = np.argsort(labs)
    
    fig, (heat, lines) = plt.subplots(1,2, figsize=(10, 5))
    cmap=heat.pcolormesh(cor_mat[slabs,:][:,slabs], cmap='RdYlBu_r')
    lines.scatter(labs[slabs], range(len(labs)), s=1)
    lines.set_ylim(0, len(labs))
    lines.set_xlim(min(labs)-0.5, max(labs)+0.5)
    
    plt.colorbar(cmap, label='PCC')
    plt.title("clustered cor mat")
    #plt.grid()

    plt.savefig(f'{name}_10_cormat_HDBSCAN.png', bbox_inches='tight')
    plt.show()

    clustergrid = sns.clustermap(cor_mat, method='ward', cmap='RdYlBu_r', figsize=(7, 7))
    row_order = clustergrid.dendrogram_row.reordered_ind
    col_order = clustergrid.dendrogram_col.reordered_ind
    
    # Reorder the original data
    reordered_cor_mat = cor_mat[row_order][:, col_order]
    plt.title("hclust")
    #plt.grid()
    plt.savefig(f'{name}_11_cormat_hcluster.png', bbox_inches='tight')
    plt.show()

    sns.heatmap(UM.hdist_all(pl.Series(candidate)[slabs]), vmax=4)
    
    plt.title("Hamming HDBSCAN")
    plt.savefig(f'{name}_12_hdist_HDBSCAN.png', bbox_inches='tight')
    plt.show()
    
    sns.heatmap(UM.hdist_all(pl.Series(candidate)[row_order]), vmax=4)
    plt.title("Hamming hclust")
    plt.savefig(f'{name}_13_hdist_hclust.png', bbox_inches='tight')
    plt.show()
    
    return cluster_labs
