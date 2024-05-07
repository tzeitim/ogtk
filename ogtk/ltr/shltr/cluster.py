import numpy as np
import seaborn as sns
import polars as pl
from scipy.cluster.vq import vq, kmeans, whiten
import umap
from sklearn.preprocessing import StandardScaler
import hdbscan
import fastcluster
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import ogtk
from ogtk.utils import db as db 

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()

#__all__ = [
#    "cluster_ibars",
#]

pl.Config.set_fmt_str_lengths(150)

@ogtk.utils.log.call
def load_parquet_for_analysis(df, plot = True, min_cov_log10=2):
    '''
    Computes the marginal coverage of each integration across samples and
    filters out integrations with low coverage. The marginal coverage is the
    sum of the number of molecules across all samples. min_cov_log10 is the
    minimum log10 of the marginal coverage to keep the integration.

    If plot is True, it plots the distribution of marginal coverage and the
    minimum coverage threshold.

    Returns a polars data frame with the marginal coverage and the log10 of the
    marginal coverage.
    '''

    logger.critical("Warning: assuming that molecules provided come from uncut samples") #pyright: ignore

    df = (
        df
        .lazy()
        # TODO: this is a dirty way of keeping sane spacers 
        # it works for now since there are plenty of uncut alleles in general
        .filter(pl.col('wt'))
        # generate a compsite ibar-spacer 
        .with_columns(raw_ibar=pl.col('raw_ibar')+'-'+pl.col('spacer'))

        .groupby(['cbc','raw_ibar','seq', 'sample_id'])
         .agg([
            pl.col('umi_dom_reads').sum().alias('allele_reads'),
            pl.count().alias('allele_umis'),
         ])
        .sort('allele_umis', descending=True)
        .groupby(['cbc','raw_ibar','seq'], maintain_order=True)
        .first()
        .filter(pl.col('allele_reads')>1)
        #.filter(    #pl.int_range(0, pl.count()).shuffle().over("sample_id")) < molecules_per_sample
        .with_columns(pl.col('cbc')+"_"+pl.col('sample_id'))
        .groupby(['cbc', 'sample_id', 'raw_ibar'])
        .agg(pl.col('allele_umis').sum())
        .sort('allele_umis')
        .with_columns(pl.col('allele_umis').sum().over(['raw_ibar', 'sample_id']).alias('ibar_marginal'))
        .with_columns(pl.col('ibar_marginal').log10().alias('ibar_marginal_log'))
        .collect()
    )

    if plot:
        fg = sns.displot(
            data=df.select('raw_ibar', 'ibar_marginal_log', 'sample_id').unique().to_pandas(),
            x='ibar_marginal_log', 
            bins=200, 
            aspect=3, 
            height=3,
            )
        fg.ax.set_yscale('log')
        fg.ax.axvline(x=min_cov_log10, c='red')
        plt.show()
        

    return df.filter(pl.col('ibar_marginal_log')>=min_cov_log10)


@ogtk.utils.log.call
def cluster_ibars(
        mols_parquet: None| str=None, 
        df: None| pl.DataFrame= None,
        plot=True, 
        min_cov_log10=2, 
        umap:bool=False, 
        out_parquet:str|None=None,
        )->pl.DataFrame:
    '''
    Runs HDBSCAN in order to cluster integrations based in the correlation
    matrix of individual integrations and their molecule counts. 

    Input: Either a parquet file path or a polars data frame of single cell
    data, at the molecule level. Setting umap to True can be very slow for
    large datasets so it is not recommended. Output: Returns 

    out_parquet = '/local/users/polivar/src/artnilet/workdir/scv2/ibar_clusters.parquet',
    '''

    assert not (mols_parquet is None and df is None), "You need to specify mols_parque or df"

    if df is None:
        df=pl.scan_parquet(mols_parquet)

    df = load_parquet_for_analysis(df, plot, min_cov_log10)

    zz = (df
          .pivot(
              values='allele_umis', 
              aggregate_function='sum',
              index=['sample_id', 'cbc'], 
              columns='raw_ibar')
          .fill_null(0)
         )

    z = zz.drop(['sample_id', 'cbc'])
    
    z = z.select(np.array(z.columns)[np.sum(z.to_numpy(), axis=0)>0])
    zn = z.transpose().fill_null(0).to_numpy()/1.0
    zn = zn[np.nansum(zn, axis=1)>0]

    candidate = zz.drop(['sample_id', 'cbc']).columns
    cor_mat = np.corrcoef(zn)
    
    clusterer = hdbscan.HDBSCAN(core_dist_n_jobs=50)
    clusterer.fit(cor_mat)

    labs = clusterer.labels_
    yyy = df.join(pl.DataFrame(dict(ibar=candidate, cluster=labs)), left_on='raw_ibar', right_on='ibar')


    final_df = (
            yyy
            .groupby(['sample_id', 'raw_ibar', 'cluster'])
            .agg(pl.count())
            .drop('count')
            .with_columns(pl.col('sample_id').str.replace('e.?$', '')).rename({'sample_id':'clone'})
            # name clusters based on size. -1 remains as noise and the largest cluster starts with 1
            .with_columns(pl.count().over('cluster').rank('dense', descending=True).alias('rank'))
            .with_columns(cluster=pl.when(pl.col('cluster')==-1).then('cluster').otherwise(pl.col('rank')))
            .drop('rank')
            # split composite ibar-spacer into two individual fields
            .with_columns(
                raw_ibar=pl.col('raw_ibar').str.split('-').list.get(0), 
                spacer=pl.col('raw_ibar').str.split('-').list.get(1)
                )
        )

    if out_parquet is not None:
        final_df.write_parquet(out_parquet)
        logger.io(f'written file {out_parquet}')

    if plot:
        plot_heatmap(zz)
        plot_cormat(cor_mat)
        plot_cluster_qc(clusterer, cor_mat, df, candidate, yyy)

    if umap:
        visualize_clusters(df, zz, zn, clusterer)

    # TODO in some cases it is desirable to remove the 'clone' and unique the results

    return final_df


def plot_heatmap(zz):
    from matplotlib.colors import LogNorm
    with plt.rc_context({'figure.figsize':(30,10), 'figure.dpi':200}):
        sns.heatmap(
                zz
                .sort('sample_id')
                .drop('cbc')
                .to_pandas()
                .set_index(['sample_id']), 
                norm=LogNorm()
                )
    plt.show()

def plot_cormat(cor_mat):
    fig, heat = plt.subplots(1,1, figsize=(5.5, 5))
    map = heat.pcolormesh(cor_mat, cmap='RdYlBu_r')
    fig.colorbar(map)
    plt.show()

def plot_cluster_qc(clusterer, cor_mat, df, candidate, yyy):
    fig, axes = plt.subplots(1,1)
    axes.scatter(clusterer.probabilities_, clusterer.labels_)

    plt.rcParams['figure.dpi'] = 300

    labs = clusterer.labels_
    slabs = np.argsort(labs)

    fig, (heat, lines) = plt.subplots(1,2, figsize=(10, 5))
    heat.pcolormesh(cor_mat[slabs,:][:,slabs], cmap='RdYlBu_r')
    lines.scatter(labs[slabs], range(len(labs)), s=1)
    lines.set_ylim(0, len(labs))
    lines.set_xlim(min(labs)-0.5, max(labs)+0.5)

    plt.show()


    # ibars across samples
    fg = sns.catplot(data=yyy
                .groupby(['sample_id', 'cluster'])
                .agg(pl.col('raw_ibar').n_unique().alias('ibars'))
                .sort('sample_id')
                .to_pandas(), 
                x='cluster', 
                y='ibars', 
                col='sample_id',
                aspect=0.8, 
                kind='bar')

    for ax in fg.axes_dict.values():
        ax.grid()
    plt.show()

    # fraction of sample covered by cluster
    fg = sns.catplot(data=yyy
                .groupby(['sample_id', 'cluster'])
                .agg((pl.col('raw_ibar').n_unique()).alias('ibars'))
                .sort('sample_id')
                .with_columns(pl.col('ibars')/pl.col('ibars').sum().over('sample_id'))
                .to_pandas(), 
                x='cluster', 
                y='ibars', 
                col='sample_id',
                aspect=0.8, 
                kind='bar')

    for ax in fg.axes_dict.values():
        ax.grid()
    fg.set(ylim=(0, 1))
    plt.show()

def visualize_clusters(df, zz, zn, clusterer):
    from scipy.cluster.vq import vq 
    from matplotlib.colors import from_levels_and_colors
    import umap

    reducer = umap.UMAP()

    cm_embedding = reducer.fit_transform(zz.drop(['sample_id', 'cbc']).to_numpy())

    # coloring objects
    cy = plt.rcParams['axes.prop_cycle'].by_key()['color']

    levels= df['sample_id'].unique().to_list()
    cmap, norm = from_levels_and_colors(range(len(levels)), cy[0:len(levels)-1], extend="neither")
    labels_c = {-1:'black', 0:'red', 1:'blue', 2:'orange', 3:'hotpink', 4:'purple'}
    sample_c = dict(zip(df['sample_id'].unique(), range(len(df['sample_id'].unique()))))


    fig, ax = plt.subplots()
    scatter = ax.scatter(
            cm_embedding[:,0], 
            cm_embedding[:,1], 
            s=0.01, cmap=cmap,
            c=[sample_c[i] for i in zz['sample_id']], 
            label=zz['sample_id'].to_list()
    )

    legend1 = ax.legend(
            handles=scatter.legend_elements()[0], 
            labels=levels,
            loc="lower left",
            title="sample", 
            bbox_to_anchor=(1, 0.5))

    ax.add_artist(legend1)


    ib_reducer = umap.UMAP(min_dist=0.45)
    ib_embedding = ib_reducer.fit_transform(zn)

    fig, ax = plt.subplots(figsize=(5,5), dpi=100)
    scatter = ax.scatter(
        ib_embedding[:,0], 
        ib_embedding[:,1], 
        cmap="Paired",
        c=clusterer.labels_)
    ax.legend(*scatter.legend_elements(), loc='center left', bbox_to_anchor=(1, 0.5))

