from logging import warning
from typing import Sequence,Optional
from sys import prefix
import pysam
import regex
import ogtk 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl
from . import ibars as ib
from . import plot as pt


def lala():
    pt.plot_sibling_noise()

def reads_to_molecules(sample_id: str,
           parquet_ifn: str, 
           corr_dict_fn: str | None= None,
           valid_ibars: Sequence | None=None, 
           min_reads: int=1, 
           max_reads: int=int(1e6),
           downsample: int | None=None,
           min_cells_per_ibar: int | None=1000,
           clone: str | None=None,
           columns_ns: str | Sequence='tsos',
           cache_dir: str='/local/users/polivar/src/artnilet/cache',
           use_cache: bool=True
           ) -> pl.DataFrame:
    ''' Analytical routine to process UMIfied shRNA from single-cell assays (parquet format).
        If `clone` is not provided it assumes that can be found as the first element of a dash-splitted `sample_id` string. 
        Off-targets are defined by reads without an ibar pattern match.
        UMIs are filtered based on the number of reads
    '''
    
    # takes some time but is light-weight
    # TODO add cache point
    cache_out = f'{cache_dir}/{sample_id}_r2mols.parquet'
    import os

    if not use_cache or not os.path.exists(cache_out):
        rdf = ib.extract_read_grammar(parquet_ifn = parquet_ifn, batch = sample_id, sample = downsample)
        tot_reads = rdf.shape[0]

        rdf = rdf.with_columns(pl.col('cbc')+"-1")

        # correct cbc if path to dictionary is provided
        if corr_dict_fn is not None:
            print('correcting with dic')
            rdf = ogtk.utils.cd.correct_cbc_pl(rdf, ogtk.utils.cd.load_corr_dict(corr_dict_fn))

        # CPU-intensive
        # TODO add thread control
        rdf = ib.ibar_reads_to_molecules(rdf, modality='single-cell')
        rdf = ib.count_essential_patterns(rdf)

        # guess from the first field of the sample_id the clone of origin
        # if `clone` is not specified
        if clone is None:
            clone = sample_id.split('_')[0]

        ib.noise_spectrum(sample_id, rdf, index = ['dss'], columns=columns_ns)    

        # determine valid_ibars in data
        rdf = ib.filter_ibars(rdf, valid_ibars, min_cells_per_ibar=min_cells_per_ibar)

        ## create mask for wt
        rdf = ib.mask_wt(rdf)

        plot = False
        if plot:
            pt.plot_sibling_noise(rdf)

        # TODO a small number (56 cases in 5M) of molecules ar the 
        # ['cbc', 'umi','raw_ibar'] level map to more than one seq

        #gather stats
        pc_offt = rdf.filter(~pl.col("valid_ibar")).shape[0]/rdf.shape[0]
        tot_umis = rdf.select('umi').n_unique()

        print(f'{pc_offt=:.2%}')
        print(f'{tot_umis=}')
        print(f'{tot_reads=}')

        rdf = (rdf
               .filter(pl.col('raw_ibar').is_not_null())
               .filter(pl.col('umi_dom_reads')>=min_reads)
               .filter(pl.col('umi_dom_reads')<=max_reads)
               .with_columns(pl.col('umi').n_unique().over(['cbc', 'raw_ibar', 'seq', 'wt']).alias('umis_allele'))
                             # ^- this feels too high for my taste??
               .with_columns(pl.col('seq').n_unique().over(['cbc', 'raw_ibar']).alias('n_sib'))
               )
        # consolidates QC metrics into df
        rdf = qc_stats(rdf, sample_id, clone, pc_offt, tot_umis, tot_reads)

        # normalize umi counts based on the size of the cell and levels of expression of the ibar
        rdf = normalize(rdf)
        rdf.write_parquet(cache_out)
    else:
        print(f'loading from {cache_out}')
        rdf = pl.read_parquet(cache_out)
    return(rdf)    

def allele_calling(
        mols: pl.DataFrame, 
        min_umis_allele: int=2,
        by='umis_allele',
        descending=True,
        )-> pl.DataFrame:
    ''' Given a mol-level data frame, returns the top allele for individual ibar-cell data points
    umis_allele is computed going over ['cbc', 'raw_ibar', 'seq', 'wt']  from `reads_to_molecules()`

    The top allele is determined by ``by``, 'umis_allele' by default but could also be:
        - 'cs_norm_umis_allele'
        - 'ib_norm_umis_allele'
        - 'db_norm_umis_allele'  
    '''
    # major collapse event where the top ranking sequence as the final allele is selected.

    df = (
        mols
        #.filter(pl.col('wt').is_not())
        .sort(by, descending=descending) # we give priority to non-wt
        #.sort('db_norm_umis_allele', descending=False) # we give priority to non-wt
        #.sort(['umis_allele', 'wt'], [True, False]) # we give priority to non-wt
        #.filter(pl.col('umis_allele')>=2)
        .groupby(['cbc', 'raw_ibar'], maintain_order=True)
        # https://github.com/pola-rs/polars/issues/10054
        .head(1) # <- this is it! # change to .first() or .top_k()
        .select(['cbc', 'raw_ibar', 'seq', 'wt', 'umi_reads',\
                'umis_allele', 'umis_cell', 'umis_ibar', 'cs_norm_umis_allele', 'ib_norm_umis_allele', 'db_norm_umis_allele',\
                'cluster', 'sample_id', 'valid_ibar'])
        .with_columns(pl.col('raw_ibar').n_unique().over('cbc').alias('cov'))
    )

    return(df)
    df =  (
            df
            #.filter(pl.col('wt').is_not())
            .sort('umis_allele', True) # we give priority to non-wt
            #.sort(['umis_allele', 'wt'], [True, False]) # we give priority to non-wt
            .filter(pl.col('umis_allele')>=min_umis_allele)
            .groupby(['cbc', 'raw_ibar'], maintain_order=True)
            .head(1) # <- this is it!
    )

def to_matlin(df, expr: None | pl.Expr, sort_by='cluster', cells=100):
    ''' Returns a matlin object from an allele-called (e.g. `allele_calling()`)

        If `expr` is `None` then the `top` cells are selected. `top` cells are defined by `cells` multiplied by the total number if ibars. 
    '''
    matl = ogtk.ltr.matlin()
    #subset = ['kalhor_id', 'nspeed', 'raw_ibar']
    tot_ibars = df['raw_ibar'].n_unique()
    top = cells * tot_ibars
    top = cells 

    if expr is None:
        matl.ingest_ibars_pl(df.sample(top), sort_by=sort_by)
    else:
        matl.ingest_ibars_pl(df.filter(expr), sort_by=sort_by)

    return matl


def to_compare_top_alleles(df):
    ''' Returns a data frame with the number of umis supporting wt and non-wt states of a given ibar-cell
    Fields:
    true = wt
    false = non-wt
    lt = log wt
    lf = log non-wt
    '''
    df = (df
       .filter(pl.col('valid_cell'))
       .filter(pl.col('valid_ibar'))
       #.select(['cbc', 'raw_ibar', 'wt', 'umis_allele'])
       .pivot(index=['cbc', 'raw_ibar'], columns='wt', values='umis_allele')
    )
    df = (df
          .fill_null(0)
          .with_columns(
              [(1+pl.col('true')).log10().alias('lt'), 
               (1+pl.col('false')).log10().alias('lf')]
          ))
    return(df)

def normalize(df):
    ''' normalize umi counts based on the size of the cell and levels of expression of the ibar
    '''
    return(df
        .with_columns(pl.col('umi').n_unique().over('cbc').alias('umis_cell'))
        .with_columns(pl.col('umi').n_unique().over('raw_ibar').alias('umis_ibar'))
        .with_columns((pl.col('umis_allele')/pl.col('umis_cell')).prefix('cs_norm_'))
        .with_columns((pl.col('umis_allele')/pl.col('umis_ibar')).prefix('ib_norm_'))
        .with_columns((pl.col('umis_allele')/pl.col('umis_cell')/pl.col('umis_ibar')).prefix('db_norm_'))
    )

def qc_stats(df, sample_id, clone, pc_offt, tot_umis, tot_reads):
    ''' helper function that consolidates QC metrics into df
    '''
    return(
            df
           .with_columns(pl.lit(sample_id).alias('sample_id'))
           .with_columns(pl.lit(clone).alias('clone'))
           .with_columns(pl.lit(pc_offt).alias('qc_pc_offt'))
           .with_columns(pl.lit(tot_umis).alias('qc_tot_umis'))
           .with_columns(pl.lit(tot_reads).alias('qc_tot_reads'))
           )


def basure():
    rich.print(f"%wt {data['wt'].mean() *100:.2f}")
    rich.print(f"%wt {data.filter(pl.col('umis_seq')>0)['wt'].mean()*100 :.2f} F")

    if plot:
        plt.figure()
        fg = sns.displot(data=data.to_pandas(), x='umis_seq', aspect=2, log_scale=10)
        fg.ax.axvline(1.5, c='r')
        fg.ax.set_title(batch)
        fg.ax.set_yscale('log')
        fg.ax.grid()

        plt.figure()
        fg = sns.displot(data=data.to_pandas(), x='umis_seq', aspect=2, log_scale=10)
        fg.ax.axvline(1.5, c='r')
        fg.ax.set_title(batch)
        fg.ax.grid()

    # actually filter umis seen once
    if not unfiltered:
        data = data.filter(pl.col('umis_seq')>1)
    if plot:
        plt.figure()
        fg = sns.displot(
            data=data.groupby('raw_ibar').agg(pl.col('cbc').n_unique()).to_pandas(),
            x = 'cbc', 
            aspect = 1.5,
            kind='kde', 
            log_scale=10,
            )
        fg.ax.set_title(batch)
        fg.ax.grid()


        fg = sns.displot(data=data.groupby('cbc').count().to_pandas(), x='count', kind='ecdf', aspect=1.0)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} ibars per cell')

        fg = sns.displot(data=data.groupby('raw_ibar').count().to_pandas(), x='count', kind='ecdf', aspect=1.0)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} cells per ibar')

        fg = sns.catplot(data=data.groupby('cbc').count().to_pandas(), y='count', kind='boxen', aspect=0.5)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} ibars per cell')
        fg.ax.set_ylim([0, fg.ax.get_ylim()[1]])

        fg = sns.catplot(data=data.groupby('raw_ibar').count().to_pandas(), y='count', kind='boxen', aspect=0.5)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} cells per ibar')
        fg.ax.set_ylim([0, fg.ax.get_ylim()[1]])

        fg = sns.catplot(data=data.groupby('cbc').count().with_columns(pl.col('count')/len(valid_ibars)).to_pandas(), y='count', kind='boxen', aspect=0.5)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} ibars per cell')
        fg.ax.set_ylim([0, fg.ax.get_ylim()[1]])

        fg = sns.catplot(data=data.groupby('raw_ibar').count().with_columns(pl.col('count')/total_cells).to_pandas(), y='count', kind='boxen', aspect=0.5)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} cells per ibar')
        fg.ax.set_ylim([0, fg.ax.get_ylim()[1]])

    #rich.print(f"median ibars·per·cell {data.groupby('cbc').count().with_columns(pl.col('count'))['count'].median()}")
    #rich.print(f"fraction of ibars covered {data.groupby('cbc').count().with_columns(pl.col('count')/len(valid_ibars))['count'].median()*100:.2f}%")
    #rich.print(f"median cells·per·ibar {data.groupby('raw_ibar').count().with_columns(pl.col('count'))['count'].median()}")
    #rich.print(f"fraction of cells covered {data.groupby('raw_ibar').count().with_columns(pl.col('count')/total_cells)['count'].median()*100:.2f}%")

    
    if False:
        mat = data.pivot(columns='raw_ibar', index='cbc', values='umis_seq')
        dd = dict(zip(*data.select(['raw_ibar', 'nspeed']).unique()))
        speed_sort = np.argsort([dd[i] for i  in mat.columns if i in dd.keys()])
        np.array([dd[i] for i  in mat.columns if i in dd.keys()])[speed_sort]
        plt.figure()
        mapp = plt.pcolormesh(mat.drop('cbc').select(pl.all().log10()).to_numpy()[:,speed_sort])
        plt.colorbar(mapp)

    if do_plot:
        plt.figure()
        fg = sns.catplot(
                data=data
                        .groupby(['raw_ibar', 'speed'])
                        .agg(pl.col('cbc').count())
                        .to_pandas(),
                y='cbc', 
                x='speed', 
                kind='boxen')
        fg.ax.set_ylim((0, total_cells))
        fg.ax.set_title(batch)
        fg.ax.grid()
        
        plt.figure()
        fg = sns.catplot(
                data=data
                        .groupby(['raw_ibar', 'speed'])
                        .agg(pl.col('cbc').count()/total_cells)
                        .to_pandas(),
                y='cbc', 
                x='speed', 
                kind='boxen')
        fg.ax.set_ylim((0,1))
        fg.ax.set_title(batch)
        fg.ax.grid()

        plt.figure()
        fg = sns.catplot(
                data=data
                        .groupby(['raw_ibar', 'speed'])
                        .agg(pl.col('cbc').count()/total_cells)
                        .to_pandas(),
                y='cbc', 
                x='speed', 
                alpha=0.8)
        fg.ax.set_ylim((0,1))
        fg.ax.set_title(batch)
        fg.ax.grid()


        fg = sns.catplot(
                data=data
                        .groupby(['raw_ibar', 'speed', 'wt'])
                        .agg([pl.col('cbc').count(), pl.col('umis_seq').mean()])
                        .to_pandas(),
                y='umis_seq', 
                x='speed', 
                hue='wt',
                kind='boxen')
        fg.ax.grid()
        fg.ax.set_title(batch)

    #fg.ax.set_ylim((0,1))
        # open questions
        # are we accounting for the G+?
        # - no: [···TSO···]G[···WT···]
        # - no: [···TSO···]GTGTAACTTAACACTGAGTG[···CNSCFL···] is non WT when it should
        # - probably this doesn't matter since we are going for the top seq, unless these cases are the dominant allele
        # analysis at the integration level (lineage tracing)
        # plot the fraction of the pool that a to top allele shows, instead of the raw umis_seq
        # annotate with wl
        # cbc and umi corrections
        # single cell stats (ibars recovered)
        # cells per ibar as ecdf
        # cells per ibar as ecdf/boxen but normalized
        # keep static list of ibars


    rich.print(':vampire:')
    data = (data
       .with_columns(pl.lit(batch).alias('sample_id'))
       .with_columns(pl.col('speed').cast(pl.Utf8))
       .with_columns(pl.col('raw_ibar').cast(pl.Utf8))
       .with_columns(pl.col('cbc').cast(pl.Utf8))
       .with_columns(pl.col('sample_id').cast(pl.Utf8))
       .with_columns(pl.col('diversity').cast(pl.Utf8))
        )
    return(data)

def compute_clonal_composition(
       df: pl.DataFrame, 
       clone_dict: dict|None = None,
       normalize_cluster_size=False,
       )->pl.DataFrame:
    ''' Provided a chimera (cells from more than one cell line) ibar mols data frame, assign
        (per-cell) the normalized|total molecule counts mapping to cell line-specific
        integrations. 

        df: is an ibar mol-level pl data frame (usually exp.mols)
        the normalization corresponds to the size of the ibar cluster
    '''

    
    #normalize counts per cluster size - not implemented
    #cell_size = 'ib_norm_umis_allele' if normalize_cluster_size else 'umis_cell'
    #values = "ncounts" if normalize_cluster_size else 'counts'
    values = 'counts'

    if clone_dict is not None:
        df = df.with_columns(pl.col('cluster').map_dict(clone_dict))

    index = ['cbc', 'umis_cell']
    per_cell_clonal_composition = (
        df
         .filter(pl.col('cluster').is_not_null())
         .with_columns(pl.col('raw_ibar').n_unique().over(['cluster']).alias('cluster_size'))
         .groupby(['cbc', 'cluster_size', 'umis_cell'])
            .agg(pl.col('cluster').value_counts()).explode('cluster').unnest('cluster').rename({"cluster":"clone_score"})
         #normalize counts per cluster size - not implemented
         #.with_columns((pl.col('counts')/pl.col('cluster_size')).prefix('n'))
         #.pivot(index=['cbc', 'umis_cell'], columns='clone_score', values=values, aggregate_function='sum')
         .pivot(index=index, columns='clone_score', values=values)
         .sort('cbc')   
         .fill_nan(0)
         .fill_null(0)
        )
    # Add ibar_cluster prefix to the pivoted columns
    prefix = "ibc_"
    per_cell_clonal_composition.columns = [prefix + col if col not in index else col 
                                           for col in per_cell_clonal_composition.columns]

    return(per_cell_clonal_composition)
