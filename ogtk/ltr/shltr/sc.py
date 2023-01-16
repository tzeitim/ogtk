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
           valid_ibars: Sequence | None=None, 
           min_reads: int=1, 
           max_reads: int=int(1e6),
           down_sample: int | None =None,
           min_cells_per_ibar: int | None=1000,
           clone: str | None=None,
           ) -> pl.DataFrame:
    ''' Analytical routine to process UMIfied shRNA from bulk assays.
        If `clone` is not provided it assumes that can be found as the first element of a dash-splitted `sample_id` string. 
        Off-targets are defined by reads without an ibar pattern match.
        UMIs are filtered based on the number of reads
    '''
    
    # takes some time but is light-weight
    # TODO add cache point
    rdf = ib.extract_read_grammar_new(parquet_ifn = parquet_ifn, batch = sample_id, sample=down_sample)
    tot_reads = rdf.shape[0]

    # CPU-intensive
    # TODO add thread control
    rdf = ib.ibar_reads_to_molecules(rdf, modality='single-cell')
    rdf = ib.count_essential_patterns(rdf)

    # guess from the first field of the sample_id the clone of origin
    # if `clone` is not specified
    if clone is None:
        clone = sample_id.split('_')[0]

    ib.noise_spectrum(sample_id, rdf, index = ['dss'], columns='tsos')    

    # determine valid_ibars in data
    rdf = ib.filter_ibars(rdf, valid_ibars, min_cells_per_ibar=min_cells_per_ibar)

    ## create mask for wt
    rdf = ib.mask_wt(rdf)

    plot = False
    if plot:
        pt.plot_sibling_noise(rdf)

    #gather stats
    pc_offt = rdf.filter(~pl.col("valid_ibar")).shape[0]/rdf.shape[0]
    tot_umis = rdf.select('umi').n_unique()

    print(f'{pc_offt=:.2%}')
    print(f'{tot_umis=}')
    print(f'{tot_reads=}')

    return(rdf
           .filter(pl.col('raw_ibar').is_not_null())
#           .filter(pl.col('raw_ibar').is_in(valid_ibars))
           .filter(pl.col('umi_dom_reads')>=min_reads)
           .filter(pl.col('umi_dom_reads')<=max_reads)
           .with_columns(pl.lit(sample_id).alias('sample_id'))
           .with_columns(pl.lit(clone).alias('clone'))
           .with_columns(pl.lit(pc_offt).alias('qc_pc_offt'))
           .with_columns(pl.lit(tot_umis).alias('qc_tot_umis'))
           .with_columns(pl.lit(tot_reads).alias('qc_tot_reads'))
           .with_column(pl.col('umi').n_unique().over(['cbc', 'raw_ibar', 'seq', 'wt']).alias('umis_allele'))
                         # ^- this feels too high for my taste??
            .with_columns(pl.col('seq').n_unique().over(['cbc', 'raw_ibar']).alias('n_sib'))
           # normalize umi counts based on the unfiltered size of the cell
           .with_column(pl.col('umi').n_unique().over('cbc').alias('umis_cell'))
           .with_column((pl.col('umis_allele')/pl.col('umis_cell')).prefix('norm_'))
          )

def allele_calling(
        df: pl.DataFrame, 
        min_umis_allele: int=2
        )-> pl.DataFrame:
    ''' Given a mol-level data frame, returns the top allele for individual ibar-cell data points
    umis_allele is computed going over ['cbc', 'raw_ibar', 'seq', 'wt']  from `reads_to_molecules()`
    '''
    # major collapse event where the top ranking sequence as the final allele is selected.
    df =  (
            df
            #.filter(pl.col('wt').is_not())
            .sort('umis_allele', True) # we give priority to non-wt
            #.sort(['umis_allele', 'wt'], [True, False]) # we give priority to non-wt
            .filter(pl.col('umis_allele')>=min_umis_allele)
            .groupby(['cbc', 'raw_ibar'], maintain_order=True)
            .head(1) # <- this is it!
    )
    return(df)

def to_matlin(df, subset='cluster', cells=100):
    ''' Returns a matlin object from an allele-called (e.g. `allele_calling()`)
    '''
    matl = ogtk.ltr.matlin()
    #subset = ['kalhor_id', 'nspeed', 'raw_ibar']
    matl.ingest_ibars_pl(df.head(cells* df['raw_ibar'].n_unique()), subset=subset)

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

        fg = sns.catplot(data=data.groupby('cbc').count().with_column(pl.col('count')/len(valid_ibars)).to_pandas(), y='count', kind='boxen', aspect=0.5)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} ibars per cell')
        fg.ax.set_ylim([0, fg.ax.get_ylim()[1]])

        fg = sns.catplot(data=data.groupby('raw_ibar').count().with_column(pl.col('count')/total_cells).to_pandas(), y='count', kind='boxen', aspect=0.5)
        fg.ax.grid()
        fg.ax.set_title(f'{batch} cells per ibar')
        fg.ax.set_ylim([0, fg.ax.get_ylim()[1]])

    #rich.print(f"median ibars·per·cell {data.groupby('cbc').count().with_column(pl.col('count'))['count'].median()}")
    #rich.print(f"fraction of ibars covered {data.groupby('cbc').count().with_column(pl.col('count')/len(valid_ibars))['count'].median()*100:.2f}%")
    #rich.print(f"median cells·per·ibar {data.groupby('raw_ibar').count().with_column(pl.col('count'))['count'].median()}")
    #rich.print(f"fraction of cells covered {data.groupby('raw_ibar').count().with_column(pl.col('count')/total_cells)['count'].median()*100:.2f}%")

    
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
       .with_column(pl.lit(batch).alias('sample_id'))
       .with_column(pl.col('speed').cast(pl.Utf8))
       .with_column(pl.col('raw_ibar').cast(pl.Utf8))
       .with_column(pl.col('cbc').cast(pl.Utf8))
       .with_column(pl.col('sample_id').cast(pl.Utf8))
       .with_column(pl.col('diversity').cast(pl.Utf8))
        )
    return(data)


