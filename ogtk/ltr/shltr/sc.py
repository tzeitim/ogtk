from typing import Sequence,Optional
import os
import ogtk 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl
import rich as rich
from . import ibars as ib
from . import plot as pt

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()


@ogtk.utils.log.call
def reads_to_molecules(sample_id: str,
           parquet_ifn: str, 
           corr_dict_fn: str | None= None,
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
    
    spacers = ib.load_wl(True)['spacer']

    cache_out = f'{cache_dir}/{sample_id}_r2mols.parquet'
    encoded_reads_parquet = f'{cache_dir}/{sample_id}_encoded_reads.parquet'

    if use_cache and os.path.exists(cache_out):
        logger.info(f'loading from cache {cache_out}')
        rdf = pl.read_parquet(cache_out)
    else:
        # generate a reads data frame 
        if os.path.exists(encoded_reads_parquet) and use_cache:
            logger.info(f"loading {encoded_reads_parquet=}")
            rdf = pl.read_parquet(encoded_reads_parquet)
        else:
            rdf = ib.encode_reads_with_dummies(parquet_ifn = parquet_ifn, batch = sample_id, sample = downsample)
            rdf.write_parquet(encoded_reads_parquet)

        tot_reads = rdf.shape[0]
        rdf = rdf.with_columns(pl.col('cbc')+"-1")

        # correct cbc if path to dictionary is provided
        if corr_dict_fn is not None:
            logger.info('correcting with dic')
            rdf = ogtk.utils.cd.correct_cbc_pl(rdf, ogtk.utils.cd.load_corr_dict(corr_dict_fn))

        # CPU-intensive
        # TODO add thread control
        
        rdf = ib.ibar_reads_to_molecules(rdf, modality='single-cell')

        rdf = ib.extract_spacer(rdf, spacers) #type: ignore
        
        rdf = ib.encode_wt(rdf, spacers) #type: ignore
        
        rdf = ib.count_essential_patterns(rdf)

        # guess from the first field of the sample_id the clone of origin
        # if `clone` is not specified
        if clone is None:
            clone = sample_id.split('_')[0]

        ib.noise_spectrum(sample_id, rdf, index = ['dss'], columns=columns_ns)    

        # determine valid_ibars in data
        #rdf = ib.mask_valid_ibars(rdf, valid_ibars, min_cells_per_ibar=min_cells_per_ibar)

        ## create mask for wt
        rdf = ib.mask_wt(rdf)

        plot = False
        if plot:
            pt.plot_sibling_noise(rdf)

        # TODO a small number (56 cases in 5M) of molecules ar the 
        # ['cbc', 'umi','raw_ibar'] level map to more than one seq

        #gather stats
        # how to determine pc_offt TODO
        #pc_offt = rdf.filter(~pl.col("valid_ibar")).shape[0]/rdf.shape[0]
        tot_umis = rdf.select('umi').n_unique()

        #logger.info(f'{pc_offt=:.2%}')
        logger.info(f'{tot_umis=}')
        logger.info(f'{tot_reads=}')

        rdf = (rdf
               .filter(pl.col('raw_ibar').is_not_null())
               .filter(pl.col('umi_dom_reads')>=min_reads)
               .filter(pl.col('umi_dom_reads')<=max_reads)
               .with_columns(pl.col('umi').n_unique().over(['cbc', 'raw_ibar', 'seq', 'wt']).alias('umis_allele'))
                             # ^- this feels too high for my taste??
               .with_columns(pl.col('seq').n_unique().over(['cbc', 'raw_ibar']).alias('n_sib'))
               )
        # consolidates QC metrics into df
        rdf = qc_stats(rdf, sample_id, clone, tot_umis, tot_reads)

        # normalize umi counts based on the size of the cell and levels of expression of the ibar
        #rdf = normalize(rdf)

        rdf.write_parquet(cache_out)
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
    alleles = (mols
        .lazy()
        #.filter(pl.col('doublet_prediction')=='singlet')
        .with_columns(pl.count().over('cbc', 'raw_ibar', 'seq').alias('umis_allele'))
        .group_by([ 'doublet_prediction', 'cbc', 'raw_ibar',])
        .agg(
            ties = ((pl.col('umis_allele') == pl.col('umis_allele').max()).sum()/pl.col('umis_allele').max())-1,
        
            umis_per_ibar = pl.count(),
            umis_top_allele = pl.col('umis_allele').max(),
            siblings = pl.col('seq').n_unique(),
            
            norm_counts = pl.col('cat_db_norm_umis_allele').min(),
            norm_counts_acc_raw = pl.col('cat_db_norm_umis_allele').gather(pl.col('umis_allele').arg_max()),
            
            raw_counts = pl.col('umis_allele').max(),
            
            top_allele_raw = pl.col('seq').gather(pl.col('umis_allele').arg_max()),
            top_allele_norm = pl.col('seq').gather(pl.col('cat_db_norm_umis_allele').arg_min()),
            #seq_raw = pl.col('seq').gather(pl.col('umis_allele').arg_max()),
                          
            umi_dom_reads = pl.col('umi_dom_reads').max(),
            umi_reads = pl.col('umi_reads').max(),

        )
        .with_columns(pl.col('raw_ibar').n_unique().over('cbc').alias('cells_per_ibar'))
        .explode('top_allele_raw')
        .explode('top_allele_norm')
        .explode('norm_counts_acc_raw')
        #.explode('seq_raw')
        .collect()
    )

    return(alleles)

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

@ogtk.utils.log.call
def normalize(df:pl.DataFrame, over_cells:List=['cbc'], over_ibars:List=['raw_ibar'], prefix:str='', expr_cells:None|pl.Expr=None):
    ''' normalize umi counts based on the size of the cell and levels of expression of the ibar
    '''
    if expr_cells is None:
        expr_cells = pl.lit(True)
    return(df
        .with_columns(pl.count().over(over_cells).alias(prefix+'umis_cell'))
       # .with_columns(pl.count().over(over_ibars).alias(prefix+'umis_ibar'))
           .with_columns(pl.when(
                            expr_cells is None
                         )
                         .then(
                            pl.count().over(over_ibars)
                         )
                         .otherwise(
                            pl.col('raw_ibar').count()/pl.col('cbc').filter(expr_cells).count()                             
                         )
                        .alias(prefix+'umis_ibar'))
        .with_columns((pl.col('umis_allele')/pl.col(prefix+'umis_cell')).name.prefix(prefix+'cs_norm_'))
        .with_columns((pl.col('umis_allele')/pl.col(prefix+'umis_ibar')).name.prefix(prefix+'ib_norm_'))
        .with_columns((pl.col('umis_allele')/(pl.col(prefix+'umis_cell')*pl.col(prefix+'umis_ibar'))).name.prefix(prefix+'db_norm_'))
           .with_columns(pl.col(prefix+'cs_norm_umis_allele').log10() * -1 )
           .with_columns(pl.col(prefix+'ib_norm_umis_allele').log10() * -1 )
           .with_columns(pl.col(prefix+'db_norm_umis_allele').log10() * -1 )
    )
@ogtk.utils.log.call
def qc_stats(df, sample_id, clone, tot_umis, tot_reads):
    ''' helper function that consolidates QC metrics into df
    '''
    return(
            df
           .with_columns(pl.lit(sample_id).alias('sample_id'))
           .with_columns(pl.lit(clone).alias('clone'))
           #.with_columns(pl.lit(pc_offt).alias('qc_pc_offt'))
           .with_columns(pl.lit(tot_umis).alias('qc_tot_umis'))
           .with_columns(pl.lit(tot_reads).alias('qc_tot_reads'))
           )

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
    #values = "ncounts" if normalize_cluster_size else 'count'
    values = 'count'

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
         #.with_columns((pl.col('count')/pl.col('cluster_size')).prefix('n'))
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
