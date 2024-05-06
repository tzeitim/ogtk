from typing import Sequence,Optional, List
import os
import ogtk 
import polars as pl
#import polars.selectors as cs
import rich as rich

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()


# function that takes a mol-level data frame and returns a cell-level data frame
# with the top allele for each cell
@ogtk.utils.log.call
def allele_calling(
        mols: pl.DataFrame, 
        excluded_cells_field:str='doublet_prediction',
        )-> pl.DataFrame:
    ''' Given a mol-level data frame, returns the top allele for individual ibar-cell data points
    umis_allele is computed going over ['cbc', 'raw_ibar', 'seq']  from `reads_to_molecules()`

    The top allele is determined by ``by``, 'umis_allele' by default but could also be:
        - 'cs_norm_umis_allele'
        - 'ib_norm_umis_allele'
        - 'db_norm_umis_allele'  
    '''
    if excluded_cells_field not in mols.columns:
        raise ValueError(f'{excluded_cells_field} not in {mols.columns}. Pleae provide a valid field name such as "doublet_prediction" or "excluded_cell"')

    alleles = (mols
        .lazy()
        #.filter(pl.col('doublet_prediction')=='singlet')
        .with_columns(pl.count().over('cbc', 'raw_ibar', 'seq').alias('umis_allele'))
               .group_by([excluded_cells_field, 'cbc', 'raw_ibar']) #type: ignore
        .agg(
            ties = (
                (pl.col('umis_allele') == pl.col('umis_allele').max()).sum()
                /pl.col('umis_allele').max()
                )-1,
        
            umis_per_ibar = pl.count(),
            umis_top_allele = pl.col('umis_allele').max(),
            siblings = pl.col('seq').n_unique(),
            
            norm_counts = pl.col('norm_counts').max(),
            norm_counts_acc_raw = pl.col('norm_counts').gather(pl.col('umis_allele').arg_max()), #type: ignore
            
            raw_counts = pl.col('umis_allele').max(),
            
            top_allele_raw = pl.col('seq').gather(pl.col('umis_allele').arg_max()), #type: ignore
            top_allele_norm = pl.col('seq').gather(pl.col('norm_counts').arg_max()), #type: ignore
                          
            umi_dom_reads = pl.col('umi_dom_reads').max(),
            umi_reads = pl.col('umi_reads').max(),
        )
        .explode('top_allele_raw')
        .explode('top_allele_norm')
        .explode('norm_counts_acc_raw')
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
def normalize_mol_counts(
        df:pl.DataFrame,
        over_cells:Sequence[str]=['cbc'],
        over_ibars:Sequence[str]=['raw_ibar'],
        bg_cells_expr:None|pl.Expr=None,
        )->pl.DataFrame:
    ''' normalize umi counts based on the size of the cell and levels of expression of the ibar
        cells can be grouped using `over_cells` and ibars using `over_ibars`
        `bg_cells_expr` is a filter expression to select cells to be used for normalization
        if `bg_cells_expr` is `None` then groupby `over_ibars` is used to compute the number of cells
        
        The normalized counts are -log10 transformed
    '''

    if bg_cells_expr is None:
        bg_cells_expr = pl.lit(True)

    model_from = pl.col('doublet_prediction')=="singlet"
    return (df
            .filter(pl.col('doublet_prediction')=='singlet')
            .with_columns(pl.col('umi').count().over(['cbc', 'doublet_prediction']).alias('umis_per_cell'))
            .with_columns(pl.col('umi').count().over(['doublet_prediction']).alias('total_umis'))
            .with_columns(pl.col('umi').count().over(['raw_ibar', 'doublet_prediction']).alias('total_umis_ibar'))
            .with_columns(pl.col('umi').filter(bg_cells_expr).count().over('raw_ibar').alias('umis_per_ibar_model'))
            .with_columns(pl.col('umi').count().over(['cbc', 'raw_ibar', 'seq']).alias('umis_allele'))
            .with_columns(p_cell = pl.col('umis_per_cell')/pl.col('total_umis'))
            .with_columns(p_ibar = pl.col('total_umis_ibar')/pl.col('total_umis'))
            .with_columns(p_join = pl.col('p_cell') * pl.col('p_ibar'))
            .with_columns(weighted_counts =pl.col('umis_allele')*pl.col('p_join'))
           )


@ogtk.utils.log.call
def qc_stats(df, sample_id, clone, tot_umis, tot_reads):
    ''' helper function that consolidates QC metrics into df
    '''
    return(
            df
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
