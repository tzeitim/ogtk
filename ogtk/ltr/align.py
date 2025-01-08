import polars as pl
import pyseq_align
from pyseq_align import NeedlemanWunsch, SmithWaterman
from functools import lru_cache
import ngs_tools

from typing import *
from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()

@lru_cache
def lambda_needlemanw(seq, foc_ref, aligner):
    ''' 
    '''
    alignment = aligner.align(foc_ref, seq)

    return {
        'cigar_str': ngs_tools.sequence.alignment_to_cigar(alignment.result_a, alignment.result_b),
        'aligned_ref':  alignment.result_a,
        'aligned_seq':  alignment.result_b,
        'alignment_score':  alignment.score,
        'start_ref': alignment.pos_a,
        'start_seq': alignment.pos_b,
        }

@lru_cache
def lambda_smithwaterman(seq, foc_ref, aligner):
    ''' 
    '''
    alignment = aligner.align(foc_ref, seq, 1)[0]

    return {
        'cigar_str': ngs_tools.sequence.alignment_to_cigar(alignment.result_a, alignment.result_b),
        'aligned_ref':  alignment.result_a,
        'aligned_seq':  alignment.result_b,
        'alignment_score':  alignment.score,
        'start_ref': alignment.pos_a,
        'start_seq': alignment.pos_b,
        }

def default_aligner():
    return NeedlemanWunsch(
            gap_open=-20,
            gap_extend=-1,
            no_start_gap_penalty=True,
        )

def return_aligned_alleles(
        df:pl.DataFrame,
        ref_str:str,
        seq_trim_field:str,
        aligner:Union["pyseq_align.NeedlemanWunsch", "pyseq_align.SmithWaterman"], #pyright: ignore
        intbc_field:str='intbc',
        min_group_size:int=100,
        )->pl.DataFrame:
    '''
    Requires a data frame that:
    - belongs to a single integration 
    - sequences have been trimmed to their minimally informative interval

    E.g.
    df_alg = (
            df
            .sample(int(1e5))
            .group_by('intbc')
            .map_groups(al.return_aligned_alleles)
    )

    Returns the original data frame with additional fields that contain a an aligned version of the original sequence
     = position-specific code:
         - -1 is for gap
         - 0 match
         - >0 insertion length

    Alignment procedure:
     - Trim the original sequence based on fuzzy 5- and 3-prime anchor
       sequences. In the case of RNA-based detection the TSO fulfills this
       role. For DNA-based (zombie) it's the 3' end of the U6 promoter 
     - 
     # TODO watch out for stammering on alignments - how do they look? search for XXX in the seq
     # TODO keep a record of sequences that can't be aligned (perhaps upstream this function)

    '''

    oschema = {"aseq":str, "alg":pl.List(pl.Int64), "aweight":pl.Int64}
    merged_schema = df.schema.copy()
    merged_schema.update(oschema)

    query_schema = {seq_trim_field:pl.Utf8}#, 'oseq':pl.Utf8}

    # alignment.schema=OrderedDict([('oseq', String), ('cigar', String), ('alg_ref', String), ('alg_seq', String), ('alg_score', Int64)])

    alignment_schema = {
           'cigar_str':pl.Utf8,
           'aligned_ref': pl.Utf8, 
           'aligned_seq': pl.Utf8, 
           'alignment_score': pl.Int64,
           'start_ref': pl.Int64,
           'start_seq': pl.Int64,
           }
    
    merged_schema = query_schema.copy()
    merged_schema.update(alignment_schema)

    #assert df['sample_id'].n_unique() == 1, \
    #        "Please provide a data frame pre-filtered with a single sample_id. \
    #        Best when this function is called within a group_by"

    assert df[intbc_field].n_unique() == 1, \
            "Please provide a data frame pre-filtered with a single integration. \
            Best when this function is called within a group_by"

    #foc_sample = df['sample_id'][0]
    foc_int = df[intbc_field][0]

    if df.shape[0]<min_group_size:
        logger.debug(f"No entries found in group/data frame ")
        return pl.DataFrame(schema=merged_schema)

    ### delete?
    if len(df)==0:
        logger.debug(f"No entries found to generate a fasta file {foc_int} ")
        return pl.DataFrame(schema=merged_schema)


    # TODO: make this optional
    query = df.select(seq_trim_field).unique()

    # perform the alignment
    alignment_udf = lambda_needlemanw if isinstance(aligner, NeedlemanWunsch) else lambda_smithwaterman

    query =(
            query
            .with_columns(
                alignment= pl.col(seq_trim_field)
                   .map_elements(
                       lambda x: alignment_udf(x, ref_str, aligner), 
                                return_dtype=pl.Struct(alignment_schema)
                    )
            )
            .select(seq_trim_field, 'alignment')
            )

    alignment = (
            query.with_columns(
                cigar_str=pl.col('alignment').struct.field('cigar_str'), 
                aligned_ref=pl.col('alignment').struct.field('aligned_ref'), 
                aligned_seq=pl.col('alignment').struct.field('aligned_seq'), 
                alignment_score=pl.col('alignment').struct.field('alignment_score'), 
                start_ref=pl.col('alignment').struct.field('start_ref'),
                start_seq=pl.col('alignment').struct.field('start_seq'),
            )
            .drop('alignment')
        )
    query = None

    alignment = alignment.select(list(merged_schema.keys()))
    return alignment

def compute_coverage(df, ref_str, seq_field='r2_seq', max_range:int|None=None)-> pl.DataFrame:
    dfa = return_aligned_alleles(df, ref_str=ref_str, seq_trim_field=seq_field, aligner=SmithWaterman())

    result = (dfa
        .with_columns(
            pl.struct(['start_ref', seq_field])
            .map_elements(lambda x: list(range(x['start_ref'], x['start_ref'] + len(x[seq_field]))), return_dtype=pl.List(pl.Int64))
            .alias('covered_positions')
        )
        .explode('covered_positions')
        .group_by('covered_positions')
        .len()
        .sort('covered_positions')
    )

    # add missing positions with zero counts by lending a 1
    if max_range is None:
        max_range = result['covered_positions'].max()

    dummy_positions = pl.DataFrame([(i, 1) for i in range(0, max_range+1)], orient='row', schema=result.schema)

    result = (
        pl.concat( [result, dummy_positions])
        .with_columns(pl.col('len'))
        .group_by('covered_positions')
        .agg(pl.col('len').sum()-1) # we get back the 1 we lend
        .sort('covered_positions')
        .rename({'len':'coverage'})
     )
    return result
