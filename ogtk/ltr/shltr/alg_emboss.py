from typing import *
from colorhash import ColorHash
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
import os
import regex
import seaborn as sns
import tempfile
from pyseq_align import NeedlemanWunsch#, SmithWaterman
import ngs_tools
from functools import lru_cache

import ogtk
import ogtk.utils as ut
import ogtk.utils.log as log

from ogtk.ltr import ltr_utils
from . import plot as pt

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()

import alg_lib as alg

from ogtk.ltr.shltr.ibars import return_feature_string, fuzzy_match_str

def return_aligned_alleles_emboss(
    df,
    lim=50,
    correct_tss_coord: None|int = 18,
    keep_intermediate=True,
    min_group_size=100,
    verbose=True,
    )->pl.DataFrame:
    '''
    Requires a data frame that:
    - belongs to a single ibar
    - belongs to a single sample
    - is annotated with a 'kalhor_id' 

    E.g.
    df_alg = (
            df
            .sample(int(1e5))
            .group_by('raw_ibar','sample_id')
            .map_groups(ogtk.shltr.ibars.return_aligned_alleles)
    )

    Returns the original data frame with additional fields that contain a an aligned version of the original sequence
     - position-specific code:
         - -1 is for gap
         - 0 match
         - >0 insertion length
    
    Alignment fields: "alg_ref", "alg_seq", "alg_score"

    '''
    ref_str = "TTTCTTATATGGG[SPACER]GGGTTAGAGCTAGA[IBAR]AATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGGTACT"
    schema = {"sample_id":str, "ibar":str, "id":str, "aseq":str, "sweight":pl.Int64, "alg":pl.List(pl.Int64)}
    oschema = {"aseq":str, "alg":pl.List(pl.Int64), "aweight":pl.Int64}
    merged_schema = df.schema.copy()
    merged_schema.update(oschema)
    
    assert correct_tss_coord is None or isinstance(correct_tss_coord, int), "correct_tss_coords must be None or an integer that defines the correction coordinate"

    assert df['sample_id'].n_unique() == 1, "Please provide a data frame pre-filtered with a single sample_id. Best when this function is called within a group_by"

    assert df['raw_ibar'].n_unique() == 1, "Please provide a data frame pre-filtered with a single ibar. Best when this function is called within a group_by"

    assert 'kalhor_id' in df.columns, "Missing field: kalhor_id"

    foc_sample = df['sample_id'][0]
    foc_ibar= df['raw_ibar'][0]


    if df.shape[0]<min_group_size:
        return pl.DataFrame(schema=merged_schema)

    # we reduce the size of the tso sequence to reduce the specificity of the match
    tso_string = return_feature_string("tso")
    fuzzy_tso = fuzzy_match_str(tso_string[5:], wildcard=".{0,1}")
    raw_fasta_entries = (
            df
            .with_columns(
                # a parenthesis is needed for the fuzzy pattern to include the wildcards at the beginning
                pl.col('oseq').str.replace(f'^.+?({fuzzy_tso})', tso_string)
                .alias('seq_trim')
            )
            .with_columns(
                pl.len().over(['wt', 'seq_trim'])
                .alias('sweight')
            )
            .with_columns(
                [pl.col(i).rank(method="dense").cast(pl.Int64).name.suffix('_encoded') for i in ['seq_trim',]]
            )
            .with_row_index(name='id')
            .with_columns(
                (
                 "WT_"+pl.col('wt')\
                 +"_kalhor_"+pl.col('kalhor_id').cast(pl.Utf8).fill_null('NA')\
                 +"_sweight_"+pl.col('sweight').cast(pl.Utf8)\
                 +"_id_"+pl.col('id').cast(pl.Utf8)\
                 +"_seqid_"+pl.col('seq_trim_encoded').cast(pl.Utf8)\
                 )
                 .alias('id')
            )
            .with_columns(
                (">"+pl.col('id')\
                 +"\n"+pl.col('seq_trim')\
                 +"\n")
                 .alias('fasta')
            )
            .drop('seq_trim_encoded')
        )

    ### delete?
    if len(raw_fasta_entries)==0:
        print("No entries found to generate a fasta file {foc_ibar} {foc_sample}")
        return pl.DataFrame(schema=merged_schema)

    foc_spacer = (
            raw_fasta_entries
            .filter(pl.col('wt'))['spacer']
            .value_counts()
            .sort('count', descending=True)
            )
    print(foc_spacer)

    #TODO change to use a defined spacer db
    if foc_spacer.shape[0]==0:
        print(f'No WT allele found for {foc_ibar} {foc_sample}')
        return pl.DataFrame(schema=merged_schema)

    foc_spacer = foc_spacer['spacer'][0]
    foc_ref = ref_str.replace('[SPACER]', foc_spacer).replace('[IBAR]', foc_ibar)

    # ibars could show different sizes as spacer are varible
    # we use the downstream ibar sequence to anchor the ibar position
    foc_ref_match = regex.search(f'{foc_ibar}{return_feature_string("dsibar")}', foc_ref)
    lim = lim if foc_ref_match is None else foc_ref_match.start()

    coord_bins = range(lim)

    # perform the alignment
    alignment_tuples = alignpw_ibar_emboss(
            raw_fasta_entries=raw_fasta_entries,
            foc_sample = foc_sample,
            foc_ibar =  foc_ibar,
            foc_ref =  foc_ref,
            verbose = verbose,
            keep_intermediate =  keep_intermediate)

    
    # create alignment data frame
    alg_df = []

    for ((ref_name, read_name), (ref_seq, aligned_seq)) in alignment_tuples:
        if correct_tss_coord is not None:
            ref_seqv = np.array(list(ref_seq))
            aligned_seqv = np.array(list(aligned_seq))
            c_idx = range(correct_tss_coord) 
            aligned_seqv[c_idx]=ref_seqv[c_idx]
            aligned_seq = ''.join(aligned_seqv)

        expansion_factor = regex.search(".+sweight_(.+)_id", read_name)
        expansion_factor = int(expansion_factor.groups()[0])

        # when reads don't have insertions
        if "-" not in ref_seq:
            alg_mask = alg.return_del_mask(aligned_seq, 200)
            aseq =  aligned_seq[0:lim]
            #{"sample_id":str, "ibar":str, "id":str, "aseq":str, "sweight":pl.Int64, "alg":pl.List(pl.Int64)}
            alg_df.append((foc_sample, foc_ibar, read_name, aseq, expansion_factor, alg_mask))
        # cases where there are insertions
        # add an additional field iseq to concatenate a string of integrations, 
        if "-" in ref_seq:
            alg_mask = alg.return_ins_mask(ref_seq, 200)
            # TODO check the policy for trimming insertions
            aseq =  aligned_seq[0:lim]
            alg_df.append((foc_sample, foc_ibar, read_name, aligned_seq, expansion_factor, alg_mask))

    alg_df = to_interval_df(alg_df, schema=schema, sort_by='sweight') 

    raw_fasta_entries = (
        raw_fasta_entries
        .join(alg_df, left_on='id', right_on='id', how='left')
        .with_columns(pl.count().over('aseq').alias('aweight'))
        .select(merged_schema.keys())
        .with_columns(pl.col('aweight').cast(pl.Int64))
    )
    return raw_fasta_entries

def alignpw_ibar_emboss(
        raw_fasta_entries:pl.DataFrame,
        foc_ibar:str,
        foc_sample:int,
        foc_ref:str,
        keep_intermediate:bool=False,
        verbose:bool=False,
        conda_prefix:str = "conda run -n emboss")-> Sequence:
    '''
    Assumes a conda env named emboss which contains the emboss suite

    To create one:
        conda create -n emboss -y -c bioconda  emboss
        
    Returns: a tuple of tuples in the form: 
            ( (query_name, ref_name), (query_seq,ref_seq)) 
        from an alignment from a fasta file
    '''
    with tempfile.TemporaryDirectory() as temp_dir:
        run_id = f'{foc_sample}_{foc_ibar}'
        reads_fa = f'{temp_dir}/{run_id}_in_pair_algn.fa'
        ref_fa = f'{temp_dir}/{run_id}_ref_pair_algn.fa' 
        final_alignment = f'{run_id}_pair_algn.fa'

        with open(reads_fa, 'wt') as reads_out, open(ref_fa, 'wt') as ref_out:
            ref_out.write(f">ref\n{foc_ref}\n")
            for i in raw_fasta_entries['fasta'].to_numpy():
                if i is not None:
                    reads_out.write(i)

        ltr_utils.mltbc_align_reads_to_ref(
                name=run_id, 
                fa_ofn=final_alignment, 
                ref_path=ref_fa, 
                ref_name='ref', 
                conda_prefix=conda_prefix,
                verbose=verbose,
                wd=temp_dir)

        #TODO clean the messy paths related to the temp_dir
        alignment_tuples = ltr_utils.return_alignment_tuples(f'{temp_dir}/{final_alignment}') 
        if keep_intermediate:
            os.system(f'cp -r {temp_dir} ./{run_id}')

    return alignment_tuples

def to_interval_df(alignment_data, schema, sort_by):
    '''
        Converts a list of tuples containing intervals data (``alignment_data``) into  a data frame following a ``schema``.
    '''
    df = pl.DataFrame(alignment_data, schema=schema).sort(sort_by, descending=True)
    return df
