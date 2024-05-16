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

def return_del_mask(aseq, lim): 
    '''
    '''
    alg_mask = np.zeros(lim, dtype=int)
    cur = 0
    for match in regex.finditer("-+", aseq):
        operation_len = np.diff(match.span())[0]
        alg_mask[range(match.start(), match.end())] = -1
        cur = operation_len-1
    return pl.Series(alg_mask)

def return_ins_mask(ref_seq, lim): 
    '''
    '''
    alg_mask = np.zeros(lim, dtype=int)
    cur = 0
    for match in regex.finditer("-+", ref_seq):
        operation_len = np.diff(match.span())[0]
        alg_mask[match.start()-cur] = operation_len
        cur = operation_len-1
    return pl.Series(alg_mask)


def to_intervals(pos_array, sample_id, lim=50):
    xy = (
            pl.Series(pos_array, dtype=pl.Float64)
            .append(pl.Series(range(lim)).cast(pl.Float64))
            .value_counts()
            .with_columns(pl.col('count')-1)
            .rename({'':'pos', "counts":sample_id})
            .with_columns(pl.col("pos").cast(pl.Float64))
            .with_columns(pl.col(sample_id)/len(pos_array))
            )
    return xy

def explode_pos(df):
    dsource = (
            df
            .drop(['alg', 'aseq', 'weight'])
            .melt(id_vars=['sample_id', 'ibar', 'id'], variable_name='pos')
            .with_columns(pl.col('pos').str.replace('pos_', '').cast(pl.Int64))
        )
    return dsource

def return_exploded_positions(df, lim=50): 
    '''
    '''
    fdata = (
            df
            .with_columns([pl.col('alg').list.get(i).alias(f"pos_{i}") for i in range(lim)])
            .with_columns( pl.col("^pos.+$") * pl.col('aweight'))
        )
    return fdata


