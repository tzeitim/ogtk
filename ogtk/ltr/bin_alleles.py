import ogtk
import pickle
import subprocess
import pyaml
import itertools
import pyfaidx
import os
import multiprocessing
import itertools
import regex
import numpy as np
import pandas as pd
import pdb
from functools import wraps
from .ltr_mtl import mltb_bin_alleles 


class Bin_Alleles():
    def __init__(self, name, intid, config_card_dir, outdir, fqm, bint_db_ifn, kmer_correct = True,  min_reads_per_umi = 4, intid2_R2_strand = None, threads = 100, end = 5000, ranked_min_cov = 5, consensus = False, umi_errors = 1, debug = False, jobs_corr=10, alg_gapopen = 20, alg_gapextend =1, correction_dict_path = None, correction_dict = None, trimming_pattern = None, use_cache = True):
        self.name = name
        self.intid = intid
        self.config_card_dir = config_card_dir
        self.outdir = outdir
        self.fqm = fqm
        self.bint_db_ifn = bint_db_ifn


