#import ogtk
#import pickle
#import subprocess
#import pyaml
#import itertools
#import pyfaidx
#import os
#import multiprocessing
#import itertools
#import regex
#import numpy as np
#import pandas as pd
#import pdb
#
class Ltr_run():
    def __init__(self, name, intid, config_card_dir, outdir, fqm, load_defaults = None, **kwargs):
        self.name = name
        self.intid = intid
        self.config_card_dir = config_card_dir
        self.outdir = outdir
        self.fqm = fqm
        self.__dict__.update(kwargs)

        if load_defaults is not None:
            self.__dict__.update(self.list_of_options(load_defaults))


    def list_of_options(self, recipe):
        ''' returns default arguments of parameters according to recipe '''
        if recipe not in ['sc_bin', 'bulk_bin', 'sc_shrna']:
            raise ValueError("recipe type not supported")

        if recipe == 'sc_bin':
            # TODO nos verified 
            arguments = {
            'bint_db_ifn' : None,\
            'kmer_correct' : True,\
            'min_reads_per_umi' : 4,\
            'intid2_R2_strand' : None,\
            'threads' : 100,\
            'end' : 5000,\
            'ranked_min_cov' : 5,\
            'consensus' : False,\
            'umi_errors' : 1,\
            'debug' : False,\
            'jobs_corr' : 10,\
            'alg_gapopen' : 20,\
            'alg_gapextend' : 1,\
            'correction_dict_path' : None,\
            'correction_dict' : None,\
            'trimming_pattern' : None,\
            'use_cache': True}

        if recipe == 'bulk_bin':
            # TODO nos verified 
            arguments = {
            'bint_db_ifn' : None,\
            'kmer_correct' : True,\
            'min_reads_per_umi' : 4,\
            'intid2_R2_strand' : None,\
            'threads' : 100,\
            'end' : 5000,\
            'ranked_min_cov' : 5,\
            'consensus' : False,\
            'umi_errors' : 1,\
            'debug' : False,\
            'jobs_corr' : 10,\
            'alg_gapopen' : 20,\
            'alg_gapextend' : 1,\
            'correction_dict_path' : None,\
            'correction_dict' : None,\
            'trimming_pattern' : None,\
            'use_cache': True}

        if 'sc_shrna':
            arguments = {
            'bint_db_ifn' : None,\
            'kmer_correct' : True,\
            'min_reads_per_umi' : 4,\
            'intid2_R2_strand' : None,\
            'threads' : 100,\
            'end' : 5000,\
            'ranked_min_cov' : 5,\
            'consensus' : False,\
            'umi_errors' : 1,\
            'debug' : False,\
            'jobs_corr' : 10,\
            'alg_gapopen' : 20,\
            'alg_gapextend' : 1,\
            'correction_dict_path' : None,\
            'correction_dict' : None,\
            'trimming_pattern' : None,\
            'use_cache': True}

        return(arguments)
