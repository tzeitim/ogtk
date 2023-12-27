# Rename to : bulk multi-channel (bmch_binalleles)
import ogtk
import pickle
import subprocess
import pyaml
import os
import itertools
import os
import pandas as pd
import numpy as np
import regex
import itertools
import pdb
import time, random
import multiprocessing
import gzip 

from .ltr_utils import *

def bulk_bin_alleles(name, intid, intid2_R2_strand, 
                     config_card_dir, 
                     outdir, 
                     fq1, 
                     fq2, 
                     bint_db_ifn, 
                     min_reads_per_umi = 4, 
                     threads = 100, 
                     bbthreads = 8, 
                     end = 5000, 
                     ranked_min_cov = 5, 
                     umi_len = 17, 
                     umi_errors = 1,
                     trimming_pattern = None, 
                     alg_gapopen = 20, 
                     alg_gapextend =1,
                     process_unmerged = True,
                     verbose = False,
                     use_cache = True):
    '''Once bulk fastq files (R1 and R2) have been split (RPI>rev_id>intid) process them into a table of barcodes 
    for downstream analysis.
    Two modalities to choose for a representative allelle: ranked (consensus = False) 
    or consensus (consensus = True)
    '''

    # TODO add tadpole-based correction as an option
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
        print("Created {}".format(outdir))

    if not os.path.isdir(config_card_dir):
        os.makedirs(config_card_dir, exist_ok=True)
        print("Created {}".format(config_card_dir))
    
    out_prefix = f'{outdir}/{name}'
    if not os.path.isdir(out_prefix):
        os.makedirs(out_prefix, exist_ok=True)
        print("Created {}".format(out_prefix))

    if intid2_R2_strand != None:
        intid2 = ogtk.UM.rev_comp(intid2_R2_strand)# same strand as initd1
    else:
        intid2 = len(intid) * "N"
   
    bint_db = pyaml.yaml.load(open(bint_db_ifn), Loader=pyaml.yaml.FullLoader)

    conf_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
    conf = {'name':name, 'desc':{}, 'outdir':outdir, 'intid1':intid, 'intid2': intid2, 'modality':'bulk'}

    ref_seq = bint_db['ref_seq'].replace('{INTID1}', intid).replace('{INTID2}', intid2).upper()
    ref_name = f'{bint_db["name_plasmid"]}_intid_{intid}_{intid2}'
    ref_path  = f'{out_prefix}/{ref_name}.fa'
    rcref_path = f'{out_prefix}/rc_{ref_name}.fa'

    #### config
    conf = conf_dict_write(
        conf,
        level = 'inputs',
        fq1 = fq1,
        fq2 = fq2
    )
    ### bbmerge
    conf = conf_dict_write(
        conf,
        level = 'bbmerge',
        fastq_merged = f'{out_prefix}_merged.fastq.gz',
        fastq_umerged1 = f'{out_prefix}_unmerged_R1.fastq.gz',
        fastq_umerged2 = f'{out_prefix}_unmerged_R2.fastq.gz',
        ihist =f'{out_prefix}_ihist.txt',
        k = 60,
        threads=bbthreads,
        Xmx = '2g',
        modality = 'ecct', # error-correct with Tadpole.sh
        log = f'{out_prefix}_merge_stats.txt',
    )
    cmd_bbmerge = f"bbmerge-auto.sh in1={fq1} in2={fq2} out={conf['bbmerge']['fastq_merged']} \
        k={conf['bbmerge']['k']} \
        outu={conf['bbmerge']['fastq_umerged1']} \
        outu2={conf['bbmerge']['fastq_umerged2']} \
        ihist={conf['bbmerge']['ihist']} \
        t={conf['bbmerge']['threads']} \
        {conf['bbmerge']['modality']} \
        extend2=20 iterations=5 interleaved=false -Xmx{conf['bbmerge']['Xmx']}"
 
    conf = conf_dict_write(
        conf,
        cmd = cmd_bbmerge
    )

    ### filters
    conf = conf_dict_write(
        conf,
        level = 'filters',
        reads_sampled = end,
        ranked_min_cov = ranked_min_cov,
        min_reads_per_umi = min_reads_per_umi,
        umi_errors = umi_errors,
        trimming_pattern = trimming_pattern
        )
#### pair-wise alignment
    conf = conf_dict_write(conf, level='desc', 
        alignment =  
        "Parameters for the pairwise alignment of every recovered allele "+
        "and the reference. The corrected files fix the issue of duplicated"+ 
        "fasta entry names"
        )
    
    conf = conf_dict_write(
        conf, 
        level='alignment', 
        fa_correctedm =f'{out_prefix}_corrected_merged.fa',
        fa_corrected1 =f'{out_prefix}_corrected_R1.fa',
        fa_corrected2 =f'{out_prefix}_corrected_R2.fa',
        alg_mode = 'needleman',
        alg_gapopen = alg_gapopen,
        alg_gapextend = alg_gapextend,
        ref_path  = ref_path,
        ref_seq  = ref_seq,
        verbose = verbose,
        bint_db_ifn = bint_db_ifn,
        rcref_path = rcref_path,
        ref_name = ref_name,
        )
               
   #### cache
    conf = conf_dict_write(
        conf, 
        level = 'cache',
        merged_readset = f'{out_prefix}_mreadset.corr.pickle',
        readset1 = f'{out_prefix}_readset1.corr.pickle',
        readset2 = f'{out_prefix}_readset2.corr.pickle',
        cseqs = f'{out_prefix}_readset.kmer_cseqs.pickle'
        )
   #### main output
    conf = conf_dict_write(
        conf, 
        level='desc', 
        lineage = "Main output of the preprocessing scripts. Tabulated files that contain the bint barcode string",
        )

    conf = conf_dict_write(
        conf,
        level ='lineage',
        consensus =    f'{out_prefix}_consensus_by_cell.txt',
        allele_reps =  f'{out_prefix}_allele_reps.pickle', 
        merged_tab =   f'{out_prefix}_merged_binned.txt',
        merged_full =  f'{out_prefix}_binned_full.txt',
        paired_tab =   f'{out_prefix}_paired_binned.txt',
        paired_full =  f'{out_prefix}_paired_binned_full.txt', # automatically defined by mlt_compute_barcode_matrix_paired
        )
   
    # Molecule data
    conf = conf_dict_write(
        conf,
        level = 'mols',
        umi_len = umi_len,
        counts =  f'{out_prefix}_umi_counts.txt',
        cov_pck = f'{out_prefix}_umi_cov.pickle'
        )
  
    # merge paired-end reads if possible
    log_merge =conf['bbmerge']['log'] 
    pp = subprocess.run(cmd_bbmerge.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(f"merge log\n{log_merge}")
    with open(log_merge, 'wb') as fout:
        fout.write(pp.stdout)
        fout.write(pp.stderr)
    
    _bulk_bin_alleles(conf_fn, conf, use_cache = use_cache, threads = threads, process_unmerged = process_unmerged)
    return(conf_fn)

def _bulk_bin_alleles(conf_fn, conf, **kwargs):
    ''' main function for binning alleles from bulk data. For reference see bulk_bin_alleles
    '''
    outdir = conf['outdir']
    name = conf['name']
    intid = conf['intid1']
    intid2 = conf['intid2']
    intid2_rc = ogtk.UM.rev_comp(conf['intid2'])
    pickled_merged_readset = conf['cache']['merged_readset']
    pickled_readset1 = conf['cache']['readset1']
    pickled_readset2 = conf['cache']['readset2']
    pickled_ceqes = conf['cache']['cseqs']
    allele_reps_fn = conf['lineage']['allele_reps']
    consensus_tab_out = conf['lineage']['consensus']
    umi_counts_ofn = conf['mols']['counts']
    umi_cov_ofn = conf['mols']['cov_pck']
    merged_tab_out = conf['lineage']['merged_tab']
    end = conf['filters']['reads_sampled']
    ranked_min_cov = conf['filters']['ranked_min_cov']
    trimming_pattern = conf['filters']['trimming_pattern']
    umi_errors = conf['filters']['umi_errors']
    ref_seq = conf['alignment']['ref_seq']
    ref_name = conf['alignment']['ref_name']
    ref_path = conf['alignment']['ref_path']
    rcref_path = conf['alignment']['rcref_path']
    # from kwargs
    threads = kwargs['threads']
    use_cache =kwargs['use_cache'] 
    process_unmerged = kwargs['process_unmerged']

    if end != None:
        print("Warning: end is not None; You are not processing all the data", end)

    rs_paths = [pickled_merged_readset, pickled_readset1, pickled_readset2]

    # do a pass on the raw fastqs and group reads by UMI
    if all([os.path.exists(i) for i in rs_paths]) and use_cache:
        rsm = pickle.load(open(pickled_merged_readset, 'rb'))
        rs1 = pickle.load(open(pickled_readset1, 'rb'))
        rs2 = pickle.load(open(pickled_readset2, 'rb'))
        print(f"!!! loaded cached ReadSets!!!: {'_'.join(rs_paths)}", end = '\n\n')
        print(f"merged total umis: {len(rsm.umis)}")
        print(f"R1 total umis: {len(rs1.umis)}")
        print(f"R2 total umis: {len(rs2.umis)}")
    else:  
    ## for merged (rsm) and unmerged (rs1, rs2)
        rsm = ogtk.UM.fastq_collapse_UMI(
            conf['bbmerge']['fastq_merged'],
            umi_len = conf['mols']['umi_len'],
            end = end, 
            keep_rid = True,
            min_reads_per_umi = conf['filters']['min_reads_per_umi'],
            trimming_pattern = trimming_pattern)

        print(f"loaded rs with trimming={trimming_pattern != None}. total umis: {len(rsm.umis)}")

        rs1, rs2 = ogtk.UM.pfastq_collapse_UMI(
            conf['bbmerge']['fastq_umerged1'],
            conf['bbmerge']['fastq_umerged2'],
            umi_len = conf['mols']['umi_len'],
            min_reads_per_umi = conf['filters']['min_reads_per_umi'],
            end = end)

    # remove umis that exist on the merged set out of the unmerged
    unmerged_by_error = set(rsm.umi_list()).intersection(set(rs1.umi_list()))
    rs1.delete(unmerged_by_error)
    rs2.delete(unmerged_by_error)

    #TODO the saturation stats should be outsourced somehow and reported for all three cases not just the merged (rsm) 
    ### TODO for the refactored version too

    if umi_errors >0:
        print(f"Correcting umis with a hdist of {umi_errors} using {threads} cores")
        rsm.correct_umis(errors = umi_errors , silent = True, jobs = threads)
        rs1.correct_umis(errors = umi_errors , silent = True, jobs = threads)
        rs2.correct_umis(errors = umi_errors , silent = True, jobs = threads)

    else:
        print("Not correcting umis")
    
    rsm.consensus = ogtk.UM.do_fastq_pileup(rsm, min_cov = ranked_min_cov, threads = threads, trim_by = None)
    rs1.consensus = ogtk.UM.do_fastq_pileup(rs1, min_cov = ranked_min_cov, threads = threads, trim_by = None)
    rs2.consensus = ogtk.UM.do_fastq_pileup(rs2, min_cov = ranked_min_cov, threads = threads, trim_by = None)

    # align allele to reference. merged and unmerged separately
    ## make sure to include the intid on the reference to help the alignment
    mlt_create_fasta_ref(ref_name, ref_seq, ref_path)
    ## also on rev comp (not sure if this is really needed)
    mlt_create_fasta_ref(ref_name+"RC", ogtk.UM.rev_comp(ref_seq), rcref_path)

    alg_mode = 'needleman'

    alg_mode = conf['alignment']['alg_mode']
    alg_gapopen = conf['alignment']['alg_gapopen']
    alg_gapextend = conf['alignment']['alg_gapextend']
    fa_correctedm = conf['alignment']['fa_correctedm']
    fa_corrected1 = conf['alignment']['fa_corrected1']
    fa_corrected2 = conf['alignment']['fa_corrected2']
    alg_verbose = conf['alignment']['verbose'] 
    bint_db_ifn = conf['alignment']['bint_db_ifn']

    # badly covered umis have been internally removed from the readsets
    if len(rsm.umis)!=0:
        mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"m", 
                                 fa_ofn = fa_correctedm, 
                                 rs = rsm, ref_path = ref_path, 
                                 verbose = alg_verbose,
                                 ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                                 gapextend = alg_gapextend)
        compute_barcode_matrix_merged(fa_ifn = fa_correctedm, tab_out= merged_tab_out, bint_db_ifn = bint_db_ifn, do_rc=False)

    if len(rs1.umis)!=0:
        mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"1", 
                                 fa_ofn = fa_corrected1, 
                                 rs = rs1, ref_path = ref_path, 
                                 verbose = alg_verbose,
                                 ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                                 gapextend = alg_gapextend)

    if len(rs2.umis)!=0:
        mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"2", fa_ofn = fa_corrected2, 
                                 rs = rs2, ref_path = rcref_path, 
                                 verbose = alg_verbose,
                                 ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                                 gapextend = alg_gapextend)


    # TODO PO verify if we need a list of corrected umis of whether they had been deleted 
    rs1c = list(rs1.umis.keys())
    rs2c = list(rs2.umis.keys())

    umi_whitelist = list(set(rs1c).intersection(set(rs2c)))

    print("whitelist:", len(umi_whitelist))

    if len(umi_whitelist)>0 and process_unmerged and (len(rs1.umis)!=0 and len(rs2.umis)!=0):
        unmerged_tab_out = conf['lineage']['paired_tab'] 

        print(f'exploring {fa_corrected1}')
        print(f'exploring {fa_corrected2}')

        mlt_compute_barcode_matrix_paired(fa_ifn = fa_corrected1, 
                                           fa_ifn2 = fa_corrected2, 
                                           umi_whitelist = umi_whitelist, 
                                           umi_blacklist = list(rsm.umis.keys()), 
                                           tab_out = unmerged_tab_out, 
                                           bint_db_ifn = bint_db_ifn, 
                                           do_rc = True)

    ### cache: save binaries
    for rs, pickled_readset in zip([rsm, rs1, rs2], rs_paths):
        print(f'pickling rs as {pickled_readset}')
        with open(pickled_readset, 'wb') as pcklout:
            pickle.dump(rs, pcklout)


    # save the count information for the UMIs 
    with open(umi_counts_ofn, 'w') as fh:
        fh.write(f'count\tumi\n')
        for i in (rsm, rs1):
            counts = i.umi_counts()
            for count,umi in zip(counts, counts.index):
                fh.write(f'{count}\t{umi}\n')
    
    # save run settings 
    with open(conf_fn, 'w') as yaml_stream:
        pyaml.yaml.dump(conf, yaml_stream)
        print("Run settings saved to config card:\n{}".format(conf_fn))
    return(conf_fn)
