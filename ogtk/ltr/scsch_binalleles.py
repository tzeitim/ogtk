## single-cell single channel bin alleles 
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
#import os
#import glob
#import pandas as pd
#import numpy as np
#import regex
#import itertools
#import pdb
#import time, random
#import multiprocessing
#import gzip 
#import anndata as an
#import ltr.ltr_utils as utils
##import .ibars


def scsch_bin_alleles(name, intid, 
    config_card_dir,
    outdir, 
    fqm, 
    h5ad_ifn,
    bint_db_ifn, 
    kmer_correct = True,  
    min_reads_per_umi = 4, 
    intid2_R2_strand = None, 
    threads = 100, 
    end = 5000, 
    ranked_min_cov = 5, 
    consensus = False, 
    umi_errors = 1,
    debug = False, 
    jobs_corr=10, 
    alg_gapopen = 20,
    alg_gapextend =1, 
    correction_dict_path = None, 
    correction_dict = None, 
    trimming_pattern = None, 
    use_cache = True):
    
    '''Once 10x fastqs (fqm) have been split by integration, merged into a single read and appended by the UMI,
    process them into a table of barcodes for downstream analysis
    Two modalities to choose for a representative allelle: ranked (consensus = False) 
    or consensus (consensus = True)
    '''
    utils.create_workdir(outdir, config_card_dir) 

    out_prefix = f'{outdir}/{name}'

    bint_db = pyaml.yaml.load(open(bint_db_ifn), Loader=pyaml.yaml.FullLoader)

    conf_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
    conf = {'name':name, 'desc':{}, 'outdir':outdir, 'intid1':intid, 'modality':'sc'}

    ref_seq = bint_db['ref_seq'].replace('{INTID1}', intid).upper()
    ref_seq = ref_seq.replace('PROTOSPACER', 'N' * 21) #TODO implement ibar db?
    ref_name = f'{bint_db["name_plasmid"]}_intid_{intid}'
    ref_path = f'{outdir}/{ref_name}.fa'
 
    anchor2 = bint_db['anchor_plasmid2']
    rxanch2 = regex.compile("ANCHOR2{e<=3}".replace("ANCHOR2", anchor2))
################ we were here
    adata = an.read_h5ad(h5ad_ifn).obs.index
    #r1 = '/local/users/polivar/src/artnilet/datain/20211217_scv2/direct_B3_shRNA_S11_R1_001.fastq.gz'
    sorted_tab = ogtk.utils.tabulate_10x_fastqs(r1 = fqm, cbc_len=16, 
    umi_len=10, end = 1000)
    n_cells = 1000
    ibars.rs_ibars_pc(name = name, 
                    tbx_ifn = sorted_tab, 
                    cell_bcs = anns[key], 
                    ibar_wl = ibar_wl, 
                    end = n_cells)
    return(0) 
    ## for 10x 5' kit
    umi_len = 26

    #### config
    conf = conf_dict_write(
        conf,
        level = 'inputs',
        fastq_merged = fqm
    )
    conf = conf_dict_write(
        conf,
        level = 'filters',
        reads_sampled = end,
        ranked_min_cov = ranked_min_cov,
        min_reads_per_umi = min_reads_per_umi,
        umi_errors = umi_errors,
        trimming_pattern = trimming_pattern,
        kmer_correct = kmer_correct
        )

    conf = conf_dict_write(conf, level='desc', 
        alignment =  
        "Parameters for the pairwise alignment of every recovered allele "+
        "and the reference. The corrected files fix the issue of duplicated"+ 
        "fasta entry names"
        )
    #### pair-wise alignment
    conf = conf_dict_write(
        conf, 
        level='alignment', 
        fa_correctedm =f'{out_prefix}_corrected_merged.fa',
        alg_mode = 'needleman',
        alg_gapopen = alg_gapopen,
        alg_gapextend = alg_gapextend,
        ref_path  = f'{outdir}/{ref_name}.fa',
        bint_db_ifn = bint_db_ifn,
        ref_name = ref_name,
        )
                
   #### cache
    conf = conf_dict_write(
        conf, 
        level = 'cache',
        readset = f'{out_prefix}_readset.corr.pickle',
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
        consensus =    f'{out_prefix}_consensus_by_cell_10x.txt',
        allele_reps =  f'{out_prefix}_allele_reps.pickle', 
        merged_tab =   f'{out_prefix}_binned_barcodes_10x.txt',
        merged_full =  f'{out_prefix}_binned_barcodes_10x_full.txt'
        )
   
    # Molecule data
    conf = conf_dict_write(
        conf,
        level = 'mols',
        umi_len = umi_len,
        counts =  outdir + f'{out_prefix}_umi_counts.txt',
        cov_pck = outdir + f'{out_prefix}_umi_cov.pickle'
        )

    print("creating fasta ref")
    create_fasta_ref(conf['alignment']['ref_name'], ref_seq, conf['alignment']['ref_path'])


    #>>>
    _sc_bin_alleles(conf_fn, conf, use_cache = use_cache, threads = threads, correction_dict_path = correction_dict_path, correction_dict = correction_dict)


def _sc_bin_alleles(conf_fn, conf, **kwargs):
    '''
        main function for binning alleles. For reference see sc_bin_alleles 
    '''
    intid = conf['intid1']
    pickled_readset = conf['cache']['readset']
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
    kmer_correct = conf['filters']['kmer_correct']

    # from kwargs
    threads = kwargs['threads']
    use_cache =kwargs['use_cache'] 
    correction_dict_path = kwargs['correction_dict_path']
    correction_dict = kwargs['correction_dict']
    
    if end != None:
        print("Warning: end is not None; You are not processing all the data", end)

    # do a pass on the raw fastqs and group reads by UMI
    if os.path.exists(pickled_readset) and use_cache:
        rssc = pickle.load(open(pickled_readset, 'rb'))
        print(f"!!! loaded cached ReadSet !!! {pickled_readset}. total umis: {len(rssc.umis)}", end = '\n\n')
    else:    
        rssc = ogtk.UM.fastq_collapse_UMI(
                        conf['inputs']['fastq_merged'], 
                        umi_len = conf['mols']['umi_len'], 
                        end = end, 
                        keep_rid = True, 
                        min_reads_per_umi = conf['filters']['min_reads_per_umi'],
                        trimming_pattern = trimming_pattern)

        print(f"loaded rs with trimming={trimming_pattern != None}. total umis: {len(rssc.umis)}")
        if umi_errors>0:
            print("Correcting umis with a hdist of {}".format(umi_errors))
            rssc.correct_umis(errors = umi_errors , silent = False, jobs = jobs_corr)
            conf['mols']['saturation']['unmerged']['corrected'] = rssc.saturation()
        else:
            print("Not correcting umis")
        
        if kmer_correct:
            kmer = 50
            print('starting kmer-based correction k={kmer}')
            umi_list = rssc.umi_counts()
            umi_list = umi_list[umi_list > 1]
            umi_list = rssc.umi_counts().index
            cseqs = ogtk.UM.monitor_kmer_corr_umi_list(umi_list, rssc, pickle_ofn=pickled_cseqs, k = kmer)
            for corr_umi in cseqs.keys():
                quals = ['F'*len(seq) for seq in cseqs[corr_umi]]
                rssc.umis[corr_umi].seqs = cseqs[corr_umi]
                rssc.umis[corr_umi].quals = quals

    conf = store_molecule_saturation_stats(conf, rssc)

    # The first step to call an allele consensus is to get a representative sequence at the molecule level
    print(f'getting rank1 seq per UMI')
    ogtk.UM.do_fastq_pileup(rssc, 
                    min_cov = ranked_min_cov, 
                    threads = threads, 
                    trim_by = None, 
                    by_alignment = False)

    print(f'determining valid cells')
    cell_umi_dict = ogtk.UM.get_list_of_10xcells(
                    rssc, 
                    correction_dict_path = correction_dict_path, 
                    correction_dict = correction_dict)

    print(f'collecting cell allele candidates', end = '...')
    cell_alleles = cell_dominant_allele(rssc, cell_umi_dict, allele_reps_fn)

    # TODO qc plot
    if False: #or plot_qcs
        def dist_rep_allele(x):
            top_seqs = x['dominant_seq'].value_counts()
            return(top_seqs.to_list()[0])

        plt.close('all')
        yd = cell_alleles.groupby('cell').apply(lambda x: dist_rep_allele(x)).value_counts(normalize = True)
        xd = yd.index
        plt.scatter(xd, np.cumsum(yd.to_list()))

    print(f'determining best candidate')
    cell_consensus = cell_alleles.groupby('cell').apply(lambda x: get_rep_allele(x, min_mol = 0, min_domix = 0.0))
    print(f'a total of discarded entries {np.sum(cell_consensus == False)}/{len(cell_consensus)} {np.sum(cell_consensus == False)/len(cell_consensus):.2f}')
    cell_consensus = cell_consensus[cell_consensus != False]

    cell_consensus.to_csv(consensus_tab_out, sep ='\t')
    #>>>>>>>>>
    #print(f"called consensuses {len(rssc.consensus)}")
    #cov_df = rssc.compute_coverage_table(cell_umis)

    # as a control one can recompute the reads back from the UMI cover 
    #rep_read = [(len(rssc.umis[i].seqs), rssc.consensus[i][1]) for i in rssc.consensus.keys()]
    #import matplotlib.pyplot as plt
    #plt.scatter(*zip(*rep_read))
    #reads_rank1 = pd.Series([v[2] for i,v in rssc.consensus.items()]).value_counts(normalize= False)

    # determine cells based on the 10x bc
    # we might need to filter at this level for cells that have more than one molecule ... but what about reads?
    #pdb.set_trace()
    
    if False:
        debug_iter = [i for i in job_iter if len(i[2])>3 ]
        print("There are {} cells and {} with more than 3 representatives".format(len(job_iter), len(debug_iter)))
        for i in debug_iter:
            ogtk.UM.mafft_consensus(i)

    if len(cell_consensus) == 0 :
        conf['success'] = False
        print(f'Failed to generate representative sequences for each cell')
        return(None)
    else:
        print(f'Success')
        conf['success'] = True
 
    ### cache: save binaries
    print(f'pickling rs as {pickled_readset}')
    with open(pickled_readset, 'wb') as pcklout:
        pickle.dump(rssc, pcklout)

    #print(f'pickling umi coverage info as {umi_cov_ofn}')
    #cov_df.to_pickle(umi_cov_ofn)

    print("aligning cell alleles to ref")
    align_reads_to_ref(
            name = conf['outdir']+"/"+intid, 
            fa_ofn = conf['alignment']['fa_correctedm'], 
            consensus_dict = cell_consensus.to_dict(), 
            ref_name = conf['alignment']['ref_name'],
            ref_path = conf['alignment']['ref_path'],
            mode = conf['alignment']['alg_mode'],
            gapopen = conf['alignment']['alg_gapopen'], 
            gapextend = conf['alignment']['alg_gapextend'])

    print(f"featurizing sequences to {merged_tab_out}")
    compute_barcode_matrix_merged(
                fa_ifn =conf['alignment']['fa_correctedm'], 
                tab_out = merged_tab_out, 
                bint_db_ifn = conf['alignment']['bint_db_ifn'], 
                do_rc = False)
   
    # save the count information for the UMIs 
    with open(umi_counts_ofn, 'w') as fh:
        fh.write(f'count\tumi\n')
        counts = rssc.umi_counts()
        for count,umi in zip(counts, counts.index):
            fh.write(f'{count}\t{umi}\n')
                
    # save run settings 
    with open(conf_fn, 'w') as yaml_stream:
        pyaml.yaml.dump(conf, yaml_stream)
    print("Run settings saved to config card:{}".format(conf_fn))
    return(conf_fn)
  
