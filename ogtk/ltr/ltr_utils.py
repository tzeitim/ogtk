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
import polars as pl
import os
import glob
import pandas as pd
import numpy as np
import regex
import itertools
import pdb
import time, random
import multiprocessing
import gzip 
#from .bin_alleles import Bin_Alleles



def extract_intid_worker(args):
    ''' input precompiled regex and sequences'''
    rxc, seqs = args
    hits = []
    for i in seqs:
        match = rxc.search(i)
        if match:
            hits.append(match.group('intid1'))
    return(hits)

def extract_intid_pair_worker(args):
    ''' input precompiled regex and sequences'''
    rxcs, seqs = args
    rxc, rxc2 = rxcs
    hits = []
    for i in seqs:
        match = rxc.search(i)
        if match:
            match2 = rxc2.search(i)
            if match2:
                hits.append(match.group('intid1') +":"+ match2.group('intid2'))
    return(hits)


def detect_integrations_dual(ifn, n_intid = None, intid_len = 8, 
    anchor1 = "CTGTTCCTGTAGAAAT", error1 = 3, \
    anchor2 = "CCGGACTCAGATCTCGAGCT", error2 = 2, \
    anchor3 = "CGAGCGCTATGAGCGACTATGGGA", error3 = 3, \
    limit = None, cores = 50, verbose = False):
    ''' process fastq/fastq.gz/bam file and detects number of integrations 
    by means of anchoring seqs (anchor1, anchor2, anchor3)
    Supports fuzzy matching via the errors argument.
    Number cores can be adjusted
    '''
    # Determine the iterator in order to fetch sequences
    if ifn.endswith("fastq"):
        with open(ifn) as ff:
            it = itertools.islice(ff, 1, limit, 4)
            seqs = [i for i in it]
           
    if ifn.endswith("fastq.gz"):
        import gzip
        with gzip.open(ifn, 'rt') as ff:
            it = itertools.islice(ff, 1, limit, 4)
            seqs = [i for i in it]

    if ifn.endswith("bam"):
        import pysam
        it= pysam.AlignmentFile(ifn)
        it = [i.seq for i in it]
        seqs = [i for i in it]
        
    # we trim the first 100bp to reduce the memory footprint
    #seqs = [i[0:100] for i in it]
    # not possible for paired mode
    # TODO this might introduce a bug 
    chunks = np.array_split(seqs, cores)
    
    rxc = regex.compile(".*({}){{e<={}}}(?P<intid1>.{{{}}}).*({}){{e<={}}}".format(anchor1, error1, intid_len, anchor2, error2))
    rxc2 = regex.compile(".*(?P<intid2>.{{{}}}).*({}){{e<={}}}".format(intid_len, anchor3, error3))

    pool = multiprocessing.Pool(cores)

    hits = np.array(pool.map(extract_intid_worker,  itertools.zip_longest([rxc], chunks, fillvalue=rxc)), dtype=object)
    hits_paired = np.array(pool.map(extract_intid_pair_worker,  itertools.zip_longest([[rxc, rxc2]], chunks, fillvalue=[rxc, rxc2])), dtype=object)

    pool.close()
    pool.join()
    filtered_hits_s = detect_integrations_filter_hits(hits, n_intid, verbose = verbose)
    filtered_hits_d = detect_integrations_filter_hits(hits_paired, n_intid, verbose = verbose)
    # TODO how to deal with both modalities in one?
    return(filtered_hits_d)

def detect_integrations_export_ntop_intids_json(intid_series, n_top, ofn):
    import pyaml
    import pandas as pd
    df = pd.DataFrame({'intid_code':[f'intid_{i:02d}' for i in range(n_top)], 'intid':intid_series.head(n_top).index})
    df.to_json(ofn)    
    
def detect_integrations_single(ifn, n_intid = None, intid_len = 8,\
        anchor1 = "CTGTTCCTGTAGAAAT", error1 = 3,\
        anchor2 = "CCGGACTCAGATCTCGAGCT", error2 = 2,\
        limit = None,\
        cores = 50,\
        do_rc = False, verbose = False):
    ''' For barcoded integrations with a single identifier.
        Process fastq/fastq.gz/bam file and detects number of integrations by means of anchoring seqs (anchor1, anchor2)
        Supports fuzzy matching via the errors argument.
        Number cores can be adjusted
    '''
    # Determine the iterator in order to fetch sequences
    if ifn.endswith("fastq"):
        with open(ifn) as ff:
            it = itertools.islice(ff, 1, limit, 4)
            seqs = [ogtk.UM.rev_comp(i.strip())+'\n' if do_rc else i for i in it]
           
    if ifn.endswith("fastq.gz"):
        import gzip
        with gzip.open(ifn, 'rt') as ff:
            it = itertools.islice(ff, 1, limit, 4)
            #seqs = [i for i in it]
            seqs = [ogtk.UM.rev_comp(i.strip())+'\n' if do_rc else i for i in it]

    if ifn.endswith("bam"):
        import pysam
        it= pysam.AlignmentFile(ifn)
        it = [i.seq for i in it]
        #seqs = [i for i in it]
        seqs = [ogtk.UM.rev_comp(i.strip())+'\n' if do_rc else i for i in it]
        
    # we trim the first 100bp to reduce the memory footprint
    #seqs = [i[0:100] for i in it]
    # not possible for paired mode
    # TODO this might introduce a bug 
    chunks = np.array_split(seqs, cores)
    
    rxc = regex.compile(".*({}){{e<={}}}(?P<intid1>.{{{}}}).*({}){{e<={}}}".format(anchor1, error1, intid_len, anchor2, error2))
    #rxc2 = regex.compile(".*(?P<intid2>.{{{}}}).*({}){{e<={}}}".format(intid_len, anchor3, error3))

    pool = multiprocessing.Pool(cores)

    hits = np.array(pool.map(extract_intid_worker,  itertools.zip_longest([rxc], chunks, fillvalue=rxc)), dtype=object)
    #hits_paired = np.array(pool.map(extract_intid_pair_worker,  itertools.zip_longest([[rxc, rxc2]], chunks, fillvalue=[rxc, rxc2])))

    pool.close()
    pool.join()
    filtered_hits = detect_integrations_filter_hits(hits, n_intid)
    return(filtered_hits)


#    hits = [item for sublist in hits for item in sublist]
#    x = pd.Series(hits).value_counts()
#    xx = pd.Series([int(np.log10(i)) for i in x], index = x.index)
#    
#    if n_intid != None:
#        return(x.head(n_intid))
#    else:
#        print(f'No expected number of integrations provided, returning an guess based on the top ranking intid counts. ')
#        valid_ints = [i for i in x[xx>(xx[0]-1)].index]
#        print("Found {} valid integrations".format(len(valid_ints)))
#        print(valid_ints)
#        return(valid_ints)
def detect_integrations_filter_hits(hits, n_intid = None, verbose=False):
    ''' Helper function to filter out valid hits from an integration number detection stepi
    '''
    hits = [item for sublist in hits for item in sublist]
    hit_counts = pd.Series(hits).value_counts()
    hit_lcounts = pd.Series([int(np.log10(i)) for i in hit_counts], index = hit_counts.index)

    if n_intid != None:
        return(hit_counts.head(n_intid))
    else:
        if verbose:
            print(f'No expected number of integrations provided, returning an guess based on the top ranking intid counts.')
        if len(hit_lcounts) == 0 :
            if verbose:
                print("No valid integrations found")
            return([])
        valid_hits = hit_counts[hit_lcounts>(hit_lcounts[0]-1)]
        valid_ints = [i for i in valid_hits.index]
        if verbose:
            print("Found {} valid integrations".format(len(valid_ints)))
        return(valid_ints)

def preprocessing_split_10xfastq_by_intid(valid_intids, intid_len, fq1, fq2, anchor1, anchor2, outdir, limit = None):
    '''Given a list of valid integration ids (5') it splits the fastq files accordingly. It also appends the 10x barcode to the split read'''
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    import gzip 

    f1 = gzip.open(fq1, 'rt') if fq1.endswith('fastq.gz') else open(fq1) 
    f2 = gzip.open(fq2, 'rt') if fq2.endswith('fastq.gz') else open(fq2) 

    rx1 = regex.compile("(?P<umi>.{28})")
    rx2_str = ".*(ANCHOR1){e<=3}(?P<intid1>.{INTIDLEN}).*(ANCHOR2){e<=3}"
    rx2_str = rx2_str.replace('ANCHOR1', anchor1)
    rx2_str = rx2_str.replace('ANCHOR2', anchor2)
    rx2_str = rx2_str.replace('INTIDLEN', str(intid_len))
    rx2 = regex.compile(rx2_str)
    limit = None
    undet = 0
    if limit != None:
        print("Warning this is a subset")
    i1 = itertools.islice(f1, 0, limit, 1)
    i2 = itertools.islice(f2, 0, limit, 1)
    
    out_files_dic2 = {}
    
    for intid in valid_intids:
        out_files_dic2[intid] = open("{}/intid_{}_indexed-raw_R1.fastq".format(outdir, intid), 'w')
    
    id =[]
    seq = []
    plus = []
    qual = []
    
    seen = {}
    for i,v in enumerate(zip(i1,i2)):
            if i%100000 == 0:
                print(f'\rundet/tot {undet}/{i}', end = '')
            l1 = v[0].strip()
            l2 = v[1].strip()
    
            # this corresponds to read id
            if i%4 == 0:
                id = [l1, l2]
            # this corresponds to read seq
            if i%4 == 1:
                seq = [l1, l2]
            # this corresponds to read '+'
            if i%4 == 2:
                plus = [l1, l2]
            # this corresponds to read qual
            if i%4 == 3:
                qual = [l1, l2]
    
            if qual:
                m2 = rx2.search(seq[1])
                if m2:
                    intid1 = m2.group('intid1')
                    if intid1 in valid_intids:
                        # TODO consider correcting the intids using hamming distance of 1 or even better, using the distribution of distances of the valid_intids
                        # TODO this might not be perfect
                        for valid in valid_intids:
                            dist = ogtk.UM.string_hamming_distance((intid1, valid))
                            if dist <= 3:
                                intid1 = valid
                        out_files_dic2[intid1].write("\n".join([id[1], seq[0]+seq[1], plus[1], qual[0]+qual[1]]) +'\n')
                    else:
                        undet +=1
    
                id = []
                seq = []
                plus = []
                qual = []
    
    print(f'\nundet {undet}/{i}')
    for i in valid_intids:
        out_files_dic2[i].close()
###########
########### BEGINNING
###########
def split_fastq_by_intid(valid_pairs_tab, fq1, fq2, dname, outdir, limit = None, len_sample_id = 6):
    valid_pairs = open(valid_pairs_tab).readlines()
    intid1_wl = {}
    intid2_wl = {}

    valids = {}
    for i,line in enumerate(valid_pairs):
        #pdb.set_trace()
        i1, i2 = line.strip().split('\t')
        key = i1+"_"+i2
        valids[key] = key
        #intid1_wl[i1] = i
        #intid2_wl[i2] = i
    #print(intid1_wl)
    #print(intid2_wl)
    #if len(intid1_wl.keys()) != len(intid2_wl.keys()):
    #    print("error, didn't find the same whitelists")
    #    return(0)

    rx1 = return_default_regex()[0]
    rx2 = return_default_regex()[1]
    
    f1 = gzip.open(fq1, 'rt') if fq1.endswith("fastq.gz") else open(fq1)
    f2 = gzip.open(fq2, 'rt') if fq2.endswith("fastq.gz") else open(fq2)

    if limit != None:
        print("Warning this is a subset")
    i1 = itertools.islice(f1, 0, limit, 1)
    i2 = itertools.islice(f2, 0, limit, 1)

    out_files_dic1 = {}
    out_files_dic2 = {}

    #print("intid1_wl, intid2_wl")
    #print(intid1_wl, intid2_wl)

    for ipair in valids.keys():#zip(intid1_wl.keys(), intid2_wl.keys()):
        out_files_dic1[ipair] = open(f"{outdir}/splitraw_{dname}_intid_{ipair}_R1.fastq".format(outdir, dname, ipair), 'w')
        out_files_dic2[ipair] = open(f"{outdir}/splitraw_{dname}_intid_{ipair}_R2.fastq".format(outdir, dname, ipair), 'w')

    id =[]
    seq = []
    plus = []
    qual = []
    chimeric = 0

    for i,v in enumerate(zip(i1,i2)):
            if i%100000 == 0:
                print("\r"+str(i), end='')
            l1 = v[0].strip()
            l2 = v[1].strip()

            # this corresponds to read id
            if i%4 == 0:
                id = [l1, l2]
            # this corresponds to read seq
            if i%4 == 1:
                seq = [l1, l2]
            # this corresponds to read '+'
            if i%4 == 2:
                plus = [l1, l2]
            # this corresponds to read qual
            if i%4 == 3:
                qual = [l1, l2]

            if qual:
                m1 = rx1.search(seq[0])
                m2 = rx2.search(seq[1])
                if m1 and m2:
                    # if both integration ids are found and valid, select the reads
                    #if m1.group('intid1') in intid1_wl.keys() and m2.group('intid2') in intid2_wl.keys():
                    i1 = m1.group('intid1')
                    i2 = m2.group('intid2')
                    key = i1+"_"+i2
                    if key in valids.keys():
                        out_files_dic1[key].write("\n".join([id[0], seq[0], plus[0], qual[0]]) +'\n')
                        out_files_dic2[key].write("\n".join([id[1], seq[1][len_sample_id-1:], plus[1], qual[1][len_sample_id-1:]]) +'\n')
                        #print(intid1, intid2, m1.group('intid1'), m2.group('intid2'))
                    else:
                        chimeric += 1
                            #if intid1 != intid2:
                            #print('mismatch between integration ids')
                            #raise AssertionError
                       # else:
                       # filter_good+=1
                    #else:
                    #    inv1.write("\n".join([id[0]+ " sample_id:%s"%(sample_id), seq[0], plus[0], qual[0]]) +'\n')
                    #    inv2.write("\n".join([id[1]+ " sample_id:%s"%(sample_id), seq[1], plus[1], qual[1]]) +'\n')
                    #    filter_bad+=1

                id = []
                seq = []
                plus = []
                qual = []
    total = i/4
    print(f"{fq1}: out of {total} reads, {chimeric} were chimeric ({chimeric/total:.2%})")
    # close the file handles
    for i in intid1_wl.values():
       out_files_dic1[i].close()
       out_files_dic2[i].close()


def return_default_regex(umi_len = 19, intid_len =8):
    rx1 = regex.compile("(?P<umi>.{"+str(umi_len)+"}).*(CTGTTCCTGTAGAAAT){e<=2}(?P<intid1>.{"+str(intid_len)+"}).*(CCGGACTCAGATCTCGA){e<=2}")
    rx2 = regex.compile(".*(GCTCATAGCGCTCG){e<=3}(?P<intid2>.{"+str(intid_len)+"}).*(CGACTCCATAGTCGAC){e<=2}")
    return(rx1, rx2)

def guess_valid_intid_pairs(outdir, f1_fn, f2_fn, sample_name, sampling_depth, umi_len = 19, intid_len = 8, rx1 = None, rx2 = None, orderm = 1.5):
    if rx1 == None:
        rx1 = return_default_regex()[0]
    if rx2 == None:
        rx2 = return_default_regex()[1]

    ## Trace the number and combination of valid integrations
    f1 = gzip.open(f1_fn, 'rt') if f1_fn.endswith("fastq.gz") else open(f1_fn)
    f2 = gzip.open(f2_fn, 'rt') if f2_fn.endswith("fastq.gz") else open(f2_fn)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print(f"dumping the integration pair table at the read level at:\n{outdir} using a stringency of base10({orderm})")
    random_start = random.randrange(1e4)
    random_start = random_start- random_start%4
    i1 = itertools.islice(f1, random_start+1, int(sampling_depth), 4)
    i2 = itertools.islice(f2, random_start+1, int(sampling_depth), 4)
    pair_tab = '{}/int_pair_table.txt'.format(outdir)
    o1 = open(pair_tab, 'w')
    for r1,r2 in zip(i1,i2):
        m1 = rx1.search(r1)
        m2 = rx2.search(r2)
        if m1 and m2:
            o1.write("\t".join([m1.group('umi'), m1.group('intid1'), m2.group('intid2')])+'\n')
    o1.close()
    f1.close()
    f2.close()
    pairs = open(pair_tab)
    
    pp = []
    for line in pairs:
        int1, int2 = line.strip().split('\t')[1:3]
        pp.append(int1+"_"+int2)
    
    pairs.close()

    x = pd.Series(pp).value_counts()
    xx = pd.Series([int(np.log10(i)) for i in x], index = x.index)
    
    valid_pair_path = "{}/valid_pairs".format(outdir)
    valid_pair_tab = open(valid_pair_path, 'w')
    if len(xx)==0:
        return(valid_pair_path)
    # keep only the most abundant pairs
    # orderm defines the number of orders of magnitude base 10 for which to filter out
    for i in x[xx>(xx[0]-orderm)].index:
        intid1, intid2 = i.split("_")
        #valid_pair_tab.write("\t".join([intid1, ogtk.UM.rev_comp(intid2)])+"\n")
        valid_pair_tab.write("\t".join([intid1, intid2])+"\n")
        print("\t".join([intid1, intid2]))
    valid_pair_tab.close()
    return(valid_pair_path)

def demult_fastq_pair(args):
    ''' wrapper for guessing and demultiplexing for pool.map()'''
    outdir, f1_fn, f2_fn, sample_name, sampling_depth, umi_len, intid_len, rx1, rx2, orderm = args
    print(f1_fn)

    valid_pairs_tab = guess_valid_intid_pairs(outdir = outdir, f1_fn = f1_fn, f2_fn = f2_fn, 
                            sample_name = sample_name, sampling_depth = sampling_depth,  
                            umi_len = umi_len, intid_len = intid_len, rx1 = rx1, rx2 = rx2, orderm = orderm)

    demult_outdir  = "/".join(valid_pairs_tab.split("/")[0:-1])
    log = f'[{time.strftime("%c")}]\nDemultiplexing {sample_name} {valid_pairs_tab} outdir = {demult_outdir}'

    print(log)
    split_fastq_by_intid(valid_pairs_tab = valid_pairs_tab, fq1 = f1_fn, fq2 = f2_fn, dname = sample_name, 
                    outdir = demult_outdir, limit = None, len_sample_id = 6)
     
def demult_per_illumina_index(rootdir, umi_len = 19, intid_len = 8, orderm = 1.5, sampling_depth = 1e5, threads = None):
    ''' Iterates over rootdir, finds R1 fastqs, sampling them {sampling_depth} to guess {orderm} the number of valid pairs - then proceeds to split the fastq files'''
    # TODO add support for gz
    import glob
    import time
                                           
    files_fq = glob.glob(rootdir+'/*R1*fastq*')
    print(f'[{time.strftime("%c")}]\nProcessing {rootdir} umi_len = {umi_len} intid_len = {intid_len}')
    rx1 = return_default_regex()[0]
    rx2 = return_default_regex()[1]
    
    if threads == None:
        threads = len(files_fq)
    pool = multiprocessing.Pool(threads)
    pool_cmds = []
    ## Trace the number and combination of valid integrations
    for f1_fn in files_fq:
        f2_fn = f1_fn.replace('R1', 'R2')
        sample_name = f1_fn.split("/")[-1].split(".")[0].split("_S")[0]
        outdir = f"{rootdir}/by_intid/{sample_name}"
        pool_cmds.append((outdir, f1_fn, f2_fn, sample_name, sampling_depth, umi_len, intid_len, rx1, rx2, orderm))

    pool.map(demult_fastq_pair, iter(pool_cmds) )
    pool.close()


def mlt_create_fasta_ref(ref_name, ref_seq, ref_path):
    fout = open(ref_path, 'w')
    fout.write('>{}\n{}\n'.format(ref_name, ref_seq))
    fout.close()

def mltbc_make_unique_fa_ref_entries(fa_ifn, fa_ofn, ref_name = 'hspdrv7_scgstl', wd='.'):
    '''
    Relabels reference entries in a FASTA file to create unique identifiers by appending the UMI (Unique Molecular Identifier).

    This function is designed to work on pair-wise sequence alignments in FASTA format. The function iterates over the input 
    FASTA file line by line, identifying header lines (starting with '>'). If a header line is identified, it modifies the 
    line by appending the UMI to create a unique reference name, before writing the line to the output FASTA file. 

    Parameters:
    fa_ifn: str
        Path to the input FASTA file. The file should contain pair-wise sequence alignments.

    fa_ofn: str
        Path where the output FASTA file will be saved. The output file will be a copy of the input file, but with modified 
        header lines for reference entries.

    ref_name: str, optional
        A string present in the header line of reference sequences. If a header line contains this string, the UMI will be 
        appended to the line. The default value is 'hspdrv7_scgstl'.

    Returns:
    None. The function writes to an output file and does not return any value.
    '''

    FF = open(fa_ifn)
    entries = []
    whole_fa = []
    for i,line in enumerate(FF):
        if line[0] == ">":
             entries.append(i)
        whole_fa.append(line)
    # if we don't set dtype to object, then the appended string is trimmed byt numpy
    whole_fa = np.array(whole_fa, dtype=object)
    start  =  1 if ref_name in whole_fa[0] else 0
    offset = -1 if ref_name in whole_fa[0] else 1
    umis = itertools.islice(entries, start, None, 2)
    refs = itertools.islice(entries, start+offset, None, 2)
    for umi,ref in zip(umis, refs):
        umi = whole_fa[umi][1:].strip()
        new_ref = "{}_{}\n".format(whole_fa[ref].strip(), umi)
        whole_fa[ref] = new_ref
    
    ofa = open(f'{wd}/{fa_ofn}', 'w')
    for i in whole_fa:
        ofa.write(i)
    ofa.close()
   

def mltbc_align_reads_to_ref(
        name, 
        fa_ofn,  
        ref_path, 
        input_fasta=None,
        ref_name = 'hspdrv7_scgstl', 
        mode = 'needleman', 
        gapopen = 20, 
        gapextend = 1, 
        rs=None,
        wd = '.',
        conda_prefix = None,
        verbose = False):
    '''
    Performs pair-wise alignment of reads against a reference, and relabels reference entries to make them unique.

    Supports a custom conda environment via the conda_prefix argument which should be a string of the form:
    ```conda run -n name_of_env```
    Assumes a conda environment with the EMBOSS suite installed:
        conda create -n emboss -y -c bioconda  emboss
    Parameters:
    name: str
        The base name for input and output files.

    fa_ofn: str
        Path where the output FASTA file with unique reference names will be saved.

    ref_path: str
        Path to the reference sequence in FASTA format.

    ref_name: str, optional
        A string present in the header line of reference sequences to be relabeled. Default is 'hspdrv7_scgstl'.

    mode: str, optional
        The alignment mode to be used. Options are 'needleman' and 'waterman'. Default is 'needleman'.

    gapopen: int, optional
        Penalty for opening a gap. Default is 20.

    gapextend: int, optional
        Penalty for extending a gap. Default is 1.

    rs: readset, optional
        A readset to writes out the consensus sequences as a FASTA file.
        Default is None.

    verbose: bool, optional
        If True, the function will print detailed messages during execution. Default is False.

    Returns:
    None. The function writes results to output files and does not return any value.

    The function first performs pair-wise alignment of reads against a reference using the 'needleall' or 'water' algorithm
    of the EMBOSS suite, depending on the 'mode' argument. The alignment is executed with subprocess.run, and the stdout and
    stderr outputs are captured and written to separate log files for each alignment format.

    The function then relabels the reference entries in the resulting alignment file by appending a unique identifier,
    using the helper function mltbc_make_unique_fa_ref_entries. The relabeled alignment is saved to the output FASTA file.
    '''

    if input_fasta is None:
        pwalg_in = f'{wd}/{name}_in_pair_algn.fa'
    else:
        pwalg_in = input_fasta
    pwalg_out = f'{wd}/{name}_out_pair_algn.fa'

    if rs is not None:
        rs.write_consensuses_fasta(pwalg_in)

    if verbose:
        print('Running pair-wise alignment')
    if mode == 'waterman':
        'water -gapextend 1 -gapopen 20 -datafile EDNAFULL -awidth3=100 -aformat3 fasta  -asequence ptol2-hspdrv7_scgstl.fa -bsequence msa_in.fa  -aaccshow3 yes -outfile ww.fasta -snucleotide2  -snucleotide1' 
    if mode == 'needleman':
        cmd_template = f'needleall -gapextend {gapextend} -gapopen {gapopen} -datafile EDNAFULL -awidth3=100  -minscore 90 -bsequence {ref_path} -asequence {pwalg_in} -aaccshow3 yes'

        if conda_prefix is not None:
            cmd_template = f"{conda_prefix} {cmd_template}"

        cmd_needleman = f'{cmd_template} -aformat3 fasta -outfile {pwalg_out}'
        # TODO change the three invokations for seqret instead ?
        # needleall -asequence seq1.fasta -bsequence seq2.fasta -outfile alignment.needle

        if verbose: print(cmd_needleman) 

        oo = subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)

        if verbose: print(oo.stdout, oo.stderr)

        with open(f"{pwalg_out}_fasta.olog", 'wb') as logo, open(f"{pwalg_out}_fasta.elog", 'wb') as loge:
            logo.write(oo.stdout)
            loge.write(oo.stderr)

        cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "score", pwalg_out.replace('fa', 'score'))

        if verbose: print(cmd_needleman) 
        oo = subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)

        with open(f"{pwalg_out}_score.olog", 'wb') as logo, open(f"{pwalg_out}_score.elog", 'wb') as loge:
            logo.write(oo.stdout)
            loge.write(oo.stderr)

        #cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "srs", pwalg_out.replace('fa', 'srs'))
        #print(cmd_needleman)
        #subprocess.run(cmd_needleman.split())

        cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "simple", pwalg_out.replace('fa', 'simple'))
        if verbose: print(cmd_needleman) 
        oo = subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)

        with open(f"{pwalg_out}_simple.olog", 'wb') as logo, open(f"{pwalg_out}_simple.elog", 'wb') as loge:
            logo.write(oo.stdout)
            loge.write(oo.stderr)
    
    mltbc_make_unique_fa_ref_entries(fa_ifn = pwalg_out, fa_ofn = fa_ofn, ref_name = ref_name, wd=wd)


def mlt_compute_barcode_matrix_paired(fa_ifn, fa_ifn2, umi_whitelist, umi_blacklist, tab_out, bint_db_ifn, tab_full_out = None, do_rc = False):

    '''Writes to disk a tabulated file with the corresponfing bint strings.
    Requires a list of valid UMIs (whitelist:union set from R1 and R2), a list of
    invalid UMIs (sometimes the same UMI cannot be merged and thus appears on both,
    merged (gets priority and defines blacklist) and unmerged fractions,  and a
    fasta file that provides the pairwise aligment between a lineage allele and the
    reference'''

    FF = pyfaidx.Fasta(fa_ifn, as_raw=True)
    FF2 = pyfaidx.Fasta(fa_ifn2, as_raw=True)
    
    #TODO: make more elegant by subtracting intersection of white and blacklists
    fa_entries = umi_whitelist#[i for i in FF.keys()]
    outf = open(tab_out, 'w')
    if tab_full_out == None:
        tab_full_out = tab_out.replace('.txt', '_full.txt')
    full_out = open(tab_full_out, 'w')

    for i in range(0, int(len(fa_entries)/2), 2):
        #if fa_entries[i] in umi_whitelist and fa_entries[i] not in umi_blacklist:
    #TODO: make more elegant by subtracting intersection of white and blacklists
        if fa_entries[i] not in umi_blacklist:
            ref_seq  = FF[i+1][:].upper() #if not do_rc else ogtk.UM.rev_comp(FF[i+1][:].upper()) 
            read_seq = FF[i][:].upper() #if not do_rc else ogtk.UM.rev_comp(FF[i][:].upper()) 

            ref_seq2  = FF2[i+1][:].upper() if not do_rc else ogtk.UM.rev_comp(FF2[i+1][:].upper()) 
            read_seq2 = FF2[i][:].upper() if not do_rc else ogtk.UM.rev_comp(FF2[i][:].upper()) 

            lineage_vector = get_lineage_vector((ref_seq, read_seq, bint_db_ifn, 'read1'))
            lineage_vector2 = get_lineage_vector((ref_seq2, read_seq2, bint_db_ifn, 'read2'))
            if lineage_vector != "trash" and lineage_vector2 != "trash":
                #full_out.write('{}\n'.format('\t'.join(["ref",  fa_entries[i], ref_seq, flineage_vector])))
                #full_out.write('{}\n'.format('\t'.join(["read", fa_entries[i], read_seq, flineage_vector])))
                flineage_vector = '\t'.join(lineage_vector)
                flineage_vector2 = '\t'.join(lineage_vector2)


                full_out.write('{}\n'.format('\t'.join(["ref",  fa_entries[i], ref_seq, flineage_vector])))
                full_out.write('{}\n'.format('\t'.join(["read1",  fa_entries[i], read_seq, flineage_vector])))
                full_out.write('{}\n'.format('\t'.join(["read2",  fa_entries[i], read_seq, flineage_vector2])))

                #full_out.write('{}\n'.format(fa_entries[i]))
                #full_out.write('{}\n'.format('\n'.join([ref_seq, read_seq])))
                #full_out.write('{}\n'.format('\n'.join([read_seq2, ref_seq2])))
                #full_out.write('\t'.join(lineage_vector)+'\n')
                #full_out.write('\t'.join(lineage_vector2)+'\n')
            
                outf.write('{}\t{}\t{}\n'.format(fa_entries[i], '\t'.join(lineage_vector), '\t'.join(lineage_vector2)))
    outf.close()
    full_out.close()

###########
########### END
###########
    
  
def align_reads_to_ref(name, fa_ofn, consensus_dict, ref_path, ref_name = 'hspdrv7_scgstl', mode='needleman', gapopen = 20, gapextend = 1, verbose = False):
    pwalg_in = '{}_in_pair_algn.fa'.format(name)
    pwalg_out = '{}_out_pair_algn.fa'.format(name)
    outf = open(pwalg_in, 'w')
    for k,v in zip(consensus_dict.keys(), consensus_dict.values()):
        outf.write(">{}\n{}\n".format(k,v))
    outf.close()

    if mode == 'waterman':
        'water -gapextend 1 -gapopen 20 -datafile EDNAFULL -awidth3=100 -aformat3 fasta  -asequence ptol2-hspdrv7_scgstl.fa -bsequence msa_in.fa  -aaccshow3 yes -outfile ww.fasta -snucleotide2  -snucleotide1' 
    if mode == 'needleman':
        cmd_template = 'needleall -gapextend {} -gapopen {} -datafile EDNAFULL -awidth3=100  -minscore 90 -bsequence {} -asequence {} -aaccshow3 yes'.format(gapextend, gapopen, ref_path, pwalg_in)
    
        cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "fasta", pwalg_out)
        if verbose: print(cmd_needleman) 
        oo= subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)

        cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "score", pwalg_out.replace('fa', 'score'))
        if verbose: print(cmd_needleman) 
        oo= subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)

        #cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "srs", pwalg_out.replace('fa', 'srs'))
        #print(cmd_needleman)
        #subprocess.run(cmd_needleman.split())

        cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "simple", pwalg_out.replace('fa', 'simple'))
        if verbose: print(cmd_needleman) 
        oo= subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)
    
    make_unique_fa_ref_entries(fa_ifn = pwalg_out, fa_ofn = fa_ofn, ref_name = ref_name)

def return_alignment_tuples(fasta_ifn):
    ''' Returns a tuple of tuples in the form: 
            ( (query_name, ref_name), (query_seq,ref_seq)) 
        from an alignment from a fasta file
    '''
    FF = pyfaidx.Fasta(fasta_ifn, as_raw=True)
    fa_entries = [i for i in FF.keys()]
    return (((fa_entries[i+1], fa_entries[i]), (FF[i+1][:].upper(), FF[i][:].upper())) for i in range(0, len(fa_entries), 2))


def compute_barcode_matrix_merged(fa_ifn, tab_out,  bint_db_ifn, tab_full_out = None, do_rc = False):
    '''
    Generates a barcode matrix by performing a pairwise alignment between a lineage allele and a reference sequence.
    The output is written to disk as a tabulated file with corresponding bint strings.

    Parameters:
    fa_ifn: str
        The input FASTA file providing the pairwise alignment between a lineage allele and the reference.

    tab_out: str
        The output file path where the tabulated file with the corresponding bint strings will be saved.

    bint_db_ifn: str
        The input file providing the binary identification (bint) strings.

    tab_full_out: str, optional
        The output file path for the full tabulated data. If not provided, '_full.txt' will be appended to the name of tab_out file.

    do_rc: bool, optional
        If True, the function performs reverse complement of the sequences. Default is False.

    The function first reads the input FASTA file and initializes output files. For each pair of sequences (reference and read) 
    in the FASTA file, it computes a lineage vector using the get_lineage_vector function. 

    If the lineage vector does not contain "trash", the function writes the vector to the output files and keeps track of the 
    number of "non-trash" hits. If the lineage vector contains "trash", the function records it but does not write it to the output.

    Finally, if there are more than one "trash" vectors, the function prints their counts and proportions, and the proportion of "non-trash" hits.
    '''

    FF = pyfaidx.Fasta(fa_ifn, as_raw=True)
    fa_entries = [i for i in FF.keys()]
    outf = open(tab_out, 'w')
    if tab_full_out is None:
        tab_full_out = tab_out.replace('.txt', "_full.txt")
    full_out = open(tab_full_out, 'w')

    hits_trash = []
    trashed = []
    for i in range(0, len(fa_entries), 2):
            ref_seq  = FF[i+1][:].upper() if not do_rc else ogtk.UM.rev_comp(FF[i+1][:].upper()) 
            read_seq = FF[i][:].upper() if not do_rc else ogtk.UM.rev_comp(FF[i][:].upper()) 


            lineage_vector = get_lineage_vector((ref_seq, read_seq, bint_db_ifn, 'merged'))

            if "trash" not in lineage_vector:
                hits_trash.append(1)
                flineage_vector = '\t'.join(lineage_vector)

                outf.write('{}\t{}\n'.format(fa_entries[i], flineage_vector))
                full_out.write('{}\n'.format('\t'.join(["ref",  fa_entries[i], ref_seq, flineage_vector])))
                full_out.write('{}\n'.format('\t'.join(["read", fa_entries[i], read_seq, flineage_vector])))
            else:
                trashed.append(lineage_vector)
                hits_trash.append(0)

    if len(trashed) >1:
        trashed = pd.Series(trashed)
        print(trashed.value_counts(normalize=True))
        print(trashed.value_counts(normalize=False))
        print(f'not trashed {sum(hits_trash)}/{len(hits_trash)} {sum(hits_trash)/len(hits_trash):.2f}')
    full_out.close()
    outf.close()

def make_unique_fa_ref_entries(fa_ifn, fa_ofn, ref_name = 'hspdrv7_scgstl'):
    '''Small helper function that relabels reference entries on a pair-wise sequence alignments.
        It appends the UMI to the ref name'''
    FF = open(fa_ifn)
    entries = []
    whole_fa = []
    for i,line in enumerate(FF):
        if line[0] == ">":
             entries.append(i)
        whole_fa.append(line)
    # if we don't set dtype to object, then the appended string is trimmed byt numpy
    whole_fa = np.array(whole_fa, dtype=object)
    start  =  1 if ref_name in whole_fa[0] else 0
    offset = -1 if ref_name in whole_fa[0] else 1
    umis = itertools.islice(entries, start, None, 2)
    refs = itertools.islice(entries, start+offset, None, 2)
    for umi,ref in zip(umis, refs):
        umi = whole_fa[umi][1:].strip()
        new_ref = "{}_{}\n".format(whole_fa[ref].strip(), umi)
        whole_fa[ref] = new_ref
    
    ofa = open(fa_ofn, 'w')
    for i in whole_fa:
        ofa.write(i)
    ofa.close()
   
def get_lineage_vector(args):
    ''' with a given reference and red sequence pair, returns a lineage vector making use of a barcode interval db (bintdb)
    Alleles are encoded in the format [mismatches, insertions, deletions] e.g. "mm.ins.dele" 
    '''
    if len(args) <5:
        ref_seq, read_seq, bint_db_ifn, read_end, = args
        debug = False
    if len(args) == 5:
        ref_seq, read_seq, bint_db_ifn, read_end, debug = args
    # mismatch filtering threshold
    max_mm = 50 
    bint_db = pyaml.yaml.load(open(bint_db_ifn), Loader=pyaml.yaml.FullLoader)
    

    last_bint_len = np.diff(list(bint_db['intervs'].values())[-1])[0]

    ins = 0
    weird = 0
    mis = 0
    de = 0
    bint_dic = {}
    bint_keys = [i for i in bint_db['intervs'].keys()]

    # find if trimming the read is possible and do it if it is
    anchor1 = regex.compile("({}){{e<=3}}".format(bint_db['anchor_plasmid'].upper()))
    anchor2 = regex.compile("({}){{e<=3}}".format(bint_db['anchor_plasmid2'].upper()))
    
    match1 = anchor1.search(read_seq)
    match2 = anchor2.search(read_seq)

    ## TODO How to filter out initial strage integrations derived from the transcriptome?
    if match1:
        #trim
        read_seq = read_seq[match1.span()[0]:]
        ref_seq = ref_seq[match1.span()[0]:]

    if not (match1 or match2):
        return('trash-unmatched')
    
    for i,bint in enumerate(bint_keys):
        bint_dic[bint] = Bint(name=bint, index=i)

    # Since ref and read seqs have been aligned we can iterate through them in parallel
    # for each operation found (mis, ins, del) create a record ["type", position, "character"]
    for i, (ref_chr, read_chr) in enumerate(zip(ref_seq, read_seq)):
        ref_chr = ref_chr.upper()
        read_chr = read_chr.upper()
        ref_chr = ref_chr.replace('N', '-') 
        read_chr = read_chr.replace('N', '-') 
        cur_bint = return_bint_index(bint_l=bint_db['intervs'], current_position = i, last_bint_len = last_bint_len)
        if cur_bint != 0: # 0 means that it couldn't find an overlapping bint
            if ref_chr == "-" and read_chr == "-":
                weird+=1
            else:
                ########## insertions
                if ref_chr == "-":
                    ins+=1
                    # update the bint db by pushing down the coords by 1
                    for j in range(bint_keys.index(cur_bint), len(bint_keys)):
                        next_bint = bint_keys[j]
                        bint_db['intervs'][next_bint][0]+= 1
                    # save integration
                    bint_dic[cur_bint].ins.append(['i', i, read_chr])
                ########## deletions
                if read_chr == "-":
                    de+=1
                    bint_dic[cur_bint].dele.append(['d', i, ref_chr])
                ########## mismatches
                if read_chr != ref_chr and ref_chr != '-' and read_chr != '-':
                    mis+=1
                    bint_dic[cur_bint].mis.append([ref_chr, i, read_chr])
    silent = True
    if silent == False:
        print("ins:", ins, "del:", de, "weird:", weird, "mis:", mis)

    # At this point all Bint objects inside the dictionary should have been filled up
    # Start and end class atributes intialization 
    # - update coordinate data inside the class
    #   the start coordinate is specified using the recorded 'intervs' list using each ith entry
    #   the end coordinate is specified using the recorded 'intervs' list using i+1th for ith
    for bint in bint_dic.keys():
        bint_dic[bint].start = bint_db['intervs'][bint][0]
        if bint_keys.index(bint) == len(bint_keys)-1:
            # last case
            bint_dic[bint].end = bint_dic[bint].start + (last_bint_len -1)
        else:
            next_bint = bint_keys[bint_keys.index(bint) + 1]
            bint_dic[bint].end = bint_db['intervs'][next_bint][0]-1

    # Update coordinates of operations (ins, del, mis) relative to bint.start and bint.end:
    for bint in bint_dic.values():
        bint.transfer_coords_to_operation_records()

    # TODO create class that stores Bints and internalize this functionality 
    # when sequencing was paired end and reads couldn't be merged
    ## discard any non covered bint with the appropiate flag 
    ## considering that insertions have already been accounted for by the updated coords

    if read_end == 'read1' or read_end == 'read2':
        # otherwise read_end == 'merged'
        bint_keys = [i for i in bint_dic.keys()]
        critical_pos = bint_db['read_length'][read_end] if read_end == 'read1' else len(ref_seq)-bint_db['read_length'][read_end]
        discard_index = bint_keys.index( return_bint_index(bint_db['intervs'], critical_pos, last_bint_len))#TODO change previus line to something more elegant
        discarded_range = range(discard_index, len(bint_dic.keys())) if read_end == 'read1' else range(0, discard_index+1)
        if 'read2' == read_end and False: # TODO remove debug
            print(ref_seq)
            print(read_seq)
            print(discard_index)
            print([bint_keys[i] for i in discarded_range])
            pp()
        for i in discarded_range:
            discard_bint = bint_keys[i]
            bint_dic[discard_bint].dele = []#[('NA','0', 'NA')] 
            bint_dic[discard_bint].mis  = []#[('NA','0', 'NA')] 
            bint_dic[discard_bint].ins = [('NA')] 
    if debug:
        print("First")
        print([bint.generate_barcode_string() for bint in bint_dic.values()])
   
    for bint in bint_dic.values():
    # flag as NA bints that surpass a mismatch threshold
        if len(bint.mis) > max_mm:
            bint.dele = []#[('NA','0', 'NA')] 
            bint.mis  = []#[('NA','0', 'NA')] 
            bint.ins = [('NA')] 
            # filter out clear excisions/NA
        if len(bint.dele)-1 == bint.end - bint.start:
            bint.dele = [('e','x','x')] 
    if debug:
        print("flagged as NA bints that surpass a mismatch threshold")
        print([bint.generate_barcode_string() for bint in bint_dic.values()])

    # return "trash" if homopolymers are found on any bint
    bint_strings = [bint.generate_barcode_string() for bint in bint_dic.values()] 

    if debug:
        print("filtered out clear excisions/NA")
        print(bint_strings)

    if any(["trash" in i for i in bint_strings]):
        return("trash-gen") 

    ## TODO How to get rid of uncovered bints that look like excisions?
    ## the logic is to assign as NAs any stretch of excisions for which we have no evidence of it's end, e.g. it reaches the last bint 
    if not match2:
        # this means that we have no evidence for a real excision, therefore turn into NAs the last observed excisions AND
        # the next to last barcode
        corrected = False
        for i in range(len(bint_strings))[::-1]:
            if bint_strings[i] == '..exx':
                bint_strings[i] = '.NA.'
            if bint_strings[i] == '..':
                corrected = True
            elif not corrected:
                bint_strings[i] = '.NA.'
                corrected = True  
    if debug:
        print(bint_strings)
        print(match1, match2)
    
    return(bint_strings)

def return_bint_index(bint_l, current_position, last_bint_len):
    '''Helper function that returns the key for the bint (barcode interval) that covers a given position on a pairwise alignment'''
    bint_index = None
    keys = [i for i in bint_l.keys()]
    for i in range(len(keys)):
        key  = keys[i]
        if i == (len(keys)-1):
            if current_position >= bint_l[keys[i]][0] and current_position <= bint_l[keys[i]][0] + last_bint_len:
                return(key)
        else:
            if current_position >= bint_l[keys[i]][0] and current_position <=bint_l[keys[i+1]][0]:
                bint_index = key
                return(bint_index)
    return(0) # means that it couldn't find an overlapping bint

def consecutive_to_string(vec):
    '''Returns a concatenated string for consecutive elements 
    regarding the site-specific substitution.
    For example it generates barcode chains; a string like "50ACT52.60C" out of [50A, 51C, 52T, 60C]
    '''
    #vec = [ [type, pos, char], ['i', '10', 'A'] ]
    if len(vec)==1:
        return(''.join(vec[0]))
    
    hits = []
    # TODO check why we get stings and not integers for i[1]
    positions = [int(i[1]) for i in vec]
    nucleotides = [i[2] for i in vec]
    bc_chains = [] 
    for i in range(len(positions)-1):
        delta = positions[i+1] - positions[i]
        if delta > 1:
            hits.append(i)
            bc_chain = str(positions[hits[0]]) + ''.join(nucleotides[hits[0]:(1+hits[-1])]) + str(positions[hits[-1]])
            bc_chains.append(bc_chain)
            hits = []
        elif delta == 1:
            hits.append(i)
        elif delta < 1:
            raise ValueError("downstream coordinates must be bigger than upstream ones")
        # we are now eval the next to last 
        if i == (len(positions)-2):
            if delta == 1:
                hits.append(i+1)
                bc_chain = str(positions[hits[0]]) + ''.join(nucleotides[hits[0]:(1+hits[-1])]) + str(positions[hits[-1]])
                bc_chains.append(bc_chain)
            elif delta > 1:
                hits.append(i)
                bc_chain = str(positions[hits[0]]) + ''.join(nucleotides[hits[0]:(1+hits[-1])]) + str(positions[hits[-1]])
                bc_chains.append(bc_chain)

                hits.append(i+1)
                bc_chain = str(positions[i+1]) + ''.join(nucleotides[i+1]) + str(positions[i+1])
                bc_chains.append(bc_chain)
                #chain = []
                #chain.append(x[i+1])
                #print('break', chain[0], chain[-1])
    #print('break', chain[0], chain[-1])
    polyA = regex.compile("AAAAAAAA")
    polyG = regex.compile("GGGGGGGG")
    for i in range(len(bc_chains)):
        if polyA.search(bc_chains[i]) or polyG.search(bc_chains[i]):
            return("trash-poly")
    return(''.join(bc_chains))    


class Bint():
    def __init__(self, name, index):
        self.name = name
        self.mis = []
        self.dele = []
        self.ins = []
        self.index = index

        # start and end are defined once a locus has been screened and they store the real start end of the barcode interval
        self.start = None
        self.end = None
    def print_bc_vectors(self):
        print(self.name)
        print(self.mis)
        print(self.dele)
        print(self.ins)

    def transfer_coords_to_operation_records(self):
        """ transfer the coordinates from .start and .end into the raw coordintaes of each operation """
        # mismatches
        for i in range(len(self.mis)):
            self.mis[i][1] = str(self.mis[i][1] - self.start)
        # deletions
        for i in range(len(self.dele)):
            self.dele[i][1] = str(self.dele[i][1] - self.start)
        # insertions
        for i in range(len(self.ins)):
            self.ins[i][1]= str(self.ins[i][1] - self.start)

    def generate_barcode_string(self):
        ''' '''
        if self.start == None:
            raise ValueError('Uncorrected coordinates for indels. Update bints first and transfer the correct coordinate to the operations')
        mismatches = self.gen_mismatch_string()
        #mismatches = ''
        insertions = self.gen_insertion_string()
        deletions =  self.gen_deletion_string()
        allele_string = ".".join([mismatches, insertions, deletions])
        return(allele_string)
        #TODO return the bint index as part of the barcode
        #return(f"b{self.index}@" + allele_string)

    def gen_mismatch_string(self):
        if len(self.mis)>1:
            return(''.join([''.join(i) for i in self.mis]))
        else:
            return("")

    def gen_insertion_string(self):
        if len(self.ins) >0:
            return(consecutive_to_string(self.ins))
        else:
            return("")

    def gen_deletion_string(self):
        if len(self.dele) >0:
            return(consecutive_to_string(self.dele))
        else:
            return("")
# ========= end of Bint

def create_fasta_ref(ref_name, ref_seq, ref_path):
    fout = open(ref_path, 'w')
    fout.write('>{}\n{}\n'.format(ref_name, ref_seq))
    fout.close()

  
def compute_dist_mat(mat, bcs, bias, cores = 50, dont_include = 3):
    cores = int(cores)
    #print(f"Warning: FORCED for  using {cores} cores")
    print(f'got {len(bcs)} barcodes')
    #print(np.__config__.show())
    if cores == 0:
        print("zero cores")
        cc = similarity_score((mat, bcs, bias, dont_include)) 
        return(cc)
    else:
        print("pooling cores")
        chunks = np.array_split(bcs, cores)
        it = iter([(mat, chunk, bias, dont_include) for chunk in chunks])

        print(f"Firing up {cores} workers")
        pool = multiprocessing.Pool(int(cores))
        cc = pool.map(similarity_score, it)
        pool.close()
        pool.join()
        print("")
        cum_mat = cc[0]
        for i in cc[1:]:
            cum_mat = cum_mat + i
        return(cum_mat)

def similarity_score(args):
    from scipy import sparse
    mat, bcs, bias, dont_include = args
    first = True
    bias = np.array([float(i) for i in bias])

    pct_i = int(len(bcs)/100)
    if pct_i == 0:
        pct_i = 1e-6

    for i,bc in enumerate(bcs):
        if bc > int(dont_include):
            # TODO fully understand the transition between the matrix types for the following three lines
            mm = mat == bc
            mm = sparse.csc_matrix(mm *bias)
            ee = mm.dot(mm.transpose())
            if first:
                first = False
                cum_mat = ee
            else:
                if i%pct_i==0:
                    print(f'{100*i/len(bcs):.2f}%', end = '\r')
                cum_mat = cum_mat + (ee.todense()*1)
    return(cum_mat)

def to_numpy(mat, ofn):
    import pickle
    with open(ofn, 'wb') as out:
        pickle.dump(mat, out)


def cell_dominant_allele(rs, cell_umi_dict, allele_reps_fn = None):
    '''
    returns [and writes] a cell-level table of dominant seqs and their stats
    '''
    cell_alleles = [["cell", "dominant_seq", "dominant_reads", "reads_per_umi", "pool_size"]]
    for cell in cell_umi_dict.keys():
        for umi in cell_umi_dict[cell]: 
            umol = rs.umis[umi]
            if umol.consensus is not None:
                cons = umol.consensus
                entry = [cell, cons['dominant_seq'], cons['dominant_reads'], cons['reads_per_umi'], cons['pool_size']]
                cell_alleles.append(entry)
    cell_alleles = pd.DataFrame(cell_alleles[1:], columns=cell_alleles[0])   
    if allele_reps_fn is not None:
        cell_alleles.to_pickle(allele_reps_fn)
        print(f'saved cell_alleles {allele_reps_fn}')
    return(cell_alleles)

def get_rep_allele(x, min_mol = 3, min_domix = 0.75):
    '''
    Helper function for cell consensus calling.
    To determine the allele present on a given cell, one must select from a pool of candidate sequences (each belonging to a prefiltered UMI)
    '''
    top_seqs = x['dominant_seq'].value_counts()
    domix = x['dominant_seq'].value_counts(normalize = True).to_list()[0]

    if len(top_seqs) >= min_mol and domix >= min_domix:
        return(top_seqs.index[0])
    else:
        return(False)

def conf_dict_write(conf, level = None, **kargs):
    if level is not None:  
        if level not in conf.keys():
            conf[level] = {}
        for k,v in kargs.items():
            conf[level][k] = v
        return(conf)
    else:
        for k,v in kargs.items():
            conf[k] = v
        return(conf)
def conf_dict_write_iter(conf, iterable):
    print('memelo')

def store_molecule_saturation_stats(conf, rs):
    '''
    Helper function to record the saturation stats.
    '''
    # TODO improve the objectives of this routine
    conf = conf_dict_write(
            conf, 
            level='mols', 
            nreads = {'unmerged':rs.nreads},
            desc = "number of umis whose rank1 sequence is a >= [1, 2, 3, 4, 5, 10, 100, 1000]",
            saturation = {'unmerged':{'uncorrected':rs.allele_saturation()}}
            )
    return(conf)

def create_workdir(outdir, config_card_dir):
    import os
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
        print("Created {}".format(outdir))

    if not os.path.isdir(config_card_dir):
        os.makedirs(config_card_dir, exist_ok=True)
        print("Created {}".format(config_card_dir))
        
def return_fq_fields():
   return ['readid', 'seq', 'plus', 'qual']

def load_fastq_to_polars(ifn):
    ''' Helper function to load an uncompressed fastq file and tabulates into polars
    Not lazy
    '''
    fq_fields = return_fq_fields()
    return (
         pl.scan_csv(ifn, separator='.', has_header=False)
        .with_columns((pl.arange(0, pl.count()) // 4)
        .alias("step"))
        .groupby("step", maintain_order=True)
        .agg([  pl.col("column_1").take(i).alias(name) for i, name in enumerate(fq_fields) ])
        .drop('step')
        .drop('plus')
        #.collect(streaming=True)
    )


def split_fastq_by_intid_pl(
    outdir, 
    ifn,
    anchor1,
    anchor2,
    min_cov = 0.1,
    dont_filter = False,
    ):
    '''

    For Nora's line these are the recommended values for the anchors
    anchor1 = 'GTGGCAGG'
    anchor2 = 'TCCTGTAG'
    '''
    ifn2 = ifn.replace('R1', 'R2')

    fq_fields = return_fq_fields()
    fields_r2 = {k:k+'2' for k in fq_fields}
    del fields_r2['plus']

    res = (pl.concat(
            [
                load_fastq_to_polars(ifn),
                load_fastq_to_polars(ifn2).rename(fields_r2)
             ],
                how="horizontal",
                rechunk=True,
            )
            .with_columns(
                [pl.col('seq').str.extract(f"{anchor1}",1).alias('intid1'),
                 pl.col('seq2').str.extract(f"{anchor2}",1).alias('intid2')
                ]
            )
            .with_columns(pl.count().over(['intid1', 'intid2']).alias('counts'))
            .with_columns((pl.col('counts')/pl.col('counts').max()).alias('ncounts'))
           )

    if dont_filter:
        return res

    res = (res
            .filter(pl.col('ncounts')>=min_cov)
    )

    valid_groups = res.select(['intid1', 'intid2']).unique()

    for g1,g2 in zip(valid_groups['intid1'], valid_groups['intid2']):
        print(f'{g1}, {g2}')
        (res
            .filter((pl.col('intid1')==g1) & (pl.col('intid2')==g2))
            .with_columns(pl.lit('+').alias('plus'))
            .select(['readid', 'seq', 'plus', 'qual'])
            .write_csv(f'{outdir}/{g1}-{g2}_R1.fastq', separator='\n', has_header=False)
         )

        (res
        .filter((pl.col('intid1')==g1) & (pl.col('intid2')==g2))
        .with_columns(pl.lit('+').alias('plus'))
        .select(['readid2', 'seq2', 'plus', 'qual2'])
        .write_csv(f'{outdir}/{g1}-{g2}_R2.fastq', separator='\n', has_header=False)
     )

