#import ogtk
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
import ogtk.bbin_alleles as bbin
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

def mltbc_make_unique_fa_ref_entries(fa_ifn, fa_ofn, ref_name = 'hspdrv7_scgstl'):
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
   

def mltbc_align_reads_to_ref(name, fa_ofn, rs, ref_path, 
    ref_name = 'hspdrv7_scgstl', 
    mode='needleman', 
    gapopen = 20, 
    gapextend = 1, 
    verbose = False):
    pwalg_in = '{}_in_pair_algn.fa'.format(name)
    pwalg_out = '{}_out_pair_algn.fa'.format(name)

    rs.write_consensuses_fasta(pwalg_in)

    if verbose:
        print('Running pair-wise alignment')
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
    
    mltbc_make_unique_fa_ref_entries(fa_ifn = pwalg_out, fa_ofn = fa_ofn, ref_name = ref_name)


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
    
def bulk_bin_alleles(name, intid, intid2_R2_strand, 
    config_card_dir, 
    outdir, 
    fq1, 
    fq2, 
    bint_db_ifn, 
    min_reads_per_umi = 4, 
    threads = 100, 
    end = 5000, 
    ranked_min_cov = 5, 
    umi_len = 17, 
    umi_errors = 1,
    trimming_pattern = None, 
    alg_gapopen = 20, 
    alg_gapextend =1,
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

        
    # TODO clean this mess redirectin to who knows where
    # TODO do we neeed thiss?
    if intid2_R2_strand != None:
        intid2 = ogtk.UM.rev_comp(intid2_R2_strand)# same strand as initd1
    else:
        intid2 = len(intid) * "N"
   
    ####bint_db_ifn = '/local/users/polivar/src/projects/mltracer/conf/gstlt_barcode_intervsx10.yaml'
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
        Xmx = '2g',
        modality = 'ecct', # error-correct with Tadpole.sh
        log = f'{out_prefix}_merge_stats.txt',
    )
    cmd_bbmerge = f"bbmerge-auto.sh in1={fq1} in2={fq2} out={conf['bbmerge']['fastq_merged']} \
        k={conf['bbmerge']['k']} \
        outu={conf['bbmerge']['fastq_umerged1']} \
        outu2={conf['bbmerge']['fastq_umerged2']} \
        ihist={conf['bbmerge']['ihist']} \
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
  


####    yaml_out_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
####    yaml_out = {'name':name, 'desc':{}, 'intid1':intid, 'intid2':intid2, 'modality':'bulk'}
    
    # for intid in dataset
    # merge
####    fqm =       outdir + '/{}_merged.fastq'.format(name)
####    fqum1 =     outdir + '/{}_unmerged_R1.fastq'.format(name)
####    fqum2 =     outdir + '/{}_unmerged_R2.fastq'.format(name)
####    ihist =     outdir + '/{}_ihist.txt'.format(name)
####    log_merge = outdir + '/{}_merge_stats.txt'.format(name)
####    cmd_merge = "bbmerge-auto.sh in1={} in2={} out={} outu={} outu2={} ihist={} ecct extend2=20 iterations=5 interleaved=false -Xmx2g ".format(fq1, fq2, fqm, fqum1, fqum2, ihist)
    
####    yaml_out['desc']['read_merge'] = "Parameters used for the merging of read 1 and read 2"
####
####    umi_counts_ofn = outdir + '/{}_umi_counts.txt'.format(name)
####    yaml_out['mols'] = {}
####    yaml_out['mols']['umi_len'] = umi_len
####    yaml_out['mols']['umi_correction'] = umi_errors
####    yaml_out['mols']['counts'] = umi_counts_ofn 
####    
####    yaml_out['read_merge'] = {}
####    yaml_out['read_merge']['fqm'] = fqm
####    yaml_out['read_merge']['fqum1'] = fqum1
####    yaml_out['read_merge']['fqum2'] = fqum2
####    yaml_out['read_merge']['ihist'] = ihist
####    yaml_out['read_merge']['log_merge'] = log_merge
####
    # align
####    fa_correctedm = outdir+'/{}_corrected_merged.fa'.format(name)
####    fa_corrected1 = outdir+'/{}_corrected_R1.fa'.format(name)
####    fa_corrected2 = outdir+'/{}_corrected_R2.fa'.format(name)
####
####    yaml_out['desc']['alignment'] = "Parameters for the pairwise alignment of every recovered allele and the reference. The corrected files fix the issue of duplicated fasta entry names"
####
####    yaml_out['alignment'] = {}
####    yaml_out['alignment']['fa_correctedm'] = fa_correctedm
####    yaml_out['alignment']['fa_corrected1'] = fa_corrected1
####    yaml_out['alignment']['fa_corrected2'] = fa_corrected2
####
####    # more outfiles
####    merged_tab_out =   outdir + '/{}_barcodes_merged.txt'.format(name)
#### 
####    yaml_out['desc']['lineage'] = "Main output of the pre-processing scripts. Tabulated files that contain the bint barcode string"
####    yaml_out['lineage'] = {}
####    yaml_out['lineage']['merged_tab'] = merged_tab_out
####    yaml_out['lineage']['merged_full'] = merged_tab_out.replace('.txt','_full.txt')
####    yaml_out['alignment']['alg_mode'] = alg_mode
####    yaml_out['alignment']['alg_gapopen'] = alg_gapopen
####    yaml_out['alignment']['alg_gapextend'] = alg_gapextend


####

    # merge paired-end reads if possible
    log_merge =conf['bbmerge']['log'] 
    pp = subprocess.run(cmd_bbmerge.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(f"merge log\n{log_merge}")
    with open(log_merge, 'wb') as fout:
        fout.write(pp.stdout)
        fout.write(pp.stderr)
    
    # >>>
    _bulk_bin_alleles(conf_fn, conf, use_cache = use_cache, threads = threads)
    return(None)

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

    if end != None:
        print("Warning: end is not None; You are not processing all the data", end)
    rs_paths = [pickled_merged_readset, pickled_readset1, pickled_readset2]

    # do a pass on the raw fastqs and group reads by UMI
    if all([os.path.exists(i) for i in rs_paths]) and use_cache:
        rsm = pickle.load(open(pickled_merged_readset, 'rb'))
        rs1 = pickle.load(open(pickled_readset1, 'rb'))
        rs2 = pickle.load(open(pickled_readset2, 'rb'))
        print(f"loaded cached ReadSets: {'_'.join(rs_paths)}")
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
            end = end)

    # remove umis that exist on the merged set out of the unmerged
    unmerged_by_error = set(rsm.umi_list()).intersection(set(rs1.umi_list()))
    rs1.delete(unmerged_by_error)
    rs2.delete(unmerged_by_error)

    #TODO the saturation stats should be outsourced somehow and reported for all three cases not just the merged (rsm) 
    ### TODO for the refactored version too
    ###yaml_out['mols']['saturation'] = {}
    ###yaml_out['mols']['saturation']['merged'] = {}
    ###yaml_out['mols']['saturation']['unmerged'] = {}
    ###yaml_out['mols']['nreads'] = {}
    ###yaml_out['mols']['nreads']['total'] = sum([rsm.nreads, rs1.nreads])
    ###yaml_out['mols']['nreads']['merged'] = rsm.nreads
    ###yaml_out['mols']['nreads']['umerged'] = rs1.nreads
    ###yaml_out['mols']['desc'] = "number of umis whose rank1 sequence is a >= [1, 2, 3, 4, 5, 10, 100, 1000]"
    ###yaml_out['mols']['saturation']['merged']['uncorrected'] = rsm.allele_saturation()
    ###yaml_out['mols']['saturation']['unmerged']['uncorrected1'] = rs1.allele_saturation()
    ###yaml_out['mols']['saturation']['unmerged']['uncorrected2'] = rs2.allele_saturation()

    if umi_errors >0:
        print(f"Correcting umis with a hdist of {umi_errors} using {threads} cores")
        rsm.correct_umis(errors = umi_errors , silent = True, jobs = threads)
        rs1.correct_umis(errors = umi_errors , silent = True, jobs = threads)
        rs2.correct_umis(errors = umi_errors , silent = True, jobs = threads)
        ###yaml_out['mols']['saturation']['merged']['corrected'] = rsm.saturation()
        ###yaml_out['mols']['saturation']['unmerged']['corrected1'] = rs1.saturation()
        ###yaml_out['mols']['saturation']['unmerged']['corrected2'] = rs2.saturation()

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
    bint_db_ifn = conf['alignment']['bint_db_ifn']


# PO removed in refactoring
#    if rsm.consensus.keys() is None:
#        yaml_out['success'] = False
#        return(None)
#    else:
#        yaml_out['success'] = True
# 

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"m", 
                            fa_ofn = fa_correctedm, 
                            rs = rsm, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"1", 
                            fa_ofn = fa_corrected1, 
                            rs = rs1, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"2", fa_ofn = fa_corrected2, 
                            rs = rs2, ref_path = rcref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    
    compute_barcode_matrix_merged(fa_ifn = fa_correctedm, tab_out= merged_tab_out, bint_db_ifn = bint_db_ifn, do_rc=False)
# TODO PO verify if we need a list of corrected umis of whether they had been deleted 
    rs1c = list(rs1.umis.keys())
    rs2c = list(rs2.umis.keys())
    umi_whitelist = list(set(rs1c).intersection(set(rs2c)))

    print("whitelist:", len(umi_whitelist))

    if len(umi_whitelist)>0:
        unmerged_tab_out = conf['lineage']['paired_tab'] 

        mlt_compute_barcode_matrix_paired(fa_ifn = fa_corrected1, 
                                           fa_ifn2 = fa_corrected2, 
                                           umi_whitelist = umi_whitelist, 
                                           umi_blacklist = list(rsm.umis.keys()), 
                                           tab_out = unmerged_tab_out, 
                                           bint_db_ifn = bint_db_ifn, 
                                           do_rc = True)

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


def sc_bin_alleles(name, intid, 
    config_card_dir,
    outdir, 
    fqm, 
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
    # TODO add tadpole-based correction as an option
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
        print("Created {}".format(outdir))

    if not os.path.isdir(config_card_dir):
        os.makedirs(config_card_dir, exist_ok=True)
        print("Created {}".format(config_card_dir))
        
    out_prefix = f'{outdir}/{name}'

    # TODO clean this mess redirectin to who knows where
    # TODO do we neeed thiss?
    if intid2_R2_strand != None:
        intid2 = ogtk.UM.rev_comp(intid2_R2_strand)# same strand as initd1
    else:
        intid2 = len(intid) * "N"
   

    bint_db = pyaml.yaml.load(open(bint_db_ifn), Loader=pyaml.yaml.FullLoader)

    conf_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
    conf = {'name':name, 'desc':{}, 'outdir':outdir, 'intid1':intid, 'modality':'sc'}

    ref_seq = bint_db['ref_seq'].replace('{INTID1}', intid).upper()
    ref_name = f'{bint_db["name_plasmid"]}_intid_{intid}'
    ref_path = outdir+f'/{ref_name}.fa'
 
    anchor2 = bint_db['anchor_plasmid2']
    rxanch2 = regex.compile("ANCHOR2{e<=3}".replace("ANCHOR2", anchor2))

    ## for 10x
    umi_len = 28

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
        print(f"loaded cached ReadSet {pickled_readset}. total umis: {len(rssc.umis)}")
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
    

def compute_barcode_matrix_merged(fa_ifn, tab_out,  bint_db_ifn, tab_full_out = None, do_rc = False):
    '''
    Writes to disk a tabulated file with the corresponfing bint strings.
    Requires a fasta file that provides the pairwise aligment between a lineage
    allele and the reference
    '''

    FF = pyfaidx.Fasta(fa_ifn, as_raw=True)
    fa_entries = [i for i in FF.keys()]
    outf = open(tab_out, 'w')
    if tab_full_out is None:
        tab_full_out = tab_out.replace('.txt', "_full.txt")
    full_out = open(tab_full_out, 'w')
    #for i in range(0, int(len(fa_entries)/2), 2):
    #### !!!!!
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
    ''' with a given reference and red sequence pair, returns a lineage vector making use of a barcode interval db (bintdb)'''
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

