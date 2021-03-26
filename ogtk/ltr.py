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


def detect_integrations_dual(ifn, intid_len =8, anchor1 = "CTGTTCCTGTAGAAAT", error1 = 3, anchor2 = "CCGGACTCAGATCTCGAGCT", error2 = 2, anchor3="CGAGCGCTATGAGCGACTATGGGA", error3=3, limit = None, cores = 50):
    ''' process fastq/fastq.gz/bam file and detects number of integrations by means of anchoring seqs (anchor1, anchor2, anchor3)
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
    hits_paired = np.array(pool.map(extract_intid_pair_worker,  itertools.zip_longest([[rxc, rxc2]], chunks, fillvalue=[rxc, rxc2])))

    pool.close()
    pool.join()

    hits = [item for sublist in hits for item in sublist]
    x = pd.Series(hits).value_counts()
    xx = pd.Series([int(np.log10(i)) for i in x], index = x.index)

    valid_ints = [i for i in x[xx>(xx[0]-1)].index]
    print("Found {} valid integrations".format(len(valid_ints)))
    if len(hits_paired) >0:
        print("Paired ints found")
        x = pd.Series([item for sublist in hits_paired for item in sublist]).value_counts()
        xx = pd.Series([int(np.log10(i)) for i in x], index = x.index)
        valid_ints = [i for i in x[xx>(xx[0]-1)].index]
    print(valid_ints)
    return(valid_ints)

def export_ntop_intids_json(intid_series, n_top, ofn):
    import pyaml
    import pandas as pd
    df = pd.DataFrame({'intid_code':[f'intid_{i:02d}' for i in range(n_top)], 'intid':intid_series.head(n_top).index})
    df.to_json(ofn)    
    
def detect_integrations_single(ifn, n_intid = None, intid_len =8, anchor1 = "CTGTTCCTGTAGAAAT", error1 = 3, anchor2 = "CCGGACTCAGATCTCGAGCT", error2 = 2, limit = None, cores = 50):
    ''' For barcoded integrations with a single identifier.
    Process fastq/fastq.gz/bam file and detects number of integrations by means of anchoring seqs (anchor1, anchor2)
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
    #rxc2 = regex.compile(".*(?P<intid2>.{{{}}}).*({}){{e<={}}}".format(intid_len, anchor3, error3))

    pool = multiprocessing.Pool(cores)

    hits = np.array(pool.map(extract_intid_worker,  itertools.zip_longest([rxc], chunks, fillvalue=rxc)), dtype=object)
    #hits_paired = np.array(pool.map(extract_intid_pair_worker,  itertools.zip_longest([[rxc, rxc2]], chunks, fillvalue=[rxc, rxc2])))

    pool.close()
    pool.join()

    hits = [item for sublist in hits for item in sublist]
    x = pd.Series(hits).value_counts()
    xx = pd.Series([int(np.log10(i)) for i in x], index = x.index)
    
    if n_intid != None:
        return(x.head(n_intid))
    else:
        print(f'No expected number of integrations provided, returning an guess based on the top ranking intid counts. ')
        valid_ints = [i for i in x[xx>(xx[0]-1)].index]
        print("Found {} valid integrations".format(len(valid_ints)))
        print(valid_ints)
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
    
                id = []
                seq = []
                plus = []
                qual = []
    
    for i in valid_intids:
        out_files_dic2[i].close()
    
def sc_bin_alleles(name, intid, config_card_dir, outdir, fqm, bint_db_ifn, intid2_R2_strand = None, threads = 100, end = 5000, ranked_min_cov = 5, consensus = False, umi_errors = 1, debug = False, jobs_corr=10, alg_gapopen = 20, alg_gapextend =1, correction_dict_path = None, correction_dict = None, trimming_pattern = None):
    
    ''' Once 10x fastqs (fqm) have been merged and split by integration,
        process them into a table of barcodes for downstream analysis on R
        Two modalities to choose for a representative allelle: ranked (consensus = False) or consensus (consensus = True)'''
    # TODO add tadpole-based correction as an option
    
    bint_db = pyaml.yaml.load(open(bint_db_ifn), Loader=pyaml.yaml.FullLoader)
    ref_seq = bint_db['ref_seq'].replace('{INTID1}', intid).upper()
    ref_name = f'{bint_db["name_plasmid"]}_intid_{intid}'
    ref_path = outdir+f'/{ref_name}.fa'
 
    anchor2 = bint_db['anchor_plasmid2']
    rxanch2 = regex.compile("ANCHOR2{e<=3}".replace("ANCHOR2", anchor2))

    ## for 10x
    umi_len = 28
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
        print("Created {}".format(outdir))

    if not os.path.isdir(config_card_dir):
        os.makedirs(config_card_dir, exist_ok=True)
        print("Created {}".format(config_card_dir))
        
    # TODO clean this mess redirectin to who knows where
    # TODO do we neeed thiss?
    if intid2_R2_strand != None:
        intid2 = ogtk.UM.rev_comp(intid2_R2_strand)# same strand as initd1
    else:
        intid2 = len(intid) * "N"
   
    yaml_out_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
    yaml_out = {'name':name, 'desc':{}, 'intid1':intid, 'modality':'sc'}
    
    # align
    fa_correctedm = outdir+'/{}_corrected_merged.fa'.format(name)

    yaml_out['desc']['alignment'] = "Parameters for the pairwise alignment of every recovered allele and the reference. The corrected files fix the issue of duplicated fasta entry names"

    
    alg_mode = 'needleman'
    yaml_out['alignment'] = {}
    yaml_out['alignment']['fa_correctedm'] = fa_correctedm
    yaml_out['alignment']['alg_mode'] = alg_mode
    yaml_out['alignment']['alg_gapopen'] = alg_gapopen
    yaml_out['alignment']['alg_gapextend'] = alg_gapextend

    # more outfiles
    merged_tab_out =   outdir + '/{}_binned_barcodes_10x.txt'.format(name)
    consensus_tab_out =   outdir + '/{}_consensus_by_cell_10x.txt'.format(name)
    allele_reps_fn =   outdir + '/{}_allele_reps.txt'.format(name)
    pickled_readset= outdir + f'{name}_readset.corr.pickle'
 
    yaml_out['desc']['lineage'] = "Main output of the pre-processing scripts. Tabulated files that contain the bint barcode string"
    yaml_out['lineage'] = {}
    yaml_out['lineage']['consensus'] = consensus_tab_out
    yaml_out['lineage']['allele_reps'] = allele_reps_fn
    yaml_out['lineage']['merged_tab'] = merged_tab_out
    yaml_out['lineage']['merged_full'] = merged_tab_out.replace('.txt','_full.txt')
   
    # Molecule data
    umi_counts_ofn = outdir + '/{}_umi_counts.txt'.format(name)
    yaml_out['mols'] = {}
    yaml_out['mols']['umi_len'] = umi_len
    yaml_out['mols']['umi_correction'] = umi_errors
    yaml_out['mols']['counts'] = umi_counts_ofn 
    
    if end != None:
        print("Warning: end is not None; You are not processing all the data", end)

    # do a pass on the raw fastqs and group reads by UMI

    rssc = ogtk.UM.fastq_collapse_UMI(fqm, umi_len = umi_len, end = end, keep_rid = True, trimming_pattern = trimming_pattern)
    
    print(f"loaded rs with trimming. total umis: {len(rssc.umis)}")

    # store the first round of molecule stats
    #TODO the saturation stats should be outsourced somehow
    yaml_out['mols']['saturation'] = {}
    yaml_out['mols']['saturation']['unmerged'] = {}
    yaml_out['mols']['nreads'] = {}
    yaml_out['mols']['nreads']['umerged'] = rssc.nreads
    yaml_out['mols']['desc'] = "number of umis whose rank1 sequence is a >= [1, 2, 3, 4, 5, 10, 100, 1000]"
    yaml_out['mols']['saturation']['unmerged']['uncorrected'] = rssc.allele_saturation()


    if umi_errors >0:
        print("Correcting umis with a hdist of {}".format(umi_errors))
        rssc.correct_umis(errors = umi_errors , silent = False, jobs = jobs_corr)
        yaml_out['mols']['saturation']['unmerged']['corrected'] = rssc.saturation()
    else:
        print("Not correcting umis")
    
    # TODO this is really important for 10x, long reads
    # The first step to call an allele consensus is to get a representative sequence for a given UMI
    # thus by_alignment == False to get a ranked based representative
    rssc.consensus = ogtk.UM.do_fastq_pileup(rssc, min_cov = ranked_min_cov, threads = threads, trim_by = None, by_alignment = False)
    # save binary

    print(f"called consensuses {len(rssc.consensus)}")

    with open(pickled_readset, 'wb') as pcklout:
        pickle.dump(rssc, pcklout)
    # determine cells based on the 10x bc
 
    # the second step is to call a consensus *of the representative molecules*
    celldb = {} #stores the candidate representative sequences per cell
    cell_umis = {} #stores a list of umis that correspond to a given cell

    # correction via correction dictionaries
    if correction_dict_path != None or correction_dict != None:
        hits_dic = []
        print(f'Using cellranger dictionary', end = "...")
        if correction_dict_path != None or correction_dict == None:
            print(f'Loading correction dictionaries {correction_dict_path}')
            crd = pickle.load(open(correction_dict_path, 'rb'))
        elif correction_dict != None:
            print("using provided")
            crd = correction_dict

        for tenex_umi in rssc.consensus.keys():
            foc_cell = tenex_umi[0:16]
            if foc_cell in crd['cell_barcodes'].keys():
                hits_dic.append(1)
                # convert cells to a cellrager-corrected cell barcode
                foc_cell = crd['cell_barcodes'][foc_cell]
                if celldb.get(foc_cell, 0) == 0:
                    celldb[foc_cell] = []
                    cell_umis[foc_cell] = []
                cell_umis[foc_cell].append(tenex_umi) ## remove?
                celldb[foc_cell].append(rssc.consensus[tenex_umi]) 
            else:
                hits_dic.append(0)
        print(f'correction stats {sum(hits_dic)/len(hits_dic)}')
    else:
        for tenex_umi in rssc.consensus.keys():
            foc_cell = tenex_umi[0:16]
            if not any([i == "N" for i in foc_cell]):
                if celldb.get(foc_cell, 0) == 0:
                    celldb[foc_cell] = []
                    cell_umis[foc_cell] = []
                # expand reads
                cell_umis[foc_cell].append(tenex_umi) ## remove?
                celldb[foc_cell].append(rssc.consensus[tenex_umi]) 

    # This might change
    msa_jobs = 50
    naive = False
 
    print(f"filled celldb {len(celldb)}")
   
    #pool = multiprocessing.Pool(threads)
    print("Calling allele consensus by aligning different UMI-controlled molecules per each cell. ")
    # fname, name, seqs, ranked_min_cov, jobs = args
    job_iter = []

    
    # TODO how can we filter good candidates?
    
    allele_reps = open(allele_reps_fn, "w")
    cc = []
    with open(allele_reps_fn, "w") as allele_reps:
        for cell in celldb.keys():
            repseqs_count_list = celldb[cell]
            tenex_umis = cell_umis[cell]
            if len(repseqs_count_list)>1:
                df = pd.DataFrame(repseqs_count_list, columns = ('seqs','counts'))
                df['raw_umi'] = tenex_umis
                df['cell'] = cell
                # TODO how to treat ties?
                top = df['seqs'].value_counts().head(1)
                cc.append((cell, str(top.index.to_list()[0])))
                #get rep seq per cell and its corresponding counts.
                #rep_seq, counts = repseqs_count_list[0]
            else:
                rep_seq, counts = repseqs_count_list[0]
                cc.append((cell, rep_seq))

            #cell_name = cell+"_reps"+str(len(tenex_umis))
            cell_name = cell #+"_reps"+str(len(tenex_umis))
            # 1 - keep a record on disk of the different candidates and their counts
            # 2 - append to cc the "best" candidate for a given cell 
            ##for allele, (umi, counts) in zip(rep_seqs, tenex_umis):
            ##    #TODO add a real count metric
            ##    #counts = int(umi_counts.split("_")[1])
            ##    umi = umi[0][16:]
            ##    reps = [allele for i in range(counts)]
            ##    job_iter.append((outdir, cell_name, reps, ranked_min_cov, msa_jobs, naive))
            ##    #allele_reps.write("{}\n".format("\t".join([cell, str(counts), umi, allele])))
            ##    allele_reps.write("\t".join([cell, str(counts), umi, allele])+'\n')
        
    # used for debugging
    if False:
        debug_iter = [i for i in job_iter if len(i[2])>3 ]
        print("There are {} cells and {} with more than 3 representatives".format(len(job_iter), len(debug_iter)))
        for i in debug_iter:
            ogtk.UM.mafft_consensus(i)

    #cc = pool.map(ogtk.UM.mafft_consensus, iter(job_iter))
    #pool.close()
    #pool.join()
 
    print(f"filled cc {len(cc)}")
    
    cell_consensus = dict(cc)
    cons_out = open(consensus_tab_out, 'w') 
    for cell in cell_consensus.keys():
        cell_cons = cell_consensus[cell]
        cons_out.write("\t".join([cell, cell_cons])+'\n')
    cons_out.close() 
   
    if len(cell_consensus.keys()) == 0 :
        yaml_out['success'] = False
        print(f'Failed to generate representative sequences for each cell')
        return(None)
    else:
        print(f'Success')
        yaml_out['success'] = True
 
    ## make sure to include the intid on the reference to help the alignment

    print("creating fasta ref")
    create_fasta_ref(ref_name, ref_seq, ref_path)


    print("aligning reads tp ref")
    align_reads_to_ref(name = outdir+"/"+intid, fa_ofn = fa_correctedm, 
                            consensus_dict = cell_consensus, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)


    print(f"featurizing sequences to {merged_tab_out}")
    compute_barcode_matrix_merged(fa_ifn = fa_correctedm, tab_out = merged_tab_out, bint_db_ifn = bint_db_ifn, do_rc = False)
   
    # save the count information for the UMIs 
    with open(umi_counts_ofn, 'w') as fh:
        fh.write(f'count\tumi\n')
        counts = rssc.umi_counts()
        for count,umi in zip(counts, counts.index):
            fh.write(f'{count}\t{umi}\n')
                
    # save run settings 
    yaml_stream = open(yaml_out_fn, 'w')
    pyaml.yaml.dump(yaml_out, yaml_stream)
    yaml_stream.close() 
    print("Run settings saved to config card:{}".format(yaml_out_fn))

    
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
    '''Writes to disk a tabulated file with the corresponfing bint strings.
    Requires a fasta file that provides the pairwise aligment between a lineage
    allele and the reference'''

    FF = pyfaidx.Fasta(fa_ifn, as_raw=True)
    fa_entries = [i for i in FF.keys()]
    outf = open(tab_out, 'w')
    if tab_full_out == None:
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


def create_fasta_ref(ref_name, ref_seq, ref_path):
    fout = open(ref_path, 'w')
    fout.write('>{}\n{}\n'.format(ref_name, ref_seq))
    fout.close()

def mltb_bin_alleles(name, intid, intid2_R2_strand, config_card_dir, outdir, fq1, fq2, threads = 100, end = 5000, ranked_min_cov = 2, umi_len = 17, umi_errors = 1, alg_gapopen = 20, alg_gapextend =1):
    '''Once fastq files have been split by: RPI/rev_id/intid process them into a table of barcodes for downstream analysis on R; ranked_min_cov determines the threshold to call a representative rank1 sequence as valid'''
    intid2 = ogtk.UM.rev_comp(intid2_R2_strand)# same strand as initd1
    # TODO consider generating a conf file out of the pre-processing step
    # constant vars
    ref_path = outdir+'/hspdrv7_scgstl_{}.fa'.format(name)
    rcref_path = outdir+'/rc_hspdrv7_scgstl_{}.fa'.format(name)
    ref_name = 'hspdrv7_scgstl_{}'.format(name)
    #ref_seq = 'gcgccaccacctgttcctgtagaaat{}ccggactcagatctcgagctcaagcttcggacagcagtatcatcgactaTGGagtcgagagcgcgctcgtcgactaTGGagtcgtcagcagtactactgacgaTGGagtcgacagcagtgtgtgagtctaTGGagtcgagagcatagacatcgagtaTGGagtcgactacagtcgctacgactaTGGagtcgacagagatatcatgcagtaTGGagtcgacagcagtatctgctgtcaTGGagtcgactgcacgacagtcgactaTGGAGTCG{}cgagcgctatgagcgactatgc'.format(intid, intid2).upper()
    ref_seq = 'gccaccacctgttcctgtagaaat{}ccggactcagatctcgagctcaagcttcggacagcagtatcatcgactaTGGagtcgagagcgcgctcgtcgactaTGGagtcgtcagcagtactactgacgaTGGagtcgacagcagtgtgtgagtctaTGGagtcgagagcatagacatcgagtaTGGagtcgactacagtcgctacgactaTGGagtcgacagagatatcatgcagtaTGGagtcgacagcagtatctgctgtcaTGGagtcgactgcacgacagtcgactaTGGAGTCG{}cgagcgctatgagcgactatgc'.format(intid, intid2).upper()
    bint_db_ifn = '/local/users/polivar/src/projects/mltracer/conf/gstlt_barcode_intervsx10.yaml'

    yaml_out_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
    yaml_out = {'name':name, 'desc':{}, 'intid1':intid, 'intid2':intid2, 'modality':'bulk'}
    
    # for intid in dataset
    # merge
    fqm =       outdir + '/{}_merged.fastq'.format(name)
    fqum1 =     outdir + '/{}_unmerged_R1.fastq'.format(name)
    fqum2 =     outdir + '/{}_unmerged_R2.fastq'.format(name)
    ihist =     outdir + '/{}_ihist.txt'.format(name)
    log_merge = outdir + '/{}_merge_stats.txt'.format(name)
    cmd_merge = "bbmerge-auto.sh in1={} in2={} out={} outu={} outu2={} ihist={} ecct extend2=20 iterations=5 interleaved=false -Xmx2g ".format(fq1, fq2, fqm, fqum1, fqum2, ihist)
    
    yaml_out['desc']['read_merge'] = "Parameters used for the merging of read 1 and read 2"

    umi_counts_ofn = outdir + '/{}_umi_counts.txt'.format(name)
    yaml_out['mols'] = {}
    yaml_out['mols']['umi_len'] = umi_len
    yaml_out['mols']['umi_correction'] = umi_errors
    yaml_out['mols']['counts'] = umi_counts_ofn 
    
    yaml_out['read_merge'] = {}
    yaml_out['read_merge']['fqm'] = fqm
    yaml_out['read_merge']['fqum1'] = fqum1
    yaml_out['read_merge']['fqum2'] = fqum2
    yaml_out['read_merge']['ihist'] = ihist
    yaml_out['read_merge']['log_merge'] = log_merge

    # align
    fa_correctedm = outdir+'/{}_corrected_merged.fa'.format(name)
    fa_corrected1 = outdir+'/{}_corrected_R1.fa'.format(name)
    fa_corrected2 = outdir+'/{}_corrected_R2.fa'.format(name)

    yaml_out['desc']['alignment'] = "Parameters for the pairwise alignment of every recovered allele and the reference. The corrected files fix the issue of duplicated fasta entry names"

    yaml_out['alignment'] = {}
    yaml_out['alignment']['fa_correctedm'] = fa_correctedm
    yaml_out['alignment']['fa_corrected1'] = fa_corrected1
    yaml_out['alignment']['fa_corrected2'] = fa_corrected2

    # more outfiles
    merged_tab_out =   outdir + '/{}_barcodes_merged.txt'.format(name)
 
    yaml_out['desc']['lineage'] = "Main output of the pre-processing scripts. Tabulated files that contain the bint barcode string"
    yaml_out['lineage'] = {}
    yaml_out['lineage']['merged_tab'] = merged_tab_out
    yaml_out['lineage']['merged_full'] = merged_tab_out.replace('.txt','_full.txt')

   
    # merge paired end reads if possible

    pp = subprocess.run(cmd_merge.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    fout = open(log_merge, 'wb')
    fout.write(pp.stdout)
    fout.write(pp.stderr)
    fout.close()

    if end != None:
        print("Warning: end is not None; You are not processing all the data", end)

    # do a pass on the raw fastqs and group reads by UMI
    ## for merged (rsm) and unmerged (rs1, rs2)
    rsm =           ogtk.UM.fastq_collapse_UMI(fqm, umi_len = umi_len, end = end, keep_rid = True)
    rs1, rs2 =      ogtk.UM.pfastq_collapse_UMI(fqum1, fqum2, umi_len = umi_len, end = end)

    # remove umis that exist on the merged set out of the unmerged
    unmerged_by_error = set(rsm.umi_list()).intersection(set(rs1.umi_list()))
    rs1.delete(unmerged_by_error)
    rs2.delete(unmerged_by_error)

    #TODO the saturation stats should be outsourced somehow and reported for all three cases not just the merged (rsm) 
    yaml_out['mols']['saturation'] = {}
    yaml_out['mols']['saturation']['merged'] = {}
    yaml_out['mols']['saturation']['unmerged'] = {}
    yaml_out['mols']['nreads'] = {}
    yaml_out['mols']['nreads']['total'] = sum([rsm.nreads, rs1.nreads])
    yaml_out['mols']['nreads']['merged'] = rsm.nreads
    yaml_out['mols']['nreads']['umerged'] = rs1.nreads
    yaml_out['mols']['desc'] = "number of umis whose rank1 sequence is a >= [1, 2, 3, 4, 5, 10, 100, 1000]"
    yaml_out['mols']['saturation']['merged']['uncorrected'] = rsm.allele_saturation()
    yaml_out['mols']['saturation']['unmerged']['uncorrected1'] = rs1.allele_saturation()
    yaml_out['mols']['saturation']['unmerged']['uncorrected2'] = rs2.allele_saturation()

    if umi_errors >0:
        print("Correcting umis with a hdist of {}".format(umi_errors))
        rsm.correct_umis(errors = umi_errors , silent = True, jobs = 80)
        rs1.correct_umis(errors = umi_errors , silent = True, jobs = 80)
        rs2.correct_umis(errors = umi_errors , silent = True, jobs = 80)
        yaml_out['mols']['saturation']['merged']['corrected'] = rsm.saturation()
        yaml_out['mols']['saturation']['unmerged']['corrected1'] = rs1.saturation()
        yaml_out['mols']['saturation']['unmerged']['corrected2'] = rs2.saturation()

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
    #alg_gapopen = 20
    #alg_gapextend = 1

    yaml_out['alignment']['alg_mode'] = alg_mode
    yaml_out['alignment']['alg_gapopen'] = alg_gapopen
    yaml_out['alignment']['alg_gapextend'] = alg_gapextend


    if len(rsm.consensus.keys()) == 0 :
        yaml_out['success'] = False
        return(None)
    else:
        yaml_out['success'] = True
 

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"m", fa_ofn = fa_correctedm, 
                            consensus_dict = rsm.consensus, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"1", fa_ofn = fa_corrected1, 
                            consensus_dict = rs1.consensus, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"2", fa_ofn = fa_corrected2, 
                            consensus_dict = rs2.consensus, ref_path = rcref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    
    mlt_compute_barcode_matrix_merged(fa_ifn = fa_correctedm, tab_out= merged_tab_out, bint_db_ifn = bint_db_ifn, do_rc=False)

    rs1c = list(rs1.consensus.keys())
    rs2c = list(rs2.consensus.keys())
    umi_whitelist = list(set(rs1c).intersection(set(rs2c)))

    print("whitelist:", len(umi_whitelist))

    if len(umi_whitelist)>0:
        unmerged_tab_out = outdir + '/{}_barcodes_unmerged.txt'.format(name)
        yaml_out['lineage']['paired_tab'] = unmerged_tab_out 
        yaml_out['lineage']['paired_full'] = unmerged_tab_out.replace('.txt','_full.txt')

        mlt_compute_barcode_matrix_paired(fa_ifn = fa_corrected1, fa_ifn2 = fa_corrected2, umi_whitelist = umi_whitelist, 
                                           umi_blacklist = list(rsm.consensus.keys()), tab_out = unmerged_tab_out, bint_db_ifn = bint_db_ifn, do_rc = True)

    # save the count information for the UMIs 
    with open(umi_counts_ofn, 'w') as fh:
        fh.write(f'count\tumi\n')
        for i in (rsm, rs1):
            counts = i.umi_counts()
            for count,umi in zip(counts, counts.index):
                fh.write(f'{count}\t{umi}\n')
    
    # save run settings 
    yaml_stream = open(yaml_out_fn, 'w')
    yaml.dump(yaml_out, yaml_stream)
    yaml_stream.close() 
    print("Run settings saved to config card:\n{}".format(yaml_out_fn))

   
def compute_dist_mat(mat, bcs, bias, cores = 50, dont_include = 3):
    cores = int(cores)
    #print(f"Warning: FORCED for  using {cores} cores")
    print(f'got {len(bcs)} barcodes')
    if cores == 0:
        print("zero cores")
        cc = similarity_score((mat, bcs, bias, dont_include)) 
        return(cc)
    else:
        print("pooling cores")
        chunks = np.array_split(bcs, cores)
        it = iter([(mat, chunk, bias, dont_include) for chunk in chunks])
        pool = multiprocessing.Pool(int(cores))
        cc = pool.map(similarity_score, it)
        pool.close()
        pool.join()
        cum_mat = cc[0]
        for i in cc[1:]:
            cum_mat = cum_mat + i
        return(cum_mat)

def similarity_score(args):
    from scipy import sparse
    mat, bcs, bias, dont_include = args
    first = True
    #bias = np.array([1.25, 1.25, 1.25, 1.25, 0.25, 0.25, 0.25 , 0.25, 0.25, 0.25])
    #bias = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    w = 0.1
    #bias = np.array([w, w, w, w, w, w, w, w, w, w])
    #bias = np.array([5*w, 5*w, 5*w, 5*w, 5*w, w, w, w, w, w])
    #bias = np.array([w/5, w/5, w/5, w/5, w/5, w, w, w, w, w])

    bias = np.array([float(i) for i in bias])
    bias = bias/len(bias)

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
                if i%500==0:
                    print(len(bcs)/i)
                cum_mat = cum_mat + (ee.todense()*1)
    return(cum_mat)

def to_numpy(mat, ofn):
    import pickle
    with open(ofn, 'wb') as out:
        pickle.dump(mat, out)


