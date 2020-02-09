import os
import pdb
import pysam
from pyfaidx import Fasta
import yaml
import subprocess
import operator
import copy
import regex
import itertools
import multiprocessing
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from .aln_to_mat import Alg as Alg


class Read_Set:
    '''
    Hosts a list UM instances and provides a number of functions to process them
    '''
    def __init__(self, name = ''):
        # gets filled with UM instances
        self.umis = {}
        self.name = name
        self.ori_umis = None
        self.corr_stages = []

        # RID to UMI dic
        self.rid_umi = {} # holds a k,v for every readid, UMI
        ## mapped attributes
        self.original_header = None
        self.header = None

        self.consensus = None
        self.conf = None
        self.prefix = ''

    def set_header(self, header):
        self.original_header = header
        self.header = header.as_dict()

    def write_fasta(self, outfn  = None, balance = False):
        outfn = self.prefix+ '.cons.fa'
        consensuses_fa = open(outfn, 'w')
        # balance the alignment using the full reference 100x times
        if balance:
            fa_ref = Fasta(self.conf['ref_fa'])
            for i in range(50):
                consensuses_fa.write(">%s\n%s\n"%(fa_ref[0].name+"_balancer_"+str(i),fa_ref[0][:]))
        for entry in self.consensus.keys():
            if len(self.consensus[entry]) >0:
                consensuses_fa.write(">%s\n%s\n"%(entry, self.consensus[entry].replace('\n', '')))
        print('save consensuses to %s' % (outfn))
        consensuses_fa.close()

    def write_tabular(self):
        outfn = self.prefix+ '.cons.txt'
        consensuses_tab = open(outfn, 'w')
        for entry in self.consensus.keys():
            
            hits, umi, cigar = entry.split("_") 
            if len(self.consensus[entry]) >0:
                consensuses_tab.write("%s\t%s\n"%("\t".join(entry.split("_")), self.consensus[entry].replace('\n', '')))
        print('save consensuses to %s' % (outfn))
        consensuses_tab.close()


    def add_umi(self, umi, seq, alignment=None):
        if self.umis.get(umi, 0) == 0:
            self.umis[umi] = UM(self.name, umi)
        self.umis[umi].seqs.append(seq)
        if alignment.__class__.__name__ == 'AlignedSegment':
            if self.umis[umi].sams.get(alignment.cigarstring, 0) == 0:
                self.umis[umi].sams[alignment.cigarstring] = []
            self.umis[umi].sams[alignment.cigarstring].append(alignment)
            self.umis[umi].is_mapped = True

    def add_rid(self, umi, rid):
        # TODO add some sanity checks here
        self.rid_umi[rid] = umi

    def umi_list(self):
        return([i for i in self.umi_counts().index])
    def umi_counts(self):
        by_counts = pd.Series([len(self.umis[i].seqs) for i in self.umis.keys()], index=self.umis.keys()).sort_values(ascending=False)
        by_counts = pd.DataFrame({"counts":by_counts, "umi":by_counts.index}, index=by_counts.index)
        by_counts = by_counts.sort_values(by=["counts", "umi"], ascending = False)
        return(by_counts['counts'])

    def umi_countss(self):
        by_counts = sorted([(len(self.umis[i].seqs), i) for i in self.umis.keys()], key = lambda x: x[1])
        
        return(by_counts)

    def correct_umis(self, errors = 1, mode = "dist", jobs =10, silent = True):
        if self.ori_umis == None:
            self.ori_umis = [i for i in self.umis.keys()]
        self.corr_stages.append(self.umi_counts())
        seqs = [len(self.umis[i].seqs) for i in self.umis.keys()]
        by_counts = pd.Series(seqs, index = [i for i in self.umis.keys()]).sort_values()
        if not silent:
            print("start", len(seqs))

        old_to_new = merge_all(self.umi_list(), errors = errors, mode = mode, jobs = jobs)
        for oldk,newk in old_to_new[::-1]:
            if oldk != newk:
                for i in self.umis[oldk].seqs:
                    self.umis[newk].seqs.append(copy.deepcopy(i))
                del self.umis[oldk]
        
        corrected_count = sum([1 for i,v in old_to_new if i !=v])
        if not silent:
            print("Found %s corrections" %(corrected_count))
            print("end", len([i for i in self.umis.keys()]))
        
        #if corrected_count>0:
        #    self.correct_umis(errors = errors, mode = mode)

    def do_msa(self, plot = False):
        assert self.consensus, "consensus has not been generated: run .do_pileup to generate it"
        if len(self.consensus) == 0:
            return()
        prefix = self.prefix
        conf_clustalo = self.conf['clustalo']
        clust = {'ifn': prefix+".cons.fa",
                 'ofn': prefix+".msa.fa",
                 'tree': prefix+".msa.tree",
                 'dm': prefix+".msa.dm"
                }
        #clust_args = (clust['ifn'], conf_clustalo['flags'], str(conf_clustalo['threads']), clust['ofn'], clust['tree'], clust['dm'])
        clust_args = (clust['ifn'], conf_clustalo['flags'], str(conf_clustalo['threads']), clust['ofn'])
        cmd_clustalo = "clustalo -i %s %s --threads %s --outfile %s  " % clust_args
        #cmd_clustalo = "clustalo -i %s %s --threads %s --outfile %s --guidetree-out %s --distmat-out %s" % clust_args
        print(cmd_clustalo)
        subprocess.run(cmd_clustalo.split())
        clean_fasta(clust['ofn'], prefix)
        
        ## post processing
        alignment = Alg(fastafn = clust['ofn'], freqfn = prefix + ".msa.freqs", colorfn = prefix + ".msa.colors")
        self.msa = alignment
        return(alignment)

    def do_pileup(self, fa_ref = 'refs/rsa_plasmid.fa', start= 0 , end = None, region = 'chr_rsa:1898-2036', threads = 50, balance = True, min_cov = 5, prefix = None):
        consensuses = do_pileup(self, fa_ref = fa_ref, start= start , end = end, region = region, threads = threads, balance = balance, min_cov = min_cov, prefix = prefix)
        return(consensuses)

    def plot_hdist(self, outpng):
        plot_hdist(self, outpng)
    def get_hdists(self):
        seqs = self.umi_list()
        dists = hdist_all(seqs)
        return(dists)
       

class UM:
    def __init__(self, sample, umi):
        self.sample = sample
        self.umi = umi
        self.seqs = []
        self.ori_umis = None
        self.is_mapped = False
     
        # stores sam entries using a cigar string as key
        self.sams = {}
        
    def add_seq(self, seq):
        self.seqs.append(seq)

    def export_to_fastq(self, fasta_ofn = "/home/polivar/test.fa"):
        if not self.is_mapped:
            cc = pd.Series(self.seqs).value_counts()
            write_to_fasta(fasta_ofn, [[str(i)+self.umi, v] for i,v in zip(cc, cc.index)])
        else:
            seqs = [ [i, ] for i in self.sams.keys()]

def pfastq_collapse_UMI(fastq_ifn1, fastq_ifn2, umi_start=0, umi_len=12, end=None):
    '''Constructs a paired Readset by transfering the UMI from R1 to R2 via the
    read id. TODO and converts R2 to the same strand''' 
    rset1 = fastq_collapse_UMI(fastq_ifn1, umi_len = umi_len, keep_rid = True, end = end)
    rset2 = fastq_collapse_UMI(fastq_ifn2, umi_len = umi_len, rid_umi = rset1.rid_umi, do_rc = False, end = end)
    
    return([rset1, rset2])
    
def fastq_collapse_UMI(fastq_ifn, name = None, umi_start=0, umi_len=12, end=None, keep_rid=False, rid_umi=None, do_rc = False):
    '''
    Process a fastq file and collapses reads by UMI making use of their consensus seq
    Keeping the read id (rid) helps to match the UMI from read1 to read2
    If a rid_umi dic is provided, it uses it instead of the sequence information to call a UMI
    '''
    if end != None:
        end = int(end)
    if name == None:
        name = fastq_ifn
    readset = Read_Set(name=fastq_ifn)
    
    fq = open(fastq_ifn).readlines()
    reads =     itertools.islice(fq, 1, end, 4)
    rids =      itertools.islice(fq, 0, end, 4)
    for i, (read, rid) in enumerate(zip(reads, rids)):
        read = read.strip()
        if rid_umi == None:
            umi = read[umi_start:umi_len] 
            seq = read[umi_len:]
        else:
            umi = rid_umi[rid.split(" ")[0]]
            seq = read.strip()
        if do_rc:
            seq = rev_comp(seq) 

        readset.add_umi(umi, seq)

        if keep_rid:
            readset.add_rid(umi, rid.split(" ")[0])

#
#####
#    if keep_rid: # flag to store read id usign the UMI as a key
#        rids = itertools.islice(fq, 0, end, 4)
#        for i, (read, rid) in enumerate(zip(it, rids)):
#            read = read.strip()
#            rid = rid.split(" ")[0] 
#            umi = read[umi_start:umi_len] 
#            seq = read[umi_len:]
#            readset.add_umi(umi, seq)
#            readset.add_rid(umi, rid)
#    else:
#        for i,read in enumerate(it):
#            read = read.strip()
#            if rid_umi == None:
#                umi = read[umi_start:umi_len] 
#            else:   
#                umi = rid_umi[rid]
#            seq = read[umi_len:]
#            readset.add_umi(umi, seq)
    print(i, "reads were processed")
    return(readset)

def bam_collapse_UMI(bam_ifn, umi_start=0, umi_len=12, end=None):
    '''
    Process a fastq file and collapse reads by UMI making use of their consensus seq
    '''
    print("processing %s" % bam_ifn)
    bamfile = pysam.pysam.AlignmentFile(bam_ifn)  
    readset = Read_Set(name = bam_ifn)
    readset.set_header(bamfile.header.copy())
    if end != None:
        end = int(end)
    bam_it = itertools.islice(bamfile, 0, end)
    for i,bam in enumerate(bam_it):
        print("\r%s"%i, end='')
        read = bam.seq
        umi = read[umi_start:umi_len] 
        seq = read[umi_len-1:]
        readset.add_umi(umi, seq, bam)
        
    print(i, "reads were processed")
    return(readset)

def fit_header_to_umis(readset):
   self = readset
   self.header = self.original_header.as_dict()
   ln = self.original_header['SQ'][0]['LN']
   for umi in self.umis.keys():
    for cigar in self.umis[umi].sams.keys():
       new_ref = return_new_ref(self.umis[umi], cigar)
       self.header['SQ'].append({'SN': new_ref, 'LN': ln})

def return_new_ref(umi, cigar):
    self = umi
    """
    Returns a code that matches a new contig for a given umi
    input is an instance of a UM class
    """
    n_reads = str(len(self.sams[cigar]))
    return("_".join([n_reads, umi.umi, cigar]))

def return_new_ref_cigar(readset, umi, cigarstring):
    self = readset
    """
    Returns a code that matches a new contig for a given umi
    """
    n_cigars = str(pd.value_counts([i.cigarstring for i in self.umis[umi]])[cigarstring])
    #return("_".join([self.mmc[0], self.hash, umi, n_cigars, cigarstring]))
    return("_".join([self.mmc[0], self.hash, umi, n_cigars]))


 
def to_sam(readset, ln = 3497, prefix= ''):
   self = readset
   """
   generates a sam file (with right headers) with upated contig fields
   corresponding to sequence variant (cigar) found for each umi 
   """
   consensus_sam = prefix+'.consensus'+".sam"
   fit_header_to_umis(self)
   sam_out = []
   header = self.header

   for umi in self.umis.keys():
       for cigar in self.umis[umi].sams.keys():
         for sam in self.umis[umi].sams[cigar]:
           new_sam = sam.to_string().split('\t')
           #new_sam[2] = self.return_new_ref(umi, sam.cigarstring)
           new_sam[2] = return_new_ref(self.umis[umi], cigar)
           sam_out.append(new_sam)
   sout = pysam.AlignmentFile(consensus_sam, 'w', header = header)
   sout.close()
   sout = open(consensus_sam, "a")

   for i in sam_out:
       sout.write("\t".join(i)+"\n")
   sout.close()
   return(consensus_sam)

def mafft_consensus(args):
    # TODO implement a quality check
    ''' Generate consensus sequences of a set of reads. Returns a (key, consensus) tuple'''
    fname, name, seqs, min_cov, jobs  = args
        
    #TODO trim by regex?
    #consensus = ()
    if len(seqs) > 1:
        umi_counts = pd.Series(seqs).value_counts()
        # TODO - do we really want to do this? What would happen if we have 5 reads with the real sequence of all seen only once?
        umi_counts = umi_counts[umi_counts >= min_cov]

        if len(umi_counts) > 1:
            mafft_in = '{}_{}_mafft_in.fa'.format(fname, name)
            mafft_out = '{}_{}_mafft_out.fa'.format(fname, name)

            fasta_entry_names = ["_".join([name, str(i), str(v)]) for i,v in enumerate(umi_counts.to_list())]
            write_to_fasta(mafft_in, zip(fasta_entry_names, umi_counts.index.to_list()))
            
            # run maff
            cmd_mafft = "mafft --maxiterate 1000 --genafpair {} ".format(mafft_in, mafft_out)
            
            # for pyhton 3.7 universal_newlines=True == text=True
            # not needed if writing using bytes
            # this is the case since the Fasta class cannot read from a stream and must open a file 
            #pp = subprocess.run(cmd_mafft.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            
            pp = subprocess.run(cmd_mafft.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # TODO maybe process the output here already?
            fout = open(mafft_out, 'wb')
            fout.write(pp.stdout)
            fout.close()

            ferr = open(mafft_out+'_err', 'wb')
            ferr.write(pp.stderr)
            ferr.close()
            
            msa_fa = Fasta(mafft_out)
            # concatenate aligned sequences (all the same size) 
            read_stream = ""
            for read in msa_fa.keys():
                read_stream = read_stream + str(msa_fa[read][:])

            # chop the stream into a matrix and trasfer it to the pandas universe
            msa_mat = [i for i in read_stream]
            msa_mat = np.array(msa_mat) 
            msa_mat = pd.DataFrame(msa_mat.reshape((len(msa_fa.keys()),-1)))

            # extract the consensus
            read_len = msa_mat.shape[1]
            
            read_counts = [int(i.split("_")[2]) for i in msa_fa.keys()]
            
            umi_consensus = ''
            for i in range(read_len):
                variants = msa_mat[i]
                # multiply each nt by the times the read was found, concatenate into string and transform into a vector
                column_expansion = [i for i in ''.join([var*count for var,count in zip(variants, read_counts)])]
                column_counts = pd.Series(column_expansion).value_counts()
                umi_consensus = umi_consensus + column_counts.head(1).index.to_list()[0].upper()
            #umi_consensus = ''.join([msa_mat[i].value_counts().head(1).index.to_list()[0] for i in range(read_len)]).upper()

            consensus = (name, umi_consensus)
        else:
            consensus = ("drop-not-cov", name)
    else:
        consensus = ("drop-one-seq", name)
             
    return(consensus) 


def do_fastq_pileup(readset, min_cov = 2, threads =  100, threads_alg = 1, trim_by = None):
    ''' Generate consensus sequences of sets of reads, grouped by UMI '''
    pool = multiprocessing.Pool(threads)
    # define the list of pool arguments
    # name, seqs, min_cov, jobs = args
    it = iter([(readset.name.replace('.fastq', ''), umi, readset.umis[umi].seqs, min_cov, threads_alg) for umi in readset.umis.keys()])
    # collapsing reads by consensus seems to be incorrect, deprecating the use of such (e.g. mafft alignment) in favor for just getting the most frequent sequence
    #cc = pool.map(mafft_consensus, it)
    cc = []
    for umi in readset.umis.keys():
        counts = pd.Series(readset.umis[umi].seqs).value_counts()
        if counts[0] >= min_cov:
            cc.append((umi, counts.head(1).index.to_list()[0]))
    #cc = [i for i in cc if "drop" not in i[0]]
    consensuses = trim_by_regex(dict(cc) , trim_by, span_index=1) if trim_by != None else dict(cc) 
    
    return(consensuses) 

 
def do_pileup(readset, fa_ref = 'refs/rsa_plasmid.fa', start= 0 , end = None, region = 'chr_rsa:1898-2036', threads = 50, balance = True, min_cov = 5, prefix = None):
    if prefix == None:
        prefix = readset.prefix
    self = readset
    #if self.consensus != None:
    #    return(None)
    # create SAM file that maps to Umi-based contigs
    print("Creating SAM")
    sam_tmp = to_sam(self, prefix = prefix)
    bam_tmp = sam_tmp.replace('sam', 'bam')

    cmd_sort = "samtools sort %s -o %s" % (sam_tmp, bam_tmp)
    cmd_index = "samtools index %s" % (bam_tmp)
    print(cmd_sort)
    print(cmd_index)

    # TODO : update region

    
    print("samtools")
    subprocess.run(cmd_sort.split())
    subprocess.run(cmd_index.split())
   
    print("opening %s" %(bam_tmp)) 
    cons_bam = pysam.AlignmentFile(bam_tmp, 'rb')
    #self.consensus = {}
    
    print("consensusing")

    sorted_entries = sorted([i['SN'] for i in cons_bam.header.as_dict()['SQ']][1:], key = lambda x: int(x.split("_")[0]))
    print("new ref bam contigs: ", len(sorted_entries)) #[i['SN'] for i in cons_bam.header.as_dict()['SQ']]))
    
   
    consensuses = {}
    cc = [] 
    pool = multiprocessing.Pool(threads)
    #entry, cons_bam, min_cov = params
    it = iter([[i, bam_tmp, min_cov] for i in sorted_entries])
    cc = pool.map(return_cons, it)

    for k,v in cc:
        if k != None:
            consensuses[k] = v
    print("kept %s hits" % len(consensuses))
    return(consensuses)
    
def iterator_collapse_UMI(name, field, it, umi_start=0, umi_len=12, end=None):
    readset = Read_Set(name=name)
    for i,read in enumerate(it):
        read = read.split('\t')[field]
        read = read.strip()
        umi = read[umi_start:umi_len] 
        seq = read[umi_len-1:]
        umis.add_umi(umi, seq)
    print(i, "reads")
    return(umis)    

def split_fastq(sample_name, out_prefix, f1_fn, f2_fn, vsamples, anchor, nreads = None, keep_id = False ):
    ''' Splits fastqs given an anchor on R1 and a sample index (usually introduced via REV primer) on R2 '''
    filter_good = 0
    filter_bad = 0
    if not os.path.exists(out_prefix):
        os.mkdir(out_prefix)
    inv1 = open('%s%s_invalid_R1.fastq'%(out_prefix, sample_name), 'w')
    inv2 = open('%s%s_invalid_R2.fastq'%(out_prefix, sample_name), 'w')
    
    dic_files1 = {}
    dic_files2 = {}
    
    for v in vsamples:
        dic_files1[v] = open(out_prefix + '%s_%s_R1.fastq'% (sample_name, v), 'w')
        dic_files2[v] = open(out_prefix + '%s_%s_R2.fastq'% (sample_name, v), 'w')
    
    print("processing files:\n%s\n%s" %(f1_fn, f2_fn))
    f1 = open(f1_fn)
    f2 = open(f2_fn)
    
    if nreads != None:
        nreads = nreads *4
    
    i1 = itertools.islice(f1, 0, nreads, 1)
    i2 = itertools.islice(f2, 0, nreads, 1)
    
    is_same = True
    id =[]
    seq = []
    plus = []
    qual = []
    
    
    for i,v in enumerate(zip(i1,i2)):
        l1 = v[0].strip()
        l2 = v[1].strip() 
               
        if i%4 == 0:
            # this corresponds to read id
            id = [l1, l2]
        if i%4 == 1:
            # this corresponds to read seq
            seq = [l1, l2]   
            sample_id = l2[0:6]
        if i%4 == 2:
            # this corresponds to read '+'
            plus = [l1, l2]
        if i%4 == 3:
            # this corresponds to read qual
            qual = [l1, l2]
    
        if qual:
            match = anchor.search(seq[0])
            if sample_id in vsamples and match:
                if keep_id:
                    sid_string = ''
                else:
                    sid_string = " sample_id:{}".format(sample_id)
                dic_files1[sample_id].write("\n".join([id[0]+sid_string, seq[0], plus[0], qual[0]]) +'\n')
                dic_files2[sample_id].write("\n".join([id[1]+sid_string, seq[1], plus[1], qual[1]]) +'\n')
                filter_good+=1
            else:
                inv1.write("\n".join([id[0]+ " sample_id:%s"%(sample_id), seq[0], plus[0], qual[0]]) +'\n')
                inv2.write("\n".join([id[1]+ " sample_id:%s"%(sample_id), seq[1], plus[1], qual[1]]) +'\n')
                filter_bad+=1
            id =[]
            seq = []
            plus = []
            qual = []

    
    inv1.close()
    inv2.close()
   
    for v in vsamples:
        dic_files1[v].close()
        dic_files2[v].close()

    total = filter_good + filter_bad
    print(sample_name, i/4)
    print("collected %.2fM (%.2f) valid reads. Filtered %.2fM (%.2f) reads.\ntotal %s"%(filter_good/1e6, filter_good/total, filter_bad/1e6, filter_bad/total, total))
 
 
def string_hamming_distance(args):
    str1 = args[0]
    str2 = args[1]
    #https://github.com/indrops/indrops/blob/2ad4669b72ea6ba794f962ed95b778560bdfde4a/indrops.py
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.
    eg "karolin" and "kathrin" is 3.
    """
    return sum(map(operator.ne, str1, str2))

def merge_umi_to_pool(args):
    '''
    Returns a pair that maps the old umi to the new one
    Expects UMIs sorted by ocurrence
    '''
    seq = args[0]
    rank = args[1][0]
    pool = args[1][1]
    errors = args[1][2]
    mode = args[1][3]
   
    
    #TODO think on how one can one reduce the number of cases (comparions, e.g. subset of pool)
    for i,best in enumerate(pool):
        if mode == "regex":
            match = regex.search(seq+'{e<='+str(errors)+'}', best)
            if match and rank < 1:
                return([seq, best])
        if mode == "dist":
            dist = string_hamming_distance((seq, best))
            if dist <= errors and i < rank:
                return([seq, best])
    return([seq, seq])

def compare_umi_to_pool(args):
    '''
    Returns the hammin distance between a given UMI to the rest of the pool
    '''
    seq = args[0]
    index = args[1][0]
    pool = args[1][1]
   
    dists = [] 
    #TODO think on how one can one reduce the number of cases (comparions, e.g. subset of pool)
    for i,best in enumerate(pool):
        dists.append(string_hamming_distance((seq, best)))
    return(dists)

def merge_all(seqs, jobs = 10, errors = 1, mode = "regex"):
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [[i, seqs, errors, mode] for i,v in enumerate(seqs)], fillvalue=seqs)
    return(pool.map(merge_umi_to_pool, it))

def hdist_all(seqs, jobs = 10):
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [[i, seqs] for i,v in enumerate(seqs)], fillvalue=seqs)
    dists = pool.map(compare_umi_to_pool, it)
    return(dists)

def plot_hdist(readset, outpng, seqs = None):
    if seqs == None:
        seqs = readset.umi_list()
    dists = hdist_all(seqs)
    chunks_iterator = itertools.zip_longest(*[iter(dists)]*len(seqs), fillvalue=None)
    plt.ioff()
    matplotlib.rcParams['figure.figsize'] = [20,20]
    hdist_mat = np.array([i for i in chunks_iterator])

    #TODO understand why is hdist_mat has an additional layer
    plt.imshow(hdist_mat[0],cmap="plasma" )
    plt.colorbar()
    plt.savefig(outpng, dpi=300)
    plt.clf()

def chain_keys(k, umis_dict):
    if k != umis_dict[k]:
        return(chain_keys(umis_dict[k], umis_dict))
    else:
        return(umis_dict[k])
    print("updating dict") 


def seq_to_fasta(seq, name=None):
    return(">%s\n%s\n"%(name,seq))

def write_to_fasta(fasta_ofn, seqs):
    with open(fasta_ofn, 'w') as fa:
        for name,seq in seqs:
            fa.write(seq_to_fasta(seq, name)) 


## STAR 

def star_build_ref(conf_fn, force = False, prepare_workspace = False):
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    if prepare_workspace:
        conf_init_workspace(conf)
    genomedir = conf['star']['genomedir'] 
    ref_fa = conf['input_files']['ref_fa']
    if not os.path.exists(genomedir) or force:
        os.makedirs(genomedir)
        ref_cmd = "STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN 100 --genomeSAindexNbases 2" %(genomedir, ref_fa)
        print(ref_cmd)
        subprocess.run(ref_cmd.split())
    else:
        print("Found pre-existing genomedir: %s" % (genomedir))

def star_map(conf_fn, fq_reads, prefix = '',  force = False, prepare_workspace = False):
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader) 
    if prepare_workspace:
        conf_init_workspace(conf)
    star = conf['star']
    genomedir = star['genomedir'] 

    if not os.path.exists(genomedir):
        star_build_ref(conf_fn, force)

    star_args = (genomedir, fq_reads, star['threads'], prefix, star['options'] )
    star_cmd = "STAR --genomeDir %s --readFilesIn %s --runThreadN %s --outFileNamePrefix %s %s" % star_args
    sort_cmd = "samtools sort %s -o %s" % (prefix+'Aligned.out.bam', prefix+"sorted.bam")
    index_cmd = "samtools index %s" % (prefix+"sorted.bam")
    
    print("searching for %s" %(prefix+"Aligned.out.bam"))
    print(os.path.exists(prefix+"Aligned.out.bam"))
    if (not os.path.exists(prefix+"Aligned.out.bam")) or force:
        print(star_cmd)
        subprocess.run(star_cmd.split())
        print(sort_cmd)
        subprocess.run(sort_cmd.split())
        print(index_cmd)
        subprocess.run(index_cmd.split())
    else:
        print("Found pre-existing aligment: %s" %(prefix+"Aligned.out.bam"))

 
def conf_return_path(conf, ftype = 'bam'):
    return("/".join([conf['dataset']['workdir'], conf['input_files'][ftype]]))

def conf_init_workspace(conf):
    workdir = conf['dataset']['workdir']
    if not os.path.exists(workdir):
        os.makedirs(workdir)

def clean_fasta(fn, prefix):
    fa = Fasta(fn)
    tmp_fa_fn = prefix+'tmp.fa'
    tmp_fa = open(tmp_fa_fn, 'w')
    for entry in fa.keys():
        tmp_fa.write(">%s\n%s\n"%(entry, fa[entry][:].seq.replace('\n', '')))
    fa.close()
    tmp_fa.close()
    cmd = "mv %s %s"%(tmp_fa_fn, fn)
    subprocess.run(cmd.split())



def fa_to_tabular(ifn, oufn, start, end):
    fa = Fasta(ifn)
    fout = open(oufn, 'w')
    fout.write("\t".join(["mm","hash","umi", "counts","seq"]) + '\n')
    for i in fa:
        name = i.name.split("_")
        mm = name[0]
        hash_key = name[1]
        umi = name[2]
        counts = name [3]
        seq = i[:].seq[start:end]
        fout.write("\t".join([mm, hash_key, umi, counts, seq])+'\n')
    fout.close()
    fa.close()



def fa_to_tabular(ifn, oufn, start, end):
    fa = Fasta(ifn)
    fout = open(oufn, 'w')
    fout.write("\t".join(["mm","hash","umi", "counts","seq"]) + '\n')
    for i in fa:
        name = i.name.split("_")
        mm = name[0]
        hash_key = name[1]
        umi = name[2]
        counts = name [3]
        seq = i[:].seq[start:end]
        fout.write("\t".join([mm, hash_key, umi, counts, seq])+'\n')
    fout.close()
    fa.close()


    
def filter_qual(query, tolerance = 0, min_qual = "@"):
    bad_nts = sum([i for i in map(lambda x: x.encode("ascii") <= min_qual.encode("ascii"), query.strip())])
    return(bad_nts < tolerance)

def trim_by_regex(dic, pattern, end = None, span_index = 0):
    import regex
    dic = copy.deepcopy(dic)
    edge_cases = []
    mean_start = []
    for i in dic.keys():
        match = regex.search(pattern, dic[i])
        if match:
            sp = match.span()
            if end == None:
                dic[i] = dic[i][sp[span_index]:]
            else:
                dic[i] = dic[i][sp[span_index]:sp[span_index]+end]
            mean_start.append(sp[span_index])
        else:
            edge_cases.append(i)
    if len(mean_start) >0:
        mean_start = int(np.mean(mean_start)) 
        for i in edge_cases:
            if end == None:
                dic[i] = dic[i][mean_start:]
            else:
                dic[i] = dic[i][mean_start:mean_start+end]
    return(dic)

    
def return_cons(params):
    import pysam
    entry, bam_tmp, min_cov = params
    fields = entry.split("_")
    umi = fields[1]
    cigar = fields[2]
    hits = int(fields[0])
    cons_bam = pysam.AlignmentFile(bam_tmp, 'rb')
    indel_pattern = regex.compile("([-\+])[0-9]+([ACGTNacgtn]+)")
    if hits >= min_cov:
        # at this point nts_per_pos is a vector of vectors (columns of the pileup) which in turn hold the list of nucleotides found at each point 
        nts_per_pos = [pos.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True) for pos in cons_bam.pileup(contig=entry, min_base_quality=0)]
        nts_per_pos = [[decide(j, indel_pattern) for j in i] for i in nts_per_pos]
        if len(nts_per_pos) >0:
            nts_per_pos = [pd.value_counts(i) for i in nts_per_pos]
            # danger: we are almost blindly assuming that the most frequent is the right one
            # here is where the consensus sequence is determined
            consensus = "".join([i.index[0] for i in nts_per_pos if len(i.index) >0])
            return([entry, consensus])
        else:
            return([None, None])
 
    else:
        return([None, None])


def decide(pattern, ins_pattern):
    match = ins_pattern.search(pattern)
    convs = {'*':'-', '>':'-', '<':'-'}
    if match:
        if match.groups()[0] == "-":
            return('')
        else:
            return(match.groups()[1])
    else:
        return(convs.get(pattern, pattern))


class Bint():
    def __init__(self, name):
        self.name = name
        self.mis = []
        self.dele = []
        self.ins = []

    def report(self):
        print(self.name)
        print(self.mis)
        print(self.dele)
        print(self.ins)


def return_bint_index(bint_l, current_position):
    '''Helper function that returns the key for the bint (barcode interval) that covers a given position on a pairwise alignment'''
    bint_index = None
    keys = [i for i in bint_l.keys()]
    for i in range(len(keys)):
        key  = keys[i]
        # Magical number warning = 27bp TODO:think twice
        if i == (len(keys)-1):
            if current_position >= bint_l[keys[i]][0] and current_position <= bint_l[keys[i]][0] + 27:
                return(key)
        else:
            if current_position >= bint_l[keys[i]][0] and current_position <= bint_l[keys[i+1]][0]:
                bint_index = key
                return(bint_index)
    return(0) # means that it couldn't find an overlapping bint

def get_lineage_vector(args):
    ref_seq, read_seq = args
    bint_db = yaml.load(open('/local/users/polivar/src/projects/mltracer/conf/gstlt_barcode_intervs.yaml'), Loader=yaml.FullLoader)
    ins = 0
    weird = 0
    mis = 0
    de = 0
    bint_dic = {}
    bint_keys = [i for i in bint_db['intervs'].keys()]
    for bint in bint_keys:
        bint_dic[bint] = Bint(name=bint)

    for i, (ref, read) in enumerate(zip(ref_seq, read_seq)):
        ref = ref.upper()
        read = read.upper()
        cur_bint = return_bint_index(bint_db['intervs'], i)
        if cur_bint != 0:
            if ref == "-" and read == "-":
                weird+=1
            else:
                if ref == "-":
                    ins+=1
                    # update the bint db by pushing down the coords by 1
                    #bint_db['intervs'][cur_bint][1] += 1
                    for j in range(bint_keys.index(cur_bint), len(bint_keys)):
                        next_bint = bint_keys[j]
                        bint_db['intervs'][next_bint][0]+= 1
                    # save integration
                    bint_dic[cur_bint].ins.append("i"+str(i)+read)
                if read == "-":
                    de+=1
                    bint_dic[cur_bint].dele.append("d"+str(i)+ref)
                if read != ref and ref != '-' and read != '-':
                    # we have a mismatch
                     mis+=1
                     print(cur_bint)
                     bint_dic[cur_bint].mis.append(ref+str(i)+read)
             #print(ref+"="+read, end = ' ')
    print("ins:",ins, "del:", de, "weird:", weird, "mis:", mis)
    for i in bint_dic.values():
        i.report()

def rev_comp(seq):
    " rev-comp a dna sequence with UIPAC characters ; courtesy of Max's crispor"
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M',
    "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n", "V" : "B", "v":"b", 
    "B" : "V", "b": "v", "W" : "W", "w" : "w", "-":"-"}

    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)



# for i in range(0, 50*(2), 2):
#     ref_seq  = FF[i+1][:].upper()
#     read_seq = FF[i][:].upper()
#     get_lineage_vector((ref_seq, read_seq))
#     if i == 26:
#         print(i,ref_seq)
#         print(i,read_seq)
#         pdb.set_trace()


