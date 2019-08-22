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


class Read_Set:
    def __init__(self, name = ''):
        self.umis = {}
        self.name = name
        self.ori_umis = None
        self.corr_stages = []

        ## mapped attributes
        self.original_header = None
        self.header = None

        self.consensus = None
        self.conf = None

    def set_header(self, header):
        self.original_header = header
        self.header = header.as_dict()


    def write_fasta(self, outfn  = None):
        prefix = self.conf['dataset']['workdir']+ "/" 
        if outfn == None:
            outfn = prefix+ '.cons.fa'
        consensuses_fa = open(outfn, 'w')
        for entry in self.consensus.keys():
            if len(self.consensus[entry]) >0:
                consensuses_fa.write(">%s\n%s\n"%(entry, self.consensus[entry].replace('\n', '')))
        print('save consensuses to %s' % (outfn))
        consensuses_fa.close()

    def add_umi(self, umi, seq, alignment=None):
        if self.umis.get(umi, 0) == 0:
            self.umis[umi] = UM(self.name, umi)
        self.umis[umi].seqs.append(seq)
        if alignment != None:
            if self.umis[umi].sams.get(alignment.cigarstring, 0) == 0:
                self.umis[umi].sams[alignment.cigarstring] = []
            self.umis[umi].sams[alignment.cigarstring].append(alignment)
            self.umis[umi].is_mapped = True

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

    def correct_umis(self, errors = 1, mode = "dist"):
        if self.ori_umis == None:
            self.ori_umis = [i for i in self.umis.keys()]
        self.corr_stages.append(self.umi_counts())
        seqs = [len(self.umis[i].seqs) for i in self.umis.keys()]
        by_counts = pd.Series(seqs, index = [i for i in self.umis.keys()]).sort_values()
        print("start", len(seqs))

        old_to_new = merge_all(self.umi_list(), errors = errors, mode = mode)
        for oldk,newk in old_to_new[::-1]:
            if oldk != newk:
                for i in self.umis[oldk].seqs:
                    self.umis[newk].seqs.append(copy.deepcopy(i))
                del self.umis[oldk]
        
        corrected_count = sum([1 for i,v in old_to_new if i !=v])
        print("Found %s corrections" %(corrected_count))
        
        print("end", len([i for i in self.umis.keys()]))
        if corrected_count>0:
            self.correct_umis(errors = errors, mode = mode)

    def do_msa(self):
        assert self.consensus, "consensus has not been generated: run .do_pileup to generate it"
        prefix = self.conf['dataset']['workdir']+ "/"  
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

def fastq_collapse_UMI(fastq_ifn, umi_start=0, umi_len=12, end=None):
    '''
    Process a fastq file and collapse reads by UMI making use of their consensus seq
    '''
    fq = open(fastq_ifn)
    it = itertools.islice(fq, 1, end, 4)
    readset = Read_Set(name=fastq_ifn)
    for i,read in enumerate(it):
        read = read.strip()
        umi = read[umi_start:umi_len] 
        seq = read[umi_len:]
        readset.add_umi(umi, seq)
    print(i, "reads were processed")
    return(readset)

def bam_collapse_UMI(bam_ifn, umi_start=0, umi_len=12, end=None):
    '''
    Process a fastq file and collapse reads by UMI making use of their consensus seq
    '''
    bamfile = pysam.pysam.AlignmentFile(bam_ifn)  
    readset = Read_Set(name = bam_ifn)
    readset.set_header(bamfile.header.copy())
    if end != None:
        end = int(end)
    bam_it = itertools.islice(bamfile, 0 ,end)
    for i,bam in enumerate(bam_it):
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

 
 
def do_pileup(readset, fa_ref = 'refs/rsa_plasmid.fa', start= 0 , end = None, region = 'chr_rsa:1898-2036', threads = 50, prefix = ''):
    self = readset
    if self.consensus != None:
        return(None)
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
    self.consensus = {}
    
    print("consensusing")
    print("this new bam file contains these chroms", len([i['SN'] for i in cons_bam.header.as_dict()['SQ']]))

    sorted_entries = sorted([i['SN'] for i in cons_bam.header.as_dict()['SQ']][1:], key = lambda x: x.split("_")[0])
    
    for i,entry in enumerate(sorted_entries):
        if i%100 ==0:
            print('\r', i, end ='')
        fields = entry.split("_")
        umi = fields[1]
        cigar = fields[2]
        hits = int(fields[0])
        #print(entry)
        if hits > 1:
            nts_per_pos = [pos.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True) for pos in cons_bam.pileup(contig=entry)]
            #print(len(nts_per_pos))
            cuini = [1 for i in nts_per_pos if regex.search("-|\+", "".join(i)) ]

            #if sum(cuini) >0:
            #    print(sum(cuini))
            #    pdb.set_trace()
            convs = {'*':'-', '>':'-', '<':'-'}
            indel_pattern = regex.compile("([-\+])[0-9]+([ACGTNacgtn]+)")

            def decide(pattern, ins_pattern):
                match = ins_pattern.search(pattern)
                if match:
                    if match.groups()[0] == "-":
                        return('')
                    else:
                        return(match.groups()[1])
                else:
                    return(convs.get(pattern, pattern))
                
            nts_per_pos = [[decide(j, indel_pattern) for j in i] for i in nts_per_pos]
            if(entry == "AAACGTTCACAC_15S18M4I99M15S_1"):
                pdb.set_trace()
            #print(entry)

            #pdb.set_trace()
            #print(len(nts_per_pos))
            if len(nts_per_pos) >0:
                nts_per_pos = [pd.value_counts(i) for i in nts_per_pos]
                #print(nts_per_pos)
                # danger: we are almost blindly assuming that the most frequent is the right one
                consensus = "".join([i.index[0] for i in nts_per_pos])
                #consensus = "".join([i[0] for i in nts_per_pos])
                #print(consensus)
                self.consensus[entry] = consensus[11:] 
        else:
            consensus = readset.umis[umi].sams[cigar][0].seq[11:]
            self.consensus[entry] = consensus 

    self.write_fasta()
    return(None)



    
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

def split_fastq(sample_name, out_prefix, f1_fn, f2_fn, vsamples, anchor, nreads = None):
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
                dic_files1[sample_id].write("\n".join([id[0]+ " sample_id:%s"%(sample_id), seq[0], plus[0], qual[0]]) +'\n')
                dic_files2[sample_id].write("\n".join([id[1]+ " sample_id:%s"%(sample_id), seq[1], plus[1], qual[1]]) +'\n')
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

def merge_all(seqs, jobs = 100, errors = 1, mode = "regex"):
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [[i, seqs, errors, mode] for i,v in enumerate(seqs)], fillvalue=seqs)
    return(pool.map(merge_umi_to_pool, it))

def hdist_all(seqs, jobs = 100):
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [[i, seqs] for i,v in enumerate(seqs)], fillvalue=seqs)
    dists = pool.map(compare_umi_to_pool, it)
    return(dists)

def plot_hdist(readset):
    seqs = readset.umi_list()
    dists = hdist_all(seqs)
    chunks_iterator = itertools.zip_longest(*[iter(dists)]*len(seqs), fillvalue=None)
    hdist_mat = np.array([i for i in chunks_iterator])
    matplotlib.rcParams['figure.figsize'] = [10,10]
    #TODO understand why is hdist_mat has an additional layer
    plt.imshow(hdist_mat[0],cmap="plasma" )
    plt.colorbar()

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


def star_build_ref(conf_fn, force = False):
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
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

def star_map(conf_fn, force = False):
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader) 
    conf_init_workspace(conf)
    name = conf['dataset']['name']
    bam_ofn = conf_return_path(conf, 'bam')
    star = conf['star']
    genomedir = star['genomedir'] 
    prefix = conf['dataset']['workdir']+"/" 

    if not os.path.exists(genomedir):
        star_build_ref(conf_fn, force)

    star_args = (genomedir, conf['input_files']['raw_fastq'], star['threads'], prefix, star['options'] )
    star_cmd = "STAR --genomeDir %s --readFilesIn %s --runThreadN %s --outFileNamePrefix %s %s" % star_args
    sort_cmd = "samtools sort %s -o %s" % (prefix+'Aligned.out.bam', prefix+"sorted.bam")
    index_cmd = "samtools index %s" % (prefix+"sorted.bam")

    if not os.path.exists(prefix+"Aligned.out.bam") or force:
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


