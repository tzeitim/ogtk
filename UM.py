import os
import operator
import copy
import regex
import itertools
import multiprocessing
import pandas as pd

class Read_Set:
    def __init__(self, name = ''):
        self.umis = {}
        self.name = name
        self.ori_umis = None
        self.corr_stages = []

    def add_umi(self, umi, seq):
        if self.umis.get(umi, 0) == 0:
            self.umis[umi] = UM(self.name, umi)
        self.umis[umi].seqs.append(seq)
    def umi_list(self):
        return([i[1] for i in self.umi_counts()])
    def umi_counts(self):
        by_counts = pd.Series([len(self.umis[i].seqs) for i in self.umis.keys()], index=self.umis.keys()).sort_values(ascending=False)
        by_counts = pd.DataFrame({"counts":by_counts, "umi":by_counts.index}, index=by_counts.index)
        by_counts = by_counts.sort_values(by=["counts", "umi"], ascending = False)
        return(by_counts['counts'])

    def umi_countss(self):
        by_counts = sorted([(len(self.umis[i].seqs), i) for i in self.umis.keys()], key = lambda x: x[1])
        
        return(by_counts)

    def correct_umis(self, errors = 1):
        if self.ori_umis == None:
            self.ori_umis = [i for i in self.umis.keys()]
        self.corr_stages.append([i for i in self.umis.keys()])
        seqs = [len(self.umis[i].seqs) for i in self.umis.keys()]
        by_counts = pd.Series(seqs, index = [i for i in self.umis.keys()]).sort_values()
        print("start", len(seqs))

        new_to_old = merge_all()
        import pdb
        pdb.set_trace()
        corrected_count = 0
        #umis_corr = {}
        #for i,seq in enumerate(by_counts.index):
        #    uncorrected = True
        #    print('\r%.4f\t%s' %(i/len(by_counts), seq), end='')
        #    for best in by_counts.index[::-1]:
        #        if uncorrected and by_counts[seq]<by_counts[best]:
        #            #match = regex.search("(%s){s<=%s}" % (seq, errors), best)
        #            dist = string_hamming_distance(seq, best)
        #            #if match:
        #            if dist <= errors:
        #                #if sum(match.fuzzy_counts) > 0:
        #                umis_corr[seq] = best
        #                corrected_count+=1
        #                uncorrected = False
        #    if uncorrected:
        #        umis_corr[seq] = seq

        #def chain_keys(k, umis_dict):
        #    if k != umis_dict[k]:
        #        return(chain_keys(umis_dict[k], umis_dict))
        #    else:
        #        return(umis_dict[k])
        #print("updating dict") 

        #for old_key in by_counts.index:
        #    new_key = chain_keys(old_key, umis_corr)
        #    if old_key != new_key:
        #        for i in self.umis[old_key].seqs:
        #            self.umis[new_key].seqs.append(copy.deepcopy(i))
        #        del self.umis[old_key]
        #print("...") 

        print("end", len([i for i in self.umis.keys()]))
        if corrected_count>0:
            self.correct_umis(errors = errors)


        
 
class UM:
    def __init__(self, sample, umi):
        self.sample = sample
        self.umi = umi
        self.seqs = []
        self.ori_umis = None

    def add_seq(self, seq):
        self.seqs.append(seq)



def fastq_collapse_UMI(fastq_ifn, umi_start=0, umi_len=12, end=None):
    '''
    Process a fastq file and collapse reads by UMI making use of their consensus seq
    '''
    fq = open(fastq_ifn)
    ifq = itertools.islice(fq, 1, end, 4)
    umis = Read_Set(name=fastq_ifn)
    for i,read in enumerate(ifq):
        read = read.strip()
        umi = read[umi_start:umi_len] 
        seq = read[umi_len-1:]
        umis.add_umi(umi, seq)
    print(i, "reads")
    return(umis)    
    
def bam_collapse_UMI():
    print("TODO")

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
                inv1.write("".join([id[0]+ " sample_id:%s"%(sample_id), seq[0], plus[0], qual[0]]) +'\n')
                inv2.write("".join([id[1]+ " sample_id:%s"%(sample_id), seq[1], plus[1], qual[1]]) +'\n')
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
    pool = args[1]
    errors = 1
    
    #TODO think on how one can reduce the number of cases (comparions, e.g. subset of pool)
    for best in pool:
        dist = string_hamming_distance((seq, best))
        if dist <= errors:
            return([seq, best])
    return([seq, best])

def merge_all(seqs, jobs = 100):
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [seqs], fillvalue=seqs)
    return(pool.map(merge_umi_to_pool, it))

