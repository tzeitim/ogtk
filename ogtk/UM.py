import os
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
import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from .aln_to_mat import Alg as Alg
import random
from .utils import *
from .mp import *
import gzip


class Read_Set:
    '''
    Hosts a list UM instances and provides a number of functions to process them
    '''
    def __init__(self, name = '', scurve_levels = None, scurve_step = None):
        # gets filled with UM instances
        self.umis = {}
        self.name = name
        self.ori_umis = None
        self.corr_stages = []
        self.nreads = 0

        # RID to UMI dic
        self.rid_umi = {} # holds a k,v for every readid, UMI
        ## mapped attributes
        self.original_header = None
        self.header = None

        self.consensus = None
        self.conf = None
        self.prefix = ''

        if scurve_levels is None:
            scurve_levels = [0, 2, 4, 8, 16]
        if scurve_step is None:
            scurve_step = 10000

        self.scurve_levels = scurve_levels
        self.scurve = {}
        for i in self.scurve_levels:
            self.scurve[i] = [(0, 0)]
        if scurve_step >= 1e4:
            self.scurve_step = int(1e4 * int(scurve_step/1e4))
        else:
            self.scurve_step = int(scurve_step)


    def log_scurve(self):
        counts = self.umi_counts()
        for i in self.scurve_levels:
            self.scurve[i].append((self.nreads, np.sum(counts>=i)))


    def saturation(self):
        return(' '.join([ str(sum(self.umi_counts() >= mincov)) for mincov in [1, 2, 3, 4, 5, 10, 100, 1000] ]))

    def true_allele_cov(self):
        ''' returns the counts for the top allele'''
        cc = []
        for umi in self.umis.keys():
            counts = pd.Series(self.umis[umi].seqs).value_counts()
            cc.append(counts[0])
        return(cc)

    def allele_saturation(self):
        return(' '.join([ str(sum([i >= mincov for i in self.true_allele_cov()])) for mincov in [1, 2, 3, 4, 5, 10, 100, 1000] ]))

    def set_header(self, header):
        self.original_header = header
        self.header = header.as_dict()

    def return_rep_reads(self, mincov = 1, unique = False):
        '''return the most frequent read for every umi in the readset'''
        rep_reads = [i.return_most_freq_read(mincov = mincov) for i in self.umis.values()]
        if not unique:
            return([i for i in rep_reads if i])
        else:
            return(list(set([i for i in rep_reads if i])))
        

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
    def write_consensuses_fasta(self, fasta_fn):
        ''' write to disk a tabulated file for each umi and its dominant sequence
        '''
        # TODO add assertion of correction
        with open(fasta_fn, 'w') as outf:
            for umi_str in self.umis.keys():
                umi = self.umis[umi_str]
                umi_seq = umi.consensus['dominant_seq']
                outf.write(">{}\n{}\n".format(umi_str, umi_seq))

    def return_consensuses_tuple(self):
        ''' returns tuples for each umi and its dominant sequence
            (umi_str, umi_seq, umi_dom_reads, umi_reads)
        '''
        # TODO add assertion of correction
        umi_tuple = []
        for umi_str in self.umis.keys():
            umi = self.umis[umi_str]
            umi_seq = umi.consensus['dominant_seq']
            umi_dom_reads = umi.consensus['dominant_reads']
            umi_reads = umi.consensus['reads_per_umi']
            umi_tuple.append((umi_str, umi_seq, umi_dom_reads, umi_reads))
        return(umi_tuple)

    def add_umi(self, umi, seq, alignment = None, qual = None):
        if self.umis.get(umi, 0) == 0:
            self.umis[umi] = UM(self.name, umi)
        self.umis[umi].seqs.append(seq)
        if qual != None:
            self.umis[umi].quals.append(qual.strip())
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

    def umi_counts(self, mincov = 1):
        by_counts = pd.Series([len(self.umis[i].seqs) for i in self.umis.keys()], index=self.umis.keys()).sort_values(ascending=False)
        by_counts = pd.DataFrame({"counts":by_counts, "umi":by_counts.index}, index=by_counts.index)
        by_counts = by_counts.sort_values(by=["counts", "umi"], ascending = False)
        if(mincov >1):
            return(by_counts['counts'][by_counts['counts']>=mincov])
        else:
            return(by_counts['counts'])

    def umi_countss(self, mincov = 1):
        by_counts = sorted([(len(self.umis[i].seqs), i) for i in self.umis.keys()], key = lambda x: x[1])
        return(by_counts)

    def correct_umis(self, errors = 1, mode = "dist", jobs =100, silent = True):
        jobs = int(jobs)
        if self.ori_umis == None:
            self.ori_umis = [i for i in self.umis.keys()]
        self.corr_stages.append(self.umi_counts())
        seqs = [len(self.umis[i].seqs) for i in self.umis.keys()]
        by_counts = pd.Series(seqs, index = [i for i in self.umis.keys()]).sort_values()
        if not silent:
            print(f"start {len(seqs)} = {(len(seqs)**2)/1e6:0.2f}M comparisons")

        old_to_new = merge_all(self.umi_list(), errors = errors, mode = mode, jobs = jobs)
        for oldk,newk in old_to_new[::-1]:
            if oldk != newk:
                for i in self.umis[oldk].seqs:
                    self.umis[newk].seqs.append(copy.deepcopy(i))
                del self.umis[oldk]
        
        corrected_count = sum([1 for i,v in old_to_new if i !=v])
        if not silent:
            print("Found %s corrections" %(corrected_count))
            print(f"Total corrected UMIs {len(list(self.umis.keys()))}")
        
        #if corrected_count>0:
        #    self.correct_umis(errors = errors, mode = mode)
    def delete(self, delete_candidates):
        '''
        deletes the keys (UMIS) specified on the argument
        '''
        for item in delete_candidates:
            del self.umis[item]

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
       
    def assemble_umis(self, threads= 10):
        ''' parallel call to tadpole to assemble contigs for each UMI'''
        umis =list(self.umis.keys())
        with multiprocessing.Pool(threads) as pool:
            print(f"parallel tadpole invocation {threads}")
            cc = pool.map(wrap_assemble, iter([self.umis[i] for i in umis]))
            print("yacabe")
            #pool.close()
            #pool.join()
            return(cc)
            #for k,contigs in zip(umis,cc):
            #    self.umis[k].contigs = contigs

    def export_to_fastq(self):
        return(print(' Not implemented yet, transfer the standalone function here'))
    
    def compute_coverage_table(self, cell_umis):
        '''
        Basic data frame with  coverage stats per UMI:
        fields returned: \t ["cell", "umi", "reads", "seqs", "r1_reads"]
        '''
        assert self.consensus, "consensus has not been generated: e.g. run .do_pileup to generate it"

        if len(self.consensus) == 0:
            return()

        consensuses_dict = self.consensus
        cov_df = [["cell", "umi", "reads", "seqs", "r1_reads"]]
        for cell in list(cell_umis.keys()):
            for umi in cell_umis[cell]:
                if umi in consensuses_dict.keys():
                    reads = len(self.umis[umi].seqs)
                    seqs = len(np.unique(self.umis[umi].seqs))
                    r1_reads = consensuses_dict[umi][2]  
                    cov_df.append([cell, umi, reads, seqs, r1_reads])

        cov_df = pd.DataFrame(cov_df[1:], columns=cov_df[0])
        print('pipi')
        return(cov_df)


def rs_export_to_fastq(rs, outfn):
    with gzip.open(outfn, 'wt') as fqout:
        for umi in rs.umis.keys():
            oumi = rs.umis[umi]
            for (i,seq), qual in zip(enumerate(oumi.seqs), oumi.quals):
                str_out = "\n".join([f'@{umi}_{i}', seq, "+", qual])
                fqout.write(str_out + '\n')

def create_tmp_dir_per_user(category = 'tadpole'):
    '''creates and returns a tmpr in a user-specific manner '''
    username = os.environ['USER']
    tmp_dir = f"/tmp/{category}_{username}/"
    if not os.path.exists(tmp_dir):
        print(f"Created tmp dir {tmp_dir}")
        os.makedirs(tmp_dir)
    subprocess.run(f"chmod -R g+rw {tmp_dir}".split())
    return(tmp_dir)

def wrap_kmer_correct(name, umi, k = 60, mincontig='auto', mincountseed=1, verbose = False, keep_tmp = True, sge_jobs = False, temp_prefix = None):
    username = os.environ['USER']
    temp_dir = create_tmp_dir_per_user(category = 'tadpole')
    if temp_prefix is None:
        temp_prefix = f'{temp_dir}/tadpolecorr_{name}_'

    fasta_ofn = temp_prefix+"_"+umi.umi+"_temp.fastq"
    umi.export_to_fastq(fasta_ofn = fasta_ofn)

    if sge_jobs:
        cmds =  tadpole_correct_fastqs(fasta_ofn, out=fasta_ofn.replace('temp', 'contig').replace('fastq', 'fa'), return_seqs=True, k = k, verbose = verbose, mincontig=mincontig, mincountseed=mincountseed, return_cmd = True)
        return(cmds)
    else:
        cseqs = tadpole_correct_fastqs(fasta_ofn, out=fasta_ofn.replace('temp', 'contig').replace('fastq', 'fa'), return_seqs=True, k = k, verbose = verbose, mincontig=mincontig, mincountseed=mincountseed)
    if not keep_tmp:
        print(f'Removing temp dir {temp_dir}')
        subprocess.run("rm {temp_dir}".split())
    return(cseqs)

def wrap_assemble(name, umi, k = 31, mincontig='auto', mincountseed=1, verbose = False, keep_tmp = True):
    username = os.environ['USER']
    temp_dir= create_tmp_dir_per_user(category = 'tadpole')
    temp_prefix= f'{temp_dir}/tadpole_{name}_'
    if len(set(umi.seqs))>1:
       fasta_ofn = temp_prefix+"_"+umi.umi+"_temp.fastq"
       umi.export_to_fastq(fasta_ofn = fasta_ofn)
       contigs = tadpole_fastqs(fasta_ofn, out=fasta_ofn.replace('temp', 'contig').replace('fastq', 'fa'), return_contigs=True, k = k, verbose = verbose, mincontig=mincontig, mincountseed=mincountseed)
    else:
       contigs = False
    if not keep_tmp:
        print(f'Removing temp dir {temp_dir}')
        subprocess.run("rm {temp_dir}".split())
    return(contigs)

def kmer_corr_umi_list_sge(umi_list, rs, ofn, k = 60, monitor_rate = 100, verbose = False, temp_prefix = None):
    '''Provided a umi list and a readset it attempts to assemble the umis while keeping a pickled record '''
    # cseqs = corrected sequences
    cseqs = {}
    for i,umi_str in enumerate(umi_list):
        umi = rs.umis[umi_str]
        if i%50 == 0 or i ==0 or i ==5 :
            print(f'\r{100*i/len(umi_list):.4f}% out of {i}/{len(umi_list)}', end = '')
        cseqs[umi.umi] = wrap_kmer_correct(umi.umi, umi, k = k, verbose = verbose, sge_jobs = True, temp_prefix = temp_prefix)
    return(cseqs)


def monitor_kmer_corr_umi_list(umi_list, rs, pickle_ofn, k = 80, monitor_rate = 100, verbose = False, temp_prefix = None):
    '''Provided a umi list and a readset it attempts to assemble the umis while keeping a pickled record '''
    # cseqs = corrected sequences
    cseqs = {}
    for i,umi_str in enumerate(umi_list):
        umi = rs.umis[umi_str]
        cseqs[umi.umi] = wrap_kmer_correct(umi.umi, umi, k = k, verbose = verbose, sge_jobs = False, temp_prefix = temp_prefix)
        if len(cseqs) % monitor_rate == 0 or len(cseqs) == 5 or len(cseqs) == 1:
            counts = [len(i) for i in cseqs.values() if i]
            print(pd.Series(counts).value_counts())
            with open(pickle_ofn, 'wb') as pcklout:
                pickle.dump(cseqs, pcklout)
    return(cseqs)


def monitor_assembly_umi_list(umi_list, rs, pickle_ofn, k = 80, monitor_rate = 100, verbose = False):
    '''Provided a umi list and a readset it attempts to assemble the umis while keeping a pickled record '''
    contigs = {}
    for umi_str in umi_list:
        umi = rs.umis[umi_str]
        contigs[umi.umi] = wrap_assemble(umi.umi, umi, k = k, verbose = verbose)
        if len(contigs) % monitor_rate == 0 or len(contigs) == 5 or len(contigs) == 1:
            counts = [len(i) for i in contigs.values() if i]
            print(pd.Series(counts).value_counts())
            with open(pickle_ofn, 'wb') as pcklout:
                pickle.dump(contigs, pcklout)

def assemble_pickled(pickle_in, k = 80, verbose = False):
    print(pickle_in)

    rs = pickle.load(open(pickle_in, 'rb'))
    if rs.__class__.__name__ != "Read_Set":
        if len(rs)==2:
            rs = rs[1]

    ucounts = rs.umi_counts()
    ucounts = ucounts[ucounts >1]
    pickle_out = pickle_in.replace('readset', 'contigs')
    print(pickle_out)
    monitor_assembly_umi_list(umi_list = ucounts.index , rs = rs,
       pickle_ofn = pickle_out,  k = k, monitor_rate = 100, verbose = verbose)

    
def qc_compare_raw_contigs(raw_bam_path, contig_bam_path, region, left, right, patches_bed = None, png_out = None):

    df_raw = bam_coverage(raw_bam_path, region)
    df_contig = bam_coverage(contig_bam_path, region)
    
    f, ax = plt.subplots(dpi=100, figsize=(20,5))
   
    raw = ax.plot(df_raw['start'], df_raw['cov']/np.sum(df_raw['cov']), linewidth=1)
    cont = ax.plot(df_contig['start'], df_contig['cov']/np.sum(df_contig['cov']), linewidth=1)
    ax.set_xlim(left = left, right = right)
    
    if patches_bed != None:
        from matplotlib import patches
        features = pd.read_table(patches_bed, header = None, names= ["chrom", "start", "end", "name"])
        features['i'] = range(len(features.index))
        feat_patches = features[['i','start', 'end']].apply(lambda x: alternating_patches(x, baseline=-0.0001), axis =1)

        for rect in feat_patches:
          ax.add_patch(rect)


    ax.title.set_text(os.path.basename(raw_bam_path.replace('bam','')))
#    f.suptitle(os.path.basename(raw_bam_path.replace('bam','')), y = 1.05, fontsize=12)
    ax.grid(True, linewidth=0.51, linestyle=':')
    ax.legend(['raw', 'assembled'], bbox_to_anchor = (1.0, 0.5))
    box = ax.get_position()
    ax.set_position([0.3, 1.0, box.width*0.3, box.height])
    ax.set_xlabel('chrom coord')
    ax.set_ylabel('fraction')


    if png_out != None:
         plt.savefig(png_out, bbox_inches='tight', dpi=150) #, loc=(1.1, 0.5))
         plt.close()
    else:
        return((f, ax, raw, cont))
    
def qc_assembly_pickle(pickle_ifn, png_out = None):
    '''Perform QCs for an assembly pickled dictionary'''
    from scipy.stats import gaussian_kde


    df = import_assembly_dict(pickle_ifn)
    bc = os.path.basename(pickle_ifn)
    
    if len(df.index) <=1:
        print(f'file {pickle_ifn} seems to contain no data. its data frame only has {len(df.index)} rows')
        return(None)
    qc_plot = plt.subplots(1, 2, figsize=(15,5))
    fig, (ax1, ax2) = qc_plot
    fig.suptitle(bc)

    # Calculate the point density
    xy = np.vstack([df['cov'], df['length']])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()

    ax1.title.set_text('Contig length histogram')
    ax1.set_xlabel('Contig length (bp)')
    ax1.set_ylabel('Counts')
    ax1.grid(True, linewidth=0.51, linestyle=':')

    ax2.title.set_text('Cov vs Length')
    ax2.set_xlabel('Contig length (bp)')
    ax2.set_ylabel('(log) Coverage')
    ax2.grid(True, linewidth=0.51, linestyle=':')

    ax1.hist(sorted(df['length']))
    ax2.set_yscale('log')
    ax2.scatter(df['length'][idx], df['cov'][idx], c=z[idx], s=100, edgecolor=None)

    if png_out != None:
        plt.savefig(png_out, bbox_inches='tight', dpi=100) #, loc=(1.1, 0.5))
        plt.close()
    else:
        return(qc_plot)

def map_assembly_pickle(pickle_ifn, ref):
    '''loads a pickled contig dict, exports it to fa and maps it using minimap2 '''
    def contigdf_to_fasta(x):
        umi, i, seq = x
        return('\n'.join([f'>{umi}_{i}', seq]))

    fa_out_path = pickle_ifn.replace('pickle', 'fa')
    bam_out_path = pickle_ifn.replace('pickle', 'bam')

    print(fa_out_path)
    print(bam_out_path)

    df = import_assembly_dict(pickle_ifn)
    if len(df.index) >1:
        with open(fa_out_path, 'wt') as outfa: 
            for i in df[['umi','i','seq']].apply(lambda x: contigdf_to_fasta(x), axis =1 ).to_list():
                outfa.write(i+'\n')
        mp2_map(ref = ref, query = fa_out_path, options = '-a --cs --splice --sam-hit-only', outfn = bam_out_path, force = True)

def import_assembly_dict(pickle_ifn):
    pcklhdl = open(pickle_ifn, 'rb')
    contigs = pickle.load(pcklhdl)
    records = []
    for k,v in contigs.items():
        if v:
            for ii, i in enumerate(v):
                length = float(i[0].split(',')[1].split("=")[1])
                cov = float(i[0].split(',')[2].split("=")[1])
                gc = float(i[0].split(',')[3].split("=")[1])
                seq = i[1]
                records.append([k, ii, length , cov, gc, seq])
    df = pd.DataFrame.from_records(records, columns=["umi", "i", "length", "cov", "gc", "seq"])
    return(df)
                      
def subrs(name,umilist, rs, keep_intact = False):
    '''Returns a new readset object from a given readset with the option of keeping the original readset intact or to remove the subset '''        
    #TODO check for a more efficient way to not need to copy and delete when not keeping intact.
    subrs = Read_Set(name = name)
    for umi in umilist:
        oumi = rs.umis[umi]
        subrs.umis[oumi.umi] = copy.deepcopy(oumi)
        if not keep_intact:
            del oumi
    print(f'extracted {len(subrs.umis)}')
    return(subrs)
        

class UM:
    def __init__(self, sample, umi):
        self.sample = sample
        self.umi = umi
        self.seqs = []
        self.quals = []
        self.ori_umis = None
        self.is_mapped = False
        self.contigs = False
     
        # stores sam entries using a cigar string as key
        self.sams = {}

        self.consensus = None
        
    def add_seq(self, seq):
        self.seqs.append(seq)

    def export_to_fasta(self, fasta_ofn = "/home/polivar/test.fa"):
        if not self.is_mapped:
            with open(fasta_ofn, 'w') as out_fa:
                cc = pd.Series(self.seqs).value_counts()
                index = 0
                for count,seq in zip(cc, cc.index):
                    for ii in range(count): 
                        out_fa.write(">"+str(index)+"_"+self.umi+"\n" + seq+"\n")
                        index+=1
        else:
            seqs = [ [i, ] for i in self.sams.keys()]


    def export_to_fastq(self, fasta_ofn = "/home/polivar/test.fastq", verbose = False):
        if verbose:
            print(f'Saving to {fasta_ofn}')

        if not self.is_mapped or True: # TODO remove this True to the else can work - > implement the else
            with open(fasta_ofn, 'w') as out_fa:
                merge_identical = False #TODO: decide how to deal with the quality for identical sequences
                if merge_identical:
                    cc = pd.Series(self.seqs).value_counts()
                    index = 0
                    for count, seq in zip(cc, cc.index):
                        for ii in range(count): 
                            out_fa.write("@"+str(index)+"_"+self.umi+"\n" + seq+"\n+\n" + qual + "\n")
                            index+=1
                else:
                    cc = pd.Series(self.seqs).value_counts()
                    for index, (seq, qual) in enumerate(zip(self.seqs, self.quals)):
                        out_fa.write("@"+str(index)+"_"+self.umi+"\n" + seq + "\n+\n" + qual + "\n")
        else:
            seqs = [ [i, ] for i in self.sams.keys()]

    def return_most_freq_read(self, mincov = 1):
        if not self.is_mapped:
            return(False) # TODO change this assertion to raise an error
        if len(self.sams.keys()) == 1 and len(self.seqs) >= mincov :
            return(list(self.sams.keys())[0])
        else:
            return(False) # TODO change this assertion to raise an error
        #counted_reads = pd.Series(self.seqs).value_counts()
        #counted_reads = counted_reads[counted_reads >= mincov]
        #if len(counted_reads) == 0:
        #    return(False)
        #top_read = counted_reads.index[0]
        #return(top_read)
    def return_kmer_counts(self, k):
        kmer_counts = {}
        for seq in self.seqs:
            for mer in build_kmers(seq, k):
                if mer not in kmer_counts.keys():
                    kmer_counts[mer] = 0
                kmer_counts[mer]+=1
        return(kmer_counts)
            

def count_umis_fastq(ifn_fastq):
    ''' memory efficient readset that only keeps record of umis without storing any other information of the read'''    
def count_umis_bam10x(ifn_bam):
    ''' memory efficient readset that only keeps record of umis without storing any other information of the read'''    



def pfastq_collapse_UMI(fastq_ifn1, fastq_ifn2, umi_start=0, umi_len=17, end=None, fuse = False, do_rc = False, min_reads_per_umi=1, **kwargs):
    '''Constructs a paired Readset by transfering the UMI from R1 to R2 via the
    read id. TODO and converts R2 to the same strand''' 
    
    rset1 = fastq_collapse_UMI(fastq_ifn1, umi_len = umi_len, keep_rid = True, end = end, min_reads_per_umi= min_reads_per_umi, scurve_levels = kwargs.get('scurve_levels', None), scurve_step = kwargs.get('scurve_step', None))

    rset2 = fastq_collapse_UMI(fastq_ifn2, umi_len = umi_len, rid_umi = rset1.rid_umi, do_rc = do_rc, end = end, min_reads_per_umi= min_reads_per_umi, scurve_levels = kwargs.get('scurve_levels', None), scurve_step = kwargs.get('scurve_step', None))

    #TODO add support for controlling the strandednes of the pooling method, for now this just dumps read2s in to read1s
    if fuse: 
        for k, v in rset1.umis.items():
            for read, qual in zip(rset2.umis[k].seqs, rset2.umis[k].quals):
                rset1.add_umi(k, seq = read, qual=qual)
        return(rset1)

    
    return([rset1, rset2])
    
def fastq_collapse_UMI(fastq_ifn, name = None, umi_start=0, umi_len=12, min_reads_per_umi =1, end=None, keep_rid=False, rid_umi=None, do_rc = False, downsample = False, threshold = 1, verbose = False, trimming_pattern = None, wl=None, **kwargs):

    '''
    Process a fastq file and collapses reads by UMI making use of their consensus seq
    Keeping the read id (rid) helps to match the UMI from read1 to read2
    If a rid_umi dic is provided, it uses it instead of the sequence information to call a UMI
    '''
    if end != None:
        end = int(end)
    if name == None:
        name = fastq_ifn
    umi_len = int(umi_len)
    readset = Read_Set(name=fastq_ifn, scurve_levels =  kwargs.get('scurve_levels', None), scurve_step =  kwargs.get('scurve_step', None))
    
    if end != None:
        end = end * 4

    fq = gzip.open(fastq_ifn, 'rt') if fastq_ifn.endswith("fastq.gz") else open(fastq_ifn)#.readlines()
    reads =     itertools.islice(fq, 1, end, 4)
    fq2 = gzip.open(fastq_ifn, 'rt') if fastq_ifn.endswith("fastq.gz") else open(fastq_ifn)#.readlines()
    rids =      itertools.islice(fq2, 0, end, 4)
    fqq = gzip.open(fastq_ifn, 'rt') if fastq_ifn.endswith("fastq.gz") else open(fastq_ifn)#.readlines()
    quals =      itertools.islice(fqq, 3, end, 4)

    # TODO improve this routine, now it is hard coded for single patterns
    if trimming_pattern is not None:
        if not isinstance(trimming_pattern, dict):
            raise TypeError('Trimming patterns must be entered as a dict for each pattern to match (value) and group name to keep (key), for example `{"read":"ANCHOR1{e<=3}(?P<read>.+)ANCHOR2{e<=3}"}`\ntrimming_pattern = {"read":"(ANCHOR1){e<=3}(?P<read>.+)(ANCHOR2){e<=3}"}')
        keep_group = list(trimming_pattern.keys())[0] # e.g keeps the 'read' group from the match object
        trim_re = regex.compile(list(trimming_pattern.values())[0])
        trim_matches = []

    if downsample:
        reads_rids_it = [(x, y, z) for x,y,z in zip(reads, rids, quals) if random.random() < threshold]
    else:
        reads_rids_it = [(x, y, z) for x,y,z in zip(reads, rids, quals)]

    if wl is not None:
        wl_skipped=0
        wl = set(wl)

    if 'match_pattern' in kwargs.keys():
        do_match = True
        match_pattern = regex.compile(kwargs['match_pattern'])
    else:
        do_match = False

    for i, (read, rid, qual) in enumerate(reads_rids_it):
        if do_match:
            if not match_pattern.search(read):
                continue
        readset.nreads += 1
        if readset.nreads%readset.scurve_step ==0:
            #readset.scurve.append((readset.nreads, len(readset.umis.keys()) ))
            readset.log_scurve()
        read = read.strip()
        if rid_umi == None:
            umi = read[umi_start:umi_len] 
            seq = read[umi_len:]
            qual = qual[umi_len:]
        else:
            umi = rid_umi[rid.split(" ")[0]]
            seq = read.strip()
            qual = qual.strip()
        if do_rc:
            seq = rev_comp(seq) 
            qual = qual[::-1]

        if trimming_pattern != None:
            match = trim_re.search(seq)
            if match:
                # we append the different anatomical components of the read and their qualities
                spans = [match.span(i) for i in match.groupdict().keys()]
                seq = "".join([seq[i:ii] for i,ii in spans])
                qual = "".join([qual[i:ii] for i,ii in spans])
                trim_matches.append(1)
            else:
                trim_matches.append(0)

        if wl is not None:
            if umi in wl:
                readset.add_umi(umi, seq, qual = qual)
            else:
                wl_skipped+=1
            continue

        readset.add_umi(umi, seq, qual = qual)

        if keep_rid:
            readset.add_rid(umi, rid.split(" ")[0])


    raw_rs_size = len(readset.umis)
    badly_covered = [umi.umi for umi in readset.umis.values() if len(umi.seqs) < min_reads_per_umi ]
    readset.delete(badly_covered)

    print(f'A total of {len(badly_covered)}::{100*len(badly_covered)/raw_rs_size:.2f}% umis were removed due to poor coverage')

    ucounts = readset.umi_counts()
    top10 = (100*ucounts[0:9]/readset.nreads).tolist()
    top10 = " ".join([i for i in map( lambda x: f'{x:.2f}%', top10[0:9])]) 
    sat = readset.saturation()
    #print(sat[1]/sat[0], sat[0]/sat[1])
    if verbose:
        print(f"A total of {readset.nreads} reads were processed with {len(ucounts)} umis. top10 {top10} ")
    if trimming_pattern != None:
        print(f'trimming stats {sum(trim_matches)}/{len(trim_matches)} {sum(trim_matches)/len(trim_matches):.2f}')
    if wl is not None:
        print(f"A total of {readset.nreads} reads were processed but {wl_skipped} were omitted final size ={len(ucounts)}")

    return(readset)

def bam_collapse_UMI(bam_ifn, umi_start=0, umi_len=12, end=None):
    '''
    Process a BAM file and collapse reads by UMI making use of their consensus seq
    '''
    print("processing %s" % bam_ifn)
    umi_len = int(umi_len)
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
        readset.add_umi(umi, seq, alignment=bam, qual = bam.qual)
        
    print(f"\r{i} reads were processed")
    return(readset)

def bam_from_merged_collapse_UMI(bam_ifn, end=None, region = None):
    '''
    Process merged BAM-file from a 10x-style read pairs that have been merged before with merge_10x_fastqs. Looks for ":umi:" as a separator. 
    '''
    print("processing %s" % bam_ifn)
    umi_len = 28
    bamfile = pysam.pysam.AlignmentFile(bam_ifn)  
    readset = Read_Set(name = bam_ifn)
    readset.set_header(bamfile.header.copy())

    if end != None:
        end = int(end)

    if region != None:
        bam_it = bamfile.fetch(region)
    else:
        bam_it = itertools.islice(bamfile, 0, end)

    for i,bam in enumerate(bam_it):
        print("\r%s"%i, end='')
        read = bam.seq
        umi = bam.query_name.split(":umi:")[1] 
        seq = bam.seq
        readset.add_umi(umi, seq, qual = bam.qual)
        
    print(f"\r{i} reads were processed")
    return(readset)

def bam10x_collapse_UMI(bam_ifn, umi_start=0, umi_len=28, end=None, use_uncorrected = False):
    '''
    Process a BAM file and collapse reads by UMI making use of their consensus seq
    '''
    print("processing %s" % bam_ifn)
    umi_len = int(umi_len)
    bamfile = pysam.pysam.AlignmentFile(bam_ifn)  
    readset = Read_Set(name = bam_ifn)
    readset.set_header(bamfile.header.copy())
    if end != None:
        end = int(end)
    bam_it = itertools.islice(bamfile, 0, end)

    for i,bam in enumerate(bam_it):
        if len(readset.umis)%10000 == 0:
            print("\r%s"%i, end='')

        if bam.has_tag("UB") and bam.has_tag("CB"):
            read = bam.seq
            seq = bam.seq
            umi = bam.get_tag("CB") + bam.get_tag("UB")
            readset.add_umi(umi, seq, alignment=bam, qual = bam.qual)

        elif use_uncorrected:
            if bam.has_tag("UR") and bam.has_tag("CR"):
                read = bam.seq
                seq = bam.seq
                umi = bam.get_tag("CR") + bam.get_tag("UR")
                readset.add_umi(umi, seq, alignment=bam, qual = bam.qual)

    print(f"\r{i} reads were processed")
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
    fname, name, seqs, min_cov, jobs, naive  = args
        
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
            read_len = len(msa_fa[0][:])
            read_counts = [int(i.split("_")[2]) for i in msa_fa.keys()]
            # extract the consensus
            ## naively
            if False:
                # concatenate aligned sequences (all the same size) 
                read_stream = ""
                for read in msa_fa.keys():
                    read_stream = read_stream + str(msa_fa[read][:])

                # chop the stream into a matrix and trasfer it to the pandas universe
                msa_mat = [i for i in read_stream]
                msa_mat = np.array(msa_mat) 
                msa_mat = pd.DataFrame(msa_mat.reshape((len(msa_fa.keys()),-1)))
                
                umi_consensus = ''
                for i in range(read_len):
                    variants = msa_mat[i]
                    # multiply each nt by the times the read was found, concatenate into string and transform into a vector
                    column_expansion = [i for i in ''.join([var*count for var,count in zip(variants, read_counts)])]
                    column_counts = pd.Series(column_expansion).value_counts()
                    umi_consensus = umi_consensus + column_counts[0].index.to_list()[0].upper()
                #umi_consensus = ''.join([msa_mat[i].value_counts().head(1).index.to_list()[0] for i in range(read_len)]).upper()
            else:
                expand_by_counts = True
                cons_in =  '{}_{}_cons_in.fa'.format(fname, name)
                cons_out = '{}_{}_cons_out.fa'.format(fname, name)
                cmd_cons = 'cons -datafile EDNAFULL -sequence {} -outseq {} -name {}'.format(cons_in, cons_out, name)

                cons_in_Fa = open(cons_in, 'w')
                for entry, times in zip(msa_fa.keys(), umi_counts.to_list()):
                    seq = msa_fa[entry][:]
                    if not expand_by_counts:
                        times = 1
                    for ins in range(times):
                        cons_in_Fa.write(">{}_{}\n{}\n".format(entry, times, seq)) 
                cons_in_Fa.close()
                # invoque emboss cons 
                pp = subprocess.run(cmd_cons.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
                cons_Fa = Fasta(cons_out, as_raw=True)
                umi_consensus = cons_Fa[0][:]
                umi_consensus = umi_consensus.replace('N', '')
                umi_consensus = umi_consensus.replace('n', '')

            consensus = (name, umi_consensus)
        else:
            if umi_counts[0] >= min_cov:
                consensus = (name, seqs[0])
            else: 
                #print("got only one allele and doesn't reach the min_cov {} with counts ".format(min_cov, umi_counts[0]))
                consensus = ('drop-no_cov', name)
    else:
        #print("got only one rep")
        #if int(name.split("_")[1]) > 1:
        consensus = (name, seqs[0])
        #else:
        #    consensus = ("drop-one-seq", name)
             
    return(consensus) 

# TODO: rename this function to something like return_consensuses
def do_fastq_pileup(readset, min_cov = 2, threads =  100, threads_alg = 1, trim_by = None, by_alignment = False):
    ''' Generate consensus sequences of sets of reads, grouped by UMI \n\n\n
    By aligment:\n\n
    collapsing reads by consensus seems to be incorrect, deprecating the use of such (e.g. mafft alignment) 
    in favor for just getting the most frequent sequence

    By rank:\n\n
    this option represents a ranked based approach where only the most common sequence is used as the representative\n
    it returns a dictionary of the form (umi, (dominant_seq, reads_per_umi, reads_rank1))

    '''
    # define the list of pool arguments
    # fname, name, seqs, min_cov, jobs  = args
    it = iter([(readset.name.replace('.fastq', ''), umi, readset.umis[umi].seqs, min_cov, threads_alg) for umi in readset.umis.keys()])
    if by_alignment:
        pool = multiprocessing.Pool(threads)
        print("Aligning read with the same UMI")
        cc = pool.map(mafft_consensus, it)
        pool.close()
        pool.join()
        #cc = [i for i in cc if "drop" not in i[0]]
        # TODO adapt to incorporate dictionary-based data aggregation
        return(consensuses) 
    else:
    # it returns a dictionary of the form (umi, (dominant_seq, reads_per_umi, reads_rank1))
        cc = []
        for umi in readset.umis.keys():
            counts = pd.Series(readset.umis[umi].seqs).value_counts()
            readset.umis[umi].consensus = {
                "type": "rank",
                "dominant_seq":counts.index.to_list()[0],
                "dominant_reads": counts[0],
                "reads_per_umi":len(readset.umis[umi].seqs),
                "pool_size": len(counts),
            }

   #         if counts[0] >= min_cov:
   #             dominant_seq = counts.index.to_list()[0]
   #             reads_per_umi = len(readset.umis[umi].seqs)
   #             reads_rank1 = counts[0]
   #             foc = (umi, (dominant_seq, reads_per_umi, reads_rank1))
   #  
    # legacy code - this function used to return a dictionary 
    #consensuses = trimdict_by_regex(dict(cc) , trim_by, span_index=1) if trim_by != None else dict(cc) 

 
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
    args = (seq, rank, pool, errors, mode)
    mode = "regex" "dist"
    '''
    seq = args[0]
    rank = args[1][0]
    pool = args[1][1]
    errors = args[1][2]
    mode = args[1][3]
   
    
    rc = regex.compile(seq+'{e<='+str(errors)+'}')
    #TODO think on how one can one reduce the number of cases (comparions, e.g. subset of pool)
    for i,best in enumerate(pool[0:rank+1]):
        if mode == "regex":
            match = rc.search(best)
            if match and i < rank:
                return([seq, best])
        if mode == "dist":
            dist = string_hamming_distance((seq, best))
            if dist <= errors and i < rank:
                return([seq, best])
    return([seq, seq])

def compare_umi_to_pool(args):
    '''
    Returns the hammin distance between a given UMI to the rest of the pool
    args = (seq, (index, pool))
    '''
    seq = args[0]
    index = args[1][0]
    pool = args[1][1]
   
    dists = [] 
    #TODO think on how one can one reduce the number of cases (comparions, e.g. subset of pool)
    for i,best in enumerate(pool):
        dists.append(string_hamming_distance((seq, best)))
    return(dists)

def merge_all(seqs, jobs = 10, errors = 1, mode = "dist"):
    ''' 
    Parallel wrapper for merge_umi_to_pool - expects seqs sorted by ocurrence 
    '''
    jobs = int(jobs)
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [[i, seqs, errors, mode] for i,v in enumerate(seqs)], fillvalue=seqs)
    ret = pool.map(merge_umi_to_pool, it)
    pool.close()
    pool.join()
    return(ret)

def hdist_all(seqs, jobs = 10):
    '''
    Pair-wise comparison of set of sequences
    '''
    jobs = int(jobs)
    pool = multiprocessing.Pool(jobs)
    it = itertools.zip_longest(seqs, [[i, seqs] for i,v in enumerate(seqs)], fillvalue=seqs)
    dists = pool.map(compare_umi_to_pool, it)
    pool.close()
    pool.join()
    chunks_iterator = itertools.zip_longest(*[iter(dists)]*len(seqs), fillvalue=None)
    hdist_mat = np.array([i for i in chunks_iterator])
    #TODO understand why is hdist_mat has an additional layer
    #TODO should a flag for as_matrix be added?
    return(hdist_mat[0])

def plot_hdist(readset = None, outpng = None, seqs = None):
    if outpng == None:
        raise ValueError('an outpng must be provided')
    if seqs == None and readset == None:
        raise ValueError('a list of sequences or a readset must be provided')
    if seqs == None:
        seqs = readset.umi_list()
    hdist_mat = hdist_all(seqs)
    plt.ioff()
    matplotlib.rcParams['figure.figsize'] = [20,20]

    #plt.imshow(hdist_mat[0],cmap="plasma" )
    plt.imshow(hdist_mat,cmap="plasma" )

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

def trimdict_by_regex(dic, pattern, end = None, span_index = 0):
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
    '''
    Helper function that returns the key for the bint (barcode interval) that covers a given position on a pairwise alignment
    '''
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


def alternating_patches(x, baseline =-0.0002, delta= 0.0001):
    ''' helper function to draw rectangles after applying it on a df usually via lambda'''
    from matplotlib import patches
    index, left, right = x
    sign = 1 if index %2 ==0 else -1
    delta = delta * sign    
    
    xcoords = [ [left, right] , [right, left]]
    ycoords = [ [baseline + delta, baseline + delta], [baseline, baseline]]
 
    xcoords = np.array(xcoords).flatten()
    ycoords = np.array(ycoords).flatten()
    coords = [i for i in zip(xcoords, ycoords)]
    return(patches.Polygon(np.array(coords), fill = False))


def get_list_of_10xcells(rs, correction_dict_path = None, correction_dict = None):
    '''
    Returns a CBC-> list(UMIs) dict of 10x-style cell barcodes. It can use a correction dictionary if provided, else it just strips the full UMI according to the range that corresponds to the CBC
    '''
    umi_list = rs.umis.keys()
    valid_cells = []
    valid_pairs = []
    if correction_dict_path is not None or correction_dict is not None:
        # UMI correction via correction dictionaries
        hits_dic = []
        print(f'Using cellranger dictionary', end = "...")
        if correction_dict_path is not None:
            print(f'Loading correction dictionaries {correction_dict_path}')
            correction_dict = pickle.load(open(correction_dict_path, 'rb'))
        if correction_dict is not None:
            print("using provided")

        for umi in umi_list:
            candidate_cell = umi[0:16]
            if candidate_cell in correction_dict['cell_barcodes'].keys():
                candidate_cell = correction_dict['cell_barcodes'][candidate_cell]
                hits_dic.append(1)
                # correct raw cellbarcode with the cellranger-corrected cellbarcode
                valid_cells.append(candidate_cell)
                valid_pairs.append((candidate_cell, umi))
            else:
                hits_dic.append(0)
        print(f'correction stats {sum(hits_dic)/len(hits_dic)}')
    else:
    # without correction
        for umi in umi_list:
            canditate_cell = umi[0:16]
            if not any([i == "N" for i in candidate_cell]):
                candidate_cell = candidate_cell +"-1"
                valid_cells.append(candidate_cell)
                valid_pairs.append((candidate_cell, umi))

    cell_umi_dict = {}

    for cell,umi in valid_pairs:
        if cell not in cell_umi_dict.keys():
            cell_umi_dict[cell] = []
        cell_umi_dict[cell].append(umi)

    return(cell_umi_dict)

