import shutil
import os
import subprocess
import glob
import itertools

#https://docs.python.org/3/library/itertools.html#itertools.zip_longest
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def bulk_merge_dir(raw_path, out_path, clean = True, force = False, errors =0, verbose = False):
    ''' This wrapper makes use of bbmerge to merge fastqs locate in a given directory '''

    if clean and os.path.exists(out_path):
        shutil.rmtree(out_path, ignore_errors=True)
        os.makedirs(out_path)
        os.makedirs(out_path+"/merged")
        os.makedirs(out_path+"/umerged")

    if not os.path.exists(out_path):
        os.makedirs(out_path)
        os.makedirs(out_path+"/merged")
        os.makedirs(out_path+"/umerged")

    ###
    returns = []
    for i in glob.glob(raw_path+"/*R1*fastq"):
        f1 = i
        f2 = i.replace('R1', 'R2')
        basename = f1.split('/')[-1]
        basename = basename.split('_R1')[0]

        out_merged = "{}/{}_merged.fastq".format(out_path, basename)
        out_unmerged1 = "{}/{}_unmerged_R1.fastq".format(out_path, basename)
        out_unmerged2 = "{}/{}_unmerged_R2.fastq".format(out_path, basename)
        unmerged_samples = []
        
        if not os.path.exists(out_merged) or force:
            ret = bbmerge_fastqs(f1 = f1, f2 = f2, out = out_merged, outu = out_unmerged1 , outu2 = out_unmerged2, verbose = verbose)
            returns.append(ret)
        if ret:
            unmerged_samples.append(basename)
    return(returns)

def bbmerge_fastqs(f1, f2, out, outu, outu2, verbose = False, bbmerge_bin = 'bbmerge.sh'):
    ''' use bbtools to merge amplicon derived reads 
    docs: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbmerge.sh
    ref: https://doi.org/10.1371/journal.pone.0185056
    If no present, install through conda bbtools
    '''
    cmd = "{} in1={} in2={} out={} outu={} outu2={} pfilter=0 mismatches=5 efilter=-1".format(bbmerge_bin, f1, f2, out, outu, outu2)
    pp = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if "RuntimeException" in pp.stderr.decode():
        #means the command failed
        return(False)
    else:
        total = ''.join([i for i in pp.stderr.decode().split('\n') if "Pairs" in i])
        total = float(total.split('\t')[-1].replace("%", ''))
        joined = ''.join([i for i in pp.stderr.decode().split('\n') if "Joined" in i])
        joined = float(joined.split('\t')[-1].replace("%", ''))/100
        log = '\n'.join([i for i in pp.stderr.decode().split('\n')])
        log_ofn =out.replace('fastq', 'log') 
        if verbose:
            print(pp.stderr.decode())
            print("bbmerge log %s" % (log_ofn))


        with open(log_ofn, 'w') as logout:
            for i in log:
                logout.write(i)

        return((total, joined, log_ofn))

def print_open_children():
    ''' prints the currently open sub-processes on ipython'''
    print(subprocess.Popen('ps -u polivar | grep ipython | wc -l', shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode())


def tadpole_fastqs(f1, out, verbose = False, k=66, threads =1, tadpole_bin = 'tadpole.sh', bm1=1, bm2=1, mincontig = "auto", mincountseed = 100, return_contigs=False):
    ''' use tadpole from bbtools to assemble a cloud of sequences controlled by a UMI
    '''
    #tadpole.sh in=tmp/f1_ACTTCGCCAGAGTTGG_GTGCGAGAGGGTA.fastq out=mini k=66 overwrite=True bm1=1 bm2=1 mincountseed=4
    cmd = f"{tadpole_bin} in={f1} out={out} k={k} overwrite=True bm1={bm1} bm2={bm2} t={threads} -Xmx6g mincontig={mincontig} mincountseed={mincountseed} rcomp=f"
    if verbose:
        print(cmd)
    pp = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if "RuntimeException" in pp.stderr.decode():
        #means the command failed
        print("RunTimeEx")
        return(False)
    contigs = pp.stderr.decode().split("Contigs generated:")[1].split('\n')[0].split('\t')[1]
    contigs = int(contigs)

    if contigs == 0:
        return(False)

    else:
        #total = ''.join([i for i in pp.stderr.decode().split('\n') if "Pairs" in i])
        #total = float(total.split('\t')[-1].replace("%", ''))
        #joined = ''.join([i for i in pp.stderr.decode().split('\n') if "Joined" in i])
        #joined = float(joined.split('\t')[-1].replace("%", ''))/100
        log = '\n'.join([i for i in pp.stderr.decode().split('\n')])
        if "contig" in out:
            log_ofn =out.replace('contig', 'log') 
        else:
            log_ofn = out+'.log'
        if verbose:
            print(pp.stderr.decode())
            print("tadpole log %s" % (log_ofn))


        with open(log_ofn, 'w') as logout:
            for i in log:
                logout.write(i)
        if return_contigs:
            from pyfaidx import Fasta
            contigs = Fasta(out, as_raw=True)
            contigs = [ (k, contigs[k][:]) for k in contigs.keys()]
            return(contigs)

        return(True)

def tadpole_correct_fastqs(f1, out, verbose = False, k=66, threads =1, tadpole_bin = 'tadpole.sh', bm1=1, bm2=1, mincontig = "auto", mincountseed = 100, return_seqs=True):
    ''' use tadpole from bbtools to error-correct a cloud of sequences controlled by a UMI
    '''
    #tadpole.sh in=tmp/f1_ACTTCGCCAGAGTTGG_GTGCGAGAGGGTA.fastq out=mini k=66 overwrite=True bm1=1 bm2=1 mincountseed=4
    cmd = f"{tadpole_bin} in={f1} out={out} k={k} mode=correct overwrite=True bm1={bm1} bm2={bm2} t={threads} -Xmx6g mincontig={mincontig} mincountseed={mincountseed} rcomp=f"
    if verbose:
        print(cmd)
    pp = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if "RuntimeException" in pp.stderr.decode():
        #means the command failed
        print("RunTimeEx")
        return(False)

    # corrected sequences - cseqs
    cseqs = pp.stderr.decode().split("Reads fully corrected:")[1].split('\n')[0].split('\t')[1]
    cseqs = int(cseqs)

    if cseqs == 0:
        return(False)

    else:
        log = '\n'.join([i for i in pp.stderr.decode().split('\n')])
        if "contig" in out:
            log_ofn =out.replace('contig', 'log') 
        else:
            log_ofn = out+'.log'
        if verbose:
            print(pp.stderr.decode())
            print("tadpole log %s" % (log_ofn))


        with open(log_ofn, 'w') as logout:
            for i in log:
                logout.write(i)

        if return_seqs:
            from pyfaidx import Fasta
            cseqs = Fasta(out, as_raw=True)
            cseqs = [ (k, cseqs[k][:]) for k in cseqs.keys()]
            return(cseqs)

        return(True)



def bam_coverage(bam_path, region):
    ''' returns a panda dataframe for a given region of a given bam file
        The hardcoded -d = 0 means to not cap the coverage (this makes this function very slow).
    '''
    import pysam
    import pandas as pd
    cached_fn =bam_path+".cached_cov.pd.pickle" 
    if os.path.exists(cached_fn):
        print(f'loading from cache {cached_fn}')
        df = pd.read_pickle(cached_fn)
        return(df)
        
    cov = pysam.samtools.depth("-a", "-d", "0", "-r", region, bam_path)
    df = pd.DataFrame([i.split('\t') for i in cov.split('\n')], columns = ['chr', 'start', 'cov'])
    df['cov'] = pd.to_numeric(df['cov'])
    df['start'] = pd.to_numeric(df['start'])

    df.to_pickle(cached_fn) 
    return(df)


def sfind(path, pattern, options = "-iname"):
    import subprocess

    sys_call = subprocess.run(f'find {path} {options} {pattern}'.split(), stdout =subprocess.PIPE, text = True)
    ls_call = [i for i in sys_call.stdout.split()] 
    return(ls_call)

def merge_10x_fastqs(indir, r1, end = None):
    ''' merge 10x fastqs into a single one (AKA R3). Uses ":umi:" as separator 
    It is recommended to run annotate_bam after mapping these files to transfer the CR and UR tags
    '''
    import gzip
    r2 = r1.replace("R1", "R2")
    r3 = r1.replace("R1", "R3")

    with gzip.open(indir + r1, 'rt') as R1, gzip.open(indir+ r2, 'rt') as R2, gzip.open(indir+r3, 'wt') as R3:
        i1 = itertools.islice(grouper(R1, 4), 0, end)
        i2 = itertools.islice(grouper(R2, 4), 0, end)
        for read1,read2 in zip(i1, i2):
            read3 = list(read2)
            umi_str = read1[1][0:28]
            read3[0] = read3[0].split()[0] + ":umi:"+ umi_str + "\n"#read1[1][0:28]
            for i in read3:
                R3.write(i)
            #print(f"read1{read1}\nread2{read2}\nread3{read3}")
    print(f"merged fastq {r3}")
    return(r3)

def annotate_bam(ibam):
    ''' Used to transfer the CR and UR tags on bam files generated by mapping merged fastqs'''
    import pysam
    obam_fn = ibam.replace(".bam", "_ann.bam")
    ibam = pysam.AlignmentFile(ibam)
    obam = pysam.AlignmentFile(obam_fn, 'wb', template = ibam)
    for i in ibam:
        molid = i.query_name.split(":umi:")[1]
        i.set_tag("UR", molid[16:], "Z")
        i.set_tag("CR", molid[0:16], "Z")
        i.set_tag("ZI", i.cigarstring, "Z")
        obam.write(i)
    ibam.close()
    obam.close()
    print(f'saved annotated bam file\n{obam_fn}')
    # call samtools for indexing
    cmd_index = f'samtools index {obam_fn}'
    subprocess.run(cmd_index.split())


def merge_10x_fastq_dir(indir, force=False, read1_pattern = "_R1_"):
    import os
    ''' merges all fastq pairs in a given directory'''
    for r1 in os.listdir(indir):
        if read1_pattern in r1 and "fastq.gz":
            print("/".join([indir,r1]))
            # define a mappe
            #mbam_path = indir + r1.replace("R1", "R3").replace("fastq.gz", "bam")
            
            if not os.path.exists((indir+r1).replace(read1_pattern, "_R3_")) or force:
                print(f"merging r1 and r2 as merged read didn't exist")
                r3 = merge_10x_fastqs(indir, r1)
            
            #if not os.path.exists(mbam_path):
            #    print(f"mapping to {mbam_path}")
            #    ogtk.mp.mp2_map(ref="/local/users/polivar/src/ivltr/refs/transgenes/trans.fa ", 
            #                query = indir + r3, 
            #                outfn = mbam_path, verbose = True, options = "-a --cs --splice")
            
            #ogtk.UM.annotate_bam(mbam_path)
            
            print(r3)
            #print(mbam_path)
     
def generate_correction_dictionaries(pickle_ofn, path_to_bam, force = False, verbose = False):
    '''
    Generates a cellbarcode correction dictionary based on the cellranger barcoded BAM tags
    TODO: there is a small fraction of uncorrected cell barcodes that seem to map to more than one corrected cell barcode
    '''
    import pysam
    import itertools
    import pickle 
    import pyaml
    import os

    # load a pre-computed dictionary of found
    if os.path.exists(pickle_ofn) and not force:
        dd =  pickle.load(open(pickle_ofn, 'rb')) 
        print(dd['qc'])
        return(dd)
    else:
    # make it if forced or not found
        bam = pysam.AlignmentFile(path_to_bam)
        cbc_corr = []
        umi_corr = []

        if verbose:
            print(f'Opening bam file {path_to_bam} to screen for correction pairs', end = '....')

        for read in bam:
            if read.has_tag('CR') and read.has_tag('CB'):
                cbc_corr.append((read.get_tag('CR'), read.get_tag('CB')))
            if read.has_tag('UR') and read.has_tag('UB'):
                umi_corr.append((read.get_tag('UR'), read.get_tag('UB')))

        if verbose:
            print(f'done')

        dict_cbc = dict(cbc_corr)
        dict_umi = dict(umi_corr)

        # quick QCs
        before_dict_cbc = set([i[1] for i in cbc_corr])
        before_dict_umi = set([i[1] for i in umi_corr])

        qc_str = f'Missed barcode cases = {len(before_dict_cbc) - len(set(dict_cbc.values()))}'
        qc_str = qc_str + f'Missed umi cases = {len(before_dict_umi) - len(set(dict_umi.values()))}'

        print(qc_str)

        dd = {'cell_barcodes':dict_cbc, 'umi_barcodes':dict_umi, 'qc':qc_str}
        with open(pickle_ofn, 'wb') as handle:
            pickle.dump(dd, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return(dd)


