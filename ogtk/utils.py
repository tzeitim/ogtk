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


def tadpole_fastqs(f1, out, verbose = False, k=66, tadpole_bin = 'tadpole.sh', bm1=1, bm2=1, mincontig = "auto", mincountseed = 100, return_contigs=False):
    ''' use tadpole from bbtools to assemble a cloud of sequences controlled by a UMI
    '''
    #tadpole.sh in=tmp/f1_ACTTCGCCAGAGTTGG_GTGCGAGAGGGTA.fastq out=mini k=66 overwrite=True bm1=1 bm2=1 mincountseed=4
    cmd = f"{tadpole_bin} in={f1} out={out} k={k} overwrite=True bm1={bm1} bm2={bm2} t=5 -Xmx6g mincontig={mincontig} mincountseed={mincountseed} rcomp=f"
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


