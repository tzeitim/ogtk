import shutil
import os
import subprocess
import glob

def bulk_merge_dir(raw_path, out_path, do_merge, clean = True, force = False, errors =0):
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
            ret = genotype_merge_fastqs(f1 = f1, f2 = f2, out = out_merged, outu = out_unmerged1 , outu2 = out_unmerged2)
        if ret:
            unmerged_samples.append(basename)

def bbmerge_fastqs(f1, f2, out, outu, outu2, silent = True, bbmerge_bin = 'bbmerge.sh'):
    ''' use bbtools to merge amplicon derived reads 
    docs: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbmerge.sh
    ref: https://doi.org/10.1371/journal.pone.0185056
    If no present, install through conda bbtools
    '''
    cmd = "{} in1={} in2={} out={} outu={} outu2={} pfilter=0 mismatches=5 efilter=-1".format(bbmerge_bin, f1, f2, out, outu, outu2)
    pp = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if not silent:
        print(pp.stderr.decode())

    if "RuntimeException" in pp.stderr.decode():
        #means the command failed
        return(False)
    else:
        result = ''.join([i for i in pp.stderr.decode().split('\n') if "Joined" in i])
        print(result)
        return(True) 
