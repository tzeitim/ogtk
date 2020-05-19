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


## BWA
def bwa_build_ref(conf_fn, force = False, prepare_workspace = False):
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    ref_fa = conf['input_files']['ref_fa']
    flag_file = ref_fa + ".bwt"
    if not os.path.exists(flag_file) or force:
        os.makedirs(flag_file)
        ref_cmd = "bwa index %s "%(ref_fa)
        print(ref_cmd)
        subprocess.run(ref_cmd.split())
    else:
        print("Found pre-existing index: %s" % (flag_file))

def bwa_map(conf_fn, fq_reads, prefix,  force = False, prepare_workspace = False):
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader) 
    if prepare_workspace:
        conf_init_workspace(conf)
    bwa = conf['bwa']
    bwa_args = bwa.get('args', '')
    ref_fa = conf['input_files']['ref_fa']
    outbam = prefix + '.bam'

    bwa_build_ref(conf_fn, force)

    bwa_cmd = 'bwa bwasw -t %s %s %s %s' % (bwa['threads'], bwa_args, ref_fa, fq_reads)
    tobam_cmd = 'samtools view -b'
    index_cmd = 'samtools index %s' % (prefix + 'sorted.bam')

    if (not os.path.exists(outbam)) or force:
        print(bwa_cmd, "|", tobam_cmd , ">" , outbam)
        c1 = subprocess.Popen(bwa_cmd.split(), stdout=subprocess.PIPE)
        c2 = subprocess.Popen(tobam_cmd.split(), stdin=c1.stdout, stdout=open(outbam, "w"))
        c1.stdout.close()
        c2.communicate()

    else:
        print("Found pre-existing aligment: %s" %outbam)
 

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

 
