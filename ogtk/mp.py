# mapping routines

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

    bwa_cmd = 'bwa %s -t %s %s %s %s' % (bwa['alg'], bwa['threads'], bwa_args, ref_fa, fq_reads)
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

def star_map(conf_fn, fq_reads, prefix = '',  force = False, prepare_workspace = False, modality = None):
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

def mp2_map(ref, query, options, outfn, verbose = False, force = False):
        '''
        Wrapper for minimap2  and samtools commands to sort, index and binarize output
        '''
        if not os.path.exists(outfn) or force:
            import re
            rebam = re.compile('\.bam|\.sam|\.BAM|\.SAM')
            mp2_cmd = f"minimap2 {options} {ref} {query} -o {outfn}"
            subprocess.run(mp2_cmd.split())

            if rebam.search(outfn):
                import uuid
                unique_filename = "rm_me_"+str(uuid.uuid4())
                bam_cmd = f'samtools view -h -b {outfn} -o {unique_filename}'
                sort_cmd = f'samtools sort {unique_filename} -o {outfn}'
                index_cmd = f'samtools index {outfn}'

                if verbose:
                    print(bam_cmd)
                subprocess.run(bam_cmd.split())

                if verbose:
                    print(sort_cmd)
                subprocess.run(sort_cmd.split())

                if verbose:
                    print(index_cmd)
                subprocess.run(index_cmd.split())

                subprocess.run(f'rm {unique_filename}*'.split())
            else:
                print('Unsupported feature. output filename not sam/bam')
                return(None)

 
def embos_align_reads_to_ref(name, fa_ifn, fa_ofn, ref_path, ref_name = 'hspdrv7_scgstl', mode='needleman', gapopen = 20, gapextend = 1, verbose = False, thorough = False, trim_start = None, trim_end = None, omit_ref = True, out_format = 'sam'):
    pwalg_out = '{}_out_pair_algn.fa'.format(name)

    if verbose:
        print(f"gapopen\t{gapopen}")
        print(f"gapextend\t{gapextend}")

    if mode == 'waterman':
        # The emboss header iw wrong
        #'water -gapextend 1 -gapopen 20 -datafile EDNAFULL -awidth3=100 -aformat3 fasta  -asequence ptol2-hspdrv7_scgstl.fa -bsequence msa_in.fa  -aaccshow3 yes -outfile ww.fasta -snucleotide2  -snucleotide1' 
        cmd_template = f'water -asequence {ref_path} -sformat1 fasta -bsequence {fa_ifn} -sformat2 fasta -gapopen {gapopen} -gapextend {gapextend}'

        cmd_waterman = f'{cmd_template} -outfile {pwalg_out} -aformat3 {out_format}'
        if verbose: print(cmd_waterman) 
        oo = subprocess.run(cmd_waterman.split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)
        
        if thorough:
            cmd_waterman = f'{cmd_template} -outfile {pwalg_out} -aformat3 pair'
            if verbose: print(cmd_waterman) 
            oo = subprocess.run(cmd_waterman.split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            if verbose: print(oo.stdout, oo.stderr)


        subprocess.run(f'mv {pwalg_out} {fa_ofn}'.split())


    if mode == 'needleman':
        cmd_template = 'needleall -gapextend {} -gapopen {} -datafile EDNAFULL -awidth3=100  -minscore 90 -bsequence {} -asequence {} -aaccshow3 yes'.format(gapextend, gapopen, ref_path, fa_ifn) 
    
        cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "fasta", pwalg_out)
        if verbose: print(cmd_needleman) 
        oo= subprocess.run(cmd_needleman.split(), stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        if verbose: print(oo.stdout, oo.stderr)
        
        make_unique_fa_ref_entries(fa_ifn = pwalg_out, fa_ofn = fa_ofn, pattern = ref_name, omit_ref= omit_ref, trim_start = trim_start, trim_end = trim_end)

        if thorough:
            cmd_needleman = '{} -aformat3 {} -outfile {}'.format(cmd_template, "score", pwalg_out.replace('fa', 'sam'))
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
    


def make_unique_fa_ref_entries(fa_ifn, fa_ofn, pattern = 'hspdrv7_scgstl', omit_ref = True, trim_start = None, trim_end = None):
    '''Small helper function that relabels reference entries on a pair-wise sequence alignments.
        It appends the UMI to the ref name; pattern is used to identify entries that belong to the reference since emboss includes both (ref, query) in its pairwise alignmnts'''
    FF = open(fa_ifn)
    entries = []
    whole_fa = []
    for i,line in enumerate(FF):
        if line[0] == ">":
             entries.append(i)
        whole_fa.append(line)

    FF.close()
    # if we don't set dtype to object, then the appended string is trimmed byt numpy
    whole_fa = np.array(whole_fa, dtype=object)
    start  =  1 if pattern in whole_fa[0] else 0
    offset = -1 if pattern in whole_fa[0] else 1
    
    umis = itertools.islice(entries, start, None, 2)
    refs = itertools.islice(entries, start+offset, None, 2)

    for umi,ref in zip(umis, refs):
        umi = whole_fa[umi][1:].strip()
        new_ref = "{}_{}\n".format(whole_fa[ref].strip(), umi)
        whole_fa[ref] = new_ref

    if omit_ref:
        import uuid
        unique_filename = "tmpfa_"+str(uuid.uuid4())
        with open(unique_filename, 'w') as ofa:
            for i in whole_fa:
                ofa.write(i)

        FF = Fasta(unique_filename)
        with open(fa_ofn, 'w') as ofa:
            for entry in FF.keys():
                if not pattern in entry:
                    if trim_start != None and trim_end != None:
                        ofa.write(f">{entry}\n{FF[entry][trim_start:trim_end]}\n")
                        #print(f">{entry}\n{FF[entry][trim_start:trim_end]}\n")
                    else:
                        ofa.write(f">{entry}\n{FF[entry][:]}\n")
        subprocess.run(f'rm {unique_filename}'.split())
    else: 
        with open(fa_ofn, 'w') as ofa:
            for i in whole_fa:
                ofa.write(i)
     




