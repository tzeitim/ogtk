import shutil
import numpy as np
import time
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


def bulk_merge_dir(raw_path, out_path, clean = True, force = True, verbose = False, file_pattern = "fastq.gz"):
    ''' 
        This wrapper makes use of bbmerge to merge fastqs locate in a given directory. 
        file_pattern can be 'fastq' or 'fastq.gz'
        clean: removes any previous instance of `out_path`
    '''
    #TODO implement arguments options

    if clean and os.path.exists(out_path):
        shutil.rmtree(out_path, ignore_errors=True)
        os.makedirs(out_path)
        os.makedirs(out_path+"/umerged")

    if not os.path.exists(out_path):
        os.makedirs(out_path)
        os.makedirs(out_path+"/umerged")

    ###
    returns = []
    for i in glob.glob(raw_path+"/*R1*"+file_pattern):
        f1 = i
        basename = f1.split('/')[-1]
        basepath = f1.split('/')[0:-1]

        f2 = "/".join(basepath)+"/"+basename.replace('_R1_', '_R2_') 
        basename = basename.split('_R1')[0]

        out_merged = "{}/{}_merged.fastq.gz".format(out_path, basename)
        out_unmerged1 = "{}/umerged/{}_unmerged_R1.fastq.gz".format(out_path, basename)
        out_unmerged2 = "{}/umerged/{}_unmerged_R2.fastq.gz".format(out_path, basename)
        unmerged_samples = []
        
        ret = False
        if not os.path.exists(out_merged) or force:
            ret = bbmerge_fastqs(f1=f1, f2=f2, out=out_merged, outu=out_unmerged1, outu2=out_unmerged2, verbose=verbose)
            returns.append(ret)
        if ret:
            unmerged_samples.append(basename)
    return(returns)

def bbmerge_fastqs(f1, f2, out, outu, outu2, verbose = False, bbmerge_bin = 'bbmerge.sh'):
    ''' use bbtools to merge amplicon derived reads. Returns (total, joined, log_ofn) tuple per file found.
    docs: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbmerge.sh
    ref: https://doi.org/10.1371/journal.pone.0185056
    If no present, install through conda bbtools
    '''
    cmd = "{} in1={} in2={} out={} outu={} outu2={} pfilter=0 mismatches=5 efilter=-1".format(bbmerge_bin, f1, f2, out, outu, outu2)

    if verbose:
        print(cmd, end='\n')

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
        log_ofn = out.replace('.gz', '') 
        log_ofn = log_ofn.replace('fastq', 'log') 
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

def tadpole_correct_fastqs(f1, out, verbose = False, k=66, threads =1, tadpole_bin = 'tadpole.sh', bm1=1, bm2=1, mincontig = "auto", mincountseed = 100, return_seqs=True, return_cmd = False):
    ''' use tadpole from bbtools to error-correct a cloud of sequences controlled by a UMI
    '''
    #tadpole.sh in=tmp/f1_ACTTCGCCAGAGTTGG_GTGCGAGAGGGTA.fastq out=mini k=66 overwrite=True bm1=1 bm2=1 mincountseed=4
    cmd = f"{tadpole_bin} in={f1} out={out} k={k} mode=correct overwrite=True bm1={bm1} bm2={bm2} t={threads} -Xmx6g mincontig={mincontig} mincountseed={mincountseed} rcomp=f"
    if verbose:
        print(cmd)
    if return_cmd:
        return(cmd)
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
    ''' wrapper for shell's find  '''
    import subprocess

    sys_call = subprocess.run(f'find {path} {options} {pattern}'.split(), stdout =subprocess.PIPE, text = True)
    ls_call = [i for i in sys_call.stdout.split()] 
    return(ls_call)

# TODO write a generalized version for tabulate_*fastqs

def tabulate_single_umified_fastq(r1, cbc_len =16 , umi_len = 10, end = None, single_molecule = False, force = False, comparable=False, rc=False):
    ''' tabulate fastq into tabix, supports a single long read 
        |readid | start | end  | cbc | umi | seq | qual|

    5ʼkit cbc_len = 16 ; umi_len = 10
    3ʼkit cbc_len = 16 ; umi_len = 12

    if single_molecule = cbc and umi are the same
    if comparable it creates a tabix file for comparison purposes at a specified depth 1e5 by default
    previous name: tabulate_10x_fastqs
    TODO: implent full python interface otherwise bgzip might not be available in the system
    '''
    import gzip
    import bgzip

    unsorted_tab = r1.split('_R1_')[0]+'.unsorted.txt'
    sorted_tab = unsorted_tab.replace('unsorted.txt', 'sorted.txt.gz')

    # round end to closest order of magnitude in powers of ten
    if end is not None:
        end = int(10**int(np.log10(end)))

    if comparable and end is not None:
        f_end = f'{end:0.0e}'.replace('+','')
        unsorted_tab = unsorted_tab.replace('.txt', f'.{f_end}.comparable.txt')
        sorted_tab = sorted_tab.replace('.txt', f'.{f_end}.comparable.txt')

    if os.path.exists(sorted_tab) and not force:
        print(f'using pre-computed {sorted_tab}')
        return(sorted_tab)

    print(subprocess.getoutput('date'))

    with open(unsorted_tab, 'wt') as R3:
        with gzip.open(r1, 'rt') as R1:#, bgzip.BGZipWriter(R3) as zR3:
            i1 = itertools.islice(grouper(R1, 4), 0, end)
            for read1 in i1:
                read_id = read1[0].split(' ')[0]
                cbc_str = read1[1][0:cbc_len] if not single_molecule else read1[1][0:cbc_len+umi_len] 
                umi_str = read1[1][cbc_len:cbc_len+umi_len] if not single_molecule else read1[1][0:cbc_len+umi_len] 
                seq =   read1[1][cbc_len+umi_len:].strip() 
                qual =  read1[3][cbc_len+umi_len:].strip()
                if rc:
                    seq = rev_comp(seq)
                    qual = qual[::-1]
                out_str = '\t'.join([read_id, '0', '1', cbc_str, umi_str, seq, qual])+'\n'
                if len(seq)>cbc_len+umi_len: 
                    R3.write(out_str)
    cmd_sort = f'sort -k4,4 -k5,5 {unsorted_tab}'
    cmd_bgzip = 'bgzip'
    
    c1 = subprocess.Popen(cmd_sort.split(), stdout = subprocess.PIPE)
    c2 = subprocess.Popen(cmd_bgzip.split(), stdin = c1.stdout, stdout = open(sorted_tab, 'w'))
    c1.stdout.close()
    c2.communicate()

    pool = [c1, c2]
    while any(map(lambda x: x.poll() == None, pool)):
        time.sleep(1)

    map(lambda x: x.terminate(), pool)
 
    print(f"tabbed fastq {sorted_tab}")

    cmd_tabix = f"tabix -f -b 2 -e 3 -s 4 --zero-based {sorted_tab}"
    c3 = subprocess.run(cmd_tabix.split())
    del_unsorted = subprocess.run(f'rm {unsorted_tab}'.split())
    print(subprocess.getoutput('date'))
    return(sorted_tab)

def paired_umified_fastqs_to_parquet():
    ''' convert paired umified fastqs (e.g 10x) into parquet format
    '''

def tabulate_paired_umified_fastqs(r1, cbc_len =16 , umi_len = 10, end = None, single_molecule = False, force = False, comparable=False, rev_comp_r2=False, export_parquet = False):
    ''' merge umified paired fastqs (e.g. 10x fastqs) into a single one tab file with fields:
        |readid | start | end  | cbc | umi | seq | qual|

    5ʼkit cbc_len = 16 ; umi_len = 10
    3ʼkit cbc_len = 16 ; umi_len = 12

    if single_molecule = cbc and umi are the same
    if comparable it creates a tabix file for comparison purposes at a specified depth 1e5 by default
    ``rev_comp_r2`` reverse-complements R2

    if export_parquet is True an additional instance of the tabulated reads will get written as a parquet file.
    previous name: tabulate_10x_fastqs
    TODO: implent full python interface otherwise bgzip might not be available in the system
    '''
    import gzip
    import bgzip
    r2 = r1.replace("R1", "R2")

    rc = 'rc_' if rev_comp_r2 else ''
    unsorted_tab = r1.split('_R1_')[0]+f'.{rc}unsorted.txt'
    sorted_tab = unsorted_tab.replace('unsorted.txt', 'sorted.txt.gz')
    parfn = unsorted_tab.replace('unsorted.txt', '.parquet')

    # round end to closest order of magnitude in powers of ten
    if end is not None:
        end = int(10**int(np.log10(end)))

    if comparable and end is not None:
        f_end = f'{end:0.0e}'.replace('+','')
        unsorted_tab = unsorted_tab.replace('.txt', f'.{f_end}.comparable.txt')
        sorted_tab = sorted_tab.replace('.txt', f'.{f_end}.comparable.txt')

    if os.path.exists(sorted_tab) and not force:
        print(f'using pre-computed {sorted_tab}')
        if export_parquet:
            if os.path.exists(parfn) and not force:
                return(parfn)
            else:
                export_tabix_parquet(sorted_tab, parfn)
                return(parfn)
        return(sorted_tab)

    print(subprocess.getoutput('date'))

    with open(unsorted_tab, 'wt') as R3:
        #with gzip.open(r1, 'rt') as R1, gzip.open(r2, 'rt') as R2:#, bgzip.BGZipWriter(R3) as zR3:
        with gzip.open(r1, 'rt') as R1, gzip.open(r2, 'rt') as R2, bgzip.BGZipWriter(R3) as zR3:
            i1 = itertools.islice(grouper(R1, 4), 0, end)
            i2 = itertools.islice(grouper(R2, 4), 0, end)
            for read1,read2 in zip(i1, i2):
                #read2 = list(read2)
                read_id = read1[0].split(' ')[0]
                cbc_str = read1[1][0:cbc_len] if not single_molecule else read1[1][0:cbc_len+umi_len] # there's a bug here. for single-mol ? it was corrected to read1[1][0:cbc_len] 
                umi_str = read1[1][cbc_len:cbc_len+umi_len] if not single_molecule else read1[1][0:cbc_len+umi_len] # also here read1[1][0:cbc_len]
                seq = read2[1].strip() 
                qual =  read2[3].strip()
                if rev_comp_r2:
                    seq = rev_comp(seq)
                    qual = qual[::-1]
                out_str = '\t'.join([read_id, '0', '1', cbc_str, umi_str, seq, qual])+'\n'
                if len(seq)>cbc_len+umi_len: 
                    R3.write(out_str)
    cmd_sort = f'sort -k4,4 -k5,5 {unsorted_tab}'
    cmd_bgzip = 'bgzip'
    
    c1 = subprocess.Popen(cmd_sort.split(), stdout = subprocess.PIPE)
    c2 = subprocess.Popen(cmd_bgzip.split(), stdin = c1.stdout, stdout = open(sorted_tab, 'w'))
    c1.stdout.close()
    c2.communicate()

    pool = [c1, c2]
    while any(map(lambda x: x.poll() == None, pool)):
        time.sleep(1)

    map(lambda x: x.terminate(), pool)
 
    print(f"tabbed fastq {sorted_tab}")

    cmd_tabix = f"tabix -f -b 2 -e 3 -s 4 --zero-based {sorted_tab}"
    c3 = subprocess.run(cmd_tabix.split())
    del_unsorted = subprocess.run(f'rm {unsorted_tab}'.split())
    print(subprocess.getoutput('date'))
    
    if export_parquet:
        export_tabix_parquet(sorted_tab, parfn)

    return(sorted_tab)

def tabulate_umified_fastqs(r1, cbc_len =16, umi_len = 10, end = None, single_molecule = False, force = False, comparable = False):
    '''
    Generate tabix file from a single fastq. Use `tabulate_paired_umified_fastqs` instead for a paired-end mode
    if single_molecule: the usual 10x division of R1 into cbc and umi is ignored and a full-length umi is used insted

    if comparable it creates a tabix file for comparison purposes at a specified depth 1e5 by default
    '''
    import gzip
    import bgzip
    import os
    unsorted_tab = r1.split('_R1_')[0]+'.unsorted.txt'
    sorted_tab = unsorted_tab.replace('unsorted.txt', 'sorted.txt.gz')

    # round end to closest order of magnitude in powers of ten
    if end is not None:
        end = int(10**int(np.log10(end)))

    if comparable and end is not None:
        f_end = f'{end:0.0e}'.replace('+','')
        unsorted_tab = unsorted_tab.replace('.txt', f'.{f_end}.comparable.txt')
        sorted_tab = sorted_tab.replace('.txt', f'.{f_end}.comparable.txt')

    if os.path.exists(sorted_tab) and not force:
        print(f'using pre-computed {sorted_tab}')
        return(sorted_tab)

    print(subprocess.getoutput('date'))

    with open(unsorted_tab, 'wt') as R3:
        with gzip.open(r1, 'rt') as R1:#, bgzip.BGZipWriter(R3) as zR3:
            i1 = itertools.islice(grouper(R1, 4), 0, end)
            for read1 in i1:
                #read2 = list(read2)
                read_id = read1[0].split(' ')[0]
                cbc_str = read1[1][0:cbc_len] if not single_molecule else read1[1][0:cbc_len+umi_len] 
                umi_str = read1[1][cbc_len:cbc_len+umi_len] if not single_molecule else read1[1][0:cbc_len+umi_len] 
                seq = read1[1][cbc_len+umi_len:].strip() 
                qual = read1[3][cbc_len+umi_len:].strip() 
                out_str = '\t'.join([read_id, '0', '1', cbc_str, umi_str, seq, qual])+'\n'
                if len(seq)>cbc_len+umi_len: 
                    R3.write(out_str)
    cmd_sort = f'sort -k4,4 -k5,5 {unsorted_tab}'
    cmd_bgzip = 'bgzip'
    
    c1 = subprocess.Popen(cmd_sort.split(), stdout = subprocess.PIPE)
    c2 = subprocess.Popen(cmd_bgzip.split(), stdin = c1.stdout, stdout = open(sorted_tab, 'w'))
    c1.stdout.close()
    c2.communicate()

    pool = [c1, c2]
    while any(map(lambda x: x.poll() == None, pool)):
        time.sleep(1)

    map(lambda x: x.terminate(), pool)
 
    print(f"tabbed fastq {sorted_tab}")

    cmd_tabix = f"tabix -f -b 2 -e 3 -s 4 --zero-based {sorted_tab}"
    c3 = subprocess.run(cmd_tabix.split())
    del_unsorted = subprocess.run(f'rm {unsorted_tab}'.split())

    print(subprocess.getoutput('date'))
    return(sorted_tab)


def df_to_tabix():
    '''Use system calls to tabindex a dataframe
    '''


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
    ''' Transfers the CR and UR tags on bam files generated after mapping merged fastqs (using `merge_11_fastqs`)'''
    print(f"annotating {ibam}")
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


def merge_10x_fastq_dir(indir, force = False, read1_pattern = "_R1_"):
    import os
    ''' merges all fastq pairs in a given directory and annotates (transfers CR and UR tags) by defautl (controlled by annotate)'''
    for r1 in os.listdir(indir):
        if read1_pattern in r1 and "fastq.gz":
            print("/".join([indir,r1]))
            # define a mappe
            #mbam_path = indir + r1.replace("R1", "R3").replace("fastq.gz", "bam")
            
            r3 = (indir+r1).replace(read1_pattern, "_R3_")
            if not os.path.exists(r3) or force:
                print(f"merging r1 and r2 as merged read didn't exist")
                r3_name = merge_10x_fastqs(indir, r1) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                print(r3_name, r3)
            else:
                print(f"previous merged files detected, set force = True to overwrite or delete {r3}")
            #if not os.path.exists(mbam_path):
            #    print(f"mapping to {mbam_path}")
            #    ogtk.mp.mp2_map(ref="/local/users/polivar/src/ivltr/refs/transgenes/trans.fa ", 
            #                query = indir + r3, 
            #                outfn = mbam_path, verbose = True, options = "-a --cs --splice")
            
            #ogtk.UM.annotate_bam(mbam_path)
            
            #print(mbam_path)
     

def generate_cbc_correction_dictionaries(path_to_bam, force = False, verbose = False, chunk_size=1e6):
    '''
    Generates a cellbarcode correction dictionary based on the cellranger barcoded BAM tags
    TODO: there is a small fraction of uncorrected cell barcodes that seem to map to more than one corrected cell barcode
    '''
    import pysam
    import os
    import polars as pl

    ofn = path_to_bam.replace('.bam', '.parquet')
    print(f'scanning {path_to_bam}')
    # load a pre-computed dictionary of found
    if os.path.exists(ofn) and not force:
        print('pre-computed map found, returning file name. Use force to regenerate')
        return(ofn)
    else:
    # make it if forced or not found
        entries = []
        if verbose:
            print(f'Opening bam file {path_to_bam} to screen for correction pairs', end = '....')
        print('') 
        bam = pysam.AlignmentFile(path_to_bam)

        df = pl.DataFrame(columns=[('CR', pl.Utf8), ('CB', pl.Utf8)] )
        
        for read in bam:
            if read.has_tag('CR') and read.has_tag('CB'):
                entries.append((read.get_tag("CR"), read.get_tag("CB")))
            

            if len(entries)%chunk_size ==0:
                df = (df.vstack(pl.DataFrame(entries, orient='row', columns=[('CR', pl.Utf8), ('CB', pl.Utf8)] ))
                     .unique()
                    )
                entries = []

        df = df.unique().write_parquet(ofn)

        if verbose:
            print(f'done')

        return(ofn)

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


def fill_sge_template(job_name, cmd, wd, pe = 'smp', ram = '7G', run_time = "00:10:00", job_reservation = False):
    ''' Main utility to generated a sge job script
    '''
    job_out = f'{wd}/{job_name}.stdout'
    job_err = f'{wd}/{job_name}.stderr'
    args = locals()
    template = return_sge_template()
    
    for k in args.keys():
        if k == 'job_reservation':
            if args[k]:
                args[k] = 'y'
            else:
                args[k] = 'n'
        template = template.replace('{'+k.upper()+'}', str(args[k]))

    return(template)

def return_sge_template(**kargs):

    ''' Helper function, returns an empty job template
    options = ['interconnects']
    '''

    options = {'interconnects': "#$ -l {INTER}            # select ib for infiniband or opa for omni-path\n"}
    template = "\
    #!/bin/bash\n\
    \n\
    #$ -N {JOB_NAME}       # specify a job name (script name as default)\n\
    #$ -o {JOB_OUT}        #\n\
    #$ -e {JOB_ERR}        #\n\
    #$ -pe {PE}      # use parallel environments orte and 32 slots\n\
    #$ -l m_mem_free={RAM} # max memory (RAM) per slot (4G as default on all.q)\n\
    #$ -l h_rt={RUN_TIME} # max run time in hh:mm:ss (4 days as default on all.q)\n\
    #$ -R {JOB_RESERVATION}             # use job reservation, needed for MPI jobs\n\
    #$ -cwd             # start the job in current working directory\n\
    "

    guix_template='\n\
    # Set user~@~Ys default guix profile:\n\
    export GUIX_PROFILE="$HOME/.guix-profile"\n\
    # Load the guix profile:\n\
    source "$GUIX_PROFILE/etc/profile"\n\
    # Use the guix locale, not of the host system\n\
    export GUIX_LOCPATH="$GUIX_PROFILE/lib/locale"\n'

    guix_original='\
    source $GUIX_PROFILE/etc/profile\n\
    /home/polivar/.guix-profile\n\
    ### add guix environment variables\n\
    source ~/.guix-profile/etc/profile\n\
    ### stacksize unlimited, needed for large MPI jobs\n\
    ulimit -s unlimited\n'

    cmd_template = "\
    ulimit -s unlimited\n\
    ### Start the job\n\
    {CMD}\n\
    "
    if len(kargs.keys()) >0:
        for key in kargs.keys():
            if key in options.keys():
                template = template + options[key]

    template = '\n'.join([template , guix_template])
    template = '\n'.join([template , cmd_template])

    return(template)

def split_bam_by_tag(bam_ifn, tab, prefix, tag = 'CB', header = None):
    import pysam
    import pandas as pd 
    import subprocess
    print(subprocess.getoutput('date'))

    df = pd.read_csv(tab, sep='\t', header=header)
    df.columns = [tag, 'ann']
    bam = pysam.AlignmentFile(bam_ifn, 'rb')
    # map tag->ann
    tag_ann_dict = dict(zip(df[tag], df['ann']))
    anns = set(tag_ann_dict.values())
    print(f" loaded {len(tag_ann_dict)} keys for {len(anns)} annotations")

    # create dictionaty of outstreams
    outs = {}
    for ann in anns:
        outs[ann] = pysam.AlignmentFile(f'{prefix}_{ann}.bam','wb', template = bam)
    outs['undet'] = pysam.AlignmentFile(f'{prefix}_undet.bam','wb', template = bam)
 
    for read in bam.fetch():
        if read.has_tag(tag):
            foc_tag = read.get_tag(tag)
            
            if foc_tag in tag_ann_dict.keys():
                outs[tag_ann_dict[foc_tag]].write(read)
            else:
                outs['undet'].write(read)

    for ann in outs.keys():
        outs[ann].close()

    print(subprocess.getoutput('date'))


def list_to_fasta(fasta_ofn, seqs, names=None):
    '''saves a list into a fasta file format, naming the sequence according to their index if names are not provided
    '''
    if names is None:
        names = range(len(seqs))
    with open(fasta_ofn, 'w') as fa:
        for name, seq in zip(names, seqs):
            fa.write(seq_to_fasta(seq, name)) 


def seq_to_fasta(seq, name=None):
    return(">%s\n%s\n"%(name,seq))

def rev_comp(seq):
    " rev-comp a dna sequence with UIPAC characters ; courtesy of Max's crispor"
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M',
    "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n", "V" : "B", "v":"b", 
    "B" : "V", "b": "v", "W" : "W", "w" : "w", "-":"-"}

    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

def export_tabix_parquet(tbxifn, parfn)->None:
    ''' Small helper function that opens a tabix (gzipped) file and exports it back as parquet
    '''
    import polars as pl
    import subprocess

    print(subprocess.getoutput('date'))

    df = pl.read_csv(tbxifn, sep='\t', has_header=False)
    df.columns=['readid',  'start' ,'end'  , 'cbc' , 'umi' , 'seq' , 'qual']
    df.write_parquet(parfn)

    print(f'saved {parfn}')
    print(subprocess.getoutput('date'))
