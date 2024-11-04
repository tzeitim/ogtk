"""
Genomic sequence mapping utilities providing wrappers for various alignment tools.

This module provides wrapper functions for common genomic sequence alignment tools
including BWA, STAR, minimap2, and EMBOSS suite. It handles index building,
mapping, and post-processing of alignment results.

Dependencies:
    - pysam
    - pyfaidx
    - yaml
    - numpy
    - pandas
    - matplotlib
    - BWA
    - STAR
    - minimap2
    - EMBOSS suite
    - samtools
"""

import os
import pysam
from pyfaidx import Fasta
import yaml
import subprocess
import regex
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional
from .aln_to_mat import Alg

def bwa_build_ref(conf_fn: str, force: bool = False, prepare_workspace: bool = False) -> None:
    """
    Build BWA index for reference genome.

    Args:
        conf_fn: Path to YAML configuration file
        force: If True, rebuild index even if it exists
        prepare_workspace: If True, initialize workspace directories

    Raises:
        FileNotFoundError: If reference FASTA file doesn't exist
        subprocess.CalledProcessError: If BWA index building fails
    """
    conf = yaml.safe_load(open(conf_fn))
    ref_fa = conf['input_files']['ref_fa']
    flag_file = ref_fa + ".bwt"
    
    if not os.path.exists(flag_file) or force:
        os.makedirs(os.path.dirname(flag_file), exist_ok=True)
        ref_cmd = f"bwa index {ref_fa}"
        print(f"Building BWA index: {ref_cmd}")
        subprocess.run(ref_cmd.split(), check=True)
    else:
        print(f"Found pre-existing index: {flag_file}")

def bwa_map(conf_fn: str, fq_reads: str, prefix: str, force: bool = False, 
            prepare_workspace: bool = False) -> None:
    """
    Map reads to reference using BWA.

    Args:
        conf_fn: Path to YAML configuration file
        fq_reads: Path to FASTQ input file
        prefix: Output file prefix
        force: If True, remap even if output exists
        prepare_workspace: If True, initialize workspace directories

    Raises:
        subprocess.CalledProcessError: If mapping fails
    """
    conf = yaml.safe_load(open(conf_fn))
    bwa = conf['bwa']
    ref_fa = conf['input_files']['ref_fa']
    outbam = f"{prefix}.bam"

    if prepare_workspace:
        conf_init_workspace(conf)  # Note: This function needs to be defined

    bwa_build_ref(conf_fn, force)

    if not os.path.exists(outbam) or force:
        bwa_cmd = f"bwa {bwa['alg']} -t {bwa['threads']} {bwa.get('args', '')} {ref_fa} {fq_reads}"
        print(f"Mapping reads: {bwa_cmd}")
        
        with subprocess.Popen(bwa_cmd.split(), stdout=subprocess.PIPE) as c1:
            with subprocess.Popen(['samtools', 'view', '-b'], 
                                stdin=c1.stdout, 
                                stdout=open(outbam, "wb")) as c2:
                c1.stdout.close()  # Allow c1 to receive a SIGPIPE if c2 exits
                c2.communicate()
    else:
        print(f"Found pre-existing alignment: {outbam}")

def star_build_ref(conf_fn: str, force: bool = False, prepare_workspace: bool = False) -> None:
    """
    Build STAR index for reference genome.

    Args:
        conf_fn: Path to YAML configuration file
        force: If True, rebuild index even if it exists
        prepare_workspace: If True, initialize workspace directories

    Raises:
        subprocess.CalledProcessError: If STAR index building fails
    """
    conf = yaml.safe_load(open(conf_fn))
    genomedir = conf['star']['genomedir']
    ref_fa = conf['input_files']['ref_fa']
    threads = conf['star']['threads']
    genomeSAindexNbases = conf['star']['genome_generation']['genomeSAindexNbases']
    env_prefix = conf['env_prefix']
    
    if prepare_workspace:
        conf_init_workspace(conf)  # Note: This function needs to be defined

    if not os.path.exists(genomedir) or force:
        os.makedirs(genomedir, exist_ok=True)
        ref_cmd = f"STAR --runMode genomeGenerate --genomeDir {genomedir} " \
                 f"--genomeFastaFiles {ref_fa} --runThreadN {threads} " \
                 f"--genomeSAindexNbases {genomeSAindexNbases}"
        print(f"Building STAR index: {ref_cmd}")
        subprocess.run(ref_cmd.split(), check=True)
    else:
        print(f"Found pre-existing genomedir: {genomedir}")

def star_map(conf_fn: str, fq_reads: str, prefix: str = '', force: bool = False, 
            prepare_workspace: bool = False, modality: Optional[str] = None) -> None:
    """
    Map reads to reference using STAR.

    Args:
        conf_fn: Path to YAML configuration file
        fq_reads: Path to FASTQ input file
        prefix: Output file prefix
        force: If True, remap even if output exists
        prepare_workspace: If True, initialize workspace directories
        modality: Optional sequencing modality specification

    Raises:
        subprocess.CalledProcessError: If mapping fails
    """
    conf = yaml.safe_load(open(conf_fn))
    star = conf['star']
    genomedir = star['genomedir']

    if prepare_workspace:
        conf_init_workspace(conf)

    aligned_bam = f"{prefix}Aligned.out.bam"
    sorted_bam = f"{prefix}sorted.bam"

    if not os.path.exists(aligned_bam) or force:
        # Build index if needed
        if not os.path.exists(genomedir):
            star_build_ref(conf_fn, force)

        # Run STAR alignment
        star_cmd = f"STAR --genomeDir {genomedir} --readFilesIn {fq_reads} " \
                  f"--runThreadN {star['threads']} --outFileNamePrefix {prefix} {star['options']}"
        subprocess.run(star_cmd.split(), check=True)

        # Sort and index BAM
        subprocess.run(f"samtools sort {aligned_bam} -o {sorted_bam}".split(), check=True)
        subprocess.run(f"samtools index {sorted_bam}".split(), check=True)
    else:
        print(f"Found pre-existing alignment: {aligned_bam}")

def mp2_map(ref: str, query: str, options: str, outfn: str, 
            verbose: bool = False, force: bool = False) -> None:
    """
    Map reads using minimap2 with optional BAM output processing.

    Args:
        ref: Path to reference sequence
        query: Path to query sequence
        options: Minimap2 command line options
        outfn: Output filename
        verbose: If True, print additional progress information
        force: If True, remap even if output exists

    Raises:
        subprocess.CalledProcessError: If mapping or BAM processing fails
    """
    if not os.path.exists(outfn) or force:
        import re
        rebam = re.compile(r'\.bam|\.sam|\.BAM|\.SAM')
        
        # Run minimap2
        mp2_cmd = f"minimap2 {options} {ref} {query} -o {outfn}"
        subprocess.run(mp2_cmd.split(), check=True)

        # Process BAM output if needed
        if rebam.search(outfn):
            import uuid
            temp_file = f"rm_me_{uuid.uuid4()}"
            
            # Convert to BAM, sort, and index
            for cmd in [
                f'samtools view -h -b {outfn} -o {temp_file}',
                f'samtools sort {temp_file} -o {outfn}',
                f'samtools index {outfn}'
            ]:
                if verbose:
                    print(cmd)
                subprocess.run(cmd.split(), check=True)
            
            # Cleanup
            subprocess.run(f'rm {temp_file}*'.split(), check=True)

def make_unique_fa_ref_entries(fa_ifn: str, fa_ofn: str, pattern: str = 'hspdrv7_scgstl',
                             omit_ref: bool = True, trim_start: Optional[int] = None, 
                             trim_end: Optional[int] = None) -> None:
    """
    Process FASTA file to make reference entries unique by appending UMIs.

    Args:
        fa_ifn: Input FASTA file
        fa_ofn: Output FASTA file
        pattern: Pattern to identify reference entries
        omit_ref: If True, exclude reference sequences from output
        trim_start: Start position for sequence trimming
        trim_end: End position for sequence trimming

    Raises:
        FileNotFoundError: If input file doesn't exist
    """
    entries = []
    whole_fa = []
    
    with open(fa_ifn) as ff:
        for i, line in enumerate(ff):
            if line.startswith(">"):
                entries.append(i)
            whole_fa.append(line)

    whole_fa = np.array(whole_fa, dtype=object)
    start = 1 if pattern in whole_fa[0] else 0
    offset = -1 if pattern in whole_fa[0] else 1

    # Process UMIs and references
    for umi_idx, ref_idx in zip(
        itertools.islice(entries, start, None, 2),
        itertools.islice(entries, start + offset, None, 2)
    ):
        umi = whole_fa[umi_idx][1:].strip()
        whole_fa[ref_idx] = f"{whole_fa[ref_idx].strip()}_{umi}\n"

    if omit_ref:
        temp_file = f"tmpfa_{uuid.uuid4()}"
        with open(temp_file, 'w') as ofa:
            ofa.writelines(whole_fa)

        ff = Fasta(temp_file)
        with open(fa_ofn, 'w') as ofa:
            for entry in ff.keys():
                if pattern not in entry:
                    seq = ff[entry][trim_start:trim_end] if (trim_start is not None and 
                          trim_end is not None) else ff[entry][:]
                    ofa.write(f">{entry}\n{seq}\n")
        
        os.remove(temp_file)
    else:
        with open(fa_ofn, 'w') as ofa:
            ofa.writelines(whole_fa)
