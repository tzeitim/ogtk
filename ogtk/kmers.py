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
import gzip


#https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
def build_kmers(sequence, ksize, as_set = True):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    if as_set:
        return set(kmers)
    else:
        return(kmers)

#https://gist.github.com/songweizhi/8f7759b602c81eb79fb8eb01950448ed
def de_bruijn_ize(st, k, nodes = None, edges = None):
    if edges == None:
        edges = []
    if nodes == None:
        nodes = set()
    for i in range(len(st) - k + 1):
        edges.append((st[i:i + k - 1], st[i + 1:i + k]))
        nodes.add(st[i:i + k - 1])
        nodes.add(st[i + 1:i + k])
    return ({"nodes":nodes, "edges":edges})


def visualize_de_bruijn(st,k):
    nodes, edges = de_bruijn_ize(st, k)
    dot_str = 'digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str += '    %s [label = "%s"] ;\n' % (node, node)
    for src, dst in edges:
        dot_str += '    %s -> %s ;\n' % (src, dst)
    return dot_str + '}\n'

def visualize_de_gp(nodes, edges):
    dot_str = 'digraph "DeBruijn graph" {\n'
    for i,node in enumerate(nodes):
        dot_str += '    %s [label = "%s"] ;\n' % (node, i)
    for src, dst in edges:
        dot_str += '    %s -> %s ;\n' % (src, dst)
    return dot_str + '}\n'
