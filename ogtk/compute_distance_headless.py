import os
import pandas as pd
import matplotlib.pyplot as plt
import sys, argparse
import pickle
import numpy as np
import ogtk


if __name__=='__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('mat', type=str, help='path to the barcode matrixy')
    parser.add_argument('out', type=str, help='path to the output score matrix')
    parser.add_argument('-c', '--cores', type=int, help='cores to perform computation', default=10)
    parser.add_argument('-w', '--weight', type=int, help='weight per shared barcode [0]. By default it computes a weight defined as 1/ncol(mat)', default=0)
    parser.add_argument('-x', '--exclude', type=int, help='exclude up to this number, e.g 1,2,3 are non informatinve, WT, NA and excision [3]', default=3)
    parser.add_argument('-v', '--verbose', type=int, help='print additional information', default=0)

    args=parser.parse_args()

    f_mat = np.load(args.mat)
    bcs = list(set(f_mat.flatten()))
    n_cols = f_mat.shape[1]
    print(f_mat.shape)

    if args.weight == 0:
        bias =[1/n_cols for i in range(n_cols)] 
    else:
        bias = args.weight
    cores = args.cores
    dont_include = args.exclude
    verbose = args.verbose
    dist_mat = ogtk.ltr.compute_dist_mat(mat=f_mat,  bcs = bcs, bias = bias, cores = cores, dont_include = dont_include)

    #out_fn = args.mat.replace("npy", "_dist.npy") 
    out_fn = args.out

    print(f'Saved to {out_fn}')
    np.save(out_fn, dist_mat)
