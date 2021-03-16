import os
import pandas as pd
import matplotlib.pyplot as plt
import sys, argparse
import pickle
import ogtk


if __name__=='__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('in_pickle', type=str, help='path to the pickled readset for which to perform assembly')
    parser.add_argument('k', type=int, help='k parameter for building kmer graph', default=80)
    parser.add_argument('-v', '--verbose', type=int, help='k parameter for building kmer graph', default=0)
    #parser.add_argument("--verbose", type=str2bool, nargs='?',
    #                    const=True, default=False,
    #                    help="Activate nice mode.")

    args=parser.parse_args()

    verbose = args.verbose
    ogtk.UM.assemble_pickled(pickle_in = args.in_pickle, k = args.k, verbose = verbose == 1)

