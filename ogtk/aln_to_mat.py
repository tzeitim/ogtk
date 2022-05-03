from pyfaidx import Fasta
import pdb
import numpy as np
import sys, argparse

class Pos:
    def __init__(self, index):
        self.index = index
        self.freq = {"a":0, "c":0, "t":0, "g":0, '-':0, 'n':0}

class Alg:
    def __init__(self, fastafn, freqfn, colorfn):
        self.pos = []
        self.init = False
        self.size = 0
        self.fasta = Fasta(fastafn)
        self.colorfn = colorfn
        self.conta = {'n':0, '-':0, 'a':1, 'c':2, 'g':3, 't':4, '\n':'\n'}

        self.read_fasta(fastafn)
        self.write_freqs(freqfn)

    def do_plot(self, plot, names = False):
        msa = self.seqtocol(self.colorfn, names= names)
        if plot:
            return(msa)
        

    def read_fasta(self, fastafn):
        for entry in self.fasta.keys():
            seq = self.fasta[entry][:]
            if not self.init:
                # this assumes that all the entries in the fasta record are the same size. 
                # this is the default setting for clustalo 
                # TODO add an assertion ro verify so
                self.size = len(seq) 
                for i in range(0, self.size):
                    self.pos.append(Pos(i))
                self.init = True

            for nt in range(0, self.size):
                self.pos[nt].freq[seq[nt].lower()]+=1 
    
    def seqtocol(self, outfn, names=False):
        outf = open(outfn, 'w')
        colors = []
        for i,entry in enumerate(self.fasta.keys()):
            outf.write(entry+','+','.join([str(self.conta[i.lower()]) for i in self.fasta[entry][:]])+'\n')
            if names:
                colors.append(entry)
            [colors.append(self.conta[i.lower()]) for i in self.fasta[entry][:]]
        outf.close()

        # TODO thisis very weird, check why one option returns the transpose
        if names:
            #colors = np.array(colors).reshape( 1+i, 1+len(self.fasta[entry][:])) 
            colors = np.array(colors).reshape( 1+len(self.fasta[entry][:]), 1+i) 
        else:
            colors = np.array(colors).reshape(1+i, len(self.fasta[entry][:])) 
        return(colors)

    def write_freqs(self, outfn):
        outf = open(outfn, 'w')
        outf.write('\t'.join(['a','c','t','g'])+'\n')
        for j in self.pos:
            outf.write('\t'.join([str(j.freq['a']),str(j.freq['c']),str(j.freq['t']),str(j.freq['g'])])+'\n')
        outf.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Covert ClustalOmega alignments for plotting on R")

    parser.add_argument('-a', '--algn', type=str, help='alignement file from clustalo', required=True)
    parser.add_argument('-o', '--outputfreqs', type=str, help='output filename for the nucleotide frequencies for each positions (occupancy)', required=True)
    parser.add_argument('-c', '--outputcolor', type=str, help='output filename for color coloded nucleotides', required=True)

    args=parser.parse_args()

    alignment = Alg(fastafn = args.algn, freqfn = args.outputfreqs, colorfn = args.outputcolor)

