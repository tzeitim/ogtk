import pyaml
from scipy.spatial.distance import pdist
import fastcluster
from scipy.cluster import hierarchy
import ogtk
import pdb
import pandas as pd 
import numpy as np
from colorhash import ColorHash
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

class matlin():
    def __init__(self):
        # code mm,i,del
        self.non_inf = {'NA': '   NA', 'wt': '   wt', 'ex':'   ex'}
        self.df = None
        self.bc = None
        self.name = ''

        self.lanes = []
        self.filtered = None

    def indel_qc(self):
        self.bcs  = pd.Series(np.hstack(self.df.values)).value_counts(normalize = True)
        self.indels = pd.DataFrame([i.split(".") for i in self.bcs.index], columns = ['mm', 'ins' ,'del'])
        self.indels.apply(lambda x: len(x.unique())).plot(kind='bar')

    def merge_lane(self, matlin2, inplace=True, encode_mat = True):
        
        """ Appends dataframe from a matlin object. This assumes that both have been
        previously merged across all their integrations and most importantly, whitelisted """

        if len(self.lanes) == 0:
            self.lanes.append(0)
            self.df.index = [f'{i}-{self.lanes[-1]}' for i in self.df.index]
        # TODO: currently assumes that the second matlin has not been merged previously; 
        # add support for multiple-indexing
        matlin2.df.index = [f'{i}-{self.lanes[-1]+1}' for i in matlin2.df.index]

        if inplace:
            print(self.df.shape)
            print(matlin2.df.shape)
            self.df = self.df.append(matlin2.df)
            print(self.df.shape)
            self.lanes.append(self.lanes[-1]+1)
            if encode_mat:
                self.__encode_mat()
        else:
            return(self.df.append(matlin2))

    def merge_int(self, tab, encode_mat = True):
        """Applies a outer merge of integrations """
        if self.__assert_none("df"):
            self.df = self.__read_table(tab)
            self.intid_len = self.df.shape[1]
        else:
            df = self.__read_table(tab)
            print(self.df.shape, df.shape)
            print(len(set(self.df.index).intersection(df.index)))
            #pdb.set_trace()
            indi = self.df.index.union(df.index, sort=False)
            self.df = self.df.merge(df, left_index=True, right_index=True, how='outer', sort=False)
            self.df = self.df.loc[indi]
            print(len(set(self.df.index).intersection(df.index)))
            self.df.columns = [f'bint{i}' for i in range(self.df.shape[1])]
            self.df.fillna(self.non_inf['NA'], inplace = True)

        if encode_mat:
            self.__encode_mat()

    def filter_wl(self, wl):
        self.df = self.df.loc[self.df.index.intersection(wl)]

    def unique_alleles(self):
        unique_alls = ~self.bc.duplicated()
        return(unique_alls.index[unique_alls])

        return(unique_alls)#~self.bc.duplicated())
        
    def merge_int_list(self, yml_path_list, encode_mat = True):
        for path in yml_path_list:
            conf = pyaml.yaml.load(open(path), Loader=pyaml.yaml.Loader) 
            self.merge_int(conf['lineage']['merged_tab'], encode_mat = encode_mat)
   
    def plot_mat(self, r, rows = range(0, 100),  cells = None, vmax=None):
            
        color_dict = dict(map(lambda x: (x, ColorHash(x).hex), set(self.df.values.flatten()) ))
        color_dict[self.non_inf['ex']] = "#000000"
        color_dict[self.non_inf['NA']] =  "#dddddd"
        color_dict[self.non_inf['wt']] =  "#ffffff"

        z = self.df.stack().sort_values().drop_duplicates().apply(lambda x: color_dict[x]).to_list()
        # factorize indels
        cmap = ListedColormap(z)
        if self.__assert_none('bc'):
            self.__encode_mat()
        mat = self.bc
        if cells.__class__.__name__ != 'NoneType':
            im = mat.loc[cells].iloc[:,r]
        else:    
            im = mat.iloc[rows, r]
        plt.close()
        if vmax == None:
            vmax = len(z)
        fig, ax = plt.subplots()
        ax.matshow(im, cmap=cmap,  aspect='auto', vmin=1, vmax=vmax)
        ax.set_xticks([i-0.5 for i in range(0, im.shape[1], 7)])
        ax.xaxis.grid(True, which='major')
        ax.set_xticklabels('')

        print(im)

    def plot_indel_freq(self, top_n, log = True, lowest_rank = 0):
        df = self.df.apply(lambda x: x.value_counts().head(top_n), axis = 0)
        fig, ax = plt.subplots()
        if log:
            ax.matshow(np.log(df.loc[df.sum(axis=1).sort_values(ascending = False).index.to_list()].iloc[lowest_rank:,:]), aspect = 'auto')
        else:
            ax.matshow(df.loc[df.sum(axis=1).sort_values(ascending = False).index.to_list()].iloc[lowest_rank:,:], aspect = 'auto')
        ax.set_xticks([i-0.5 for i in range(0, df.shape[1], 7)])
        ax.xaxis.grid(True, which='major')
        ax.set_xticklabels('')

    def allele_distance(self, cores = 0):
        if self.__assert_none('amat'):
            self.__encode_mat()

        dist = ogtk.ltr.compute_dist_mat(
                        self.amat, 
                        bcs=self.amat.stack().drop_duplicates(), 
                        bias = [1/self.amat.shape[1] for i in range(self.amat.shape[1])], cores=cores)
        self.dist = dist

    def hclust(self, condense = True, optimal_ordering = False):
        if condense:
            dist = pdist(self.dist)
        else:
            dist = self.dist
        print('computing linkage')
        if optimal_ordering:
            Z = hierarchy.linkage(dist, optimal_ordering=optimal_ordering)
        else:
            Z = fastcluster.linkage(dist)

        zl = hierarchy.leaves_list(Z)

        self.dist = pd.DataFrame(self.dist)
        self.dist.index = self.amat.index
        self.dist.columns = self.dist.index
        hsorted_index = self.dist.index[zl]
        self.dist  = self.dist.loc[hsorted_index, hsorted_index]
        allele_size = self.ar.value_counts()
        self.expanded_alleles = np.hstack(np.array([ [i]*allele_size[i] for i in hsorted_index], dtype=object))
        pdb.set_trace()
        # we need to get a cell indexes at the cell level 

        #self.bc = self.bc.loc[self.dist.index]
        self.Z = Z
        self.zl = zl

    def plot_detection(self):
        mask = self.df != self.non_inf['NA']
        mask.apply(lambda x: sum(x)/7, axis = 1).value_counts(sort=True).plot(kind='bar')

    def intid(self, x):
        x = x * 7
        return(range(x, x+7))

    def __read_table(self, ifn, omit_mm = True):
        tt  = pd.read_table(ifn, header=None)
        # TODO add headers since the beginning
        tt.columns = np.hstack(['cbc', [f'bint{i}' for i in range(tt.shape[1]-1)]])
        tt.set_index('cbc', inplace = True, drop = True, append=False)
        f_na = tt==".NA."
        f_ex = tt==".exx"
        f_wt = tt==".."
        if omit_mm:
            print('removing mms')
            tt = tt.apply(lambda x : [self.__omit_mm(i) for i in x] ) 
        tt[f_na] = self.non_inf["NA"]
        tt[f_wt] = self.non_inf["wt"]
        tt[f_ex] = self.non_inf["ex"]
        #pdb.set_trace()
        #tt = tt.apply(lambda x: [self.annotate_bint(i, x) for i in x], axis = 0)
        return(tt)

    def __omit_mm(self, bc):
        ''' Omits a given index from a barcode of the format mm.ins.del'''
        bcs = bc.split('.')
        return(f'.{bcs[1]}.{bcs[2]}')

    def __encode_mat(self):
        bcs = self.df.stack()#.rank(method='dense')#.unstack()
        codes = bcs.rank(method='dense')
        bcs_dic = dict(zip(codes.drop_duplicates(), bcs.drop_duplicates()))
        self.bc = codes.unstack()
        self.ar = self.bc.apply(self.__concat_allele_codes, axis=1)
        # allele matrix: encoded indels matrix without duplicates
        self.amat = self.bc[~self.bc.duplicated()]
        self.amat.set_index(self.amat.apply(self.__concat_allele_codes, axis=1), inplace=True)

    def __concat_allele_codes(self, x):
        return('.'.join([str(int(i)) for i in x]))
    def __assert_none(self, attr):
        return(self.__getattribute__(attr).__class__.__name__ == "NoneType")
    def __annotate_bint(self, bc, x):
        ''' Deprecated as it is best to keep track of the same NHEJ results '''
        if bc in self.non_inf.values():
            return(bc)
        else:
            return(f"{x.name}.{bc}")




