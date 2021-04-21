import pyaml
import pandas as pd 
import numpy as np
from colorhash import ColorHash
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

class matlin():
    def __init__(self):
        self.df = None
        self.name = ''

        self.non_inf = ["   NA", "   wt","   ex"]

    def __read_table(self, ifn):
        tt  = pd.read_table(ifn, header=None)
        # TODO add headers since the beginning
        tt.columns = np.hstack(['cbc', [f'bint{i}' for i in range(tt.shape[1]-1)]])
        tt.set_index('cbc', inplace = True, drop = True, append=True)
        tt[tt == ".NA."] = "   NA"
        tt[tt == ".."] = "   wt"
        import pdb
        tt[tt == "..exx"] = "   ex"
        #tt = tt.apply(lambda x: [self.annotate_bint(i, x) for i in x], axis = 0)
        return(tt)

    def merge_tab(self, tab):
        if self.df.__class__.__name__ != "DataFrame":
            self.df = self.__read_table(tab)
        else:
            df = self.__read_table(tab)
            print(self.df.shape, df.shape)
            self.df = self.df.merge(df, on='cbc', how='outer')
            self.df.columns = [f'bint{i}' for i in range(self.df.shape[1])]
            self.df.fillna('   NA', inplace = True)
        
    def merge_list(self, yml_path_list):
        for path in yml_path_list:
            conf = pyaml.yaml.load(open(path), Loader=pyaml.yaml.Loader) 
            self.merge_tab(conf['lineage']['merged_tab'])

    def plot_mat(self, r, cells = range(0,100), vmax=None):
        color_dict = dict(map(lambda x: (x, ColorHash(x).hex), set(self.df.values.flatten()) ))
        color_dict['   ex'] = "#ff0000"
        color_dict['   NA'] =  "#000000"
        color_dict['   wt'] =  "#ffffff"
        z = self.df.stack().sort_values().drop_duplicates().apply(lambda x: color_dict[x]).to_list()
        # factorize indels
        mat = self.df.stack().rank(method='dense').unstack()
        cmap = ListedColormap(z)
        im = mat.iloc[cells,r]
        plt.close()
        if vmax == None:
            vmax = len(z)
        plt.matshow(im, cmap=cmap,  aspect='auto', vmin=1, vmax=vmax)
        print(im)

    def annotate_bint(self, bc, x):
        if bc in self.non_inf:#["_NA", "_wt", "_ex"]:
            return(bc)
        else:
            return(f"{x.name}.{bc}")

    def intid(self, x):
        x = x * 7
        return(range(x, x+7))
