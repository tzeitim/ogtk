import pyaml
from scipy.spatial.distance import pdist
import fastcluster
from scipy.cluster import hierarchy
import ogtk
import pdb
import pandas as pd 
import polars as pl
import numpy as np
from colorhash import ColorHash
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

class matlin():
    def __init__(self, bints = 7):
        # code mm,i,del
        self.non_inf = {'NA': '   NA', 'wt': '   wt', 'ex':'   ex'}
        self.df = None
        self.bc = None
        self.name = ''
        self.bins_per_intid = bints

        self.dist = None
        
        self.lanes = []
        self.filtered = None

        self.Z = None
        self.is_clustered = False


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
            self.df = self.df.append(matlin2.df)
            self.lanes.append(self.lanes[-1]+1)
            if encode_mat:
                self.__encode_mat()
        else:
            return(self.df.append(matlin2))

    def merge_int(self, tab, encode_mat = True, omit_mm = True):
        """Applies a outer merge of integrations """
        #if self.__assert_none("df"):
        if self.df is None:
            self.df = self.__read_table(tab, omit_mm = omit_mm)
            self.intid_len = self.df.shape[1]
        else:
            df = self.__read_table(tab, omit_mm = omit_mm)
            indi = self.df.index.union(df.index, sort=False)
            self.df = self.df.merge(df, left_index=True, right_index=True, how='outer', sort=False)
            self.df = self.df.loc[indi]
            self.df.columns = [f'bint{i}' for i in range(self.df.shape[1])]
            self.df = self.df.fillna(self.non_inf['NA'], inplace = False)

        if encode_mat:
            self.__encode_mat()

    def filter_wl(self, wl, copy = False):
        if copy:
            import copy
            new_instance = copy.deepcopy(self)
            new_instance.filter_wl(wl, copy = False)
            new_instance._matlin__encode_mat()
            return(new_instance)
            
        self.df = self.df.loc[self.df.index.intersection(wl)]
        self.__encode_mat()


    def unique_alleles(self):
        unique_alls = ~self.bc.duplicated()
        return(unique_alls.index[unique_alls])

        return(unique_alls)#~self.bc.duplicated())
        
    def merge_int_list(self, yml_path_list, encode_mat = True, omit_mm = True):
        if omit_mm:
            print('removing mms')

        for path in yml_path_list:
            conf = pyaml.yaml.load(open(path), Loader=pyaml.yaml.Loader) 
            self.merge_int(conf['lineage']['merged_tab'], encode_mat = encode_mat, omit_mm = omit_mm)

    def allele_distance(self, cores = 0):
        #if self.__assert_none('amat'):
        if self.amat is None:
            self.__encode_mat()

        bcs = self.amat.stack().drop_duplicates()
        dist = ogtk.ltr.compute_dist_mat(
                        self.amat, 
                        bcs=bcs,
                        bias = [1/self.amat.shape[1] for i in range(self.amat.shape[1])], cores=cores)
        
        
        dist = pd.DataFrame(dist)
        dist.index = self.amat.index
        dist.columns = dist.index
        # expand the encoded alleles by their number of children
        ii = np.hstack([ [i] * len(self.allele_children[i]) for i in self.amat.index])
        # get the original children ids
        iii = np.hstack([self.allele_children[i] for i in dist.index])
        # expand distance matrix
        dist = dist.loc[ii,:].loc[:,ii]
        # rename distance index on the cell space
        dist.index = iii
        dist.columns = iii
        
        self.dist = dist

        print(f'dist matrix computed {self.dist.shape}')

    def hclust(self, condense = True, method = 'ward', optimal_ordering = False):
        '''
        perform hierarchical clustering. By default a `fastcluster` command is invoked. if optimal_ordering is True then a slower modality using `scipy.cluster.hierarchy.linkage` is used. 
        '''
        if condense:
            dist = pdist(self.dist)
        else:
            dist = self.dist

        print('computing linkage')

        if optimal_ordering:
            print("optimal ordering can be slow")
            Z = hierarchy.linkage(dist,method=method, optimal_ordering=optimal_ordering)
        else:
            print("using fastcluster")
            Z = fastcluster.linkage(dist, method=method)

        zl = hierarchy.leaves_list(Z)

        #self.dist = pd.DataFrame(self.dist)
        #self.dist.index = self.amat.index
        hsorted_index = self.dist.index[zl]
        self.dist  = self.dist.loc[hsorted_index, hsorted_index]
        self.df = self.df.loc[hsorted_index, :]
        allele_size = self.ar['allele'].value_counts()

        #self.expanded_alleles = np.hstack(np.array([ [i]*allele_size[i] for i in hsorted_index], dtype=object))
        # we need to get a cell indexes at the cell level 
        #self.cells_hsorted = np.hstack(np.array([self.allele_children[i] for i in hsorted_index], dtype=object))

        #self.bc = self.bc.loc[self.cells_hsorted]
        self.bc = self.bc.loc[hsorted_index]
        self.Z = Z
        self.zl = zl
        self.is_clustered = True
   
    def cutting_stats(self):
        self.bc.apply(lambda x: [int(100*i)/100 for i in [np.mean(x==1), np.mean(x==2), np.mean(x==3), sum(x==3)/(len(x)-sum(x==1)-sum(x==2))]], axis = 0, result_type='reduce')

    def create_color_dict(self):
        color_dict = dict(map(lambda x: (x, ColorHash(x).hex), set(self.df.values.flatten()) ))
        color_dict[self.non_inf['ex']] = "#000000"
        color_dict[self.non_inf['NA']] =  "#dddddd"
        color_dict[self.non_inf['wt']] =  "#ffffff"
        return(color_dict)
        

    def plot_mat(self, col_range = None, rows = range(0, 100),  cells = None, vmax=None, add = False, ax = None):
            
        #color_dict = dict(map(lambda x: (x, ColorHash(x).hex), set(self.df.values.flatten()) ))
        #color_dict[self.non_inf['ex']] = "#000000"
        #color_dict[self.non_inf['NA']] =  "#dddddd"
        #color_dict[self.non_inf['wt']] =  "#ffffff"

        color_dict = self.create_color_dict()

        z = self.df.stack().sort_values().drop_duplicates().apply(lambda x: color_dict[x]).to_list()
        # factorize indels
        cmap = ListedColormap(z)

        #if self.__assert_none('bc'):
        if self.bc is None:
            self.__encode_mat()
        mat = self.bc
        #if cells.__class__.__name__ != 'NoneType':
        if cells is not None:
            #cells =  [i for i in mat.index if i in cells]
            cells = self.__intersect_cells(cells, mat.index)
            cells = [cells[i] for i in np.arange(len(cells)-1, -1, -1)]
            #return(cells)
            if col_range is not None:
                im = mat.loc[cells].iloc[:,col_range]
            else:
                im = mat.loc[cells]#.iloc[:,0:r]
        else:    
            if col_range is not None:
                im = mat.iloc[rows, col_range]
            else:
                im = mat.iloc[rows, :]

        if not add:
            plt.close()
        if vmax == None:
            vmax = len(z)
        if ax == None:
            fig, ax = plt.subplots()
        ax.matshow(im, cmap=cmap,  aspect='auto', vmin=1, vmax=None)#np.max(self.bc.to_numpy()))
        ax.set_xticks([i-0.5 for i in range(0, im.shape[1], self.bins_per_intid)])
        ax.xaxis.grid(True, which='major')
        ax.set_xticklabels('')

    def plot_indel_freq_mat(self, top_n, log = True, lowest_rank = 0):
        df = self.df.apply(lambda x: x.value_counts().head(top_n), axis = 0)
        fig, ax = plt.subplots()
        if log:
            ax.matshow(np.log(df.loc[df.sum(axis=1).sort_values(ascending = False).index.to_list()].iloc[lowest_rank:,:]), aspect = 'auto')
        else:
            ax.matshow(df.loc[df.sum(axis=1).sort_values(ascending = False).index.to_list()].iloc[lowest_rank:,:], aspect = 'auto')
        ax.set_xticks([i-0.5 for i in range(0, df.shape[1], self.bins_per_intid)])
        ax.xaxis.grid(True, which='major')
        ax.set_xticklabels('')

    def plot_indel_freq(self, cells = None, intid = None, axes = None):
        #if cells.__class__.__name__ != 'NoneType':
        if cells is not None:
            cells = self.__intersect_cells(cells, self.df.index)
        else:
            cells = self.df.index
        
        #if intid.__class__.__name__ == 'NoneType':
        if intid is None:
            intid = range(self.df.shape[1])
        else:
            intid = self.intid(intid)
    
        color = self.create_color_dict()
        aa = self.df.loc[cells].iloc[:,intid].apply(lambda x: x.value_counts())

        if axes is None:
            fig, axes = plt.subplots(2,1)

        for i in axes:
            i.set_facecolor("#fafafafa")
            i.set_axisbelow(True)

        #return(aa)
        aa.transpose().plot(kind='bar', stacked = True, legend = False, color = color, ax = axes[0])
        axes[0].grid(axis='y')
        
        aa[[not(i == self.non_inf['wt'] or i == self.non_inf['NA']) for i in  aa.index]].transpose().plot(kind='bar', legend = False, stacked = True, color = color, ax = axes[1])
        axes[1].grid(axis='y')

    def dextract(self, cells = None, expanded = True):
        ''' extracts a distance matrix for a given query '''
        #if self.__assert_none('dist'):
        if self.dist is None:
           raise ValueError("No distance matrix found. Generate one by calling .alle_distance and optionally .hclust to group alleles/cells") 
        #if cells.__class__.__name__ != 'NoneType':
        if cells is not None:
            # keep valid cells and hsort them
            ocells = self.__intersect_cells(cells, self.bc.index)

            cells = [cells[i] for i in np.arange(len(cells)-1, 0, -1)]
            # convert cells to allele indexes
            #cells = [self.ar['allele'][i] for i in ocells]
            mat = self.dist.loc[cells, cells].copy()
            mat.index = cells
            mat.columns = cells
            return(mat)
        else:
            if expanded:
                cells = self.expanded_alleles
            else:
                cells = self.dist.index
                return(self.dist.loc[cells, cells].copy())

    def plot_dist_heatmap(self, cmap = 'magma_r', cells = None, expanded = True, ax = None, aspect = 'equal'): 
        mat = self.dextract(cells = cells, expanded = expanded)
        if ax == None:
            fig, ax = plt.subplots()

        ax.matshow(mat, cmap=cmap, aspect = aspect)

    def plot_detection(self, ax = None):
        if ax == None:
            fig, ax = plt.subplots()

        mask = self.df != self.non_inf['NA']
        mask.apply(lambda x: sum(x)/self.bins_per_intid, axis = 1).value_counts(sort=True).plot(kind='bar', ax = ax)

    
    def plot_uncut_pct(self, cells = None , ax = None):
        def return_uncut_pct(x):
            y = x[x!=self.non_inf['NA']]
            y = y.value_counts(normalize = True).head(10)
            if self.non_inf['wt'] in y.keys():
                return(y[self.non_inf['wt']])
            else:
                print(y)
                return(np.nan)
        
        #if cells.__class__.__name__ != 'NoneType':
        if cells is not None:
            cells = self.__intersect_cells(cells, self.df.index)
        else:
            cells = self.df.index

        vals = []
        for i in range(0, int(self.df.shape[1]/self.bins_per_intid)):
            wt_per_int =self.df.loc[cells].iloc[:,self.intid(i)].apply(return_uncut_pct, axis = 0, result_type='reduce') 
            vals.append(wt_per_int.to_list())
        pd.DataFrame(np.array(vals)).transpose().plot(kind='bar', ax = ax)
        if ax != None:
            ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5)) 
        else:
            plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5)) 

        plt.title(f'{self.name} uncut fraction per intid', color='black')


    def intid(self, x):
        bints = self.bins_per_intid
        x = x * bints
        return(range(x, x+bints))

    def ingest_ibars_pl(self, 
            df, 
            encode_mat=True, 
            sort_by=['kalhor_id', 'nspeed', 'raw_ibar'],
            ):

        '''Direct convertion of a polars data frame into the character matrix
        of self.matlin object 
        '''
            
        df = df.sort(sort_by)
                
        df =(df
            .with_columns(pl.when(~pl.col('wt'))
                .then(pl.col('seq'))
                .otherwise("..").alias('seq'))
            .pivot(columns='raw_ibar', index='cbc', values='seq')
            )

        self.df  = df.to_pandas().set_index('cbc')
        self.df = self.df.fillna(self.non_inf['NA'], inplace = False)
        self.intid_len = self.df.shape[1]
        f_wt = self.df==".."
        self.df[f_wt] = self.non_inf["wt"]

        self.ibar_to_bint = dict([(f'bint{i}',v) for i,v in enumerate(self.df.columns)])
        self.df.columns = [f'bint{i}' for i in range(self.df.shape[1])]

        if encode_mat:
            self.__encode_mat()

    def ingest_ibars_csv(self, ifn, encode_mat=True):
        ''' CSV convertion of a polars data frame into the character matrix of self.matlin object
        '''
        self.df  = pd.read_csv(ifn).set_index('cbc')
        self.df = self.df.fillna(self.non_inf['NA'], inplace = False)
        self.intid_len = self.df.shape[1]
        f_wt = self.df==".."
        self.df[f_wt] = self.non_inf["wt"]

        self.ibar_to_bint = dict([(f'bint{i}',v) for i,v in enumerate(self.df.columns)])
        self.df.columns = [f'bint{i}' for i in range(self.df.shape[1])]

        if encode_mat:
            self.__encode_mat()

    def __read_table(self, ifn, omit_mm = True):
        tt  = pd.read_table(ifn, header=None)
        # TODO add headers to files since the pre-processing step
        tt.columns = np.hstack(['cbc', [f'bint{i}' for i in range(tt.shape[1]-1)]])
        tt.set_index('cbc', inplace = True, drop = True, append=False)
        if omit_mm:
            tt = tt.apply(lambda x : [self.__omit_mm(i) for i in x] ) 

        f_na = tt==".NA."
        f_ex = tt=="..exx"
        f_wt = tt==".."
        tt[f_na] = self.non_inf["NA"]
        tt[f_wt] = self.non_inf["wt"]
        tt[f_ex] = self.non_inf["ex"]
        #tt = tt.apply(lambda x: [self.annotate_bint(i, x) for i in x], axis = 0)
        return(tt)

    def __omit_mm(self, bc):
        ''' Omits mismatch element from a barcode of the format mm.ins.del'''
        return(F'.{bcs[1]}.{bcs[2]}')

    def __encode_mat(self):
        bcs = self.df.stack()#.rank(method='dense')#.unstack()
        codes = bcs.rank(method='dense')
        bcs_dic = dict(zip(codes.drop_duplicates(), bcs.drop_duplicates()))
        self.bc = codes.unstack()

        # allele matrix: encoded indels matrix without duplicates
        self.amat = self.bc[~self.bc.duplicated()]
        self.amat.set_index(self.amat.apply(self.__concat_allele_codes, axis = 1), inplace = True)

        # allele representative
        self.ar = self.bc.apply(self.__concat_allele_codes, axis=1).to_frame()
        self.ar['cell'] = self.ar.index
        self.ar.columns = ['allele', 'cell']
        allele_children = dict([(i,[]) for i in self.amat.index])

        x = self.ar.apply(lambda x: allele_children[x['allele']].append(x['cell']), axis = 1)
        self.allele_children = allele_children

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

    def __intersect_cells(self, cells, index):
        ''' Returns a sorted list of cells that intersect a query '''
        #if cells.__class__.__name__ != 'NoneType':
        if cells is not None:
            cellsi = set(cells)
            cells =  [i for i in index if i in cellsi]
            return(cells)
        else:
            return(None)



def plot_viral_clones(x, xv, top_n = 5):
    ''' x is a matlin object \n e.g. 
    plot_viral_clones(x = o3to6 , xv = o3to6_virus, top_n = 5)'''
    with plt.rc_context({'figure.dpi':100}):
        for i in xv['vid'].value_counts().head(top_n).index.to_list():
            cells = set(xv[xv['vid'] == i]['cbc']).intersection(x.df.index)
            x.plot_mat(range(x.bc.shape[1]), cells=cells, add = True)
            x.plot_dist_heatmap(cells = cells)



def return_matlins_per_viral_clone(matlin, xv, top_n = 5, do_plot = False):
    '''
    x is a matlin object xv the viral table
    '''
    if do_plot:
        xv['vid'].value_counts().head(top_n).plot(kind ='bar')
    split_matlins = []
    viral_labels_code = []
    viral_labels = xv['vid'].value_counts().head(top_n).index.to_list() 
    for i,v in enumerate(viral_labels):
        cells = set(xv[xv['vid'] == v]['cbc']).intersection(matlin.df.index)
        print(f'{len(cells)} cells from viral id {v}')
        split_matlins.append((i, matlin.filter_wl(wl = cells, copy = True)))
        viral_labels_code.append((i,v))
    split_matlins.append(('code', viral_labels_code))
    return(dict(split_matlins))

def plott(X, truncate_mode = None, p = 30, add_trans = False, freq = False, figsize = (15,5), dpi=150, **kwargs):
    with plt.rc_context({'figure.figsize':figsize, 'figure.dpi':dpi, 'figure.autolayout':True}):
        for x in [i for i in X.values()]:
            if isinstance(x, matlin) or x.__class__.__name__ == 'matlin':
                if not x.is_clustered:
                    x.allele_distance(cores = 0)
                    x.hclust(condense = True, optimal_ordering = True)
                if freq:
                    fig2, ax2 = plt.subplots(1,2)
                    x.plot_indel_freq(axes=ax2)
                fig, ax = plt.subplots(1,3)
                
                Z = x.Z
                #Z = Z[np.arange(Z.shape[0]-1, 0, -1), :]
                def map_cluster(x):
                    if x in inters:
                        return(mm.loc[x]['leiden'])
                    else:
                        return('::')
                

                den_ax = ax[0]
                if truncate_mode is None and add_trans:
                    mm = kwargs['mm']
                    sdata = kwargs['sdata']

                    inters = x.df.index.intersection(mm.index)
                    labels = [map_cluster(i) for i in x.df.index]
                    x.dd= hierarchy.dendrogram(Z, p = p, labels = labels, ax = den_ax, truncate_mode = truncate_mode, orientation = 'left')

                    xlbls = den_ax.get_ymajorticklabels()
                    for i,ii in enumerate(xlbls):
                        cluster  = ii.get_text()
                        if cluster == '::':
                            ii.set_text('---')
                            ii.set_color('#ffffff')
                        else:
                            ii.set_color(sdata.uns['leiden_colors'][int(cluster)])
                            ii.set_text('---')
                    den_ax.set_yticklabels(['---' for i in xlbls])
                else:
                    x.dd= hierarchy.dendrogram(Z, p = p, ax = den_ax, truncate_mode = truncate_mode, orientation = 'left')

                lims = den_ax.get_xlim()
                lims = np.abs(np.diff(lims)[0])
                print(lims)
                den_ax.set_xlim(left = lims*0.5, right=lims*-0.1)
                x.plot_mat(r=range(x.df.shape[1]), cells=x.df.index, add=True, ax = ax[1])
                x.plot_dist_heatmap(cells = x.df.index, ax = ax[2], aspect = 'auto')
                    
            else:
                print(x.__class__.__name__)
                print('no')
    

