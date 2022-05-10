import matplotlib.pyplot as plt
import colorhash   
from colorhash import ColorHash
import seaborn as sns


def reads_per_unit(unit, umi_counts, sample_id, png_prefix=None):

    with plt.rc_context({'figure.figsize':(8.5, 8.5)}):
        fig_title = f'Reads per single {unit}'
        fig, ax = plt.subplots(1,1)
        sns.ecdfplot(umi_counts.to_list(), ax = ax)
        ax.set_xlim(-30,500)
        ax.set_ylim(0,1.15)
        ax.set_title(fig_title)
        ax.grid()

        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()


def reads_per_umi_ecdf(ibars_df, sample_id, min_cov, png_prefix=None):
    ### ### ### ###
    fig_title = 'ibar Reads per UMI'
    fig, ax = plt.subplots(1,1)
    ax.grid()
    ### ### ### ###
    sns.ecdfplot(data=ibars_df, x='umi_reads', hue='expected', ax =ax)

    ax.set_ylim(0,1.15)
    ax.set_title(f'{sample_id}\n{fig_title}')
    ax.axvline(min_cov, color='r')

    if png_prefix is not None:
        fig.savefig(f'{png_prefix}/{sample_id}_ecdf_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
        plt.close()
 
def boxen_reads_per_ibar(ibars_df, sample_id, png_prefix=None):
    reads_ibar = ibars_df.groupby(['expected', 'ibar']).apply(lambda x: sum(x.umi_dom_reads))

    with plt.rc_context({'figure.figsize':(5,5)}):
        fig_title = 'reads per ibar'
        fig, ax = plt.subplots(1,1)
        sns.boxenplot(data=reads_ibar.reset_index(name='reads'), 
            x='expected', 
            y='reads', 
            orient='v',
            ax=ax)
    
        ax.set_title(f'{sample_id}\n{fig_title}')
        ax.grid()
        fig_title = fig_title.split('\n')[0]
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_boxen_reads_ibar_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()
    
def boxen_mols_per_ibar(ibars_df, sample_id, png_prefix=None):
    mol_ibar = ibars_df.groupby(['expected', 'ibar']).size().sort_values(ascending=False)

    with plt.rc_context({'figure.figsize':(5,5)}):
        fig_title = 'reads per ibar'
        fig, ax = plt.subplots(1,1)
        ax.grid()
        sns.boxenplot(data=mol_ibar.reset_index(name='molecules'), 
            x='expected', 
            y='molecules', 
            orient='v')

#        sns.boxenplot(data=[mol_ibar, omol_ibar], ax=ax)
    
        ax.set_title(f'{sample_id}\n{fig_title}')
        fig_title = fig_title.split('\n')[0]
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_boxen_reads_ibar_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()
 
 
def expression_levels_curve(ibars_df, sample_id, png_prefix=None):
    '''Formerly referred to as 'decay'
    '''
    fig, ax = plt.subplots(1,1)
    df = ibars_df.groupby(['expected', 'ibar']).size().reset_index(name='mols')
    fig_title = 'wl mols per ibar\n(wl vs orphans)'
    ax.set_xlim((-1, 70))
    ax.grid()
    art = df.groupby('expected').apply(lambda x: \
        ax.plot(x.sort_values(['mols'], ascending = False)['mols'].to_list())
        )
    fig.legend(art.keys())

    ax.set_title(f'{sample_id}\n{fig_title}')
    fig_title = fig_title.split('\n')[0]

    if png_prefix is not None:
        fig.savefig(f'{png_prefix}/{sample_id}_decay_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
        plt.close()
   

def macro_expression_levels_bar(ibars_df, sample_id, png_prefix=None):
    with plt.rc_context({'figure.figsize':(7.5 * 3, 7.5)}):
        fig, ax = plt.subplots(1,3, sharex=True, sharey=True)
        fig.text(0.5, 0.04, 'ibar', ha='center') #figx
        fig.text(0.04, 0.5, 'molecules', va='center', rotation='vertical') #figy

        expression_levels_bar(ibars_df, sample_id=sample_id, expected=None, device=(fig, ax[0]))
        expression_levels_bar(ibars_df, sample_id=sample_id, expected=True, device=(fig, ax[1]))
        expression_levels_bar(ibars_df, sample_id=sample_id, expected=False, device=(fig, ax[2]))

        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_decay_bar_macro', bbox_inches = "tight")
            plt.close()

def expression_levels_bar(ibars_df, sample_id, expected = None, png_prefix=None, device=None):
    if device is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig, ax = device

    df = ibars_df.groupby(['expected', 'ibar']).size().reset_index(name='mols')

    if expected is not None:
        wlisted = 'is' if expected else 'not'
        df = df[df.expected == expected]
    else:
        wlisted = 'agnostic'

    fig_title = f'{wlisted} wl mols per ibar'
    ax.set_xlim((-1, 70))
    ax.grid()

    x=df.sort_values('mols', ascending = False)
    color=[ColorHash(colorhash.colorhash.color_hash(i)).hex for i in x.ibar] 
    x=x.set_index('ibar')
    ###
    ax.bar(range(len(x)), x['mols'], color=color)

    ax.set_title(f'{sample_id}\n{fig_title}')
    fig_title = fig_title.split('\n')[0]

    if png_prefix is not None and device is not None:
        fig.savefig(f'{png_prefix}/{sample_id}_decay_bar_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
        plt.close()

def kde_mols_per_unit(ibars_df, unit, sample_id, log_base=10, png_prefix=None):
    df = ibars_df.groupby(['expected', 'ibar']).size().reset_index(name='mols')
    with plt.rc_context({'figure.figsize':(15.5, 15.5)}):
        fig = sns.displot(data=df, 
            kind='kde', 
            x='mols', 
            hue='expected', 
            log_scale=log_base)
        fig_title = f'log({log_base}) molecules per {unit}'
        fig.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
        fig.ax.grid()

        fig.ax.set_title(f'{sample_id}\n{fig_title}')
        fig_title = fig_title.replace(')','')
        fig_title = fig_title.replace('(','')
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_mol_kde_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()



def ibar_confidence(df, correction, sample_id, top_n=80, png_prefix=None):
    '''
    https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.axes_grid1.axes_grid.ImageGrid.html
    '''

    import pandas as pd
    from ogtk import UM
    from mpl_toolkits.axes_grid1 import ImageGrid
    from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

    fig_title = f'{sample_id} ibar confidence {"corr" if correction else "raw"}'
    #ibar len
    vmin = 0
    vmax = 6

    df = df.drop_duplicates().groupby(['ibar']).size().sort_values(ascending=False)
    df = df.head(top_n).reset_index(name='mols')

   
    ribar_list = df.ibar

    # TODO bee doo : generalized form using the following fn
    def return_cuadra(ibar_list, correction):
        if correction:
            ibar_list =  UM.merge_all(ibar_list, jobs=1, errors =1, mode='dist')
            ibar_list = [i[1] for i in ibar_list]

        dist_mat = UM.hdist_all(ibar_list)
        Z = linkage(dist_mat, 'ward')
        zl = leaves_list(Z)
    
        mats = [
                dist_mat,
                dist_mat[:, zl], 
                dist_mat[zl,:] ,
                dist_mat[zl,:][:, zl] 
                ]
        return(mats)

    # cuadra corr
    if correction:
        cibar_list =  UM.merge_all(ribar_list, jobs=1, errors =1, mode='dist')
        cibar_list = [i[1] for i in cibar_list]

        dist_corr = UM.hdist_all(cibar_list)
        Z_corr = linkage(dist_corr, 'ward')
        zl_cor = leaves_list(Z_corr)

        mats = [
                dist_corr,
                dist_corr[:, zl_cor], 
                dist_corr[zl_cor,:] ,
                dist_corr[zl_cor,:][:, zl_cor] 
                ]

    # cuadra raw
    else:
        dist_raw = UM.hdist_all(ribar_list)
        Z_raw = linkage(dist_raw, 'ward')
        zl_raw = leaves_list(Z_raw)
    
        mats = [
                dist_raw,
                dist_raw[:, zl_raw], 
                dist_raw[zl_raw,:] ,
                dist_raw[zl_raw,:][:, zl_raw] 
                ]

    titles = [f'', '', '', '']

    #sns.set(font_scale=0.65)
    with plt.rc_context({'figure.figsize':(5*2 , 5*2)}):
        # Set up figure and image grid
        fig = plt.figure()
        fig.text(0.5, 0.04, sample_id, ha='center') #figx
        fig.text(0.04, 0.5, '', va='center', rotation='vertical') #figy

        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                        nrows_ncols=(2,2),
                        axes_pad=0.251,
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="7%",
                        cbar_pad=0.15,
                        )

        # Add data to image grid
        for ax, mat, title in zip(grid, mats, titles):
            im = ax.matshow(mat, vmin=vmin, vmax=vmax, cmap='gist_stern')
            ax.set_title(title)
            ax.tick_params('x', top=False, bottom=True, labelbottom=False, labeltop=False)
            

        # Colorbar
        ax.cax.colorbar(im)
        ax.cax.toggle_label(True)

        #plt.tight_layout()    # Works, but may still require rect paramater to keep colorbar labels visible
        #plt.show()

        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_ibar_hamm_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()

#    with plt.rc_context({'figure.figsize':(15.5, 8.5)}):
#        
#        fig, ax = plt.subplots(1,3)
#        vmin = 0
#        vmax = 6
#       
#        cwl_ibar = df[df.wl].ibar.value_counts().index
#        cbl_ibar = df[~df.wl].ibar.value_counts().index
#        call_ibar = df.ibar.value_counts().index
#
#
#        cwl_ibar =  ogtk.UM.merge_all(cwl_ibar,  jobs=1, errors =1, mode='dist')
#        cbl_ibar =  ogtk.UM.merge_all(cbl_ibar,  jobs=1, errors =1, mode='dist')
#        call_ibar = ogtk.UM.merge_all(call_ibar, jobs=1, errors =1, mode='dist')
#
#        cwl_ibar = pd.Series(list([i[1] for i in cwl_ibar])).sort_values()
#        ubl_ibar = pd.Series([i[0] for i in cbl_ibar])
#        cbl_ibar = pd.Series([i[1] for i in cbl_ibar])
#
#        call_ibar = pd.Series(list([i[1] for i in call_ibar])) 
#
#        wl_mat = ogtk.UM.hdist_all(cwl_ibar)
#        nwl_mat = ogtk.UM.hdist_all(cbl_ibar)
#        all_mat = ogtk.UM.hdist_all(ubl_ibar)
#
#        #print(f'{png_prefix}/{name}_dist.csv') 
#        #pd.DataFrame(all_mat).to_csv(f'{png_prefix}/{name}_dist.csv', sep='\t')
#
#
#        from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
#        Z = linkage(nwl_mat, 'ward', optimal_ordering=True)
#        zl = leaves_list(Z)
#
#        ax[0].matshow(nwl_mat, vmin=vmin, vmax=vmax)
#        #dendrogram(Z, ax=ax[0])
#        ax[0].set_title(f'{sample_id} wl')
#        ax[1].matshow(nwl_mat[zl,:][:,zl], vmin=vmin, vmax=vmax)
#        ax[1].set_title('~wl')
#        ax[2].matshow(all_mat[zl,:][:,zl], vmin=vmin, vmax=vmax)
#        ax[2].set_title('all')
#
#        if png_prefix is not None:
#            fig.savefig(f'{png_prefix}/{sample_id}_ibar_chamm_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
#            plt.close()
#