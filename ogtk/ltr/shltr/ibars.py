import pysam
import regex
import ogtk 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def es(string, errors=0):
    ''' error string '''
    gen = "(XXX){e<=YYY}".replace("XXX", string).replace("YYY", str(errors))
    return(gen)

def ibar_regex(errors, rc = False):
    '''Returns a regex object capable of capturing ibars
    '''
    anchor_pam = "GGGTTAGAGCTAGA"
    anchor_u6  = "AGGACGAAACACC"
    anchor_pad = "AATAGCAAGTTAA"

    if not rc:
        anatomy = "(?P<spacer>.+)"+es(anchor_pam, errors)+"(?P<ibar>.{6})"+es(anchor_pad, errors)+"(.+)"
        reanatomy = regex.compile(anatomy)
        return(reanatomy)
    else:
        anchor_pam = ogtk.UM.rev_comp(anchor_pam)
        anchor_u6 =  ogtk.UM.rev_comp(anchor_u6)
        anchor_pad = ogtk.UM.rev_comp(anchor_pad)
        print(anchor_pam, anchor_u6, anchor_pad)
        anatomy = "(?P<spacer>.+)"+es(anchor_pam, errors)+"(?P<ibar>.{6})"+es(anchor_pad, errors)+"(.+)"
        reanatomy = regex.compile(anatomy)
        return(reanatomy)



def return_ibar(seq, rg, rc = True):
    ''' Helper function that returns ibars given a string and a regex object
    '''
    seq = ogtk.UM.rev_comp(seq) if rc else seq
    match = rg.match(seq)
    if match:
        return(match.group('ibar'))
    else:
        return(None)


def plot_ibars_pc(name, tbx_ifn, cell_bcs, end = 1000, rc = True):
    ''' rc = do rev comp
    '''
    import pysam

    rg = ibar_regex(0)
    tbx = pysam.TabixFile(tbx_ifn)

    ibars = []
    ibars_st = []
    missed = []
    seq_field = 5
    for cell in cell_bcs[0:end]:
        try:
            ff = [i for i in tbx.fetch(cell.replace("-1", ""), parser = pysam.asTuple())]
            ss = [return_ibar(i[seq_field], rg, rc) for i in ff]
            ss = set(ss)
            ibars.append(len(ss))
            ibars_st.append(list(ss))
        except ValueError:
            missed.append(cell)
            #pass
    with plt.rc_context({'figure.dpi':100, 'figure.figsize':(15,5)}):
        fig, axes = plt.subplots(1,3)
        sns.histplot(data=ibars, ax=axes[0])
        sns.boxenplot(data=ibars, ax=axes[1])
        axes[2].plot(pd.Series(np.hstack(ibars_st)).value_counts().tolist())
        axes[2].set_xlim((0,100))
    print(f'{len(missed)/len(cell_bcs):0.2f}')
    return((ibars, ibars_st))


def parse_tabix_fastq(tbx_lst):
    '''
        |readid | start | end  | cbc | umi | seq | qual|
    '''
    rid_f = 0
    cbc_f = 3
    seq_f = 5
    qual_f = 6
    #ibar = return_ibar(tbx_lst[seq_f], rg, rc)
    readid = tbx_lst[rid_f]
    seq = tbx_lst[seq_f]
    qual = tbx_lst[qual_f]
    NL = '\n'
    fq_str = f"{readid}{NL}{seq}{NL}+{NL}{qual}{NL}"
    return(fq_str)


def extract_ibars_alleles_from_tabixed_fastq(sample_id, tbx_ifn,  cell_bcs, ibar_wl, rc = True, richf = False):
    ''' Derive dominant sequence (allele) for each molecule while controlling for ibar.
    Fastq files must be tabixed first. rich formatting is supported but disabled by default.

    To facilitate the approach all the analysis assumes that sequences are being read in the same strand as pol iii transcription ocurrs. 
    cell/ibar/umi/read 
    rc = do rev comp  -> This should change for zombie mode

    returns a dictionary that keeps track of molecules with an ibar and orphanes ones.

        return({'df':df_out, 'orphans':orphans, 'dfo':dfo_out, 'read_missed':read_missed})
    '''
    import pysam
    import rich
    from rich.progress import Progress
    from rich.console import Console


    ibar_wl = set(ibar_wl)
    rg = ibar_regex(0)
    tbx = pysam.TabixFile(tbx_ifn)
    TSO_SEQ = 'TTTCTTATATGGG'
    anatomy = f'.+({TSO_SEQ}) '

    #(umi_str, umi_seq, umi_dom_reads, umi_reads)
    header  = [ 'cell', 'ibar', 'umi', 'seq', 'umi_dom_reads', 'umi_reads']
    df = []
    dfo = []
    imissed = [] # no ibar could be fetched
    read_missed = []# no TSO was found -> considered bad read
    tmissed = []  # rare cases where tabix misses an entry 
    seq_field = 5
    nreads = 0

    console = Console(record=False, force_interactive= False)
    _context = open('rm_me', 'wt') if not richf else Progress(
        f"[yellow]{sample_id}[/]",
                rich.progress.TimeRemainingColumn(),
                rich.progress.BarColumn(),
                "wlisted\t{task.fields[wlisted]}|",
                "orphans\t{task.fields[orphans]}|",
                "[red]missed ibar\t{task.fields[imissed]}|[/]",
                transient=True,
                refresh_per_second=1,
                console = console)
    with _context as progress:
        if richf:
            task_cells = progress.add_task('cell', total=len(cell_bcs), wlisted="", orphans = "", imissed = "")

        cellsp = 0
        for cell in cell_bcs:
            cellsp +=1
            if richf:
                progress.update(task_cells, 
                    advance = 1, 
                    wlisted=len(df), 
                    orphans=len(dfo), 
                    imissed=f'{len(imissed)} {len(read_missed)} {len(imissed)/(1e-10+len(read_missed)):0.2f}' 
                )
            
            # define the readsets per valid ibar into a dictionary
            irsets = {}
            orphans = {}

            try:
                query = [i for i in tbx.fetch(cell.replace("-1", ""), parser = pysam.asTuple())]
                for i in query:
                    nreads +=1 
                    umi = i[4] 

                    seq = ogtk.UM.rev_comp(i[seq_field]) if rc else i[seq_field]
                    qual = i[6]
                    ibar = return_ibar(i[seq_field], rg, rc)
                    #TODO what happens to reads with no detectable ibar
                    if ibar is None:
                        imissed.append(umi)
                        continue

                    anatomy = f'.+({TSO_SEQ})(?P<read>.+){ibar}.+'
                    match = regex.search(anatomy, seq)

                    if match:
                        seq = match.group('read')
                    else:
                        read_missed.append(seq)
                        continue

                    if ibar in ibar_wl:
                        if ibar not in irsets.keys():
                            #print(ibar)
                            irsets[ibar] = ogtk.UM.Read_Set(name = ibar)
                        irsets[ibar].add_umi(umi, seq, qual = qual)
                    else:
                        if ibar not in orphans.keys():
                            orphans[ibar] = ogtk.UM.Read_Set(name = ibar)
                        orphans[ibar].add_umi(umi, seq, qual = qual)

                # recover dominant sequence for each ibar
                # store as a df for both ibar groups
                igroups = [irsets, orphans]
                dfms = [df, dfo]

                for rsetg, dfoc in zip(igroups, dfms):
                    for ibar in rsetg.keys():
                        ogtk.UM.do_fastq_pileup(readset = rsetg[ibar])
                        # extract the dominant sequnce and collate mol stats 
                        # (umi_str, umi_seq, umi_dom_reads, umi_reads)
                        mol_cons = [[cell, ibar, a,b,c,d] for a,b,c,d in rsetg[ibar].return_consensuses_tuple()]

                        for ustat in mol_cons:
                            dfoc.append(ustat)

            except ValueError:
                # some times the tabix file is missing an entry specified in the cbcs
                tmissed.append(cell)

        df_out = pd.DataFrame(df, columns = header)
        dfo_out = pd.DataFrame(dfo, columns = header)

        found = ibar_wl.intersection(set(df_out['ibar']))
        print(f'{len(ibar_wl)} expected ; found {len(found)}/{len(ibar_wl)} = {100*len(found)/len(ibar_wl):0.2f}%')
        if richf:
            rich.print(f'[red] ibars missed {len(imissed)}+ reads missed {len(read_missed)}/nreads {nreads}  = {(len(imissed)+len(read_missed))/nreads*100:0.2f}%')
        return({'df':df_out, 'dfo':dfo_out})

def genotype_ibar_clone(
    sample_id, 
    clone, 
    r1, 
    single_molecule,
    single_end,
    q=0.85, 
    nmols=None, 
    min_cov=None,
    end=int(1e5), 
    force=True, 
    richf=True,
    genotyping_db=None,
    do_plots=True,
    h5ad_path=None,
    png_prefix=None):
    ''' 
    genotyping_db = '/local/users/polivar/src/artnilet/resources/results/genotyping.txt'
    samples_tab = pd.read_csv('/local/users/polivar/src/artnilet/conf/bulk_RNA_clones_samples.txt', delimiter='\t')
    
    '''
    if h5ad_path is None and not single_molecule:
        raise  ValueError('if not single-molecule (thus single cell) a h5d adata paht must be provided')
 
    if do_plots and png_prefix is None:
        richf = False

    from rich import print
   

    unit ="molecule" if single_molecule else "cell"  

    print(f'[red] indexing\n{r1}')
    print(f'single mol {single_molecule}')
    ## index original reads
    if single_end:
        print(f'single_end {single_end}')
        tabix_fn = ogtk.utils.tabulate_umified_fastqs(r1 = r1, cbc_len=0, umi_len=25, single_molecule=single_molecule, end = end, force = force)
    else:
        print(f'single_end {single_end}')
        tabix_fn = ogtk.utils.tabulate_paired_umified_fastqs(r1 = r1, cbc_len=16, umi_len=10, single_molecule=single_molecule, end = end, force = force)

    ## extract all reads and, on the fly, and count reads per umis or cells
    ## depends on what is on field 3

    #    |readid | start | end  | cbc | umi | seq | qual|
    umi_counts = pd.read_csv(tabix_fn, compression='gzip', sep = '\t', header = None)[3].value_counts(normalize= False)
    if min_cov is None:
        min_cov = np.quantile(umi_counts, q)
 
    #cell_bcs = umi_counts[umi_counts>=min_cov].index.to_list()
    if True:
        if not single_molecule:
            # import cell barcodes from cleaned andata's index
            # TODO change thi to something less embarrassing
            import anndata as ad
            import gc
            adata = ad.read_h5ad(h5ad_path)
            cell_bcs=[i.split('-')[0] for i in adata[adata.obs.batch == sample_id].obs.index.to_list()]
            del(adata)
            gc.collect()
        else:
            cell_bcs = umi_counts[umi_counts>=min_cov].index.to_list()

    print(f'quantile {q} represents a min of {min_cov} reads')
    print(f'[red]{unit}s[/red] screened {len(cell_bcs)}')

    
    if clone is not None and genotyping_db is not None:
        gen = pd.read_csv(genotyping_db, sep= '\t')
        gen_mask =gen.clone == clone.upper() 
        print(f'looking for clone {clone.upper()} in {sum(gen_mask)} rows')
        if sum(gen_mask) ==0:
            #TODO add warning
            print(f'[bold red] Warning: clone {clone.upper()} is not present in genotyping db')
            ibar_wl=['NULL']
        else:
            ibar_wl = gen[gen_mask]['ibar_str'].drop_duplicates()
    else:
        ibar_wl = ['NULL']
        
    rc = not single_end # False if single_end
    ibars_dict = extract_ibars_alleles_from_tabixed_fastq(
        sample_id, 
        tabix_fn,
        cell_bcs,
        ibar_wl, 
        richf = richf,
        rc = rc)
    

    if len(ibars_dict['df']) ==0:
        print(f'no ibar found for clone {clone} using only orphaned ibars')
        ibars_dict['df'] = ibars_dict['dfo']
        
    if not do_plots:
            return(ibars_dict)
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


    return(repe(ibars_dict, unit = unit, clone = clone, min_cov = min_cov, sample_id=sample_id, png_prefix=png_prefix))

    return(ibars_dict) 

def repe(ibars_dict, unit, clone, min_cov, sample_id, png_prefix=None):
    import colorhash   
    from colorhash import ColorHash


    ### plotting ###
    sns.set_context("talk")
    plt.rc('axes', axisbelow=True)

    plt.rcParams['figure.dpi'] = 200
    if png_prefix is not None:
        import os
        if not os.path.exists(png_prefix):
            os.makedirs(png_prefix)

    ### ### ### ###
    fig_title = 'ibar Reads per UMI'
    fig, ax = plt.subplots(1,1)
    ax.grid()
    ### ### ### ###
    sns.ecdfplot(ibars_dict['df']['umi_reads'], ax=ax)
    ax.set_ylim(0,1.15)

    if len(ibars_dict['dfo']) >0 and clone is not None:
        sns.ecdfplot(ibars_dict['dfo']['umi_reads'], ax=ax)
    ax.set_title(f'{sample_id}\n{fig_title}')
    ax.axvline(min_cov, color='r')

    if png_prefix is not None:
        fig.savefig(f'{png_prefix}/{sample_id}_ecdf_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
        plt.close()
 
    ### ### ### ###
    if len(ibars_dict['dfo']) >0:
        fig_title = 'Dominant reads per UMI\n(wl vs orphan)'
        fig, ax = plt.subplots(1,1)
        ax.grid()
        ax.set_axisbelow(True)

        sns.boxenplot( data = [
                    ibars_dict['df'].groupby(['ibar']).apply(lambda x: sum(x['umi_dom_reads'])), 
                    ibars_dict['dfo'].groupby(['ibar']).apply(lambda x: sum(x['umi_dom_reads']))], ax=ax)
        ax.set_title(f'{sample_id}\n{fig_title}')

        fig_title = fig_title.split('\n')[0]
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_boxen_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()

    reads_ibar = ibars_dict['df'].groupby(['ibar']).apply(lambda x: sum(x.umi_dom_reads))
    mol_ibar = ibars_dict['df'].groupby(['ibar']).size().sort_values(ascending=False) 
    
    oreads_ibar = ibars_dict['dfo'].groupby(['ibar']).apply(lambda x: sum(x.umi_dom_reads))
    omol_ibar = ibars_dict['dfo'].groupby(['ibar']).size().sort_values(ascending=False) 
    omol_ibar = omol_ibar.head(len(mol_ibar))
    
    ###
    with plt.rc_context({'figure.figsize':(5,5)}):
        fig_title = 'mols per ibar\n(wl vs orphans)'
        fig, ax = plt.subplots(1,1)
        ax.grid()
        sns.boxenplot(data=[mol_ibar, omol_ibar], ax=ax)
    
        ax.set_title(f'{sample_id}\n{fig_title}')
        fig_title = fig_title.split('\n')[0]
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_boxen_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()
    
        fig_title = 'mols per ibar\n(wl vs orphans)'
    
            
    df = pd.DataFrame({'wl':True, 'mols':mol_ibar})
    df = pd.concat([df, pd.DataFrame({'wl':False, 'mols':omol_ibar})])
    df = df.reset_index()
    ###

    fig, ax = plt.subplots(1,1)

    fig_title = 'wl mols per ibar\n(wl vs orphans)'
    ax.set_xlim((-10,80))
    ax.grid()
    df.groupby('wl').apply(lambda x: \
        ax.plot(x.sort_values(['mols'], ascending = False)['mols'].to_list())
        )


    ax.set_title(f'{sample_id}\n{fig_title}')
    fig_title = fig_title.split('\n')[0]
    if png_prefix is not None:
        fig.savefig(f'{png_prefix}/{sample_id}_decay_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
        plt.close()
    ###

    ###
   
    def bar_lambda(df, wlisted = None):
        fig, ax = plt.subplots(1,1)

        if wlisted is None:
            wlisted = '' if df.wl.iloc[0] else 'not'

        fig_title = f'{wlisted} wl mols per ibar'
        ax.set_xlim((-10,80))
        ax.grid()
    
        x=df.sort_values('mols', ascending = False)
        color=[ColorHash(colorhash.colorhash.color_hash(i)).hex for i in x.ibar] 
        x=x.set_index('ibar')
        ax.bar(range(len(x)), x['mols'], color=color)

        ax.set_title(f'{sample_id}\n{fig_title}')
        fig_title = fig_title.split('\n')[0]
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_decay_bar_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()

    bar_lambda(df, wlisted ='agnostic')
    df.groupby('wl').apply(bar_lambda)


    ###


    with plt.rc_context({'figure.figsize':(15.5, 15.5)}):
        fig = sns.displot(data=np.log10(df.mols), kind='kde')
        fig_title = f'molecules per {unit}'
        fig.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
        fig.fig.suptitle(f'{sample_id}\n{fig_title}')
        fig.ax.grid()

        ax.set_title(f'{sample_id}\n{fig_title}')
        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_mol_kde_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()
        
    with plt.rc_context({'figure.figsize':(15.5, 8.5)}):
        
        fig, ax = plt.subplots(1,3)

        vmin = 0
        vmax = 6


        wl_mat = ogtk.UM.hdist_all(df[df.wl].ibar)
        nwl_mat = ogtk.UM.hdist_all(df[~df.wl].ibar)
        all_mat = ogtk.UM.hdist_all(df.ibar)

        from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
        Z = linkage(wl_mat, 'ward')
        zl = leaves_list(Z)
        
        ax[0].matshow(wl_mat[zl,:][:,zl], vmin=vmin, vmax=vmax)
        ax[0].set_title(f'{sample_id} wl')
        ax[1].matshow(nwl_mat, vmin=vmin, vmax=vmax)
        ax[1].set_title('~wl')
        ax[2].matshow(all_mat, vmin=vmin, vmax=vmax)
        ax[2].set_title('all')

        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_ibar_hamm_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()

    with plt.rc_context({'figure.figsize':(15.5, 8.5)}):
        
        fig, ax = plt.subplots(1,3)
        vmin = 0
        vmax = 6
       
        cwl_ibar = df[df.wl].ibar.value_counts()

        cwl_ibar = cwl_ibar.index
        cbl_ibar = df[~df.wl].ibar.value_counts().index
        call_ibar = df.ibar.value_counts().index


        cwl_ibar =  ogtk.UM.merge_all(cwl_ibar,  jobs=1, errors =2, mode='dist')
        cbl_ibar =  ogtk.UM.merge_all(cbl_ibar,  jobs=1, errors =1, mode='dist')
        call_ibar = ogtk.UM.merge_all(call_ibar, jobs=1, errors =2, mode='dist')

        cwl_ibar = pd.Series(list([i[1] for i in cwl_ibar])).sort_values()

        ubl_ibar = pd.Series([i[0] for i in cbl_ibar])
        cbl_ibar = pd.Series([i[1] for i in cbl_ibar])

        call_ibar = pd.Series(list([i[1] for i in call_ibar])) 

        wl_mat = ogtk.UM.hdist_all(cwl_ibar)

        nwl_mat = ogtk.UM.hdist_all(cbl_ibar)
        all_mat = ogtk.UM.hdist_all(ubl_ibar)

        #print(f'{png_prefix}/{name}_dist.csv') 
        #pd.DataFrame(all_mat).to_csv(f'{png_prefix}/{name}_dist.csv', sep='\t')


        from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
        Z = linkage(nwl_mat, 'ward', optimal_ordering=True)
        zl = leaves_list(Z)

        ax[0].matshow(nwl_mat, vmin=vmin, vmax=vmax)
        #dendrogram(Z, ax=ax[0])
        ax[0].set_title(f'{sample_id} wl')
        ax[1].matshow(nwl_mat[zl,:][:,zl], vmin=vmin, vmax=vmax)
        ax[1].set_title('~wl')
        ax[2].matshow(all_mat[zl,:][:,zl], vmin=vmin, vmax=vmax)
        ax[2].set_title('all')

        if png_prefix is not None:
            fig.savefig(f'{png_prefix}/{sample_id}_ibar_chamm_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
            plt.close()

#    sns.boxenplot(data = [np.mean(all_mat[df.wl,:][], axis=0), 
#            np.mean(ogtk.UM.hdist_all(bulk_b3.ibar)[~bulk_b3.wl], axis=0)])            
    return(ibars_dict)


def df_genotype_ibar_clone(x, rootdir, quantile = 0.75, end = int(1e6), do_plots= True, png_prefix='/local/users/polivar/src/artnilet/figures/', h5ad_path = None, force = True, genotyping_db=None):
    '''Wrapper function can be mapped to a data frame  in order to automate parameter definition. 
    Ther are several expected fields in data frame: sample_id, clone, lin_bin, stage
    The stereotypical data frame would be artnilet/conf/xpdb_datain.txt

    #TODO improve logics, for example cc is now a handle for bulk when in reality it means cell culture 
    '''
    print(x)
    r1 = f'{rootdir}/datain/{x.lin_lib}'
    print(r1)
    ibars_dict = genotype_ibar_clone(
                sample_id = f'{x.sample_id}', 
                clone = x.clone,
                end=end,
                r1=r1,
                q=quantile, 
                do_plots=do_plots,
                force=force,
                single_molecule= False, # cc = cell culture
                single_end=False,
                genotyping_db=genotyping_db,
                h5ad_path = h5ad_path,
                png_prefix=png_prefix)

    return(ibars_dict)
