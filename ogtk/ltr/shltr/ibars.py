from typing import *
from colorhash import ColorHash
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
import os
import regex
import seaborn as sns
import tempfile
from pyseq_align import NeedlemanWunsch#, SmithWaterman
import ngs_tools
from functools import lru_cache

import ogtk
import ogtk.utils as ut
import ogtk.utils.log as log

from ogtk.ltr import ltr_utils
from . import plot as pt

from ogtk.utils.log import Rlogger
logger = Rlogger().get_logger()

@ogtk.utils.log.call
def reads_to_molecules(sample_id: str,
           parquet_ifn: str, 
           modality: str,
           zombie: bool=False,
           corr_dict_fn: str | None= None,
           min_reads: int=1, 
           max_reads: int=int(1e6),
           do_sampling: int | None=None,
           min_cells_per_ibar: int | None=1000,
           clone: str | None=None,
           columns_ns: str | Sequence='np_tsos',
           cache_dir: str='/local/users/polivar/src/artnilet/cache',
           force: bool=False, 
           ) -> pl.DataFrame:
    ''' Analytical routine to process UMIfied shRNA from single-cell assays (parquet format).
        If `clone` is not provided it assumes that can be found as the first element of a dash-splitted `sample_id` string. 
        Off-targets are defined by reads without an ibar pattern match.
        UMIs are filtered based on the number of reads and quality (e.g. if they contain Ns)
    '''
    
    if modality not in ['single-cell', 'single-molecule']:
        raise ValueError(f'incorrect modality. please use "single-cell" or "single-molecule"')

    spacers = load_wl(True)['spacer']

    cache_out = f'{cache_dir}/{sample_id}_r2mols.parquet'
    encoded_reads_parquet = f'{cache_dir}/{sample_id}_encoded_reads.parquet'

    if clone is None:
        clone = sample_id.split('_')[0]

    if os.path.exists(cache_out) and not force:
        logger.info(f'loading from cache {cache_out}')
        rdf = pl.read_parquet(cache_out)
        return rdf
    else:
        # generate a reads data frame 
        if os.path.exists(encoded_reads_parquet) and not force:
            logger.info(f"loading {encoded_reads_parquet=}")
            rdf = pl.read_parquet(encoded_reads_parquet)
        else:
            rdf = encode_reads_with_dummies(
                     parquet_ifn = parquet_ifn,
                     do_sampling = do_sampling,
                     zombie=zombie
                )

            rdf.write_parquet(encoded_reads_parquet)

        tot_reads = rdf.shape[0]
        if modality == 'single-cell':
            rdf = rdf.with_columns(pl.col('cbc')+"-1")

        # correct cbc if path to dictionary is provided
        #if corr_dict_fn is not None:
        #    logger.info('correcting with dic')
        #    rdf = ogtk.utils.cd.correct_cbc_pl(rdf, ogtk.utils.cd.load_corr_dict(corr_dict_fn))

        # CPU-intensive
        # TODO add thread control
        
        rdf = ibar_reads_to_molecules(rdf, modality=modality)
        rdf = extract_spacer(rdf, spacers) #type: ignore
        rdf = encode_wt(rdf, spacers) #type: ignore
        rdf = count_essential_patterns(rdf)

        # plot noise and return noise-related data frames
        noise_df, noise_df_pc = \
                noise_spectrum(
                        sample_id,
                        rdf,
                        index = ['np_dss'],
                        columns=columns_ns)    

        path_noise_out = f'{cache_dir}/{sample_id}_noise_spectrum.parquet'
        path_noisep_out = f'{cache_dir}/{sample_id}_noise_spectrum_pc.parquet'

        logger.io(f"writing {path_noise_out}")
        logger.io(f"writing {path_noisep_out}")

        noise_df.write_parquet(path_noise_out)
        noise_df_pc.write_parquet(path_noisep_out)

        logger.step('dropping noise pattern fields')
        rdf = rdf.drop(cs.by_name("^np_.+$"))

        rdf = mask_wt(rdf)

        plot = False
        if plot:
            pt.plot_sibling_noise(rdf)

        #gather stats
        # how to determine pc_offt TODO
        #pc_offt = rdf.filter(~pl.col("valid_ibar")).shape[0]/rdf.shape[0]
        tot_umis = rdf.select('umi').n_unique()

        #logger.info(f'{pc_offt=:.2%}')
        logger.info(f'{tot_umis=}')
        logger.info(f'{tot_reads=}')

        umis_allele_group = ['cbc', 'raw_ibar', 'seq'] if modality == "single-cell" else ['raw_ibar', 'seq']
        n_sib_group = ['cbc', 'raw_ibar'] if modality == "single-cell" else ['raw_ibar']

        rdf = (rdf
                .with_columns(pl.lit(sample_id).alias('sample_id'))
                .with_columns(pl.lit(clone).alias('clone'))
                .filter(pl.col('raw_ibar').is_not_null())
                .filter(pl.col('raw_ibar').str.contains('N').not_())
                .filter(pl.col('spacer').str.contains('N').not_())
                .filter(pl.col('umi').str.contains('N').not_())
                .filter(pl.col('umi_dom_reads')<=max_reads)
                .filter(pl.col('umi_dom_reads')>=min_reads)
                .with_columns(pl.col('umi').n_unique().over(umis_allele_group).alias('umis_allele'))
                .with_columns(pl.col('seq').n_unique().over(n_sib_group).alias('n_sib'))
               )
        # consolidates global QC metrics into df; these are agnostic to downsampling
        rdf = qc_stats(rdf, tot_umis, tot_reads)
        path_qc_stats= f'{cache_dir}/{sample_id}.qc_stats.parquet'

        logger.io(f"writing {path_qc_stats}") #pyright: ignore
        rdf.select(cs.by_name("^(qc_.+|sample_id)$")).unique().write_parquet(path_qc_stats)
        rdf = rdf.drop(cs.by_name("^qc_.+$"))

        rdf.write_parquet(cache_out)
        logger.io(cache_out)
    return(rdf)    

@ogtk.utils.log.call
def qc_stats(df, tot_umis, tot_reads):
    ''' helper function that consolidates global (downsampling agnostic) QC metrics into a data frame
    '''
    return(
            df
           #.with_columns(pl.lit(pc_offt).alias('qc_pc_offt'))
           .with_columns(pl.lit(tot_umis).alias('qc_tot_umis'))
           .with_columns(pl.lit(tot_reads).alias('qc_tot_reads'))
           )

def error_string(string, errors=0, error_type='e'):
    ''' es = error string 

        error_type is a string that determined whether general errors (`e`), deletions (`d`) or substitutions
        returns a string of the form(PATTERN){e<=ERRORS}
    ''' 
    if error_type not in ['e', 'd', 's']:
        raise  ValueError('errors_type must be either of: "e", "s", "d"')
    gen = "(PATTERN){ERROR_TYPE<=ERRORS}".replace("PATTERN", string).replace("ERRORS", str(errors)).replace("ERROR_TYPE", error_type)
    return(gen)

def ibar_regex(errors, rc = False, zombie=False):
    '''Returns a regex object capable of capturing ibars
        anchor_pam = "GGGTTAGAGCTAGA"
        anchor_u6  = "AGGACGAAACACC"
        anchor_u6_1  = "GGCTTTATATATC"
        anchor_pad = "AATAGCAAGTTAA"
        anatomy = "(?P<spacer>.+)"+error_string(anchor_pam, errors)+"(?P<ibar>.{6})"+error_string(anchor_pad, errors)+"(.+)"
    '''
    anchor_pam = "GGGTTAGAGCTAGA"
    anchor_u6_1  = "GGCTTTATATATC"
    anchor_u6  = "AGGACGAAACACC"
    anchor_pad = "AATAGCAAGTTAA"

    if zombie:
        # here we should find the u6 so the anatomy would look like
        anchor_pam_short = "GAGCTAGA"
        anatomy = (
                "(?P<upstream>.+)"+error_string(anchor_u6_1, errors)+
                "(?P<spacer>.+)"+error_string(anchor_pam_short, errors)+
                "(?P<ibar>.{6})"+error_string(anchor_pad, errors)+"(.+)"
        )
        #print(anatomy)
    if not rc:
        #anatomy = "(?P<spacer>.+)"+error_string(anchor_pam, errors)+"(?P<ibar>.{6})"+error_string(anchor_pad, errors)+"(.+)"
        anatomy = "(?P<spacer>.+)"+error_string(anchor_pam, errors)+"(?P<ibar>.{6})"+error_string(anchor_pad, errors)+"(.+)"
        #print(anatomy)
        reanatomy = regex.compile(anatomy)
        return(reanatomy)
    else:
        anchor_pam = ogtk.UM.rev_comp(anchor_pam)
        anchor_u6 =  ogtk.UM.rev_comp(anchor_u6)
        anchor_pad = ogtk.UM.rev_comp(anchor_pad)
        print(anchor_pam, anchor_u6, anchor_pad)
        #anatomy = "(?P<spacer>.+)"+error_string(anchor_pam, errors)+"(?P<ibar>.{6})"+error_string(anchor_pad, errors)+"(.+)" <- this might be wrong 25.06.22
        anatomy = "(?P<spacer>.+)"+error_string(anchor_pad, errors)+"(?P<ibar>.{6})"+error_string(anchor_pam, errors)+"(.+)"
        #print(anatomy)
        reanatomy = regex.compile(anatomy)
        return(reanatomy)



def return_ibar(seq, rg, rc = True):
    '''
        Helper function that returns an ibar string     
        given a read sequence and a regex object of the type returned by `ibar_regex`
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


def extract_ibars_alleles_from_tabixed_fastq(sample_id, tbx_ifn, cell_bcs, ibar_wl, anchor_ok, rg = ibar_regex(0), rc = True, richf = False):
    
    ''' Derive dominant sequence (allele) for each molecule while controlling for ibar.
    Fastq files must be tabixed first. rich formatting is supported but disabled by default.

    To facilitate the approach all the analysis assumes that sequences are being read in the same strand as pol iii transcription ocurrs. 
    cell/ibar/umi/read 
    rc = do rev comp  -> This should change for zombie mode

    returns a data frame

        return({'expected':exp_out, 'orphans':orphans, 'unexpected':unex_out, 'read_missed':read_missed})
    '''
    import pysam
    import rich
    from rich.progress import Progress
    from rich.console import Console

    ibar_wl = set(ibar_wl)
    tbx = pysam.TabixFile(tbx_ifn)

    #(umi_str, umi_seq, umi_dom_reads, umi_reads)
    header  = [ 'cell', 'ibar', 'umi', 'seq', 'umi_dom_reads', 'umi_reads']
    expected = []
    unexpected = []
    no_ibars = []
    ibar_missed = [] # no ibar could be fetched
    read_missed = []# no TSO was found -> considered bad read
    tmissed = []  # rare cases where tabix misses an entry 
    seq_field = 5
    nreads = 0

    debug_messages =0
    #console = Console(record=False, force_interactive= False)
    console = Console(record=True, force_interactive= True)
    if len(cell_bcs) < 10:
        print(f'the number of cells provided looks too small {len(cell_bcs)=}')
    _context = open('rm_me', 'wt') if not richf else Progress(
                f"[yellow]{sample_id}[/]",
                rich.progress.TimeRemainingColumn(),
                rich.progress.BarColumn(),
                "wlisted\t{task.fields[wlisted]}|",
                "orphans\t{task.fields[orphans]}|",
                "[red]missed ibar\t{task.fields[ibar_missed]}|[/]",
                transient=True,
                refresh_per_second=1,
                console = console)
    with _context as progress:
        if richf:
            task_cells = progress.add_task('cell', total=len(cell_bcs), wlisted="", orphans = "", ibar_missed = "")

        cellsp = 0
        for cell in cell_bcs:
            cellsp +=1
            if richf:
                progress.update(task_cells, 
                    advance = 1, 
                    wlisted=len(expected), 
                    orphans=len(unexpected), 
                    ibar_missed=f'{len(ibar_missed)} {len(read_missed)} {len(ibar_missed)/(1e-10+len(read_missed)):0.2f}' 
                )
            
            # define the readsets per valid ibar into a dictionary
            irsets = {}
            orphans = {}
            no_ibars_rsets={}
            no_ibars_rsets['no_ibar'] = ogtk.UM.Read_Set(name = 'no_ibar')

            try:
                query = [i for i in tbx.fetch(cell.replace("-1", ""), parser = pysam.asTuple())]
                # traverse all reads for a given cell
                for tabix_line in query:
                    nreads +=1 
                    umi = tabix_line[4] 
                    read_str = tabix_line[seq_field]
                    seq = ogtk.UM.rev_comp(read_str) if rc else read_str
                    qual = tabix_line[6]
                    ibar = return_ibar(read_str, rg, rc)
                    #TODO what happens to reads with no detectable ibar
                    if ibar is None:
                        ibar_missed.append(umi)
                        no_ibars_rsets['no_ibar'].add_umi(umi, seq, qual = qual)
                        continue

                    anatomy = f'.+({anchor_ok})(?P<read>.+)'+error_string(f'?P<ibar>{ibar}', 1, 'd')+'.+'
                    if debug_messages<10:
                        print(anatomy)
                        debug_messages+=1
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
                # with all reads processed for a given cell
                # recover dominant sequence for each ibar
                # store as a df for both ibar groups

                # create iterable for the different ibar groups
                # and their pointers to data frames
                igroups = [irsets, orphans, no_ibars_rsets]
                dfms = [expected, unexpected, no_ibars]

                for rsetg, dfoc in zip(igroups, dfms):
                    for ibar in rsetg.keys():
                        rset=rsetg[ibar] 
                        # Generate consensus sequences of sets of reads, grouped by UMI 
                        # (wrongly named function do_fastq_pileup) 
                        ogtk.UM.do_fastq_pileup(readset=rset)
                        # extract the dominant sequnce and collate mol stats 
                        # (umi_str, umi_seq, umi_dom_reads, umi_reads)
                        mol_cons = [[cell, ibar, a,b,c,d] for a,b,c,d in rset.return_consensuses_tuple()]

                        for ustat in mol_cons:
                            dfoc.append(ustat)

            except ValueError:
                # some times the tabix file is missing an entry specified in the cbcs
                tmissed.append(cell)

        # print final stats
        # TODO: save them somewhere?

        if richf:
            progress.update(task_cells, 
                advance = 1, 
                wlisted=len(expected), 
                orphans=len(unexpected), 
                ibar_missed=f'{len(ibar_missed)} {len(read_missed)} {len(ibar_missed)/(1e-10+len(read_missed)):0.2f}' 
            )
        (
            print(f'{len(expected)=}\n'
                +f'{len(unexpected)=}\n'
                +f'{len(ibar_missed)=}\n'
                +f'{len(read_missed)=}\n'
                +f'missed ibars {len(ibar_missed)/(1e-10+nreads):0.2f}\n'
                #+f'reads missed/nreads {len(read_missed)}/{nreads}  = {(len(ibar_missed)+len(read_missed))/nreads*100:0.2f}%'
            )
        )
        exp_out = pd.DataFrame(expected, columns = header)
        unexp_out = pd.DataFrame(unexpected, columns = header)
        no_ibars_out = pd.DataFrame(no_ibars, columns = header)

        exp_out['sample_id'] = sample_id
        unexp_out['sample_id'] = sample_id
        no_ibars_out['sample_id'] = sample_id

        exp_out['expected'] = True
        unexp_out['expected'] = False
        no_ibars_out['expected'] = False


        found = ibar_wl.intersection(set(exp_out['ibar']))
        print(f'{len(ibar_wl)} expected ; found {len(found)}/{len(ibar_wl)} = {100*len(found)/len(ibar_wl):0.2f}%')
        if richf:
            rich.print(f'[red] ibars missed {len(ibar_missed)} + reads missed/nreads {len(read_missed)}/{nreads}  = {(len(ibar_missed)+len(read_missed))/nreads*100:0.2f}%')
        #return({'expected':exp_out, 'unexpected':unexp_out})
        return(pd.concat([exp_out, unexp_out, no_ibars_out]))

def genotype_ibar_clone(
    sample_id, 
    clone, 
    r1, 
    single_molecule,
    single_end,
    quantile=0.85, 
    nmols=None, 
    min_cov=None,
    end=int(1e5), 
    force=True, 
    richf=True,
    genotyping_db=None,
    do_plots=True,
    h5ad_path=None,
    png_prefix=None, **kwargs):
    ''' 
    Returns a UMI-controlled table for  
    It is a wrapper function of `extract_ibars_alleles_from_tabixed_fastq`
    genotyping_db = '/local/users/polivar/src/artnilet/resources/results/genotyping.txt'
    samples_tab = pd.read_csv('/local/users/polivar/src/artnilet/conf/bulk_RNA_clones_samples.txt', delimiter='\t')
    
    '''
    if h5ad_path is None and not single_molecule:
        raise  ValueError('if not single-molecule (and thus, single cell) a h5d adata path must be provided')
 
    if do_plots and png_prefix is None:
        print('automatically disabled rich formatting')
        richf = False

    from rich import print

    unit ="molecule" if single_molecule else "cell"  

    print(f'[red] indexing\n{r1}')
    print(f'single mol {single_molecule}')

    ## index a small sample of original reads
    if single_end:
        print(f'single_end {single_end}')
        tabix_fn = ut.tabulate_umified_fastqs(r1 = r1,
                     cbc_len=0,
                     umi_len=25,
                     single_molecule=single_molecule,
                     end=end,
                     force=force,
                     comparable=True)
    else:
        print(f'single_end {single_end}')
        tabix_fn = ut.tabulate_paired_umified_fastqs(r1 = r1,
                     cbc_len=16,
                     umi_len=10,
                     single_molecule=single_molecule,
                     end=end,
                     force=force,
                     comparable=True)
    
    if not single_molecule:
        # import cell barcodes from cleaned andata's index
        # TODO change thi to something less embarrassing
        import anndata as ad
        adata_obs = ad.read_h5ad(h5ad_path).obs
        cell_bcs=[i.split('-')[0] 
                  for i in adata_obs[adata_obs['batch'] == sample_id].index.to_list()
                  ]
    else:
        # TODO add support for plotting what it means a given cutoff in the 
        # read count per umi distribution
        # e.g by calling pt.reads_per_unit()
        if min_cov is None:
            cell_bcs, min_cov = filter_umis_quantile_tabix(tabix_fn=tabix_fn, quantile=quantile)
            print(f'quantile {quantile} represents a min of {min_cov} reads')
            print(f'[red]{unit}s[/red] screened {len(cell_bcs)}')
        else:
            cell_bcs = filter_umis_cov_tabix(tabix_fn=tabix_fn, min_cov=min_cov)
   
    if clone is not None and genotyping_db is not None:
        gen = pd.read_csv(genotyping_db, sep= '\t')
        gen_mask =gen.clone == clone.upper() 
        print(f'looking for clone {clone.upper()} in {sum(gen_mask)} rows')
        if sum(gen_mask)==0:
            #TODO add warning
            print(f'[bold red] Warning: clone {clone.upper()} is not present in genotyping db')
            ibar_wl=['NULL']
        else:
            ibar_wl = gen[gen_mask]['ibar_str'].drop_duplicates()
    else:
        ibar_wl = ['NULL']
        
    rc = not single_end # False if single_end

    if 'zombie' in kwargs:
        zombie = kwargs['zombie']
    else:
        zombie = False

    if zombie:
        rc = False
        print(":zombie:")
    else:
        print(":dna:")

    print(f'{rc=}')
    rg = ibar_regex(errors = 0, zombie=zombie)
    # validaton sequences
    # These are sequences to grant that the origin of the read is bona fide (e.g. real TSO event)
    # TODO: verify that this still holds true for zombie mode. For guides it makes sense since the TSO
    # generally covered on the R2
    TSO_SEQ = 'TTTCTTATATGGG'
    ZOM_SEQ = 'AACTTGAAAGTAT'

    if len(cell_bcs) <10:
        import warnings
        warnings.warning("The number of cells detected is too low for {sample_id}. Returning None", DeprecationWarning)
        return(pd.DataFrame())
    ibars_df = extract_ibars_alleles_from_tabixed_fastq(
        sample_id=sample_id, 
        tbx_ifn=tabix_fn,
        cell_bcs=cell_bcs,
        ibar_wl=ibar_wl, 
        anchor_ok= TSO_SEQ if not zombie else ZOM_SEQ,
        richf = richf,
        rg = rg,
        rc = rc)
    
    return(ibars_df)

def generate_qc_plots(ibars_df, unit,  min_cov, sample_id, png_prefix=None):
    import colorhash   

    #df = pd.DataFrame({'wl':True, 'mols':mol_ibar})
    #df = pd.concat([df, pd.DataFrame({'wl':False, 'mols':omol_ibar})])
    #df = df.reset_index()
    ####


    ### plotting ###
    sns.set_context("talk")
    plt.rc('axes', axisbelow=True)

    plt.rcParams['figure.dpi'] = 150

    if png_prefix is not None:
        import os
        if not os.path.exists(png_prefix):
            os.makedirs(png_prefix)

    ###
    ###

    ###
    pt.reads_per_umi_ecdf(ibars_df, sample_id=sample_id, min_cov = min_cov, png_prefix=png_prefix)

    ###
    pt.boxen_reads_per_ibar(ibars_df, sample_id=sample_id, png_prefix=png_prefix)
 
    ###
    pt.boxen_mols_per_ibar(ibars_df, sample_id=sample_id, png_prefix=png_prefix)
            
    ###
    pt.expression_levels_curve(ibars_df, sample_id=sample_id, png_prefix=png_prefix)
    ###
    pt.macro_expression_levels_bar(ibars_df, sample_id=sample_id, png_prefix=png_prefix)
    #pt.expression_levels_bar(ibars_df, sample_id=sample_id, expected=None, png_prefix=png_prefix)
    #pt.expression_levels_bar(ibars_df, sample_id=sample_id, expected=True, png_prefix=png_prefix)
    #pt.expression_levels_bar(ibars_df, sample_id=sample_id, expected=False, png_prefix=png_prefix)
    ###


    ###

    pt.kde_mols_per_unit(ibars_df, unit=unit, sample_id=sample_id, png_prefix=png_prefix)
    pt.ibar_confidence(ibars_df, correction = True, sample_id=sample_id, png_prefix=png_prefix)
    pt.ibar_confidence(ibars_df, correction = False, sample_id=sample_id, png_prefix=png_prefix)
#    with plt.rc_context({'figure.figsize':(15.5, 8.5)}):
#        
#        fig, ax = plt.subplots(1,3)
#
#        vmin = 0
#        vmax = 6
#
#
#        wl_mat = ogtk.UM.hdist_all(df[df.wl].ibar)
#        nwl_mat = ogtk.UM.hdist_all(df[~df.wl].ibar)
#        all_mat = ogtk.UM.hdist_all(df.ibar)
#
#        from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
#        Z = linkage(wl_mat, 'ward')
#        zl = leaves_list(Z)
#        
#        ax[0].matshow(wl_mat[zl,:][:,zl], vmin=vmin, vmax=vmax)
#        ax[0].set_title(f'{sample_id} wl')
#
#        ax[1].matshow(nwl_mat, vmin=vmin, vmax=vmax)
#        ax[1].set_title('~wl')
#
#        ax[2].matshow(all_mat, vmin=vmin, vmax=vmax)
#        ax[2].set_title('all')
#
#        if png_prefix is not None:
#            fig.savefig(f'{png_prefix}/{sample_id}_ibar_hamm_{fig_title.lower().replace(" ", "_")}', bbox_inches = "tight")
#            plt.close()
#

##    sns.boxenplot(data = [np.mean(all_mat[df.wl,:][], axis=0), 
##            np.mean(ogtk.UM.hdist_all(bulk_b3.ibar)[~bulk_b3.wl], axis=0)])            


def df_genotype_ibar_clone(x, outdir, rootdir, quantile = 0.75, end = int(1e6), do_plots= True, png_prefix='/local/users/polivar/src/artnilet/figures/', h5ad_path = None, force = True, force_df = True,  genotyping_db=None):
    '''Wrapper function can be mapped to a data frame in order to automate parameter definition. 
    Ther are several expected fields in data frame: sample_id, clone, lin_bin, stage
    The stereotypical data frame would be artnilet/conf/xpdb_datain.txt

    force = force the generation of a tabix file for the reads
    #TODO improve logics, for example cc is now a handle for bulk when in reality it means cell culture 
    '''
    import os
    print(x)
    r1 = f'{rootdir}/datain/{x.lin_lib}'
    out_tab = f'{outdir}/{x.sample_id}_sc_ibars.hdf'

    if os.path.exists(out_tab) and not force_df:
        print(f'loading pre-computed\n{out_tab}')
        #ibars_df = pd.read_csv(out_tab, sep='\t')
        ibars_df = pd.read_hdf(out_tab, key='df')
    else:
        print(r1)
        ibars_df = genotype_ibar_clone(
                    sample_id = f'{x.sample_id}', 
                    clone = x.clone,
                    end=end,
                    r1=r1,
                    quantile=quantile, 
                    do_plots=do_plots,
                    force=force,
                    single_molecule=False,
                    single_end=False,
                    genotyping_db=genotyping_db,
                    h5ad_path = h5ad_path,
                    png_prefix=png_prefix)
        #ibars_df.to_csv(out_tab, sep='\t', index=False )
        ibars_df.to_hdf(out_tab, key='df', mode='w')

    if not do_plots:
            return(ibars_df)
    generate_qc_plots(ibars_df, unit = 'cell', min_cov = 1, sample_id=x.sample_id, png_prefix=png_prefix)

    return(ibars_df)

def filter_umis_quantile_tabix(tabix_fn, quantile):
    ## extract all reads and, on the fly, and count reads per umis or cells
    ## depends on what is on field 3

    #    |readid | start | end  | cbc | umi | seq | qual|
    import pandas as pd
    umi_counts = pd.read_csv(tabix_fn, compression='gzip', sep = '\t', header = None)[3].value_counts(normalize= False)
    min_cov = np.quantile(umi_counts, quantile)
    filtered_umis = umi_counts[umi_counts>=min_cov].index.to_list()
    
    # plot filtering QCs
    pt.reads_per_unit(unit, umi_counts, sample_id, png_prefix)
 
    return((min_cov, filtered_umis))

def filter_umis_cov_tabix(tabix_fn, min_cov):
    #    |readid | start | end  | cbc | umi | seq | qual|
    import pandas as pd
    umi_counts = (
        pd.read_csv(tabix_fn, compression='gzip', sep = '\t', header = None)[3]
        .value_counts(normalize= False)
    )   
    filtered_umis = umi_counts[umi_counts>=min_cov].index.to_list()
    return((min_cov, filtered_umis))

def guess_ibars_mx(mat, verbose: bool= False)-> Sequence:
    ''' Corrects detected ibars using the ranking method
    '''
    import pdb
    col_sum = np.sum(mat, axis=0)
    row_sum = np.sum(mat, axis=1)

    top_cols = np.flip(np.argsort(col_sum))
    top_rows = np.flip(np.argsort(row_sum))    

    iraw = col_sum.sort_values(ascending=False).index 
    
    #icor_map = ogtk.UM.merge_all(iraw, errors=1, mode='dist')
    icor_map2 = ogtk.UM.merge_all(iraw, errors=1, mode='dist')
    icor2 = [i[1] for i in icor_map2]
    
    #cind = col_sum[top_cols].index.isin(icor)
    #cind2 = col_sum[top_cols].index.isin(icor2)
    cind2 = col_sum.index.isin(icor2)
    
    #good= col_sum[top_cols][cind2].index
    good= col_sum[cind2].index
    good = [i for i in good if i!='null']
    good = [i for i in good if 'N' not in i]

    #pdb.set_trace()
    return(good)

def convert_df_to_mat(df: pl.DataFrame, ibar_field: str ='raw_ibar') -> pd.DataFrame:
    '''
        Returns a pandas data frame representing a matrix of umis across cell barcodes (rows) and ibars (columns)
    '''
    mat = (
        df
        .group_by(['cbc', ibar_field])
        .agg(pl.col('umi').n_unique().alias("umis"))
        .pivot(values='umis', columns=ibar_field, index='cbc')
        .drop('cbc')
        .to_pandas()
    )
    return(mat)


def guess_ibars_df(df: pl.DataFrame,
                   ibar_field: str='raw_ibar',
                   min_cells_per_ibar: int = int(1e3),
                   min_mols_per_ibar: int=1,
                   verbose: bool=False,
                   )-> Sequence:
    ''' By grouping a dataframe by cbc and ibar into a matrix of umis, a correction by rank is applied
    '''
    print(f'{min_cells_per_ibar=}')
    print(f'{min_mols_per_ibar=}')

    mat = convert_df_to_mat(df, ibar_field)
    print(f'{mat.shape=}')

    cells_per_ibar = np.nansum((mat.fillna(0)>0).to_numpy(), axis = 0)
    cmask_ibars = cells_per_ibar > min_cells_per_ibar
    
    mols_per_ibar = np.nansum(mat.values, axis = 0)
    mmask_ibars = mols_per_ibar > min_mols_per_ibar
    mask_ibars = np.logical_and(mmask_ibars,cmask_ibars)

    mat = mat.iloc[:, mask_ibars]
    verbose = True
    if verbose:
        import rich
        rich.print(f'{min_cells_per_ibar=}')
        rich.print(f'{len(cells_per_ibar)=}')
        rich.print(f'{cmask_ibars.sum()=}')
        rich.print(f'{mols_per_ibar.sum()=}')
        rich.print(f'{mmask_ibars.sum()=}')
        rich.print(f'{mask_ibars.sum()=}')
        rich.print(f'{mat.shape=}')


    valid_ibars= guess_ibars_mx(mat)
    return(valid_ibars)

def return_compiled_corrector(): 
    ''' returns a regex object capable of capturing
        the spacer region of a read while accounting for TSS variability. It assumes an
        anchor sequence (GGGTTAGA) at the beginning of the scaffold.  
        '(G{1,2}T)(?P<spacer>.+)(?P<rest>GGGTTAGA)'
    ''' 
    import regex
    spacer_corrector = regex.compile('(G{1,2}T)(?P<spacer>.+)(?P<rest>GGGTTAGA)')
    return(spacer_corrector)

def return_corrected_spacer(x, spacer_corrector, correction_pad='GGT'):
    match = spacer_corrector.search(x)
    return(f"{correction_pad}{match.group('spacer')}{match.group('rest')}"  if match else x)

def load_wl(as_pl: bool=False, turn_to_cat: bool= False):
    ''' load metadata table
    '''

    wl = pd.read_csv('/local/users/polivar/src/artnilet/conf/protospacers_singlecell.csv')
    wl.index = wl['kalhor_id']
    #cat_type = CategoricalDtype(categories=["b", "c", "d"], ordered=True)
    if turn_to_cat:
        wl['speed']= wl.speed.astype(pd.CategoricalDtype(categories=["fast", "mid", "slow"], ordered=True))
        wl['diversity'] = wl['diversity'].astype(pd.CategoricalDtype(ordered=True))
    if as_pl:
        return(pl.DataFrame(wl))
    return(wl)

def load_mol_ibarspl(sample_id, min_reads_umi=2, min_dom_cov=2):
    wl = load_wl()
    print(f'loading {sample_id}')
    mol_ibars = (
            pl.scan_csv('/local/users/polivar/src/artnilet/workdir/scv2/ibar_all_filtered.csv', separator='\t')
            .filter(pl.col('ibar')!='no_ibar')
            .filter(pl.col('sample_id')==sample_id)
            .filter(pl.col('umi_reads')>=min_reads_umi)
            .filter(pl.col('umi_dom_reads')>=min_dom_cov)
            .with_columns(pl.col('seq').str.extract("(G{1,2}T)(.+)(GGGTTAGA.+)", 2).str.replace("^", "GGT").alias('cseq'))
            .with_columns(pl.col("cseq").str.replace("GGGTTAGA", "").is_in(wl.spacer.to_list()).alias("wt"))
            .collect()
            )
    # calleles = corrected alleles
    alleles = (
        mol_ibars
        .group_by(['sample_id', 'cell', 'ibar'])
        .agg(
            [
                pl.col('cseq').n_unique().alias('n_calleles'),
                pl.col('seq').n_unique().alias('n_alleles'), 
            ])
        )

    mol_ibars = (
                mol_ibars.join(alleles, 
                            left_on=['sample_id', 'cell', 'ibar'], 
                            right_on=['sample_id', 'cell', 'ibar'])
                )
    return(mol_ibars)

def annotate_ibars(sample_id='h1e1', min_reads_umi=2, min_dom_cov=2):
    ''' returns a dictionary that maps ibars to a kalhor_id
    '''
    wl = pl.DataFrame(load_wl())
    spacer_id = dict([(i,ii) for i,ii in zip(wl.spacer, wl.kalhor_id)])
    spacer_speed = dict([(i,ii) for i,ii in zip(wl.spacer, wl.speed)])

    uncut = load_mol_ibarspl(sample_id, min_reads_umi, min_dom_cov)
    suncut = uncut.with_columns([
           pl.col('cseq').is_in(wl.spacer.to_list()).alias('wt'), 
           pl.col('cseq').map_groups(lambda x: spacer_id.get(x, None)).alias('kalhor'),
           pl.col('cseq').map_groups(lambda x: spacer_speed.get(x, None)).alias('speed')
        ])
    uncut = uncut.join(wl, left_on='cseq', right_on='spacer', how='left')

    ss = (
        uncut
        .select(['ibar', 'kalhor_id', 'cseq'])
        .filter(pl.col('kalhor_id').is_not_null())
        .group_by('ibar')
        .agg(pl.col('cseq').value_counts(sort=True).head(2))
        #.explode('cseq').unnest('cseq').rename({'':'cseq'}) # remove rename
        .explode('cseq').unnest('cseq') # remove rename
        .sort(['ibar','count'], reverse=True)
        .join(wl, left_on='cseq', right_on='spacer', how='left')
    )

    print([format(i, '.2f') for i in np.quantile(ss.counts.to_list(), np.arange(0,1, 0.1,))])
    break_point = np.quantile(ss.counts.to_list(), 0.55)
    ss = ss.filter(pl.col('count')>break_point)
    return(dict(zip(*ss.select(['ibar', 'kalhor_id']))))

def compute_ibar_table_from_tabix_zombie(tbx_ifn, valid_cells):
    valid_cells = pl.Series(valid_cells).str.replace("-1.*", "").to_list()
    u6='AGGACGAAACACC'
    anch1='GCTAGA'
    anch2='AATA'
    can_spacer="GGGTTA"

    wl = ogtk.shltr.ibars.load_wl(True)

    fuzzy_u6 =   ogtk.shltr.ibars.fuzzy_match_str(u6)
    fuzzy_anch1 = ogtk.shltr.ibars.fuzzy_match_str(anch1)
    fuzzy_anch2 = ogtk.shltr.ibars.fuzzy_match_str(anch2)
    fuzzy_canspa= ogtk.shltr.ibars.fuzzy_match_str(can_spacer)

    anchp = f'.+({fuzzy_anch1})(.{"{6}"})({fuzzy_anch2})'
    u6p = f'.+({fuzzy_u6})(.+)({fuzzy_anch2})'

    df= pl.read_csv(tbx_ifn, sep='\t', has_header=False)
    df.columns=['readid',  'start' ,'end'  , 'cbc' , 'umi' , 'seq' , 'qual']

    df = df.filter((pl.col('cbc')).is_in(valid_cells))

    dff = (
            df
            .with_columns([
                            pl.col('seq').str.extract(u6p, 2).alias('seq'),
                            pl.col('seq').str.extract(anchp, 2).alias('ibar'),
                            ])
            
            .group_by(['cbc', 'umi', 'ibar'], maintain_order=False)
                .agg([pl.col('seq').value_counts(sort=True).head(1), pl.col('seq').count().alias('umi_reads')])
                #.explode('seq').unnest('seq').rename({'':'seq', 'counts':'umi_dom_reads'}) # remove rename
                .explode('seq').unnest('seq').rename({'counts':'umi_dom_reads'}) # remove rename
    )

    valid_ibars = guess_ibars_df(dff)
    data = (
        dff
        .with_columns(pl.when(pl.col('ibar').is_in(valid_ibars)).then(pl.col('ibar')).otherwise('no_ibar'))
        .with_columns(pl.col('seq').str.extract(f"(G+T)(.+?)({fuzzy_canspa})", 2).str.replace('^',"GGT").alias('can_spacer'))
        .with_columns(pl.col('can_spacer').is_in(wl.spacer.to_list()).alias('wt'))
    )

    al = (
        data
        .group_by(['cbc', 'ibar'])
        .agg(
            [
                pl.col('seq').n_unique().alias('n_calleles'),
            ])
        )

    cfs =['cbc', 'ibar'] 
    data =  data.join(al, left_on=cfs, right_on=cfs, how='left')# <<<  annotate corrected alleles 
    return(data)


def get_ibar_table_from_tabix(tbx_ifn, valid_cells, out_fn, name, force=False, slow=False):
    import os
    if not os.path.exists(out_fn) or force:
        df = compute_ibar_table_from_tabix(tbx_ifn, valid_cells, name, slow=slow)
        df.write_parquet(out_fn)
    else:
        print(f'loading pre-computed {out_fn}')
        df = pl.read_parquet(out_fn)
    return(df)

def annotate_ibars_db(df: pl.DataFrame):
    '''Annotates data frame using a pre-determined data base
    '''
    wl = load_wl(True)
    df = df.join(wl, left_on='spacer', right_on='spacer', how='left') # <<< kalhor annotate
    return(df)

def export_ibar_mols_to_matlin(rin, sample_id ='h1e11'):
    ''' ``rin`` should the output of a cured ibar table e.g. ``compute_ibar_table()`` 
    '''
    out_fn = f'/local/users/polivar/src/artnilet/workdir/scv2/imatlin.{sample_id}'
    print(out_fn) 
    #tt = (rin.filter((pl.col("sample_id")==sample_id) & (pl.col('n_calleles')==1))
    tt = (rin.filter((pl.col('n_calleles')==1)& (pl.col('batch')==sample_id))
        .with_columns(pl.when(~pl.col('wt'))
                    .then(pl.col('seq'))
                    .otherwise("..").alias('seq'))
        .sort(['nspeed','kalhor_id'])
        )
    tt.write_csv(out_fn, has_header=True)
    return(tt)

def fuzzy_match_str(string, wildcard=".{0,1}", include_original = True, sep='|'):
    ''' returns a pattern string to be used as a fuzzy match for ``string``
    '''
    fuzz = []
    if include_original:
        fuzz.append(string)
    for ii in range(len(string)):
        s2=''
        for i,c in enumerate(string):
            if i==ii:
                s2 += wildcard
            else:
                s2 += c
        fuzz.append(s2)
    fuzz = sep.join(fuzz)
    return(fuzz)

def plot_imols_matrix(imols):
    import numpy as np
    import matplotlib.pyplot as plt
    mat= convert_df_to_mat(imols)
    cells_per_ibar = np.nansum((mat.fillna(0)>0).to_numpy(), axis = 0)
    ibars_per_cell = np.nansum((mat.fillna(0)>0).to_numpy(), axis = 1)
    ibars_per_cell = np.argsort(ibars_per_cell)[::-1]
    icells_per_ibar = np.argsort(cells_per_ibar, )[::-1]
    i100 = np.nansum((mat.fillna(0)>0).to_numpy()[:, icells_per_ibar], axis=0)
    
    plt.pcolormesh(mat.to_numpy()[:,icells_per_ibar][:, i100>1][ibars_per_cell,:])
    plt.plot(i100[i100>1])


def return_pibar(imols, ext_pattern = '(.{3})(.{4})(.+)', position=2, min_molecules=1e4):
    ''' a pibar is generated by pre-pending an arbitrary, but reasonably stable, string from the spacer.
    The idea is to provide additional specificity for filtering purposes.
    The ext pattern is a regular expression for which the pre-pended string will be generated. Position 2 is the one extracted
    e.g '(.{3})(.{4})(.+)' will capture the 4nt after the first 3.
    the final integration barcode would be ibar.ext
    '''    
    data = (
        imols
        .filter(pl.col('wt'))    # and uncut
        .with_columns(pl.col('raw_ibar')+ "." + pl.col('can_spacer').str.extract(ext_pattern, position)) # perform the actual encoding of ibar+ext
        .group_by(['raw_ibar', 'can_spacer'])
            .agg(pl.col('can_spacer').count().alias('molecules'))
        .sort('molecules', reverse=True)
    )
    valid = data.filter(pl.col('molecules')>min_molecules)['raw_ibar'].to_list()
    #sns.displot(data['molecules'].to_pandas(), kind='ecdf', aspect=3, log_scale=10)
    #sns.displot(data['molecules'].to_pandas(), kind='kde', aspect=3, log_scale=10)
    hvalids = (
            data.filter(pl.col('molecules')>min_molecules)
            .with_columns(pl.col('raw_ibar')+ "." + pl.col('can_spacer').str.extract(ext_pattern, position))
            .select('raw_ibar').to_series()
                .map_groups(lambda x: ogtk.UM.compare_umi_to_pool((x, (0, valid)))).to_list()
    )
    hvalids = np.array(hvalids)
    #sns.clustermap(hvalids, method='ward', figsize=(10,10), row_cluster=False, col_cluster=False, vmax=3, )

def stutter():
    x = fuzzy_match_str("GTGGGGTTAGA", ".") # scaffold stutter
    x = fuzzy_match_str("CACCTTCA", ".{0,1}")#36
    x = fuzzy_match_str("TTGGGTAC", ".{0,1}") # no 33 -> there is a homologous stretch  
    x = fuzzy_match_str("GCGAGACGTG", ".{0,1}") # 39 -> similar to the scaffold mm=1
    x = fuzzy_match_str("TAAATTGC", ".{0,1}") # 55 also hits in the scaffold mm=1; twice since it's on the hairpin
    x = fuzzy_match_str("CTTAACACT", ".{0,1}") # #kalhor 6 
    x = fuzzy_match_str("AAGCTGCG", ".{0,1}") # #kalhor#22 
    x = fuzzy_match_str("GTGGGGTTAGA", ".")
    px = f'.*?({x}).*?({x})'



   

def pl_tabix_to_df(tbx_ifn, do_sampling = False):
    ''' load tabix fastq as polars df
    '''
    df=pl.read_csv(tbx_ifn, sep='\t', has_header=False)
    df.columns=['readid',  'start' ,'end'  , 'cbc' , 'umi' , 'seq' , 'qual']
    return(df)

def pl_fastq_to_df(fastq_ifn, rc=False, end = None, export = False):
    ''' unfinished since the real need is about paired fastqs
    '''
    step = 4
    if end is not None:
        end = int(end*4)
    df = pl.read_csv(fastq_ifn, has_header=False, n_rows = end, sep='\n')
    fastq_fields = ['rid', 'seq', 'plus', 'qual']
    df = (
            df.with_columns(
            (pl.arange(0, pl.count()) // step).alias("step")
            ).group_by("step", maintain_order=True)
        .agg([
        pl.col("column_1").take(i).alias(name) for i, name in enumerate(fastq_fields)
        ])
        .drop(['step', 'plus'])
        )
    df = df.lazy().with_columns(pl.col('rid').str.replace("^@", "")).collect()
    if rc:
        df = (df
                .lazy()
                .with_columns([
                    pl.col('seq')
                        .str.split('')
                        .arr.reverse()
                        .arr.join('')
                        .str.replace_all('A', 't')
                        .str.replace_all('C', 'g')
                        .str.replace_all('G', 'c')
                        .str.replace_all('T', 'a')
                        .str.to_uppercase(),
                    pl.col('qual')
                        .str.split('')
                        .arr.reverse()
                        .arr.join('')
                    ])
                .collect()
            )
    if export:
        parquet_fn = fastq_ifn+'.parquet'
        print(parquet_fn)
        df.write_parquet(parquet_ifn)
    else:
        return(df)

def pl_10xfastq_to_df(fastq_ifn, end = None, export = False):
    ''' 
    DON'T USE THIS - IS IS INSANELY BAD IN TERMS OF MEMORY

    Converts a pair of 10x-derived fastq files into a polars data frame
    unfinished since the real need is about paired fastqs
    
    This version keeps read1 as is while reverse complementing read2
    
    --
    Performance notes:
    The difference in performance among concat vs join is negligible at 1e6 reads
     28.7 s for 'concat'
     28.6 s for 'join'
     since join guarantees that the reads are in the right order it is preferable
    --
    '''
    read1 = pl_fastq_to_df(fastq_ifn, end=end, export = False)
    read2 = pl_fastq_to_df(fastq_ifn.replace('_R1_', '_R2_'), rc=True, end=end, export = False)

    #read1 = read1.join(read2.with_columns(pl.col('rid').str.replace(" 2:N:", " 1:N:")), left_on='rid', right_on='rid')
    read1 = pl.concat([read1, read2])
    #read1 = read1.rename({'seq':'r1_seq', "qual":"r1_qual", "seq_right":"r2_seq", "qual_right":"r2_qual"})

    if export:
        parquet_fn = fastq_ifn+'.parquet'
        print(parquet_fn)
        read1.write_parquet(parquet_ifn)
    else:
        return(read1)

def return_feature_string(key: str):
    ''' For a provided key return the default string  
    Supported features: 

     - 'u6mx'
     - 'tso'
     - 'u6bk'
     - 'usibar'
     - 'dsibar'
     - 'canscaf'
     - 'primer'
     - 'scaf2'
    ''' 
    features_dict = {
    'u6mx':   'GACGAAACACC', # 'GACGAAACACC'
    'tso':    'TTTCTTATATGGG',
    'u6bk':   'GGAAAGAAACACCG', # broken u6
    'usibar': 'GCTAGA', # upstream ibar
    'dsibar': 'AATAGCAA', # downstream ibar
    'canscaf':'GGGTTAGAG',
    'primer': 'GGCTAGTCCGTTATCAACTTG', # applies only to RNA-based capture
    'scaf2':  'GTTAACCTAA'
    }
    assert key in features_dict, f"{key} is not a valid feature. Look at {features_dict.keys()}"
    return(features_dict[key])

@ogtk.utils.log.call
def encode_reads_with_dummies(
            parquet_ifn: str| None = None,
            df: pl.DataFrame | None = None,
            zombie = False,
            do_sampling = None,
            ) -> pl.DataFrame:
    ''' Encodes raw reads based on regular expressions. It doesn't aggregate results
    '''

    wl = load_wl(True)
    wts = "|".join(wl['spacer'])

    u6mx =     return_feature_string('u6mx')
    tso =      return_feature_string('tso')
    u6bk =     return_feature_string('u6bk')
    usibar=    return_feature_string('usibar')
    dsibar=    return_feature_string('dsibar')
    canscaf =  return_feature_string('canscaf')
    primer =   return_feature_string('primer')
    scaf2 =    return_feature_string('scaf2')

    fuzzy_bu6 =    fuzzy_match_str(u6bk, '.{0,1}')
    fuzzy_u6 =     fuzzy_match_str(u6mx, '.{0,1}')
    fuzzy_usibar = fuzzy_match_str(usibar, wildcard=".{0,1}")
    fuzzy_dsibar = fuzzy_match_str(dsibar, wildcard=".{0,1}")
    fuzzy_tso =    fuzzy_match_str(tso, wildcard=".{0,1}")  # account for sequencing errors
    fuzzy_scaff2 = fuzzy_match_str(scaf2, wildcard=".{0,1}")
    fuzzy_stammering = fuzzy_match_str("GTGGGGTTA", ".") # scaffold stutter
    fuzzy_primer=  fuzzy_match_str(primer, wildcard=".{0,1}?") # ??

    if not zombie:
        fuzzy_canspa=  fuzzy_match_str(canscaf, wildcard="") # ??
    else:
        fuzzy_canspa=  fuzzy_match_str(canscaf, wildcard="") # ?? used to be a wildcard='.'

    ####
    if all(i is None for i in [df, parquet_ifn]):
        raise ValueError("you need to provide parquet_ifn or a pl")


    def pfn(path:str):
        ''' Parse file name
        '''
        return path.split('/')[-1]
            

    # TODO change name of parquet_ifn
    # TODO warning: current downsampling only `.head()`
    # depending on whether `parquet_ifn` is a list or not, read the files accordingly
    use_pyarrow=True
    
    drop = ['readid', 'qual' , 'start', 'end']

    @ogtk.utils.log.call
    def return_parquets_scan(parquet_ifns, do_sampling=False):
        '''
        '''
        if not isinstance(parquet_ifns, list):
            parquet_ifns = [parquet_ifns]

        if do_sampling:
            return (
                    pl.scan_parquet(i)
                    .head(do_sampling)
                    .drop(drop)
                    .with_columns(pl.lit(pfn(i)).alias('ifn'))
                    .collect()  
                    for i in parquet_ifns
                    )
        else:
            return (
                    pl.scan_parquet(i)
                    .drop(drop)
                    .with_columns(pl.lit(pfn(i)).alias('ifn'))
                    .collect() 
                    for i in parquet_ifns
                    )
    # TODO check whether this splitting of the function is still needed 
    # https://github.com/pola-rs/polars/issues/8338
    @ogtk.utils.log.call
    def return_parquets_read(parquet_ifns, do_sampling=False, use_pyarrow=True):
        '''
        '''
        if not isinstance(parquet_ifns, list):
            parquet_ifns = [parquet_ifns]
        if do_sampling:
            return (
                    pl.read_parquet(i, use_pyarrow=use_pyarrow)
                    .head(do_sampling)
                    .drop(drop)
                    .with_columns(pl.lit(pfn(i)).alias('ifn'))
                    for i in parquet_ifns
                    )
        else:
            return (
                    pl.read_parquet(i, use_pyarrow=use_pyarrow)
                    .drop(drop)
                    .with_columns(pl.lit(pfn(i)).alias('ifn'))
                    for i in parquet_ifns
                    )

    if df is None:
            fs = return_parquets_read(parquet_ifn, do_sampling, use_pyarrow)
            df = pl.concat(fs)

    logger.info(f'total_reads={df.shape[0]}')
    

    if "cbc" in df.columns:
        total_cells =df["cbc"].n_unique() 
        #print(f'{total_cells=}')
    else:
        total_umis =df["umi"].n_unique() 
        #print(f'{total_umis=}')


    logger.step(f'[red]fuzzy encoding[/]', extra={"markup": True})
    # It is faster to use fuzzy matching only once, so we encode first and then count
    if zombie:
        fprime_dummie = '[···U6···]'
    df = (
        df
        .lazy()
        .with_columns(pl.col('seq').alias('oseq'))
        #.with_columns(pl.col('seq').str.replace_all(f'({wts})', f'[···WT···]'))
        .with_columns(pl.col('seq').str.replace_all(f'.*?({fuzzy_tso})', '[···TSO···]'))
        .with_columns(pl.col('seq').str.replace_all(f'.*({fuzzy_u6})', '[···U6···]'))
        .with_columns(pl.col('seq').str.replace_all(f'.*({fuzzy_bu6})', '[···bU6···]'))
        .with_columns(pl.col('seq').str.extract(f'CTAGA(.{"{6}"})({fuzzy_dsibar})', 1).alias('raw_ibar'))
        .with_columns(pl.col('seq').str.replace_all(f'({fuzzy_dsibar})', '[···SCF1···]'))
        .with_columns(pl.col('seq').str.replace_all(f'({fuzzy_scaff2})', '[···SCF2···]'))
        .with_columns(pl.col('seq').str.replace_all(f'({canscaf})', '[···CNSCFL···]'))
        .with_columns(pl.col('seq').str.replace_all(f'({fuzzy_primer})', '[···LIB···]'))
        .with_columns(pl.col('seq').str.replace_all(r'\[···SCF2···\].+\[···LIB···\]', '[···SCF2···][···LIB···]'))
        #.with_columns(pl.col('seq').str.extract(f'([A-Z]{"{6}"})\[···SCF1···\]',1).alias('raw_ibar'))

        #######.with_columns(pl.col('seq').str.replace_all(f'(\]{"CTAGA"})', '][···SCF0···]'))
        .with_columns(pl.col('seq').str.replace_all(r'\[···SCF2···\]\[···LIB···\].+$', ''))
        .with_columns(pl.col('seq').str.replace_all(f'{fuzzy_stammering}', '[XXX]'))
        .with_columns(pl.col('seq').str.replace_all(r'\[···U6···\]', '[···TSO···]'))
        .with_columns(pl.col('seq').str.replace_all(r'\[···TSO···\]\[···WT···\]', '[···TSO···][···WT···]'))
        .with_columns(pl.col('seq').str.replace_all(r'\[···SCF1···\]\[···SCF2···\].*', '[···SCF1···]')) # trims ends of zombie libraries that do not show the primer sequence


        .drop([i for i in ['qual',  'readid', 'start', 'end'] if i in df.columns])
        #.filter(pl.col('seq').str.contains(r'[···SCF1···][END]'))
        .collect()
    )
    logger.step(f'[green] done[/]', extra={"markup": True})

    return(mask_wt(df))

def empirical_metadata_annotation(
        df: pl.DataFrame,
        top_n:int=1,
        drop_counts: bool=True)->pl.DataFrame:
    ''' 
    Determines the metadata annotation to a given ibar-spacer pair. This
    function should be run on samples that have not been induced. Optionally
    (``drop_counts``) show the evidence for a given ibar-spacer pair.

    Returns a data frame that can be used for annotation (join).
    '''

    wl = ogtk.shltr.ibars.load_wl(as_pl=True).drop('clone')

    # determine the top spacer per raw_ibar
    tsp_df = (
        df
        .group_by(['sample_id', 'raw_ibar', 'spacer'])
        .len()
        .filter(
            pl.col('len').rank(method='dense', descending=True).over(['sample_id', 'raw_ibar']) <= top_n
        )
        .with_columns(nlen=pl.col('len')/pl.col('len').sum().over('sample_id') * 1e4)
        .join(wl, right_on='spacer', left_on='spacer', how='inner') #pyright:ignore

    )

    if drop_counts:
        tsp_df = tsp_df.drop('len')

    if 'clone' in tsp_df.columns:
        tsp_df = tsp_df.rename({'clone':'clone_ann'})
    return(tsp_df)



@ogtk.utils.log.call
def count_essential_patterns(df: pl.DataFrame, seq_col: str = 'seq') -> pl.DataFrame:
    ''' Count appearances of expected patterns in encoded seqs, for example TSOs, U6s, cannonical scaffolds, uncut matches
    '''
    # count matches for a series of patterns
    df = (
        
        df
        .lazy()
        .with_columns([
            pl.col(seq_col).str.count_match('TSO').alias('np_tsos'),
            pl.col(seq_col).str.count_match('SCF1').alias('np_dss'),
            pl.col(seq_col).str.count_match('XXX').alias('np_sts'),
            pl.col(seq_col).str.count_match('·U6·').alias('np_u6s'),
            pl.col(seq_col).str.count_match('bU6').alias('np_bu6s'),
            pl.col(seq_col).str.count_match('CNSCFL').alias('np_can'),
            pl.col(seq_col).str.count_match('WT').alias('np_wts'),
            pl.col('spacer'),
        ])
        .collect()
     )
    return(df)

def noise_properties(df: pl.DataFrame) -> pl.DataFrame:
    ''' Extract some metrics of systematic errors e.g. low complexity or length of spacer
    '''
    df = (
        df
        .lazy()
        .with_columns(pl.col('seq').str.extract('(\[···TSO···\])(.+?)(\[···CNSCFL···\])',2).str.n_chars().alias('spacer_len'))
        .with_columns(pl.col('seq').str.extract('(\[···TSO···\])(.+?)(\[···CNSCFL···\])',2).str.count_match('TG').alias('tg_cmp'))
        .with_columns((pl.col('tg_cmp')*2/pl.col('spacer_len')).alias('tg_cmp'))
        .with_columns(pl.when(pl.col('spacer').is_null()).then(pl.col('spacer_len')).otherwise(pl.col('spacer').str.n_chars()).alias('spacer_len'))
        .collect()
    )
    return df

def plot_noise_properties(df):
    '''
    '''
    from matplotlib.colors import LogNorm

    sns.displot(df,
                x='tg_cmp',
                y='spacer_len',
                cmap='plasma',
                bins=(10,30),
                cbar=True,
                cbar_kws=dict(norm=LogNorm())
                )

@ogtk.utils.log.call
def ibar_reads_to_molecules(
        df: pl.DataFrame,
        modality: str ='single-cell',
        )-> pl.DataFrame:
    ''' Returns a data frame at the level of molecules keep count of the reads behind each UMI as `umi_reads`
        Additionally, top representative encoded sequences and original sequences are kept
    '''
    # collapse reads into molecules
    groups =['cbc', 'umi', 'raw_ibar', 'seq'] 

    if modality =='single-molecule':
        groups = groups[1:]

    df = (
        df
        .lazy()
        .filter(pl.col('raw_ibar').is_not_null())
        .with_columns(pl.count('seq').over(groups).alias('seq_counts'))
        .group_by(groups[:-1]) #type: ignore
            .agg(
                umi_reads = pl.count(),
                has_ties = (pl.col('seq_counts') == pl.col('seq_counts').max()).sum() # perhaps can be simplified with .unique_counts ?
                            != pl.col('seq_counts').max(),
                umi_dom_reads =  pl.col('seq_counts').max(),
                seq =  pl.col('seq').take(pl.col('seq_counts').arg_max()),
                oseq = pl.col('oseq').take(pl.col('seq_counts').arg_max())
            )
        .collect()
        .explode('seq')
        .explode('oseq')
    )

    logger.info(f'filtering {df["has_ties"].sum()} ties {df["has_ties"].mean():.2%}')

    return df.filter(pl.col('has_ties').is_not()).drop('has_ties')

@ogtk.utils.log.call
def noise_spectrum(batch: str,
                   df: pl.DataFrame,
                   columns: str | Sequence = 'wts',
                   index: Sequence= ['np_sts', 'np_u6s', 'np_can', 'np_dss'],
                   do_plot:bool=False,
                   out_fn: str| None= None):
    '''
    Generate the complex heatmaps that encode noise according to features
    '''
    #TODO improve this function
    df = count_essential_patterns(df)
    xxx = df.pivot(index=index, 
                   columns=columns, 
                   values='umi', 
                   aggregate_function='count', 
                   sort_columns=True).sort(index)

    rename_d =dict([(str(i), f'{columns}_{i}') for i in range(xxx.shape[1]-len(index))]) 

    xxxp = (
        xxx.rename(rename_d)
        .to_pandas()
        .set_index(index)
        .apply(lambda x: 100*x/(df.shape[0]/1.0))
        .fillna(0)
        )

    logger.info(xxxp.apply(sum))
    xxx = (
        xxx.rename(rename_d)
        .to_pandas()
        .set_index(index)
        .apply(lambda x: np.log10(x))
        .fillna(0)
    )    
    

    if do_plot:
        fig, axes = plt.subplots(2,1, dpi=60, figsize=(12, 6*len(index)))
        
        sns.heatmap(xxxp, annot=True, fmt='.1f', ax=axes[0], annot_kws={"fontsize":16}, cmap="RdYlBu_r")
        axes[0].set_title(batch)
        axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
        
        sns.heatmap(xxx, annot=True, fmt='.1f', ax=axes[1], annot_kws={"fontsize":16}, cmap="RdYlBu_r")
        axes[1].set_title(batch)
        axes[1].set_yticklabels(axes[1].get_yticklabels(), rotation=0)
        
        plt.tight_layout()
        if out_fn is not None:
            fig.savefig(out_fn)
        
        plt.show()

    return pl.DataFrame(xxx), pl.DataFrame(xxxp)


@ogtk.utils.log.call
def mask_valid_ibars(
        df: pl.DataFrame,
        valid_ibars: Sequence | None=None,
        min_mols_per_ibar: int=1000,
        min_cells_per_ibar: int | None =1000,
        fraction_cells_ibar: float= 0.2
        )->pl.DataFrame:
    ''' 
    '''
    if min_cells_per_ibar is None and 'valid_cells' not in df.columns:
            raise  ValueError('in order to estimate a fraction of cells, a "valid_cell" boolean mask must exist on the df or an estimation must be provided (min_cells_per_ibar)')
    # flag ibars 
    if valid_ibars is None:
        # determine valid ibars automatically
        logger.info('determining number of valid ibars')
        if min_cells_per_ibar is None and "valid_cell" in df.columns:
            total_cells =df.filter(pl.col('valid_cell'))["cbc"].n_unique() 
            logger.info(f'{total_cells=}')
            min_cells_per_ibar = int(total_cells*fraction_cells_ibar)
            logger.info(f'{min_cells_per_ibar=} ({fraction_cells_ibar:%})')
        else:
            logger.info(f'{min_cells_per_ibar=} ')

        valid_ibars = guess_ibars_df(
                df, 
                min_cells_per_ibar = min_cells_per_ibar, #type: ignore
                min_mols_per_ibar=min_mols_per_ibar)
    else:
        logger.info('using provided valid ibars')
    
    logger.info(f'{len(valid_ibars)=}')

    df = df.with_columns(pl.col('raw_ibar').is_in(valid_ibars).alias('valid_ibar'))
    return(df)

@ogtk.utils.log.call 
def extract_spacer(df: pl.DataFrame, spacers:List|pl.Series) -> pl.DataFrame:
    ''' Adds 'spacer' field to data frame based on a list of spacers
    '''
    spacers_regex = '|'.join(spacers)
    return (
        df
        .with_columns(pl.col('seq').str.extract(f'({spacers_regex})', 1).alias('spacer'))
        #purify spacer sequence    
        .with_columns(pl.col('seq').str.replace(r'CTAGA.+?$', '').alias('xspacer'))
        .with_columns(pl.col('xspacer').str.replace(r'\[···TSO···\]', ''))
        .with_columns(pl.col('xspacer').str.replace(r'\[···CNSCFL···\]', ''))
        .with_columns(pl.col('xspacer').str.replace(r'\[···bU6···\]', ''))
        .with_columns(pl.when(pl.col('spacer').is_null()).then(pl.col('xspacer')).otherwise(pl.col('spacer')).alias('spacer'))
        .drop('xspacer')
        #.with_columns(pl.col('seq').str.replace_all(f'({wts})', f'[···WT···]'))
    )

@ogtk.utils.log.call
def encode_wt(df: pl.DataFrame, spacers:List) -> pl.DataFrame:
    ''' Simplifies the 'seq' field by replacing an uncut spacer for the string [···WT···]
    '''
    spacers_regex = '|'.join(spacers)
    return (
        df
        .with_columns(pl.col('seq').str.replace_all(f'({spacers_regex})', f'[···WT···]'))
    )
@ogtk.utils.log.call
def mask_wt(df: pl.DataFrame) -> pl.DataFrame:
    ''' Simple routine to add a boolean field that defines uncut entries
    '''
    wt_str = r'\[···WT···\]\[···CNSCFL···\]'

    df= (df
         .with_columns(
             ((pl.col('seq').str.count_match('WT')==1) 
                       & (pl.col('seq').str.contains(wt_str)))
        .alias('wt'))
     )
    return df

def create_rectangle_patches(
        data: pl.DataFrame,
        rects: List[patches.Rectangle],
        x_delta: float = 0.0) -> List[patches.Rectangle]:
    """
    Create rectangle patches to be plotted based on the input DataFrame.
    
    Args:
        data (pl.DataFrame): Input DataFrame containing the necessary information to create rectangles.
        rects (List[patches.Rectangle]): List of existing rectangle patches.
        x_delta (float, optional): Shift in x-coordinate for the rectangles. Defaults to 0.0.

    Returns:
        List[patches.Rectangle]: List of rectangle patches to be plotted.
    """
    for ii, i in enumerate(data.iter_rows(named=True)):
        ii = ii * x_delta
        width = 0.75 
        xy =  (ii+i['x'] - width/2, i['y'])  # adjust x to be centered on x-tick
        height = i['freq']
        facecolor = i['color'] if i['wts']<1 else '#dddddd'
        rect = patches.Rectangle(xy=xy, width=width, height=height, facecolor=facecolor)
        rects.append(rect)
    return rects

def plot_stacked_rectangles(
        rects: List[patches.Rectangle],
        title: str,
        groups: List[str]) -> None:
    """
    Plot stacked rectangles on a plot.
    
    Args:
        rects (List[patches.Rectangle]): List of rectangle patches to be plotted.
        title (str): Title of the plot.
        groups (List[str]): List of group names for the x-ticks.
    """
    fig, ax = plt.subplots(figsize=(len(groups), 5))
    for rect in rects:
        ax.add_patch(rect)
        
    ax.set_xlim((-0.5, len(groups) - 0.5))  # set x limit to ensure all rectangles are displayed
    ax.set_title(title)
    ax.set_xticks(range(len(groups))) 
    ax.set_xticklabels(groups, rotation=90)
    plt.show()

def prepare_stacked_df(
        df: pl.DataFrame,
        element: str,
        grouping_field: str,
        grouping_value: str,
        n_top: int,
        element_field: str='raw_ibar',
        x: float =0.5,
        bc='seq') -> pl.DataFrame:
    """
    Prepare a DataFrame to be plotted as stacked rectangles.
    
    Args:
        df (pl.DataFrame): Input DataFrame.
        element (str): Element value to filter on.
        grouping_field (str): Field to group the DataFrame by.
        grouping_value (str): Value to filter the DataFrame on based on the grouping field.
        n_top (int): Number of top rows to select after sorting.
        element_field (str, optional): Field to filter based on element. Defaults to 'raw_ibar'.
        x (float, optional): Value to assign to the 'x' column of the DataFrame. Defaults to 0.5.
        bc (str, optional): either 'seq' or 'spacer', represents the source of the barcode. Defaults to 'seq'.
        
    Returns:
        pl.DataFrame: DataFrame ready to be plotted.
    
    Raises:
        AssertionError: If bc is not 'seq' or 'spacer'
    """

    assert bc in ['seq', 'spacer'], "bc must be equal to 'seq' or 'spacer' as it represents the source of the barcode"
    data= (
        df
        .filter(pl.col(element_field)==element)
        .filter(pl.col(grouping_field)==grouping_value)
        .group_by([element_field, bc, 'wts']).count()
        .sort(['count', bc], descending=True)
        .head(n_top)
    )

    if data.shape[0] == 0:
        return(pl.DataFrame({'x':x, 'y':1, 'freq':1, 'wts':1 }))

    data = (
        data
        .with_columns(pl.col(bc).apply(lambda x: ColorHash(x).hex, skip_nulls=True).alias('color'))
        .with_columns((pl.col('count')/pl.col('count').sum()).alias('freq'))
        .with_columns(pl.col('freq').shift(1).fill_null(0).alias('y'))
        .with_columns(pl.col('y').cumsum())
        .with_columns(pl.lit(x).alias('x'))        
    )
    return data

def plot_indel_stacked_fractions(element: str,
        df: pl.DataFrame,
        groups: List[str],
        element_field: str='raw_ibar',
        grouping_field: str ='sample_id',
        n_top: int =10,
        bc: str='seq') -> None:
    """
    Plot indel stacked fractions based on the input DataFrame and groups.
    
    Args:
        element (str): Element value to filter on.
        df (pl.DataFrame): Input DataFrame.
        groups (List[str]): List of groups to be plotted.
        element_field (str, optional): Field to filter based on element. Defaults to 'raw_ibar'.
        grouping_field (str, optional): Field to group the DataFrame by. Defaults to 'sample_id'.
        n_top (int, optional): Number of top rows to select after sorting. Defaults to 10.
        bc (str, optional): either 'seq' or 'spacer', represents the source of the barcode. Defaults to 'seq'.
    """
    rects = []
    for idx, grouping_value in enumerate(groups):
        data = prepare_stacked_df(
               df,
               element,
               element_field=element_field,
               grouping_field=grouping_field,
               grouping_value=grouping_value,
               n_top=n_top,
               x=idx,
               bc=bc)
        rects = create_rectangle_patches(data, rects)

    plot_stacked_rectangles(rects, title=f'{element_field} {element} {bc}', groups=groups)
# TODO support passing a dictionary of ibar -> uncut spacers

@ogtk.utils.log.call
def return_aligned_alleles(
    df,
    aligner:Union["pyseq_align.NeedlemanWunsch", "pyseq_align.SmithWaterman"], #pyright: ignore
    ref_spacer:str|None=None,
    lim:int=50,
    correct_tss_coord:None|int=18,
    min_group_size:int=100,
    tso_pad:str|None="",
    )->pl.DataFrame:
    '''
    Requires a data frame that:
    - belongs to a single ibar

    E.g.
    df_alg = (
            df
            .sample(int(1e5))
            .group_by('raw_ibar')
            .map_groups(ogtk.shltr.ibars.return_aligned_alleles)
    )

    Returns the original data frame with additional fields that contain a an aligned version of the original sequence
     = position-specific code:
         - -1 is for gap
         - 0 match
         - >0 insertion length

    Alignment procedure:
     - Trim the original sequence based on fuzzy 5- and 3-prime anchor
       sequences. In the case of RNA-based detection the TSO fulfills this
       role. For DNA-based (zombie) it's the 3' end of the U6 promoter 
     - 
     # TODO watch out for stammering on alignments - how do they look? search for XXX in the seq


    '''

    #         [    TSO     ]
    ref_str = "TTTCTTATATGGG[SPACER]GGGTTAGAGCTAGA[IBAR]AATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTG"
    ref_str = "TTTCTTATATGGG[SPACER]GGGTTAGAGCTAGA[IBAR]AATAGCAAGTTAACCTAAGGCTAGTCC" #GTTATCAACTTGGTACT"
    schema = {"sample_id":str, "ibar":str, "id":str, "aseq":str, "sweight":pl.Int64, "alg":pl.List(pl.Int64)}
    oschema = {"aseq":str, "alg":pl.List(pl.Int64), "aweight":pl.Int64}
    merged_schema = df.schema.copy()
    merged_schema.update(oschema)

    query_schema = {'seq_trim':pl.Utf8, 'oseq':pl.Utf8}

    # alignment.schema=OrderedDict([('oseq', String), ('cigar', String), ('alg_ref', String), ('alg_seq', String), ('alg_score', Int64)])

    alignment_schema = {
           'cigar_str':pl.Utf8,
           'aligned_ref': pl.Utf8, 
           'aligned_seq': pl.Utf8, 
           'alignment_score': pl.Int64,
           }
    
    merged_schema = query_schema.copy()
    merged_schema.update(alignment_schema)

    assert correct_tss_coord is None or isinstance(correct_tss_coord, int), \
            "correct_tss_coords must be None or an integer that defines the correction coordinate"

    #assert df['sample_id'].n_unique() == 1, \
    #        "Please provide a data frame pre-filtered with a single sample_id. \
    #        Best when this function is called within a group_by"

    assert df['raw_ibar'].n_unique() == 1, \
            "Please provide a data frame pre-filtered with a single ibar. \
            Best when this function is called within a group_by"

    #foc_sample = df['sample_id'][0]
    foc_ibar= df['raw_ibar'][0]

    if df.shape[0]<min_group_size:
        logger.info("group too small")
        return pl.DataFrame(schema=merged_schema)

    # reduce the size of the TSO sequence to decrease the specificity of the match
    tso_string = return_feature_string("tso")
    tso_pad = tso_string if tso_pad is None else tso_pad

    fuzzy_tso = fuzzy_match_str(tso_string[5:], wildcard=".{0,1}")
    fuzzy_u6 = fuzzy_match_str(return_feature_string('u6mx')[5:], wildcard=".{0,1}")

    fuzzy_5pend = f'{fuzzy_tso}|{fuzzy_u6}'
    #                              [Common region]
    fuzzy_3pend = fuzzy_match_str("CCTAAGGCTAGTCC", wildcard=".{0,1}")

    
    # filter out sequences that can't be controlled
    # TODO keep a record
    df = (
            df
            .lazy()
            .filter(
                pl.col('oseq')
                .str.contains(f'{fuzzy_5pend}'))
            .filter(
                pl.col('oseq')
                .str.contains(f'{fuzzy_3pend}'))
            .collect()
        )
    # compute seq_trim
    df = (
            df
            .lazy()
            .with_columns(
                # a parenthesis is needed for the fuzzy pattern to include the wildcards at the beginning
                pl.col('oseq')
                .str.replace(f'^.+?({fuzzy_5pend})', tso_pad)
                .str.replace(f'{fuzzy_3pend}', '')
                .alias('seq_trim')
            )
            .collect()
        )

    ### delete?
    if len(df)==0:
        logger.debug(f"No entries found to generate a fasta file {foc_ibar} ")
        return pl.DataFrame(schema=merged_schema)

    if ref_spacer is None:
        # attempts to guess a spacer based on the most frequent wt allele
        guess_spacer = (
                df
                .filter(pl.col('wt'))['spacer']
                .value_counts()
                .sort('count', descending=True)
                )

        if guess_spacer.shape[0]==0:
            logger.debug(f'No WT allele found for {foc_ibar} ')
            return pl.DataFrame(schema=merged_schema)

        ref_spacer = guess_spacer['spacer'][0]
        guess_spacer = None

    assert isinstance(ref_spacer, str), "ref_spacer should be a string"

    foc_ref = ref_str.replace('[SPACER]', ref_spacer).replace('[IBAR]', foc_ibar) 
    foc_ref = foc_ref.replace(return_feature_string("tso"), tso_pad)

    logger.debug(f'{foc_ref=}')

    # ibars could show different sizes as spacer are varible
    # we use the downstream ibar sequence to anchor the ibar position
    foc_ref_match = regex.search(f'{foc_ibar}{return_feature_string("dsibar")}', foc_ref)
    lim = lim if foc_ref_match is None else foc_ref_match.start()

    coord_bins = range(lim)
    query = df.select('seq_trim', 'oseq').unique()

    # perform the alignment
    query =(
            query
            .with_columns(
                alignment= pl.col('seq_trim')
                   .map_elements(lambda x: lambda_needlemanw(x, foc_ref, aligner), 
                            return_dtype=pl.Struct(alignment_schema))
            )
            .select('oseq', 'alignment', 'seq_trim')
            )

    alignment = (
            query.with_columns(
                cigar_str=pl.col('alignment').struct.field('cigar_str'), 
                aligned_ref=pl.col('alignment').struct.field('aligned_ref'), 
                aligned_seq=pl.col('alignment').struct.field('aligned_seq'), 
                alignment_score=pl.col('alignment').struct.field('alignment_score'), 
            )
            .drop('alignment')
        )
    query = None

    alignment = alignment.select(list(merged_schema.keys()))
    return alignment

def compare_ibar_hdistance(seqs:Sequence, mols, select:None|Sequence=None, title='', jobs = 10):
    # how to find ibar shadows
    mm = ogtk.UM.hdist_all(seqs, jobs=jobs)
    assert (mm==mm.T).all(), "Hamming distance matrix is not symmetric"

    mm = np.triu(mm, k=0)

    mp = pl.DataFrame(mm, schema=list(seqs)).with_columns(raw_ibar=seqs)
    mp = mp.melt(id_vars='raw_ibar', variable_name="raw_ibar2", value_name="hdist")
    mp.shape

    dselect = ['raw_ibar', 'ibar_umisg']

    select = np.hstack([dselect, select]) if select is not None else dselect

    mp  = (
        mp
        .join(mols.select(select).unique(), 
              left_on='raw_ibar', 
              right_on='raw_ibar', 
              how='left')
        .join(mols.select(select).unique(), 
              left_on='raw_ibar2',
              right_on='raw_ibar', 
              how='left', 
              suffix="2")
        .filter((pl.col("hdist")>0) | (pl.col('raw_ibar')==pl.col('raw_ibar2')))
    )
    
    return mp

    #g = sns.relplot(data=mp
    #            .filter(pl.col('cluster')=="1")
    #            .sort("hdist", descending=True)
    #                .to_pandas(),
    #            x="ibar_umisg",
    #            y="ibar_umisg2",
    #            col="cluster2",
    #            row="hdist",
    #           )
    #
    #g.set(xscale='log')
    #g.set(yscale='log')

    grata = (mp.group_by("cluster", "cluster2", "hdist").agg(
        pear=pl.corr("ibar_umisg", "ibar_umisg2", method="pearson"),
        spear=pl.corr("ibar_umisg", "ibar_umisg2", method="spearman"),
        ).sort('cluster')
    )

    g = sns.catplot(data=grata.to_pandas(), y="pear", x="hdist", col="cluster", kind="box")
    g.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
    g.fig.suptitle(f'Pearson {title}')
    g.set(ylim=(-1, 1))
    plt.show()
    
    g = sns.catplot(data=grata.to_pandas(), y="spear", x="hdist", col="cluster", kind="box")
    g.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
    g.fig.suptitle(f'Spearman {title}')
    g.set(ylim=(-1, 1))
    plt.show()

@lru_cache
def lambda_needlemanw(seq, foc_ref, aligner):
    ''' 
    '''
    alignment = aligner.align(foc_ref, seq)

    return {
        'cigar_str': ngs_tools.sequence.alignment_to_cigar(alignment.result_a, alignment.result_b),
        'aligned_ref':  alignment.result_a,
        'aligned_seq':  alignment.result_b,
        'alignment_score':  alignment.score
        }



def default_aligner():
    return NeedlemanWunsch(
            gap_open=-20,
            gap_extend=-1,
            no_start_gap_penalty=True,
        )


def align_alleles(df,
                  ref_df:None|pl.DataFrame=None,
                  alignment_groups=['raw_ibar', 'sample_id'],
                  show_progress=True,
                  alg_kwargs:None|Dict[str,Union[None, Any]]=None,
                  ) -> pl.DataFrame:
    '''
    Aligns the alleles to their corresponding references in place.
    Adds the following fields:
    - 'cigar_str',
    - 'aligned_ref',
    - 'aligned_seq',
    - 'alignment_score'
    show_progress controls the execution with a progress bar present
    '''
    from rich.progress import Progress

    assert all([i in df.columns for i in ['raw_ibar', 'sample_id', 'oseq']]), \
        "df must include raw_ibar, sample_id and oseq in its fields"

    groups = df.group_by(*alignment_groups)
    
    if alg_kwargs is None:
        alg_kwargs = {'aligner':default_aligner(), 'min_group_size':1}

    assert alg_kwargs is not None, "a set of arguments needs to be passed to the alignment"

    def return_spacer(ref_df, group):
        ''' '''
        if ref_df is None:
            return None

        foc_ibar = group['raw_ibar'].unique()
        assert len(foc_ibar)==1, "single ibar groups should be passed"
        foc_ibar = foc_ibar[0]

        foc_spacer=ref_df.filter(pl.col('raw_ibar')==foc_ibar)['spacer']
        assert len(foc_spacer)>0, f"ibar {foc_ibar} is not present in ref_df" 

        foc_spacer = foc_spacer[0]
        logger.debug(f'{group.shape=}\n{foc_ibar=} {foc_spacer=}')
        return foc_spacer

    if show_progress:
        with Progress() as progress:
            total = groups.len().shape[0]
            alignment_task = progress.add_task("[red]Aligning...", total=total)

            def process_group_pro(group):
                alg_kwargs['ref_spacer'] = return_spacer(ref_df, group) #pyright:ignore
                result = return_aligned_alleles(group, **alg_kwargs) #pyright:ignore
                progress.update(alignment_task, advance=1)
                return result

            df_alg = groups.map_groups(process_group_pro)
    else:
        def process_group(group):
            alg_kwargs['ref_spacer'] = return_spacer(ref_df, group) #pyright:ignore
            return return_aligned_alleles(group, **alg_kwargs) #pyright:ignore

        df_alg = groups.map_groups(process_group)

    return df.join(other=df_alg.unique(), left_on='oseq', right_on='oseq', how='left')


def parse_cigar(df, block_dels:bool=False):
    import rogtk as rr
    # import polars_hash as plh
    return (
            df
            .with_columns(
                cigop=rr.parse_cigar('cigar_str', block_dels) #pyright: ignore
            )
            .with_columns(pl.col('cigop').str.split('|')).explode('cigop')
            .with_columns(pl.col('cigop').str.split(','))
            .with_columns(
                operation=pl.col('cigop').list.get(0),
                coord=pl.col('cigop').list.get(1).cast(pl.Int32),
                oper_len=pl.col('cigop').list.get(2).cast(pl.Int32),
            )
            .drop('cigop')
        )

def explode_alignment_by_coord(df, 
       max_coord=40,
       statistic = 'cigar_str', 
       grouping_fields= [
               'operation',
               'coord',
               'raw_ibar',
               "sample_id",
               "speed",
               "kalhor_id",
               'clone'],
    ):
    '''
    '''
    gf = set(grouping_fields)
    gf.add(statistic)

    assert 'cigar_str' in df.columns, \
            "An aligned data frame must be provided. Expecting `cigar_str` field."

    mm = parse_cigar(df).filter(pl.col('coord')<=max_coord)

    mmm = len_all_combos(
            df=mm.select(list(gf)),
            key=list(gf),
            fill_value=0
    )

    gf.remove(statistic)
    mmmm= (
            mmm.group_by(*list(gf))
            .agg(sum=pl.col('len').sum(), mean=pl.col('len').mean()*1e4)
            )
    return mmmm

def all_keys(df_keys):
    '''
    Create data frame of all unique combinations of values. \n
    Courtesy of mcrumiller.
    '''
    columns = df_keys.columns
    df = df_keys.select(columns[0]).unique()
    for column in columns[1:]:
        df = df.join(df_keys.select(column).unique(), how="cross")
    return df


def len_all_combos(df, key, fill_value=None):
    '''
    Perform group-by with op and include all possible key combinations. \n
    Courtesy of mcrumiller.
    '''
    keys = all_keys(df.select(key).unique())
    #df_grouped = df.group_by(key).len()
    df_grouped = df.select(key).unique().group_by(key).len()


    df_out = keys.join(df_grouped, on=key, how="left")

    if fill_value is not None:
        df_out = df_out.with_columns(pl.col('len').fill_null(fill_value))

    return df_out
