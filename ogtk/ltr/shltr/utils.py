import pdb
from typing import Optional
import itertools
import regex
import gzip
import pandas as pd
import pdb
import IPython
import numpy as np
import matplotlib.pyplot as plt
import ogtk
import pysam
from importlib import reload

from .ibars import *

def return_bg_model(ifn, q = 0.5, end = int(1e4)):
    bg = []
    with zipped(ifn) as fastq, zipped(ifn) as fastq2:
        if end == None:
            end = None
        else:
            end = int(end) *4
        rit = itertools.islice(fastq, 1, end, 4)
        qit = itertools.islice(fastq2, 3, end, 4)
        it = zip(rit, qit)
        for read, qual in it:
            bg.append([ord(i)-33 for i in qual.strip()])

    bgg = np.array(bg)
    return(np.quantile(bgg, q = q, axis=0))


def zipped(ifn, mode = 'rt'):
    ''' returns the right handle for a file depending on
        whether it is gzipped or not
    '''
    if ifn.endswith('.gz'):
        print(f'Opening compressed file {ifn}')
        return(gzip.open(ifn, mode))
    else:
        print(f'Opening txt file {ifn}')
        return(open(ifn, mode))

def annotate_process_fastq(ifn, wl, wd, errors = 0, sample_end = int(5e6), model_end = int(1e4), model_quant = 0.5, ntop=100, min_reads_per_ibar=1000, sample_prob = 2/3):
    ''' Randomized exploration of reads with a sample_end depth aiming to recover up to reads_per_ibar reads for the ntop ibars. It returns a cloud of read ids beloning to a the topn ibars in order to annotate chains 
    '''
    # there are some concerns around this approach as ntop is somehow arbitrary, we might be adding too many false positives
    # with the same weight as true ones. The disjoint nature of the analysis (4c pipeline and this script on the other).
    cs1 = "GCTTTAAGGCCGGTCCTAGCAA"
    cs2 = "GCTCACCTATTAGCGGCTAAGG"
    css = {cs1:"cs1", cs2:"cs2"}
    anchor_pam = "GGGTTAGAGCTAGA"
    anchor_u6 = "AGGACGAAACACC"
    anchor_pad = "AATAGCAAGTTAA"
    anatomy = es(anchor_u6, errors)+"(?P<spacer>.+)"+es(anchor_pam, errors)+"(?P<ibar>.{6})"+es(anchor_pad, errors)+"(.+)"
    reanatomy = regex.compile(anatomy)
    #return(reanatomy)
    step = 10000

    qc_stats = [(0,0)]
    
    bgg = return_bg_model(ifn, q = model_quant, end = model_end)
    #print(bgg)

    limit_cov = ntop*min_reads_per_ibar
    limit_hits = limit_cov*3
    limit_delta = step*0.1
    
    with zipped(ifn) as fastq, zipped(ifn) as fastq1, zipped(ifn) as fastq2 :
        ibars = {cs1:[], cs2:[]}
        spacers = {cs1:[], cs2:[]}
        entries = {}    
        ii = []
        cases = nowl= hits = nohits = nomatch = noqual = no_cs = 0
        
        if sample_end == None:
            sample_end = None
        else:
            sample_end = sample_end * 4
            
        idit = itertools.islice(fastq, 0, sample_end, 4)
        rit = itertools.islice(fastq1, 1, sample_end, 4)
        qit = itertools.islice(fastq2, 3, sample_end, 4)
        it = zip(idit, rit, qit)
        
        debug = 0
        for rid, read, qual in it:
            if 'N' not in read:
                match = reanatomy.search(read)
                if match and np.random.uniform() > sample_prob:
                    cases+=1
                    spacer = match.group('spacer')
                    ibar = match.group('ibar')
                    quals = [ord(i)-33 for i in qual.strip()]
                    if ibar and spacer:
                        if spacer in wl['pspacer'].to_list():
                            is0, is1 = match.span('ibar')
                            ib_rq = quals[is0:is1]
                            ib_mq = bgg[is0:is1]
                            succ_ib = [ r>=m for r,m in zip(ib_rq, ib_mq) ]
                            
                            sp0, sp1 = match.span('spacer')
                            sp_rq = quals[sp0:sp1]
                            sp_mq = bgg[sp0:sp1]
                            succ_sp = [ r>=m for r,m in zip(sp_rq, sp_mq) ]
                            
                            if debug and False:
                                pdb.set_trace()
                            if debug <=3:
                                print(f'{np.mean(succ_ib):.2f}', f'{np.mean(succ_sp):.2f}')
                                debug+=1
                            
                            if np.mean(succ_ib) >= 0.009 and np.mean(succ_sp) >= 0.005:
                                #ibars[cs].append(ibar)
                                #spacers[cs].append(spacer)
                                hits+=1
                                skey = f'{ibar}_{spacer}'
 
                                if skey not in entries.keys():
                                    entries[skey] = []
                                if len(entries[skey]) <= min_reads_per_ibar:
                                        rid = rid.replace("@","").split(" ")[0]
                                        entries[skey].append(rid)
                                if hits%step ==0:
                                #if hits%int(sample_end*0.1) ==0:
                                    top_cov = np.sum(sorted([len(i) for i in entries.values()][0:ntop]))
                                    delta =top_cov-qc_stats[-1][0]
                                    qc_x = (cases, hits, top_cov, limit_cov, hits, limit_hits)
                                    #qc_stats.append(qc_x)
                                    qc_stats.append((top_cov, delta, hits))
                                    if top_cov > limit_cov or hits > limit_hits or delta < limit_delta:
                                        break
                            else:
                                noqual+=1
                        else:
                            nowl+=1
                    else:
                        no_cs+=1
                else:
                    nomatch+=1
            else:
                nohits+=1
    print(f'hits\t{hits}\t{4*hits/sample_end:.2f}\n'
          f'nowl\t{nowl}\t{4*nowl/sample_end:.2f}\n'
          f'nohits\t{nohits}\t{4*nohits/sample_end:.2f}\n'
          f'bad cs\t{no_cs}\t{4*no_cs/sample_end:.2f}\n'
          f'noqual\t{noqual}\t{4*noqual/sample_end:.2f}\n'
        f'nomatch\t{nomatch}\t{4*nomatch/sample_end:.2f}')

    seq_to_id = dict([(i,"id"+str(ii)) for i,ii in zip(wl['pspacer'], wl['id'])])
   
    df = pd.DataFrame([(k,len(v)) for k,v in entries.items()])
    df.columns = ['ibar', 'count']
    df = df.sort_values(['count', 'ibar'], ascending=[False, True])
    df = df[df['count'] >= np.max(df['count']) * 0.95]
    
    def write_ibar_reads(ibar_reads_fn, reads):
        return(pd.Series(reads).to_csv(ibar_reads_fn, index=False, header=False)
               ) 

    def parse_line(i):
        f_ = i.split('\t')
        out = []
        for i in f_[2:]:
            i = i.split(':')
            i = f'{i[1]}\t{i[2]}'
            out.append(f'{f_[0]}\t{i}')
        return(out)

    df_out = []
    for ibar in df['ibar'].to_list():
        ibar_reads_fn = f'{wd}/ibar_{ibar}.txt'
        write_ibar_reads(ibar_reads_fn, entries[ibar])
        stdout = call_grep_reads(wd = wd, reads_fn = ibar_reads_fn)

        dff = [parse_line(i) for i in stdout]
        dff = pd.DataFrame([i.split('\t') for i in np.hstack(dff) if 'chrLin' not in i])
        dff.columns = ['rid', 'chrom', 'start']
        dff.start = dff.start.map(int)

        per = dff['chrom'].value_counts(normalize = True)
        top = int(np.round(1/per[0]))
        which_max = per.head(top).index.map(str)
        per = per.head(top)
        dff = dff[dff.chrom.apply(lambda x: x in which_max)]
        dff['ibar'] = ibar
        dff = dff.groupby(['ibar', 'chrom'])['start'].apply(lambda x: \
                                        np.median(x.drop_duplicates())).reset_index()
        df_out.append(dff.loc[:,['ibar','chrom','start']])#(ibar, which_max, np.median(dff['start'])))

    df_out = pd.concat(df_out)
    #df_out.columns = ['ibar', 'chrom', 'start']
    df_out.start = df_out.start.map(int)
    df_out = df_out.sort_values(['chrom', 'start'])
    df_out['spacer'] = df_out.ibar.apply(lambda  x: x.split('_')[1])
    df_out['ibar'] = df_out.ibar.apply(lambda  x: x.split('_')[0])
    df_out.to_csv(f'{wd}/ibar_int_coords.txt', sep='\t', index=False, header = True)


def call_grep_reads(wd, reads_fn):
    import subprocess
    cmd_grep = f'grep -F -f {reads_fn} {wd}/chain'
    cmd1 = subprocess.Popen(cmd_grep.split(), stdout=subprocess.PIPE, text=True)
    cmd1_stdout, cmd1_stderr = cmd1.communicate()
    return(cmd1_stdout.split('\n'))


def genotype_process_fastq(ifn, wl, errors = 0, sample_end = int(1e4), model_end = int(1e4), model_quant = 0.5):
    cs1 = "GCTTTAAGGCCGGTCCTAGCAA"
    cs2 = "GCTCACCTATTAGCGGCTAAGG"
    css = {cs1:"cs1", cs2:"cs2"}
    anchor1 = "GGGTTAGAGCTAGA"
    anchoru6 = "AGGACGAAACACC"
    anchorovlp = "TTATCAACTTGA" #delimits cs 5'
    anchor2 = "AAAAGTGGCACCG"
    anatomy = es(anchoru6, errors)+"(?P<spacer>.+)"+es(anchor1, errors)+"(?P<ibar>.{6}).+"+es(anchorovlp, errors)+"(?P<cs>.+)"+es(anchor2, errors)
    reanatomy = regex.compile(anatomy)
    
    bgg = return_bg_model(ifn, q = model_quant, end = model_end)
    #print(bgg)
    
    with zipped(ifn) as fastq:
        ibars = {cs1:[], cs2:[]}
        spacers = {cs1:[], cs2:[]}
        entries = []    
        ii = []
        nohits = 0
        noqual =0
        no_cs = 0
        
        if sample_end == None:
            sample_end = None
        else:
            sample_end = sample_end * 4
            
        rit = itertools.islice(fastq, 1, sample_end, 4)
        qit = itertools.islice(fastq, 3, sample_end, 4)
        it = zip(rit, qit)
        
        debug = 0
        for read, qual in it:
            if 'N' not in read:
                match = reanatomy.search(read)
                if match:
                    spacer = match.group('spacer')
                    ibar = match.group('ibar')
                    cs = css.get(match.group('cs'), False)
                    quals = [ord(i)-33 for i in qual.strip()]
                    if cs:
                        is0, is1 = match.span('ibar')
                        ib_rq = quals[is0:is1]
                        ib_mq = bgg[is0:is1]
                        succ_ib = [ r>=m for r,m in zip(ib_rq, ib_mq) ]
                        
                        sp0, sp1 = match.span('spacer')
                        sp_rq = quals[sp0:sp1]
                        sp_mq = bgg[sp0:sp1]
                        succ_sp = [ r>=m for r,m in zip(sp_rq, sp_mq) ]
                        
                        if debug and False:
                            pdb.set_trace()
                        if debug <=3:
                            print(f'{np.mean(succ_ib):.2f}', f'{np.mean(succ_sp):.2f}')
                            debug+=1
                        
                        if np.mean(succ_ib) >= 0.009 and np.mean(succ_sp) >= 0.005:
                            #ibars[cs].append(ibar)
                            #spacers[cs].append(spacer)
                            entries.append([cs,spacer,ibar])
                        else:
                            noqual+=1
                    else:
                        no_cs+=1
                else:
                    nohits+=1
            else:
                nohits+=1
    print(f'nohits\t{nohits}\nbad cs\t{no_cs}\nnoqual\t{noqual}')
    #print(f'nohits\t{nohits/end}\nbad cs\t{lalal/end}\nnoqual\t{noqual/end}')
    #ibars[cs1] = pd.Series(ibars[cs1])
    #ibars[cs2] = pd.Series(ibars[cs2])
    #spacers[cs1] = pd.Series(spacers[cs1])
    #spacers[cs2] = pd.Series(spacers[cs2])
    #plt.close()
    #plt.plot(spacers[cs1].value_counts().to_list())
    #plt.plot(spacers[cs2].value_counts().to_list())
    seq_to_id = dict([(i,"id"+str(ii)) for i,ii in zip(wl['pspacer'], wl['id'])])
    
    df = pd.DataFrame(entries, columns=["cs_seq", "spacer_seq", "ibar"])
    
    df['ibar'] = df['ibar'].astype(pd.CategoricalDtype(categories = set(df['ibar']), ordered=False))
    df['spacer'] = df['spacer_seq'].apply(lambda x: seq_to_id.get(x, np.nan))
    df['spacer_seq'] = df['spacer_seq'].astype(pd.CategoricalDtype(categories = wl['pspacer']))
    return(df)

def genotype_process_entries(df, top_n = 10):
    #df.set_index(["cs","spacer"], inplace=True)
    aa = df.groupby(['cs_seq', 'spacer'])
    mm = aa.apply(lambda x: x['ibar'].value_counts(normalize=False).head(top_n)).unstack(fill_value=0)
    mask = np.nansum(mm, axis=1) > 20
    mask_ibar = np.nansum(mm, axis=0) >2
    fig, ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle(f'Top {top_n} iBARs per proto-spacer', fontsize=16)
    
    #ax[0].matshow(mm, aspect= 'auto')
    ax.matshow(mm.iloc[mask,:].iloc[:, mask_ibar], aspect= 'auto')
    idx = [ "\n".join(i) for i in mm.iloc[mask,:].iloc[:, mask_ibar].index.to_list()]
    
    ax.set_yticks(range(len(idx)))
    ax.set_yticklabels(idx)
    #[i for i in ax[0].get_yticks()])
    #mm.iloc[mask,:].iloc[:, mask_ibar].plot(kind='bar', figsize=(10,2)).legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
    
    

def zuang(x):
    '''guess number of cases based on fraction of top 
    '''
    per = x['ibar'].value_counts(normalize = True)
    top = int(np.round(1/per[0]))
    per = per.head(top)
    return(per)

#https://docs.python.org/3/library/itertools.html#itertools.zip_longest
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

    
def transfer_intids_fastqs(fq1):
    ''' Transfer the integration information from R1 to the id of R2. 
    '''
    import os
    import itertools
    errors = 0
    anchor_pam = "GGGTTAGAGCTAGA"
    anchor_u6 = "AGGACGAAACACC"
    anchor_pad = "AATAGCAAGTTAA"
    anatomy = es(anchor_u6, errors)+"(?P<spacer>.+)"+es(anchor_pam, errors)+"(?P<ibar>.{6})"+es(anchor_pad, errors)+"(.+)"
    reanatomy = regex.compile(anatomy)

    fq1s = fq1.split('/')
    dout = '/'.join(fq1s[0:-1])+'/merged'

    if not os.path.exists(dout):
        os.makedirs(dout)

    outfn = dout + f'/{fq1s[-1]}'
    print(outfn)

    fq2 = fq1.replace("_R1_", "_R2_")

#    it1 = grouper(itertools.islice(zipped(fq1), 0, 10*4), 4, '')
#    it2 = grouper(itertools.islice(zipped(fq2), 0, 10*4), 4, '')

    it1 = grouper(zipped(fq1), 4, '')
    it2 = grouper(zipped(fq2), 4, '')

    it = zip(it1, it2)
    with open(outfn, 'wt') as fout:
        for i,ii in it:
            i1, r1, p1, q1 = i
            i2, r2, p2, q2 = ii
            match = reanatomy.search(r1)
            if match:
                intid = f'{match.group("ibar")}_{match.group("spacer")}'
            else:
                intid = 'naibar_spacer'

            i2 = i2.replace(" ", f":{intid} ")
            fout.write("".join([i2, r2, p2, q2]))


        
def chain_to_df(chain_ifn: Optional[str] = 'workdir/umi4c/sample/chain', 
                end: Optional[int] = None) -> None:
    ''' Writes a tab separated file with the mapping points annotated by ibar+spacer, using the umi4c chain file 
    '''
    df = []
    if end is not None:
        end = int(end)
    with open(chain_ifn) as chain:
        for i in itertools.islice(chain,0, end):
            f_ = i.split('\t')
            ibar = f_[0].split(":")[-2]
            rid = ":".join(f_[0].split(":")[0:-2])
            if 'naibar' not in ibar:
                for f in f_[2:]:
                    if "Lin" not in f:
                        f = f.split(':')
                        umi, chrom, coord = f[0:3]
                        df.append((rid, ibar,chrom, coord))
    df = pd.DataFrame(df)
    df.columns = ['rid', 'ibar', 'chrom', 'start']
    df = df.sort_values(['ibar', 'chrom', 'start'], ascending = [True, True, True])
    df['id'] = df['ibar'].astype('category').cat.codes
    df

    dff = df.drop(columns='rid').drop_duplicates()

    bg = dff.id.value_counts(sort=False) > 30
    bg_ = bg[bg].index
    dff = dff[dff.id.apply(lambda x: x in bg_)]
    dff['start'] = dff['start'].map(int) +1
    dff['end'] = dff['start']+1
    dff = dff.loc[:, ['chrom','start', 'end', 'id','ibar',]]
    dff.to_csv(ifn.replace('chain', 'ibar_ann_intervs'), index=False, sep='\t')

def cbra_convert_adj_coords_to_cis(
                adj_ifn: Optional[str] = 'workdir/umi4c/sample/adj.full.coord', 
                end: Optional[int] = None) -> None:
    ''' Writes a tab separated file with the mapping points annotated by ibar+spacer, using the umi4c chain file 
    '''
    # TODO assertion for filenames
    df = pd.read_csv(adj_ifn, sep='\t')
    f_lin = df.chr1 == 'Lin'
    f_cis = df.chr1 == df.chr2
    f_ = f_lin & ~f_cis
    df.loc[f_, 'start1'] = df.loc[f_, 'start2']
    df.loc[f_, 'end1'] = df.loc[f_, 'start1']+int(1)
    df.loc[f_, 'end2'] = df.loc[f_, 'start1']+1
    
    df.loc[f_, 'chrom1'] = df.loc[f_, 'chr2'].apply(lambda x: 'chr'+x)
    df.loc[f_, 'chrom2'] = df.loc[f_, 'chr2'].apply(lambda x: 'chr'+x)
    
    df = df.loc[f_,['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'mol']]
    df.end1 = df.end1.map(int)
    df.end2 = df.end2.map(int)
    df.to_csv(adj_ifn.replace('adj.full.coord','contacts.txt'), sep='\t', index=False, header = True)

def clone_purity_prs(ifn, end = int(5e6), min_reads_per_umi = 4):
    rs= ogtk.UM.pfastq_collapse_UMI(
                        fastq_ifn1 = ifn, 
                        fastq_ifn2 = ifn.replace("_R1_","_R2_"), 
                        umi_len = 26, 
                        min_reads_per_umi = min_reads_per_umi,
                        end = None if end is None else int(end), 
                        fuse = True, 
                        do_rc = True)

    return(rs)

def clone_purity(rs, wl = None):
    ''' returns [pairs_m, mm, mm_d, Z, hcl]
    '''
    pairs = clone_purity_pairs(rs, wl = wl)
    pairsd = clone_purity_code(pairs)
    pairs_m = clone_purity_crosstab(pairsd)
    result = clone_purity_plot(pairs_m) # returs [pairs_m, mm, Z, hcl]
    return(result)

def clone_purity_pairs(rs, tabulate = False, wl=None):
    ''' based on a readset, screen for molecules and fetch their ibars; by default it generates a list of tuples (pairs)
    used to generate a graph; optionally it returns a dataframe of the cell, umi, ibar order
    '''
    rg = ibar_regex(0)
    cells = {}
    for i in rs.umis.keys():
        cell_bc = i[0:16]
        if cell_bc not in cells.keys():
            cells[cell_bc] = []
        cells[cell_bc].append(i)

    cell_size = pd.DataFrame([(k, len(v)) for k,v in cells.items()])
    cell_size.columns = ['cell', 'size']
    cell_size = cell_size.sort_values('size', ascending = False)
    cell_size = cell_size.set_index('cell')

    if wl is not None:
        cell_pool = set(cell_size.index).intersection(set(wl))
        print(f'found {len(cell_pool)} on inter whitelist ({len(wl)}) {len(cell_pool)/len(wl):0.2f}')
        print(f"cells with min obs {cell_size.loc[list(cell_pool)[-1]]}")
    else:
        cells_sampled = int(min(2e4, len(cells)))
        print(cells_sampled)
        cell_pool = cell_size.index.to_list()[0:cells_sampled]
        print(f"traversing {cells_sampled} cells with min obs {cell_size.iloc[cells_sampled-1]}")

    # extract iBAR for each molecule
    outs = []
    for cell in cell_pool:
        for umi in cells[cell]:
            #seq = pd.Series(rs.umis[umi].seqs).value_counts().head(1).index.to_list()[0]
            seqs = rs.umis[umi].seqs
            for seq in seqs[0:min([30, len(seqs)])]:
                match = rg.match(seq)
                if match:
                 outs.append((cell, umi, match.group('ibar')))
    print("done")
    df = pd.DataFrame(outs)
    df.columns = ['cell', 'umi', 'ibar']

    if tabulate:
        return(df)

    # compute a set of pair-wise combinations for all ibars found in each cell
    def ibar_pairs(x):
        return([i for i in itertools.combinations(x['ibar'].value_counts().index, 2)])

    topn = len(wl) if wl is not None else int(2e4)
    print(f'keeping {topn} top cells')
    topc = df.cell.value_counts().head(topn).index
#    rata = df.set_index('cell').loc[topc,:].groupby('cell').apply(lambda x: [i for i in itertools.combinations(x['ibar'].value_counts().index, 2)])
    raw_pairs = df.set_index('cell').loc[topc,:].groupby('cell').apply(ibar_pairs)

    # unpack the list of iBAR pairs
    pairs = []
    for i in raw_pairs.to_list():
        if len(i)>0:
            for k in i:
                pairs.append(k)

    return(pairs)

def clone_purity_code(pairs):
    pairsd = pd.DataFrame(pairs)
    pairsd.columns = ['i', 'ii']
    cats = set(np.hstack(pairsd.to_numpy()))

    pairsd['ai'] = pd.Categorical(pairsd['i'], categories=cats)
    pairsd['aii'] = pd.Categorical(pairsd['ii'], categories=cats)

    return(pairsd)

def clone_purity_crosstab(pairsd):
    pairs_m = pd.crosstab(pairsd['i'], pairsd['ii'])
    pairs_m = pairs_m.reindex(pairs_m.index, axis=1, fill_value=0)
    pairs_m = np.tril(pairs_m) + np.triu(pairs_m, 1).T
    pairs_m = np.tril(pairs_m) +  np.triu(pairs_m.T, 1)
    return(pairs_m)

def clone_purity_mix_pairs(v1, v2):
    mix = []

    for i in [v1, v2]:
        for ii in i:
            mix.append(ii)
    return(mix)

def clone_purity_plot(pairs_m):
    ''' returns [pairs_m, mm, mm_d, Z, hcl]
    '''
    from scipy.cluster import hierarchy
    mm_d = clone_purity_denoise(pairs_m)
    mm = 1/(1+mm_d)
    Z = hierarchy.linkage(mm, 'ward')
    hcl = hierarchy.leaves_list(Z)

    plt.matshow(mm[hcl,:][:, hcl], cmap='RdYlBu')
    plt.matshow(np.log(1+mm[hcl,:][:, hcl]), cmap='rainbow_r')
    return([pairs_m, mm, mm_d, Z, hcl])

def clone_purity_denoise(pairs_m):
    marginal_sum =np.sum(pairs_m, axis=1) 
    #v, cts = dict(zip(*numpy.unique(marginal_sum, return_counts=True)))
    #his = dict(zip(*numpy.unique(marginal_sum, return_counts=True)))

    his = pd.Series(marginal_sum).value_counts(sort=False)
    ii = sorted(his.index)
    cs = np.cumsum(his[ii].to_numpy())
    trend = cs[1:]-cs[0:-1]
    mfilter = cs[np.argmin(trend)]

    #mfilter = np.quantile(marginal_sum, 0.5)
    #filter_mask = np.where(marginal_sum <= mfilter)[0]
    #filter_mask = np.where(marginal_sum <= mfilter)[0]
    filter_mask = marginal_sum >= mfilter
    pairs_m = pairs_m[filter_mask,:][:, filter_mask]

    return(pairs_m)


