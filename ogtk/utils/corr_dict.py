import polars as pl
from logging import warning
from typing import Sequence,Optional

def generate_cbc_correction_dictionaries(path_to_bam, force = False, verbose = False, chunk_size=1e6):
    '''
    Generates a cellbarcode correction dictionary based on the cellranger barcoded BAM tags
    TODO: there is a small fraction of uncorrected cell barcodes that seem to map to more than one corrected cell barcode
    '''
    import pysam
    import os
    import polars as pl
    import rich

    ofn = path_to_bam.replace('.bam', '.parquet')
    print(f'scanning {path_to_bam}')
    # load a pre-computed dictionary of found
    if os.path.exists(ofn) and not force:
        rich.print('pre-computed map found, returning file name. Use force to regenerate')
        return(ofn)
    else:
    # make it if forced or not found
        entries = []
        if verbose:
            rich.print(f'Opening bam file {path_to_bam} to screen for correction pairs', end = '....')
        print('') 
        bam = pysam.AlignmentFile(path_to_bam)

        df = pl.DataFrame(columns=[('CR', pl.Utf8), ('CB', pl.Utf8)] )
        
        for read in bam:
            if read.has_tag('CR') and read.has_tag('CB'):
                entries.append((read.get_tag("CR"), read.get_tag("CB")))
            

            if len(entries)%chunk_size ==0:
                df = (df.vstack(pl.DataFrame(entries, orient='row', columns=[('CR', pl.Utf8), ('CB', pl.Utf8)] ))
                     .unique()
                    )
                rich.print(df.shape)
                entries = []

        df = df.unique().write_parquet(ofn)

        if verbose:
            print(f'done')

        return(ofn)
def correct_cbc_pl(df: pl.DataFrame, corr_dic_df: pl.DataFrame)-> pl.DataFrame:
    ''' Provided a polars dataframe and a correction dictionary polars dataframe, return a corrected version of the input data frame
    '''
    df = df.join(corr_dic_df, left_on='cbc', right_on='CR', how='left')
    cbc_cr = df.with_columns(pl.col("cbc")==pl.col("CB"))["cbc"].mean()
    print(f'cbc ==>CB {cbc_cr}')
    df = df.drop('cbc').rename({'CB':'cbc'})

    # rather convoluted way of not using np.hstack()
    col_i = pl.Series([pl.Series([-1]), pl.Series(range(0, len(df.columns)-1))]).explode()
    columns  = pl.Series(df.columns)[col_i].to_list()

    return(df.select(columns))

def load_corr_dict(ifn):
    ''' load a parquet file and sanitize cbc names (remove "-1")
    '''
    pp=pl.scan_parquet(ifn).with_columns(pl.col('CB').str.replace('-1', '')).collect()
    return(pp)

def generate_correction_dictionaries(pickle_ofn, path_to_bam, force = False, verbose = False):
    '''
    Generates a cellbarcode correction dictionary based on the cellranger barcoded BAM tags
    TODO: there is a small fraction of uncorrected cell barcodes that seem to map to more than one corrected cell barcode
    '''
    import pysam
    import itertools
    import pickle 
    import pyaml
    import os

    # load a pre-computed dictionary of found
    if os.path.exists(pickle_ofn) and not force:
        dd =  pickle.load(open(pickle_ofn, 'rb')) 
        print(dd['qc'])
        return(dd)
    else:
    # make it if forced or not found
        bam = pysam.AlignmentFile(path_to_bam)
        cbc_corr = []
        umi_corr = []

        if verbose:
            print(f'Opening bam file {path_to_bam} to screen for correction pairs', end = '....')

        for read in bam:
            if read.has_tag('CR') and read.has_tag('CB'):
                cbc_corr.append((read.get_tag('CR'), read.get_tag('CB')))
            if read.has_tag('UR') and read.has_tag('UB'):
                umi_corr.append((read.get_tag('UR'), read.get_tag('UB')))

        if verbose:
            print(f'done')

        dict_cbc = dict(cbc_corr)
        dict_umi = dict(umi_corr)

        # quick QCs
        before_dict_cbc = set([i[1] for i in cbc_corr])
        before_dict_umi = set([i[1] for i in umi_corr])

        qc_str = f'Missed barcode cases = {len(before_dict_cbc) - len(set(dict_cbc.values()))}'
        qc_str = qc_str + f'Missed umi cases = {len(before_dict_umi) - len(set(dict_umi.values()))}'

        print(qc_str)

        dd = {'cell_barcodes':dict_cbc, 'umi_barcodes':dict_umi, 'qc':qc_str}
        with open(pickle_ofn, 'wb') as handle:
            pickle.dump(dd, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return(dd)

