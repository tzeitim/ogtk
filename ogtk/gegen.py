import itertools
#import mappy # consider using minimap's api instead of a system invocation?
import pysam
import re
import copy
import pdb
import IPython

#TODO: - manual coordinates ; - add chr to the mini ref
def sitedb_generate_miniref(sitedb_yaml):
    '''
    Using the sitedb generate a fasta file (mini reference) for mapping reads from endogenous targets, for example.
    Depends on a yaml file that contains a ['loci'] level with 'chr' 'gene_start' 'gene_end' keys.
    '''
    import pyaml
    import uuid
    import pyfaidx
    import subprocess

    tmp_coords = f"/tmp/{uuid.uuid4()}.temp_coords"
    tmp_fa = tmp_coords.replace('coords', 'fa')
    print(tmp_coords)

    sitedb = pyaml.yaml.load(open(sitedb_yaml), Loader = pyaml.yaml.Loader)
    loci = sitedb['loci']
    genes = list(loci.keys())
    path_2bit = sitedb['ref']['path_2bit']
    out_fa =  sitedb['ref']['mini_ref']

    # generate fasta file with the small reference sequence
    with open(tmp_coords, 'w') as coords:
        for gene in genes:
            chrom = loci[gene]['chr'].replace('','')
            start = loci[gene]['gene_start']
            end =   loci[gene]['gene_end']
            coords.write(f'{chrom}:{start}-{end}'+'\n')

    # Call twobit2fasta 
    twobit_cmd = f'twoBitToFa -noMask -seqList={tmp_coords} {path_2bit} {tmp_fa}'
    #print(twobit_cmd)
    subprocess.run(twobit_cmd.split())

    # provide names for mini reference
    FF = pyfaidx.Fasta(tmp_fa) 
    with open(out_fa, 'w') as outfa:
        for entry in FF.keys():
            c_chrom = entry.split(':')[0]
            c_start = int(entry.split(':')[1].split('-')[0])
            for gene in genes:
                if loci[gene]['gene_start'] == c_start and loci[gene]['chr'] == c_chrom:
                    outfa.write(f'>{loci[gene]["name"]}\n{FF[entry][:]}\n')

    # clean the mess
    subprocess.run(f'rm {tmp_fa} {tmp_coords}'.split())
    print(f"miniref created at:\n{out_fa}")


def extract_reads(sitedb_yaml, use_unmapped = True, FORCE = False, with_chr = True):
    '''
    Expects a sitedb yaml file that provides hits reference in FASTA format.
 
    Required YAML fields:
        xp.workdir
        xp.name 
        xp.tenex_bam
        xp.all_hit_bam
        ref.mini_ref
 
    Steps:
        1. Get unmapped reads into fastq file retaining CR/UR tags
        2. Extract reads around genes of interest into fastq with CR/UR tags 
        3. Concatenate hits and unmapped
        4. Map using minimap2
        5. Profit
    '''
    import os
    import time
    import subprocess
    import pyaml
    import ogtk

    # globals
    sitedb = pyaml.yaml.load(open(sitedb_yaml), Loader = pyaml.yaml.Loader)
    workdir = sitedb['xp']['workdir'] 
    xp_name = sitedb['xp']['name']
    tenex_bam = sitedb['xp']['tenex_bam']
    all_hits_bam = sitedb['xp']['all_hits_bam'] 

    umapped_fastq = f"{workdir}/{xp_name}_umapped.fastq.gz"
    hits_genes_fastq = f"{workdir}/{xp_name}_hits_genes.fastq.gz"
    mini_ref =  sitedb['ref']['mini_ref']


    if not os.path.exists(workdir):
        os.makedirs(workdir)
                       
    # step 1
    umap_cmd = f'samtools view -h -b -f 4 {tenex_bam}'
    umap_fastq_cmd = f'samtools fastq -T CR,UR'
    umap_format_cmd = f'sed -e s/\\t/_/g'
    umap_gzip = 'gzip'
    
    if use_unmapped:

        if os.path.exists(umapped_fastq): 
            print(f'A previous instance of {umapped_fastq} was found! skipping its computation. Use the FORCE = True argument to override this assertion...', end  = '')
            if FORCE:
                print(f'Removing it since FORCE = {FORCE}.')
                subprocess.run(f'rm {umapped_fastq}'.split())
            else:
                print(f'Use the FORCE = True argument to override this assertion')

        if os.path.exists(umapped_fastq) == False or FORCE == True: 
            print(umap_cmd, "|", umap_fastq_cmd , "|", umap_format_cmd , "| gzip >", umapped_fastq)
            c1 = subprocess.Popen(umap_cmd.split(), stdout = subprocess.PIPE)
            c2 = subprocess.Popen(umap_fastq_cmd.split(), stdin=c1.stdout, stdout=subprocess.PIPE)
            c3 = subprocess.Popen(umap_format_cmd.split(), stdin=c2.stdout, stdout=subprocess.PIPE)
            c4 = subprocess.Popen(umap_gzip.split(), stdin=c3.stdout, stdout=open(umapped_fastq, 'w'))
            c1.stdout.close()
            c2.stdout.close()
            c3.stdout.close()
            c4.communicate()
        
            while any(map(lambda x: x.poll() == None, [c1,c2,c3,c4])):
                time.sleep(1)

            map(lambda x: x.terminate(), [c1,c2,c3,c4])
    else:
        cmd_touch = f'touch {umapped_fastq}'
        subprocess.run(cmd_touch.split())

    # step 2
    # step 2.1 get all reads that map to genes of interest
    coords_str = sitedb_return_coord_str(sitedb_yaml, with_chr = with_chr)
    hits_cmd = f'samtools view -h -b {tenex_bam} {coords_str}'
    
    c1 = subprocess.Popen(hits_cmd.split(), stdout = subprocess.PIPE)
    c2 = subprocess.Popen(umap_fastq_cmd.split(),stdin=c1.stdout, stdout = subprocess.PIPE)
    c3 = subprocess.Popen(umap_format_cmd.split(), stdin=c2.stdout, stdout=subprocess.PIPE)
    c4 = subprocess.Popen(umap_gzip.split(), stdin=c3.stdout, stdout=open(hits_genes_fastq, 'w'))
    c1.stdout.close()
    c2.stdout.close()
    c3.stdout.close()
    c4.communicate()
   
    while any(map(lambda x: x.poll() == None, [c1,c2,c3,c4])):
        time.sleep(1)

    map(lambda x: x.terminate(), [c1,c2,c3,c4])

    # step 3
    # not true but almost
    #'minimap2 -a --cs --splice --sam-hit-only  ./mini_ref.fa umapped.fastq.gz hits_genes.fastq.gz | samtools sort | samtools view -b > fffyyy.bam'
    # added --sr --splice-flank=no
    ogtk.mp.mp2_map(ref = mini_ref, query = " ".join([hits_genes_fastq, umapped_fastq]), options = '-a --cs --splice --splice-flank=no --sam-hit-only --sr', outfn = all_hits_bam, verbose = True)

       

   
    
def sitedb_return_coord_str(path_to_yaml, with_chr = True):

    import pyaml
    
    sitedb = pyaml.yaml.load(open(path_to_yaml), Loader = pyaml.yaml.Loader)
    genes = sitedb['loci'].keys()
    coords = []
    loci = sitedb['loci']
    for gene in genes:
        chrom = loci[gene]['chr']  
        chrom = chrom if with_chr else chrom.replace('chr', '')
        start = loci[gene]['gene_start']
        end   = loci[gene]['gene_end']
        coord = f"{chrom}:{start}-{end}"
        coords.append(coord)
    coords_str = ' '.join(coords)
    return(f'{coords_str}')

    # step 2
    

def generate_correction_dictionaries(sitedb_yaml, FORCE = False, verbose = False):
    '''
    Generates a cellbarcode correction dictionary based on the cellranger barcoded BAM tags
    TODO: there is a small fraction of uncorrected cell barcodes that seem to map to more than one corrected cell barcode
    '''
    import pysam
    import itertools
    import pickle 
    import pyaml
    import os

    sitesdb = pyaml.yaml.load(open(sitedb_yaml), Loader = pyaml.yaml.Loader)

    workdir = sitesdb['xp']['workdir']
    xpname  = sitesdb['xp']['name']

    pickle_ofn = f'{workdir}/{xpname}.corrdic.pickle'
    
    # load a pre-computed dictionary of found
    if os.path.exists(pickle_ofn) and not FORCE:
        dd =  pickle.load(open(pickle_ofn, 'rb')) 
        print(dd['qc'])
        return(dd)
    else:
    # make it if forced or not found
        path_to_bam = sitesdb['xp']['tenex_bam']
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

def call_alleles(sitedb_yaml, correction_dictionaries, window = 120/2, n_reads = None):
    '''
    ifn_bam = output from minimap of a concatenated set of reads of hits on targets 
    and unmapped reads. Reads are expected to include the CR and UR tags from the 
    10x bam, e.g.:
    
    @NS500648:472:HF7WWAFX2:1:11102:17555:11614_CR:Z:AAAAAAAAAAAACAAT_UR:Z:AAAAAAAAAAAA
    
    Parameters:
        n_reads: number of reads to consume from the bam; None means all.
    
    Real call example:
    _call_alleles(ifn_bam, ofn, loci, correction_dictionaries)
    '''
    import pyaml
    sitedb = pyaml.yaml.load(open(sitedb_yaml), Loader=pyaml.yaml.Loader)
    ifn_bam = sitedb['xp']['all_hits_bam'] 
    ofn = sitedb['xp']['allele_tab']
    loci = sitedb['loci']
    
    # this is just a wrapper so we call here the real function
    _call_alleles(ifn_bam, ofn, loci, correction_dictionaries, n_reads, window)

def _call_alleles(ifn_bam, ofn, loci, correction_dictionaries = None, n_reads = None, window = 120/2, region = None):
    weird = []
    bam = pysam.AlignmentFile(ifn_bam)
    if region != None:
        bam_i = bam.fetch(region)
    else:
        bam_i = itertools.islice(bam, 0 , n_reads)

    out = open(ofn, 'w')

    with open(ofn, 'w') as out:
        out.write('\t'.join(['chr', 'start', 'cr', 'cb', 'ur','ub', 'molid', 'cs', 'cigar', 'iv', 'dv', 'hallele', 'sallele', 'seq'])+'\n')
        for read in bam:
            if not read.has_tag('cs'):
                weird.append(read.query_name)
            else:
                cs_str = read.get_tag('cs')

                foc_chrom = read.reference_name
                chrom_str = foc_chrom.replace('chr', '') 
                if chrom_str in loci.keys():
                    # here we need to add code to clip out the 'chr' that willbe added in the future #TODO
                    guide_mini = loci[chrom_str]['start'] - loci[chrom_str]['gene_start'] 
                    foc_start = guide_mini - window 
                    foc_end =   guide_mini + window  

                    if read.reference_start >= foc_start and read.reference_start <= foc_end:
                        import pdb
                        #pdb.set_trace()
                        #print(chrom_str, loci[chrom_str]['start'], foc_start, '<',read.reference_start, '<', foc_end )
                        allele = cs_str.split(':')
                        mm_al = []
                        old_allele = copy.copy(allele)
                        current = 0
                        for i,v in enumerate(allele):
                            if i == 0:
                                allele[i] = read.reference_start
                                current = allele[i]
                            else:
                                if "-" in v:
                                    match = re.search("(?P<numeric_part>[0-9]+)-(?P<deletion>.+)", v)
                                    if match:
                                        start_step = int(match.group("numeric_part"))
                                        deletion = match.group('deletion')
                                        if "*" in deletion:
                                            mms = deletion.split('*')
                                            allele[i] =f"{current+start_step}-{mms[0]}"
                                            current = current + start_step + len(mms[0])
                                            for mm in mms[1:]:
                                                mm_al.append(f'{current}*{mm}')
                                                current += 1
                                        # change this to replace 
  #                                      print(f"{current+start_step}-{deletion}")
                                        else:
                                            allele[i] =f"{current+start_step}-{deletion}"
                                            current = current + start_step + len(deletion)
                                elif "+" in v:
                                    match = re.search("(?P<numeric_part>[0-9]+)\+(?P<insertion>.+)", v)
                                    if match:
                                        start_step = int(match.group("numeric_part"))
                                        insertion = match.group('insertion')
                                        if "*" in insertion:
                                            mms = insertion.split('*')
                                            allele[i] = f'{current+start_step}+{mms[0]}'
                                            current = current + start_step
                                            for mm in mms[1:]:
                                                mm_al.append(f'{current}*{mm}')
                                                current += 1 
                                        else:
 #                                       print(f"{current+start_step}+{insertion}")
                                            allele[i] = f"{current+start_step}+{insertion}"
                                            current = current + start_step
                                elif "~" in v:
                                    match = re.search("(?P<numeric_part>[0-9]+)~(?P<splice_donor>[actgn]{2})(?P<intron>[0-9]+)(?P<splice_acc>.+)", v)
                                    if match:
                                        start_step = int(match.group('numeric_part'))
                                        intron = match.group('intron')
                                        splice_acc = match.group('splice_acc')
                                        if "*" in splice_acc:
                                            mms = splice_acc.split('*')
                                            intron_str = match.group('splice_donor') + match.group('intron') + mms[0]
                                            allele[i] = f'{current+start_step}ii{intron_str}'
                                            current = current + start_step + int(intron)
                                            for mm in mms[1:]:
                                                mm_al.append(f'{current}*{mm}')
                                                current += 1 
                                        else:
                                            intron_str = match.group('splice_donor') + match.group('intron') + match.group('splice_acc')
                                            allele[i] = f'{current+start_step}ii{intron_str}' # we use i no keep a record that this allle was derived from an intron type
                                            current = current + start_step + int(intron)
                                elif "*" in v:
                                    match = re.search("(?P<numeric_part>[0-9]+)\*(?P<mismatch>.+)",v)
                                    if match:
                                        start_step = int(match.group("numeric_part"))
                                        mismatch = match.group('mismatch')
#                                        print(f"{current+start_step}*{mismatch}")
                                        mm_al.append(f"{current+start_step}*{mismatch}")
                                        current = current + start_step + 1
                                        allele[i] = 'MM'

                                else:
                                    match = re.search("(?P<numeric_part>[0-9]+)",v)
                                    if match:
                                        start_step = int(match.group("numeric_part"))
                                        current = current + start_step
                                        allele[i] = f"{current}"
                        allele = [i for i in allele if i != 'MM']
                        al_cs = ':'.join(allele[1:-1])
                        mm_al_cs = ':'.join(mm_al)
                        
                        iv = []
                        dv = []
                        for op,ll in read.cigartuples:
                            if op == 1:
                                iv.append(str(ll))
                            if op == 2:
                                dv.append(str(ll))

                        if len(iv) == 0:
                            iv = '0'
                        if len(dv) == 0:
                            dv = '0'

                        iv = ",".join(iv)
                        dv = ','.join(dv)
                        if al_cs == '':
                            al_cs = "wildtype"
                        if mm_al_cs == '' or mm_al_cs == None:
                            mm_al_cs = 'NA'

                        chrom = read.reference_name
                        if chrom is None:
                            pdb.set_trace()
                        map_start = str(read.reference_start)

                        if read.has_tag("CR"):
                            cr = read.get_tag("CR")
                        else:
                            cr = read.query_name.split(':Z:')[1].split("_")[0]
                        if read.has_tag("UR"):
                            ur = read.get_tag("UR")
                        else:
                            ur = read.query_name.split(':Z:')[2]

                        if correction_dictionaries != None:
                            if cr in correction_dictionaries['cell_barcodes'].keys():
                                cb = correction_dictionaries['cell_barcodes'][cr]
                            else:
                                cb = "NA"
                            if ur in correction_dictionaries['umi_barcodes'].keys():
                                ub = correction_dictionaries['umi_barcodes'][ur]
                            else:
                                ub = "NA"
                        else:
                            ub = 'NA'
                            cb = 'NA'

                        molid = cr + ur
                        cs = cs_str
                        cigar = read.cigarstring
                        hallele = al_cs
                        sallele = mm_al_cs
                        seq = read.seq
                        if seq is None:
                            weird.append(read)
                        else:
                            str_out = '\t'.join([chrom, map_start, cr, cb, ur, ub, molid, cs, cigar, iv, dv, hallele, sallele, seq])
                            #print(str_out)
                            out.write(str_out+'\n')
                            #print(read.query_name, '\n', read.cigarstring, '\n', read.cigartuples, '\n', allele, '\n', old_allele, '\n' , mm_al , '\n This is the allele cs', al_cs, '\n This is the mismatch cs', mm_al_cs, '\n', read.seq)
        if len(weird)>0:
            print(f'Warning: a total of {len(weird)} reads without CS tag were found! or that do not have a seq (?)')


