##1 pool3
#basedir = "/local/users/polivar/src/projects/mltracer/datain/20201229_bulk_rna_pilot2/GSTLT_bulk2/"
#sbcs = ['AGACTC', 'AGCTTC', 'CATGAG', 'CAGATC', 'TCACAG', 'GTCTAG', 'GTTGCA', 'GTGACA']
#pool1_sbcs = ['AGACTC', 'AGCTTC', 'CATGAG', 'CAGATC', 'TCACAG']
#pool2_sbcs = ['GTCTAG', 'GTTGCA', 'GTGACA', 'AGACTC', 'AGCTTC']
#
#end = None
#prs = ogtk.UM.pfastq_collapse_UMI(fastq_ifn1=basedir+'pool2_S2_R1_001.fastq.gz',
#                                  fastq_ifn2=basedir+'pool2_S2_R2_001.fastq.gz',
#                                  umi_len=28,
#                                  end= end)
##2 pool2
#with open('pool2_all.txt', 'w') as outf:
#    tot_umis = len(prs[1].umis)
#    total = 0
#
#    for k,umi in prs[1].umis.items():
#        sample = k[-6:]
#
#        if sample in sbcs:
#            total+=1
#            print("\r"+str(total/tot_umis), end='')
#            megastr = ''.join(umi.seqs)
#            if "aTGGagtcg".upper() in megastr or 'accacctgttcctgtagaaa'.upper() in megastr or "CGAGCGCTATGAGCGACTAT" in
# megastr:
#                contigs = ogtk.UM.wrap_assemble(umi, k = 29)
#                if contigs:
#                    for con in contigs:
#                        outf.write('\t'.join([sample, umi.umi, "contig", con[1].strip()])+'\n')
#                for s,seq in enumerate(umi.seqs):
#                    outf.write('\t'.join([sample, umi.umi, str(s), seq])+'\n')
#
#
#
#
#

def fastq_to_contig(name, ifn1, ifn2, end = None, verbose = False):
    #1 pool1
    #basedir = "/local/users/polivar/src/projects/mltracer/datain/20201229_bulk_rna_pilot2/GSTLT_bulk2/"
    import ogtk
    import pandas as pd
    sbcs = ['AGACTC', 'AGCTTC', 'CATGAG', 'CAGATC', 'TCACAG', 'GTCTAG', 'GTTGCA', 'GTGACA']
    #pool1_sbcs = ['AGACTC', 'AGCTTC', 'CATGAG', 'CAGATC', 'TCACAG']
    #pool2_sbcs = ['GTCTAG', 'GTTGCA', 'GTGACA', 'AGACTC', 'AGCTTC']
    
    prs = ogtk.UM.pfastq_collapse_UMI(
                fastq_ifn1=ifn1,
                fastq_ifn2=ifn2,
                umi_len=28,
                end = end)
    
    
    #2 pool1
    with open(f'{name}_all.txt', 'w') as outf:
        tot_umis = len(prs[1].umis)
        total = 0
    
        for k,umi in prs[1].umis.items():
            sample = k[-6:]
    
            if sample in sbcs:
                total+=1
                print(f"\r{100*(total/tot_umis):.4f}", end='')
                megastr = ''.join(umi.seqs)
                if "aTGGagtcg".upper() in megastr or 'accacctgttcctgtagaaa'.upper() in megastr or "CGAGCGCTATGAGCGACTAT" in megastr:
                    contigs = ogtk.UM.wrap_assemble(name, umi, k = 29, verbose = verbose)
                    if contigs:
                        for con in contigs:
                            contig_str = con[1].strip()
                            outf.write('\t'.join([sample, umi.umi, "contig", str(len(contig_str)), contig_str])+'\n')
                    else:
                        counts = pd.Series(umi.seqs).value_counts()
                        seq = counts.index.to_list()[0]
                        counts = counts.to_list()[0]
                        outf.write('\t'.join([sample, umi.umi, "read", str(len(seq)), seq])+'\n')
                    #for s,seq in enumerate(umi.seqs):
                    #    outf.write('\t'.join([sample, umi.umi, str(s), str(len(seq)), seq])+'\n')
    
    
def fastq10x_to_contig(name, outdir, ifn1, ifn2, K = 80, end = None, just_set = False, do_rc = False, verbose = False, mincontig='auto', mincountseed=100, fuse_readsets = False):
    import ogtk
    import pandas as pd
    import re
    import os

    if not os.path.exists(outdir):
        print(f'Creating output dir')
        os.makedirs(outdir)
    print(f'Dumping results to {outdir}')
    K = int(K) 
    if mincontig != 'auto':
        mincontig = int(mincontig)
    mincountseed = int(mincountseed)
    polyT = re.compile("TTTTTTTTTTTTTTTTTTTTTT")
    
    prs = ogtk.UM.pfastq_collapse_UMI(
                fastq_ifn1=ifn1,
                fastq_ifn2=ifn2,
                umi_len=28,
                end = end,
                do_rc = do_rc,
                fuse = fuse_readsets)

    if just_set:
        if fuse_readsets:
            prs = prs[1]
        return(prs)
 
    if not fuse_readsets:
        prs = prs[1]

   
    with open(f'{outdir}/{name}_all.txt', 'w') as outf, open(f'{outdir}/{name}_contigs.fa', 'w') as fout, open(f'{outdir}/{name}_contigs.fastq', 'w') as fqout:
        tot_umis = len(prs.umis)
        total = 0
    
        sorted_umis = prs.umi_counts().index.to_list()
        
        for k in sorted_umis:
            umi = prs.umis[k]
            sample = k[-6:]
    
            total+=1
            print(f"\r{100*(total/tot_umis):.4f}", end='')
            megastr = ''.join(umi.seqs)
            counts = pd.Series(umi.seqs).value_counts()
            contigs = ogtk.UM.wrap_assemble(name= name, umi = umi, 
                    k = K, verbose = verbose, mincontig=mincontig, mincountseed=mincountseed)
            cbc_umi = umi.umi
            cbc = cbc_umi[0:16]
            umi_str = cbc_umi[16:]

            if contigs:
                for con in contigs:
                    contig_str = con[1].strip()
                    if polyT.search(contig_str):
                        contig_str = ogtk.UM.rev_comp(contig_str)
                    contig_cov = con[0].split('=')[2].split(',')[0]
                    contig_len = len(contig_str)
                    outf.write('\t'.join([cbc, umi_str, cbc_umi, "contig", contig_cov, str(len(umi.seqs)) ,str(contig_len), contig_str])+'\n')
                    fout.write(">"+'_'.join([cbc, umi_str, cbc_umi, "contig", contig_cov, str(len(umi.seqs)), str(contig_len), "_MU:Z:"+ sample])+'\n'+contig_str+'\n')
                    fqout.write("@"+'_'.join([cbc, umi_str, cbc_umi, "contig", str(contig_len), "_MU:Z:"+ sample])+'\n'+contig_str+'\n+\n'+str(contig_len*"E")+"\n")
            else:
                seq = counts.index.to_list()[0]
                counts = counts.to_list()[0]
                #outf.write('\t'.join([sample, umi.umi, "read", str(len(seq)), seq])+'\n')
                outf.write('\t'.join([cbc, umi_str, cbc_umi, "read", str(counts), str(len(umi.seqs)), str(len(seq)), seq])+'\n')
            outf.flush()
            #for s,seq in enumerate(umi.seqs):
            #    outf.write('\t'.join([cbc, umi_str, cbc_umi, str(s), "NA", str(len(umi.seqs)), str(len(seq)), seq])+'\n')
    
 


def bam10x_to_contig(name, ifn,  K = 80, end = None, just_set = False, verbose = False, mincontig = 100, mincountseed=1, try_contig = True, report_rep_read = False):
    import ogtk
    import pandas as pd

    K = int(K) 
    mincontig = int(mincontig)
    mincountseed = int(mincountseed)
    prs = ogtk.UM.bam10x_collapse_UMI(
                bam_ifn = ifn,
                umi_len=28,
                end = end)
    if just_set:
        return(prs)
    
    with open(f'{name}_all.txt', 'w') as outf, open(f'{name}_contigs.fa', 'w') as fout, open(f'{name}_contigs.fastq', 'w') as fqout:
        tot_umis = len(prs.umis)
        total = 0
    
        sorted_umis = prs.umi_counts().index.to_list()
        
        for k in sorted_umis:
            umi = prs.umis[k]
    
            total+=1
            print(f"\r{100*(total/tot_umis):.4f}", end='')
            megastr = ''.join(umi.seqs)
            #contigs = ogtk.UM.wrap_assemble(name= name, umi = umi, k = K, verbose = verbose, mincontig=mincontig, mincountseed=mincountseed)
            contigs = False
            cbc_umi = umi.umi
            cbc = cbc_umi[0:16]
            umi_str = cbc_umi[16:]
            reads = len(umi.seqs)

            if contigs and try_contig:
                for con in contigs:
                    contig_str = con[1].strip()
                    contig_cov = con[0].split('=')[2].split(',')[0]
                    contig_len = len(contig_str)
                    outf.write('\t'.join([cbc, umi_str, cbc_umi, "contig", contig_cov, str(reads) ,str(contig_len), contig_str])+'\n')
                    fout.write(">"+'_'.join([cbc, umi_str, cbc_umi, "contig", contig_cov, str(reads), str(contig_len)])+'\n'+contig_str+'\n')
                    fqout.write("@"+'_'.join([cbc, umi_str, cbc_umi, "contig", str(contig_len)])+'\n'+contig_str+'\n+\n'+str(contig_len*"E")+"\n")
            if report_rep_read:
                counts = pd.Series(umi.seqs).value_counts()
                seq = counts.index.to_list()[0]
                counts = counts.to_list()[0]
                #outf.write('\t'.join([sample, umi.umi, "read", str(len(seq)), seq])+'\n')
                outf.write('\t'.join([cbc, umi_str, cbc_umi, "read", str(counts), str(reads), str(len(seq)), seq])+'\n')



            #for s,seq in enumerate(umi.seqs):
            #    outf.write('\t'.join([cbc, umi_str, cbc_umi, str(s), "NA", str(len(umi.seqs)), str(len(seq)), seq])+'\n')
    
def format_alleles(name, states_ifn, fa_ofn, do_rc, anchor1 = "GATCCACCAGGCCCTGAAG", ref_path = "/local/db/danrer11/trackdb/seq/Brainbow2.0L.fa", length = 90, ref_name = 'brainbow', verbose = False, trim_start = 1180, trim_end = 1270, keep_tmp = False, force = False, gapopen = 50, gapextend = 0.5, mode = ["water", "mini2", "needleman"][0]):
    from pyfaidx import Fasta
    import re
    import pdb
    import ogtk
    import subprocess
    import pysam

    trim_start = int(trim_start)
    trim_end = int(trim_end)

    length = int(length)
    anchor2 = "TBD"
    #fa_ofn = 'fclean_out.fa'
    fa_1 = fa_ofn+'.fstates.tmp'
    rr = re.compile(anchor1)
    filtered = []

    # export records to fasta format
    with open(states_ifn) as ifn, open(fa_1, 'w') as outf:
        for line in ifn.readlines():
            f = line.split('\t')
            is_contig = f[3]
            seq = ogtk.UM.rev_comp(f[7]) if do_rc else f[7]

            #if is_contig == "contig" or force:
            if False:
                match = rr.search(seq)
                if match:
                    start = match.span()[0]
                    end = start + length
                    state = seq[start:end].strip()
                    if len(state) >= length:
                        ostr1 = ">"+"_".join([i for i in f[0:7]])+"\n"
                        ostr2  = ostr1+state+'\n'
                        outf.write(ostr2)
                    else:
                        filtered.append(state)
            else:
                state = seq
                if len(state) >= length:
                    ostr1 = ">"+"_".join([i for i in f[0:7]])+"\n"
                    ostr2  = ostr1+state+'\n'
                    outf.write(ostr2)


    if mode == "water":
        #water -asequence /local/db/danrer11/trackdb/seq/Brainbow2.0L.fa  -sformat1 fasta -bsequence fbli_clean.fa.fstates.tmp -sformat2 fasta -gapopen 55  -gapextend 0.02  -outfile sam -aformat3 sam
        print(">>>",fa_ofn)
        ogtk.mp.embos_align_reads_to_ref(
                        name=name, fa_ifn=fa_1, fa_ofn = fa_ofn, 
                        mode = "waterman",
                        ref_path = ref_path, ref_name = ref_name, 
                        trim_start = trim_start, trim_end = trim_end,  
                        gapopen=gapopen, gapextend=gapextend, 
                        verbose = verbose, thorough = True)
   
        print(">>>>", fa_1, fa_ofn)
        format_alleles_from_sam(ref_path, fa_ofn) 



    if mode == "needleman":
        ogtk.mp.embos_align_reads_to_ref(
                        name=name, fa_ifn=fa_1, fa_ofn = fa_ofn, 
                        mode = 'needleman',
                        ref_path = ref_path, ref_name = ref_name, 
                        trim_start = trim_start, trim_end = trim_end,  
                        gapopen=gapopen, gapextend=gapextend, 
                        verbose = verbose, thorough = True)

        with open(fa_ofn.replace('fa', 'concise'), 'w') as tab:
            tab.write("\t".join(['cbc', 'umi', 'molid', 'contig_str', 'contig_cov', 'reads', 'contig_len', 'allele']) + '\n')

            FF = Fasta(fa_ofn, as_raw=True)
            for entry in FF.keys():
                cbc, umi, molid, contig_str, contig_cov, reads, contig_len = entry.split('_')
                tab.write("\t".join([cbc, umi, molid, contig_str, contig_cov, reads, contig_len, FF[entry][:]]) + '\n')


    #ogtk.mp.embos_align_reads_to_ref(name=name, fa_ifn=fa_1, fa_ofn = fa_ofn, ref_path = ref_path, ref_name = ref_name,  gapopen=20, gapextend=0.55, verbose = verbose)
    #ogtk.mp.embos_align_reads_to_ref(name='fclean', fa_ifn='fclean.states', fa_ofn = fa_ofn, ref_path = ref_path, ref_name = 'brainbow', omit_ref = False)
    if mode == "mini2":
        with open(fa_ofn.replace('fa', 'concise'), 'w') as tab:
            tab.write("\t".join(['cbc', 'umi', 'molid', 'contig_str', 'contig_cov', 'reads', 'contig_len', 'seq', 'cs', 'allele', 'cigar']) + '\n')
            # TODO add to ogtk
            bam_fn = fa_ofn.replace('fa', 'concise.bam')
            ogtk.mp.mp2_map(ref=ref_path, query = fa_1, outfn = bam_fn, options= "-a --cs -A 5" )
            bam = pysam.AlignmentFile(bam_fn) 
            cs = re.compile("(?P<offset>[0-9]+)(?P<op>\*[a-z][a-z]+|[=\+\-][A-Za-z]+)|(?P<matches>[0-9]+)")
            for i in bam:
                cbc, umi, molid, contig_str, contig_cov, reads, contig_len = i.query_name.split('_')

                if molid == "CTCAGAAGTAGTTAGACCAGCGAAAGTC":
                    target = True
                else:
                    target = False

                if i.has_tag('cs') and i.seq != None:
                    cs_str = i.get_tag('cs')
                    start = i.reference_start
                    current = start
                    if target:
                        print(cs_str, current, start, i.seq, i.cigarstring)
                    #print(i.seq, start, i.get_tag('cs'))
                    if i.seq == None:
                        print("culerosa")
                        pdb.set_trace()
                    ops = [j for j in cs_str.split(":")]
                    allele = []
                    for j in ops:
                        match = cs.search(j)
                        if match:
                            if match.groups()[2] == None:
                                offset = int(match.group('offset'))
                                op =match.group('op')
                                op_step = len(op) - 1
                                if target:
                                    print(f"j{j}, {current}, {offset}, {current+offset}, {current+offset+op_step}")
                                current = current + offset
                                if current < 1515:

                                    allele.append("".join([str(current), op[0],str(op_step)]))
                                    #allele.append("".join([str(current), match.group('op')]))
                                else:
                                    if target:
                                        print(current, j)
                                current = current + op_step
                    if len(allele) == 0:
                        allele = ['maybe_wt', "wt"] 
        
                    allele = ":".join(allele)
                    if target and False:
                        print(allele)
                        pdb.set_trace()

                    if i.cigarstring == None:
                        cigar= "None"
                    else:
                        cigar = i.cigarstring
                    #print("\t".join([cbc, umi, molid, contig_str, contig_cov, reads, contig_len, i.seq, i.get_tag('cs'), allele, cigar]) + '\n')
                    tab.write("\t".join([cbc, umi, molid, contig_str, contig_cov, reads, contig_len, i.seq, i.get_tag('cs'), allele, cigar]) + '\n')


    if not keep_tmp:
        subprocess.run(f'rm {fa_1}'.split())



def extract_10xcbc_whitelist(barcodes_tsv_ifn, output_filename = None):
    import gzip
    #"/data/junker/users/jmintch/mapped/20200812/Hr6hpiII_v4-Dr11/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

    if output_filename != None:
        with gzip.open(barcodes_tsv_ifn, 'rt') as bcs, open(output_filename, 'w') as outf:
            for i in bcs.readlines():
                    outf.write(i.strip()[0:-2:]+'\n')
    else:
        import numpy as np
        outv = []
        with gzip.open(barcodes_tsv_ifn, 'rt') as bcs:
            for i in bcs.readlines():
                    outv.append(i.strip()[0:-2:])
            return(outv)



def format_alleles_from_sam(ref_path, fa_ofn):
    # we add the right header #TODO improve
    import uuid
    from pyfaidx import Fasta
    import subprocess
    import pysam

    unique_filename = "tmpfa_"+str(uuid.uuid4())

    with open(unique_filename, 'w') as ofa, open(fa_ofn) as infile:
        FF = Fasta(ref_path)
        entries = [(f, len(FF[f][:])) for f in FF.keys()]
        for entry in entries:
            ofa.write(f'@SQ\tSN:{entry[0]}\tLN:{entry[1]}\n')
            print(entry)

        for line in infile.readlines():
            #print(line)
            ofa.write(line)

    subprocess.run(f'mv {unique_filename} {fa_ofn}'.split())
    
    with open(fa_ofn.replace('fa', 'concise'), 'w') as tab:
        tab.write("\t".join(['cbc', 'umi', 'molid', 'contig_str', 'contig_cov', 'reads', 'contig_len','cigar', 'allele', 'pruneda']) + '\n')

        for read in pysam.AlignmentFile(fa_ofn):
            start = read.reference_start
            trim_start = 1111
            trim_end = 1515
            refi = 0
            quei = 0
            #if start<trim_start:
            current = start
            allele = []
            for i in read.cigartuples:
                if i[0] == 0: #match
                    quei += i[1]
                    if len(allele) == 0:
                        allele.append(f'{current+i[1]-trim_start}M')
                    elif len(allele) == 100*len(read.cigartuples)-1:
                        allele.append(f'{trim_end -current+i[1]}M')
                    else:
                        allele.append(f'{current+i[1]}M')
                    current = current + i[1]
                if i[0] == 1: #insertion
                    allele.append(f'{current+i[1]}I')
                    current = current + i[1]
                if i[0] == 2: #deletion
                    allele.append(f'{current+i[1]}D')
                    current = current + i[1]

            cbc, umi, molid, contig_str, contig_cov, reads, contig_len = read.query_name.split('_')
            pruneda = ":".join(allele[1:-1])
            if len(pruneda) == 0 or len(allele) == 0:
                pruneda = 'NA'
            tab.write("\t".join([cbc, umi, molid, contig_str, contig_cov, reads, contig_len, read.cigarstring, ':'.join(allele), pruneda]) + '\n')
    



