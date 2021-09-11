
def mltb_bin_alleles(name, intid, intid2_R2_strand, config_card_dir, outdir, fq1, fq2, threads = 100, end = 5000, ranked_min_cov = 2, umi_len = 17, umi_errors = 1, alg_gapopen = 20, alg_gapextend =1):
    '''Once fastq files have been split by: RPI/rev_id/intid process them into a table of barcodes for downstream analysis on R; ranked_min_cov determines the threshold to call a representative rank1 sequence as valid'''
    intid2 = ogtk.UM.rev_comp(intid2_R2_strand)# same strand as initd1
    # TODO consider generating a conf file out of the pre-processing step
    # constant vars
    ref_path = outdir+'/hspdrv7_scgstl_{}.fa'.format(name)
    rcref_path = outdir+'/rc_hspdrv7_scgstl_{}.fa'.format(name)
    ref_name = 'hspdrv7_scgstl_{}'.format(name)
    #ref_seq = 'gcgccaccacctgttcctgtagaaat{}ccggactcagatctcgagctcaagcttcggacagcagtatcatcgactaTGGagtcgagagcgcgctcgtcgactaTGGagtcgtcagcagtactactgacgaTGGagtcgacagcagtgtgtgagtctaTGGagtcgagagcatagacatcgagtaTGGagtcgactacagtcgctacgactaTGGagtcgacagagatatcatgcagtaTGGagtcgacagcagtatctgctgtcaTGGagtcgactgcacgacagtcgactaTGGAGTCG{}cgagcgctatgagcgactatgc'.format(intid, intid2).upper()
    ref_seq = 'gccaccacctgttcctgtagaaat{}ccggactcagatctcgagctcaagcttcggacagcagtatcatcgactaTGGagtcgagagcgcgctcgtcgactaTGGagtcgtcagcagtactactgacgaTGGagtcgacagcagtgtgtgagtctaTGGagtcgagagcatagacatcgagtaTGGagtcgactacagtcgctacgactaTGGagtcgacagagatatcatgcagtaTGGagtcgacagcagtatctgctgtcaTGGagtcgactgcacgacagtcgactaTGGAGTCG{}cgagcgctatgagcgactatgc'.format(intid, intid2).upper()
    bint_db_ifn = '/local/users/polivar/src/projects/mltracer/conf/gstlt_barcode_intervsx10.yaml'

    yaml_out_fn = config_card_dir + '/config_card_{}.yaml'.format(name)
    yaml_out = {'name':name, 'desc':{}, 'intid1':intid, 'intid2':intid2, 'modality':'bulk'}
    
    # for intid in dataset
    # merge
    fqm =       outdir + '/{}_merged.fastq'.format(name)
    fqum1 =     outdir + '/{}_unmerged_R1.fastq'.format(name)
    fqum2 =     outdir + '/{}_unmerged_R2.fastq'.format(name)
    ihist =     outdir + '/{}_ihist.txt'.format(name)
    log_merge = outdir + '/{}_merge_stats.txt'.format(name)
    cmd_merge = "bbmerge-auto.sh in1={} in2={} out={} outu={} outu2={} ihist={} ecct extend2=20 iterations=5 interleaved=false -Xmx2g ".format(fq1, fq2, fqm, fqum1, fqum2, ihist)
    
    yaml_out['desc']['read_merge'] = "Parameters used for the merging of read 1 and read 2"

    umi_counts_ofn = outdir + '/{}_umi_counts.txt'.format(name)
    yaml_out['mols'] = {}
    yaml_out['mols']['umi_len'] = umi_len
    yaml_out['mols']['umi_correction'] = umi_errors
    yaml_out['mols']['counts'] = umi_counts_ofn 
    
    yaml_out['read_merge'] = {}
    yaml_out['read_merge']['fqm'] = fqm
    yaml_out['read_merge']['fqum1'] = fqum1
    yaml_out['read_merge']['fqum2'] = fqum2
    yaml_out['read_merge']['ihist'] = ihist
    yaml_out['read_merge']['log_merge'] = log_merge

    # align
    fa_correctedm = outdir+'/{}_corrected_merged.fa'.format(name)
    fa_corrected1 = outdir+'/{}_corrected_R1.fa'.format(name)
    fa_corrected2 = outdir+'/{}_corrected_R2.fa'.format(name)

    yaml_out['desc']['alignment'] = "Parameters for the pairwise alignment of every recovered allele and the reference. The corrected files fix the issue of duplicated fasta entry names"

    yaml_out['alignment'] = {}
    yaml_out['alignment']['fa_correctedm'] = fa_correctedm
    yaml_out['alignment']['fa_corrected1'] = fa_corrected1
    yaml_out['alignment']['fa_corrected2'] = fa_corrected2

    # more outfiles
    merged_tab_out =   outdir + '/{}_barcodes_merged.txt'.format(name)
 
    yaml_out['desc']['lineage'] = "Main output of the pre-processing scripts. Tabulated files that contain the bint barcode string"
    yaml_out['lineage'] = {}
    yaml_out['lineage']['merged_tab'] = merged_tab_out
    yaml_out['lineage']['merged_full'] = merged_tab_out.replace('.txt','_full.txt')

   
    # merge paired end reads if possible

    pp = subprocess.run(cmd_merge.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    fout = open(log_merge, 'wb')
    fout.write(pp.stdout)
    fout.write(pp.stderr)
    fout.close()

    if end != None:
        print("Warning: end is not None; You are not processing all the data", end)

    # do a pass on the raw fastqs and group reads by UMI
    ## for merged (rsm) and unmerged (rs1, rs2)
    rsm =           ogtk.UM.fastq_collapse_UMI(fqm, umi_len = umi_len, end = end, keep_rid = True)
    rs1, rs2 =      ogtk.UM.pfastq_collapse_UMI(fqum1, fqum2, umi_len = umi_len, end = end)

    # remove umis that exist on the merged set out of the unmerged
    unmerged_by_error = set(rsm.umi_list()).intersection(set(rs1.umi_list()))
    rs1.delete(unmerged_by_error)
    rs2.delete(unmerged_by_error)

    #TODO the saturation stats should be outsourced somehow and reported for all three cases not just the merged (rsm) 
    yaml_out['mols']['saturation'] = {}
    yaml_out['mols']['saturation']['merged'] = {}
    yaml_out['mols']['saturation']['unmerged'] = {}
    yaml_out['mols']['nreads'] = {}
    yaml_out['mols']['nreads']['total'] = sum([rsm.nreads, rs1.nreads])
    yaml_out['mols']['nreads']['merged'] = rsm.nreads
    yaml_out['mols']['nreads']['umerged'] = rs1.nreads
    yaml_out['mols']['desc'] = "number of umis whose rank1 sequence is a >= [1, 2, 3, 4, 5, 10, 100, 1000]"
    yaml_out['mols']['saturation']['merged']['uncorrected'] = rsm.allele_saturation()
    yaml_out['mols']['saturation']['unmerged']['uncorrected1'] = rs1.allele_saturation()
    yaml_out['mols']['saturation']['unmerged']['uncorrected2'] = rs2.allele_saturation()

    if umi_errors >0:
        print("Correcting umis with a hdist of {}".format(umi_errors))
        rsm.correct_umis(errors = umi_errors , silent = True, jobs = 80)
        rs1.correct_umis(errors = umi_errors , silent = True, jobs = 80)
        rs2.correct_umis(errors = umi_errors , silent = True, jobs = 80)
        yaml_out['mols']['saturation']['merged']['corrected'] = rsm.saturation()
        yaml_out['mols']['saturation']['unmerged']['corrected1'] = rs1.saturation()
        yaml_out['mols']['saturation']['unmerged']['corrected2'] = rs2.saturation()

    else:
        print("Not correcting umis")

    
    rsm.consensus = ogtk.UM.do_fastq_pileup(rsm, min_cov = ranked_min_cov, threads = threads, trim_by = None)
    rs1.consensus = ogtk.UM.do_fastq_pileup(rs1, min_cov = ranked_min_cov, threads = threads, trim_by = None)
    rs2.consensus = ogtk.UM.do_fastq_pileup(rs2, min_cov = ranked_min_cov, threads = threads, trim_by = None)

    # align allele to reference. merged and unmerged separately
    ## make sure to include the intid on the reference to help the alignment
    mlt_create_fasta_ref(ref_name, ref_seq, ref_path)
    ## also on rev comp (not sure if this is really needed)
    mlt_create_fasta_ref(ref_name+"RC", ogtk.UM.rev_comp(ref_seq), rcref_path)

    alg_mode = 'needleman'
    #alg_gapopen = 20
    #alg_gapextend = 1

    yaml_out['alignment']['alg_mode'] = alg_mode
    yaml_out['alignment']['alg_gapopen'] = alg_gapopen
    yaml_out['alignment']['alg_gapextend'] = alg_gapextend


    if len(rsm.consensus.keys()) == 0 :
        yaml_out['success'] = False
        return(None)
    else:
        yaml_out['success'] = True
 

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"m", fa_ofn = fa_correctedm, 
                            consensus_dict = rsm.consensus, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"1", fa_ofn = fa_corrected1, 
                            consensus_dict = rs1.consensus, ref_path = ref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    mltbc_align_reads_to_ref(name = outdir + "/"+name+"_"+"2", fa_ofn = fa_corrected2, 
                            consensus_dict = rs2.consensus, ref_path = rcref_path, 
                            ref_name = ref_name, mode=alg_mode, gapopen = alg_gapopen, 
                            gapextend = alg_gapextend)

    
    mlt_compute_barcode_matrix_merged(fa_ifn = fa_correctedm, tab_out= merged_tab_out, bint_db_ifn = bint_db_ifn, do_rc=False)

    rs1c = list(rs1.consensus.keys())
    rs2c = list(rs2.consensus.keys())
    umi_whitelist = list(set(rs1c).intersection(set(rs2c)))

    print("whitelist:", len(umi_whitelist))

    if len(umi_whitelist)>0:
        unmerged_tab_out = outdir + '/{}_barcodes_unmerged.txt'.format(name)
        yaml_out['lineage']['paired_tab'] = unmerged_tab_out 
        yaml_out['lineage']['paired_full'] = unmerged_tab_out.replace('.txt','_full.txt')

        mlt_compute_barcode_matrix_paired(fa_ifn = fa_corrected1, fa_ifn2 = fa_corrected2, umi_whitelist = umi_whitelist, 
                                           umi_blacklist = list(rsm.consensus.keys()), tab_out = unmerged_tab_out, bint_db_ifn = bint_db_ifn, do_rc = True)

    # save the count information for the UMIs 
    with open(umi_counts_ofn, 'w') as fh:
        fh.write(f'count\tumi\n')
        for i in (rsm, rs1):
            counts = i.umi_counts()
            for count,umi in zip(counts, counts.index):
                fh.write(f'{count}\t{umi}\n')
    
    # save run settings 
    yaml_stream = open(yaml_out_fn, 'w')
    yaml.dump(yaml_out, yaml_stream)
    yaml_stream.close() 
    print("Run settings saved to config card:\n{}".format(yaml_out_fn))

 
