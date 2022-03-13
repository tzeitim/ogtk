import yaml
import pdb
import regex

import sys
#sys.path.append('/local/users/polivar/src/projects/ogtk')
from ogtk import UM as UM
#import UM

def split_from_yaml(conf_fn = sys.argv[1]):
    conf_fn = sys.argv[1]
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    errors = conf['errors']
    keep_id = conf.get('keep_id', False)
    out_prefix = conf['sample_dir']+"/"+"by_sample/"
    xp_name = conf['xp_name']

    for sample in conf['samples']:
        anchor = regex.compile(conf['samples'][sample].get("anchor", "A")+"{e<="+str(errors)+"}")
        f1_fn = conf['sample_dir'] + "/" + conf['samples'][sample]['f1_fn']
        f2_fn = conf['sample_dir'] + "/" + conf['samples'][sample]['f2_fn']
        vsamples = conf['samples'][sample]['vsamples']
        UM.split_fastq(sample_name = sample, out_prefix = out_prefix, f1_fn = f1_fn, f2_fn = f2_fn, vsamples = vsamples, anchor = anchor, keep_id = keep_id)
    
if __name__=="__main__":
    import sys, argparse

    use_txt = """    This script is able to demultiplex fastq files by using two main
    elements. 1) an anchor sequence on R1 and 2) a sample id (specified through a
    yaml conf file).  The anchor supports fuzzy matching (errors) to account for
    small sequencing errors. Index samples are exact.

    An example of a yaml conf file follows:

            xp_name: nope_nanog
            sample_dir: /local/users/polivar/tmp/RM_20191128/data
            errors: 3
            samples:
             64-cells_nanog:
              f1_fn: 64_cell_S5_R1_001.fastq
              f2_fn: 64_cell_S5_R2_001.fastq
              anchor: \"AGCAGTCGAGA\"
              vsamples:
               - TAGATC
               - CTCTCT
               - TATCCT
               - AGAGTA
               - GTAAGG
               - ACTGCA
"""
    parser=argparse.ArgumentParser(description=use_txt, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('conf', type=str, help='path to YAML config file')
    #parser.add_argument('command', type=str, choices=['all', 'map', 'call'], help='defines the main mode of the script')
    #parser.add_argument('-c', '--conf', type=str, help='[map] path to YAML config file')
    #parser.add_argument('-w', '--wildcard', type=str, help='use this parameter to pass any str; used in development')
    #parser.add_argument('-m', '--mincov', type=int, help='override config file mincov value', default=None)
    #parser.add_argument('-n', '--number', type=int, help='[map] subsample how many reads from original fastq', default=None)
    args=parser.parse_args()
 
    split_from_yaml(conf_fn = args.conf)
