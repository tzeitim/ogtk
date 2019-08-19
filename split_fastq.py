import yaml
import pdb
import regex

import sys
sys.path.append('/local/users/polivar/src/projects/ogtk')
from UM import UM as UM
import UM

def split_from_yaml(conf_fn = sys.argv[1]):
    conf_fn = sys.argv[1]
    conf = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    errors = conf['errors']


    out_prefix = conf['sample_dir']+"/"+"by_sample/"
    xp_name = conf['xp_name']

    for sample in conf['samples']:
        anchor = regex.compile(conf['samples'][sample].get("anchor", "A")+"{e<="+str(errors)+"}")
        f1_fn = conf['sample_dir'] + "/" + conf['samples'][sample]['f1_fn']
        f2_fn = conf['sample_dir'] + "/" + conf['samples'][sample]['f2_fn']
        vsamples = conf['samples'][sample]['vsamples']
        UM.split_fastq(sample, out_prefix, f1_fn, f2_fn, vsamples, anchor)
    
if __name__=="__main__":
    split_from_yaml()
