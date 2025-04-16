import rich
from typing import List, Dict, TypedDict, Optional, Any
import os
from pyaml import yaml
from rich.console import Console
from rich.text import Text
import polars as pl

from functools import wraps

from ogtk.utils.log import CustomLogger, Rlogger, call
logger = Rlogger().get_logger()

class SystemConfig(TypedDict):
    prefixes: Dict[str, str]
    default: str

def init_logger(self):
    #from import Rlogger
    self.logger = Rlogger().get_logger()
    self.rlogger = Rlogger()  # Keep a reference to the Rlogger instance
    #logger.set_level("DEBUG")


def run_bcl2fq(xp, force=False, dry=False, **args):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess

    if not xp.consolidated:
        raise ValueError("Check the configuration be consolidated")

    # extract pre-processing config
    pp = xp.pp #db.xp(conf_dict=getattr(xp, 'pp'))

    # extract bcl2fq config
    b2fq = pp['b2fq']

    # populate xp-specific information
    b2fq['outdir'] = xp.wd_fastq

    #print(b2fq)
    b2fq_g = yaml.load(open(b2fq['template']), Loader=yaml.FullLoader)

    if args is not None:
        for k in args.keys():
            b2fq[k] = args[k]
            
    for k in b2fq.keys() & b2fq_g.keys():
        b2fq_g[k] = b2fq[k]
        logger.debug(f'setting {k}\t-->\t{b2fq_g[k]}')

    # populate command
    # sanitize types
    for k,v in b2fq_g.items():
        b2fq_g[k] = str(v)

    logger.debug(b2fq)
    logger.debug(b2fq_g)

    #bcl2fastq  --processing-threads THREADS --no-lane-splitting --barcode-mismatches BARCODEMM --sample-sheet SAMPLE_SHEET --runfolder-dir RUNDIR --output-dir OUTDIR OPTIONS'}
    cmd = (
        b2fq_g['cmd']
            .replace('THREADS', b2fq_g['threads'])
            .replace('BARCODEMM', b2fq_g['barcode_mismatches'])
            .replace('SAMPLE_SHEET', b2fq_g['samplesheet'])
            .replace('RUNDIR', b2fq_g['rundir'])
            .replace('OUTDIR', b2fq_g['outdir'])
            .replace('OPTIONS', b2fq_g['options'])
        )
    cmd = f"conda run -n {b2fq_g['conda_env']} {cmd}"
    done_token = f"{xp.wd_logs}/.bcl2fq_done"

    logger.debug(cmd)
    logger.debug(f"{xp.wd_logs}/bcl2fastq.out\n{xp.wd_logs}/bcl2fastq.err")

    if not dry:
        if os.path.exists(done_token) and not force:
            return(0)

        p1 = subprocess.run(cmd.split(),
                            shell=False,
                            stdout=open(f'{xp.wd_logs}/bcl2fastq.out', 'w'),
                            stderr=open(f'{xp.wd_logs}/bcl2fastq.err', 'w'))

        if p1.returncode == 0:
            subprocess.run(f'touch {done_token}'.split())

        return(p1)
    else:
       return(0)
import os
import subprocess
import hashlib



# changes here might break the invokation in the pipeline since the will be missing arguments
@call
def tabulate_xp(xp, modality, cbc_len, umi_len, force=False)->List:
    ''' 
    Tabulates paired fastq of umified reads (R1:UMI:CB, R2:RNA) into the
    parquet format. The xp configuration requires a `tabulate` field which in
    turn needs a list of prefixes bound to a boolean variable that will
    determine whether read2 is reverse-complemented.
    For example (yaml format):
    ```
     tabulate:
       shrna: true
       zshrna: false
    ```

    Returns: list of paths of tabulated files
    '''
    import ogtk.utils as ut

    if 'tabulate' in vars(xp):
        for suffix in xp.valid_tab_suffixes():
            logger.debug(f"{suffix}")

            path_to_reads = f'{xp.wd_fastq}/{suffix}'
            rev_comp_r2 = xp.tabulate[suffix]

            logger.debug("path to reads:")
            logger.debug(path_to_reads)

            pattern =f'{xp.target_sample}*R1*.fastq.gz'
            logger.debug(f"pattern={pattern}")
            logger.debug(f"reverse complement r2 ={rev_comp_r2}")


            r1_input_files = [
                            i for i in 
                                ut.sfind(path_to_reads, pattern = "*_R1_*fastq.gz") 
                            if not i.split('/')[-1].startswith("Undetermined") 
                            ]

            logger.debug(r1_input_files)

            if len(r1_input_files)==0:
                raise ValueError(f'No files found under the pattern {pattern}')

            if not isinstance(r1_input_files, list):
                r1_input_files = [r1_input_files]

            tabulated_files = []
            for found in r1_input_files:
                logger.debug(f"tabbing {found}")
                outdir = f'{xp.wd_xp}/{suffix}'
                out_fn = f"{outdir}/{found.split('/')[-1]}".replace('.fastq.gz', '.mols.parquet')
                logger.io(out_fn)

                if not os.path.exists(out_fn) or force:
                    ut.tabulate_paired_10x_fastqs_rs(
                        file_path=found,
                        cbc_len=cbc_len,
                        umi_len=umi_len,
                        modality=modality,
                        out_fn=out_fn,
                        force=force,
                        do_rev_comp=rev_comp_r2,
                        )
                else:
                    logger.io(f"loading from cache:\n{out_fn}")

                tabulated_files.append(out_fn)
            return tabulated_files

    else:
        raise ValueError('No "tabulate" attribute in xp. When specified, add an additional prefix field bound to a boolean variable that will determine the reverse-complementarity of read2. yaml example:\ntabulate:\n  shrna: true\n  zshrna: false\n')

def print_template(conf_fn: str = '/home/polivar/src/artnilet/conf/xps/template.yaml'):
    ''' pretty-print experiment template to ease manual completion
    '''
    conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    rich.print(conf_dict)


class Xp():
    ''' Imports a yaml file into an instance of the class xp (experiment). Attributes are directly assigned via the yaml file
    '''
    # Required attributes
    system: SystemConfig
    xp_datain: str
    xp_template: str
    prefix: str
    consolidated: bool
    special_patterns: List[str]
    
    # Optional attributes that may come from template or config
    project: Optional[str]
    target_sample: Optional[str]
    gene_conf_fn: Optional[str]
    steps: Optional[List[str]]
    
    # Internal attributes
    console: Console
    quiet: bool
    logger: CustomLogger
    rlogger: Any 
    conf_fn: Optional[str]

    def __init__(self, conf_fn=None, conf_dict=None, quiet=True):
        self.conf_fn = conf_fn
        self.consolidated = False
        self.console = Console(width=800)
        self.quiet = quiet

        self.logger = Rlogger().get_logger()
        self.rlogger = Rlogger()  # Keep a reference to the Rlogger instance

        # populate the conf dir via a file or directly from an argument
        if conf_fn is not None:
            conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
        
        if conf_dict is not None:
            for k,v in conf_dict.items():
                setattr(self, k, v)
       
        # resolve prefix
        self._resolve_system_prefix()

        # resolve type of experiment using xp_template
        self.consolidate_conf()

        if 'gene_conf_fn' in vars(self):
            self.__init_genes()
    
    def __str__(self):
        return '\n'.join([f'{i}:\t{ii}' for i,ii in self.__rich_repr__()])

    def __rich_repr__(self):
        for k,v in vars(self).items():
            yield k,v

    def __init_genes(self):
        ''' import a pre-defined set of genes and gene_conf
        '''
        conf_dict = yaml.load(open(self.gene_conf_fn), Loader=yaml.FullLoader) #pyright: ignore
        for k,v in conf_dict.items():
            setattr(self, k, v)
    
    def _resolve_system_prefix(self):
        '''
        '''
        # check if prefix is in the environment
        env_prefix = os.environ.get('OGTK_SYSPREFIX')
        if env_prefix:
            self.logger.info(f"Getting system prefix from environment:\n{env_prefix}")
            setattr(self, 'prefix', env_prefix)
            return 


        # if not found as environ var get it from the conf 
        if "system" not in vars(self):
            raise ValueError(self._prefix_help())

        if "prefixes" not in self.system and "default" not in self.system:
            raise ValueError(self._prefix_help())

        prefix = self.system['prefixes'][self.system['default']]
        self.logger.info(f"Using system prefix from config file:\n{prefix}")

        setattr(self, 'prefix', prefix)


    def _populate_special(self, dic: Dict, var_pref: str, update=False):
        ''' '''
        for k,v in dic.items():
            if k.startswith(var_pref):
                if k not in vars(self) or update: 
                    setattr(self, k, eval(v))
                else:
                    logger.debug(f'kept {k} from experiment conf instead of template:\n{getattr(self, k)}')

    @staticmethod
    def _prefix_help():
        return "A system prefix configuration needs to be provided in the form \n"  \
                            "system:\n "\
                            " prefixes:\n"\
                            "    host1: '/home/user/path/to/project/'\n"\
                            "    host2: '/mount/server/user/path/to/project/'\n"\
                            "    laptop: '/Volumes/mount/user/'\n"\
                            " default: 'host1'  # Default prefix to use\n"
                            
    def logger_set_level(self, level):
        '''
        '''
        self.rlogger.set_level(level)

    def to_pl(self): 
        ''' Manually sanitize the vars() dictionary for direct coversion to a polars object
        '''
        attrs = {}
        for k,v in vars(self).items():
             attrs[k] =[v]
        return pl.DataFrame(attrs)

    def consolidate_conf(self, update=False):
        ''' The self-referencing pointers in the configuration are evaluated '''

        if not hasattr(self, 'xp_template'):
            raise ValueError("An experiment template must be provided")
        else:
            # import information from experiment template

            self.xp_template = self.xp_template.replace("${prefix}", self.prefix)
            xp_template = yaml.load(open(self.xp_template), Loader=yaml.FullLoader) #pyright:ignore

            if "special_patterns" not in xp_template:
                raise ValueError("special_patterns must be an attribute of an xp_template. e.g.:\nspecial_patterns: [xp_, pro_, sample_]")

            # some variables are special and need to be evaluated
            # following the hierachy: xp_ -> pro_ -> sample_
            self.special_patterns = xp_template['special_patterns']
            for var_pref in self.special_patterns:
                self._populate_special(xp_template, var_pref, update)

            # match special patterns to variables
            self.special_vars = []
            for k,v in vars(self).items():
                for pattern in self.special_patterns:
                    if pattern in k:
                        self.special_vars.append(k)

            # direct assignment of non-special variables
            for k,v in xp_template.items():
                if k not in vars(self) and k not in self.special_vars:
                    setattr(self, k, v)
                else:
                    logger.debug(f'kept {k} from experiment conf instead of template:\n{getattr(self, k)}')
            
            if 'tabulate' in vars(self):
                assert isinstance(self.tabulate, dict), "The .tabulate attribute of an expriment must be a dictionary" #pyright:ignore

                for suffix in self.tabulate.keys(): #pyright:ignore
                    allowed_suffixes = [v for k,v in xp_template.items() if k.endswith('_suffix')]

                    if suffix in allowed_suffixes:
                        xp_template[f'wd_{suffix}'] = "f'{self.wd_xp}/{suffix}'"
                    else:
                        logger.critical("The provided tabulation suffixes do not match the experiment template")

                        
            setattr(self, "consolidated", True)

    def init_workdir(self):
        ''' create dir tree for a given experiment
        '''
        if not self.consolidated:
            raise ValueError("An experiment object needs to be consolidated first, for this an `xp_template is needed`")

        for i in [i for i in vars(self) if i in self.special_vars]: 
            wd_dir = getattr(self, i)
            if not wd_dir:
                raise ValueError(f"Empty directory path for special variable '{i}'. Special variables must have non-empty values for directory creation.")
                
            if not os.path.exists(wd_dir):
                os.system(f"mkdir -p {wd_dir}")
                logger.debug(f":construction:\t{wd_dir}", extra={"markup": True})
            else:
                logger.debug(f":mag:\t{wd_dir}", extra={"markup": True})

    def valid_tab_suffixes(self)->List|None:
        if 'tabulate' in vars(self):
            return [k for k,v in self.tabulate.items()]
        else:
            return None

    def reset_done(self, pattern='*'):
        ''' patterns can be: "cr", "bcl2fq", ""
        '''
        cmd = f'rm {self.wd_logs}/.{pattern}_done'
        logger.debug(cmd)
        os.system(cmd)

    def print(self, text, style="bold magenta", force=False, *args, **kwargs):
        '''
        '''
        #text = Text(text)
        #text.stylize(style)
        if not self.quiet or force:
            self.console.print(text, style=style, *args, **kwargs)

    def export_xpconf(self, xp_conf_keys = None, out_fn = None, out_dir=None):
        ''' Saves current instance of an experiment to the sample_wd directory as default
        '''
        if out_fn is None:
            out_fn = f'{self.target_sample}_xpconf.yaml'
        if out_dir is None:
            out_dir = f'{self.wd_xp}'

        out_fn = f'{out_dir}/{out_fn}'
        logger.io(f'Saving xp conf to {out_fn}')

        if xp_conf_keys is None:
            xp_conf_keys = vars(self).keys()

        ignored_keys = ['console']
        
        conf_dict = dict([(i, vars(self)[i]) for i in xp_conf_keys if i not in ignored_keys])
    
        with open(out_fn, 'w') as outfile:
            yaml.dump(conf_dict, outfile)

    @call
    @wraps(run_bcl2fq)
    def demux(self, *args, **kwargs):
        ''' demultiplex
        '''
        run_bcl2fq(self, *args, **kwargs)


def return_file_name(target_sample, field):
    '''field can be:
        [ge_lib, lin_lib, tabix, rc_tabix]
    '''
    rootdir = '/local/users/polivar/src/artnilet/'
    import pandas as pd

    xps = (
        pl.read_csv('/local/users/polivar/src/artnilet/conf/xpdb_datain.txt', separator='\t')
        .filter(pl.col('target_sample')==target_sample)
    )

    lin_lib = xps['lin_lib'].to_list()[0]
    if 'tabix' in field:
        value = xps['lin_lib'].str.replace(f'_R1.+.fastq.gz', '.sorted.txt.gz').to_list()[0]
        #value = lin_lib.replace(f'S.{"{1,}"}_R1.+.fastq.gz', '.sorted.txt.gz')
        if 'rc' in field:
            value = value.replace('sorted', 'rc_sorted')

    return_value = f'{rootdir}/datain/{value}'
    if not os.path.exists(return_value):
        print(f'There was an issue while trying to open file {return_value}')
    return(return_value)


def load_db(rootdir: str='/local/users/polivar/src/artnilet')-> None:
    rich.print(f'loaded {rootdir}')

def create_db(rootdir: str | None= None, fields = ['target_sample', 'ge_lib', 'lin_lib'])-> None:
    '''
    '''
    if rootdir is None:
        raise ValueError("A root dir must be provided") 

    if not os.path.exists(rootdir):
        raise ValueError("path doesn't exist. A root directory must be manually created")

    else:
        print("not implemented")

def run_cranger(xp, force=False, dry=False, **args):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess

    if not xp.consolidated:
        raise ValueError("Check the configuration be consolidated")

    # extract pre-processing config
    pp = xp.pp #db.xp(conf_dict=getattr(xp, 'pp'))

    # extract bcl2fq config
    cr = pp['cr']

    # populate xp-specific information
    cr['outdir'] = xp.wd_fastq

    #print(b2fq)
    cr_g = yaml.load(open(cr['template']), Loader=yaml.FullLoader)
    cr_g = cr_g[cr['version']]

    if args is not None:
        for k in args.keys():
            cr[k] = args[k]
            
    for k in cr.keys() & cr_g.keys():
        cr_g[k] = cr[k]
        logger.debug(f'setting {k}\t-->\t{cr_g[k]}')

    # populate command
    # sanitize types

    for k,v in cr.items():
        cr_g[k] = str(v)

    for k,v in cr_g.items():
        cr_g[k] = str(v)

    logger.debug(cr)
    logger.debug(cr_g)
    #cmd: BIN count --uiport=UIPORT --id=SAMPLE_ID --fastqs=FASTQ_DIR --sample=FASTQ_SAMPLE_STR --transcriptome=TRANSCRIPTOME --localcores=LOCAL_CORES --localmem=LOCAL_MEM OPTIONS
    cmd = (
            cr_g['cmd']
            .replace('BIN', cr_g['bin'])
            .replace('UIPORT', cr_g['uiport'])
            .replace('SAMPLE_ID', xp.target_sample)
            .replace('FASTQ_DIR', f'{xp.wd_fastq}/{xp.scrna_suffix}')
            .replace('FASTQ_SAMPLE_STR', f'{xp.target_sample}_{xp.scrna_suffix}')
            .replace('TRANSCRIPTOME', cr_g['transcriptome'])
            .replace('LOCAL_CORES', cr_g['localcores'])
            .replace('LOCAL_MEM', cr_g['localmem'])
            .replace('OPTIONS', cr_g['options'])
            )

    done_token = f"{xp.wd_logs}/.cr_done"

    logger.debug(cmd)
    logger.debug(f"{xp.wd_logs}/cr.out\n{xp.wd_logs}/cr.err")
    # TODO add a cr entry for the cmd ran

    if not dry:
        if os.path.exists(done_token) and not force:
            return(0)

        p1 = subprocess.run(cmd.split(), 
                            shell=False, 
                            stdout=open(f'{xp.wd_logs}/cr.out', 'w'), 
                            stderr=open(f'{xp.wd_logs}/cr.err', 'w'), 
                            cwd=xp.wd_scrna)

        if p1.returncode == 0:
            subprocess.run(f'touch {done_token}'.split())
        else:
            print(f"something went wrong have a look at:\n{xp.wd_logs}/cr.out\n{xp.wd_logs}/cr.err")

        return(p1)
    else:
        return(0)

