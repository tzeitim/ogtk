import rich
import os
from pyaml import yaml
from rich.console import Console
from rich.text import Text
import polars as pl

class Xp():
    ''' Imports a yaml file into an instance of the class xp (experiment). Attributes are directly assigned via the yaml file
    '''
    def __init__(self, conf_fn=None, conf_dict=None, quiet=True):
        self.conf_fn = conf_fn
        self.consolidated = False
        self.console = Console(width=800)
        self.quiet = quiet

        if conf_fn is not None:
            conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)

        if conf_dict is not None:
            for k,v in conf_dict.items():
                setattr(self, k, v)

        if "xp_template" in vars(self):
            self.consolidate_conf()
    
    def __str__(self):
        return '\n'.join([f'{i}:\t{ii}' for i,ii in self.__rich_repr__()])

    def __rich_repr__(self):
        for k,v in vars(self).items():
            yield k,v
        
    def to_pl(self): 
        ''' Manually sanitize the vars() dictionary for direct coversion to a polars object
        '''
        attrs = {}
        for k,v in vars(self).items():
             attrs[k] =[v]
        return pl.DataFrame(attrs)

    def consolidate_conf(self):
        ''' The self-referencing pointers in the configuration are evaluated
        '''
        # import information from experiment template
        if self.xp_template is not None:
            xp_template = yaml.load(open(self.xp_template), Loader=yaml.FullLoader)

            # still undecided on whether to import everything or be more selective.
            # for now we import all
            #self.workdir =  xp_template['workdir']
            self.wd_datain = xp_template['datain']

            for k in xp_template.keys():
                if k.startswith('wd_'):
                    # some variables are special and need to be evaluated
                    if k not in vars(self):
                        setattr(self, k, eval(xp_template[k]))
                    else:
                        self.print(f'kept {k} from experiment conf instead of template:\n{getattr(self, k)}')
                else:
                    # direct assignment
                    if k not in vars(self):
                        setattr(self, k, xp_template[k])
                    else:
                        self.print(f'kept {k} from experiment conf instead of template:\n{getattr(self, k)}')


            setattr(self, "consolidated", True)
        else:
            print('an experiment template is needed')

    def init_workdir(self):
        ''' create dir tree for a given experiment
        '''
        if not self.consolidated:
            raise ValueError("An experiment object needs to be consolidated first, for this an `xp_template is needed`")

        for i in [i for i in vars(self) if i.startswith('wd_')]: 
            wd_dir = getattr(self, i)
            if not os.path.exists(wd_dir):
                os.system(f"mkdir -p {wd_dir}")
                self.print(f":construction:\t{wd_dir}")
            else:
                self.print(f":mag:\t{wd_dir}")

    def reset_done(self, pattern='*'):
        ''' patterns can be: "cr", "bcl2fq", ""
        '''
        cmd = f'rm {self.wd_logs}/.{pattern}_done'
        print(cmd)
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
            out_fn = f'{self.sample_id}_xpconf.yaml'
        if out_dir is None:
            out_dir = f'{self.wd_samplewd}'

        out_fn = f'{out_dir}/{out_fn}'
        self.print(f'Saving xp conf to {out_fn}')

        if xp_conf_keys is None:
            xp_conf_keys = vars(self).keys()

        ignored_keys = ['console']
        
        conf_dict = dict([(i, vars(self)[i]) for i in xp_conf_keys if i not in ignored_keys])
    
        with open(out_fn, 'w') as outfile:
            yaml.dump(conf_dict, outfile)



def return_file_name(sample_id, field):
    '''field can be:
        [ge_lib, lin_lib, tabix, rc_tabix]
    '''
    rootdir = '/local/users/polivar/src/artnilet/'
    import pandas as pd

    xps = (
        pl.read_csv('/local/users/polivar/src/artnilet/conf/xpdb_datain.txt', separator='\t')
        .filter(pl.col('sample_id')==sample_id)
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

def create_db(rootdir: str | None= None, fields = ['sample_id', 'ge_lib', 'lin_lib'])-> None:
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
        xp.print(f'setting {k}\t-->\t{cr_g[k]}', 'yellow')

    # populate command
    # sanitize types

    for k,v in cr.items():
        cr_g[k] = str(v)

    for k,v in cr_g.items():
        cr_g[k] = str(v)

    xp.print(cr)
    xp.print(cr_g)
    #cmd: BIN count --uiport=UIPORT --id=SAMPLE_ID --fastqs=FASTQ_DIR --sample=FASTQ_SAMPLE_STR --transcriptome=TRANSCRIPTOME --localcores=LOCAL_CORES --localmem=LOCAL_MEM OPTIONS
    cmd = (
            cr_g['cmd']
            .replace('BIN', cr_g['bin'])
            .replace('UIPORT', cr_g['uiport'])
            .replace('SAMPLE_ID', xp.sample_id)
            .replace('FASTQ_DIR', f'{xp.wd_fastq}/{xp.scrna_suffix}')
            .replace('FASTQ_SAMPLE_STR', f'{xp.sample_id}_{xp.scrna_suffix}')
            .replace('TRANSCRIPTOME', cr_g['transcriptome'])
            .replace('LOCAL_CORES', cr_g['localcores'])
            .replace('LOCAL_MEM', cr_g['localmem'])
            .replace('OPTIONS', cr_g['options'])
            )

    done_token = f"{xp.wd_logs}/.cr_done"

    xp.print(cmd)
    xp.print(f"{xp.wd_logs}/cr.out\n{xp.wd_logs}/cr.err", "green")
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
        xp.print(f'setting {k}\t-->\t{b2fq_g[k]}', 'yellow')

    # populate command
    # sanitize types
    for k,v in b2fq_g.items():
        b2fq_g[k] = str(v)

    xp.print(b2fq)
    xp.print(b2fq_g)

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

    xp.print(cmd)
    xp.print(f"{xp.wd_logs}/bcl2fastq.out\n{xp.wd_logs}/bcl2fastq.err", "green")

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
    
def tabulate_xp(xp, force=False):
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
    '''
    import ogtk.utils as ut

    if 'tabulate' in vars(xp):
        for suffix in xp.tabulate.keys():
            xp.print(f"{suffix}", 'bold red')

            path_to_reads = f'{xp.wd_fastq}/{suffix}'
            rev_comp_r2 = xp.tabulate[suffix]

            xp.print("path to reads:")
            xp.print(path_to_reads, 'bold white')

            pattern =f'{xp.sample_id}*R1*'
            xp.print(f"pattern={pattern}", 'bold green')
            xp.print(f"reverse complement r2 ={rev_comp_r2}", 'bold green')

            r1 = ut.sfind(path_to_reads, pattern=pattern)
            print(r1)

            if len(r1)==0:
                raise ValueError(f'No files found under the pattern {pattern}')

            if not isinstance(r1, list):
                r1 = [r1]

            for found in r1:
                xp.print(f"tabbing {found}")
                outdir = f'{xp.wd_samplewd}/{suffix}'
                ut.tabulate_paired_umified_fastqs(r1=found, 
                                                  force=force, 
                                                  rev_comp_r2=rev_comp_r2, 
                                                  export_parquet=True, 
                                                  outdir=outdir)

    else:
        raise ValueError('No "tabulate" attribute in xp. When specified, add an additional prefix field bound to a boolean variable that will determine the reverse-complementarity of read2. yaml example:\ntabulate:\n  shrna: true\n  zshrna: false\n')

def print_template(conf_fn: str = '/home/polivar/src/artnilet/conf/xps/template.yaml'):
    ''' pretty-print experiment template to ease manual completion
    '''
    conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    rich.print(conf_dict)





