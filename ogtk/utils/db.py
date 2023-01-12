import rich
from pyaml import yaml

class Xp():
    ''' Imports a yaml file into an instance of the class xp (experiment). Attributes are directly assigned via the yaml file
    '''
    def __init__(self, conf_fn=None, conf_dict=None):
        self.conf_fn = conf_fn
        self.consolidated = False

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
        import polars as pl
        attrs = {}
        for k,v in vars(self).items():
             attrs[k] =[v]
        return pl.DataFrame(attrs)

    def consolidate_conf(self):
        ''' The self-referencing pointers in the configuration are evaluated
        '''
        if self.xp_template is not None:
            xp_template = yaml.load(open(self.xp_template), Loader=yaml.FullLoader)

            self.workdir =  xp_template['workdir']
            self.wd_datain = xp_template['datain']

            for k in [i for i in xp_template.keys() if i.startswith('wd_')]:
                setattr(self, k, eval(xp_template[k]))

            setattr(self, "consolidated", True)
            rich.print(":vampire:")
        else:
            print('an experiment template is needed')

    def init_workdir(self):
        ''' create dir tree for a given experiment
        '''
        import os
        if not self.consolidated:
            raise ValueError("An experiment object needs to be consolidated first, for this an `xp_template is needed`")

        for i in [i for i in vars(self) if i.startswith('wd_')]: 
            wd_dir = getattr(self, i)
            if not os.path.exists(wd_dir):
                os.system(f"mkdir -p {wd_dir}")
                print(f"created\t{wd_dir}")
            else:
                print(f"found\t{wd_dir}")

    def reset_done(self, pattern='*'):
        '''
        '''
        import os
        cmd = f'rm {self.wd_logs}/.{pattern}_done'
        print(cmd)
        os.system(cmd)


def return_file_name(sample_id, field):
    '''field can be:
        [ge_lib, lin_lib, tabix, rc_tabix]
    '''
    rootdir = '/local/users/polivar/src/artnilet/'
    import polars as pl
    import pandas as pd
    import os

    xps = (
        pl.read_csv('/local/users/polivar/src/artnilet/conf/xpdb_datain.txt', sep='\t')
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
    import rich
    import polars as pl

    rich.print(f'loaded {rootdir}')

def create_db(rootdir: str | None= None, fields = ['sample_id', 'ge_lib', 'lin_lib'])-> None:
    '''
    '''
    import os
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
    import rich
    import os
    from pyaml import yaml

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
        print(f'setting {k}-->\t{cr_g[k]}')

    # populate command
    # sanitize types

    for k,v in cr.items():
        cr_g[k] = str(v)

    for k,v in cr_g.items():
        cr_g[k] = str(v)

    rich.print(cr)
    rich.print(cr_g)
    #cmd: BIN count --uiport=UIPORT --id=SAMPLE_ID --fastqs=FASTQ_DIR --sample=FASTQ_SAMPLE_STR --transcriptome=TRANSCRIPTOME --localcores=LOCAL_CORES --localmem=LOCAL_MEM OPTIONS
    cmd = (
            cr_g['cmd']
            .replace('BIN', cr_g['bin'])
            .replace('UIPORT', cr_g['uiport'])
            .replace('SAMPLE_ID', xp.sample_id)
            .replace('FASTQ_DIR', f'{xp.wd_fastq}/ge')
            .replace('FASTQ_SAMPLE_STR', xp.sample_id+'_ge')
            .replace('TRANSCRIPTOME', cr_g['transcriptome'])
            .replace('LOCAL_CORES', cr_g['localcores'])
            .replace('LOCAL_MEM', cr_g['localmem'])
            .replace('OPTIONS', cr_g['options'])
            )

    done_token = f"{xp.wd_logs}/.cr_done"

    if not dry:
        if os.path.exists(done_token) and not force:
            return(0)

        p1 = subprocess.run(cmd.split(), shell=False, stdout=open(f'{xp.wd_logs}/cr.out', 'w'), stderr=open(f'{xp.wd_logs}/cr.err', 'w'), cwd=xp.wd_scrna)

        if p1.returncode == 0:
            subprocess.run(f'touch {done_token}'.split())
        else:
            print(f"something went wrong have a look at:\n{xp.wd_logs}/cr.out\n{xp.wd_logs}/cr.err")

        return(p1)
    else:
        return(cmd)

def run_bcl2fq(xp, force=False, dry=False, **args):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess
    import rich
    import os
    from pyaml import yaml

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
        print(f'setting {k}-->\t{b2fq_g[k]}')

    # populate command
    # sanitize types
    for k,v in b2fq_g.items():
        b2fq_g[k] = str(v)

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

    if not dry:
        if os.path.exists(done_token) and not force:
            return(0)

        p1 = subprocess.run(cmd.split(), shell=False, stdout=open(f'{xp.wd_logs}/bcl2fastq.out', 'w'), stderr=open(f'{xp.wd_logs}/bcl2fastq.err', 'w'))

        if p1.returncode == 0:
            subprocess.run(f'touch {done_token}'.split())

        return(p1)
    else:
        return(cmd)
    
