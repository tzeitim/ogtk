import rich
from pyaml import yaml

class Xp():
    ''' Imports a yaml file into an instance of the class xp (experiment). Attributes are directly assigned via the yaml file
    '''
    def __init__(self, conf_fn=None, conf_dict=None):
        if conf_fn is not None:
            conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)

        for k,v in conf_dict.items():
            setattr(self, k, v)

        setattr(self, 'conf_fn', conf_fn)
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
        xp_template = yaml.load(open(self.xp_template), Loader=yaml.FullLoader)
        setattr(self, "workdir", xp_template['workdir'])
        for k in [i for i in xp_template.keys() if i.startswith('wd_')]:
            setattr(self, k, eval(xp_template[k]))

        setattr(self, "consolidated", True)
        rich.print(":vampire:")

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

def run_cr(xp_conf):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess
    import subprocess
    import rich
    # load xp parameters
    xp = db.xp(conf_fn=xp_conf)
    # extract pre-processing config
    pp = db.xp(conf_dict=getattr(xp, 'pp'))
    # extract cell ranger config
    crg = db.xp(conf_fn=pp.cr_template)
    cr = db.xp(conf_dict=getattr(crg, pp.cr_version))
    
    rich.print(cr)
    rich.print(xp)
    subprocess.run(
        [
        cr.bin, 
        cr.cmd, 
        f'--uiport={cr.uiport}',
        f'--id={xp.sample_id}',
        f'--fastqs={xp.raw_reads}',
        f'--sample={xp.sample_id}',
        f'--transcriptome={cr.transcriptome}',
        f'--localcores={cr.localcores}',  
        f'--localmem={cr.localmem}',
        ],
        cwd=crg.scwd)# stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #cd $SC_OUTPUT_DIR 
    #/local/users/polivar/bin/cellranger-7.1.0/cellranger count \
    #--uiport=7776 \
    #--id=$CR_RUN_ID \
    #--fastqs=$FASTQ_DIR \
    #--sample=$FASTQ_SAMPLE_STR \
    #--transcriptome=$CR_TRANSCRIPTOME \
    #--localcores=$LOCAL_CORES \
    #--localmem=$LOCAL_MEM \
    #$OPTIONS

def run_bcl2fq(xp_conf, **args):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess
    import rich
    from pyaml import yaml

    # load xp parameters
    xp = Xp(conf_fn=xp_conf)
    if not xp.consolidated:
        raise ValueError("Check the configuration to be automatically consolidated")
    # extract pre-processing config
    pp = xp.pp #db.xp(conf_dict=getattr(xp, 'pp'))
    # extract bcl2fq config
    b2fq = pp['b2fq']
    # populate xp-specific information
    b2fq['outdir'] = xp.raw_reads
    #print(b2fq)
    b2fq_g = yaml.load(open(b2fq['template']), Loader=yaml.FullLoader)

    if args is not None:
        for k in args.keys():# & b2fq.keys():
            b2fq[k] = args[k]
            
    for k in b2fq.keys() & b2fq_g.keys():
        b2fq_g[k] = b2fq[k]
        print(f'setting {k}-->{b2fq_g[k]}')
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
    subprocess.run(cmd.split(), shell=False)#, stdout=open(, 'w'), stderr=open(, 'w'))
    return(cmd)
    
