
class xp():
    ''' Imports a yaml file into an instance of the class xp (experiment). Attributes are directly assigned via the yaml file
    '''
    def __init__(self, conf_fn=None, conf_dict=None):
        if conf_fn is not None:
            from pyaml import yaml
            conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)

        for k,v in conf_dict.items():
            setattr(self, k, v)

        setattr(self, 'conf_fn', conf_fn)
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
    # load xp parameters
    xp = db.xp(conf_fn=xp_conf)
    # extract pre-processing config
    pp = db.xp(conf_dict=getattr(xp, 'pp'))
    # extract bcl2fq config
    b2fq = db.xp(conf_dict=pp.b2fq)
    #print(b2fq)
    b2fq_g = db.xp(conf_fn=b2fq.template)
    args = vars(b2fq)
    if args is not None:
        for k,v in args.items():
            setattr(b2fq, k, v)
            
    for i in set(vars(b2fq).keys()).intersection(vars(b2fq_g).keys()):
        print(f'-->{i}')
    
    return(b2fq_g)
    import subprocess
    import rich
    # load xp parameters
    xp = db.xp(conf_fn=xp_conf)
    bcl = db.xp(conf_dict=xp.pp['bcl2fq'])
    # extract pre-processing config
    rich.print(bcl)
    cmd = f"conda run -n {bcl.conda_env} {bcl.cmd}"
    print(cmd)
    
    
    subprocess.run(cmd.split(), shell=False)
    return(cmd)
    cmd = f'bash -i -c "conda activate {bcl.conda_env}; {bcl.cmd}"'
    
    cmd = subprocess.run(cmd.split(), shell=True)
    
    #subprocess.run(bcl.cmd.split())
    
