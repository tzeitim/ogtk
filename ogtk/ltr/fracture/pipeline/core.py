from pathlib import Path 
from functools import wraps
from typing import Any, Callable, NamedTuple
from enum import Enum
import polars as pl

from ogtk.utils.log import Rlogger, call
from ogtk.utils import sfind, tabulate_paired_10x_fastqs_rs
from .plotting import PlotDB
from .types import FractureXp

class StepResults(NamedTuple):
    """Container for data to be passed to plotting"""
    results: dict  
    
class PipelineStep(Enum):
    """Enum defining available pipeline steps"""
    PARQUET = {
        'required_params': {'umi_len', 'rev_comp', 'modality'},
        'description': "Convert fastq reads to parquet format"
    }
    
    PREPROCESS = {
        'required_params': {'umi_len', 'anchor_ont'},
        'description': "Preprocess the loaded data"
    }
    
    FRACTURE = {
        'required_params': {'umi_len', 'anchor_ont'},
        'description': "Assemble short reads into contigs"
    }
    
    MAKE_TEST = {
        'required_params': {'target_sample'},
        'description': "Makes a downsampled fastq file for testing purposes"
    }
    
    @classmethod
    def from_string(cls, step_name: str) -> "PipelineStep":
        """Convert string to PipelineStep enum"""
        try:
            return cls[step_name.upper()]
        except KeyError:
            raise ValueError(f"Invalid step name: {step_name}. Valid steps are: {[s.name for s in cls]}")


def pipeline_step(step: PipelineStep):
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        # Pipeline instance = ppi
        def wrapper(ppi: Any, *args: Any, **kwargs: Any) -> Any:
            if not ppi.should_run_step(step):
                ppi.logger.io(f"Skipping {step.name.lower()}")
                return None
            
            try:
                # Validate parameters
                missing_params = [
                    param for param in step.value['required_params']
                    if not hasattr(ppi.xp, param)
                ]
                
                if missing_params:
                    error_msg = f"Missing required parameters for {step.name}: {', '.join(missing_params)}"
                    ppi.logger.error(error_msg)
                    raise ValueError(error_msg)

                # actually run the Step
                results = func(ppi, *args, **kwargs)
                
                if getattr(ppi.xp, 'do_plot', False):
                    try:
                        ppi.logger.step(f'Plotting {step.name.lower()} results')

                        # makes use of PlotDB to fetch the method
                        plot_method = getattr(ppi.plotdb, f"plot_{step.name.lower()}", None)

                        if plot_method:
                            ppi.logger.debug(f"Generating plots for {step.name.lower()}")
                            plot_method(ppi, results)

                    except Exception as e:
                        ppi.logger.error(f"Plot generation failed: {str(e)}")
                
                return results
            except Exception as e:
                ppi.logger.error(f"Step {step.name} failed: {str(e)}")
                raise
        return wrapper
    return decorator

class Pipeline:
    plotdb: PlotDB
    """Main pipeline class for data processing using Xp configuration"""
    dry: bool
    def __init__(self, xp: FractureXp):
        self.xp = xp
        self.logger = Rlogger().get_logger()
        self.plotdb = PlotDB()
        
        # Create output directory if it doesn't exist
        self.xp.init_workdir()

        # Initialize plotting dependencies if do_plot is True
        if getattr(self.xp, 'do_plot', False):
            try:
                import seaborn as sns
                import matplotlib.pyplot as plt
                self.sns = sns
                self.plt = plt
                self.logger.debug("Initialized plotting dependencies")
            except ImportError as e:
                self.logger.error(f"Failed to import plotting dependencies: {str(e)}")
                self.xp.do_plot = False 


    def should_run_step(self, step: PipelineStep) -> bool:
        """Check if a step should be run based on configuration"""
        if hasattr(self.xp, 'steps'):
            return step.name.lower() in [s.lower() for s in self.xp.steps] 
        return True  # Run all steps if not specified
        
    @call
    @pipeline_step(PipelineStep.PARQUET)
    def to_parquet(self) -> None:
        """Convert fastq reads to parquet format"""
        try:
            self.logger.info(f"Loading data from {self.xp.pro_datain}")
            input_files = sfind(f"{self.xp.pro_datain}", "*R1*.fastq.gz") # added R1

            
            if len(input_files) == 0:
                self.logger.error(f"No files found in {self.xp.pro_datain}")

            self.logger.io(f"found {input_files}")
            
            required_params = ['umi_len', 'rev_comp', 'modality']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")

            sample_to_file = self.xp.organize_files_by_sample(input_files, self.xp.samples, max_files=1)
            
            if not self.xp.dry:
                for sample_id, file in sample_to_file.items():
                    file = file[0]
                    sample_dir = f'{self.xp.pro_workdir}/{sample_id}/' 
                    Path(sample_dir).mkdir(parents=True, exist_ok=True)

                    tabulate_paired_10x_fastqs_rs(
                            file_path=file, 
                            out_fn=f'{sample_dir}/parsed_reads.parquet', 
                            out_dir=f'{sample_dir}',
                            modality=self.xp.modality,     
                            umi_len=self.xp.umi_len,      
                            do_rev_comp=self.xp.rev_comp,  
                            force=True) # add to step interface
            pass

        except Exception as e:
            self.logger.error(f"Failed to load data: {str(e)}")
            raise
            
    @call  
    @pipeline_step(PipelineStep.PREPROCESS)
    def preprocess(self) -> StepResults:
        """Preprocess the loaded data"""
        try:
            self.logger.step("Parsing reads and collecting molecule stats")

            required_params = ['umi_len', 'anchor_ont']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")
            in_file = f"{self.xp.pro_datain}/{self.xp.target_sample}_R1_001.merged.parquet"
            out_file = f"{self.xp.sample_wd}/parsed_reads.parquet" #pyright:ignore
            self.logger.io(f"exporting reads to {out_file}")

            if not self.xp.dry:
                (
                    pl
                    .scan_parquet(in_file)
                    .pp.parse_reads(umi_len=self.xp.umi_len, anchor_ont=self.xp.anchor_ont) #pyright: ignore
                    .collect()
                    .write_parquet(out_file)
                )
                return StepResults(
                        results={'xp': self.xp,
                                 'parsed_reads':out_file}
                        )
            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}")
            raise
            
    @call  
    @pipeline_step(PipelineStep.FRACTURE)
    def fracture(self) -> None:
        """Asseble short reads into contigs"""
            
        try:
            self.logger.step("Assemblying molecules")

            required_params = ['umi_len', 'anchor_ont']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")
            in_file = f"{self.xp.pro_datain}/{self.xp.target_sample}_R1_001.merged.parquet"
            out_file = f"{self.xp.sample_wd}/parsed_reads.parquet" #pyright:ignore
            self.logger.io(f"exporting reads to {out_file}")

            if not self.xp.dry:
                (
                    pl
                    .scan_parquet(in_file)
                    .pp.parse_reads(umi_len=self.xp.umi_len, anchor_ont=self.xp.anchor_ont) #pyright: ignore 
                    .collect()
                    .write_parquet(out_file)
                )
            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}")
            raise
            
    @call  
    @pipeline_step(PipelineStep.MAKE_TEST)
    def make_test(self) -> None:
        """Makes a downsampled fastq file for testing purposes"""
        try:
            self.logger.step("Downsampling parquet to generate synthetic FASTQ")

            input_files = sfind(f"{self.xp.pro_datain}", f"{self.xp.target_sample}*merged*parquet")
            sample_to_file = self.xp.organize_files_by_sample(files=input_files, samples=self.xp.samples, max_files=1)
            umi_len = self.xp.umi_len

            
            in_file = sample_to_file[self.xp.target_sample][0]
            
            out_file1= Path(f"{self.xp.pro_datain}/TEST_{self.xp.target_sample}_R1_001.fastq.gz").as_posix()
            out_file2= Path(f"{self.xp.pro_datain}/TEST_{self.xp.target_sample}_R2_001.fastq.gz").as_posix() 

            self.logger.io(f"exporting reads to:\n{out_file1}\n{out_file2}")

            if not self.xp.dry:



                df = (
                    pl.scan_parquet(in_file)
                    .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
                    .with_columns(pl.len().over('umi').alias('reads'))
                    .with_columns(
                        pl.col('reads')
                            .qcut(
                                quantiles  = [0.0, 0.1,0.49, 0.50, 0.52, 0.99],
                                labels = ['0','10%', '<50%', '50%', '>50%', '99%', '>99%'])
                               .alias('qreads')
                       )
                       .filter(pl.col('qreads').is_in(['10%', '50%', '99%']))
                       .filter(
                           pl.col('umi').is_in(
                               pl.col('umi').filter(pl.int_range(pl.len()).shuffle().over("qreads") < 10)
                               )
                           )
                         .collect()
           #         .write_parquet(out_file)
                )
                print(df.group_by('qreads').agg(pl.col('umi').n_unique()))

                #Schema([('read_id', String),
                #('r1_seq', String),
                #('r1_qual', String),
                #('r2_seq', String),
                #('r2_qual', String)])
    #def to_fastq(self, read_id_col: str, read_qual_col: str, read_col: str)-> pl.DataFrame:
                import gzip

                with gzip.open(out_file1, 'wb') as r1gz:
                    (
                        df.dna.to_fastq(read_id_col="read_id", read_qual_col="r1_qual", read_col='r1_seq')
                        .select('r1_seq_fastq')
                        .write_csv(r1gz, 
                                   include_header=False,
                                   )
                    )
                with gzip.open(out_file2, 'wb') as r2gz:
                    (
                        df.dna.to_fastq(read_id_col="read_id", read_qual_col="r2_qual", read_col='r2_seq')
                        .select('r2_seq_fastq')
                        .write_csv(r2gz, 
                                   include_header=False,
                                   )
                    )
            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}")
            raise

    def run(self) -> bool:
        """Run the complete pipeline"""
        if 'make_test' in self.xp.steps:
            try:
                self.logger.info("=== TEST MODE ===")
                self.make_test()

                return True

            except Exception as e:
                self.logger.error(f"Pipeline tests failed {str(e)}")
                return False

        try:
            self.to_parquet()
            self.preprocess()
            #self.analyze_data()
            #self.save_results()
            
            self.logger.info("Pipeline completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False

