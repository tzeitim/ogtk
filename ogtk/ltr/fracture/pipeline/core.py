from pathlib import Path 
from functools import wraps
from typing import Any, Callable, NamedTuple, Dict, Set
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
    
    TEST = {
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


def log_invocation_params(logger: Any, step: PipelineStep, xp: Any, kwargs: Dict = None) -> None:
    """
    Log the invocation parameters for a pipeline step to both the logger and a file.
    
    Args:
        logger: Logger instance to use for logging
        step: The pipeline step being executed
        xp: The experiment configuration object
        kwargs: Additional parameters passed to the step function
    """
    from datetime import datetime
    
    # Helper function to recursively extract nested parameters
    def extract_nested_params(obj, prefix=""):
        params = {}
        if isinstance(obj, dict):
            for k, v in obj.items():
                key = f"{prefix}.{k}" if prefix else k
                if isinstance(v, (dict, list)):
                    params.update(extract_nested_params(v, key))
                else:
                    params[key] = v
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                key = f"{prefix}[{i}]"
                if isinstance(item, (dict, list)):
                    params.update(extract_nested_params(item, key))
                else:
                    params[key] = item
        return params
    
    # Log required parameters to logger
    logger.info(f"=== {step.name} Parameters ===")
    for param in step.value['required_params']:
        if hasattr(xp, param):
            param_value = getattr(xp, param)
            logger.info(f"  {param}: {param_value}")
    
    # Log step-specific nested configurations
    step_specific_configs = {
        PipelineStep.FRACTURE: ['fracture'],
        PipelineStep.PARQUET: ['force_tab'],
        # Add other step-specific params as needed
    }
    
    if step in step_specific_configs:
        for config_attr in step_specific_configs[step]:
            if hasattr(xp, config_attr):
                config_value = getattr(xp, config_attr)
                if isinstance(config_value, dict):
                    logger.info(f"  {config_attr} configuration:")
                    for k, v in config_value.items():
                        logger.info(f"    {k}: {v}")
                else:
                    logger.info(f"  {config_attr}: {config_value}")
    
    # Log additional parameters from kwargs
    if kwargs:
        logger.info("  Additional parameters:")
        for k, v in kwargs.items():
            logger.info(f"    {k}: {v}")
    
    # Record parameters to a separate file for easy retrieval
    step_log_dir = Path(f"{xp.pro_workdir}/{xp.target_sample}/logs")
    # Create logs directory if it doesn't exist
    step_log_dir.mkdir(parents=True, exist_ok=True)
    param_file = step_log_dir / f"{step.name.lower()}_params.log"
    
    with open(param_file, 'w') as f:
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        f.write(f"=== {step.name} Parameters - {timestamp} ===\n")
        
        # Write required parameters
        for param in step.value['required_params']:
            if hasattr(xp, param):
                param_value = getattr(xp, param)
                f.write(f"{param}: {param_value}\n")
        
        # Write step-specific nested configurations
        if step in step_specific_configs:
            f.write(f"\n=== {step.name} Specific Configuration ===\n")
            for config_attr in step_specific_configs[step]:
                if hasattr(xp, config_attr):
                    config_value = getattr(xp, config_attr)
                    if isinstance(config_value, (dict, list)):
                        nested_params = extract_nested_params(config_value, config_attr)
                        for k, v in nested_params.items():
                            f.write(f"{k}: {v}\n")
                    else:
                        f.write(f"{config_attr}: {config_value}\n")
        
        # Include other important configuration parameters
        important_params = ['target_sample', 'pro_workdir', 'pro_datain']
        f.write("\n=== Additional Configuration ===\n")
        for param in important_params:
            if hasattr(xp, param):
                param_value = getattr(xp, param)
                f.write(f"{param}: {param_value}\n")
        
        if kwargs:
            f.write("\n=== Function Arguments ===\n")
            for k, v in kwargs.items():
                f.write(f"{k}: {v}\n")


def pipeline_step(step: PipelineStep):
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(ppi: Any, *args: Any, **kwargs: Any) -> Any:
            if not ppi.should_run_step(step):
                ppi.logger.io(f"Skipping {step.name.lower()}")
                return None
            
            try:
                # Set up step-specific logging
                step_log_dir = Path(f"{ppi.xp.pro_workdir}/{ppi.xp.target_sample}/logs")
                # Create logs directory if it doesn't exist
                step_log_dir.mkdir(parents=True, exist_ok=True)
                step_log_file = step_log_dir / f"{step.name.lower()}.log"
                
                # Enable step-specific file logging
                logger = Rlogger()
                logger.enable_file_logging(
                    filepath=step_log_file,
                    level=ppi.xp.log_level if hasattr(ppi.xp, 'log_level') else "INFO",
                    format_string="%(asctime)s - %(levelname)s - %(message)s",
                    mode='w'  # Start fresh log for each step run
                )
                
                # Validate parameters
                missing_params = [
                    param for param in step.value['required_params']
                    if not hasattr(ppi.xp, param)
                ]
                
                if missing_params:
                    error_msg = f"Missing required parameters for {step.name}: {', '.join(missing_params)}"
                    ppi.logger.error(error_msg)
                    raise ValueError(error_msg)

                # Run the step
                ppi.logger.info(f"Starting {step.name} step")
                
                # Log invocation parameters
                log_invocation_params(ppi.logger, step, ppi.xp, kwargs)
                
                # Execute the step function
                results = func(ppi, *args, **kwargs)
                
                if getattr(ppi.xp, 'do_plot', False):
                    try:
                        ppi.logger.step(f'Plotting {step.name.lower()} results')
                        plot_method = getattr(ppi.plotdb, f"plot_{step.name.lower()}", None)
                        
                        if plot_method:
                            ppi.logger.debug(f"Generating plots for {step.name.lower()}")
                            plot_method(ppi, results)
                    except Exception as e:
                        ppi.logger.error(f"Plot generation failed: {str(e)}")
                
                ppi.logger.info(f"Completed {step.name} step")
                
                # Disable step-specific logging
                logger.disable_file_logging()
                
                return results
                
            except Exception as e:
                ppi.logger.error(f"Step {step.name} failed: {str(e)}")
                # Ensure we disable file logging even if there's an error
                logger.disable_file_logging()
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
            self.logger.info(f"Conversion to parquet:\nloading data from {self.xp.pro_datain}")
            input_files = sfind(f"{self.xp.pro_datain}", f"{self.xp.target_sample}*R1*.fastq.gz") # added R1
            if len(input_files) == 0:
                raise ValueError(f"No files found in {self.xp.pro_datain}")

            self.logger.io(f"found {input_files}")
            
            required_params = ['umi_len', 'rev_comp', 'modality']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")

            sample_to_file = self.xp.organize_files_by_sample(input_files, [{'id':self.xp.target_sample}], max_files=1)
            
            # TODO cache
            if not self.xp.dry:
                for sample_id, file in sample_to_file.items():
                    file = file[0]
                    sample_dir = f'{self.xp.pro_workdir}/{sample_id}' 
                    Path(sample_dir).mkdir(parents=True, exist_ok=True)
                    if hasattr(self.xp, "force_tab"):
                        force_tab = self.xp.force_tab
                    else:
                        force_tab = False

                    tabulate_paired_10x_fastqs_rs(
                            file_path=file, 
                            out_fn=f'{sample_dir}/parsed_reads.parquet', 
                            out_dir=f'{sample_dir}',
                            merged_fn = f'{sample_dir}/merged_reads.parquet',
                            modality=self.xp.modality,     
                            umi_len=self.xp.umi_len,      
                            do_rev_comp=self.xp.rev_comp,  
                            force=force_tab) # add to step interface
            pass

        except Exception as e:
            self.logger.error(f"Failed to convert {self.xp.target_sample} to parquet: {str(e)}")
            raise
            
    @call  
    @pipeline_step(PipelineStep.PREPROCESS)
    def preprocess(self) -> StepResults|None:
        """Preprocess the loaded data"""
        try:
            self.logger.step("Parsing reads and collecting molecule stats")

            required_params = ['umi_len', 'anchor_ont']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")

            in_file = f"{self.xp.pro_workdir}/{self.xp.target_sample}/merged_reads.parquet"
            out_file = f"{self.xp.sample_wd}/parsed_reads.parquet" #pyright:ignore
            self.logger.io(f"reading from {in_file}")

            if not self.xp.dry:
                (
                    pl
                    .scan_parquet(in_file)
                    .pp.parse_reads(umi_len=self.xp.umi_len,  
                                    anchor_ont=self.xp.anchor_ont,
                                    intbc_5prime=self.xp.intbc_5prime)
                    .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                    .collect()
                    .write_parquet(out_file)
                )
                self.logger.info(f"exported parsed reads to {out_file}")

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
    def fracture(self) -> StepResults|None:
        """Asseble short reads into contigs"""
            
        try:
            self.logger.step("Assemblying molecules")

            required_params = ['umi_len', 'anchor_ont']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")
            in_file = f"{self.xp.sample_wd}/parsed_reads.parquet" #pyright:ignore

            if not hasattr(self.xp, 'fracture'):
                self.logger.warning("using hardcoded arguments for assembly as no configuration was found.")
                setattr(self.xp, 'fracture', {})
                self.xp.fracture['start_min_coverage'] = 25
                self.xp.fracture['start_k'] = 17
                self.xp.fracture['min_reads'] = 100 

            if not self.xp.dry:
                self.logger.info(f'Reading {in_file}')
                for priority in [True, False]:
                    out_file = f"{self.xp.sample_wd}/contigs_pl_{priority}.parquet" #pyright:ignore
                    self.logger.info(f"exporting assembled contigs to {out_file}")
                    
                    chi = (
                        pl.read_parquet(in_file)
                        .pp.assembly_with_opt( #pyright: ignore
                            start_k=self.xp.fracture['start_k'], 
                            start_min_coverage=self.xp.fracture['start_min_coverage'],
                            min_reads=self.xp.fracture['min_reads'], 
                            prioritize_length=priority,
                            method=self.xp.fracture['assembly_method'],
                            start_anchor=self.xp.start_anchor,
                            end_anchor=self.xp.end_anchor,
                            )
                        .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                    )
                    chi.write_parquet(out_file)

                    self.logger.critical(f"{(chi.get_column('length')==0).mean():.2%} failed")
                
                return StepResults(
                        results={'xp': self.xp,
                                 'contigs':out_file}
                        )

            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}")
            raise
            
    @call  
    @pipeline_step(PipelineStep.TEST)
    def test(self) -> StepResults|None:
        """Makes a downsampled fastq file for testing purposes"""
        try:
            import os
            self.logger.step("Downsampling parquet to generate synthetic FASTQ")

            #input file is the merged parquet from the original FASTQ
            original_sample = self.xp.target_sample.replace('TEST_', '') #pyright: ignore
            pattern = f"{original_sample}*merged*parquet"
            input_files = sfind(f"{self.xp.pro_workdir}/{original_sample}/", pattern)

            if len(input_files)==0:
                raise ValueError(f"\nNo files matched:\n'{pattern}'\nin\n{self.xp.pro_datain}")

            sample_to_file = self.xp.organize_files_by_sample(input_files, [{'id':original_sample}], max_files=1)
            umi_len = self.xp.umi_len
            
            # with sample_to_file we get the _merged_parquet of the original sample
            # the parsed file is not useful since we will recreate the fastq 
            in_file = sample_to_file[original_sample][0]
            
            out_file= Path(f"{self.xp.sample_wd}/{self.xp.target_sample}.test.parquet").as_posix()

            out_file1= Path(f"{self.xp.pro_datain}/{self.xp.target_sample}_R1_001.fastq.gz").as_posix()
            out_file2= Path(f"{self.xp.pro_datain}/{self.xp.target_sample}_R2_001.fastq.gz").as_posix() 

            self.logger.io(f"exporting reads to:\n{out_file1}\n{out_file2}")

            if os.path.exists(out_file1) and os.path.exists(out_file2):
                self.logger.info(f"Using previous test files use --clean to remove all generated files")
                return StepResults(
                        results={'xp': self.xp}
                        )
                pass

            umis_per_group = 30
            if not self.xp.dry:
                df = (
                    pl.scan_parquet(in_file)
                    .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
                    .with_columns(pl.len().over('umi').alias('reads'))
                    .with_columns(
                        pl.col('reads')
                            .qcut(
                                quantiles  = [0.0, 0.1,0.49, 0.50, 0.52, 0.90],
                                allow_duplicates=True,
                                labels = ['0','10%', '<50%', '50%', '>50%', '90%', '>90%'])
                               .alias('qreads')
                       )
                       .filter(pl.col('qreads').is_in(['10%', '50%', '90%']))
                       .filter(
                           pl.col('umi').is_in(
                               pl.col('umi').filter(pl.int_range(pl.len()).shuffle().over("qreads") <= umis_per_group)
                               )
                           )
                         .collect()
                )
                #QC
                #print(df.group_by('qreads').agg(pl.col('umi').n_unique()))
                df.write_parquet(out_file)
                import gzip

                with gzip.open(out_file1, 'wb') as r1gz:
                    (
                        df.dna.to_fastq(read_id_col="read_id", read_qual_col="r1_qual", read_col='r1_seq')
                        .select('r1_seq_fastq')
                        .write_csv(r1gz, 
                                   quote_style='never',
                                   include_header=False,
                                   )
                    )
                    # TODO: there is a bug here since we are not rev_complementing back the sequences/nor reversing the quality
                with gzip.open(out_file2, 'wb') as r2gz:
                    (
                        df
                        .with_columns(pl.col('read_id').str.replace("1:N:0:", "2:N:0:"))
                        .dna.to_fastq(read_id_col="read_id", read_qual_col="r2_qual", read_col='r2_seq')
                        .select('r2_seq_fastq')
                        .write_csv(r2gz, 
                                   quote_style='never',
                                   include_header=False,
                                   )
                    )
                
                return StepResults(
                        results={'xp': self.xp}
                        )
            pass
        except Exception as e:
            self.logger.error(f"Failed generating tests: {str(e)}")
            raise

    def clean_test_outputs(self) -> None:
        """Remove test-related outputs from previous runs"""
        try:
            import shutil

            self.logger.step("Cleaning test outputs")
            
            # Clean test fastq files
            test_pattern = "TEST_*_R[12]_001.fastq.gz"
            test_files = Path(self.xp.pro_datain).glob(test_pattern)
            for file in test_files:
                self.logger.io(f"Removing {file}")
                file.unlink(missing_ok=True)
                
            # Clean test working directories
            if hasattr(self.xp, 'samples'):
                for sample in self.xp.samples:
                    if str(sample['id']).startswith('TEST_'):
                        test_dir = Path(f"{self.xp.pro_workdir}/{sample['id']}")
                        if test_dir.exists():
                            self.logger.io(f"Removing directory {test_dir}")
                            shutil.rmtree(test_dir)

            self.logger.info("Test outputs cleaned successfully")
            
        except Exception as e:
            self.logger.error(f"Failed to clean test outputs: {str(e)}")
            raise

    def clean_all(self) -> None:
        """Remove all pipeline outputs from previous runs"""
        try:
            import shutil

            self.logger.step("Cleaning all pipeline outputs")
            
            # First clean test outputs
            self.clean_test_outputs()
            
            # Clean working directories for all samples
            if hasattr(self.xp, 'samples'):
                for sample in self.xp.samples:
                    sample_dir = Path(f"{self.xp.pro_workdir}/{sample['id']}")
                    if sample_dir.exists():
                        self.logger.io(f"Removing directory {sample_dir}")
                        shutil.rmtree(sample_dir)
                    
            self.logger.info("All pipeline outputs cleaned successfully")
            
        except Exception as e:
            self.logger.error(f"Failed to clean pipeline outputs: {str(e)}")
            raise

    def run(self) -> bool:
        """Run the complete pipeline"""
        if self.xp.make_test:
            try:
                self.logger.info("=== TEST MODE ===")
                self.xp.target_sample = f'TEST_{self.xp.target_sample}'
                self.xp.samples.append({'id':self.xp.target_sample})
                self.xp.consolidate_conf(update=True)
                self.xp.init_workdir()

                self.logger.info(f"Created test alias {self.xp.target_sample} and initialized workdir")
                self.test()
                self.run_all()

                self.logger.info("Pipeline (TEST) completed successfully")
                return True

            except Exception as e:
                self.logger.error(f"Pipeline (TEST) tests failed {str(e)}")
                return False

        try:
            self.run_all()
            
            self.logger.info("Pipeline completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False

    def run_all(self):
        ''' Run all steps'''
        self.logger.info("==== Running pipeline ====")

        self.to_parquet()
        self.preprocess()
        self.fracture()
        #self.save_results()
