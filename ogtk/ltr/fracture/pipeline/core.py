from logging import lastResort
from pathlib import Path 
from functools import wraps
from typing import Any, Callable, NamedTuple, Dict, Optional
from enum import Enum
import polars as pl

from ogtk.utils.log import Rlogger, call
from ogtk.utils import sfind, tabulate_paired_10x_fastqs_rs
from .plotting import PlotDB
from .types import FractureXp
from . import qc


class StepResults(NamedTuple):
    #TODO is it necessary so many levels? (results, metrics)
    """Container for data to be passed to plotting or QCs
        - results: is meant to contain elements or objects related to the pipeline progression
        - metrics: keeps track of values that are used to generate the final summary and compute stats 
    """
    results: dict  
    metrics: dict = {}

    def has_metrics(self):
        """Check if a step returned any metrics """
        return len(self.metrics)>0

    
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
    logger.debug(f"=== {step.name} Parameters ===")
    for param in step.value['required_params']:
        if hasattr(xp, param):
            param_value = getattr(xp, param)
            logger.debug(f"  {param}: {param_value}")
    
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
                # why is this not the same logger as the one in the class??
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
                    ppi.logger.error(error_msg, with_traceback=True)
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
                        ppi.logger.error(f"Plot generation failed: {str(e)}", with_traceback=True)

                # run QCs relevant for the step
                ppi.logger.step(f'Running QCS {step.name}')
                
                ppi.logger.info(f"Completed {step.name} step")
                ppi.update_pipeline_summary(step, results)
                logger.disable_file_logging()
                
                return results
                
            except Exception as e:
                ppi.logger.error(f"Step {step.name} failed: {str(e)}", with_traceback=True)
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
                self.logger.error(f"Failed to import plotting dependencies: {str(e)}", with_traceback=True)
                self.xp.do_plot = False 


    def should_run_step(self, step: PipelineStep) -> bool:
        """Check if a step should be run based on configuration"""
        if hasattr(self.xp, 'steps'):
            return step.name.lower() in [s.lower() for s in self.xp.steps] 
        return True  # Run all steps if not specified
        
    @call
    @pipeline_step(PipelineStep.PARQUET)
    def to_parquet(self) -> StepResults|None:
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

            max_files = self.xp.allow_wildcards if hasattr(self.xp, "allow_wildcards") else 1  
            sample_to_file = self.xp.organize_files_by_sample(
                                input_files, 
                                [{'id':self.xp.target_sample}], 
                                max_files=max_files
                            )
            
            # TODO cache
            if not self.xp.dry:
                for sample_id, files in sample_to_file.items():
                    sample_dir = f'{self.xp.pro_workdir}/{sample_id}' 
                    Path(sample_dir).mkdir(parents=True, exist_ok=True)
                    force_tab = getattr(self.xp, "force_tab", False)
                    out_fn=f'{sample_dir}/parsed_reads.parquet'
                    merged_fn = f'{sample_dir}/merged_reads.parquet'
                    
                    # TODO the Rust implementation doesn't seem to respect the limit argument
                    limit = getattr(self.xp, "limit",  None)

                    if limit is not None:
                        self.logger.warning(f"Only loading {limit} reads")

                    # Check if multiple files should be processed
                    if len(files) > 1 and max_files > 1:
                        self.logger.info(f"Processing {len(files)} input files for sample {sample_id}")
                        temp_parsed_files = []
                        temp_merged_files = []
                        
                        for i, file in enumerate(files):
                            temp_out_fn = f'{sample_dir}/temp_parsed_reads_{i}.parquet'
                            temp_merged_fn = f'{sample_dir}/temp_merged_reads_{i}.parquet'
                            
                            self.logger.info(f"Processing file {i+1}/{len(files)}: {file}")
                            tabulate_paired_10x_fastqs_rs(
                                file_path=file, 
                                out_fn=temp_out_fn,
                                merged_fn=temp_merged_fn,
                                out_dir=sample_dir,
                                modality=self.xp.modality,     
                                umi_len=self.xp.umi_len,      
                                do_rev_comp=self.xp.rev_comp,  
                                limit=limit,
                                force=force_tab)
                            
                            temp_parsed_files.append(temp_out_fn)
                            temp_merged_files.append(temp_merged_fn)
                        
                        self.logger.info(f"Concatenating results from {len(files)} files")
                        (pl.concat([pl.scan_parquet(f) for f in temp_parsed_files]).sink_parquet(out_fn))
                        (pl.concat([pl.scan_parquet(f) for f in temp_merged_files]).sink_parquet(merged_fn))

                        #TODO if limit is set, downsample the concatenation too
                        
                        for f in temp_parsed_files + temp_merged_files:
                            Path(f).unlink(missing_ok=True)
                            
                        self.logger.info(f"Successfully concatenated results")
                    else:
                        file = files[0]
                        self.logger.info(f"Processing single file: {file}")
                        tabulate_paired_10x_fastqs_rs(
                            file_path=file, 
                            out_fn=out_fn,
                            merged_fn=merged_fn,
                            out_dir=sample_dir,
                            modality=self.xp.modality,     
                            umi_len=self.xp.umi_len,      
                            do_rev_comp=self.xp.rev_comp,  
                            force=force_tab)# expose to cli?

                # TODO why are merged files and parsed files so similar? can we get read of merged?
                return StepResults(
                        results={
                        'parsed_fn': out_fn,
                        },
                        metrics = {
                        'total_reads': (
                                pl.scan_parquet(out_fn)
                                .collect()
                                .height),
                        'total_umis': (
                                pl.scan_parquet(out_fn)
                                .select('umi')
                                .unique()
                                .collect()
                                .height),
                        'downsampled': limit
                        }
                        )

            pass

        except Exception as e:
            self.logger.error(f"Failed to convert {self.xp.target_sample} to parquet: {str(e)}", with_traceback=True)
            raise
            
    @call  
    @pipeline_step(PipelineStep.PREPROCESS)
    def preprocess(self) -> StepResults|None:
        """Preprocess the loaded data"""
        try:
            self.logger.step("Parsing reads and collecting molecule stats")

            required_params = ['umi_len', 'anchor_ont', 'sbc_len']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")

            in_file = f"{self.xp.pro_workdir}/{self.xp.target_sample}/merged_reads.parquet"
            out_file = f"{self.xp.sample_wd}/parsed_reads.parquet" #pyright:ignore
            out_file_qc = out_file.replace('.parquet', '_qc.parquet')
            self.logger.io(f"reading from {in_file}")

            if not self.xp.dry:
                ldf = (
                    pl
                    .scan_parquet(in_file)
                    .pp.parse_reads(umi_len=self.xp.umi_len,  
                                    sbc_len=self.xp.sbc_len,
                                    anchor_ont=self.xp.anchor_ont,
                                    )
                    .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                )

                # save without gc fields
                if not self.xp.parse_read1:
                    ldf.drop('^_.+$').sink_parquet(out_file)
                else:
                    # if r1 is long and informative: e.g long paired-end
                    self.logger.warning(f"extracting juice from R1")
                    (
                        pl.concat(
                            [
                                ldf, 
                                ldf.pp.parse_read1(anchor_ont=self.xp.anchor_ont)
                            ]
                        ).sink_parquet(out_file)
                     )

                # extract QC fields
                metrics_df = ldf.drop('^[^_].+$').unique().collect()
                
                metrics = {i[1:]:metrics_df[i][0] for i in metrics_df.columns}

                (
                    metrics_df.
                    rename({i:i[1:] for i in metrics_df.columns})
                    .write_parquet(out_file_qc)
                )

                self.logger.info(f"exported parsed reads to {out_file}")
                self.logger.info(f"exported parsed QC metrics to {out_file_qc}")

                return StepResults(
                        results={'xp': self.xp,
                                 'parsed_reads':out_file,
                                 'preprocess_qc': out_file_qc,
                                 },
                        metrics=metrics,
                        )
            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}", with_traceback=True)
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
                setattr(self.xp, 'fracture', {})
                self.xp.fracture['start_min_coverage'] = 25
                self.xp.fracture['start_k'] = 17
                self.xp.fracture['min_reads'] = 10 #qc.find_read_count_threshold(in_file, method="kmeans")
            else:
                self.xp.fracture['min_reads'] = 10 #qc.find_read_count_threshold(in_file, 
                                                   #              method=self.xp.fracture['method_th'])


            if not self.xp.dry:
                self.logger.info(f'Reading {in_file}')
                out_file = f"{self.xp.sample_wd}/contigs_pl_direct.parquet" #pyright:ignore
                self.logger.info(f"exporting assembled contigs to {out_file}")
                
                filter_expr= pl.col('reads')>=self.xp.fracture['min_reads']

                df_contigs = (
                        pl.read_parquet(in_file)
                        .filter(filter_expr)
                        .pp.assemble_umis( #pyright: ignore
                          k=self.xp.fracture['start_k'], 
                          min_coverage=self.xp.fracture['start_min_coverage'],
                          start_anchor=self.xp.start_anchor,
                          end_anchor=self.xp.end_anchor,
                          min_length=None,
                          auto_k=False,
                          export_graphs=False,
                          only_largest=True,
                          method=self.xp.fracture['assembly_method'],
                          )
                        .with_columns(pl.col('contig').str.len_chars().alias('length'))
                        .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                        .join(
                            pl.scan_parquet(in_file)
                              .select('umi','reads')
                              .filter(filter_expr)
                              .unique()
                              .collect(), 
                           left_on='umi', right_on='umi', how='left')
                        )

                df_contigs.write_parquet(out_file)

                success_rate = (df_contigs.get_column('length')>0).mean()
                total_assembled = (df_contigs.get_column('length')>0).sum()
               
                self.logger.critical(f"{total_assembled} assembled ({success_rate:.2%} success)")
                
                return StepResults(
                        results={'xp': self.xp, 
                                 'contigs':out_file},
                        metrics={'success_rate':success_rate,
                                 'total_assembled':total_assembled,
                                 'mean_contig_length': df_contigs.filter(pl.col('length')>0).get_column('length').mean(),
                                 'median_contig_length': df_contigs.filter(pl.col('length')>0).get_column('length').median(),
                                 },
                        )

            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}", with_traceback=True)
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
                        results={'xp': self.xp},
                        metrics={}
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
                        results={'xp': self.xp},
                        )
            pass
        except Exception as e:
            self.logger.error(f"Failed generating tests: {str(e)}", with_traceback=True)
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
            self.logger.error(f"Failed to clean test outputs: {str(e)}", with_traceback=True)
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
            self.logger.error(f"Failed to clean pipeline outputs: {str(e)}", with_traceback=True)
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
                self.logger.error(f"Pipeline (TEST) tests failed {str(e)}", with_traceback=True)
                return False

        try:
            self.run_all()
            
            self.logger.info("Pipeline completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}", with_traceback=True)
            return False

    def run_all(self):
        ''' Run all steps'''
        self.logger.info("==== Running pipeline ====")

        self.to_parquet()
        self.preprocess()
        self.fracture()
        #self.run_qcs()

    def update_pipeline_summary(self, step: PipelineStep, results: StepResults) -> None:
        """Update the pipeline summary with results from a specific step."""
        import json
        import datetime
        from pathlib import Path
        
        # TODO include at the same level of the steps in the JSON file a
        # summary of the sample, like name, date, input material

        # Path to the summary file
        summary_path = Path(f"{self.xp.pro_workdir}/{self.xp.target_sample}/pipeline_summary.json")
        
        # Load existing summary if available
        if summary_path.exists():
            try:
                with open(summary_path, 'r') as f:
                    summary = json.load(f)
            except json.JSONDecodeError:
                summary = {}
        else:
            summary = {}
        
        # Get timestamp
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        # Create step summary
        step_summary = {
            'sample_id': self.xp.target_sample,
            'timestamp': timestamp,
            'params_file': f"logs/{step.name.lower()}_params.log",
            'log_file': f"logs/{step.name.lower()}.log",
        }
        
        # Add step-specific metrics
        if results and results.results and results.has_metrics():
            metrics = results.metrics
            if step == PipelineStep.PARQUET:
                step_summary['metrics'] = metrics
            elif step == PipelineStep.PREPROCESS:
                step_summary['metrics'] = metrics
            elif step == PipelineStep.FRACTURE:
                step_summary['metrics'] = metrics
        
        # Add figures information
        figure_mapping = {
            PipelineStep.PARQUET: [],
            PipelineStep.PREPROCESS: [
                f"{self.xp.target_sample}_coverage.png",
                f"{self.xp.target_sample}_sat-coverage.png",
                f"{self.xp.target_sample}_anchors.png",
                f"{self.xp.target_sample}_feasible.png"
            ],
            PipelineStep.FRACTURE: []
        }
        
        if step in figure_mapping and hasattr(self.xp, 'sample_figs'):
            figure_files = figure_mapping[step]
            if figure_files:
                # Check which figures actually exist
                figures_dir = Path(self.xp.sample_figs)
                if figures_dir.exists():
                    existing_figures = []
                    for figure in figure_files:
                        figure_path = figures_dir / figure
                        if figure_path.exists():
                            existing_figures.append(f"../figures/{figure}")
                            
                    if existing_figures:
                        step_summary['figures'] = existing_figures
        
        # Check for output files
        if results and results.results:
            output_files = []
            if step == PipelineStep.PARQUET and 'parsed_fn' in results.results:
                output_files.append(results.results['parsed_fn'])
            elif step == PipelineStep.PREPROCESS and 'parsed_reads' in results.results:
                output_files.append(results.results['parsed_reads'])
            elif step == PipelineStep.FRACTURE and 'contigs' in results.results:
                output_files.append(results.results['contigs'])
                
            if output_files:
                step_summary['output_files'] = output_files
                
        # Update summary
        summary[step.name.lower()] = step_summary
        
        # Write updated summary
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Generate markdown report
        self._generate_markdown_report(summary)
    def _generate_markdown_report(self, summary: dict) -> None:
        """
        Generate a markdown report from the summary JSON.
        Includes a summary table of all metrics across all steps and figures generated during the pipeline.
        """
        from pathlib import Path
        import datetime
        
        report_path = Path(f"{self.xp.pro_workdir}/{self.xp.target_sample}/pipeline_summary.md")
        
        with open(report_path, 'w') as f:
            # Generate report header
            current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            f.write(f"# Pipeline Summary for {self.xp.target_sample}\n\n")
            f.write(f"*Last updated: {current_time}*\n\n")
            
            # Create a consolidated metrics table
            f.write("## Consolidated Metrics\n\n")
            
            # Collect all unique metrics across all steps
            all_metrics = set()
            for step_info in summary.values():
                if 'metrics' in step_info:
                    all_metrics.update(step_info['metrics'].keys())
            
            # Sort metrics for consistent display
            all_metrics = sorted(all_metrics)
            
            if all_metrics:
                # Generate table header
                f.write("| Metric | " + " | ".join(step.capitalize() for step in summary.keys()) + " |\n")
                f.write("|--------|" + "|".join(["------" for _ in summary]) + "|\n")
                
                # Generate table rows
                for metric in all_metrics:
                    metric_display = metric.replace('_', ' ').title()
                    row = f"| {metric_display} |"
                    
                    for step in summary.keys():
                        step_info = summary[step]
                        if 'metrics' in step_info and metric in step_info['metrics']:
                            value = step_info['metrics'][metric]
                            # Format numeric values with commas and round floats
                            if isinstance(value, int):
                                formatted_value = f"{value:,}"
                            elif isinstance(value, float):
                                formatted_value = f"{value:.2f}"
                            else:
                                formatted_value = str(value)
                            row += f" {formatted_value} |"
                        else:
                            row += " - |"
                    
                    f.write(row + "\n")
                
                f.write("\n")
            else:
                f.write("*No metrics available yet*\n\n")
            
            # Timeline of step execution
            f.write("## Pipeline Timeline\n\n")
            f.write("| Step | Timestamp | Status |\n")
            f.write("|------|-----------|--------|\n")
            
            for step_name, step_info in summary.items():
                timestamp = step_info.get('timestamp', 'Unknown')
                
                # Determine status based on metrics
                if 'metrics' in step_info and step_info['metrics']:
                    status = "✅ Completed"
                else:
                    status = "⚠️ Incomplete"
                    
                f.write(f"| {step_name.capitalize()} | {timestamp} | {status} |\n")
            
            f.write("\n")
            
            # Add figures section at the top level
            f.write("## Generated Figures\n\n")
            figure_path = Path(f"{self.xp.sample_figs}")
            
            if figure_path.exists():
                figures = list(figure_path.glob(f"{self.xp.target_sample}*.png"))
                if figures:
                    # Group figures by category based on filename patterns
                    figure_categories = {
                        "Coverage Analysis": ["coverage.png"],
                        "Saturation Coverage": ["sat-coverage.png"],
                        "Anchor Analysis": ["anchors.png", "feasible.png"],
                        "Threshold Detection": ["kmeans", "kneedle"],
                        "Other": []
                    }
                    
                    for category, patterns in figure_categories.items():
                        category_figures = []
                        for fig in figures:
                            if any(pattern in fig.name for pattern in patterns):
                                category_figures.append(fig)
                                
                        if category_figures:
                            f.write(f"### {category}\n\n")
                            for fig in category_figures:
                                rel_path = Path("../figures") / fig.name
                                f.write(f"![{fig.stem}]({rel_path})\n\n")
                            
                            f.write("\n")
                else:
                    f.write("*No figures generated yet*\n\n")
            else:
                f.write("*Figures directory not found*\n\n")
            
            # Detailed section for each step
            f.write("## Step Details\n\n")
            
            for step_name, step_info in summary.items():
                f.write(f"### {step_name.capitalize()}\n\n")
                f.write(f"**Run Time:** {step_info.get('timestamp', 'Unknown')}\n\n")
                
                # Show links to log files
                if 'params_file' in step_info:
                    f.write(f"**Parameter File:** [{step_name}_params.log]({step_info['params_file']})\n\n")
                
                if 'log_file' in step_info:
                    f.write(f"**Log File:** [{step_name}.log]({step_info['log_file']})\n\n")
                
                # Show step-specific figures from summary
                if 'figures' in step_info and step_info['figures']:
                    f.write("**Figures:**\n\n")
                    for fig_path in step_info['figures']:
                        fig_name = Path(fig_path).name
                        f.write(f"- [{fig_name}]({fig_path})\n")
                    f.write("\n")
                
                # Show output files
                if 'output_files' in step_info and step_info['output_files']:
                    f.write("**Output Files:**\n\n")
                    for file_path in step_info['output_files']:
                        file_name = Path(file_path).name
                        f.write(f"- {file_name}\n")
                    f.write("\n")
                
                # Show detailed metrics
                if 'metrics' in step_info and step_info['metrics']:
                    f.write("**Metrics:**\n\n")
                    f.write("| Metric | Value |\n")
                    f.write("|--------|-------|\n")
                    for metric, value in step_info['metrics'].items():
                        metric_display = metric.replace('_', ' ').title()
                        
                        # Format numeric values
                        if isinstance(value, int):
                            formatted_value = f"{value:,}"
                        elif isinstance(value, float):
                            formatted_value = f"{value:.2f}"
                        else:
                            formatted_value = str(value)
                            
                        f.write(f"| {metric_display} | {formatted_value} |\n")
                
                f.write("\n---\n\n")

