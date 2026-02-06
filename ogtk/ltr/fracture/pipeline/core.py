from pathlib import Path 
from functools import wraps
from typing import Any, Callable, Dict
from enum import Enum
import polars as pl

from ogtk.utils.log import Rlogger, call
from ogtk.utils.general import sfind, tabulate_paired_10x_fastqs_rs
from .plotting import PlotDB
from .types import FractureXp,StepResults

class PipelineStep(Enum):
    """Enum defining available pipeline steps"""

    DORADO = {
        'required_params': {'target_sample'},
        'description': "Run dorado basecaller on POD5 files"
    }
    
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
        PipelineStep.DORADO: ['pp'],
        PipelineStep.FRACTURE: ['fracture'],
        PipelineStep.PARQUET: ['force_tab'],
        # Add other step-specific params as needed
    }
    
    if step in step_specific_configs:
        for config_attr in step_specific_configs[step]:
            if hasattr(xp, config_attr):
                config_value = getattr(xp, config_attr)
                if isinstance(config_value, dict):
                    logger.debug(f"  {config_attr} configuration:")
                    for k, v in config_value.items():
                        logger.debug(f"    {k}: {v}")
                else:
                    logger.debug(f"  {config_attr}: {config_value}")
    
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
                
                # Validate pipeline parameters
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
                
                if isinstance(results, StepResults):
                    results = [results]

                for result in results:
                    if getattr(ppi.xp, 'do_plot', False):
                        try:
                            ppi.logger.step(f'Plotting {step.name.lower()} results')
                            plot_method = getattr(ppi.plotdb, f"plot_{step.name.lower()}", None)
                            
                            if plot_method:
                                ppi.logger.debug(f"Generating plots for {step.name.lower()}")
                                plot_method(ppi, result)
                        except Exception as e:
                            ppi.logger.error(f"Plot generation failed: {str(e)}", with_traceback=True)


                # run QCs relevant for the step
                ppi.logger.step(f'Running QCS {step.name}')
                
                ppi.logger.info(f"Completed {step.name} step")
                for result in results:
                    ppi.update_pipeline_summary(step, result)
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

        self.enabled_extensions = getattr(self.xp, 'extensions', [])

    def _par_detect_parquet_input_format(self):
        """Detect input format and return (format, files)"""
        fastq_files = sfind(f"{self.xp.pro_datain}", f"{self.xp.target_sample}*R1*.fastq.gz")
        bam_files = sfind(f"{self.xp.pro_datain}", f"{self.xp.target_sample}*.bam")
        
        if bam_files and fastq_files:
            self.logger.warning("Found both BAM and FASTQ files, prioritizing BAM")
            return 'bam', bam_files
        elif bam_files:
            return 'bam', bam_files  
        elif fastq_files:
            return 'fastq', fastq_files
        else:
            raise ValueError("No FASTQ or BAM files found")
        
    def _par_validate_format_params(self, input_format: str) -> None:
        """Validate required parameters for specific input format"""
        format_requirements = {
            'fastq': ['rev_comp'],
            'bam': ['anchor_ont', 'sbc_len', 'anchor_orient']
        }
        
        if input_format not in format_requirements:
            return  # No additional requirements
            
        required_params = format_requirements[input_format]
        missing_params = [param for param in required_params if not hasattr(self.xp, param)]
        
        if missing_params:
            param_list = ', '.join(missing_params)
            raise ValueError(f"Missing required parameter(s) for {input_format.upper()}: {param_list}")

    def _par_process_input_files(self, files, sample_dir, out_fn, merged_fn, input_format):
        """Converts FASTQ or BAM files into parquet"""
        
        force_tab = getattr(self.xp, "force_tab", False)
        limit = getattr(self.xp, "limit", None)
        
        if limit is not None:
            self.logger.warning(f"Only loading {limit} reads")

        max_files = getattr(self.xp, 'allow_wildcards', 1)
        
        if len(files) > 1 and max_files > 1:
            self.logger.info(f"Processing {len(files)} {input_format.upper()} files for sample {self.xp.target_sample}")
            temp_parsed_files = []
            temp_merged_files = []
            
            for i, file in enumerate(files):
                temp_out_fn = f'{sample_dir}/temp_parsed_reads_{i}.parquet'
                temp_merged_fn = f'{sample_dir}/temp_merged_reads_{i}.parquet'
                
                self.logger.info(f"Processing {input_format.upper()} file {i+1}/{len(files)}: {file}")

                self._par_process_single_file(file, temp_out_fn, temp_merged_fn, sample_dir, input_format, force_tab, limit)
                
                temp_parsed_files.append(temp_out_fn)
                temp_merged_files.append(temp_merged_fn)
            
            # Combine multiple parquet files (same for both formats)
            self.logger.info(f"Combining {len(temp_parsed_files)} {input_format.upper()}-derived parquet files")
            
            self._par_combine_parquet_files(temp_parsed_files, out_fn)
            self._par_combine_parquet_files(temp_merged_files, merged_fn)
                    
        else:
            # Single file processing
            file = files[0]
            self.logger.info(f"Processing single {input_format.upper()} file: {file}")
            
            # Format-specific processing
            self._par_process_single_file(file, out_fn, merged_fn, sample_dir, input_format, force_tab, limit)

    def _par_process_single_file(self, file, out_fn, merged_fn, sample_dir, input_format, force_tab, limit):
        """Process a single input file - format-specific logic"""
        
        if input_format == 'fastq':
            tabulate_paired_10x_fastqs_rs(
                file_path=file, 
                out_fn=out_fn,
                merged_fn=merged_fn,
                out_dir=sample_dir,
                modality=self.xp.modality,
                umi_len=self.xp.umi_len,
                do_rev_comp=self.xp.rev_comp,
                force=force_tab,
                limit=limit
            )
            
        elif input_format == 'bam':
            # Store raw BAM parquet in intermediate directory
            intermediate_dir = Path(sample_dir) / 'intermediate'
            intermediate_dir.mkdir(exist_ok=True)
            raw_bam_fn = str(intermediate_dir / 'raw_bam.parquet')

            import rogtk

            if not Path(raw_bam_fn).exists() or force_tab:
                rogtk.bam_to_parquet(
                    bam_path=file,
                    parquet_path=raw_bam_fn,
                    include_sequence=True,
                    include_quality=True,
                    compression='snappy',
                    limit=limit
                )
            (
                pl.scan_parquet(raw_bam_fn)
                .pp.ont_to_paired_format(
                    umi_len=self.xp.umi_len,
                    sbc_len=self.xp.sbc_len,
                    anchor_orient=self.xp.anchor_orient,
                    anchor_ont=self.xp.anchor_ont,
                    modality=self.xp.modality,
                    cbc_len=getattr(self.xp, 'cbc_len', 16),
                    do_rev_comp=self.xp.rev_comp,
                )
                .sink_parquet(merged_fn)
            )

            # For BAM, parsed = merged (no additional processing)
            import shutil
            shutil.copy2(merged_fn, out_fn)

            # Clean up intermediate file unless save_intermediate_files is set
            if not getattr(self.xp, 'save_intermediate_files', False):
                Path(raw_bam_fn).unlink(missing_ok=True)

    def _par_combine_parquet_files(self, file_list, output_file, cleanup=True):
        """Combine multiple parquet files into one (lazy/streaming)"""
        if file_list:
            # Use lazy evaluation to avoid loading all files into RAM
            pl.concat([pl.scan_parquet(f) for f in file_list]).sink_parquet(output_file)

            # Clean up temporary files
            if cleanup:
                for temp_file in file_list:
                    Path(temp_file).unlink(missing_ok=True)

    def _par_calculate_parquet_metrics(self, out_fn):
        """Calculate metrics from the processed parquet file"""
        limit = getattr(self.xp, "limit", None)
        
        return {
            'total_reads': (
                pl.scan_parquet(out_fn)
                .collect()
                .height
            ),
            # TODO  there is a bug here where 'umi' is not found
            # this wa probably introduced when importing from bam_path
            # the real issue is that the progression of file names is not correct
            # what is the difference between parsed and merged and why do parsed gets overwritten ?
            'total_umis': (
                pl.scan_parquet(out_fn)
                #.select('umi')
                #.unique()
                #
                #.select(pl.len())
                .collect()
                .height
                #.item()
            ),
            'downsampled': limit
        }

    @call  
    @pipeline_step(PipelineStep.DORADO)
    def dorado_basecall(self) -> StepResults|None:
        """Run dorado basecaller on POD5 files"""
        try:
            from ogtk.utils.db import run_dorado
            
            self.logger.info(f"Running dorado basecalling from {self.xp.pro_datain}")
            
            # Check if this is an Xp object with preprocessing config
            if not hasattr(self.xp, 'pp') or 'dorado' not in self.xp.pp:
                raise ValueError("Missing dorado configuration in xp.pp['dorado']")
            
            # Run dorado basecalling - save BAMs to datain/fracture directory
            result = run_dorado(self.xp, 
                              force=getattr(self.xp, 'force_dorado', False), 
                              dry=self.xp.dry,
                              bam_output_dir=self.xp.pro_datain)
            
            if result != 0 and not self.xp.dry:
                raise RuntimeError(f"Dorado basecalling failed with exit code: {result}")
            
            # Count generated BAM files
            from pathlib import Path
            output_dir = Path(self.xp.pro_datain)  # BAMs saved to datain/fracture
            bam_files = list(output_dir.glob("*.bam"))
            
            metrics = {
                'bam_files_generated': len(bam_files),
                'output_directory': str(output_dir)
            }
            
            self.logger.info(f"Dorado basecalling completed. Generated {len(bam_files)} BAM files")
            
            return StepResults(
                results={'bam_output_dir': str(output_dir), 'bam_files': [str(f) for f in bam_files]},
                metrics=metrics
            )
            
        except Exception as e:
            self.logger.error(f"Failed to run dorado basecalling: {str(e)}", with_traceback=True)
            raise

    @call  
    @pipeline_step(PipelineStep.PARQUET)
    def to_parquet(self) -> StepResults|None:
        """Convert input reads to parquet format (auto-detects FASTQ vs BAM)"""
        try:
            self.logger.info(f"Converting reads to parquet from {self.xp.pro_datain}")
            
            required_params = ['umi_len', 'modality']
            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")
            
            input_format, input_files = self._par_detect_parquet_input_format()
            self.logger.info(f"Detected input format: {input_format}")
            
            self._par_validate_format_params(input_format)
            
            sample_to_file = self.xp.organize_files_by_sample(
                input_files, 
                [{'id': self.xp.target_sample}], 
                max_files=getattr(self.xp, 'allow_wildcards', 1)
            )
            
            # Process each sample
            if not self.xp.dry:
                for sample_id, files in sample_to_file.items():
                    sample_dir = f'{self.xp.pro_workdir}/{sample_id}'
                    Path(sample_dir).mkdir(parents=True, exist_ok=True)
                    
                    out_fn = f'{sample_dir}/parsed_reads.parquet'
                    merged_fn = f'{sample_dir}/merged_reads.parquet'
                    
                    # ** main invocation ** #
                    self._par_process_input_files(files, sample_dir, out_fn, merged_fn, input_format)
                    
                    metrics = self._par_calculate_parquet_metrics(out_fn)
                    
            return StepResults(
                results={'parsed_fn': out_fn}, 
                metrics=metrics
            )
            
        except Exception as e:
            self.logger.error(f"Failed to convert to parquet: {str(e)}", with_traceback=True)
            raise

    def tto_parquet(self) -> StepResults|None:
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
            out_file_inv = f"{self.xp.sample_wd}/parsed_reads_invalid.parquet" #pyright:ignore

            out_file_qc = out_file.replace('.parquet', '_qc.parquet')
            total_reads = self.get_metric_from_summary("parquet", 'total_reads') 

            self.logger.info(f"reading from {in_file} with {total_reads/1e6:0.2f} M reads")

             
            if not self.xp.dry:
                ldf = (
                    pl
                    .scan_parquet(in_file)
                )

                self.logger.step(f"ldf")

                ldf = (
                    ldf
                    .pp.parse_reads(umi_len=self.xp.umi_len,  #pyright:ignore
                                    sbc_len=self.xp.sbc_len,
                                    anchor_ont=self.xp.anchor_ont,
                                    modality=self.xp.modality,
                                    cbc_len=getattr(self.xp, 'cbc_len', 16),
                                    )
                    .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                 )

                # save without gc fields
                if not self.xp.parse_read1:
                    self.logger.info(f"Exporting lazily to {out_file}")
                    (
                            ldf
                            .filter(pl.col('valid_umi'))
                            .filter(pl.col('ont'))
                            .sink_parquet(out_file)
                    )

                    self.logger.info(f"Exporting lazily to {out_file_inv}")
                    (
                            ldf
                            .filter((~pl.col('valid_umi')) | (~pl.col('ont')))
                            .sink_parquet(out_file_inv)
                    )
                else:
                    # TODO change the logics to an adaptive expression instead of boiler plate
                    # if r1 is long and informative: e.g long paired-end
                    self.logger.warning(f"extracting juice from R1")
                    (
                        pl.concat(
                            [
                                ldf
                                .filter(pl.col('valid_umi'))
                                .filter(pl.col('ont')),
                                ldf
                                .filter(pl.col('valid_umi'))
                                .filter(pl.col('ont'))
                                .pp.parse_read1(anchor_ont=self.xp.anchor_ont),
                            ]
                        ).sink_parquet(out_file)
                     )

                    (
                        pl.concat(
                            [
                                ldf
                                .filter((~pl.col('valid_umi')) | (~pl.col('ont'))),
                                ldf
                                .filter((~pl.col('valid_umi')) | (~pl.col('ont')))
                                .pp.parse_read1(anchor_ont=self.xp.anchor_ont),
                            ]
                        ).sink_parquet(out_file_inv)
                     )

                # extract QC fields and compute UMI-level metrics
                metrics_df = (
                        ldf
                        .select('umi', 'valid_umi', '^metric_.*$')
                        .unique()
                        .with_columns(metric_fraction_valid_umis = pl.col('valid_umi').mean())
                        .with_columns(metric_n_valid_umis = pl.col('valid_umi').sum())
                        .with_columns(metric_n_invalid_umis = pl.col('valid_umi').not_().sum())
                        .select('^metric_.*$')
                        .unique()
                        .collect()
                )

                metrics_df = metrics_df
                
                metrics = {f.replace('metric_', ''):metrics_df[f][0] 
                           for f in metrics_df.columns}

                (
                    metrics_df.
                    rename({f:f.replace('metric_', '') for f in metrics_df.columns})
                    .write_parquet(out_file_qc)
                )

                self.logger.info(f"exported parsed reads to {out_file}")
                self.logger.info(f"exported parsed QC metrics to {out_file_qc}")

                return StepResults(
                        results={'xp': self.xp,
                                 'parsed_reads': out_file,
                                 'parsed_reads_invalid': out_file_inv,
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

            if not hasattr(self.xp, 'fracture'):
                setattr(self.xp, 'fracture', {})
                self.xp.fracture['start_min_coverage'] = 25
                self.xp.fracture['start_k'] = 17
                self.xp.fracture['min_reads'] = 10 #qc.find_read_count_threshold(in_file, method="kmeans")
            #else:
            #    self.xp.fracture['min_reads'] = 10 #qc.find_read_count_threshold(in_file,
                                                   #              method=self.xp.fracture['method_th'])

            # Check if segmented assembly should be used
            use_segmentation = self.xp.fracture.get('use_segmentation', False)
            if use_segmentation:
                # Load METAs from features CSV or dedicated metas file
                metas_csv = self.xp.fracture.get('metas_csv') or getattr(self.xp, 'features_csv', None)
                if not metas_csv:
                    raise ValueError("Segmented assembly requires 'metas_csv' or 'features_csv' configuration")
                self.logger.info(f"Segmented assembly enabled, loading METAs from {metas_csv}")
                metas_df = pl.read_csv(metas_csv)

            in_files = {
                    'valid': f"{self.xp.sample_wd}/parsed_reads.parquet", #pyright:ignore
                    'invalid':f"{self.xp.sample_wd}/parsed_reads_invalid.parquet", #pyright:ignore
                    }

            for key,in_file in in_files.items():
                if not self.xp.dry:
                    self.logger.info(f'Reading {in_file}')

                    ldf = pl.scan_parquet(in_file)

                    # Apply masking if configured (skip if using segmentation - they're alternative strategies)
                    if hasattr(self.xp, 'features_csv') and self.xp.features_csv and not use_segmentation:
                        self.logger.info(f"Masking repetitive sequences using {self.xp.features_csv}")
                        n_reads = ldf.select(pl.len()).collect().item()

                        # Apply reverse complement if configured
                        if getattr(self.xp, 'mask_reverse_complement', False):
                            self.logger.info("Applying reverse complement to r2_seq before masking")
                            ldf = ldf.with_columns(
                                pl.col('r2_seq').dna.reverse_complement().alias('r2_seq') #pyright:ignore
                            )

                        # Flank-based masking works for both edited and unedited cassettes
                        ldf = ldf.pp.mask_flanks(
                            features_csv=self.xp.features_csv,
                            column_name='r2_seq',
                            fuzzy_pattern=getattr(self.xp, 'mask_fuzzy_pattern', True),
                            fuzzy_kwargs=getattr(self.xp, 'mask_fuzzy_kwargs', None)
                        )

                        self.logger.info(f"Applied masking to {n_reads} reads")

                        if getattr(self.xp, 'save_intermediate_files', False):
                            intermediate_dir = Path(self.xp.sample_wd) / 'intermediate'
                            intermediate_dir.mkdir(exist_ok=True)
                            masked_file = intermediate_dir / f"masked_reads_{key}.parquet"
                            self.logger.info(f"Saving masked reads to {masked_file}")
                            ldf.sink_parquet(str(masked_file))
                            ldf = pl.scan_parquet(str(masked_file))

                    filter_expr= pl.col('reads')>=self.xp.fracture['min_reads']

                    # Determine strategy name for output file
                    if use_segmentation:
                        strategy = 'segmented'
                    elif hasattr(self.xp, 'features_csv') and self.xp.features_csv:
                        strategy = 'masked'
                    else:
                        strategy = 'direct'

                    out_file = f"{self.xp.sample_wd}/contigs_{strategy}_{key}.parquet"
                    self.logger.info(f"Exporting assembled contigs ({strategy}) to {out_file}")

                    if use_segmentation:
                        debug_path = None
                        if getattr(self.xp, 'save_intermediate_files', False):
                            intermediate_dir = Path(self.xp.sample_wd) / 'intermediate'
                            intermediate_dir.mkdir(exist_ok=True)
                            debug_path = str(intermediate_dir / "segments_debug.parquet")

                        assembly_method = self.xp.fracture.get('assembly_method', 'compression')

                        # Get heterogeneity threshold from config (default 0.20 = 20% of dominant)
                        heterogeneity_threshold = self.xp.fracture.get('heterogeneity_threshold', 0.20)

                        df_contigs = (
                            ldf
                            .filter(filter_expr)
                            .pp.assemble_segmented(
                                metas=metas_df,
                                cassette_start_anchor=self.xp.start_anchor,
                                cassette_end_anchor=self.xp.end_anchor,
                                k=self.xp.fracture['start_k'],
                                min_coverage=self.xp.fracture['start_min_coverage'],
                                debug_path=debug_path,
                                method=assembly_method,
                                heterogeneity_threshold=heterogeneity_threshold,
                            )
                            .with_columns(pl.col('contig').str.len_chars().alias('length'))
                            .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                            .join(
                                pl.scan_parquet(in_file)
                                  .select('umi','reads')
                                  .filter(filter_expr)
                                  .unique(),
                               left_on='umi', right_on='umi', how='left')
                        )
                    else:
                        # Standard assembly
                        df_contigs = (
                                ldf
                                .filter(filter_expr)
                                .pp.assemble_umis( #pyright: ignore
                                  k=self.xp.fracture['start_k'],
                                  min_coverage=self.xp.fracture['start_min_coverage'],
                                  start_anchor=self.xp.start_anchor,
                                  end_anchor=self.xp.end_anchor,
                                  min_length=getattr(self.xp, 'min_contig_length', None),
                                  auto_k=False,
                                  export_graphs=False,
                                  only_largest=True,
                                  method=self.xp.fracture['assembly_method'],
                                  modality=self.xp.modality,
                                  )
                                .with_columns(pl.col('contig').str.len_chars().alias('length'))
                                .with_columns(pl.lit(self.xp.target_sample).alias('sample_id'))
                                .join(
                                    pl.scan_parquet(in_file)
                                      .select('umi','reads')
                                      .filter(filter_expr)
                                      .unique(),
                                      #.collect(),
                                   left_on='umi', right_on='umi', how='left')
                                )

                    # Unmask contigs if masking was applied (skip if using segmentation)
                    if hasattr(self.xp, 'features_csv') and self.xp.features_csv and not use_segmentation:
                        self.logger.info("Restoring original sequences (unmasking contigs)")
                        df_contigs = df_contigs.pp.unmask_flanks(
                            features_csv=self.xp.features_csv,
                            column_name='contig',
                            fuzzy_pattern=getattr(self.xp, 'mask_fuzzy_pattern', True),
                            fuzzy_kwargs=getattr(self.xp, 'mask_fuzzy_kwargs', None)
                        )
                        # Recalculate length after unmasking
                        df_contigs = df_contigs.with_columns(
                            pl.col('contig').str.len_chars().alias('length')
                        )

                    df_contigs.sink_parquet(out_file)

                    success_rate = (df_contigs.collect().get_column('length')>0).mean()
                    total_assembled = (df_contigs.collect().get_column('length')>0).sum()
                   
                    self.logger.critical(f"from {in_file} {total_assembled} assembled ({success_rate:.2%} success) {key}")
                    
                    return StepResults(
                            results={'xp': self.xp, 
                                     'contigs':out_file},
                            metrics={'success_rate':success_rate,
                                     'total_assembled':total_assembled,
                                     'mean_contig_length': df_contigs.collect().filter(pl.col('length')>0).get_column('length').mean(),
                                     'median_contig_length': df_contigs.collect().filter(pl.col('length')>0).get_column('length').median(),
                                     },
                            )

                pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}", with_traceback=True)
            raise
    @call
    def run_extensions(self) -> bool:
        """Run enabled post-processing extensions"""
        # Check if extensions are configured
        if not hasattr(self.xp, 'extensions') or not self.xp.extensions:
            self.logger.debug("No extensions configured")
            return True
        
        try:
            # Import here to avoid circular dependency
            from ..extensions import extension_registry
        except ImportError:
            self.logger.warning("Extensions module not available")
            return True
        
        # Look for contigs file (check all strategy variants)
        # on valid reads only
        key = 'valid'
        contigs_path = None
        for strategy in ['segmented', 'masked', 'direct']:
            candidate = Path(f"{self.xp.sample_wd}/contigs_{strategy}_{key}.parquet")
            if candidate.exists():
                contigs_path = candidate
                self.logger.info(f"Found contigs file: {contigs_path.name}")
                break

        if not contigs_path:
            self.logger.error("No contigs found for post-processing. Run fracture step first.")
            return False
        
        self.logger.info(f"Running {len(self.xp.extensions)} extensions")
        
        for ext_name in self.xp.extensions:
            extension = extension_registry.create_extension(ext_name, self.xp)
            if not extension:
                self.logger.error(f"Extension '{ext_name}' not found")
                self.logger.info(f"Available extensions: {extension_registry.get_available()}")
                continue

            # Set up extension-specific logging
            ext_log_dir = Path(self.xp.sample_wd) / "logs"
            ext_log_dir.mkdir(parents=True, exist_ok=True)
            ext_log_file = ext_log_dir / f"ext_{ext_name}.log"

            ext_logger = Rlogger()
            ext_logger.enable_file_logging(
                filepath=ext_log_file,
                level=self.xp.log_level if hasattr(self.xp, 'log_level') else "INFO",
                format_string="%(asctime)s - %(levelname)s - %(message)s",
            )

            try:
                self.logger.info(f"Running extension: {ext_name}")
                self.logger.info(f"Extension log file: {ext_log_file}")
                results = extension.process(contigs_path)

                # Create a synthetic step for the extension summary
                class ExtStep:
                    def __init__(self, name):
                        self.name = f"EXT_{name.upper()}"
                        self.value = {'required_params': extension.required_params}

                self.update_pipeline_summary(ExtStep(ext_name), results)
                self.logger.info(f"Extension {ext_name} completed successfully")
                ext_logger.disable_file_logging()

            except Exception as e:
                self.logger.error(f"Extension {ext_name} failed: {e}", with_traceback=True)
                ext_logger.disable_file_logging()
                return False
        
        self.logger.info("Extensions completed successfully!")
        return True        

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
    def should_run_step(self, step: PipelineStep) -> bool:
        """Check if a step should be run based on configuration"""
        if hasattr(self.xp, 'steps') and self.xp.steps is not None:
            return step.name.lower() in [s.lower() for s in self.xp.steps] 
        return True  # Run all steps if not specified

    def run(self) -> bool:
        """Run the complete pipeline and/or extensions"""
        if self.xp.make_test:
            try:
                self.logger.info("=== TEST MODE ===")
                self.xp.target_sample = f'TEST_{self.xp.target_sample}'
                self.xp.samples.append({'id' : self.xp.target_sample})
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
            has_main_steps = (hasattr(self.xp, 'steps') and 
                 self.xp.steps is not None and 
                 len(self.xp.steps) > 0)

            if has_main_steps:
                self.run_all()
            else:
                self.logger.info("No main pipeline steps specified, skipping to extensions")
            
            self.logger.info("Pipeline completed successfully")

            self.run_extensions()
            
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}", with_traceback=True)
            return False

    def run_all(self):
        ''' Run all steps (respecting step configuration)'''
        self.logger.info("==== Running pipeline ====")
        
        # Run steps in order, but only if configured to run
        if self.should_run_step(PipelineStep.DORADO):
            self.dorado_basecall()
        if self.should_run_step(PipelineStep.PARQUET):
            self.to_parquet()
        if self.should_run_step(PipelineStep.PREPROCESS):
            self.preprocess()
        if self.should_run_step(PipelineStep.FRACTURE):
            self.fracture()
        #self.run_qcs()

    def get_metric_from_summary(self, step: str, metric: str):
        import json
        summary_path = Path(f"{self.xp.pro_workdir}/{self.xp.target_sample}/pipeline_summary.json")
        if summary_path.exists():
            try:
                with open(summary_path, 'r') as f:
                    summary = json.load(f)
            except json.JSONDecodeError:
                summary = {}
        else:
            summary = {}
        return summary[step]['metrics'][metric]


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
                output_files.append(results.results['parsed_reads_invalid'])
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

