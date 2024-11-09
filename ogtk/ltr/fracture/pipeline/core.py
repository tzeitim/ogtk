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
                            out_dir='{sample_dir}/merged_reads.parquet',
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
    @pipeline_step(PipelineStep.FRACTURE)
    def test(self) -> None:
        """Makes a downsampled fastq file for testing purposes"""
            
        try:
            self.logger.step("Downsampling parquet to generate synthetic FASTQ")

            required_params = ['target_sample']

            for param in required_params:
                if param not in vars(self.xp):
                    raise ValueError(f"Missing required parameter: {param}")

            in_file = f"{self.xp.pro_datain}/{self.xp.target_sample}_R1_001.merged.parquet"
            out_file1= f"{self.xp.pro_datain}/TEST_{self.xp.target_sample}_R1_001.fastq.gz"
            out_file2= f"{self.xp.pro_datain}/TEST_{self.xp.target_sample}_R2_001.fastq.gz" 
            self.logger.io(f"exporting reads to {out_file1}\n{out_file2}")

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

    def run(self) -> bool:
        """Run the complete pipeline"""
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

