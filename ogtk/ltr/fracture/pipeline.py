from pathlib import Path
import polars as pl
import pandas as pd
from functools import wraps
import argparse
from enum import Enum, auto
from typing import Any, NamedTuple, Callable, List
from ogtk.utils.log import CustomLogger, Rlogger, call
from ogtk.utils import tabulate_paired_10x_fastqs_rs, sfind #, another_function, yet_another_function
from ogtk.utils.db import Xp
from dataclasses import dataclass
from ogtk.ltr.fracture.api_ext import PlDNA, PlPipeline, PllPipeline


class StepResults(NamedTuple):
    """Container for data to be passed to plotting"""
    results: dict  
    
class PipelineStep(Enum):
    """Enum defining available pipeline steps"""
    PARQUET = auto()
    PREPROCESS = auto()
    FRACTURE = auto()
    SAVE = auto()
    TEST = auto()
    
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
        def wrapper(pipeline_instance: Any, *args: Any, **kwargs: Any) -> Any:
            if not pipeline_instance.should_run_step(step):
                pipeline_instance.logger.io(f"Skipping {step.name.lower()}")
                return None
            
            result = func(pipeline_instance, *args, **kwargs)
            
            if getattr(pipeline_instance.xp, 'do_plot', False):
                try:
                    pipeline_instance.logger(f'{step.name.lower()}')
                    plot_method = getattr(pipeline_instance, f"plot_{step.name.lower()}", None)
                    if plot_method:
                        pipeline_instance.logger.debug(f"Generating plots for {step.name.lower()}")
                        plot_method(result)
                except Exception as e:
                    pipeline_instance.logger.error(f"Plot generation failed: {str(e)}")
            
            return result
        return wrapper
    return decorator

@call
def plot_preprocess(results):
    ''' ''' 
    xp = results['xp'] 
    ifn = results['parsed_reads']
    out_path = f'{xp.sample_figs}/{xp.target_sample}_coverage.png'
    fig = sns.displot(data=
        pl.scan_parquet(ifn)
            .select('umi', 'reads').unique()
            .collect(),
            y='reads', 
            log_scale=(10, 10), 
            kind='ecdf', 
            complementary=True, 
            stat='count')

    plt.grid()
    plt.title("Reads per UMI")

    xp.logger.io(f"saved {out_path}")
    fig.savefig(out_path)


class FractureXp(Xp):
    """Extended Xp class with pipeline-specific functionality"""
    steps: Any
    dry: bool
    logger: CustomLogger
    pro_datain: str
    modality: str
    umi_len: int
    rev_comp: bool
    anchor_ont: str
    samples: List[str]
    pro_workdir: str

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.steps = getattr(self, 'steps', None)
        self.dry = getattr(self, 'dry', False)

    def organize_files_by_sample(self, files, samples):
        # Create a dictionary to store files for each sample
        sample_files = {sample['id']: [] for sample in samples}
        
        # Sort files into appropriate sample groups
        for f in files:
            for sample_id in sample_files:
                if sample_id in f:
                    sample_files[sample_id].append(f)
                    break
        
        return sample_files

class Pipeline:
    """Main pipeline class for data processing using Xp configuration"""
    dry: bool
    def __init__(self, xp: FractureXp):
        self.xp = xp
        self.logger = Rlogger().get_logger()
        
        # Create output directory if it doesn't exist
        self.xp.init_workdir()


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

            sample_to_file = self.xp.organize_files_by_sample(self.xp.samples, input_files)

            if not self.xp.dry:
                for sample_id, file in sample_to_file.items():
                    self.logger.io(f'{sample_id} {file}')
                    sample_dir = f'{self.xp.pro_workdir}/{sample_id}/' 
                    Path(sample_dir).mkdir(parents=True, exist_ok=True)

                    tabulate_paired_10x_fastqs_rs(
                            file_path=file, 
                            out_fn=f'{sample_dir}/parsed_reads.parquet', 
                            modality=self.xp.modality,     
                            umi_len=self.xp.umi_len,      
                            do_rev_comp=self.xp.rev_comp,  
                            force=True) # add to step interface
                    print(file)

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
                    .pp.parse_reads(umi_len=self.xp.umi_len, anchor_ont=self.xp.anchor_ont) 
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
                    .pp.parse_reads(umi_len=self.xp.umi_len, anchor_ont=self.xp.anchor_ont) 
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
                    .pp.parse_reads(umi_len=self.xp.umi_len, anchor_ont=self.xp.anchor_ont) 
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

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Data processing pipeline")
    parser.add_argument("--config", required=True, help="Path to config file")
    parser.add_argument(
        "--steps",
        nargs="+",
        choices=[s.name.lower() for s in PipelineStep],
        help="Specific steps to run (overrides config file)"
    )
    parser.add_argument("--log-level", 
                       default="INFO",
                       choices=["CRITICAL", "INFO", "IO", "STEP", "DEBUG"],
                       help="Logging level")
    return parser.parse_args()


def main():
    '''
    # Run all steps defined in config
    python pipeline.py --config config.yml

    # Run only specific steps (overrides config)
    python pipeline.py --config config.yml --steps load analyze save

    # Set log level
    python pipeline.py --config config.yml --log-level DEBUG
    '''

    # Parse arguments
    args = parse_args()
    
    # Initialize logger
    logger = Rlogger().get_logger()
    Rlogger().set_level(args.log_level)
    
    try:
        # Initialize Xp configuration
        xp = FractureXp(conf_fn=args.config)
        
        # Override steps if specified in command line
        if args.steps:
            xp.steps = args.steps

        if getattr(xp, 'do_plot', False):
            import seaborn as sns
            import matplotlib.pyplot as plt
        
        # Initialize and run pipeline
        
        pipeline = Pipeline(xp)
        success = pipeline.run()
        
        return 0 if success else 1
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        return 1

if __name__ == "__main__":
    exit(main())
