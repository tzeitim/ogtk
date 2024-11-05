from pathlib import Path
import polars as pl
import pandas as pd
from functools import wraps
import argparse
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Set
from dataclasses import dataclass
from ogtk.utils.log import Rlogger, call
from ogtk.utils import tabulate_paired_10x_fastqs_rs, sfind #, another_function, yet_another_function

@pl.api.register_dataframe_namespace("dna")
class PlDNA:
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df

    def to_fasta(self, read_id_col: str, read_col: str) -> pl.DataFrame:
        return  self._df.with_columns(
                (">"+pl.col(read_id_col)\
                 +"\n"+pl.col(read_col)\
                 +"\n")
                 .alias(f'{read_col}_fasta')
        )
    def to_fastq(self, read_id_col: str, read_qual_col: str, read_col: str)-> pl.DataFrame:
        return self._df.with_columns(
                ("@"+pl.col(read_id_col)\
                 +"\n"+pl.col(read_col)\
                 +"\n+"\
                 +"\n"+pl.col(read_qual_col)\
                 +"\n"\
                )
                 .alias(f'{read_col}_fastq')
        )


@pl.api.register_dataframe_namespace("pp")
class PlPipeline:
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df

@pl.api.register_lazyframe_namespace("pp")
class PllPipeline:
    def __init__(self, ldf: pl.LazyFrame) -> None:
        self._ldf = ldf

    def parse_reads(self, umi_len, anchor_ont):
        return(
                self._ldf
                    .with_columns(pl.col('r1_seq').str.slice(0, umi_len).alias('umi'))
                    .with_columns(pl.col('r1_seq').str.contains(anchor_ont).alias('ont'))
                    .filter(pl.col('ont'))
                    .with_columns(pl.len().over('umi').alias('reads'))
                    # strip out the UMI from R1
                    .with_columns(pl.col('r1_seq').str.replace(f"^.+?{anchor_ont}", anchor_ont))
                )    

class PipelineStep(Enum):
    """Enum defining available pipeline steps"""
    PARQUET = auto()
    PREPROCESS = auto()
    ANALYZE = auto()
    SAVE = auto()
    
    @classmethod
    def from_string(cls, step_name: str) -> "PipelineStep":
        """Convert string to PipelineStep enum"""
        try:
            return cls[step_name.upper()]
        except KeyError:
            raise ValueError(f"Invalid step name: {step_name}. Valid steps are: {[s.name for s in cls]}")

class Pipeline:
    """Main pipeline class for data processing using Xp configuration"""

    def __init__(self, xp: Xp):
        self.xp = xp
        self.logger = Rlogger().get_logger()
        
        # Extract dry run setting before consolidation
        self.dry_run = getattr(self.xp, 'dry', False)
        
        # Create output directory if it doesn't exist
        Path(self.xp.output_dir).mkdir(parents=True, exist_ok=True)
        if self.dry_run:
            self.logger.critical("Running in dry run mode!")

    def should_run_step(self, step: PipelineStep) -> bool:
        """Check if a step should be run based on configuration"""
        if hasattr(self.xp, 'steps'):
            return step.name.lower() in [s.lower() for s in self.xp.steps]
        return True  # Run all steps if not specified
        
    @call
    def validate_input(self) -> bool:
        """Validate input data and configuration"""
        try:
            self.logger.step("Validating input data and configuration")
            
            # Check required attributes
            required_attrs = ['input_dir', 'output_dir']
            for attr in required_attrs:
                if not hasattr(self.xp, attr):
                    self.logger.error(f"Missing required configuration: {attr}")
                    return False
                
            # Validate step dependencies
            if (hasattr(self.xp, 'steps') and 
                'analyze' in self.xp.steps and 
                'parquet' not in self.xp.steps):
                self.logger.warning("Analysis step requires data loading step")
                return False
                
            return True
        except Exception as e:
            self.logger.error(f"Validation failed: {str(e)}")
            return False
            
    @call
    def to_parquet(self) -> None:
        """Convert fastq reads to parquet format"""
        if not self.should_run_step(PipelineStep.PARQUET):
            self.logger.io("Skipping data loading step")
            return
            
        try:
            self.logger.info(f"Loading data from {self.config.input_dir}")
            input_files = sfind(f"{self.config.input_dir}", "*.fastq.gz")
            self.logger.io(f"found {input_files}")

            if not self.config.dry:
                for file in input_files:
                    tabulate_paired_10x_fastqs_rs(
                            file_path=file, 
                            out_fn=file.replace("fastq.gz", '.parquet'), 
                            modality=self.config.parameters.modality,     #pyright: ignore
                            umi_len=self.config.parameters.umi_len,       #pyright: ignore
                            do_rev_comp=self.config.parameters.rev_comp,  #pyright: ignore
                            force=True)
                    print(file)

            pass
        except Exception as e:
            self.logger.error(f"Failed to load data: {str(e)}")
            raise
            
    @call  
    def preprocess_data(self) -> None:
        """Preprocess the loaded data"""
        if not self.should_run_step(PipelineStep.PREPROCESS):
            self.logger.io("Skipping preprocessing step")
            return
            
        try:
            self.logger.step("Parsing reads and collecting molecule stats")

            parameters = self.config.parameters

            parameters = self.xp.parameters #pyright:ignore
            required_params = ['umi_len', 'anchor_ont']

            for param in required_params:
                if param not in parameters:
                    raise ValueError(f"Missing required parameter: {param}")

            out_file = f"{self.xp.output_dir}/{self.xp.sample_id}.dff.parquet" #pyright:ignore
            self.logger.io(f"exporting reads to {out_file}")

            if not getattr(self.xp, 'dry', False):
                (
                    pl
                    .scan_parquet(f"{input_dir}/{sample}_R1_001.merged.parquet")
                    .pp.parse_reads(umi_len=umi_len, anchor_ont=anchor_ont) #pyright:ignore
                    .collect()
                    .write_parquet(out_file)
                )
            pass
        except Exception as e:
            self.logger.error(f"Failed to preprocess data: {str(e)}")
            raise
            
    @call
    def analyze_data(self) -> None:
        """Perform data analysis"""
        if not self.should_run_step(PipelineStep.ANALYZE):
            self.logger.io("Skipping analysis step")
            return
            
        try:
            self.logger.step("Analyzing data")
            # Add your analysis steps here
            pass
        except Exception as e:
            self.logger.error(f"Failed to analyze data: {str(e)}")
            raise
            
    @call
    def save_results(self) -> None:
        """Save processing results"""
        if not self.should_run_step(PipelineStep.SAVE):
            self.logger.io("Skipping save step")
            return
            
        try:
            self.logger.io(f"Saving results to {self.xp.output_dir}")
            # Add your saving logic here
            pass
        except Exception as e:
            self.logger.error(f"Failed to save results: {str(e)}")
            raise
            
    def run(self) -> bool:
        """Run the complete pipeline"""
        try:
            if not self.validate_input():
                return False
                
            self.to_parquet()
            self.preprocess_data()
            self.analyze_data()
            self.save_results()
            
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
        xp = Xp(conf_fn=args.config)
        xp.consolidate_conf()  # Ensure configuration is consolidated
        
        # Override steps if specified in command line
        if args.steps:
            xp.steps = args.steps
        
        # Initialize and run pipeline
        pipeline = Pipeline(xp)
        success = pipeline.run()
        
        return 0 if success else 1
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        return 1

if __name__ == "__main__":
    exit(main())
