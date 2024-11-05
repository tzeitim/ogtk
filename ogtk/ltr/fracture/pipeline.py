from pathlib import Path
import yaml
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

@dataclass
class PipelineConfig:
    """Configuration class for pipeline parameters"""
    input_dir: Path
    output_dir: Path
    steps: Set[PipelineStep]  # Steps to run
    parameters: Dict[str, Any]
    dry: bool
    
    @classmethod
    def from_yaml(cls, config_path: str, override_steps: Optional[List[str]] = None) -> "PipelineConfig":
        """
        Load pipeline configuration from YAML file
        
        Args:
            config_path: Path to YAML config file
            override_steps: Optional list of steps to override config file
        """
        with open(config_path) as f:
            config = yaml.safe_load(f)
            
        # Get steps from config or use all steps if not specified
        config_steps = config.get("steps", [s.name for s in PipelineStep])
        
        # Override steps if provided
        steps_to_use = override_steps if override_steps is not None else config_steps
        
        # Convert step names to PipelineStep enums
        steps = {PipelineStep.from_string(step) for step in steps_to_use}
            
        return cls(
            input_dir=Path(f'{config["prefix"]}/{config["input_dir"]}'),
            output_dir=Path(f'{config["prefix"]}/{config["output_dir"]}'),
            steps=steps,
            parameters=config.get("parameters", {}),
            dry=config.get("dry", False)
        )

class Pipeline:
    """Main pipeline class for data processing"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = Rlogger().get_logger()
        self.data = None
        
        # Create output directory if it doesn't exist
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        if self.config.dry:
            self.logger.critical("Running dry!")

    def should_run_step(self, step: PipelineStep) -> bool:
        """Check if a step should be run based on configuration"""
        return step in self.config.steps
        
    @call
    def validate_input(self) -> bool:
        """Validate input data and configuration"""
        try:
            self.logger.step("Validating input data and configuration")
            
            # Validate step dependencies
            if PipelineStep.ANALYZE in self.config.steps and PipelineStep.PARQUET not in self.config.steps:
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

            input_dir= self.config.input_dir
            output_dir= self.config.output_dir
            sample = parameters['sample']  
            umi_len = parameters['umi_len'] 
            anchor_ont = parameters['anchor_ont'] 


            out_file = f"{output_dir}/{sample}.dff.parquet"
            self.logger.io(f"exporting reads to {out_file}")

            if not self.config.dry:
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
            self.logger.io(f"Saving results to {self.config.output_dir}")
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
        # Load configuration with optional step override
        config = PipelineConfig.from_yaml(
            args.config,
            override_steps=args.steps
        )
        
        # Initialize and run pipeline
        pipeline = Pipeline(config)
        success = pipeline.run()
        
        return 0 if success else 1
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        return 1

if __name__ == "__main__":
    exit(main())
