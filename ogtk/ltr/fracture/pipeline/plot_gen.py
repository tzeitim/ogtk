from pathlib import Path
from typing import Optional, List
from .core import Pipeline, PipelineStep, StepResults

class PlotRegenerator:
    """Utility class for regenerating plots from previous pipeline runs"""
    
    def __init__(self, pipeline: Pipeline):
        self.pipeline = pipeline
        self.xp = pipeline.xp
        self.logger = pipeline.logger
        
    def regenerate_plots(self, 
                        steps: Optional[List[str]] = None, 
                  ) -> None:
        """Regenerate plots for specified steps or all steps if none specified
        
        Args:
            steps: List of step names to regenerate plots for. If None, regenerates all.
            force: If True, regenerates even if output files exist
        """
        if not self.xp.do_plot:
            self.logger.info("Plotting is disabled in configuration")
            return
            
        available_steps = [step.name.lower() for step in PipelineStep]
        
        if steps:
            steps_to_run = [s.lower() for s in steps]
            invalid_steps = [s for s in steps_to_run if s not in available_steps]
            if invalid_steps:
                raise ValueError(f"Invalid step names: {invalid_steps}. Valid steps are: {available_steps}")
        else:
            steps_to_run = available_steps
            
        for step_name in steps_to_run:
            self.logger.step(f"Regenerating plots for {step_name}")
            
            # Get plot method from PlotDB
            plot_method = getattr(self.pipeline.plotdb, f"plot_{step_name}", None)
            if not plot_method:
                self.logger.info(f"No plotting method found for step {step_name}")
                continue
                
            # Get results from previous run
            results = self._load_step_results(step_name)
            if not results:
                self.logger.info(f"No previous results found for step {step_name}")
                continue
                
            try:
                plot_method(self.pipeline, results)
                self.logger.info(f"Successfully regenerated plots for {step_name}")
            except Exception as e:
                self.logger.error(f"Failed to regenerate plots for {step_name}: {str(e)}")
                
    def _load_step_results(self, step_name: str) -> Optional[StepResults]:
        """Load results from previous pipeline run for given step"""
        step = PipelineStep[step_name.upper()]
        
        # Example for preprocess step
        if step == PipelineStep.PREPROCESS:
            parsed_reads = f"{self.xp.sample_wd}/parsed_reads.parquet"
            if Path(parsed_reads).exists():
                return StepResults(results={'xp': self.xp, 'parsed_reads': parsed_reads})
                
        # Example for fracture step    
        elif step == PipelineStep.FRACTURE:
            contigs = f"{self.xp.sample_wd}/contigs.parquet"
            if Path(contigs).exists():
                return StepResults(results={'xp': self.xp, 'contigs': contigs})
                
        return None
