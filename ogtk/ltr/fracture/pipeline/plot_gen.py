from pathlib import Path
from typing import Optional, List, Set
from .core import Pipeline, PipelineStep, StepResults
from .formats import read_file


def get_available_extensions() -> Set[str]:
    """Get available extension names from registry."""
    try:
        from ..extensions import extension_registry
        return set(extension_registry.get_available())
    except ImportError:
        return set()

class SummaryRegenerator:
    """Utility class for regenerating the markdown summary from previous pipeline runs"""
    
    def __init__(self, pipeline: Pipeline):
        self.pipeline = pipeline
        self.xp = pipeline.xp
        self.logger = pipeline.logger
        
    def regenerate_summary(self) -> None:
        """Regenerate the markdown summary from existing JSON summary file.
        
        This loads the existing JSON summary file and generates a new 
        markdown report based on its contents.
        """
        from pathlib import Path
        import json
        
        # Path to the summary file
        summary_path = Path(f"{self.xp.pro_workdir}/{self.xp.target_sample}/pipeline_summary.json")
        
        if not summary_path.exists():
            self.logger.warning(f"No pipeline summary found at {summary_path}")
            return
            
        try:
            # Load the summary JSON
            with open(summary_path, 'r') as f:
                summary = json.load(f)
                
            # Generate a new markdown report
            self.logger.step(f"Regenerating markdown summary for {self.xp.target_sample}")
            self.pipeline._generate_markdown_report(summary)
            self.logger.info("Successfully regenerated markdown summary")
            
        except Exception as e:
            self.logger.error(f"Failed to regenerate markdown summary: {str(e)}", with_traceback=True)

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

        # Pipeline steps + additional plot-only steps like segmentation, checkhealth + extensions
        extensions = [f"ext_{ext}" for ext in get_available_extensions()]
        available_steps = [step.name.lower() for step in PipelineStep] + ['segmentation', 'checkhealth'] + extensions

        if steps:
            steps_to_run = [s.lower() for s in steps]
            invalid_steps = [s for s in steps_to_run if s not in available_steps]
            if invalid_steps:
                raise ValueError(f"Invalid step names: {invalid_steps}. Valid steps are: {available_steps}")
        else:
            steps_to_run = available_steps
            
        for step_name in steps_to_run:
            self.logger.step(f"Regenerating plots for {step_name}")

            # Handle extensions separately (prefixed with ext_)
            if step_name.startswith('ext_'):
                ext_name = step_name[4:]  # Remove 'ext_' prefix
                results = self._load_extension_results(ext_name)
                if results:
                    self._regenerate_extension_plots(ext_name, results)
                else:
                    self.logger.info(f"No previous results found for extension {ext_name}")
                continue

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
                self.logger.error(f"Failed to regenerate plots for {step_name}: {str(e)}", with_traceback=True)
                
    def _load_step_results(self, step_name: str) -> Optional[StepResults]:
        """Load results from previous pipeline run for given step"""
        # Handle segmentation separately (not a PipelineStep enum member)
        if step_name.lower() == 'segmentation':
            segments = f"{self.xp.sample_wd}/intermediate/segments_debug.parquet"
            assembled = f"{self.xp.sample_wd}/intermediate/assembled_debug.parquet"
            contigs = f"{self.xp.sample_wd}/contigs_segmented_valid.parquet"
            if all(Path(f).exists() for f in [segments, assembled, contigs]):
                return StepResults(results={
                    'xp': self.xp,
                    'segments': segments,
                    'assembled': assembled,
                    'contigs_segmented_valid': contigs,
                })
            return None

        # Handle checkhealth separately
        if step_name.lower() == 'checkhealth':
            raw_bam = f"{self.xp.sample_wd}/intermediate/raw_bam.parquet"
            if not Path(raw_bam).exists():
                raw_bam = f"{self.xp.sample_wd}/intermediate/raw_bam.arrow"
            if Path(raw_bam).exists():
                return StepResults(results={
                    'xp': self.xp,
                    'raw_bam': raw_bam,
                })
            return None

        step = PipelineStep[step_name.upper()]

        # Example for preprocess step
        if step == PipelineStep.PREPROCESS:
            parsed_reads = f"{self.xp.sample_wd}/parsed_reads.parquet"
            parsed_reads_invalid = f"{self.xp.sample_wd}/parsed_reads_invalid.parquet"
            if Path(parsed_reads).exists():
                return StepResults(results={
                    'xp': self.xp,
                    'parsed_reads': parsed_reads,
                    'parsed_reads_invalid': parsed_reads_invalid,
                })

        # Example for fracture step
        elif step == PipelineStep.FRACTURE:
            contigs = f"{self.xp.sample_wd}/contigs.parquet"
            if Path(contigs).exists():
                return StepResults(results={'xp': self.xp, 'contigs': contigs})

        return None

    def _load_extension_results(self, ext_name: str) -> Optional[StepResults]:
        """Load results for regenerating extension plots."""
        if ext_name == 'cassiopeia_petracer':
            intermediate_dir = Path(self.xp.sample_wd) / "intermediate"
            segments_path = intermediate_dir / "segments_debug.parquet"
            if not segments_path.exists():
                segments_path = intermediate_dir / "segments.parquet"

            refs_path = Path(self.xp.sample_wd) / "refs.parquet"

            # Find contigs file (check for both IPC and Parquet formats)
            contigs_path = None
            for pattern in ["contigs_segmented_valid.arrow", "contigs_segmented_valid.parquet",
                            "contigs_*.arrow", "contigs_*.parquet"]:
                matches = list(Path(self.xp.sample_wd).glob(pattern))
                if matches:
                    contigs_path = matches[0]
                    break

            if segments_path.exists() and refs_path.exists():
                return StepResults(results={
                    'xp': self.xp,
                    'segments': str(segments_path),
                    'refs': str(refs_path),
                    'contigs': str(contigs_path) if contigs_path else None,
                })
        return None

    def _regenerate_extension_plots(self, ext_name: str, results: StepResults) -> None:
        """Regenerate plots for an extension."""
        from ..extensions import extension_registry

        extension = extension_registry.create_extension(ext_name, self.xp)
        if not extension:
            self.logger.warning(f"Extension {ext_name} not found")
            return

        try:
            if ext_name == 'cassiopeia_petracer':
                # Load refs into extension
                extension.refs = read_file(results.results['refs'])
                extension.workdir = Path(self.xp.sample_wd)

                # Load ldf if contigs available (needed for some QC)
                if results.results.get('contigs'):
                    from .formats import scan_file
                    extension.ldf = scan_file(results.results['contigs'])

                # Call the plot regeneration step
                extension._regenerate_filtered_qc()
                self.logger.info(f"Successfully regenerated plots for extension {ext_name}")
        except Exception as e:
            self.logger.error(f"Failed to regenerate plots for extension {ext_name}: {str(e)}", with_traceback=True)
