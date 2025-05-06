import argparse
import sys
from importlib import import_module
from typing import List

# Lazy load the heavy dependencies
__all__ = [
        'main',
        'parse_args',
]

# Define the pipeline step names for argument parsing without importing PipelineStep
PIPELINE_STEP_NAMES = ['parquet', 'preprocess', 'fracture', 'test']

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Data processing pipeline")

    parser.add_argument("--config", required=True, help="Path to config file")

    parser.add_argument(
        "--steps",
        nargs="+",
        choices=PIPELINE_STEP_NAMES,
        help="Specific steps to run (overrides config file)"
        )

    parser.add_argument(
        "--make-test",
        action="store_true",
        help="Generate test data (overrides config file setting)"
        )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["CRITICAL", "INFO", "IO", "STEP", "DEBUG"],
        help="Logging level"
        )

    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove all previous pipeline outputs before running"
    )
    parser.add_argument(
        "--clean-test",
        action="store_true",
        help="Remove only test-related outputs before running"
    )

    parser.add_argument(
        "--target-sample",
        type=str,
        help="Target sample to process (overrides config file setting)"
    )

    parser.add_argument(
        "--all-all",
        action="store_true",
        help="Run all steps for all samples. No questions asked."
    )

    parser.add_argument(
        "--all-samples",
        action="store_true",
        help="Run specified steps (conf file or argument) for all samples in the configuration"
    )

    parser.add_argument(
		"--regenerate-plots",
		nargs="*",
		help="Regenerate plots for specified steps (or all if no steps provided)"
	)
    
    parser.add_argument(
        "--regenerate-summary",
        action="store_true",
        help="Regenerate markdown summary from existing JSON data"
    )

    return parser

def lazy_import():
    """Lazily import heavier modules only when needed"""
    # Import all required modules
    modules = {}
    # First import the logger, which is needed but might still be heavy
    from ogtk.utils.log import Rlogger
    modules['Rlogger'] = Rlogger
    
    # Then import the really heavy pipeline modules
    modules['PipelineStep'] = import_module('ogtk.ltr.fracture.pipeline.core').PipelineStep
    modules['Pipeline'] = import_module('ogtk.ltr.fracture.pipeline.core').Pipeline
    modules['FractureXp'] = import_module('ogtk.ltr.fracture.pipeline.types').FractureXp
    modules['PlotRegenerator'] = import_module('ogtk.ltr.fracture.pipeline.plot_gen').PlotRegenerator
    modules['SummaryRegenerator'] = import_module('ogtk.ltr.fracture.pipeline.plot_gen').SummaryRegenerator
    return modules

def main():
    """Run the Fracture pipeline with command line arguments.

    This function serves as the main entry point for the Fracture pipeline,
    processing command line arguments and executing the appropriate pipeline steps.

    Examples
    --------
    Process all samples with all steps::

        $ python pipeline.py --config config.yml --all-all

    Process specified steps for all samples::

        $ python pipeline.py --config config.yml --all-samples

    Process specific target sample::

        $ python pipeline.py --config config.yml --target-sample "sb_rna_fracture_S3"

    Combine with other arguments::

        $ python pipeline.py --config config.yml --target-sample "sb_rna_fracture_S3" --steps parquet preprocess

    Generate test data for specific sample::

        $ python pipeline.py --config config.yml --target-sample "sb_rna_fracture_S3" --make-test

    Clean and process specific sample::

        $ python pipeline.py --config config.yml --target-sample "sb_rna_fracture_S3" --clean

    Regenerate plots for all steps::

        $ python pipeline.py --config config.yml --regenerate-plots

    Regenerate plots for specific steps::

        $ python pipeline.py --config config.yml --regenerate-plots preprocess fracture

    Regenerate plots for a specific sample::

        $ python pipeline.py --config config.yml --target-sample "sample_id" --regenerate-plots
        
    Regenerate markdown summary from existing JSON data::
    
        $ python pipeline.py --config config.yml --target-sample "sample_id" --regenerate-summary
        
    Regenerate both plots and summary::
    
        $ python pipeline.py --config config.yml --target-sample "sample_id" --regenerate-plots --regenerate-summary

    Returns
    -------
    int
        0 for successful execution, 1 for failure
    """
    # Handle --help or no arguments first, before creating the parser
    # This avoids importing any heavy modules for help display
    if len(sys.argv) == 1 or '--help' in sys.argv or '-h' in sys.argv:
        parser = parse_args()
        parser.print_help()
        return 0
    
    # Parse arguments
    parser = parse_args()
    args = parser.parse_args()

    try:
        # Only import modules when we're actually running the pipeline
        modules = lazy_import()
        Rlogger = modules['Rlogger']
        
        # Initialize logger
        logger = Rlogger().get_logger()
        Rlogger().set_level(args.log_level)
        
        # Reference the heavy modules
        logger.debug("Loading pipeline modules...")
        PipelineStep = modules['PipelineStep']
        Pipeline = modules['Pipeline']
        FractureXp = modules['FractureXp']
        PlotRegenerator = modules['PlotRegenerator']
        SummaryRegenerator = modules['SummaryRegenerator']
        
        # Initialize Xp configuration
        xp = FractureXp(conf_fn=args.config)

        # Skip the pipeline if we are only regenerating plots or summary
        if args.regenerate_plots is not None or args.regenerate_summary:
            pipeline = Pipeline(xp)
            
            # Handle plot regeneration
            if args.regenerate_plots is not None:
                logger.info("Regenerating plots")
                plot_regenerator = PlotRegenerator(pipeline)
                plot_regenerator.regenerate_plots(steps=args.regenerate_plots)
            
            # Handle summary regeneration
            if args.regenerate_summary:
                logger.info("Regenerating markdown summary")
                summary_regenerator = SummaryRegenerator(pipeline)
                summary_regenerator.regenerate_summary()
                
            return 0

        # Handle --all-all flag (all samples, all steps)
        if args.all_all:
            logger.info("Processing all steps for all samples")
            # Set steps to include all available steps
            xp.steps = [step.name.lower() for step in PipelineStep]
            success = True

            # Process each sample
            for sample in xp.samples:
                logger.info(f"\nProcessing sample: {sample['id']}")
                xp.target_sample = sample['id']
                xp.consolidate_conf(update=True)
                pipeline = Pipeline(xp)

                # Handle cleaning if requested
                if args.clean:
                    pipeline.clean_all()
                elif args.clean_test:
                    pipeline.clean_test_outputs()

                sample_success = pipeline.run()
                success = success and sample_success

            return 0 if success else 1

        # Handle --all-samples flag
        elif args.all_samples:
            logger.info("Processing specified steps for all samples")
            # Use steps from command line if provided, otherwise use config
            if args.steps:
                xp.steps = args.steps

            success = True

            # Process each sample
            for sample in xp.samples:
                logger.info(f"\nProcessing sample: {sample['id']}")
                xp.target_sample = sample['id']
                xp.consolidate_conf(update=True)
                pipeline = Pipeline(xp)

                # Handle cleaning if requested
                if args.clean:
                    pipeline.clean_all()
                elif args.clean_test:
                    pipeline.clean_test_outputs()

                sample_success = pipeline.run()
                success = success and sample_success

            return 0 if success else 1

        # Handle single sample processing
        if args.target_sample:
            logger.info(f"Overriding target sample from config with: {args.target_sample}")
            # Verify sample exists in samples list
            sample_ids = [sample['id'] for sample in xp.samples]
            if args.target_sample not in sample_ids:
                raise ValueError(f"Target sample '{args.target_sample}' not found in samples list: {sample_ids}")

            xp.target_sample = args.target_sample
            xp.consolidate_conf(update=True)

        pipeline = Pipeline(xp)

        # Handle cleaning if requested
        if args.clean:
            pipeline.clean_all()
        elif args.clean_test:
            pipeline.clean_test_outputs()

        # Override make_test if specified in command line
        if args.make_test:
            xp.make_test = True
            if "test" not in [step.lower() for step in xp.steps]:
                xp.steps.append("test")

        # Override steps if specified in command line
        if args.steps:
            xp.steps = args.steps

        success = pipeline.run()

        return 0 if success else 1

    except Exception as e:
        # If it's just a command line error or config error, don't show traceback
        if isinstance(e, (ValueError, FileNotFoundError, argparse.ArgumentError)):
            logger.error(f"Error: {str(e)}")
        else:
            logger.error(f"Pipeline execution failed: {str(e)}", with_traceback=True)
        return 1

if __name__ == "__main__":
    exit(main())
