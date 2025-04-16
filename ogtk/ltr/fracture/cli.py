import argparse
from ogtk.utils.log import Rlogger

from ogtk.ltr.fracture.pipeline.core import PipelineStep, Pipeline
from ogtk.ltr.fracture.pipeline.types import FractureXp
from ogtk.ltr.fracture.pipeline.plot_gen import PlotRegenerator

__all__ = [
        'main',
        'parse_args',
]

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Data processing pipeline")

    parser.add_argument("--config", required=True, help="Path to config file")

    parser.add_argument(
        "--steps",
        nargs="+",
        choices=[name.lower() for name in PipelineStep.__members__],
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

    return parser

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

    Returns
    -------
    int
        0 for successful execution, 1 for failure
    """
    # Parse arguments
    args = parse_args().parse_args()

    # Initialize logger
    logger = Rlogger().get_logger()
    Rlogger().set_level(args.log_level)

    try:
        # Initialize Xp configuration
        xp = FractureXp(conf_fn=args.config)

		# Skip the pipeline if we are only plotting
        if args.regenerate_plots is not None:
            logger.info("Regenerating plots")
            pipeline = Pipeline(xp)
            plot_regenerator = PlotRegenerator(pipeline)
            plot_regenerator.regenerate_plots(steps=args.regenerate_plots)
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
        logger.error(f"Pipeline execution failed: {str(e)}")
        return 1

if __name__ == "__main__":
    exit(main())
