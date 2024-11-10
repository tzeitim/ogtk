import argparse
from ogtk.utils.log import Rlogger

from pipeline.core import PipelineStep, Pipeline
from pipeline.types import FractureXp


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
