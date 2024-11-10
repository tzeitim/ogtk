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

        if 'make_test' in xp.steps:
            xp.target_sample = f'TEST_{xp.target_sample}'
        # Initialize and run pipeline
        pipeline = Pipeline(xp)
        success = pipeline.run()
        
        return 0 if success else 1
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        return 1

if __name__ == "__main__":
    exit(main())
