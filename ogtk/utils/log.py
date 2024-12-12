from rich.logging import RichHandler
from rich.console import Console
from functools import wraps
import logging
import inspect
import polars as pl
import pandas as pd
from typing import Any, Optional
from pathlib import Path

# Define custom log levels
IO_LEVEL_NUM = 19
STEP_LEVEL_NUM = 18

# Define a custom logger class extending the standard Logger class
class CustomLogger(logging.Logger):
    """Custom logger with additional IO and STEP levels"""
    
    def io(self, message: Any, *args: Any, **kws: Any) -> None:
        """Log message at IO level"""
        if self.isEnabledFor(IO_LEVEL_NUM):
            self._log(IO_LEVEL_NUM, message, args, **kws)
    
    def step(self, message: Any, *args: Any, **kws: Any) -> None:
        """Log message at STEP level"""
        if self.isEnabledFor(STEP_LEVEL_NUM):
            self._log(STEP_LEVEL_NUM, message, args, **kws)

class Rlogger:
    _instance: Optional['Rlogger'] = None
    logger: CustomLogger
    file_handler: Optional[logging.FileHandler] = None

    levels = {
        "CRITICAL": logging.CRITICAL, #50
        "INFO": logging.INFO, #20
        "IO": IO_LEVEL_NUM,
        "STEP": STEP_LEVEL_NUM,
        "DEBUG": logging.DEBUG, #10
    }

    def __new__(cls) -> 'Rlogger':
        if cls._instance is None:
            cls._instance = super(Rlogger, cls).__new__(cls)
            cls._instance.setup_logger()
        return cls._instance

    def setup_logger(self) -> None:
        # Register custom levels
        logging.addLevelName(self.levels['STEP'], "STEP")
        logging.addLevelName(self.levels['IO'], "IO")

        # Set the logger class to our custom logger
        logging.setLoggerClass(CustomLogger)

        # Initialize the logger with the custom class
        self.logger = logging.getLogger("rich_logger")  # type: ignore
        self.logger.setLevel(logging.INFO)
        
        # Console handler setup
        console_handler = RichHandler(console=Console(width=250))
        console_formatter = logging.Formatter("%(message)s", datefmt="[%X]")
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)

    def enable_file_logging(self, 
                          filepath: str | Path, 
                          level: str = "INFO", 
                          mode: str = 'a',
                          format_string: str = "%(asctime)s - %(levelname)s - %(message)s") -> None:
        """
        Enable logging to a file with customizable settings.
        
        Args:
            filepath: Path to the log file
            level: Logging level for the file output
            mode: File opening mode ('w' for write, 'a' for append)
            format_string: Format string for log messages
        """
        # Create directory if it doesn't exist
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        
        # Remove existing file handler if present
        if self.file_handler is not None:
            self.logger.removeHandler(self.file_handler)
            self.file_handler.close()

        # Create new file handler
        self.file_handler = logging.FileHandler(filepath, mode=mode)
        self.file_handler.setLevel(self.levels.get(level, logging.INFO))
        
        # Create formatter and add it to the handler
        file_formatter = logging.Formatter(format_string)
        self.file_handler.setFormatter(file_formatter)
        
        # Add the handler to the logger
        self.logger.addHandler(self.file_handler)
        self.logger.debug(f"File logging enabled: {filepath}")

    def disable_file_logging(self) -> None:
        """Disable file logging and close the file handler."""
        if self.file_handler is not None:
            self.logger.debug("File logging disabled")
            self.logger.removeHandler(self.file_handler)
            self.file_handler.close()
            self.file_handler = None

    def get_logger(self) -> CustomLogger:
        """Get the configured logger instance"""
        return self.logger

    def set_level(self, level: str) -> None:
        """Set the logger level"""
        if level not in self.levels:
            raise ValueError(f"levels supported {list(self.levels.keys())}")
        self.logger.setLevel(self.levels[level])

# Your call decorator and helper functions remain unchanged
def call(func):
    """Decorator to log function calls with arguments"""
    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        logger = Rlogger().get_logger()
        logger.step(f"{func.__name__}")

        bound_args = inspect.signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        args_str = ',\n'.join([f"{k}={format_arg(v)}" for k, v in bound_args.arguments.items()])
        logger.debug(f"Calling [bold red]{func.__name__}[/] with args:\n{args_str}", extra={"markup": True})
        value = func(*args, **kwargs)
        logger.debug(f"{func.__name__} returned:\n{format_value(value)}")
        return value
    return wrapper

def format_arg(arg: Any) -> str:
    """Format argument for logging"""
    if isinstance(arg, pd.DataFrame):
        return f"pd.DataFrame(shape={arg.shape}, cols={arg.columns!r})"
    if isinstance(arg, pl.DataFrame):
        return f"pl.DataFrame(shape={arg.shape}, cols={arg.columns!r})"
    return is_iterable(arg)

def is_iterable(obj: Any) -> str:
    """Check if object is iterable and format accordingly"""
    try:
        iter(obj)
        if isinstance(obj, pl.Series):
            return f"pl.Series(shape={obj.shape}, values={obj.head(5)!r})"
        return repr(obj)[0:200]
    except TypeError:
        return repr(obj)

def format_value(value: Any) -> str:
    """Format return value for logging"""
    if isinstance(value, pd.DataFrame):
        return f"pd.DataFrame(shape={value.shape}, cols={value.columns!r})"
    if isinstance(value, pl.DataFrame):
        return f"pl.DataFrame(shape={value.shape}, cols={value.columns!r})"
    return is_iterable(value)
