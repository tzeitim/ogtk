from rich.logging import RichHandler
from rich.console import Console
from functools import wraps
import logging
import inspect
import polars as pl
import pandas as pd
from typing import Any, Optional

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
        
        handler = RichHandler(console=Console(width=250))
        formatter = logging.Formatter("%(message)s", datefmt="[%X]")
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def get_logger(self) -> CustomLogger:
        """Get the configured logger instance"""
        return self.logger

    def set_level(self, level: str) -> None:
        """Set the logger level"""
        if level not in self.levels:
            raise ValueError(f"levels supported {list(self.levels.keys())}")
        self.logger.setLevel(self.levels[level])

# Your call decorator
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

# Your helper functions
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
