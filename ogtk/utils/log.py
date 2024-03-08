from rich.logging import RichHandler
from rich.console import Console
from functools import wraps
import logging
import inspect
import polars as pl
import pandas as pd

# Define a custom logger class extending the standard Logger class
class CustomLogger(logging.Logger):
    def io(self, message, *args, **kws)->None:
        IO_LEVEL_NUM = 19
        if self.isEnabledFor(IO_LEVEL_NUM):
            self._log(IO_LEVEL_NUM, message, args, **kws)
    
    def step(self, message, *args, **kws)->None:
        STEP_LEVEL_NUM = 25
        if self.isEnabledFor(STEP_LEVEL_NUM):
            self._log(STEP_LEVEL_NUM, message, args, **kws)

# rich logger
class Rlogger:
    _instance = None

    levels = {
        "CRITICAL": logging.CRITICAL, #50
        "INFO": logging.INFO, #20
        "IO": 19,
        "STEP": 18,
        "DEBUG": logging.DEBUG, #10
    }

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Rlogger, cls).__new__(cls)
            cls._instance.setup_logger()
        return cls._instance

    def setup_logger(self):
        # Set the logger class to our custom logger
        logging.setLoggerClass(CustomLogger)

        # Add custom levels' names
        logging.addLevelName(self.levels['STEP'], "STEP")
        logging.addLevelName(self.levels['IO'], "IO")

        # Initialize the logger with the custom class
        self.logger = logging.getLogger("rich_logger")
        self.logger.setLevel(logging.INFO)
        handler = RichHandler(console=Console(width=250))
        formatter = logging.Formatter("%(message)s", datefmt="[%X]")
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)


    def get_logger(self):
        return self.logger

    def set_level(self, level):
        if level not in self.levels:
            raise ValueError(f"levels supported {list(self.levels.keys())}")
        self.logger.setLevel(self.levels[level])

def format_arg(arg):
    if isinstance(arg, pd.DataFrame):
        return f"pd.DataFrame(shape={arg.shape}, cols={arg.columns!r})"
    if isinstance(arg, pl.DataFrame):
        return f"pl.DataFrame(shape={arg.shape}, cols={arg.columns!r})"
    return is_iterable(arg)

def is_iterable(obj):
    try:
        iter(obj)
        if isinstance(obj, pl.Series):
            return f"pl.Series(shape={obj.shape}, values={obj.head(5)!r})"
        return repr(obj)[0:200]
    except TypeError:
        return repr(obj)

def format_value(value):
    if isinstance(value, pd.DataFrame):
        return f"pd.DataFrame(shape={value.shape}, cols={value.columns!r})"
    if isinstance(value, pl.DataFrame):
        return f"pl.DataFrame(shape={value.shape}, cols={value.columns!r})"
    return is_iterable(value)

def call(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
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
