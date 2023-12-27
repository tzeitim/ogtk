from rich.logging import RichHandler
from functools import wraps
import logging
import inspect
import polars as pl
import pandas as pd

class Rlogger:
    _instance = None

    levels = {
        "INFO": logging.INFO,
        "CRITICAL": logging.CRITICAL,
        "DEBUG": logging.DEBUG,
        "STEP": 25
    }

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Rlogger, cls).__new__(cls)
            cls._instance.setup_logger()
        return cls._instance

    def setup_logger(self):
        STEP_LEVEL_NUM = self.levels['STEP']
        
        def step(self, message, *args, **kws):
            if self.isEnabledFor(STEP_LEVEL_NUM):
                # Yes, logger takes its '*args' as 'args'.
                self._log(STEP_LEVEL_NUM, message, args, **kws) 

        logging.addLevelName(STEP_LEVEL_NUM, "STEP")
        setattr(logging.Logger, "step", step)

        self.logger = logging.getLogger("rich_logger")
        self.logger.setLevel(logging.INFO)
        handler = RichHandler()
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
    else:
        return repr(arg)


def format_value(value):
    if isinstance(value, pd.DataFrame):
        return f"pd.DataFrame(shape={value.shape}, cols={value.columns!r})"
    if isinstance(value, pl.DataFrame):
        return f"pl.DataFrame(shape={value.shape}, cols={value.columns!r})"
    else:
        return repr(value)

def call(func):
    def wrapper(*args, **kwargs):
        logger = Rlogger().get_logger()
        logger.step(f"{func.__name}")

        bound_args = inspect.signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        args_str = ',\n'.join([f"{k}={format_arg(v)}" for k, v in bound_args.arguments.items()])
        logger.debug(f"Calling [bold red]{func.__name__}[/] with args:\n{args_str}", extra={"markup": True})
        value = func(*args, **kwargs)
        logger.debug(f"{func.__name__} returned:\n{format_value(value)}")
        return value
    return wrapper
