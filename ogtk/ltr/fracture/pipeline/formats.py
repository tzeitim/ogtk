"""Format abstraction for Parquet/Arrow IPC dual support.

This module provides format-agnostic I/O operations for the fracture pipeline,
supporting both Parquet (for final outputs and archival) and Arrow IPC (for
faster intermediate file I/O).

Arrow IPC is 2-5x faster than Parquet for I/O operations, making it ideal
for intermediate files that are frequently read/written during pipeline
execution. Parquet remains preferred for final outputs due to better
compression and wider tool compatibility.
"""
import polars as pl
from pathlib import Path


def get_extension(use_ipc: bool, is_intermediate: bool = True) -> str:
    """Return file extension based on format settings.

    Parameters
    ----------
    use_ipc : bool
        Whether IPC format is enabled
    is_intermediate : bool
        Whether this is an intermediate file (True) or final output (False)

    Returns
    -------
    str
        File extension including the dot (e.g., '.arrow' or '.parquet')
    """
    if use_ipc and is_intermediate:
        return '.arrow'
    return '.parquet'


def scan_file(path: str | Path) -> pl.LazyFrame:
    """Scan file in appropriate format based on extension.

    Parameters
    ----------
    path : str | Path
        Path to the file to scan

    Returns
    -------
    pl.LazyFrame
        Lazy frame from the file
    """
    path_str = str(path)
    if path_str.endswith('.arrow'):
        return pl.scan_ipc(path_str)
    return pl.scan_parquet(path_str)


def sink_file(ldf: pl.LazyFrame, path: str | Path) -> None:
    """Sink LazyFrame in appropriate format based on extension.

    Parameters
    ----------
    ldf : pl.LazyFrame
        LazyFrame to sink
    path : str | Path
        Output path (extension determines format)
    """
    path_str = str(path)
    if path_str.endswith('.arrow'):
        ldf.sink_ipc(path_str)
    else:
        ldf.sink_parquet(path_str)


def read_file(path: str | Path) -> pl.DataFrame:
    """Read file in appropriate format based on extension.

    Parameters
    ----------
    path : str | Path
        Path to the file to read

    Returns
    -------
    pl.DataFrame
        DataFrame from the file
    """
    path_str = str(path)
    if path_str.endswith('.arrow'):
        return pl.read_ipc(path_str)
    return pl.read_parquet(path_str)


def write_file(df: pl.DataFrame, path: str | Path) -> None:
    """Write DataFrame in appropriate format based on extension.

    Parameters
    ----------
    df : pl.DataFrame
        DataFrame to write
    path : str | Path
        Output path (extension determines format)
    """
    path_str = str(path)
    if path_str.endswith('.arrow'):
        df.write_ipc(path_str)
    else:
        df.write_parquet(path_str)


def get_file_path(
    base_name: str,
    use_ipc: bool,
    is_intermediate: bool = True,
    directory: str | Path | None = None
) -> str:
    """Get file path with appropriate extension.

    Parameters
    ----------
    base_name : str
        Base filename without extension
    use_ipc : bool
        Whether IPC format is enabled
    is_intermediate : bool
        Whether this is an intermediate file (True) or final output (False)
    directory : str | Path | None
        Optional directory to prepend to the path

    Returns
    -------
    str
        Complete file path with extension
    """
    ext = get_extension(use_ipc, is_intermediate)
    filename = f"{base_name}{ext}"
    if directory:
        return str(Path(directory) / filename)
    return filename
