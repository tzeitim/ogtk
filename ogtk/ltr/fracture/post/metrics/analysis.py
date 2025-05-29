"""Analysis functions for PipelineMetricsCollection parquet file operations."""

from pathlib import Path
from typing import Dict, List, Optional, Callable, Any
import polars as pl
from .summary import PipelineMetricsCollection

def query_parquet_files(collection: PipelineMetricsCollection, query_func: Callable, 
                        file_patterns: Optional[Dict[str, str]] = None, 
                        sample_ids: Optional[List[str]] = None,
                        base_path: Optional[Path] = None) -> Dict[str, Any]:
    if file_patterns is None:
        file_patterns = {
            'valid': 'parsed_reads.parquet',
            'invalid': 'parsed_reads_invalid.parquet'
        }
    
    if sample_ids is None:
        sample_ids = list(collection.samples.keys())
    
    if base_path is None:
        base_path = Path("/home/projects/nyosef/pedro/projects/lt/workdir")
    
    results = {}
    
    for sample_id in sample_ids:
        if sample_id not in collection.samples:
            continue
            
        # Parse sample_id: format is date/sample (experiment names are always dates)
        date_dir_name, sample_dir = sample_id.split('/', 1)
        sample_path = base_path / date_dir_name / sample_dir
        
        # Build file paths based on patterns
        file_paths = {}
        if sample_path.exists():
            for pattern_name, pattern in file_patterns.items():
                file_path = sample_path / pattern
                if file_path.exists():
                    file_paths[pattern_name] = file_path
        
        if file_paths:
            try:
                result = query_func(file_paths)
                results[sample_id] = result
            except Exception as e:
                print(f"Error processing sample {sample_id}: {e}")
                results[sample_id] = None
    
    return results

def efficiency(collection: PipelineMetricsCollection, 
              sample_ids: Optional[List[str]] = None) -> Dict[str, Any]:

    def query_func(file_paths):
        lval = pl.scan_parquet(file_paths['valid'])
        linv = pl.scan_parquet(file_paths['invalid'])

        return pl.concat([lval, linv]).group_by('valid_umi', 'ont').len().collect()
    
    return query_parquet_files(collection, query_func, sample_ids=sample_ids)

def group_by_valid_umi_ont(collection: PipelineMetricsCollection, 
                      sample_ids: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Example analysis: group by valid_umi and ont columns, count rows.
    Equivalent to: pl.concat([lval, linv]).group_by('valid_umi', 'ont').len().collect()
    """
    def query_func(file_paths):
        dataframes = []
        for pattern_name, path in file_paths.items():
            df = pl.scan_parquet(path)
            dataframes.append(df)
        
        if dataframes:
            return pl.concat(dataframes).group_by('valid_umi', 'ont').len().collect()
        return None
    
    return query_parquet_files(collection, query_func, sample_ids=sample_ids)
