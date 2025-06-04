"""Extended collection with heavy analysis capabilities."""
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import json
from .base import SampleMetrics, StepMetrics
from ..core.collection import PipelineMetricsCollection as CoreCollection

def lazy_import_polars():
    """Import polars only when needed."""
    try:
        import polars as pl
        return pl
    except ImportError:
        raise ImportError("polars required for advanced analysis. Install with: pip install polars")

class PipelineMetricsCollection(CoreCollection):
    """Extended collection with heavy analysis capabilities."""
    
    @classmethod
    def from_summary_files(cls, file_paths=None, directory=None) -> 'PipelineMetricsCollection':
        """Create collection from summary JSON files using polars."""
        if directory and not file_paths:
            directory_path = Path(directory)
            file_paths = list(directory_path.glob("**/pipeline_summary.json"))
            if len(file_paths) == 0:
                raise ValueError("No summary files found")
        
        if not file_paths:
            raise ValueError("Either file_paths or directory must be provided")
        
        collection = cls()
        
        for file_path in file_paths:
            try:
                with open(file_path) as f:
                    data = json.load(f)
                
                sample_id = cls._extract_sample_id(file_path)
                sample_metrics = SampleMetrics(sample_id)
                
                for step_name, step_data in data.items():
                    if isinstance(step_data, dict) and 'metrics' in step_data:
                        step_metrics = StepMetrics.from_dict(step_name, step_data)
                        sample_metrics.steps[step_name] = step_metrics
                
                collection.samples[sample_id] = sample_metrics
                
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
        
        return collection
    
    def as_long_df(self):
        """Convert to long format DataFrame (multiple rows per sample)."""
        pl = lazy_import_polars()
        rows = []
        
        for sample_id, sample in self.samples.items():
            for step_name, step in sample.steps.items():
                for metric_name, value in step.metrics.items():
                    rows.append({
                        "sample_id": sample_id,
                        "step": step_name,
                        "metric": metric_name,
                        "value": value,
                        "timestamp": step.timestamp
                    })
        
        if not rows:
            return pl.DataFrame(schema={
                "sample_id": pl.Utf8,
                "step": pl.Utf8,
                "metric": pl.Utf8,
                "value": pl.Object,
                "timestamp": pl.Datetime
            })
        
        return pl.DataFrame(rows)
    
    def as_wide_df(self):
        """Convert to wide format DataFrame (one row per sample)."""
        pl = lazy_import_polars()
        
        if not self.samples:
            return pl.DataFrame(schema={"sample_id": pl.Utf8})
        
        rows = [sample.to_dict() for sample in self.samples.values()]
        return pl.DataFrame(rows)
    
    def get_metric_comparison(self, step: str, metric: str):
        """Compare a specific metric across all samples."""
        pl = lazy_import_polars()
        data = self.get_basic_comparison(step, metric)
        return pl.DataFrame(data)
    
    def calculate_read_coverage(self):
        """Calculate reads per UMI for all samples."""
        pl = lazy_import_polars()
        data = []
        
        for sample_id, sample in self.samples.items():
            total_reads = sample.get_metric("parquet", "total_reads")
            total_umis = sample.get_metric("parquet", "total_umis")
            
            if total_reads and total_umis:
                data.append({
                    "sample_id": sample_id,
                    "total_reads": total_reads,
                    "total_umis": total_umis,
                    "reads_per_umi": total_reads / total_umis
                })
        
        return pl.DataFrame(data)
    
    def get_valid_umi_stats(self):
        """Calculate valid UMI percentages."""
        pl = lazy_import_polars()
        data = []
        
        for sample_id, sample in self.samples.items():
            valid_umis = sample.get_metric("preprocess", "n_valid_umis")
            invalid_umis = sample.get_metric("preprocess", "n_invalid_umis")
            
            if valid_umis is not None and invalid_umis is not None:
                total = valid_umis + invalid_umis
                percentage = (valid_umis / total * 100) if total > 0 else 0
                
                data.append({
                    "sample_id": sample_id,
                    "valid_umis": valid_umis,
                    "invalid_umis": invalid_umis,
                    "total_umis": total,
                    "valid_percentage": percentage
                })
        
        return pl.DataFrame(data).sort("valid_percentage", descending=True)
