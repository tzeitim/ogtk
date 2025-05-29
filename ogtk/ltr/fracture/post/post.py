from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any
import polars as pl
from datetime import datetime

@dataclass
class StepMetrics:
    """Class representing metrics for a single pipeline step."""
    step_name: str
    timestamp: datetime
    metrics: Dict[str, Any]
    
    def __post_init__(self):
        # Convert timestamp string to datetime if needed
        if isinstance(self.timestamp, str):
            self.timestamp = datetime.strptime(self.timestamp, "%Y-%m-%d %H:%M:%S")
    
    def get_metric(self, name: str, default: Any = None) -> Any:
        """Get a specific metric value."""
        return self.metrics.get(name, default)
    
    def as_dict(self) -> Dict[str, Any]:
        """Convert to a flat dictionary with prefixed keys."""
        result = {}
        for k, v in self.metrics.items():
            result[f"{self.step_name}_{k}"] = v
        return result

@dataclass
class SampleMetrics:
    """Class representing all metrics for a sample across pipeline steps."""
    sample_id: str
    steps: Dict[str, StepMetrics] = field(default_factory=dict)
    
    def add_step(self, step_name: str, timestamp: str, metrics: Dict[str, Any]) -> None:
        """Add metrics for a pipeline step."""
        self.steps[step_name] = StepMetrics(step_name, timestamp, metrics)
    
    def get_step(self, step_name: str) -> Optional[StepMetrics]:
        """Get metrics for a specific step."""
        return self.steps.get(step_name)
    
    def get_metric(self, step_name: str, metric_name: str, default: Any = None) -> Any:
        """Get a specific metric from a specific step."""
        step = self.get_step(step_name)
        if step:
            return step.get_metric(metric_name, default)
        return default
    
    def as_dict(self) -> Dict[str, Any]:
        """Convert to a flat dictionary with all metrics."""
        result = {"sample_id": self.sample_id}
        for step in self.steps.values():
            result.update(step.as_dict())
        return result
    
    def as_df(self) -> pl.DataFrame:
        """Convert to a DataFrame with one row."""
        return pl.DataFrame([self.as_dict()])

@dataclass
class PipelineMetricsCollection:
    """Collection of metrics for multiple samples."""
    samples: Dict[str, SampleMetrics] = field(default_factory=dict)
    
    @classmethod
    def from_summary_files(cls, file_paths=None, directory=None) -> 'PipelineMetricsCollection':
        """Create collection from summary JSON files."""
        # Collect JSON files
        if directory and not file_paths:
            directory_path = Path(directory)
            file_paths = list(directory_path.glob("*/pipeline_summary.json"))
        
        if not file_paths:
            raise ValueError("Either file_paths or directory must be provided")
        
        collection = cls()
        
        for file_path in file_paths:
            # Extract sample ID from path
            sample_id = file_path.parent.name
            
            try:
                # Read JSON file
                df = pl.read_json(file_path)
                
                sample_metrics = SampleMetrics(sample_id)
                
                for step in df.columns:
                    step_data = df.select(pl.col(step)).to_dicts()[0][step]
                    if step_data and 'metrics' in step_data:
                        timestamp = step_data.get('timestamp', '2000-01-01 00:00:00')
                        sample_metrics.add_step(step, timestamp, step_data['metrics'])
                
                collection.samples[sample_id] = sample_metrics
            
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        return collection
    
    def get_sample(self, sample_id: str) -> Optional[SampleMetrics]:
        """Get metrics for a specific sample."""
        return self.samples.get(sample_id)
    
    def as_long_df(self) -> pl.DataFrame:
        """Convert to long format DataFrame (multiple rows per sample)."""
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
    
    def as_wide_df(self) -> pl.DataFrame:
        """Convert to wide format DataFrame (one row per sample)."""
        if not self.samples:
            return pl.DataFrame(schema={"sample_id": pl.Utf8})
        
        return pl.concat([sample.as_df() for sample in self.samples.values()])
    
    def get_metric_comparison(self, step: str, metric: str) -> pl.DataFrame:
        """Compare a specific metric across all samples."""
        data = {
            "sample_id": [],
            f"{step}_{metric}": []
        }
        
        for sample_id, sample in self.samples.items():
            value = sample.get_metric(step, metric)
            if value is not None:
                data["sample_id"].append(sample_id)
                data[f"{step}_{metric}"].append(value)
        
        return pl.DataFrame(data)
    
    def calculate_read_coverage(self) -> pl.DataFrame:
        """Calculate reads per UMI for all samples."""
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
    
    def get_valid_umi_stats(self) -> pl.DataFrame:
        """Calculate valid UMI percentages."""
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
