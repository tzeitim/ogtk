from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any
from datetime import datetime
from .base import SampleMetrics

@dataclass
class PipelineMetricsCollection:
    """Collection of metrics for multiple samples."""
    samples: Dict[str, SampleMetrics] = field(default_factory=dict)
    source_info: Dict[str, Any] = field(default_factory=dict)
        
    @classmethod
    def from_summary_files(cls, file_paths=None, directory=None) -> 'PipelineMetricsCollection':
        """Create collection from summary JSON files."""
        # Collect JSON files
        if directory and not file_paths:
            directory_path = Path(directory)
            file_paths = list(directory_path.glob("**/pipeline_summary.json"))
            if len(file_paths) == 0:
                raise ValueError("No summary files found")
        
        if not file_paths:
            raise ValueError("Either file_paths or directory must be provided")
        
        collection = cls()
        
        for file_path in file_paths:
            sample_dir = file_path.parent.name
            experiment_dir = file_path.parent.parent.name
            sample_id = f"{experiment_dir}/{sample_dir}"
            
            try:
                # Read JSON file
                import polars as pl
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

    def __repr__(self) -> str:
        n_samples = len(self.samples)
        if n_samples == 0:
            return "PipelineMetricsCollection(empty)"
        
        sample_ids = list(self.samples.keys())[:3]
        steps = set()
        for sample in self.samples.values():
            steps.update(sample.steps.keys())
        
        sample_preview = ", ".join(sample_ids)
        if n_samples > 3:
            sample_preview += f", ... (+{n_samples-3} more)"
        
        return f"PipelineMetricsCollection({n_samples} samples: {sample_preview}; steps: {sorted(steps)})"

    
    def get_sample(self, sample_id: str) -> Optional[SampleMetrics]:
        """Get metrics for a specific sample."""
        return self.samples.get(sample_id)
    
    def as_long_df(self):
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
        
        import polars as pl
        
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
        import polars as pl
        
        if not self.samples:
            return pl.DataFrame(schema={"sample_id": pl.Utf8})
        
        return pl.concat([sample.as_df() for sample in self.samples.values()])
    
    def get_metric_comparison(self, step: str, metric: str):
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
        
        import polars as pl
        return pl.DataFrame(data)
    
    def calculate_read_coverage(self):
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
        
        import polars as pl
        return pl.DataFrame(data)
    
    def get_valid_umi_stats(self):
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
        
        import polars as pl
        return pl.DataFrame(data).sort("valid_percentage", descending=True)
