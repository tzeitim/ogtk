"""Collection class with minimal dependencies."""
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import json
from .metrics import SampleMetrics, StepMetrics

class PipelineMetricsCollection:
    """Lightweight collection - only depends on stdlib + core metrics."""
    
    def __init__(self):
        self.samples: Dict[str, SampleMetrics] = {}
        self.source_info: Dict[str, Any] = {}
    
    @classmethod
    def from_summary_files(cls, file_paths: List[Path]) -> 'PipelineMetricsCollection':
        """Load from JSON files - no heavy dependencies."""
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
    
    @staticmethod
    def _extract_sample_id(file_path: Path) -> str:
        """Extract sample ID from file path."""
        sample_dir = file_path.parent.name
        experiment_dir = file_path.parent.parent.name
        return f"{experiment_dir}/{sample_dir}"
    
    def get_sample(self, sample_id: str) -> Optional[SampleMetrics]:
        """Get metrics for a specific sample."""
        return self.samples.get(sample_id)
    
    def to_simple_dict(self) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """Export to simple dict structure for lightweight processing."""
        result = {}
        for sample_id, sample in self.samples.items():
            result[sample_id] = {}
            for step_name, step in sample.steps.items():
                result[sample_id][step_name] = step.metrics
        return result
    
    def get_basic_comparison(self, step: str, metric: str) -> Dict[str, Any]:
        """Basic comparison without heavy dependencies."""
        result = {"sample_id": [], f"{step}_{metric}": []}
        
        for sample_id, sample in self.samples.items():
            value = sample.get_metric(step, metric)
            if value is not None:
                result["sample_id"].append(sample_id)
                result[f"{step}_{metric}"].append(value)
        
        return result
    
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
