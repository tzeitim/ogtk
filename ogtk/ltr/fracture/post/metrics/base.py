from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any
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
    
    def as_df(self):
        """Convert to a DataFrame with one row."""
        import polars as pl
        return pl.DataFrame([self.as_dict()])

