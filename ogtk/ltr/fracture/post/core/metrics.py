"""Lightweight data structures with no heavy dependencies."""
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from datetime import datetime
import json

@dataclass
class StepMetrics:
    """Lightweight step metrics - no external deps."""
    step_name: str
    timestamp: datetime
    metrics: Dict[str, Any]
    
    def __post_init__(self):
        if isinstance(self.timestamp, str):
            self.timestamp = datetime.fromisoformat(self.timestamp.replace(' ', 'T'))
    
    @classmethod
    def from_dict(cls, step_name: str, data: dict) -> 'StepMetrics':
        return cls(
            step_name=step_name,
            timestamp=data.get('timestamp', '2000-01-01T00:00:00'),
            metrics=data.get('metrics', {})
        )
    
    def get_metric(self, name: str, default: Any = None) -> Any:
        """Get a specific metric value."""
        return self.metrics.get(name, default)

@dataclass  
class SampleMetrics:
    """Lightweight sample metrics - stdlib only."""
    sample_id: str
    steps: Dict[str, StepMetrics] = field(default_factory=dict)
    
    def add_step(self, step_name: str, timestamp: str, metrics: Dict[str, Any]) -> None:
        """Add metrics for a pipeline step."""
        self.steps[step_name] = StepMetrics(step_name, timestamp, metrics)
    
    def get_step(self, step_name: str) -> Optional[StepMetrics]:
        """Get metrics for a specific step."""
        return self.steps.get(step_name)
    
    def get_metric(self, step: str, metric: str, default: Any = None) -> Any:
        """Get a specific metric from a specific step."""
        step_obj = self.steps.get(step)
        return step_obj.get_metric(metric, default) if step_obj else default
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to flat dictionary with all metrics."""
        result = {"sample_id": self.sample_id}
        for step in self.steps.values():
            for metric_name, value in step.metrics.items():
                result[f"{step.step_name}_{metric_name}"] = value
        return result