"""Lightweight post-processing with optional heavy features."""

# Always available (lightweight - stdlib only)
from .core.metrics import StepMetrics, SampleMetrics
from .core.collection import PipelineMetricsCollection

# Lazy imports for heavy features
def get_advanced_analysis():
    """Get advanced analysis with polars - imported on demand."""
    from .analysis.advanced import AdvancedAnalysis
    return AdvancedAnalysis

def get_viewer_app():
    """Get Textual viewer app - imported on demand.""" 
    from .viewer.app import FractureExplorer
    return FractureExplorer

# Heavy collection with advanced features
def get_extended_collection():
    """Get extended collection with polars-based analysis capabilities."""
    from .metrics.summary import PipelineMetricsCollection as ExtendedCollection
    return ExtendedCollection

__all__ = [
    'StepMetrics',
    'SampleMetrics', 
    'PipelineMetricsCollection',
    'get_advanced_analysis',
    'get_viewer_app',
    'get_extended_collection'
]