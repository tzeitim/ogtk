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

# Legacy compatibility (will import heavy deps)
def get_legacy_collection():
    """Get the old heavy collection for backward compatibility."""
    import warnings
    warnings.warn(
        "Legacy collection imports heavy dependencies. "
        "Consider using the new lightweight PipelineMetricsCollection.",
        DeprecationWarning,
        stacklevel=2
    )
    from ..post.metrics.summary import PipelineMetricsCollection as LegacyCollection
    return LegacyCollection

__all__ = [
    'StepMetrics',
    'SampleMetrics', 
    'PipelineMetricsCollection',
    'get_advanced_analysis',
    'get_viewer_app',
    'get_legacy_collection'
]