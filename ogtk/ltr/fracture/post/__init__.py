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

def get_tree_qc():
    """Get tree QC functions - imported on demand."""
    from .tree_qc import (
        tree_parsimony_qc,
        batch_tree_qc,
        compute_parsimony_score,
        generate_null_distribution,
        compute_qc_stats,
        qc_result_to_dataframe,
        get_tree_dirs,
        save_qc_result,
        collect_qc_results,
        TreeQCResult,
        ColumnQCResult,
    )
    return {
        'tree_parsimony_qc': tree_parsimony_qc,
        'batch_tree_qc': batch_tree_qc,
        'compute_parsimony_score': compute_parsimony_score,
        'generate_null_distribution': generate_null_distribution,
        'compute_qc_stats': compute_qc_stats,
        'qc_result_to_dataframe': qc_result_to_dataframe,
        'get_tree_dirs': get_tree_dirs,
        'save_qc_result': save_qc_result,
        'collect_qc_results': collect_qc_results,
        'TreeQCResult': TreeQCResult,
        'ColumnQCResult': ColumnQCResult,
    }

__all__ = [
    'StepMetrics',
    'SampleMetrics',
    'PipelineMetricsCollection',
    'get_advanced_analysis',
    'get_viewer_app',
    'get_extended_collection',
    'get_tree_qc',
]