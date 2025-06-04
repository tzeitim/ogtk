"""Advanced analysis operations - only imported when needed."""
from typing import TYPE_CHECKING, Optional, Dict, Any
from ..metrics.analysis import query_parquet_files

if TYPE_CHECKING:
    import polars as pl

def lazy_import_polars():
    """Import polars only when needed."""
    try:
        import polars as pl
        return pl
    except ImportError:
        raise ImportError("polars required for advanced analysis. Install with: pip install polars")

class AdvancedAnalysis:
    """Advanced analysis operations using polars and parquet queries."""
    
    def __init__(self, collection):
        self.collection = collection
        self._pl: Optional['pl'] = None
    
    @property
    def pl(self):
        """Lazy polars import."""
        if self._pl is None:
            self._pl = lazy_import_polars()
        return self._pl
    
    def query_parquet_files(self, query_func, file_patterns=None, sample_ids=None, base_path=None):
        """Query parquet files for advanced analysis."""
        return query_parquet_files(self.collection, query_func, file_patterns, sample_ids, base_path)
    
    def get_efficiency_stats(self, sample_ids=None):
        """Get efficiency stats from parquet files."""
        from ..metrics.analysis import efficiency
        return efficiency(self.collection, sample_ids)