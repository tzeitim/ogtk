"""Heavy analysis operations - only imported when needed."""
from typing import TYPE_CHECKING, Optional, Dict, Any

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
    """Heavy analysis operations using polars."""
    
    def __init__(self, collection):
        self.collection = collection
        self._pl: Optional['pl'] = None
    
    @property
    def pl(self):
        """Lazy polars import."""
        if self._pl is None:
            self._pl = lazy_import_polars()
        return self._pl
    
    def to_long_dataframe(self):
        """Convert to long format DataFrame (multiple rows per sample)."""
        rows = []
        
        for sample_id, sample in self.collection.samples.items():
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
            return self.pl.DataFrame(schema={
                "sample_id": self.pl.Utf8,
                "step": self.pl.Utf8,
                "metric": self.pl.Utf8,
                "value": self.pl.Object,
                "timestamp": self.pl.Datetime
            })
        
        return self.pl.DataFrame(rows)
    
    def to_wide_dataframe(self):
        """Convert to wide format DataFrame (one row per sample)."""
        if not self.collection.samples:
            return self.pl.DataFrame(schema={"sample_id": self.pl.Utf8})
        
        rows = [sample.to_dict() for sample in self.collection.samples.values()]
        return self.pl.DataFrame(rows)
    
    def calculate_read_coverage(self):
        """Calculate reads per UMI for all samples."""
        data = []
        
        for sample_id, sample in self.collection.samples.items():
            total_reads = sample.get_metric("parquet", "total_reads")
            total_umis = sample.get_metric("parquet", "total_umis")
            
            if total_reads and total_umis:
                data.append({
                    "sample_id": sample_id,
                    "total_reads": total_reads,
                    "total_umis": total_umis,
                    "reads_per_umi": total_reads / total_umis
                })
        
        return self.pl.DataFrame(data)
    
    def get_valid_umi_stats(self):
        """Calculate valid UMI percentages."""
        data = []
        
        for sample_id, sample in self.collection.samples.items():
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
        
        return self.pl.DataFrame(data).sort("valid_percentage", descending=True)