# Fracture Pipeline Implementation Details

Technical documentation for the fracture pipeline implementation.

## Component Overview

### Core Components

- `core.py`: Pipeline infrastructure
  - `PipelineStep`: Step definition and validation
  - `Pipeline`: Main pipeline class
  - Progress tracking
  - Dependency management

- `types.py`: Type definitions
  - `FractureXp`: Extended experiment configuration
  - Supporting type definitions

- `api_ext.py`: Polars API extensions
  - Custom DataFrame operations
  - DNA/RNA processing utilities

- `plotting.py`: Visualization utilities
  - Quality control plots
  - Progress visualization
  - Results plotting

## Pipeline Architecture

### Step Definition

Steps are defined using the `PipelineStep` enum with dictionary metadata:

```python
class PipelineStep(Enum):
    PARQUET = {
        'required_params': {'umi_len', 'rev_comp', 'modality'},
        'description': "Convert fastq reads to parquet format"
    }
    
    PREPROCESS = {
        'required_params': {'umi_len', 'anchor_ont'},
        'description': "Preprocess the loaded data"
    }
    
    FRACTURE = {
        'required_params': {'umi_len', 'anchor_ont'},
        'description': "Assemble short reads into contigs"
    }

    TEST = {
        'required_params': {'target_sample'},
        'description': "Makes a downsampled fastq file for testing purposes"
    }
```

### Step Execution

Steps are executed using the `@pipeline_step` decorator which handles:

- Parameter validation against `required_params`
- Logging and progress tracking
- Error handling and recovery
- Result collection and summary generation

```python
@pipeline_step(PipelineStep.PARQUET)
def to_parquet(self) -> StepResults|None:
    """Convert fastq reads to parquet format"""
    # Implementation here
    return StepResults(results=..., metrics=...)
```


## Adding New Steps

1. Define step in `PipelineStep` enum:
```python
NEW_STEP = {
    'required_params': {'param1', 'param2'},
    'description': "Step description"
}
```

2. Implement step function with decorator:
```python
@pipeline_step(PipelineStep.NEW_STEP)
def new_step(self) -> StepResults|None:
    """Implementation"""
    return StepResults(
        results={'output_file': 'path/to/output'},
        metrics={'metric1': value1, 'metric2': value2}
    )
```

## Development Guidelines

### Code Style

- Use type hints
- Document functions with docstrings
- Follow PEP 8
- Use dataclasses for structured data

### Testing

Add tests for:
- Step validation
- Parameter validation
- Progress tracking
- Result collection and metrics

### Error Handling

The pipeline uses Python's built-in exceptions. The `@pipeline_step` decorator automatically:
- Catches exceptions and logs them with traceback
- Validates required parameters before step execution
- Ensures proper cleanup of logging resources

```python
@pipeline_step(PipelineStep.EXAMPLE)
def example_step(self) -> StepResults|None:
    try:
        # Step implementation
        return StepResults(results=..., metrics=...)
    except Exception as e:
        # Decorator handles logging and re-raising
        raise
```

## Extending the Pipeline

### Adding Polars Extensions

The pipeline uses several custom Polars extensions defined in `api_ext.py`:

**Available namespaces:**
- `df.dna.*` - DNA sequence manipulation (reverse complement, translation, etc.)
- `df.pp.*` - Pipeline preprocessing (read parsing, UMI extraction)
- `df.ppp.*` - Post-processing utilities
- `df.kmer.*` - K-mer analysis

```python
@pl.api.register_dataframe_namespace("custom")
class PlCustom:
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df
        
    def custom_method(self):
        """Implementation"""
        pass
```

**Example usage:**
```python
# Parse reads using the pp namespace
df = pl.read_parquet("reads.parquet")
processed = df.pp.parse_reads(umi_len=12, sbc_len=6, anchor_ont="AGATCGGAAGAGC")

# Convert sequences to FASTA format
fasta_df = df.dna.to_fasta(seq_col="sequence", id_col="read_id")
```

### Adding Plot Types

Plots are defined in the `PlotDB` class in `plotting.py`. Each step can have associated plots:

```python
class PlotDB:
    @call
    def plot_new_step(self, ppi, results):
        """Add new plot type for a pipeline step"""
        sns = ppi.sns  # Access seaborn
        plt = ppi.plt  # Access matplotlib
        xp = ppi.xp    # Access experiment config
        
        # Generate plot using results.results and results.metrics
        fig = sns.scatterplot(data=results.metrics)
        fig.savefig(f"{xp.sample_figs}/{xp.target_sample}_new_analysis.png")
        xp.logger.info(f"Saved plot to {xp.sample_figs}")
```

**Available plot methods:**
- `plot_preprocess()` - Coverage plots, anchor analysis
- Custom plots can be added by following the naming convention `plot_{step_name}()`

## Performance Considerations

- **Use Polars LazyFrames** for memory-efficient processing of large datasets
- **Streaming operations** - Use `.scan_parquet()` and `.sink_parquet()` for large files
- **Parallel processing** - Polars automatically parallelizes operations
- **Memory monitoring** - Large datasets are processed in chunks
- **Caching** - Intermediate results are cached in parquet format

**Example of efficient processing:**
```python
# Memory-efficient processing of large files
result = (
    pl.scan_parquet("large_file.parquet")
    .filter(pl.col("quality") > 20)
    .group_by("umi")
    .agg(pl.col("sequence").count().alias("reads"))
    .sink_parquet("output.parquet")
)
```

## Debugging

The pipeline provides comprehensive logging at multiple levels:

**Enable debug logging:**
```bash
python pipeline.py --config config.yml --log-level DEBUG
```

**Log files are automatically created:**
- `{sample_workdir}/logs/{step}_params.log` - Parameters used for each step
- `{sample_workdir}/logs/{step}.log` - Step-specific execution logs
- `{sample_workdir}/pipeline_summary.json` - Machine-readable summary
- `{sample_workdir}/pipeline_summary.md` - Human-readable report

**Debugging features:**
- **Dry run mode** - Set `dry: true` in config to validate without execution
- **Test mode** - Use `make_test: true` to create downsampled test data
- **Step-specific execution** - Run only specific steps with `steps: ["parquet", "preprocess"]`
- **Progress tracking** - Built-in progress indicators and metrics collection

## Configuration and Data Types

### FractureXp Configuration

The pipeline uses the `FractureXp` class (extending `Xp`) for configuration:

```python
class FractureXp(Xp):
    # Core parameters
    steps: List[str]              # Pipeline steps to run
    dry: bool                     # Dry run mode
    target_sample: str            # Sample identifier
    pro_datain: str              # Input data directory
    pro_workdir: str             # Working directory
    
    # Sequencing parameters
    modality: str                 # Sequencing modality
    umi_len: int                  # UMI length
    sbc_len: int                  # Sample barcode length
    rev_comp: bool                # Reverse complement R2
    anchor_ont: str               # ONT anchor sequence
    
    # Processing parameters
    fracture: dict                # Assembly parameters
    force_tab: bool               # Force tabulation
    parse_read1: bool             # Parse R1 sequences
    allow_wildcards: bool         # Allow multiple input files
    
    # Extensions and plotting
    extensions: List[str]         # Post-processing extensions
    do_plot: bool                 # Generate plots
```

### StepResults Data Structure

Pipeline steps return `StepResults` containing:

```python
class StepResults(NamedTuple):
    results: dict    # File paths, objects, configurations
    metrics: dict    # Quantitative metrics for summary
    
    def has_metrics(self) -> bool:
        """Check if step returned metrics"""
        return len(self.metrics) > 0
```

## Pipeline Flow

1. **PARQUET** - Convert FASTQ to Parquet format
   - Input: `{pro_datain}/{target_sample}_R[12]_001.fastq.gz`
   - Output: `{pro_workdir}/{target_sample}/parsed_reads.parquet`
   - Metrics: `total_reads`, `total_umis`, `downsampled`

2. **PREPROCESS** - Parse reads and extract UMIs
   - Input: `merged_reads.parquet`
   - Output: `parsed_reads.parquet`, `parsed_reads_invalid.parquet`
   - Metrics: `fraction_valid_umis`, `n_valid_umis`, `n_invalid_umis`

3. **FRACTURE** - Assemble reads into contigs
   - Input: `parsed_reads.parquet`
   - Output: `contigs_{strategy}_{key}.parquet` where strategy is:
     - `segmented` - segment → assemble → stitch (for long cassettes with METAs)
     - `masked` - mask → assemble → unmask (for repetitive sequences)
     - `direct` - direct assembly (no masking/segmentation)
   - Metrics: `success_rate`, `total_assembled`, `mean_contig_length`

4. **Extensions** - Post-processing analysis (optional)
   - Input: `contigs_{strategy}_{key}.parquet` (auto-detected)
   - Output: Extension-specific results

## API Reference

### Core Classes
- `Pipeline` - Main pipeline orchestrator (core.py:228)
- `PipelineStep` - Step enumeration (core.py:12)
- `FractureXp` - Configuration class (types.py:8)
- `StepResults` - Results container (types.py:71)

### Decorators
- `@pipeline_step(step)` - Step execution decorator (core.py:148)
- `@call` - Logging decorator for methods

### Polars Extensions
- `df.dna.*` - DNA sequence operations
- `df.pp.*` - Pipeline preprocessing
- `df.kmer.*` - K-mer analysis

See individual module docstrings for detailed API documentation.
