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

Steps are defined using the `PipelineStep` enum with `StepIO` metadata:

```python
class PipelineStep(Enum):
    PARQUET = StepIO(
        required_files=[
            "${pro_datain}/{target_sample}_R1_001.fastq.gz",
            "${pro_datain}/{target_sample}_R2_001.fastq.gz"
        ],
        output_files=[
            "${pro_workdir}/{target_sample}/parsed_reads.parquet"
        ],
        description="Convert fastq reads to parquet format",
        required_params={'umi_len', 'rev_comp', 'modality'}
    )
```

### Dependency Management

Dependencies are automatically managed through file IO specifications:

```python
def get_step_dependencies(self) -> Dict[str, Set[str]]:
    """Get dependencies between pipeline steps based on file IO"""
    dependencies = {}
    step_outputs = {}
    
    # Collect all step outputs
    for step in PipelineStep:
        _, outputs = step.format_paths(self.xp)
        for output in outputs:
            step_outputs[output] = step.name
```


## Adding New Steps

1. Define step in `PipelineStep` enum:
```python
NEW_STEP = StepIO(
    required_files=[...],
    output_files=[...],
    description="Step description",
    required_params={...}
)
```

2. Implement step function with decorator:
```python
@pipeline_step(PipelineStep.NEW_STEP)
def new_step(self) -> None:
    """Implementation"""
    pass
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
- File dependencies
- Parameter validation
- Progress tracking

### Error Handling

Use custom exceptions:
```python
class PipelineError(Exception):
    """Base class for pipeline exceptions"""
    pass

class ValidationError(PipelineError):
    """Validation failed"""
    pass
```

## Extending the Pipeline

### Adding Polars Extensions

```python
@pl.api.register_dataframe_namespace("custom")
class PlCustom:
    def __init__(self, df: pl.DataFrame) -> None:
        self._df = df
        
    def custom_method(self):
        """Implementation"""
        pass
```

### Adding Plot Types

```python
@call
def plot_new_analysis(self, ppi, results):
    """Add new plot type"""
    pass
```

## Performance Considerations

- Use Polars for data processing
- Implement progress callbacks
- Enable parallel processing where possible
- Monitor memory usage

## Debugging

Enable debug logging:
```bash
python pipeline.py --config config.yml --log-level DEBUG
```

## API Reference

See individual module docstrings for detailed API documentation.
