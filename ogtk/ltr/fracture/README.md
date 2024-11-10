# Fracture Pipeline

FRagmentation and Assembly for Comprehensive Target Unit REconstruction

A pipeline for processing and assembling molecule-based sequencing data. This pipeline provides functionality for converting FASTQ files to parquet format, preprocessing data, and assembling contigs.

## Quick Start

```bash
# Set your system prefix
export OGTK_SYSPREFIX="/path/to/your/workspace"

# Run the pipeline
python pipeline.py --config config.yml --target-sample "your_sample_id"
```

## Directory Structure

```
fracture/
├── README.md           # This file
├── config.yml         # Configuration template
├── pipeline.py        # Main entry point
└── pipeline/          # Core pipeline implementation
    ├── api_ext.py    # Polars API extensions
    ├── core.py       # Pipeline core functionality
    ├── plotting.py   # Plotting utilities
    ├── steps.py      # Pipeline step definitions
    └── types.py      # Type definitions
```

## Installation

1. Set up your environment:
```bash
conda create -n fracture python=3.10
conda activate fracture
pip install -r requirements.txt
```

2. Configure your system paths:
```yaml
# config.yml
system:
  prefixes:
    your_system: "/your/path/"
  default: "your_system"
```

## Basic Usage

### Running the Pipeline

```bash
# Process a specific sample
python pipeline.py --config config.yml --target-sample "sample_id"

# Run specific steps
python pipeline.py --config config.yml --steps parquet preprocess

# Generate test data
python pipeline.py --config config.yml --make-test
```

### Available Steps

1. **parquet**: Convert FASTQ to parquet format
2. **preprocess**: Parse reads and collect molecule stats
3. **fracture**: Assemble reads into contigs
4. **test**: Generate test data

## Configuration

See `config.yml` for a template configuration file. Key parameters:

```yaml
project: "20241020"
modality: "single-molecule"
umi_len: 12
target_sample: "sample_id"
samples:
  - id: "sample_id"
```

## Environment Variables

- `OGTK_SYSPREFIX`: Set the system prefix path
```bash
export OGTK_SYSPREFIX="/path/to/workspace"
```

## Advanced Features

- Progress tracking with rich progress bars
- Automatic dependency management
- Configurable plotting
- Test data generation
- Pipeline visualization

## Contributing

Please read the [pipeline/README.md](pipeline/README.md) for technical details about the implementation.

## Troubleshooting

Common issues and solutions:

1. **Missing Files**
   ```
   Error: Missing required files for step PARQUET
   ```
   - Check file paths in configuration
   - Verify OGTK_SYSPREFIX setting

2. **Configuration Issues**
   ```
   Error: Missing required parameter: umi_len
   ```
   - Review config.yml
   - Check parameter requirements for each step


