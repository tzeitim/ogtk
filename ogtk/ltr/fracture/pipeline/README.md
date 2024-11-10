# Fracture Pipeline Tutorial

This tutorial guides you through using the fracture pipeline for processing and assembling sequencing data. The pipeline provides functionality for converting FASTQ files to parquet format, preprocessing data, and assembling contigs.

## Table of Contents
1. [Installation](#installation)
2. [Pipeline Configuration](#pipeline-configuration)
3. [Basic Usage](#basic-usage)
4. [Pipeline Steps](#pipeline-steps)
5. [Sample Management](#sample-management)
6. [Advanced Usage](#advanced-usage)
7. [Troubleshooting](#troubleshooting)

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd pipeline
```

2. Set up the environment:
```bash
conda create -n fracture python=3.10
conda activate fracture
pip install -r requirements.txt
```

## Pipeline Configuration

The pipeline uses YAML configuration files to specify parameters and paths. Here's a basic configuration template:

```yaml
system:
  prefixes:
    yoseflab1: "/home/pedro/wexac/"
    wexac: "/home/labs/nyosef/pedro/"
    mac: "/Volumes/pedro/"
  default: "wexac"  # Default prefix to use

project: "20241020" # always start with date
xp_template: "${prefix}/projects/lt/conf/fracture_xp_template.yml"
modality: "single-molecule"
umi_len: 12

target_sample: "sb_rna_fracture_S3"  # Default target sample
samples:  # List of all available samples
  - id: "sb_rna_fracture_S3"
  - id: "s4_dna_fracture_S5"
  - id: "sc_rna_fracture_S4"

dry: False  # Set to True for dry run
do_plot: True  # Enable plotting

steps:
  - parquet
  - preprocess
  - fracture
```

### Key Configuration Parameters

- `system.prefixes`: Define different system paths for different environments
- `system.default`: Select which prefix to use
- `project`: Project identifier (use date format YYYYMMDD)
- `xp_template`: Path to experiment template file
- `modality`: Data modality (e.g., "single-molecule")
- `umi_len`: Length of UMI sequences
- `target_sample`: Default sample to process
- `samples`: List of all samples in the project
- `dry`: Enable dry run mode (no file writing)
- `do_plot`: Enable generation of plots
- `steps`: List of pipeline steps to execute

## Basic Usage

### Running the Complete Pipeline

```bash
# Process default target sample
python pipeline.py --config config.yml

# Process specific sample
python pipeline.py --config config.yml --target-sample "s4_dna_fracture_S5"
```

### Running Specific Steps

```bash
# Run steps on default sample
python pipeline.py --config config.yml --steps parquet preprocess

# Run steps on specific sample
python pipeline.py --config config.yml --target-sample "sc_rna_fracture_S4" --steps parquet preprocess
```

### Generating Test Data

```bash
# Generate test data for default sample
python pipeline.py --config config.yml --make-test

# Generate test data for specific sample
python pipeline.py --config config.yml --target-sample "s4_dna_fracture_S5" --make-test
```

### Cleaning Previous Outputs

```bash
# Clean all outputs
python pipeline.py --config config.yml --clean

# Clean only test outputs
python pipeline.py --config config.yml --clean-test

# Clean and process specific sample
python pipeline.py --config config.yml --clean --target-sample "sc_rna_fracture_S4"
```

### Setting Log Level

```bash
python pipeline.py --config config.yml --log-level DEBUG
```

## Pipeline Steps

### 1. Parquet Conversion (`parquet`)
Converts FASTQ files to parquet format for efficient processing.

Required parameters:
- `umi_len`
- `rev_comp`
- `modality`

```bash
python pipeline.py --config config.yml --steps parquet
```

### 2. Preprocessing (`preprocess`)
Preprocesses the data, including read parsing and molecule statistics.

Required parameters:
- `umi_len`
- `anchor_ont`

```bash
python pipeline.py --config config.yml --steps preprocess
```

### 3. Fracture Assembly (`fracture`)
Assembles short reads into contigs.

Required parameters:
- `umi_len`
- `anchor_ont`

```bash
python pipeline.py --config config.yml --steps fracture
```

### 4. Test Data Generation (`test`)
Creates downsampled data for testing.

Required parameters:
- `target_sample`

```bash
python pipeline.py --config config.yml --make-test
```

## Sample Management

### Working with Multiple Samples

The pipeline supports processing multiple samples defined in the config file. You can:

1. Define samples in config:
```yaml
samples:
  - id: "sb_rna_fracture_S3"
  - id: "s4_dna_fracture_S5"
  - id: "sc_rna_fracture_S4"
```

2. Select target sample via command line:
```bash
python pipeline.py --config config.yml --target-sample "s4_dna_fracture_S5"
```

### Sample Validation

The pipeline validates that the specified target sample exists in the samples list:
- Prevents processing undefined samples
- Ensures consistent sample naming
- Maintains data organization

### Sample-specific Operations

```bash
# Generate test data for specific sample
python pipeline.py --config config.yml --target-sample "s4_dna_fracture_S5" --make-test --clean-test

# Process specific sample with custom steps
python pipeline.py --config config.yml --target-sample "sc_rna_fracture_S4" --steps parquet preprocess

# Clean and process specific sample
python pipeline.py --config config.yml --target-sample "sb_rna_fracture_S3" --clean --steps fracture
```

## Advanced Usage

### Working with Test Data

1. Generate test data for specific sample:
```bash
python pipeline.py --config config.yml --target-sample "s4_dna_fracture_S5" --make-test --clean-test
```

2. Run pipeline on test data:
```bash
python pipeline.py --config config.yml --target-sample "TEST_s4_dna_fracture_S5" --steps parquet preprocess fracture
```

### Monitoring Progress

The pipeline includes progress bars for tracking:
- Overall step progress
- Individual operations within steps
- Assembly progress for each UMI

Example output:
```
Assembling contigs [######################] 100%
├── Validating input
├── Loading and preprocessing data
└── Running assembly
    └── Assembling molecules by UMI [#####################] 1000/1000
```

### Using Different Log Levels

Available log levels:
- `CRITICAL`: Only critical messages
- `INFO`: General information
- `IO`: Input/output operations
- `STEP`: Step-level progress
- `DEBUG`: Detailed debugging information

```bash
python pipeline.py --config config.yml --log-level DEBUG
```

## Troubleshooting

### Common Issues

1. Missing Configuration Parameters
```
Error: Missing required parameter: umi_len
```
Solution: Check your configuration file includes all required parameters for the step.

2. Invalid Sample Specification
```
Error: Target sample 'invalid_sample' not found in samples list
```
Solution: Verify the sample ID exists in the config file's samples list.

3. File Not Found
```
Error: No files found in /path/to/data
```
Solution: Verify file paths in configuration and check file permissions.

4. Memory Issues
```
Error: Memory error during assembly
```
Solution: Reduce batch size or use test data to verify pipeline.

### Debug Mode

Run in debug mode for detailed logging:
```bash
python pipeline.py --config config.yml --log-level DEBUG
```

### Cleaning Up

If you encounter issues, try cleaning previous outputs:
```bash
python pipeline.py --config config.yml --clean
```

### Getting Help

Run with --help to see all available options:
```bash
python pipeline.py --help
```

## Best Practices

1. Always use version control for configuration files
2. Start with test data before running on full dataset
3. Use dry runs to verify configuration
4. Monitor log files for errors and warnings
5. Keep backup copies of important data
6. Use descriptive sample names and consistent naming conventions
7. Verify sample IDs before processing
8. Clean outputs between different sample runs

## Next Steps

After successfully running the pipeline:
1. Check output quality metrics
2. Verify assembly success rates
3. Review generated plots if enabled
4. Archive results and configurations
5. Document any custom parameters used
6. Compare results across different samples

For more details, refer to the codebase documentation or reach out to the development team.
