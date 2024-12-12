# Olivares Genomic toolkit
# Installation

## Basic Installation

### Linux/Windows
```bash
pip install ogtk
```

### macOS
macOS requires additional system dependencies before installation:

1. Install LLVM and OpenMP support:
```bash
brew install llvm
```

2. Set up the environment (add this to your ~/.zshrc or ~/.bash_profile):
```bash
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
```

3. Reload your shell configuration:
```bash
source ~/.zshrc  # or source ~/.bash_profile
```

4. Install the package:
```bash
pip install ogtk
#or for all-inclusive functionality
pip install 'ogtk[sc,rs'] 
```

## With Rust Support

### Linux
```bash
pip install ogtk[rs]
```

### macOS
1. Install system dependencies:
```bash
brew install llvm patchelf
```

2. Set up the environment as described above

3. Install the package:
```bash
pip install ogtk[rs]
```
## astral uv

After having installed [astral-uv](https://docs.astral.sh/uv/)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Initialize project and add packages. 
```bash
uv init uvogtk
cd uvogtk
add ogtk
add ipython
uv venv --python 3.10
uv sync
```

## Troubleshooting

If you encounter build errors related to OpenMP on macOS:
1. Make sure LLVM is properly installed: `brew install llvm`
2. Verify environment variables are set correctly
3. Try reinstalling with: `pip install --no-cache-dir ogtk`

For other platforms, if you see OpenMP-related errors:
- Ubuntu/Debian: `sudo apt-get install libomp-dev`
- Fedora: `sudo dnf install libomp-devel`

# Functionality
Three main modules are supported.

## ltr
Comprehends lineage tracing analysis for bulk and single-cell data a various modalities. It is capable of pre-processing of multi-site lineage reporters such as GESTALT in addition to single-sites, namely shRNA.

## UM
This represents the main workframe for treating with single-molecule data. 

## Utils
Various recipes for anlayising single-cell RNA-seq data
([metacells](https://github.com/tanaylab/metacells) and
[single-cells](https://github.com/scverse/scanpy)) and implements wrappers for
other genomic toolkits such as [bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/).

---

Migration to [pola-rs](https://github.com/pola-rs) is work in progress but largely supported.


Notes:
If there are problems with missing zlib-devel not found install `pysam` via conda/mamba since zlib-devel cannot be installed using pip.

