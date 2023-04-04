# Olivares Genomic toolkit

Or something like that

Install `pysam` via conda since zlib-devel cannot be installed using pip.
# Installation
```
python3 setup.py sdist bdist_wheel
pip install setuptools_rust
pip install -e /local/users/polivar/src/projects/ogtk/
```

Three main modules are supported.

# ltr
Comprehends lineage tracing analysis for bulk and single-cell data a various modalities. It is capable of pre-processing of multi-site lineage reporters such as GESTALT in addition to single-sites, namely shRNA.

# UM
This represents the main workframe for treating with single-molecule data. 

# Utils
Various recipes for anlayising single-cell RNA-seq data (metacells and single-cell) and implements wrappers for other genomic toolkits such as bbtools.
