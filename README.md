# Olivares Genomic toolkit

# Installation
```
git clone git@github.com:tzeitim/ogtk.git
cd ogtk
pip install -e .[sc,doc,rs]
```

Three main modules are supported.

# ltr
Comprehends lineage tracing analysis for bulk and single-cell data a various modalities. It is capable of pre-processing of multi-site lineage reporters such as GESTALT in addition to single-sites, namely shRNA.

# UM
This represents the main workframe for treating with single-molecule data. 

# Utils
Various recipes for anlayising single-cell RNA-seq data
([metacells](https://github.com/tanaylab/metacells) and
[single-cells](https://github.com/scverse/scanpy)) and implements wrappers for
other genomic toolkits such as [bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/).

---

Migration to [pola-rs](https://github.com/pola-rs) is work in progress but largely supported.


Notes:
If there are problems with missing zlib-devel not found install `pysam` via conda/mamba since zlib-devel cannot be installed using pip.

