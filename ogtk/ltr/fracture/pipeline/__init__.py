"""
FRACTURE Pipeline
=========================

FRagmentation and Assembly for Comprehensive Target Unit REconstruction.

A framework for processing and assembling molecule-based sequencing data. Provides 
functionality for converting FASTQ files to parquet format, preprocessing data,
and assembling contigs.

Components
---------

Core (:mod:`.core`)
    Contains pipeline infrastructure and execution logic
    
    - :class:`.Pipeline`: Main pipeline class for orchestrating data processing
    - :class:`.PipelineStep`: Pipeline step definitions

Plotting (:mod:`.plotting`)
    QC and results visualization
    
    - :class:`.PlotDB`: Plot generation functions

Plot Generation (:mod:`.plot_gen`)
    Plot regeneration tools
    
    - :class:`.PlotRegenerator`: Regenerates plots from previous runs

Types (:mod:`.types`)
    Type definitions
    
    - :class:`.FractureXp`: Extended experiment configuration

API Extensions (:mod:`.api_ext`) 
    Polars extensions
    
    - :class:`.PlDNA`: DNA sequence operations
    - :class:`.PlPipeline`: Assembly optimization
"""

from .core import Pipeline, PipelineStep
from .plotting import PlotDB
from .plot_gen import PlotRegenerator
from .types import FractureXp
from .api_ext import PlDNA, PlPipeline#, PllPipeline

__all__ = [
    'Pipeline',
    'PipelineStep',
    'PlotDB',
    'PlotRegenerator',
    'FractureXp',
    'PlDNA', 
    'PlPipeline',
#not implemented yet    'PllPipeline'
]
