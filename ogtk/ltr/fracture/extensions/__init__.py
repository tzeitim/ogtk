"""FRACTURE Pipeline Extension System"""
from .base import PostProcessorExtension
from .registry import ExtensionRegistry, extension_registry
from . import cassiopeia_lineage  # Auto-register extensions

__all__ = ['PostProcessorExtension', 'ExtensionRegistry', 'extension_registry']
