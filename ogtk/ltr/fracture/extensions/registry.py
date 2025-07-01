from typing import Dict, Optional, List, Type
from .base import PostProcessorExtension
from ..pipeline.types import FractureXp

class ExtensionRegistry:
    """Manages available post-processing extensions"""
    
    def __init__(self):
        self._extension_classes: Dict[str, Type[PostProcessorExtension]] = {}
    
    def register(self, extension_class: Type[PostProcessorExtension]):
        """Register an extension class"""
        # Get name from class property (need to create temp instance)
        temp_instance = extension_class.__new__(extension_class)
        name = temp_instance.name
        self._extension_classes[name] = extension_class
    
    def get_available(self) -> List[str]:
        """Get list of available extension names"""
        return list(self._extension_classes.keys())
    
    def create_extension(self, name: str, xp: FractureXp) -> Optional[PostProcessorExtension]:
        """Create extension instance with xp"""
        extension_class = self._extension_classes.get(name)
        if extension_class:
            return extension_class(xp)
        return None
    
    def __contains__(self, name: str) -> bool:
        """Check if extension is registered"""
        return name in self._extension_classes

# Global registry instance
extension_registry = ExtensionRegistry()
