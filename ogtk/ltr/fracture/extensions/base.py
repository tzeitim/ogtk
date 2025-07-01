# ogtk/ltr/fracture/extensions/base.py
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Set, Type
from ..pipeline.types import StepResults, FractureXp
from .config import ExtensionConfig

class PostProcessorExtension(ABC):
    """Base class for post-processing extensions"""

    def __init__(self, xp: FractureXp):
        self.xp = xp
        self.config = self._load_config()
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Extension identifier"""
        pass
    
    @abstractmethod
    def process(self, contigs_path: Path) -> StepResults:
        """Process contigs and return results"""
        pass
    
    def should_run_step(self, step_name: str) -> bool:
        """Check if a specific extension step should run"""
        extension_steps_config = getattr(self.xp, 'extension_steps', {})
        extension_steps = extension_steps_config.get(self.name, None)
        
        if extension_steps is None:
            # If no specific steps configured, run all steps
            return True
        
        return step_name.lower() in [s.lower() for s in extension_steps]

    def _load_config(self):
        """Load and validate extension config"""
        config_class = self.get_config_class()
        if config_class:
            raw_config = self.xp.extension_config.get(self.name, {})
            return config_class.from_dict(raw_config)
        return {}
    
    @abstractmethod
    def get_config_class(self) -> Optional[Type[ExtensionConfig]]:
        """Return config class for this extension"""
        pass
    
