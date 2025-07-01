# ogtk/ltr/fracture/extensions/base.py
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any, Optional, Set
from ..pipeline.types import StepResults, FractureXp

class PostProcessorExtension(ABC):
    """Base class for post-processing extensions"""

    def __init__(self, xp: FractureXp):
        self.xp = xp
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Extension identifier"""
        pass
    
    @property 
    @abstractmethod
    def required_params(self) -> Set[str]:
        """Parameters required from experiment config"""
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
