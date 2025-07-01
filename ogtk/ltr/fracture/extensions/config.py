from abc import ABC
from dataclasses import dataclass, fields, MISSING
from typing import Dict, Any, Set

class ExtensionConfig(ABC):
    """Base class for extension configurations"""
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'ExtensionConfig':
        """Create config from dictionary with validation"""
        # Get class field names
        field_names = {f.name for f in fields(cls)} if hasattr(cls, '__dataclass_fields__') else set()
        
        # Filter config_dict to only include valid fields
        filtered_config = {k: v for k, v in config_dict.items() if k in field_names}
        
        # Check for missing required fields
        missing = cls.get_required_fields() - set(filtered_config.keys())
        if missing:
            raise ValueError(f"{cls.__name__} missing required config: {missing}")
            
        return cls(**filtered_config)
    
    @classmethod 
    def get_required_fields(cls) -> Set[str]:
        """Get fields without defaults"""
        if not hasattr(cls, '__dataclass_fields__'):
            return set()
            
        required = set()
        for field in fields(cls):
            if field.default == field.default_factory == MISSING:
                required.add(field.name)
        return required
    
    @classmethod
    def get_function_params(cls, func_name: str) -> Set[str]:
        """Get parameters needed for a specific function"""
        if not hasattr(cls, '_FUNCTION_PARAMS'):
            return set()
        return cls._FUNCTION_PARAMS.get(func_name, set())
    
    def get_function_config(self, func_name: str) -> Dict[str, Any]:
        """Get config subset for specific function"""
        needed_params = self.get_function_params(func_name)
        if not needed_params:
            # If no specific params defined, return all config
            return self.to_dict()
            
        return {k: getattr(self, k) for k in needed_params if hasattr(self, k)}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary"""
        if hasattr(self, '__dataclass_fields__'):
            return {f.name: getattr(self, f.name) for f in fields(self)}
        return vars(self)
