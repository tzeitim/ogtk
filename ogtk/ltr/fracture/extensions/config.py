from abc import ABC
from dataclasses import dataclass, fields, MISSING
from typing import Dict, Any, Set
import inspect

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
    
   
    def get_function_config(self, func_or_name) -> Dict[str, Any]:
        """Extract config automatically based on function signature"""
        
        # Get the function object
        if isinstance(func_or_name, str):
            # Try to get function from calling module's globals
            import sys
            frame = sys._getframe(1)
            caller_globals = frame.f_globals
            func = caller_globals.get(func_or_name)
            if not func:
                raise ValueError(f"Function {func_or_name} not found in calling module")
        else:
            func = func_or_name
            
        # Get function signature
        sig = inspect.signature(func)
        config_dict = self.to_dict()
        
        # Extract parameters, respecting defaults
        func_params = {}
        data_params = {'ldf', 'df', 'refs', 'self'}  # Skip these common params
        
        for param_name, param in sig.parameters.items():
            # Skip data/context parameters
            if param_name in data_params:
                continue
                
            # Check if we have this parameter in config
            if param_name in config_dict:
                config_value = config_dict[param_name]
                
                # Only include if config value is not None or if param has no default
                if config_value is not None or param.default == inspect.Parameter.empty:
                    func_params[param_name] = config_value
                    
            elif param.default == inspect.Parameter.empty:
                # Required parameter missing from config
                raise ValueError(
                    f"Required parameter '{param_name}' missing from config for function {func.__name__}"
                )
            # If param has default and not in config, let function default handle it
                
        return func_params
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary"""
        if hasattr(self, '__dataclass_fields__'):
            return {f.name: getattr(self, f.name) for f in fields(self)}
        return vars(self)
