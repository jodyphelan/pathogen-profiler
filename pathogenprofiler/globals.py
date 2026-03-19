"""Shared global storage for the package"""
from typing import Any


class GlobalStorage:
    """A simple object for storing global key-value pairs accessible throughout the package."""
    
    def __getattr__(self, name: str) -> Any:
        try:
            return self.__dict__[name]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'") from None
    
    def __setattr__(self, name: str, value: Any) -> None:
        self.__dict__[name] = value
    
    def __delattr__(self, name: str) -> None:
        try:
            del self.__dict__[name]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'") from None
    
    def __contains__(self, name: str) -> bool:
        return name in self.__dict__
    
    def __iter__(self):
        return iter(self.__dict__)
    
    def get(self, name: str, default: Any = None) -> Any:
        return self.__dict__.get(name, default)
    
    def pop(self, name: str, default: Any = None) -> Any:
        return self.__dict__.pop(name, default)
    
    def clear(self) -> None:
        self.__dict__.clear()


# Create a single global instance accessible throughout the package
g = GlobalStorage()