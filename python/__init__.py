"""CRM Catheter - Cosserat Rod Model for MRI-Actuated Robotic Catheter."""

__version__ = "1.0.0"
__author__ = "Your Name"

try:
    from . import crm_cpp
    _cpp_available = True
except ImportError:
    _cpp_available = False
    import warnings
    warnings.warn("C++ bindings not available. Build with CMake first.")

if _cpp_available:
    from .crm_wrapper import CatheterModel

from . import envs
from . import models
from . import rl

__all__ = [
    'CatheterModel',
    'envs',
    'models',
    'rl',
    'crm_cpp',
]
