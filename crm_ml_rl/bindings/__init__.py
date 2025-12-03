"""
C++ bindings for MRI-Actuated Catheter Cosserat Rod Solver.

Wraps your existing C++ implementation (CRM + CRMDYN) for use in Python ML/RL.
"""

try:
    import crm_cpp
    __all__ = ['crm_cpp']
except ImportError as e:
    raise ImportError(
        "C++ bindings not found. Please build the extension module:\n"
        "  cd crm_ml_rl/bindings\n"
        "  mkdir build && cd build\n"
        "  cmake ..\n"
        "  make -j8\n"
    ) from e
