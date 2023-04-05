import os

from . import clode_cpp_wrapper as _clode

_runtime = None
_clode_root_dir: str = os.path.join(os.path.dirname(__file__), "cpp", "")

def _get_runtime():
    global _runtime
    if _runtime is None:
        _runtime = _clode.opencl_resource()
    return _runtime

def _get_clode():
    return _clode

def print_devices():
    return _get_runtime().print_open_cl()