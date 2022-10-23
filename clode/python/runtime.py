import clode.cpp.clode_cpp_wrapper as _clode
import os

_runtime = None
_clode_root_dir: str = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                                    "clode", "cpp", "")


def _get_runtime():
    global _runtime
    if _runtime is None:
        _runtime = _clode.opencl_resource()
    return _runtime
