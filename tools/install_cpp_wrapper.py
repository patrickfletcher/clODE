from __future__ import annotations

import shutil
import sys
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CPP_DIR = REPO_ROOT / "clode" / "cpp"
WRAPPER_CANDIDATES = {
    "clode_cpp_wrapper",
    "clode_cpp_wrapper.so",
    "libclode_cpp_wrapper.so",
    "libclode_cpp_wrapper.dylib",
    "clode_cpp_wrapper.dll",
    "clode_cpp_wrapper.pyd",
}


def _extension_suffix() -> str:
    if sys.platform.startswith("win"):
        for suffix in EXTENSION_SUFFIXES:
            if suffix.endswith(".pyd"):
                return suffix
    for suffix in EXTENSION_SUFFIXES:
        if suffix.endswith(".so"):
            return suffix
    raise SystemExit("Could not determine a Python extension-module suffix for this interpreter.")


def _find_wrapper_binary() -> Path:
    candidates: list[Path] = []
    for root_name in ("bazel-bin", "bazel-out"):
        root = REPO_ROOT / root_name
        if not root.exists():
            continue
        candidates.extend(
            path
            for path in root.rglob("*")
            if path.is_file() and path.name in WRAPPER_CANDIDATES
        )

    if not candidates:
        raise SystemExit(
            "No Bazel-built clode_cpp_wrapper shared library was found. "
            "Run 'bazel build //clode/cpp:clode_cpp_wrapper' first."
        )

    candidates.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    return candidates[0]


def main() -> int:
    source = _find_wrapper_binary()
    destination = CPP_DIR / f"clode_cpp_wrapper{_extension_suffix()}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    print(destination.relative_to(REPO_ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
