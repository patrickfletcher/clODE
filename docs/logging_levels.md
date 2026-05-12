# Logging and diagnostics

clODE uses the standard Python `logging` package. Runtime logs live under the `clode` logger namespace.

The package-level convenience helpers are:

- `clode.configure_logging(...)`: install a stream handler for `clode` logs and optionally capture Python warnings.
- `clode.get_logger(name)`: retrieve `clode` or one of its child loggers.

## Basic setup

```python
import clode

clode.configure_logging(level="INFO")
logger = clode.get_logger("demo")
logger.info("Starting simulator build")
```

You can pass either an integer logging level or a standard logging level name such as `"DEBUG"`, `"INFO"`, `"WARNING"`, or `"ERROR"`.

## PyOpenCL diagnostics

`configure_logging(...)` can also enable the PyOpenCL compiler-output switch for subsequent builds:

```python
import clode

clode.configure_logging(level="DEBUG", compiler_output=True)
```

That sets `PYOPENCL_COMPILER_OUTPUT=1`, which tells PyOpenCL to surface compiler messages emitted during `Program.build()`.

If you keep `capture_warnings=True` (the default), Python warnings such as `pyopencl.CompilerWarning` are routed through logging as well.

clODE also logs whether the selected device reports PyOpenCL source-build cache support when debug logging is enabled.

PyOpenCL already owns the main build-diagnostics switches, so prefer its native controls when debugging compilation:

- `PYOPENCL_COMPILER_OUTPUT=1`: show compiler messages during builds
- `PYOPENCL_NO_CACHE=1`: disable PyOpenCL's on-disk build cache
- `PYOPENCL_BUILD_OPTIONS=...`: append extra OpenCL build options

## Explicit reports

The following are explicit report helpers, not logger-driven output:

- `clode.print_opencl()`
- `simulator.print_devices()`
- `simulator.print_status()`

They print when you call them, regardless of the configured log level.

## Advanced use

If your application already configures logging globally, you can skip `clode.configure_logging(...)` entirely and just use standard logger configuration with the `clode` namespace:

```python
import logging

logging.basicConfig(level=logging.INFO)
logging.getLogger("clode").setLevel(logging.DEBUG)
```
