# IVP API And Test Plan

Purpose: define the first-pass `InitialValueProblem` API, the delivery order for the active PR, and a deliberately slim test matrix.
Read when: implementing the IVP PR, reviewing the public problem-facing API, or deciding what to test first.
Update when: the staged plan, exposed IVP surface, or agreed test scope changes.

## Bottom line

- Treat `InitialValueProblem` as the user-facing semantic owner of RHS semantics, default state, default parameters, and basic batch shaping.
- Keep the strict OpenCL `getRHS(...)` contract internal to compilation and executor plumbing.
- Preserve simulator classes as orchestration handles that build or consume an IVP.
- Keep the curated `clode.problem` public story centered on `InitialValueProblem`, `OpenCLConverter`, and format-ingestion helpers; treat `ProblemInfo`, `RhsSource`, and low-level source preparation as derived support concepts rather than promoted user-facing API.
- Add SciPy only as a development and test dependency for this PR.
- Keep the first test pass small: focus on ownership, shaping contracts, Python-backed callability, and one or two end-to-end confirmation paths.

## Delivery order for this PR

1. Settle the IVP API and the accompanying test matrix before moving code.
2. Land IVP ownership and simulator delegation.
3. Add Python-backed IVP callability for SciPy-style `(t, y) -> dydt` use.

This sequencing is intentional. It avoids mixing semantic ownership changes with callability mechanics before the object boundary is clear.

## First-pass IVP API direction

### Core construction

The first-pass IVP should be constructible from the same problem-definition sources already supported today:

- Python RHS
- OpenCL source file
- XPP file after conversion

The object should own:

- ordered variable names
- ordered parameter names
- auxiliary variable names
- `num_noise`
- default initial-state values
- default parameter values
- RHS source metadata needed by the OpenCL build path
- remembered batch shape and ensemble size

The object should not own:

- solver parameters
- runtime selection
- observer policy
- feature or trajectory storage policy
- continuation state

### Proposed public shape

This is a planning target, not frozen API text:

```python
ivp = clode.InitialValueProblem(
    variables={"x": 1.0, "y": 0.0},
    parameters={"mu": 1.0},
    rhs_equation=get_rhs,
    aux=["energy"],
    num_noise=0,
)
```

Expected first-pass properties or read-only accessors:

- `variable_names`
- `parameter_names`
- `aux_names`
- `num_variables`
- `num_parameters`
- `num_aux`
- `num_noise`
- `default_initial_state`
- `default_parameters`
- `ensemble_shape`
- `ensemble_size`

Expected first-pass mutating or builder-style operations:

- `with_defaults(...)` or equivalent update path for default state and parameters
- `set_ensemble(...)` or equivalent batch-shaping method using the current simulator semantics as the starting point
- `set_repeat_ensemble(...)` or equivalent repeat helper
- a normalization helper that returns the `(ensemble_size, num_variables)` and `(ensemble_size, num_parameters)` arrays the executors need

Current implementation note:

- the landed branch currently exposes `get_initial_state()`, `get_parameter_values()`, `set_ensemble()`, `set_repeat_ensemble()`, `evaluate_rhs(...)`, and `__call__(...)` as the main user-facing IVP operations
- direct source-preparation and raw problem-data mutation helpers are now internal implementation details

The exact mutator names are less important than the ownership boundary. Batch shaping should move under IVP ownership even if some compatibility methods remain on `Simulator` as thin delegates for one release cycle.

### Simulator integration target

The simulator layer should move toward one of these forms:

```python
sim = clode.FeatureSimulator(ivp=ivp, observer=..., stepper=..., ...)
```

or a compatibility constructor that internally builds the IVP:

```python
sim = clode.FeatureSimulator(
    variables={...},
    parameters={...},
    rhs_equation=get_rhs,
    observer=..., stepper=..., ...
)
```

For this PR, the compatibility path is desirable to avoid a broad public break while the API is still settling.

## Python-backed callability

### Goal

Let a Python-authored IVP feel natural in a SciPy-style workflow without exposing the OpenCL boundary contract as the primary user API.

### Recommended approach

If an IVP is built from a Python RHS, preserve that Python callable and expose:

- an explicit adapter method such as `evaluate_rhs(t, y, *args, parameters=None, aux=None, wiener=None)`
- `__call__(t, y, *args)` delegating to the default SciPy-style evaluation path

The default `__call__` behavior should:

- accept a one-dimensional state vector `y`
- use the IVP's default parameter values when no overrides are supplied
- support `solve_ivp(..., args=(...))` by mapping positional `*args` onto parameter overrides in declared parameter order
- allocate temporary derivative, aux, and wiener arrays as needed
- return the derivative vector as a NumPy array

This keeps the strict converter signature internal while preserving correctness for the existing converter path.

### Non-goals for this pass

- no generic callability for OpenCL-only or XPP-only IVPs unless there is still a retained Python callable
- no attempt to reconstruct Python semantics from generated OpenCL text
- no bidirectional or round-trippable RHS IR project in this PR
- no vectorized SciPy adapter beyond the ordinary `(t, y)` contract unless it falls out cheaply later

## Test matrix

Bias toward a small, high-value matrix. The goal is to protect the new object boundary and one useful interop story, not to freeze every evolving convenience API.

### Unit and contract tests to add

1. IVP construction from a Python RHS preserves ordered names, defaults, and metadata.
2. IVP batch shaping reproduces the current size, shape, scalar-broadcast, and mixed-array rules now implemented on `Simulator`.
3. Simulator compatibility construction builds or consumes an IVP without changing core solve behavior.
4. Python-backed `ivp(t, y)` returns the same derivative as the wrapped RHS for the default parameter set.
5. Python-backed `ivp(t, y, *args)` matches the declared parameter ordering so SciPy `solve_ivp(..., args=(...))` can override parameters naturally.
6. Python-backed `ivp.evaluate_rhs(...)` supports explicit parameter overrides if that API lands in the same pass.
7. Non-Python-backed IVPs fail clearly on `__call__` with a precise error message.

### End-to-end checks to keep slim

1. One runtime/API test that constructing a simulator from an IVP reaches the current OpenCL solve path and returns expected shape-level results.
2. One SciPy-backed comparison test for a simple deterministic model, used only to validate the Python-backed callability path and basic agreement of problem specification.

### Tests explicitly not required in this PR

- broad duplication of existing numerical regression coverage
- an exhaustive matrix over every simulator subclass, observer, and stepper
- extensive API snapshot testing while the IVP surface is still settling
- tests for generic SciPy interop across OpenCL-only or XPP-only problems

## Dependency plan

- add `scipy` to development and test extras, not to required runtime dependencies
- keep SciPy-dependent tests isolated so the rest of the frontend or runtime suite does not depend on SciPy semantics

If marker separation becomes useful, a small frontend or runtime marker split for SciPy-backed IVP tests is acceptable, but avoid proliferating markers unless the suite actually needs it.

## Likely file targets

- `clode/problem/definition.py` or a neighboring problem-layer module for `InitialValueProblem`
- `clode/problem/__init__.py` and `clode/__init__.py` for exports
- `clode/simulation/base.py` and simulator subclasses for IVP ownership and delegation
- `test/` contract tests focused on IVP ownership and Python-backed callability
- `pyproject.toml` for SciPy as an optional dev/test dependency
- `docs/specifying_odes.md` once the wrapper TODO becomes real behavior
