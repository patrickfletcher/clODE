# Python Simulation Flow Reference

Status: legacy comparison reference. This note documents the wrapper-backed simulator flow that remains available only through explicit `_CLODE_BACKEND=cpp` selection in a source checkout.

## Purpose

This note documents how the current Python-side simulator wrappers actually drive the C++ and OpenCL layers, with emphasis on repeated calls and continuation.

The focus is on:

- constructor defaults
- `transient()`
- `trajectory()`
- `features()`
- how solver state is advanced on the device
- which helper calls mutate device state versus only mutate Python caches
- what repeated calls on the same simulator object actually mean today

This is based on direct code inspection plus small verification runs on platform 1 device 0.

## The Four State Layers

The current system has four relevant layers of state.

### 1. Python wrapper state

The Python classes keep local bookkeeping such as:

- `_sp`: Python-side solver parameter struct
- `_t_span`: Python-side copy of requested time span
- `_ensemble_size`, `_ensemble_shape`
- `_device_initial_state`, `_device_parameters`
- cached fetched results such as `_device_final_state`, `_device_dt`, `_device_tf`
- trajectory caches: `_device_t`, `_device_x`, `_device_dx`, `_device_aux`
- feature caches: `_device_features`
- feature observer parameters: `_op`

These fields are not the solver itself. They are Python-side mirrors or caches.

### 2. C++ host-side state

The C++ `CLODE` object keeps host vectors and metadata such as:

- `tspan`
- `x0`
- `pars`
- `xf`
- `dt`
- `tf`
- `RNGstate`
- `sp`

The trajectory subclass adds host output arrays:

- `t`, `x`, `dx`, `aux`, `nStored`

The features subclass adds:

- `F`
- `op`
- `observerInitialized`

### 3. OpenCL device buffers

The actual live solver state used by kernels sits in device buffers:

- `d_tspan`
- `d_sp`
- `d_x0`
- `d_pars`
- `d_xf`
- `d_dt`
- `d_tf`
- `d_RNGstate`

Trajectory adds device output buffers:

- `d_t`, `d_x`, `d_dx`, `d_aux`, `d_nStored`

Features adds:

- `d_odata`
- `d_op`
- `d_F`

### 4. Kernel-local state

Each kernel invocation copies device values into private locals like:

- `ti`
- `dt`
- `xi[]`
- `p[]`
- `dxi[]`
- `auxi[]`
- RNG state `rd`

The kernels integrate using those locals, then write selected results back to device buffers.

## Constructor Flow

## Base `Simulator` constructor

The Python constructor in `clode/solver.py` does the following in order.

1. Resolve the model source.
   - `src_file` is used directly unless it is an `.xpp` file, in which case it is converted.
   - `rhs_equation` is converted into a generated `clode_rhs.cl` file in the current working directory.

2. Build `ProblemInfo` from the provided variable names, parameter names, aux names, and `num_noise`.

3. Choose the OpenCL runtime with `initialize_runtime(...)`.
   - If no explicit platform/device is given, the wrapper asks for `DEVICE_TYPE_DEFAULT` and `VENDOR_ANY`.

4. Create the pybind-backed integrator object.
   - Base class: `SimulatorBase`
   - Trajectory subclass: `TrajectorySimulatorBase`
   - Feature subclass: `FeatureSimulatorBase`

5. Build the OpenCL program immediately.
   - The Python wrappers do not delay initial compilation.
   - The base simulator compiles `transient.cl` plus the model source.
   - The trajectory simulator compiles `transient.cl`, `trajectory.cl`, plus the model source.
   - The feature simulator compiles `transient.cl`, `initializeObserver.cl`, `features.cl`, plus the model source and observer build defines.

6. Set solver parameters.
   - If no explicit `SolverParams` is passed, Python constructs one from the constructor arguments.
   - Python calls `set_solver_parameters()`, which pushes the `SolverParams` struct to `d_sp`.

7. Set the time span.
   - Python stores `_t_span` and pushes the same values to `d_tspan`.

8. Set default problem data.
   - Python stores the default variable and parameter dictionaries.
   - It constructs a size-1 initial state array and size-1 parameter array.
   - It calls `_set_problem_data(...)`, which pushes both to the C++ layer.

9. C++ `setProblemData(...)` calls `setNpts(1)`.
   - This allocates all `nPts`-dependent device buffers.
   - It allocates `d_x0`, `d_pars`, `d_RNGstate`, `d_dt`, `d_tf`, `d_xf`.
   - It initializes `d_dt` from `sp.dt`.
   - It automatically calls `seedRNG()` after resizing.

10. C++ writes the initial state and parameters into `d_x0` and `d_pars`.

After this, the object is ready to simulate without any further setup.

## Subclass-specific constructor details

### `TrajectorySimulator`

The only extra constructor logic on the Python side is cache initialization:

- `_device_t = None`
- `_device_x = None`
- `_device_dx = None`
- `_device_aux = None`

The actual kernel specialization is already chosen by overriding `_create_integrator()` before the base constructor runs.

### `FeatureSimulator`

The feature simulator creates observer parameters before calling the base class constructor.

That matters because the C++ feature integrator constructor needs:

- the observer name
- the `ObserverParams`

and the initial OpenCL build needs those to determine feature layout and `ObserverData` size.

So observer choice and observer parameter defaults are effectively part of the initial compiled program configuration.

## The Atomic State Mutations

These methods are the building blocks that the higher-level simulation calls compose.

## `set_tspan(t_span)`

Python behavior:

- store `_t_span`
- call `self._integrator.set_tspan(t_span)`

C++ behavior:

- replace host `tspan`
- allocate a fresh `d_tspan`
- copy the two endpoint values to the device

This does not touch:

- `d_x0`
- `d_pars`
- `d_dt`
- `d_RNGstate`
- cached Python final results

So `set_tspan()` changes only the requested time window.

## `shift_tspan()`

Python behavior:

- call C++ `shift_tspan()`
- refresh `_t_span` from the integrator

C++ behavior:

- compute `newTspan = {tspan[1], tspan[1] + (tspan[1] - tspan[0])}`
- push that back through `setTspan(...)`

This is important:

- `shift_tspan()` uses the previously requested span
- it does not use the actual reported `tf`

So it advances the requested window by its nominal duration, not by the actual endpoint reached on the device.

### When `shift_tspan()` matters most

This distinction matters differently depending on the model class.

- For autonomous systems, the RHS depends only on state and parameters, so restarting the next kernel with the same nominal time window does not by itself change the continuous dynamics seen by the solver.
- For non-autonomous systems, the RHS can depend explicitly on `t`, so reusing the same nominal time window repeats the same forcing/time-dependent coefficients instead of continuing from the true attained time.
- It also matters when an adaptive solver does not actually end at the requested final time, because the nominal next window can drift away from the true elapsed time seen by the solver.

So `shift_tspan()` being requested-time based is mostly harmless for autonomous state evolution, but it is potentially important for non-autonomous problems and for any workflow that needs globally meaningful timestamps.

## `shift_x0()`

There is no Python wrapper method on `Simulator`, but `transient()`, `trajectory()`, and `features()` all call the C++ `shift_x0()` internally when `update_x0=True`.

C++ behavior:

- device-to-device copy: `d_x0 <- d_xf`

This is the central state-continuation primitive in the current design.

## `set_solver_parameters(...)`

Python behavior:

- update `_sp`
- push `_sp` to the integrator with `set_solver_params(...)`

C++ behavior:

- replace host `sp`
- fill host `dt` vector with `sp.dt`
- rewrite `d_sp`

Critical current detail:

- C++ `setSolverParams(...)` does not copy the refreshed host `dt` vector back into `d_dt`

That means changing `dt` through `set_solver_parameters(...)` does not reset the per-trajectory device timestep buffer on an existing ensemble.

In practice, the next kernel launch still starts from whatever `d_dt` currently holds, which is usually the final `dt` from the previous run.

This matters most for adaptive solvers.

Verified observation:

- after an adaptive transient, `get_dt()` returned about `9.54e-07`
- calling `set_solver_parameters(dt=0.2, dtmax=0.2)` left `get_dt()` unchanged
- a second run still started from that tiny persisted `d_dt`

## `seed_rng(seed=None)`

Python behavior:

- call C++ `seed_rng()` or `seed_rng(seed)`

C++ behavior:

- overwrite `RNGstate`
- copy it to `d_RNGstate`

This is the only explicit RNG reset API exposed at the Python level.

Also note:

- any `nPts` change calls `setNpts(...)`
- `setNpts(...)` automatically calls random `seedRNG()`

So changing ensemble size implicitly reseeds the stochastic state unless the user reseeds explicitly afterward.

## `set_repeat_ensemble()` and `set_ensemble()`

These are not just shape setters. They also define what counts as the next initial state.

The internal logic in `_make_problem_data(...)` is:

- if the previous ensemble size equals the new size, or the previous size is `1`, then use the current device-backed initial state and parameters as the base arrays
- otherwise, fall back to the original constructor defaults

For variables, the base array comes from `get_initial_state()`, which fetches `d_x0` if necessary.

That means:

- after a run with `update_x0=True`, the shifted final state becomes the base state for later ensemble construction
- `set_repeat_ensemble(3)` after a transient duplicates the current continued state, not the original constructor state

Verified observation:

- after one transient, `x0` became `[[1.21306121, -0.42975727]]`
- `set_repeat_ensemble(3)` replicated that continued state into all three ensemble members

Parameters follow the same persistence rule using `_device_parameters`.

## Getter behavior

### `get_initial_state()`

- lazily fetches `d_x0`
- reshapes to `(ensemble_size, num_variables)` in Fortran order

After any simulation call with `update_x0=True`, this returns the continued state because `d_x0` has already been replaced with `d_xf`.

### `get_final_state()`

- intended to lazily fetch `d_xf`
- currently broken on the second call without invalidation

Current bug:

- the local variable `final_state` is only assigned when `_device_final_state is None`
- on a second call, Python raises `UnboundLocalError`

Verified observation:

- first `get_final_state()` succeeds
- second `get_final_state()` on the same cached result raises `UnboundLocalError`

### `get_dt()`

- lazily fetches `d_dt`
- returns one value per trajectory in the ensemble

This is the actual device-side continuation timestep state.

### `get_final_time()`

- lazily fetches `d_tf`
- returns one value per trajectory in the ensemble

### `get_tspan()`

Current bug:

- Python refreshes `_t_span` from the integrator
- but never returns it

So `get_tspan()` always returns `None` today.

Verified observation:

- `sim.get_tspan()` returns `None`
- `sim.shift_tspan()` also returns `None`, which is fine for a mutator
- after shifting, the internal cached `_t_span` does update

## Default `transient()` Flow

The Python `transient()` method composes the atomic operations in this order.

1. If `t_span` argument is supplied, call `set_tspan(...)`.
2. Launch the C++ transient kernel.
3. Invalidate Python caches for:
   - `_device_final_state`
   - `_device_dt`
   - `_device_tf`
4. If `update_x0=True` (the default), call `shift_x0()` on the integrator.
5. If `update_x0=True`, invalidate `_device_initial_state`.
6. If `fetch_results=True`, fetch and return `get_final_state()`.

Kernel-side, the transient kernel:

- starts from `ti = tspan[0]`
- loads `dt` from `d_dt[i]`
- loads `x0` from `d_x0`
- loads RNG state from `d_RNGstate`
- writes back `xf`, `tf`, `d_dt`, and `d_RNGstate`
- does not write back `d_x0`

That last point is why `shift_x0()` is the continuation step.

## What repeated `transient()` means by default

With the defaults, a second call to `transient()` on the same object means:

- start from the previous `xf` because `d_x0` was shifted
- start at the same requested `tspan[0]` unless the caller changed `tspan`
- start from the previous device `d_dt`
- start from the previous device RNG state

So the default is:

- state continuation
- timestep continuation
- RNG continuation
- but not automatic time-window continuation

Verified observation:

- two consecutive default transients both reported `tf = 1.000000119...`
- the second run started from the first run's final state
- time restarted from the same requested window

This is the most important continuation fact on the base simulator.

## Default `trajectory()` Flow

The Python `trajectory()` method does this.

1. If `t_span` argument is supplied, call `set_tspan(...)`.
2. Call C++ `trajectory()`.
3. Invalidate Python trajectory caches:
   - `_device_t`
   - `_device_x`
   - `_device_dx`
   - `_device_aux`
4. Invalidate base final-result caches:
   - `_device_final_state`
   - `_device_dt`
   - `_device_tf`
5. If `update_x0=True` (default), call `shift_x0()`.
6. If `fetch_results=True` (default), fetch trajectory arrays and package them as `TrajectoryOutput`.

C++ `CLODEtrajectory::trajectory()` first calls `resizeTrajectoryVariables()`.

That means each call:

- ensures the output buffers exist for the current `nPts` and `sp.max_store`
- overwrites the trajectory output buffers for the current run
- does not append to previous outputs

The kernel stores:

- the initial point
- then every `sp.nout`-th accepted step
- plus `xf`, `tf`, `d_dt`, RNG state

## What repeated `trajectory()` means by default

It means:

- the next segment starts from the previous final state because `d_x0` was shifted
- the returned time array still starts at `tspan[0]`
- the returned trajectory is only the current segment
- the caller must manually concatenate segment outputs if a long continuation is desired

Verified observation:

- first trajectory segment ran from `t=0.0` to `t=1.000000119...`
- second trajectory segment also ran from `t=0.0` to `t=1.000000119...`
- the first stored state of the second segment equaled the final state of the first segment

So trajectory continuation is state continuation with fresh output buffers, not automatic time concatenation.

## Default `features()` Flow

The feature simulator has one extra state machine: persistent observer state.

### The observer state lives in `d_odata`

The key C++ flag is:

- `observerInitialized`

and the persistent device buffer is:

- `d_odata`

`CLODEfeatures::features()` does this:

1. `resizeFeaturesVariables()`
2. if `observerInitialized` is false, call `initializeObserver()`
3. run the feature kernel

`resizeFeaturesVariables()` sets `observerInitialized = false` only when feature buffers are reallocated, which mainly happens when `nPts` changes.

So once initialized, observer state persists across later feature runs on the same object.

### Python `features()` wrapper flow

The Python wrapper does this.

1. If `t_span` argument is supplied, call `set_tspan(...)`.
2. If `initialize_observer is None`, call the no-argument C++ `features()`.
3. If `initialize_observer` is a boolean, call the overloaded C++ `features(bool)`.
4. On the default no-argument path only, invalidate:
   - `_device_features`
   - `_device_final_state`
   - `_device_dt`
   - `_device_tf`
5. If `update_x0=True` (default), call `shift_x0()`.
6. If `fetch_results=True` (default), return `get_observer_results()`.

### How the C++ boolean overload works

`CLODEfeatures::features(bool reinitialize_observer)` does:

- `observerInitialized = !reinitialize_observer`
- then calls the normal `features()`

So:

- `features(initialize_observer=True)` forces reinitialization
- `features(initialize_observer=False)` forces continuation of existing observer state

### What `initializeObserver()` itself changes

The initialization kernel reads:

- `d_tspan`
- `d_x0`
- `d_pars`
- `d_sp`
- `d_dt`
- `d_RNGstate`

but it only writes back `d_odata`.

For one-pass observers like `basicall`, initialization is essentially a reset of observer accumulators.

For two-pass detectors, `initializeObserver()` may integrate over the current `tspan` to warm up detector thresholds, but it rewinds the local state before finalizing initialization and still only writes observer data.

So `initializeObserver()` does not advance `x0`, `tf`, `d_dt`, or RNG state on the device.

### Observer continuation uses the true attained interval

The feature kernel calls `finalizeObserverData(...)` after each run.

Across the observer implementations, this cleanup step typically computes:

```c
realtype T = *ti - tspan[0];
```

and then shifts time-based observer members by that true attained interval.

Examples:

- `basic` and `basicall` shift `t_start`
- `localmax` shifts `t_start`, `tLastMax`, `tLastMin`, and the rolling `tbuffer`
- threshold/neighbourhood observers use the same pattern for their continuation-relevant time members

This is an important correction to a simplistic reading of the Python wrapper behavior: even if Python does not advance `tspan`, the persistent observer state is not naively tied to the nominal requested interval. Its internal time base is rebased by the actual `tf - tspan[0]` attained by the kernel.

That makes observer continuation more faithful than the bare `shift_tspan()` helper.

### But this applies to continuation state, not every emitted timestamp feature

The rebasing in `finalizeObserverData(...)` is aimed at the internal observer state needed for the next run.

It does not imply that every timestamp-like feature written to `F` is turned into a globally continuous absolute time. In event-based observers, the emitted timestamp arrays are written during `finalizeFeatures(...)`, and some stored event timestamp lists are not themselves rebased afterward.

So there are two separate notions:

- internal observer time members used to continue detection/statistics across runs
- feature outputs that the user may interpret as timestamps

The first is explicitly rebased by true elapsed time. The second is not guaranteed to be globally continuous across repeated windows.

## What repeated `features()` means by default

By default, repeated `features()` calls on the same object mean:

- continue from the previous final state because `shift_x0()` runs by default
- continue from the previous device `d_dt`
- continue from the previous RNG state
- continue from the previous observer state unless reinitialized
- do not automatically shift the requested `tspan`

The last two points need to be read together:

- Python does not automatically advance the requested time window
- but the observer continuation state is rebased by the actual elapsed interval inside `finalizeObserverData(...)`

So for autonomous systems and relative-time observer logic, repeated `features()` calls are more internally consistent than the Python-level `tspan` handling alone would suggest.

This is the most continuation-heavy wrapper of the three.

So the intended kernel-side pattern for multi-window feature measurement is:

1. advance state with `update_x0=True`
2. advance the requested time window manually when absolute time matters
3. choose whether observer state should continue or reinitialize

## Continuation patterns supported by the current defaults

## Pattern 1: state-only continuation

This is what repeated calls do automatically.

If the caller simply does:

```python
sim.transient()
sim.transient()
```

or:

```python
traj.trajectory()
traj.trajectory()
```

then the second run starts from the previous final state, but with the same requested time window.

This is state continuation, not full time continuation.

For autonomous systems, that is often numerically acceptable as a continuation of the state trajectory, because absolute time does not enter the RHS.

For non-autonomous systems, it is not true absolute-time continuation because the next kernel again starts from the same nominal `tspan[0]`.

## Pattern 2: requested-window continuation

To advance both state and requested time window using the built-in helpers, the caller must do both:

- keep `update_x0=True`
- also call `shift_tspan()` or `set_tspan(...)`

Example intent:

```python
sim.transient()
sim.shift_tspan()
sim.transient()
```

Important current detail:

- `shift_tspan()` uses the nominal requested interval, not the actual reported `tf`

So this is requested-window continuation, not exact actual-time continuation.

That distinction is most important when:

- the RHS is explicitly time-dependent
- a solver does not finish exactly at the requested endpoint
- downstream code expects globally continuous timestamps

## Pattern 3: feature continuation

For features, repeated calls with defaults continue both:

- the solver state in `d_x0`
- the observer state in `d_odata`

That is the default path most aligned with cumulative feature measurement over successive windows.

And because `finalizeObserverData(...)` rebases the internal observer time members by the true attained interval, this continuation is closer to true cumulative measurement than the Python-level `tspan` handling alone would imply.

But the caller still has to manage the requested time window explicitly when absolute time itself matters, especially for non-autonomous systems or globally interpreted event timestamps.

## Pattern 4: restart from the same initial condition

To rerun from the same initial state on the same object, the caller must not rely on the defaults.

The defaults shift `x0`.

A restart-like run requires at least:

- `update_x0=False`
- and, for stochastic or adaptive solvers, explicit handling of RNG state and `d_dt`

Because even with `update_x0=False`, the kernels still mutate:

- `d_dt`
- `d_RNGstate`

So `update_x0=False` alone is not a full reset.

## Python-side cache invalidation behavior

The wrapper caches are only partially synchronized.

### Good invalidation cases

- `transient()` clears cached final state, `dt`, and `tf`
- `trajectory()` clears cached trajectory arrays plus final state, `dt`, and `tf`
- `features()` on the default no-argument path clears cached features plus final state, `dt`, and `tf`
- all three clear cached initial state when `update_x0=True`

### Weak or missing invalidation cases

- `set_tspan()` does not invalidate cached `tf`
- `set_solver_parameters()` does not invalidate cached `dt`
- `_set_problem_data()` does not invalidate cached final results
- `features(initialize_observer=...)` does not invalidate cached features or final results

So once getters have been called, later setter calls can leave Python-side caches stale unless a simulation call happens to clear them.

## Current quirks and bugs relevant to continuation

These are the main issues uncovered while tracing the default flow.

## 1. `get_tspan()` never returns the span

Current behavior:

- refreshes `_t_span`
- returns `None`

This is a Python wrapper bug.

## 2. `get_final_state()` breaks on a second cached call

Current behavior:

- first call works
- second call raises `UnboundLocalError`

This is a Python wrapper bug in cache handling.

## 3. `set_solver_parameters()` does not reset the device `d_dt` buffer

This is a C++ state propagation bug with strong continuation consequences.

It means:

- `sp.dt` changes
- host `dt` vector changes
- but kernel start `dt` does not reset on the device

So adaptive continuation is partly controlled by hidden prior device state.

## 4. `features(initialize_observer=...)` skips cache invalidation

Current behavior:

- prints `Setting initialize_observer=...` to stdout
- runs the kernel
- but does not clear `_device_features`, `_device_final_state`, `_device_dt`, or `_device_tf`

So a same-object feature rerun through the boolean path can return stale cached Python results.

## 5. `_cl_program_is_valid` is currently inert

The Python wrappers set `_cl_program_is_valid = False` in some places, such as:

- `FeatureSimulator.set_observer(...)`
- `FeatureSimulator.set_observer_parameters(...)` when `max_event_timestamps` changes

But the lazy rebuild blocks in `transient()`, `trajectory()`, and `features()` are commented out.

So marking the program invalid does not automatically rebuild it on the next run.

## 6. `shift_tspan()` is requested-time based, not actual-time based

This is not a wrapper bug by itself, but it matters for continuation.

Verified observation:

- after a run with actual `tf = 1.000000119...`
- `shift_tspan()` still produced cached `_t_span = [1.0, 2.0]`

So time continuation is tied to the requested window, not to the reported final time.

For autonomous systems, this mainly affects reported/stored time coordinates rather than the state evolution itself.

For non-autonomous systems, this can change the actual dynamics because future kernel calls see the wrong absolute time.

## 7. Ensemble resizing changes RNG behavior

If `set_ensemble()` or `set_repeat_ensemble()` changes `nPts`, then `setNpts(...)` runs and auto-calls random `seedRNG()`.

So stochastic continuation across ensemble-size changes is not seed-preserving by default.

## Summary of what repeated calls mean today

## `Simulator.transient()`

Default repeated call behavior:

- continue `x0`
- continue `d_dt`
- continue RNG state
- keep the same requested `tspan` unless changed manually

## `TrajectorySimulator.trajectory()`

Default repeated call behavior:

- continue `x0`
- continue `d_dt`
- continue RNG state
- overwrite output buffers for the new segment
- keep the same requested `tspan` unless changed manually

## `FeatureSimulator.features()`

Default repeated call behavior:

- continue `x0`
- continue `d_dt`
- continue RNG state
- continue observer state
- keep the same requested `tspan` unless changed manually

## Bottom line

The current Python API is built around persistent simulator objects whose device state survives across calls.

The main continuation primitive is not the kernel itself. It is the wrapper's post-run `shift_x0()` call.

So the default model is:

- kernels advance `xf`, `tf`, `d_dt`, RNG state, and optional output buffers
- wrappers optionally promote `xf` into the next run's `x0`
- wrappers do not automatically advance the requested time window

For feature observers, there is one more layer:

- the kernel also rebases continuation-relevant observer time members by the true attained interval in `finalizeObserverData(...)`

So the current system really has three different continuation notions:

- state continuation: yes, by default
- observer-internal relative-time continuation: yes, for features
- absolute-time continuation: only if the caller manages `tspan` appropriately

For autonomous systems, the first notion is usually the important one.

For non-autonomous systems, the third notion becomes essential, and ideally the next window should be based on the true final time rather than the nominal requested endpoint.
