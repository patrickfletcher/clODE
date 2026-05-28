# Observer Threshold-Family Background

Purpose: preserve the longer threshold-family rollout argument and early observer-family ranking that were moved out of `observer_concept_audit.md` once they stopped being the main live decision surface.
Read when: you need the historical reasoning behind the threshold-family rollout, the earlier family-ranking argument, or the more detailed explanation of why the threshold and neighborhood families were split the way they were.
Update when: this archive signpost becomes inaccurate. Do not rewrite the historical argument to reflect current priorities; update the live observer note instead.

Historical only: the current source of truth for live observer semantics and the active bundle question is `.design/reference/observer_concept_audit.md`.

## Near-term candidate family ranking at the time of the threshold-family rollout

| Candidate family | Current evidence | Near-term fit | Design read |
| --- | --- | --- | --- |
| Absolute threshold crossing in state-variable units | now landed as `Observer.threshold_crossing` with a lean event stream and runtime direction selection | landed first proof target | confirms that not every observer-family extension needs schema-changing specialization |
| Warmup-derived normalized threshold crossing | now landed as `Observer.normalized_threshold_crossing` with a lean one-stream readout and warmup-derived threshold parameterization | landed second threshold proof target | separates threshold parameterization from Schmitt topology without overloading the Schmitt family |
| Absolute Schmitt trigger | now landed as `Observer.schmitt_trigger` with a lean up/down timestamp readout and thresholds expressed directly in event-variable units | landed follow-on proof target | confirms that trigger topology and threshold parameterization can vary independently while still sharing one public config model |
| Local extremum family | now landed as `Observer.local_extremum` with `ExtremumPolarity` selection while `Observer.local_max` remains the retained heavier legacy workflow | landed adjacent proof target | confirmed that extremum polarity can stay runtime-side while the lean event schema stays fixed |
| Generalized hyperplane crossing | conceptually matches the useful part of `solve_ivp`-style scalar events | later | promising long-term trigger umbrella, but it needs a declaration for scalar trigger functions or hyperplane parameters that the current runtime path does not yet package cleanly |
| Neighborhood return triggers | now landed as `Observer.neighborhood_return` for the lean specialized path while `Observer.neighbourhood_2` remains the retained heavier legacy workflow | landed specialized follow-on | keeps `nhood2`'s better-behaved trigger semantics available without forcing the broader periodicity bundle onto the lean path; `nhood1` remains the clearest retirement candidate |

## Trigger and parameter taxonomy retained from the earlier live note

### Static trigger families

- Summary reducers and `local_max` do not need a first pass to discover internal trigger parameters.
- Their runtime knobs are mostly scalar selectors, thresholds, or counters.
- These families were the easiest place to separate trigger geometry from retained outputs or running summary groups.

### Online-discovery trigger families

- `neighbourhood_1` does not run a formal warmup pass, but it still discovers internal trigger parameters during the live pass.
- The anchor point `x0` is set from the first local minimum in `eVarIx`, and the normalization basis in `eventFunction(...)` uses trajectory ranges that are still changing over time.
- That made `nhood1` a distinct design case rather than just another one-pass detector.

### Warmup-derived trigger families

- `threshold_2` and `neighbourhood_2` both need a first pass to establish internal geometry before the main event pass starts.
- That was a good fit for a two-pass family concept, but the live code still mixed warmup-derived geometry with live-pass feature bundles and retained-output policy under one built-in name.

### Threshold topology versus threshold parameterization

- The threshold discussion became clearer once it was split along two independent axes instead of treated as one linear family ladder.
- The first axis was trigger topology: a single-boundary directional crossing versus a two-boundary Schmitt-style state machine.
- The second axis was threshold parameterization: absolute state-variable units versus warmup-derived fractional position inside the observed event-variable range.
- `Observer.threshold_crossing` occupied the absolute plus single-boundary corner.
- `Observer.schmitt_trigger` occupied the absolute plus Schmitt-style corner.
- `Observer.normalized_schmitt_trigger` occupied the warmup-derived plus Schmitt-style corner.
- `Observer.normalized_threshold_crossing` landed in the warmup-derived plus single-boundary corner with the same lean readout style as `Observer.threshold_crossing`.

### Why threshold crossing and Schmitt triggering stayed separate user-facing concepts

- Threshold crossing answers a one-boundary question: when does `x` pass this level, in which direction?
- Schmitt triggering answers a stateful two-boundary question: when does `x` enter and leave a latched upstate without counting chatter as repeated events?
- A Schmitt trigger with `xUp == xDown` and both directions can collapse to a degenerate single-boundary case, but that equivalence is better treated as an implementation detail than as the public UX model.
- Separate public concepts kept the docs legible, reduced inert parameters on the simple path, and made it easier to explain why a lean threshold family and a heavier hysteresis family naturally wanted different default readouts.

### What the threshold-family config seam proved

- Family-specific config objects were enough to improve the UX without changing the executor or kernel contracts.
- Once one config class was shared across an absolute and a warmup-derived pair, the config type could no longer determine the observer family by itself. The selected observer remained the semantic source of truth.
- Recurring knobs such as `min_amp` did not justify a giant public bundle by themselves; the better rule was to surface a knob where it was active for that family.

### Why neighborhood-return triggers remained valuable

- `nhood2` was materially better behaved than `nhood1` because the warmup pass fixed the trajectory ranges first, and the live pass only accepted `x0` after `x[eVarIx]` dropped below a threshold derived from that warmup.
- That mattered when `eVarIx` was a slow variable used to identify a meaningful phase of a bursting or otherwise multiscale oscillation.
- The family was not just a threshold surrogate. Its intended semantics were to define a small neighborhood around a good full-state anchor point and then detect departure after a return to that neighborhood.

## Additional design implications retained for historical context

### Why parameter source became a first-class design axis

- Pass count alone did not explain what kind of declaration a family needed.
- The more useful axis was where the trigger geometry got its internal parameters: nowhere, from static runtime knobs, from live-pass discovery, or from a warmup pass.

### Why `eVarIx` versus `fVarIx` had to stay explicit

- The observer catalog could not collapse event geometry and feature readout onto one implied primary variable.
- A real workflow is to trigger on a slower or more reliable reference variable while measuring period, amplitude, maxima, or other summaries on a different faster variable.

### What the minimum absolute-threshold family proved

- The landed `threshold_crossing` slice did not clone `threshold_2`'s full readout bundle into a one-pass shell. That smaller proof target was the right call.
- Trigger geometry: scalar crossing on `eVarIx` in absolute state-variable units.
- Required runtime knobs: `eVarIx`, crossing direction, one absolute threshold, `maxEventCount`, and retained-timestamp capacity.
- Retained outputs: event count plus optional crossing timestamps, with no full copied event-state snapshot.
- Running summaries: none in the first slice; the family stayed lean instead of dragging along the all-state plus all-aux bundle.

### What `threshold_2` still taught beyond the lean threshold families

- Warmup-derived normalized thresholds remained useful when the oscillation amplitude of `eVarIx` changed substantially across parameter sweeps.
- The current `threshold_2` implementation combined warmup-derived fractional placement, Schmitt-style hysteresis, and optional derivative gates.
- The derivative gates appeared most justifiable on noisy or stochastic traces where value thresholds alone still overcount shallow or ripple-driven crossings.

### Constraint from the existing stub

- The semantic `observer_threshold_crossing.clh` file became the lean proof target for absolute threshold observers, while the older `observer_threshold_1.clh` file remained only historical evidence.
- Wiring the older stub directly would have pulled clODE back toward a monolithic all-state snapshot layout instead of extending the cleaner declaration boundary established by the summary slice.

### Why `neighbourhood_1` was treated as a warning sign

- `nhood1` showed that “one pass” is not automatically the simplest design.
- Because its anchor point and normalization basis were discovered during the same pass that used them for events, the trigger geometry could move as the solve proceeded.

### Why two-pass families should not be grouped only by warmup

- `threshold_2` and `nhood2` both need a warmup pass, but they use that warmup for different reasons.
- `threshold_2` derives scalar trigger thresholds from extrema on one observed variable.
- `nhood2` derives normalization ranges from warmup, then discovers the anchor point `x0` online in the second pass only after a threshold-qualified crossing in `eVarIx`.

### Earlier observer-design split proposed in the live note

The earlier live note proposed five observer-design axes that mattered more than the monolithic built-in names:

1. trigger geometry
2. parameter source: static, live-pass discovery, or warmup-derived
3. running-summary bundle
4. retained sparse-output policy
5. selected variable scope for any all-state or all-aux summary bundles

That historical argument is preserved here because it helped motivate the later summary-family and threshold-family cleanup, even though the current live note now focuses on the accepted-step solution-buffer and observer-bundle question instead.
