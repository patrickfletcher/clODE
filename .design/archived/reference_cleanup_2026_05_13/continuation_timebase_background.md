# Continuation Time-base Background

Purpose: preserve the debugging-level continuation analysis and alternative observer-time-model discussion that no longer belongs in the slim live note.
Read when: you need bug archaeology for the landed continuation fix or want the detailed reasoning behind the solver-owned absolute-time choice.
Update when: a signpost to the current live continuation note needs correction.

Historical only. The current live summary is `.design/reference/continuation_timebase_note.md`.

## Earlier correction to the explanation

For the running-mean observers, `od->t_last` absolutely matters.

The relevant update is:

```c
realtype dt = *ti - od->t_last;
od->t_last = *ti;
realtype elapsedTime = *ti - od->t_start;
mean = runningMeanTime(mean, value, dt, elapsedTime);
```

So the earlier claim that the bug was fundamentally about rebasing `od->t_start` was incomplete. The old bug was that only a subset of the observer's continuation-relevant times was rebased.

For the basic running-mean observers, rebasing both `od->t_last` and `od->t_start` backward by the attained interval would have repaired the observer-side running-mean algebra.

## Two coherent observer-time models that were considered

### Model A: absolute-time observer state

- keep all persisted observer times in the attained absolute-time frame
- make the next kernel start from the previous attained final time
- exact continuation therefore requires advancing `t_span` between calls

This is the model the current implementation now uses.

### Model B: rebased relative-time observer state

- after each window, shift every persisted observer time left by the attained interval
- keep the next kernel in a local window frame again
- repeated `features()` calls could then continue observer state even if the nominal window is reused

For the basic or allvar running-mean observers, this would mean rebasing at least:

- `t_start`
- `t_last`

For event observers it would also mean rebasing all continuation-relevant stored times, for example:

- `tbuffer[]`
- `tLastEvent`
- `tLastMax`, `tLastMin`
- any other persisted time members that affect later updates or event detection

## Why the old behavior still was not a complete continuation model

Even if Model B were implemented correctly for observer state, it would only solve the observer-side time bookkeeping.

The kernel time passed into the RHS would still restart from the current requested `t_span[0]`. That means:

- autonomous systems could still look acceptable because only state continuation matters
- non-autonomous systems would still not be following the same absolute-time trajectory as one uninterrupted run

So the observer-time rebase model is viable as an observer-local continuation policy, but it does not by itself define the full simulator continuation semantics.

## Why Model A was chosen for the landed fix

The landed fix chose Model A because it gives one consistent absolute-time frame for:

- kernel time seen by the RHS
- observer continuation state
- attained final-time handoff between windows
- user-visible timestamps reported by observers

This is a cleaner base for non-autonomous correctness and for comparing split windows against one uninterrupted run.

## What would have been required to keep the old user expectation

If the public expectation is that repeated `features()` calls should "just continue" without any explicit `t_span` change, there were two different possible interpretations:

1. continue only the state and observer internals in a local-window sense
2. continue the actual absolute-time simulation as if one long run had been split

Only the second interpretation is compatible with non-autonomous correctness.

If the package wants repeated calls to behave like that by default, then some automatic `t_span` advancement policy is needed. A helper or explicit continuation-policy API is likely safer than more implicit rebasing rules.

## Fixed-step overshoot and attained final time

The attained final time can differ from the requested endpoint for fixed-step methods.

- In Model A, the next requested window should start from attained `tf`.
- In Model B, the observer rebase would also need to use attained `tf - tspan[0]`, not just the nominal requested duration.
