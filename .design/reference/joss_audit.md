# JOSS Audit

Purpose: assess whether the current public package story is strong enough for JOSS and identify the missing evidence or repo-surface work.
Read when: editing `paper/paper.md`, planning citation metadata, or deciding whether public claims are mature enough for publication-facing wording.
Update when: the paper story, public evidence base, citation state, or release/presentation polish changes materially.

## Bottom line

clODE looks plausibly publishable in JOSS, but not yet as a strong submission in its current public presentation. The codebase now has enough substance, age, license, docs, tests, and CI to be in scope. The main weakness is not whether this is research software. The weakness is whether the public story is current, differentiated, and well evidenced.

## What clODE can credibly claim today

- Public development history since 2019, which comfortably clears JOSS's public-history bar.
- MIT-licensed Python package with a current PyOpenCL runtime dependency and automated CI coverage.
- A distinctive workflow around very large ensembles of ODE solves on OpenCL devices, with one work-item per solve and compile-time specialization by precision, stepper, observer, and problem dimensions.
- On-device observer/feature extraction that can avoid storing full trajectories when users only need summary statistics or events.
- Multiple front doors into the solver stack: OpenCL RHS source, Python-authored RHS conversion, and XPP parsing/interoperability.
- A cleaner PyOpenCL-only architecture than the current paper draft describes.

## Strongest publication angles

1. Large parameter-sweep and ensemble workflow rather than "yet another ODE solver".
2. On-the-fly feature extraction as a memory and throughput story, not just a GPU acceleration story.
3. OpenCL and PyOpenCL portability for researchers who do not want a CUDA-only stack.
4. XPP-to-GPU bridge for existing dynamical-systems workflows.
5. A pragmatic niche between general-purpose CPU integrators and ML-oriented differentiable ODE toolchains.

## State-of-the-field argument to sharpen

The paper should not try to win on a raw claim that GPU ODE solving is otherwise unavailable. A stronger comparison is:

- SciPy and Matlab style solvers: strong general-purpose integrators, weak for massive ensemble sweeps and online feature extraction.
- DifferentialEquations.jl and diffeqpy: much broader solver coverage, but a different ergonomics and deployment story; clODE is narrower and more workflow-specialized.
- torchdiffeq, torchode, and nearby JAX or ML stacks: focused on differentiable or batched tensor ODE workflows, not the same observer and feature-extraction use case.
- XPPAUT: strong model authoring and analysis ecosystem, not a GPU ensemble engine.

The submission should explain why clODE exists alongside these tools instead of implying that they do not.

## Gaps before this is a strong JOSS contender

- The evidence base is still thin. Public materials now have an updated repo/package overview, docs home split, contributor guide, refreshed paper framing, and performance-notes scaffolding, but they still need benchmark data and comparison cases tied to real research workflows.
- The scholarly impact story needs updating. The existing draft cites older work and current applications, but the public submission package should show present-day use, citations, or concrete near-term users.
- The contributor and review signals are improved but still incomplete. CI, docs, security policy, license, and a contributor guide are present, but there is still no citation file and the public repo metadata could be tighter.
- The release story is a little rough. Git tags stop at `v0.8.1` while the package reports `0.10.0`, which makes the public maintenance story look less polished than the codebase now is.
- Core correctness work still matters. The live continuation and observer-state issues are not necessarily a JOSS blocker by themselves, but they sit underneath important future claims around chunking, continuation, and reproducibility.

## Recommended path

Short term:

- Keep the paper, README, and docs surfaces aligned around the current PyOpenCL-only package and its workflow niche.
- Build a comparison table centered on workload shape, feature extraction, portability, and model-ingestion paths.
- Turn the current performance-notes scaffolding into one or two benchmark scripts and result pages that reviewers can rerun without archaeology.
- Add citation metadata once the paper and software citation details are ready.

Medium term:

- Close the continuation and state-semantics gap so the public story is not ahead of the implementation.
- Keep the scope tight. "High-throughput ensemble ODE simulation with online observers" is a more defensible JOSS claim than a broad claim about being a fully general solver platform.

## Verdict

JOSS is realistic. The main work is not inventing a publishable angle from scratch. The main work is updating the public story so it matches the current code and demonstrates why this narrower workflow-oriented package matters.
