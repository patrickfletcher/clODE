# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Target refresh required (previous target landed)

## Current Status

- The previous active target (clean-slate semantic observer surface) is now landed in code and tests.
- This file still carried the historical rollout plan and no longer reflected active work.

## Immediate Maintainer Action

- Select one current open item from `.design/ideas.md` as the next active implementation target.
- Rewrite this file to the selected scope with narrow acceptance criteria.

## Guardrails For The Refresh

- Do not re-open completed observer-retirement rollout work as active scope.
- Keep the next target small and directly testable.
- Treat `.design/package_state.md` as the current-state source of truth when drafting the refreshed target.
