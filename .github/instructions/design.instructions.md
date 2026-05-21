---
applyTo: ".design/**"
description: Use when editing .design docs. Keep the root planning surface small, authoritative, and cheap to read.
---

# Design Doc Rules

- Root `.design/` files are the live planning surface. Keep them short enough that common work does not need `reference/`, `archived/`, or `tmp/` by default.
- `README.md` routes; `package_state.md` carries current facts and contributor routing; `next_pr.md` tracks one active target; `ideas.md` stays terse; `development_roadmap.md` keeps rationale and sequencing.
- Prefer updating an existing note over creating a near-duplicate.
- Move deep dives to `.design/reference/`, historical closeout to `.design/archived/`, and scratch material to `.design/tmp/`.
- Avoid repeating the same package fact, backlog item, or roadmap rationale across multiple root docs.
- New long-lived notes should start with `Purpose:`, `Read when:`, and `Update when:`.
- When touching a root doc, keep the opening screen high-signal and fix stale links or duplicate entries in the same pass.