# .design Maintenance

Purpose: keep maintainer docs fast to scan, current, and worth reading for both humans and coding assistants.
Read when: you are changing root `.design` docs, routing rules, or the structure of the design tree.
Update when: file roles, maintenance rules, or assistant workflow expectations change.

## Goals

- Minimize startup reads for a typical PR.
- Keep one authoritative home for current package state, the current target, the active backlog, and the longer roadmap rationale.
- Separate live docs, deep dives, history, and scratch work.
- Favor updates over proliferation.

## Token-Efficient Defaults

- The root `.design/` surface should answer common contributor questions without forcing `reference/` for ordinary work.
- Keep each root doc's first screen high-signal: current facts, current target, or current backlog only.
- One durable fact should have one authoritative home. Link instead of copying paragraphs.
- Move lookup tables, historical closeouts, and long comparative analysis out of root docs.
- If a section mainly explains landed history, archive it or compress it into one sentence.
- Prefer a short summary plus a pointer over repeating the same context across `package_state.md`, `next_pr.md`, `ideas.md`, and `development_roadmap.md`.

## Directory Contract

- Root `.design/`: live, high-signal docs assistants should consider first.
- `.design/reference/`: focused deep dives and audit notes. Task-specific, not default reading.
- `.design/archived/`: historical context only.
- `.design/tmp/`: scratch work and copied artifacts with no authority.

## Live File Contracts

- `README.md`: routing only. Keep it short and explicit about where to stop reading.
- `package_state.md`: factual current-state map and contributor routing. No speculative backlog.
- `next_pr.md`: one active implementation target with a narrow scope and explicit acceptance criteria.
- `ideas.md`: terse open backlog only. Keep items to one line when possible, and move completed lines out once their durable outcome is reflected elsewhere.
- `development_roadmap.md`: medium-lived rationale and ordering, not a task tracker or closeout log.
- `reference/README.md`: index of deep dives with clear "read when" guidance.

## Completed Backlog Items

- When an `ideas.md` item lands, first update the smallest live doc that should carry the durable result: usually `package_state.md`, `next_pr.md`, `development_roadmap.md`, or the relevant reference note.
- After that live-doc update, remove the completed line from `ideas.md` and move its original wording into a dated archive cluster under `.design/archived/`.
- Keep the archive entry historical. Do not let archived completed items become the only place that current contributors have to learn a still-relevant package fact.
- Prefer one dated closeout cluster per cleanup pass over sprinkling completed-item history across multiple root docs.

## Inbox Processing

- `ideas.md` `## Inbox` is a temporary holding area, not a second backlog.
- Use the `design-ideas-inbox` skill when rough inbox entries need to be routed into the maintained `.design` structure.
- Most inbox items should become normalized backlog lines in an existing `ideas.md` section.
- Only promote an inbox item into `next_pr.md`, `package_state.md`, `development_roadmap.md`, or `.design/reference/` when it materially belongs there.
- If an item is too vague to place, rewrite it into a smaller explicit question instead of silently dropping it.

## Design Audits

- Use the `design-doc-audit` skill when checking whether `.design` is stale, internally inconsistent, or out of sync with the live package state.
- Priority changes, blocker cleanup, and stale-item removal should be evidence-driven and conservative.
- Prefer `report-only` audits for broad priority questions; use cleanup mode only for low-risk edits.

## When To Create A New Doc

Create a new long-lived doc only if all of these are true:

- The information does not fit cleanly into an existing live or reference doc.
- The topic is likely to matter across more than one PR.
- Expanding an existing doc would make that doc materially harder to scan.

If the note is narrow and temporary, use `.design/tmp/` instead.

## Writing Rules

- Start new long-lived docs with three short lines: `Purpose:`, `Read when:`, and `Update when:`.
- Keep openings short; front-load decisions, constraints, and status.
- Prefer bullets and short sections over long prose.
- Avoid restating package facts already covered in `package_state.md`.
- Avoid restating backlog or roadmap context already covered in `ideas.md` or `development_roadmap.md`.
- Keep root docs task-oriented: current fact, active plan, or sequencing rationale.
- Link to the smallest useful supporting doc, not a broad file dump.
- Mark historical statements explicitly before moving them to `.design/archived/`.

## Update Triggers

- Package or module moves, canonical import changes, or build/runtime assumption changes: update `package_state.md`.
- Current implementation focus changes: update `next_pr.md`.
- Priority or sequencing changes: update `ideas.md` and, if the rationale changed, `development_roadmap.md`.
- Docs, paper, public-surface, or testing strategy decisions change: update the relevant file in `.design/reference/`.
- A reference note becomes active delivery scope: summarize the live decision back into `next_pr.md` or `package_state.md`.
- A note stops guiding current work: move it to `.design/archived/`.
- Scratch notes that remain useful after the PR: promote them to `.design/reference/` or move them to `.design/archived/`; otherwise delete them.

## Reference Notes

- Treat each reference note as one of three types and make that obvious near the top: stable guardrail, current design audit, or public-surface plan.
- Stable guardrail notes should capture durable boundaries, constraints, or workflows and should avoid carrying active `next_pr.md` sequencing once the root docs already cover it.
- Current design audits may name the live design question, but once a rollout lands, move the proof-path history, earlier family rankings, or superseded sequencing into `.design/archived/`.
- Keep each reference doc single-topic.
- Put recurring task keywords near the top so assistants can route quickly.
- If a reference note grows past ordinary scan length, add a first-screen `## Fast path` or equivalent stop-guidance block near the top.
- If a reference doc grows past roughly 200 lines, split it or archive stale landed-history sections unless nearly every section is still active decision support.
- When updating a long reference note, refresh the summary or fast-path section first and move stale comparative or landed-history material to `.design/archived/` instead of letting the opening screen grow.
- If a reference doc grows broad enough to answer multiple unrelated questions, split it.
- If a reference doc becomes stale because the live package changed, update it or archive it. Do not let it silently rot.

## Archived Notes

- Archive by dated folder or by a clearly named closeout cluster.
- Add or update `.design/archived/README.md` when adding a new archive cluster.
- Do not rewrite archived notes to reflect the present unless the change is a signpost such as "historical only; current source of truth is ...".

## PR Checklist

- Does the root `.design/` still contain only live docs?
- Did you update the smallest authoritative doc instead of adding a near-duplicate?
- Did newly completed `ideas.md` lines move to an archive cluster after the live docs were updated?
- Are moved docs discoverable from `README.md` or `reference/README.md`?
- Are broken links or stale file paths fixed in the same pass?
- Did you clean up `.design/tmp/` if you used it?
- If assistant workflow changed, did you update the relevant `.github/instructions/*.instructions.md` file?
