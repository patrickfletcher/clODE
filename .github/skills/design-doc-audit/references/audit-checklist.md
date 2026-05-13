# Design Doc Audit Checklist

Purpose: help the audit skill verify whether the live `.design` docs still match the current package and whether stale items can be cleaned safely.
Read when: auditing `.design` for stale items, removed blockers, duplicate backlog entries, or doc drift.
Update when: the `.design` framework or audit rules change.

## Audit Questions

### Current-state drift

- Does `.design/package_state.md` still describe the current package layout and compatibility surface?
- Does any root `.design` doc contradict `package_state.md` without evidence from the code or an explicit planned change?
- Do live docs still point at moved or archived files?

### Active target and priorities

- Does `.design/next_pr.md` still describe the current highest-priority active target?
- Do top-priority items in `.design/ideas.md` still support or follow from `next_pr.md`?
- Has a priority become obviously stale because its blocker is already resolved in the current package or in another live doc?

### Backlog hygiene

- Are there duplicate or near-duplicate items in different sections of `ideas.md`?
- Does an item marked with `depends:` or `blocks:` still have that relationship?
- Is any item clearly complete according to the live docs or current code?
- Is a backlog line carrying too much detail and better suited to `development_roadmap.md` or a single-topic reference note?

### Reference-note hygiene

- Does a note under `.design/reference/` still provide current decision support?
- Has a reference note effectively become active delivery scope and therefore need a summary in a live root doc?
- Has a reference note gone stale enough that it should be updated, split, or archived?

## Evidence Hierarchy

1. Current live code and tests.
2. `package_state.md` for current-state narrative.
3. `next_pr.md` for the active target.
4. `ideas.md` and `development_roadmap.md` for planned work.
5. `.design/reference/` notes for focused context.
6. `.design/archived/` for historical rationale only.

## Safe Automatic Cleanup

- Merge obvious duplicate backlog items.
- Remove or rewrite a blocker when a live doc or the current code clearly shows it is resolved.
- Fix moved paths, stale references, and outdated cross-links.
- Mark obviously completed housekeeping items as done or remove them if they were temporary cleanup notes.
- Tighten wording that conflicts with current live docs when the correction is factual rather than strategic.

## Report Instead Of Auto-Fixing

- Replacing the active target in `next_pr.md`.
- Broad reprioritization across sections of `ideas.md`.
- Removing an item because it seems less important without direct repo evidence.
- Changing roadmap sequencing where the tradeoff is architectural rather than factual.
