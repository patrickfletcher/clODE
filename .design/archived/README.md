# Archived Notes

Purpose: historical rationale, migration history, and bug archaeology.
Read when: the live `.design` docs are insufficient and you need past reasoning, reproduction details, or removed-architecture context.
Update when: a new archive cluster is added or an archive description becomes inaccurate.

## Rules

- Archived notes are never the source of truth for the current package.
- Use archived notes to explain why something exists or to recover a prior investigation, then verify against the live code and root `.design` docs.
- Prefer signpost edits over rewriting archived content.

## Current Folders

- `backend_migration_history_2026_05_07/`: historical backend-overhaul planning and scope maps from the transition period.
- `pyopencl_cleanup_closeout_2026_05_08/`: closeout notes from the PyOpenCL migration and wrapper cleanup.
- `pre_backend_readiness_2026_05_05/`: pre-migration bug reproductions and readiness audits.
- `post_logging_cleanup_2026_05_12/`: logging closeout plus follow-up boundary and layout notes after the logging cleanup.
- `reference_cleanup_2026_05_13/`: historical spillover moved out of live reference notes so `.design/reference/` stays current and compact.
