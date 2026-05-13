---
name: design-doc-audit
description: 'Audit whether the `.design/` docs are up to date. Use when the user asks to audit `.design`, verify priorities or blockers against package state, find stale backlog items, clean stale docs, or check whether planning docs still match the live package.'
argument-hint: '[optional: report-only, apply-cleanup, or a specific doc or section to audit]'
user-invocable: true
---

# Design Doc Audit

Audit the live `.design/` framework for staleness, contradictions, duplicate backlog items, outdated blockers, and mismatches against the current package state.

## When To Use

- The user asks whether `.design` is current.
- The user wants stale items, outdated blockers, or completed work cleaned up.
- The user wants priorities checked against the live package state.
- The user wants `.design/ideas.md`, `.design/next_pr.md`, `.design/package_state.md`, and `.design/reference/` checked for consistency.
- A large PR may have invalidated planning notes and the user wants a repo-local documentation audit.

## Modes

- `report-only`: find issues and propose edits without changing files.
- `apply-cleanup`: make low-risk cleanup edits directly and report any higher-risk judgment calls separately.
- targeted audit: focus on one doc, one section, or one topic area when the user names it.

Default to `report-only` unless the user explicitly asks for cleanup or fixes.

## Required Reads

Read these first:

1. `.design/README.md`
2. `.design/MAINTENANCE.md`
3. `.design/package_state.md`
4. `.design/next_pr.md`
5. `.design/ideas.md`
6. `.design/development_roadmap.md`
7. `.design/reference/README.md`

Then read only the focused reference notes, archived notes, or code files needed to evaluate the claimed issue.

Use the [audit checklist](./references/audit-checklist.md) when deciding what to verify and what can be cleaned automatically.

## Evidence Rules

- Treat `package_state.md` as the authoritative current-state narrative unless the live code clearly contradicts it.
- Use the current codebase and current tests as stronger evidence than stale planning language.
- Treat archived notes as historical context only.
- Do not rewrite priorities or blockers based on guesswork.
- If a supposed stale item needs broad product judgment rather than repo evidence, report it instead of silently changing it.

## Default Bias

- Prefer conservative cleanup.
- Auto-fix only low-risk issues: duplicates, clearly completed items, obviously removed blockers, broken internal paths, stale wording that conflicts with current live docs, or reference notes that now point to moved files.
- For priority changes, `next_pr.md` changes, or broader roadmap reshaping, report the recommendation unless the user explicitly asked for reprioritization.

## Procedure

1. Identify the audit scope: whole `.design`, one file, or one section.
2. Compare live planning docs against `package_state.md`, `next_pr.md`, current code, and any directly relevant reference notes.
3. Check for duplicates, stale `depends:` or `blocks:` relationships, completed items left open, and contradictions between the live docs.
4. Distinguish low-risk cleanup from high-risk judgment:
   - low-risk: duplicate backlog entries, resolved blockers with direct evidence, stale paths, outdated references, obviously completed cleanup items
   - high-risk: reprioritizing whole workstreams, replacing the active target, or removing items that are still conceptually valid but merely delayed
5. In `report-only` mode, return findings first, ordered by severity or confidence.
6. In `apply-cleanup` mode, make only the low-risk edits automatically, then report any remaining judgment calls.
7. Validate all touched markdown files after editing.

## Success Criteria

- The audit clearly separates facts, low-risk cleanup, and recommendation-level judgment.
- `.design` docs become more internally consistent and less stale.
- The inbox skill remains focused on routing rough new ideas, not repo-wide planning audits.
- The root `.design/` surface stays small and authoritative.