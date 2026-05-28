---
name: design-ideas-inbox
description: 'Process the `.design/ideas.md` Inbox. Use when the user asks to process the ideas inbox, triage new design ideas, clean up the Inbox section, or fit rough notes into the `.design` docs framework.'
argument-hint: '[optional: process the full inbox, or name a topic or item to focus on]'
user-invocable: true
---

# Design Ideas Inbox

Turn rough items in `.design/ideas.md` `## Inbox` into well-placed, normalized `.design` updates.

## When To Use

- The user asks to process the ideas inbox.
- The user wants rough design ideas triaged into the maintained `.design` structure.
- The inbox has accumulated ad hoc notes that should become backlog items, roadmap rationale, or deeper reference notes.
- A PR added new inbox items and they need to be deduplicated, clarified, or promoted.

## Required Reads

Read these first:

1. `.design/README.md`
2. `.design/MAINTENANCE.md`
3. `.design/ideas.md`

Read more only if the inbox items need them:

- `.design/package_state.md` for current facts and boundaries
- `.design/next_pr.md` for the active target
- `.design/development_roadmap.md` for broader sequencing or rationale
- `.design/reference/README.md` for topic-specific deep dives

Use the [placement guide](./references/placement-guide.md) when an item could plausibly fit more than one doc.

## Default Bias

- Most new ideas should land in the closest existing section of `.design/ideas.md`.
- Do not create a new doc for a one-line backlog item.
- Do not update `package_state.md` with speculative proposals.
- Do not rewrite `next_pr.md` unless the user is clearly reprioritizing the active target.
- If an inbox item clearly describes already-landed work, update the smallest live doc only if needed and avoid turning it back into an open backlog line.
- Use `design-doc-audit` instead if the user wants repo-wide stale-item cleanup, priority verification, or blocker auditing.

## Procedure

1. Reread the current `## Inbox` before editing. Preserve recent user edits.
2. For each inbox item, decide whether it is a backlog item, live fact, roadmap rationale, deep-dive note, historical context, or scratch work.
3. Search `ideas.md` for an existing near-duplicate before adding a new item.
4. Prefer merging into an existing `ideas.md` section. Normalize backlog entries to the repo format: `[ ] P0/P1/P2 item. depends: ... blocks: ... refs: ...`.
5. Only update another `.design` doc when the item genuinely belongs there under the placement rules.
6. After promoting or merging an inbox item, remove the original inbox line.
7. If an item is too vague to place without inventing substance, leave it in `## Inbox` but rewrite it as a smaller explicit question or placeholder.
8. Fix stale paths or references introduced by the edits, then validate the touched markdown files.

## Success Criteria

- The inbox shrinks or becomes empty.
- No idea is silently dropped; each item is merged, moved, rewritten, or intentionally left as an explicit open question.
- `ideas.md` stays terse and sectioned by domain.
- The root `.design/` surface stays limited to live docs.
- Any new long-lived doc follows the `Purpose:`, `Read when:`, `Update when:` contract.