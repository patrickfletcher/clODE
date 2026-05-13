# Design Inbox Placement Guide

Purpose: route rough inbox ideas into the right `.design` home without duplicating or bloating the live docs.
Read when: an inbox item could reasonably fit more than one `.design` file.
Update when: the `.design` structure or inbox-routing rules change.

## Default Bias

- Most proposals belong in the closest existing section of `.design/ideas.md`.
- Keep the root `.design/` docs small and high-signal.
- Promote enduring decisions, not raw brainstorming.
- Prefer updating an existing item or doc over creating a near-duplicate.

## Routing Rules

| Idea type | Default home | Use when | Avoid when |
| --- | --- | --- | --- |
| Rough proposal or follow-on work | `.design/ideas.md` | The item can stay as a concise backlog line | It needs multi-paragraph rationale or an already settled decision |
| Current active target | `.design/next_pr.md` | The user is clearly changing or confirming the active implementation target | The idea is only a backlog candidate |
| Current package or repo fact | `.design/package_state.md` | The statement is already true in the live package | The item is still a proposal |
| Broader sequencing or rationale | `.design/development_roadmap.md` | A one-line backlog entry is not enough to preserve the reasoning | The idea can be captured tersely in `ideas.md` |
| Single-topic deep dive or policy note | `.design/reference/*.md` | The topic will likely matter across multiple PRs and does not belong in a live root doc | The note is narrow, temporary, or just a one-line backlog item |
| Historical closeout or obsolete rationale | `.design/archived/` | The note is useful for archaeology but no longer drives current work | The note still guides current implementation |
| Scratch reproduction or copied artifact | `.design/tmp/` | The material is temporary and still being worked out | The content is an enduring rule or decision |

## Processing Checklist

- Search for an existing item first and merge if the new idea is a near-duplicate.
- Normalize backlog lines to the repository format only when the extra fields are meaningful.
- Keep `ideas.md` terse; move detail elsewhere only when that detail will keep paying for itself.
- Leave an item in `## Inbox` only if categorizing it would require making up intent that is not present in the note.
- When an inbox item changes the meaning of another live doc, update that doc in the same pass.
- If you create a new long-lived note, add the `Purpose:`, `Read when:`, and `Update when:` header lines.
