# AGENTS.md

## Purpose

Shared collaboration rules for Claude, Codex, and any other agent working in this repository.

## Roles

- **Claude**: primary planning and implementation agent unless the user explicitly assigns work differently
- **Codex**: review-focused secondary agent by default; may make scoped follow-up edits when asked
- **User**: final authority on scope, sequencing, naming, and merge decisions

## Planning And Implementation Boundary

- Do not move from planning into implementation unless the user explicitly says to proceed.
- Original prompts under `.github/prompts/` are reference inputs and must not be edited.
- Revised, working plans live under `plan/phaseN.md` and should stay self-contained.

## Handoff Rules

- Leave the repository in a state that another agent can understand without chat history.
- When work changes plan status, assumptions, or scope, update the active `plan/phaseN.md` before ending the session.
- Preserve prior agent notes. Append new review/response sections instead of rewriting history unless the user asks for consolidation.
- Record what was validated, what was not validated, and any blockers that remain.

## Review Thread Convention

- Reviews and responses are append-only sections in the relevant plan or working document.
- Use Roman numerals and keep the number paired across a review and its response.
- Standard pattern:
  - `## Codex Review I`
  - `## Claude Response I`
  - `## Codex Review II`
  - `## Claude Response II`
- If Claude performs the review first, invert the names:
  - `## Claude Review I`
  - `## Codex Response I`
- Do not rename or repurpose an existing numbered section after it exists.
- If a prior review needs correction or follow-up, add the next numbered review instead of editing the old one.

## Review Content Expectations

- Reviews should lead with findings ordered by severity.
- Responses should address each finding clearly: fixed, accepted, deferred, or rejected with reason.
- If a response includes code or document changes, name the affected files and any validation that was run.
- Keep each section readable on its own so a fresh agent can resume immediately.

## Ownership Defaults

- Claude owns the mainline implementation path on the active phase branch unless the user delegates otherwise.
- Codex should prefer review, validation, and documentation updates unless the user asks for implementation work.
- Neither agent should revert or rewrite the other agent's work without explicit user direction.

## Push And Merge

- Use `~/bin/push-ncar` instead of `git push`.
- No phase should move into implementation until explicitly authorized.
- No subsequent phase begins until the prior phase PR is reviewed and merged.
