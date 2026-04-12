# CLAUDE.md

## Project

TUV-x — photolysis rate constant calculator. C++ rewrite.

## Auxiliary Repository

**[TUV-x-cpp](https://github.com/NCAR/tuv-x-cpp)** (`../TUV-x-cpp`) is the companion sandbox repo containing:
- `fortran/` — original Fortran source (reference, not built)
- `plan/` — revised phase plans and review threads
- `journal/` — session narratives
- `HISTORY.md` — journal index and AI Compute Log
- `docs/`, `tutorial/`, `examples/`, `docker/`, `packaging/`, `etc/`, `tool/` — legacy infrastructure

Use TUV-x-cpp for planning, journaling, and Fortran reference. This repo (TUV-x) is the lean C++ codebase for CI/CD and PR review — reviewers need only clone this repo.

## Workflow

**Never proceed to implementation after planning unless explicitly instructed.**

**Phase lifecycle:**
1. **Plan** — write/revise plans in TUV-x-cpp `plan/phaseN.md`, discuss until satisfied
2. **Test plan** — define tests, tolerances. Part of the plan, not implementation.
3. **Generate reference data** — capture outputs as CSV in `test/reference/phaseN/` (stays in this repo)
4. **Implement** — only when explicitly told to proceed, on the phase branch
5. **PR + human review** — GitHub PR to merge the phase branch into `cpp-rewrite`. User reviews before merge.
6. **Next phase** — only after the PR is merged

No phase begins implementation until the prior phase's PR is merged.

**Session continuity**: Plans in TUV-x-cpp must contain enough context to initialize a fresh session. If implementation is interrupted, update the plan file with progress before ending.

## Agent Coordination

- Shared handoff rules live in `AGENTS.md`. Keep `CLAUDE.md` and `AGENTS.md` aligned when the collaboration protocol changes.
- Claude remains the primary coding agent unless the user explicitly assigns work differently.
- Codex is typically used for plan reviews, code reviews, validation checks, and targeted follow-up edits.

## Review Thread Naming

- Review threads in TUV-x-cpp `plan/*.md` are append-only.
- Use Roman numerals for paired review/response sections.
- Standard pattern: `## Codex Review I` / `## Claude Response I`, etc.
- Do not overwrite or rename an earlier numbered section. Add the next one instead.

## Repository Layout

- `include/tuvx/` — public C++ headers
- `src/` — C++ implementation files (.cpp)
- `test/` — C++ tests (including `test/reference/` baseline data)
- `benchmark/` — Google Benchmark files
- `data/` — NetCDF reference data
- `cmake/`, `CMakeLists.txt` — build system

## Communication Style

The user is direct and concise. Short corrections are normal workflow — not rudeness. Never apologize when corrected; just adjust and move on.

At the end of each planning session, share a Zen proverb. Record it in the journal entry (in TUV-x-cpp).

## Project History

Journals and history live in TUV-x-cpp. Update `HISTORY.md` and `journal/YYYY-MM-DD.md` there at the end of each session.

## Validation

- Reference data generated from Fortran, stored as CSV in `test/reference/phase*/`
- Committed to git — tracked and reproducible
- No bit-for-bit requirement — tolerances empirical, adjusted per test
- Capture intermediate outputs (not just final rates) to localize divergences
- Generate reference data *before* implementing each phase's C++ code

## Terminology

- **Radiator → Constituent**: The Fortran codebase uses "radiator." In radiative transfer theory the standard term is "constituent." Rename incrementally as each file is modified during each phase.

## Push

Use `~/bin/push-ncar` instead of `git push` (requires NCAR SSH key).
