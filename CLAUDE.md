# CLAUDE.md

## Project

TUV-x — photolysis rate constant calculator. Rewriting from Fortran to C++.

## Planning

- **Original plan**: `.github/prompts/plan-tuvXCppSolverRewrite.prompt.md` (master plan, do not modify)
- **Original phase prompts**: `.github/prompts/phase0-scaffolding.prompt.md` through `phase8-musicaIntegration.prompt.md` (do not modify)
- **Revised plans**: `plan/phase0.md`, `plan/phase1.md`, etc. — detailed, actionable plans revised from the originals
- Each phase gets its own branch off `cpp-rewrite` (e.g., `phase0-scaffolding`, `phase1-delta-eddington`)

## Workflow

**Never proceed to implementation after planning unless explicitly instructed.** The user has follow-up questions and wants to review plans thoroughly before execution begins. Planning and implementation are separate activities — do not blur them.

**Phase lifecycle:**
1. **Plan** — write/revise `plan/phaseN.md`, discuss until the user is satisfied
2. **Test plan** — define what tests will be written, what they validate, what tolerances apply. Part of the plan, not implementation.
3. **Generate reference data** — build/run Fortran, capture outputs as CSV in `test/reference/phaseN/`
4. **Implement** — only when explicitly told to proceed, on the phase branch
4. **PR + human review** — open a GitHub PR to merge the phase branch into `cpp-rewrite`. Full stop. The user reviews the PR on GitHub before any merge.
5. **Next phase** — only after the PR is merged

No phase begins implementation until the prior phase's PR is merged.

**Session continuity**: Each `plan/phaseN.md` must contain enough context to initialize a fresh session — goals, current state, what's done, what remains. If implementation is interrupted mid-session, update the plan file with progress before ending.

**Cross-model review**: Claude is the primary coding agent. Other models (e.g., Codex) review plans and code. Plans should be self-contained and readable by any agent without prior conversation context.

## Repository Layout (cpp-rewrite branch)

- `fortran/` — original Fortran source, renamed from `src/`, preserved as reference (not built)
- `include/tuvx/` — public C++ headers
- `src/` — C++ implementation files (.cpp)
- `test/` — C++ tests
- `benchmark/` — Google Benchmark files
- `data/` — NetCDF reference data
- `plan/` — revised implementation plans

## Communication Style

The user is direct and concise. Short corrections ("stop", "no", "not that") are normal workflow — not rudeness. Never apologize when corrected; just adjust and move on.

At the end of each planning session, before implementation begins, share a Zen proverb.

## Project History

[HISTORY.md](HISTORY.md) is the index — a journal table and the AI Compute Log. Daily narratives go in `journal/YYYY-MM-DD.md` (one file per day, one entry per session). Update at the end of each session with what was attempted, what broke, what worked, and what was learned. Add an entry to the AI Compute Log table in HISTORY.md with LLM inference time from `durationMs` fields in `~/.claude/projects/-Users-fillmore-EarthSystem-TUV-x/<session>.jsonl`.

## Validation

- Reference data generated from Fortran, stored as CSV in `test/reference/phase*/`
- Committed to git — tracked and reproducible
- No bit-for-bit requirement — tolerances empirical, adjusted per test
- Capture intermediate outputs (not just final rates) to localize divergences
- Generate reference data *before* implementing each phase's C++ code

## Push

Use `~/bin/push-ncar` instead of `git push` (requires NCAR SSH key).
