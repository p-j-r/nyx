# 20. Version control: a branch per task, never on main, a jj workspace per working session, status files as handoffs

Date: 2026-11-03

## Status

Accepted (`CLAUDE.md`, commit `a047506`; paul wants jj used, 2026-11-03).

## Context

nyx is built by many short working sessions, human and AI, on one machine, with `main` local and unpushed. Work was lost or tangled when sessions shared one working copy: a thread "switched the shared `nyx-jj` working copy mid-session" (`MISSING_FEATURES.md` section 6), uncommitted `geom` and `ai/` trees lived only in the main checkout for days (`CAPABILITIES_MEMO.md`, section 0.1).

## Decision

* One branch (jj bookmark) per task, named for the task. **Never commit to `main` directly.** Commit messages say what changed and why in one line.
* Each working session uses its **own jj workspace** (`jj workspace add ../<name>`), not the shared clone `nyx-jj`.
* `docs/*_STATUS.md` files are the handoff between sessions: read the relevant one before starting, update it before finishing.
* Nothing is left only in an uncommitted working tree or a side folder.
* `runs/`, `dist/` and results are git-ignored and never committed. `NYX_WRITE_GOLDEN=1` rewrites are committed only with a statement that the kernel is known to be right (ADR 10).
* Build artefacts: `target/` is local; `ai/` has its own `target` and lockfile.

## Consequences

* `main` is a local branch with no remote; "pushed" is not available as a safety net. Git bundles of earlier branches (`kernel-port*.bundle`, `physics-codes.bundle`) sit in `nyx-jj/target`; they are the only off-branch copies I saw.
* Two checkouts exist beside the jj clone: `nyx` (the git main checkout, with two untracked docs, R24) and `nyx-main-checkout-backup-2026-11-01`.
* This architecture work lives on the bookmark `arch-model-adrs` in the workspace `nyx-arch`.

## Code

`CLAUDE.md`, `docs/*_STATUS.md`, `.gitignore`.

## Contract

None. The rule is process, not a data shape.
