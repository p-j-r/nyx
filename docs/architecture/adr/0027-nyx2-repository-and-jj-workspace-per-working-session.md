# 27. The restructure happens in a second repository, `nyx2`, with its own history and one jj workspace per working session

Date: 2026-11-04

## Status

Accepted by the coordinator's instruction for the overnight restructure (2026-11-03). Refines ADR 20.

## Context

The restructure changes the crate layout, the mesh type, the linear algebra, the problem type and the file format. The old `nyx` repository (main `a047506`) must keep working and must not be touched. Several working sessions run in parallel and must not trample each other's working copies.

## Decision

* `nyx2` (`C:\Users\robin-adm\Documents\opomech\nyx2`) is a colocated jj repository cloned from the `arch-model-adrs` tip of `nyx` (`f74c12ea`), **keeping the full history**: the goldens, the ports' provenance, the commit hashes quoted in the ADRs and `git blame` stay valid. Its remote was removed: nothing is pushed anywhere. Its `main` is that tip; nobody commits to it directly.
* The tags `pre-restructure` and `sbs-last` both point at `f74c12ea` (code identical to `a047506`); `last-old-format` will be set in phase 8.
* **One jj workspace per working session**, next to the main one: `jj workspace add ../nyx2-<name> -r <base>`, with its own bookmark `restructure-<n>-<name>`. The base is the bookmark of the phase it builds on (`restructure-foundation` after this phase). Workspace folders are siblings of `nyx2` and share its history and `.cargo/config.toml` build directory.
* A finished phase bookmark is merged into the integration line by whoever coordinates, never by the phase's own session. Sessions do not push, do not move `main`, and update their own section of `docs/RESTRUCTURE_STATUS.md` before finishing.
* Baselines are produced by `tools/baseline.sh` (tests and the five benchmarks with fixed arguments) and recorded in the status file.

## Consequences

* The old repository is read-only for the restructure; the old `nyx-jj`, `nyx-arch` workspaces are unaffected.
* Disk: workspaces are cheap (sources only) because the target directory is shared (ADR 26).
* Conflicts between sessions are resolved by jj at rebase time; the layering table and `CODE_TO_ADR.md` are the expected one-line conflict spots.

## Code

`tools/baseline.sh`, `.cargo/config.toml`, `docs/RESTRUCTURE_STATUS.md`.

## Contract

None.
