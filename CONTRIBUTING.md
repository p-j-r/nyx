# Contributing to nyx

Thank you for your interest in contributing to `nyx`.

## Before Opening a Pull Request

1. **Read the conventions:** See `README.md`, `AGENTS.md` / `CLAUDE.md`, and `docs/architecture/ADR_INDEX.md` for crate layering and architecture decisions.
2. **Keep architecture mappings up to date:** Every new tracked file must match a row in `docs/architecture/CODE_TO_ADR.md` (`bash docs/architecture/check_code_to_adr.sh`).
3. **Do not silently edit golden test baselines** (`**/tests/golden/**`, `crates/geom/golden/**`) without explaining the numerical cause and verification in the PR description.
4. **Always test in `--release` mode** — debug builds of the FEM solvers are too slow for the test suite.

## Local Checks (Tier 1 CI Mirror)

Run the fast gate locally before pushing:

```bash
# Linux / macOS / WSL
bash tools/vvuq_ci.sh --tier 1

# Windows (PowerShell)
powershell -NoProfile -ExecutionPolicy Bypass -File .\tools\vvuq_ci.ps1 -Tier 1
```

Or run the core checks individually:

```bash
bash docs/architecture/check_code_to_adr.sh
cargo fmt --check -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine
cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
cargo test --release --workspace
```

> Note: Build `nyx-script` (`cargo build --release -p nyx-script`) and `nyx-worker` (`cargo build --release -p nyx-engine --bin nyx-worker`) before running integration suites that spawn workers.

## Pull Request Process

- Target `main` via a feature branch and pull request.
- Keep PRs focused and include tests for bug fixes and new physics/geometry features.
- External PRs require maintainer approval and passing `Tier 1: PR Fast Gate` checks on Linux and Windows.

## Licensing of Contributions

`nyx` is dual-licensed under **AGPL-3.0-only OR LicenseRef-nyx-Commercial** (see `LICENSE` and `NOTICE`), while `crates/geom` and ported solver/mesher modules in `crates/femcore` remain under BSD-2-Clause.

By submitting a pull request, you confirm that your contribution is your original work (or that you have the right to submit it) and agree that it may be distributed under the repository's dual licence (AGPL-3.0-only and commercial licences offered by the copyright holder), or BSD-2-Clause for files under BSD-2-Clause modules.
