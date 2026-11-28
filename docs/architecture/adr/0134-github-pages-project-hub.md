# 134. GitHub Pages project hub (`website/`), `site-prep` extraction, and confidentiality gate

Date: 2026-11-28

## Status

Accepted. Extends ADR 13, ADR 20, ADR 60, ADR 88, ADR 89, ADR 108, and ADR 133.

## Context

The repository already maintains machine-readable sources of truth for every public facet of `nyx`:
* Design colour roles and 8-colour CVD-safe region palettes in `docs/ux/tokens.json` (ADR 88).
* Multi-source constitutive material definitions in `materials/*.toml` (ADR 108, `format = 2`).
* Example simulation projects in `examples/*.nyx.json` (ADR 60).
* Versioned ASME V&V 20 validation datasets in `benchmarks/validation/*.json` and `docs/SQUADDS_BENCHMARK_REPORT.md` (ADR 133).
* Auto-generated engine command pages (`114` commands) in `docs/api/*.md` (ADR 89).

Publishing a public project hub on GitHub Pages requires turning these sources into a cohesive static site without introducing a Node.js/Python web-bundler toolchain (violating ADR 13), duplicating data by hand (causing drift), or leaking pre-publication/confidential substrate and material records into public build artifacts.

## Decision

* **Static site generator (Zola `0.23.6`):** Author the public site under `website/` using Zola (a single static Rust binary with built-in Tera templates, Sass compilation, and link checking) rather than JavaScript or Python static site frameworks.
* **Standalone Rust extraction CLI (`website/site-prep`):** Keep `website/site-prep` as an isolated `[workspace]` Rust package (excluded from the root workspace `Cargo.toml` and `Cargo.lock`) that parses `docs/ux/tokens.json`, `materials/*.toml`, `examples/*.nyx.json`, `benchmarks/validation/*.json`, `docs/SQUADDS_BENCHMARK_REPORT.md`, and `docs/api/*.md` into:
  * `website/sass/_tokens.scss` (`:root` / `[data-theme="light"]`, `[data-theme="dark"]`, and `@media (prefers-color-scheme: dark)` CSS custom properties for all ADR 88 roles, brand shadow tokens, 8-colour region palettes, and Swiss Lab Monograph font stacks).
  * `website/data/{materials,examples,benchmarks,api}.json`.
* **Swiss Lab Monograph aesthetic:** Pair the warm off-white page/card surfaces (`#F8F6F8` / `#FFFFFF`), dark mode (`#141218` / `#1D1A22`), and deep plum accent (`#820C4F` / `#D96AA4`) from `docs/ux/tokens.json` with editorial serif display typography (`Newsreader`), UI sans (`Inter`), monospace data tables (`IBM Plex Mono`), and `website/static/brand/nyx-mark.{png,svg}` derived from `docs/brand/nyx-mark.png`.
* **Fail-closed confidentiality gate (`site-prep --check`):** Exclude pre-publication dielectric/superconducting datasets, examples, and benchmark sub-cases from generated JSON payloads, sanitise cross-references in retained entries without modifying upstream repository sources, and scan all generated and authored `website/` files (excluding `website/public/` and `target/`) against the deny-list before `zola check` and `zola build`.
* **Automated Pages deployment (`.github/workflows/pages.yml`):** On pushes to `main` touching `website/**` or any upstream source of truth (plus `workflow_dispatch`), run `site-prep --repo . --out website --check`, `zola --root website check`, and `zola --root website build`, deploying `website/public` via `actions/upload-pages-artifact@v3` and `actions/deploy-pages@v4`.

## Consequences

* Root workspace builds, `Cargo.lock`, `crates/nyx-engine/tests/api_docs.rs`, and `docs/ux/tokens.json` remain completely untouched.
* Any drift in upstream `docs/api/*.md` command count or catalog hash, `docs/ux/tokens.json`, or any accidental confidential token anywhere in `website/` fails CI before static pages are built or deployed.
* `website/public/` and `website/site-prep/target/` are ignored in `.gitignore`.

## Code

`website/site-prep/*`, `website/sass/_tokens.scss`, `website/data/{materials,examples,benchmarks,api}.json`, `website/static/brand/nyx-mark.{png,svg}`, `website/{config.toml,content/*,templates/*,sass/*,static/*}`, `tools/gen_site_figures.*`, `.github/workflows/pages.yml`, `.gitignore`.

## Contract

`cargo run --manifest-path website/site-prep/Cargo.toml -- --repo . --out website --check` must exit `0` with zero confidentiality hits and verified `docs/api/` command count (`114`) and catalog hash before `zola --root website check` and `zola --root website build`.
