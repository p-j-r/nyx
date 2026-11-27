## Summary
<!-- What changed and why? Reference any related issue or ADR. -->

## Verification
- [ ] `bash docs/architecture/check_code_to_adr.sh` passes (new files mapped in `docs/architecture/CODE_TO_ADR.md`)
- [ ] `cargo fmt --check` passes on formatted crates (ADR 65)
- [ ] `cargo clippy --release --workspace --all-targets --no-deps -- -D warnings` passes
- [ ] Relevant `--release` tests pass (`tools/vvuq_ci.sh --tier 1` or `tools/vvuq_ci.ps1 -Tier 1`)
- [ ] No unintended changes to golden baselines (`**/tests/golden/**`, `crates/geom/golden/**`)

## Licensing
- [ ] I agree to the contribution licensing terms in `CONTRIBUTING.md`.
