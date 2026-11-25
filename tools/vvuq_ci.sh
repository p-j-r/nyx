#!/usr/bin/env bash
# VVUQ 3-tier verification runner (Linux / macOS / WSL local mirror of .github/workflows/vvuq-ci.yml).
# Usage:
#   bash tools/vvuq_ci.sh --tier 1 [--exclude-nyx]
#   bash tools/vvuq_ci.sh --tier 2
#   bash tools/vvuq_ci.sh --tier 3
set -euo pipefail

cd "$(dirname "$0")/.."

TIER=1
EXCLUDE_NYX=0

while [ $# -gt 0 ]; do
  case "$1" in
    --tier)
      TIER="${2:?--tier requires 1, 2, or 3}"
      shift 2
      ;;
    --exclude-nyx)
      EXCLUDE_NYX=1
      shift
      ;;
    -h|--help)
      sed -n '2,6p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

run_step() {
  local label="$1"
  shift
  echo "=== [VVUQ Tier ${TIER}] ${label} ==="
  "$@"
}

run_tier1() {
  run_step "ADR governance (check_code_to_adr.sh)" bash docs/architecture/check_code_to_adr.sh
  run_step "ADR 65 formatting (7 crates)" \
    cargo fmt --check -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine

  if [ "$EXCLUDE_NYX" -eq 1 ]; then
    run_step "Clippy (--release, --exclude nyx)" \
      cargo clippy --release --workspace --exclude nyx --all-targets --no-deps -- -D warnings
  else
    run_step "Clippy (--release, workspace)" \
      cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
  fi

  run_step "Build nyx-script worker & nyx-connect fake_mcp" bash -c '
    set -euo pipefail
    cargo build --release -p nyx-script
    cargo build --release -p nyx-connect --example fake_mcp
  '

  run_step "nyx-engine contracts, layering & security gates" bash -c '
    set -euo pipefail
    cargo test --release -p nyx-engine \
      --test layering --test contracts --test geomdesc --test api_docs \
      --test path_guard --test format
    case "$(uname -s)" in
      MINGW*|MSYS*|CYGWIN*)
        cargo test --release -p nyx-engine --test h5
        ;;
      *)
        cargo test --release -p nyx-engine --test h5 -- \
          --skip h5py_reads_export_results_l2_output_at_zero_ulp \
          --skip nyx_worker_has_zero_hdf5_dll_dependency_and_runs_without_anaconda_in_path
        ;;
    esac
  '

  run_step "femcore contracts" \
    cargo test --release -p femcore \
      --test contracts --test linalg_contract --test material_contract

  run_step "Workspace crate contracts (9 crates)" \
    cargo test --release \
      -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems \
      -p nyx-store -p nyx-mcp -p nyx-connect -p nyx-install \
      --test contracts

  run_step "nyx-script contract" cargo test --release -p nyx-script --test contract

  run_step "P0/P1 soundness, SSRF, proptest & finite_guard regressions" \
    bash -c '
      set -euo pipefail
      cargo test --release -p nyx-engine --lib material_fetch
      cargo test --release -p geom --test props --test gds
      cargo test --release -p nyx-sms --test finite_guard
      for crate in nyx-dfs nyx-ems; do
        if [ -f "crates/${crate}/tests/finite_guard.rs" ]; then
          cargo test --release -p "${crate}" --test finite_guard
        fi
      done
      for target in fuzz_gds fuzz_nyx_storage fuzz_script_trace fuzz_path_and_archive fuzz_mcp_framing; do
        hit=$(find crates -maxdepth 3 -name "${target}.rs" | head -n 1 || true)
        if [ -n "$hit" ]; then
          pkg=$(printf "%s" "$hit" | cut -d/ -f2)
          cargo test --release -p "$pkg" --test "$target"
        fi
      done
    '
}

run_tier2() {
  run_step "Build nyx-script worker" cargo build --release -p nyx-script
  run_step "Full workspace --release test suite (includes nyx-storage gates & examples)" \
    cargo test --release --workspace
  run_step "Differential portable-kernels (default features)" \
    cargo test --release -p femcore --test portable_kernels
  run_step "Differential portable-kernels (--features portable-kernels)" \
    cargo test --release -p femcore --features portable-kernels --test portable_kernels
  run_step "femcore linalg_contract (--features portable-kernels)" \
    cargo test --release -p femcore --features portable-kernels --test linalg_contract
}

run_tier3() {
  run_step "Build nyx-script worker" cargo build --release -p nyx-script
  run_step "Ignored deep verification tests" cargo test --release --workspace -- --ignored

  if rustup toolchain list | grep -q '^nightly'; then
    for mode in stacked tree; do
      if [ "$mode" = "stacked" ]; then
        mflags="-Zmiri-disable-isolation"
      else
        mflags="-Zmiri-tree-borrows -Zmiri-strict-provenance -Zmiri-disable-isolation"
      fi
      run_step "Miri (${mode} borrows, NYX_THREADS=2)" env NYX_THREADS=2 MIRIFLAGS="$mflags" bash -c '
        set -euo pipefail
        cargo +nightly miri test -p femcore --lib linalg::par::
        cargo +nightly miri test -p femcore --lib linalg::csr::
        cargo +nightly miri test -p femcore --lib linalg::basis::
        cargo +nightly miri test -p femcore --lib linalg::vecops::
        cargo +nightly miri test -p femcore --lib linalg::amg_elastic::
        cargo +nightly miri test -p geom --lib
        cargo +nightly miri test -p nyx-physics --lib
        cargo +nightly miri test -p nyx-store --lib
      '
    done
  else
    echo "ERROR: nightly toolchain + miri required for Tier 3 Miri checks" >&2
    exit 1
  fi

  if command -v cargo-kani >/dev/null 2>&1 || cargo kani --version >/dev/null 2>&1; then
    run_step "Kani K1 (femcore)" bash -c '
      set -euo pipefail
      cargo kani -p femcore --harness verify_csr_rows_into
      cargo kani -p femcore --harness verify_triplets_to_csr
      cargo kani -p femcore --harness verify_csr_sum_all
      cargo kani -p femcore --harness verify_raw_cols_disjoint_chunks
      cargo kani -p femcore --harness verify_raw_cols_rejects_overflow
    '
    run_step "Kani K2 (nyx-storage)" cargo kani -p nyx-storage --harness verify_read_budget_no_overflow
    run_step "Kani K3 (geom)" cargo kani -p geom --harness verify_gds_ibm_and_framing
  else
    echo "ERROR: cargo-kani required for Tier 3 Kani checks" >&2
    exit 1
  fi

  if command -v cargo-bolero >/dev/null 2>&1; then
    run_step "Bolero time-boxed libFuzzer + ASan campaigns (900s/target)" bash -c '
      set -euo pipefail
      for target in fuzz_gds fuzz_nyx_storage fuzz_script_trace fuzz_path_and_archive fuzz_mcp_framing; do
        hit=$(find crates -maxdepth 3 -name "${target}.rs" | head -n 1 || true)
        if [ -n "$hit" ]; then
          pkg=$(printf "%s" "$hit" | cut -d/ -f2)
          cargo +nightly bolero test -p "$pkg" --engine libfuzzer --sanitizer address --time 900 "$target"
        fi
      done
    '
  else
    echo "NOTE: cargo-bolero not installed; skipping time-boxed libFuzzer campaign"
  fi
}

case "$TIER" in
  1) run_tier1 ;;
  2) run_tier2 ;;
  3) run_tier3 ;;
  *)
    echo "Invalid tier '$TIER' (expected 1, 2, or 3)" >&2
    exit 2
    ;;
esac

echo "=== [VVUQ Tier ${TIER}] PASSED ==="
