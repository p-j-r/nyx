#!/usr/bin/env bash
# Records the test counts and benchmark timings that every restructure phase is compared against (plan section 3, phase 0).
# Usage:  bash tools/baseline.sh <output-dir> [examples-package]   (run from the workspace root; always --release)
#         BENCH_ONLY=1 bash tools/baseline.sh <output-dir>        skips the test run (timings only; keep the machine otherwise idle)
# Writes <output-dir>/tests.log and one log per benchmark, and prints a summary line per step.
# Benchmarks are the five named in the plan; the arguments are fixed here so a later phase can rerun exactly the same.
set -u
out=${1:?output dir}; pkg=${2:-femcore}
mkdir -p "$out"
step() { # name, command...
  local name=$1; shift
  local t0=$(date +%s%N)
  "$@" > "$out/$name.log" 2>&1; local rc=$?
  local t1=$(date +%s%N)
  local ms=$(( (t1 - t0) / 1000000 )); printf "%-14s rc=%d  %d.%03d s wall (includes any compile)
" "$name" "$rc" $((ms/1000)) $((ms%1000))
}
if [ -z "${BENCH_ONLY:-}" ]; then
  step tests cargo test --release --workspace --no-fail-fast
  grep -E '^test result|FAILED|failed' "$out/tests.log" | sort | uniq -c | sort -rn | head -40 > "$out/tests.summary"
fi
# build the examples first so the timings below are run time only
cargo build --release -p "$pkg" --examples > "$out/build_examples.log" 2>&1
# `bench`, `profile3d`, `mem_profile` and `eig_beam` drive the physics through the window bridge, so they live in nyx-engine since restructure phase 4
cargo build --release -p nyx-engine --examples >> "$out/build_examples.log" 2>&1
export NYX_PROFILE=1
step bench         cargo run --release -p nyx-engine --example bench
step profile3d     cargo run --release -p nyx-engine --example profile3d -- 16 16 16 p2
step mem_profile   cargo run --release -p nyx-engine --example mem_profile -- 12 12 12 p2 static auto
step cavity_bench  cargo run --release -p nyx-engine --example cavity_bench -- 0.012 2
step solve_scaling cargo run --release -p "$pkg" --example solve_scaling
