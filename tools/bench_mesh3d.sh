#!/usr/bin/env bash
# Best of three of the 3D meshing benchmarks (ADR 50 to 55). Usage: bash tools/bench_mesh3d.sh   (run from the workspace root, machine otherwise idle)
# For a control run, make a jj workspace at the parent commit, copy crates/femcore/examples/mesh3d_bench.rs into it (it uses nothing newer than the
# parent), run the same script there with its own CARGO_TARGET_DIR, and compare line by line. mesh3d_features_bench has no control: the parent refuses those features.
set -u
cargo build --release -p femcore --example mesh3d_bench --example mesh3d_features_bench > /dev/null 2>&1 \
  && cargo build --release -p nyx-engine --example cavity_bench > /dev/null 2>&1 || { echo "build failed"; exit 1; }
t=${CARGO_TARGET_DIR:-target}
for i in 1 2 3; do
  echo "== run $i mesh3d_bench"; "$t/release/examples/mesh3d_bench"
  echo "== run $i mesh3d_features_bench"; "$t/release/examples/mesh3d_features_bench" 2>/dev/null
  s=$(date +%s%N); "$t/release/examples/cavity_bench" 0.012 2 > /dev/null 2>&1; e=$(date +%s%N)
  echo "== run $i cavity_bench wall $(( (e - s) / 1000000 )) ms"
done
