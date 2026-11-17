#!/usr/bin/env bash
# Scripted screenshots of the window in both looks and two sizes, and a smoke test of the whole shell (ADR 88; docs/ux/UX_MERGED_PLAN.md, U0).
# Usage (from the workspace root):   bash tools/ux-shots.sh [output-dir]       default <target dir>/ux-shots
# It builds the release binary, then drives `nyx --script` for every view x light/dark x size, ALWAYS in a scratch NYX_WORKSPACE that is
# deleted afterwards: the first run of the UX review wrote into the real workspace (docs/UX_STATUS.md), so this script refuses to run without
# one. It fails when a run exits non-zero or a shot is missing or empty.
# UX_SHOTS_FAST=1 skips the views that need a solve (results, analysis). The shots show example data and the home path of the machine:
# do not commit them (docs/ux/shots holds the few that were reviewed).
set -eu
target=${CARGO_TARGET_DIR:-target}
out=${1:-$target/ux-shots}
mkdir -p "$out"
cargo build --release -p nyx
exe="$target/release/nyx"; [ -x "$exe.exe" ] && exe="$exe.exe"
[ -x "$exe" ] || { echo "no binary at $exe" >&2; exit 2; }
win() { if command -v cygpath >/dev/null 2>&1; then cygpath -m "$1"; else printf '%s' "$1"; fi; }
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
export NYX_WORKSPACE; NYX_WORKSPACE=$(win "$scratch")
[ -n "$NYX_WORKSPACE" ] || { echo "scratch workspace not set; refusing to touch the real one" >&2; exit 2; }
fail=0
for size in 1440x900 1100x700; do
  for look in light dark; do
    p="$(win "$out")/$look-$size"
    # views that need no solve: dashboard, settings, the command palette, the new-workflow sheet, a workflow's concept and build sections
    s="theme $look; page dashboard; shot $p-dashboard.png; page settings; shot $p-settings.png; page dashboard"
    s="$s; sheet palette; shot $p-palette.png; sheet new; shot $p-new.png; sheet none"
    s="$s; wf example 1; section concept; shot $p-concept.png; section build; shot $p-build.png"
    # views that show a result: wait for the solve (the script waits while the job runs)
    if [ -z "${UX_SHOTS_FAST:-}" ]; then
      s="$s; compute; shot $p-results.png; section analysis; shot $p-analysis.png"
    fi
    s="$s; quit"
    if ! "$exe" --size "$size" --script "$s"; then echo "FAIL: $look $size: nyx exited non-zero" >&2; fail=1; fi
  done
done
for f in "$out"/*.png; do
  [ -s "$f" ] || { echo "FAIL: empty or missing $f" >&2; fail=1; }
done
n=$(ls "$out"/*.png 2>/dev/null | wc -l)
[ -n "${UX_SHOTS_FAST:-}" ] && want=24 || want=32
[ "$n" -ge "$want" ] || { echo "FAIL: $n shots, expected $want" >&2; fail=1; }
echo "$n shots in $out"
exit $fail
