#!/usr/bin/env bash
# Scripted screenshots of the window in both looks and two sizes, and a smoke test of the whole shell (ADR 88; docs/ux/UX_MERGED_PLAN.md, U0).
# Usage (from the workspace root):   bash tools/ux-shots.sh [output-dir]       default <target dir>/ux-shots
# It builds the release binary, then drives `nyx --script` for every view x light/dark x size, ALWAYS in a scratch NYX_WORKSPACE that is
# deleted afterwards: the first run of the UX review wrote into the real workspace (docs/UX_STATUS.md), so this script refuses to run without
# one. It fails when a run exits non-zero or a shot is missing or empty.
# The shots of the build view include the other canvas, an isolated legend row and colouring by region (plan U4), and a job that waits with the queue popover open (plan U5), and a sweep of one parameter with its running points, queue, table and chart (plan U5).
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
    s="$s; sheet palette; shot $p-palette.png; sheet new; shot $p-new.png; sheet about; shot $p-about.png; sheet none"
    s="$s; contrast high; page settings; shot $p-settings-high.png; contrast standard; page dashboard"
    s="$s; wf example 1; section concept; shot $p-concept.png; section build; shot $p-build.png"
    # the region colours: the other canvas, one legend row isolated, and colouring by region
    s="$s; canvas theme; shot $p-build-canvas-theme.png; isolate 1; shot $p-build-isolate.png; isolate none; colour region; shot $p-build-by-region.png; colour material; canvas dark"
    # views that show a result: wait for the solve (the script waits while the job runs)
    if [ -z "${UX_SHOTS_FAST:-}" ]; then
      # a job that waits (the queue is paused) with the queue popover open; resumed, the job runs and the script waits for it again
      s="$s; wait off; queue pause; compute; shot $p-queued.png; queue open; shot $p-queue.png; queue close; queue resume; wait on"
      s="$s; compute; shot $p-results.png; canvas theme; shot $p-results-canvas-theme.png; canvas dark; section analysis; shot $p-analysis.png"
      # P6: the trust state of a result on the field page: exploratory, then checked by a check that passes, then exploratory again by one that fails (the cantilever's numbers move with the mesh)
      s="$s; wf example 1; compute; section analysis; view field; shot $p-validated-exploratory.png; check 0.5; shot $p-validated-checked.png; check 1e-9; shot $p-validated-failed.png"
      # a workflow that was saved and is marked done moves to the Earlier table of the Dashboard
      s="$s; done 0; page dashboard; shot $p-earlier.png"
      # a sweep of one parameter: the points running and waiting (the chip and the queue popover), the finished table and chart, and the runs on the Analysis page with ticked runs ringed
      s="$s; wf example 14; section build; wait off; sweep run w 4e-7 6e-7 24; shot $p-sweep-running.png; queue open; shot $p-sweep-queue.png; queue close; wait on; shot $p-sweep-done.png; section analysis; tick 1 12 24; shot $p-sweep-analysis.png"
    fi
    s="$s; quit"
    if ! "$exe" --size "$size" --script "$s"; then echo "FAIL: $look $size: nyx exited non-zero" >&2; fail=1; fi
  done
done
# U6: the interface at 150 and 200 percent (the window at its 900 x 600 minimum in logical pixels: nothing may be cut or overlap), and the canvas as
# a person with a colour-vision difference sees it (cvd.rs: protan, deutan, tritan; a simulation, not a person). Light look, scratch workspace.
for sc in 1.5 2.0; do
  case $sc in 1.5) size=1350x900;; *) size=1800x1200;; esac
  p="$(win "$out")/scale-$sc"
  s="theme light; scale $sc; page dashboard; shot $p-dashboard.png; sheet new; shot $p-new.png; sheet none; wf example 14; section build; shot $p-build.png; quit"
  if ! "$exe" --size "$size" --script "$s"; then echo "FAIL: scale $sc: nyx exited non-zero" >&2; fail=1; fi
done
for k in protan deutan tritan; do
  p="$(win "$out")/cvd-$k"
  s="theme dark; scale 1.0; wf example 14; section build; cvd $k; shot $p-build.png"
  [ -z "${UX_SHOTS_FAST:-}" ] && s="$s; compute; shot $p-results.png"
  if ! "$exe" --size 1440x900 --script "$s; quit"; then echo "FAIL: cvd $k: nyx exited non-zero" >&2; fail=1; fi
done
for f in "$out"/*.png; do
  [ -s "$f" ] || { echo "FAIL: empty or missing $f" >&2; fail=1; }
done
n=$(ls "$out"/*.png 2>/dev/null | wc -l)
[ -n "${UX_SHOTS_FAST:-}" ] && want=53 || want=108
[ "$n" -ge "$want" ] || { echo "FAIL: $n shots, expected $want" >&2; fail=1; }
echo "$n shots in $out"
exit $fail
