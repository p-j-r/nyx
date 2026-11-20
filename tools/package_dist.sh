#!/usr/bin/env bash
# Builds the release programs of nyx and puts them in dist/<name>/ with the licence files (ADR 100): the program, the MCP server, the workers
# and nyx-connect, LICENSE, NOTICE and the licence text of the font. Usage: bash tools/package_dist.sh [name]  (default: nyx-main-<short commit>).
# Set CARGO_TARGET_DIR to a folder of your own when several jj workspaces build on one machine.
set -euo pipefail
cd "$(dirname "$0")/.."
name="${1:-nyx-main-$(jj --ignore-working-copy log -r 'main' --no-graph -T 'commit_id.short(8)' 2>/dev/null || echo local)}"
target="${CARGO_TARGET_DIR:-../nyx2-target}"
cargo build --release -p nyx -p nyx-engine -p nyx-mcp -p nyx-script -p nyx-connect
out="dist/$name"
mkdir -p "$out"
# the file names of the programs that make up a nyx install (the window finds the workers and the script worker next to itself)
for exe in nyx nyx-mcp nyx-worker nyx-script-worker nyx-connect; do
    cp "$target/release/$exe.exe" "$out/" 2>/dev/null || cp "$target/release/$exe" "$out/"
done
cp LICENSE NOTICE "$out/"
cp crates/app/assets/fonts/OFL-EBGaramond.txt "$out/"
echo "packaged: $out"
