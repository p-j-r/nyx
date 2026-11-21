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
# Portability note & verification (ADR 107, QFR-7): nyx3 uses 100% pure-Rust `hdf5-pure` + `zlib-rs` (`crates/nyx-storage`),
# requiring zero external `hdf5.dll` / `zlib.dll` runtime DLLs or `HDF5_DIR` on PATH. Verify no dynamic HDF5 linkage slipped in.
checked_hdf5=0
for exe in nyx nyx-mcp nyx-worker nyx-script-worker nyx-connect; do
    bin="$out/$exe.exe"
    [ -f "$bin" ] || bin="$out/$exe"
    hit=""
    if command -v objdump >/dev/null 2>&1; then
        checked_hdf5=1
        hit="$({ objdump -p "$bin" 2>/dev/null || true; } | { grep -iE 'DLL Name:.*hdf5|NEEDED.*hdf5' || true; })"
    elif command -v ldd >/dev/null 2>&1; then
        checked_hdf5=1
        hit="$({ ldd "$bin" 2>/dev/null || true; } | { grep -iE 'hdf5\.dll|libhdf5' || true; })"
    fi
    if [ -n "$hit" ]; then
        echo "error: unexpected dynamic HDF5 linkage in $bin (violates ADR 107 pure-Rust hdf5-pure): $hit" >&2
        exit 1
    fi
done
if [ "$checked_hdf5" -eq 1 ]; then
    echo "verified zero dynamic hdf5.dll/libhdf5 dependencies (pure-Rust hdf5-pure + zlib-rs)"
fi
cp LICENSE NOTICE "$out/"
cp crates/app/assets/fonts/OFL-EBGaramond.txt "$out/"
rm -rf "$out/materials" "$out/examples" "$out/fonts"
cp -r materials examples "$out/"
cp -r crates/app/assets/fonts "$out/"
echo "packaged: $out"
