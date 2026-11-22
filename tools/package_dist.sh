#!/usr/bin/env bash
# Builds the release programs of nyx and invokes `nyx-pack` to produce dist/<name>/ (split layout),
# dist/<name>.zip, and the self-extracting installer dist/<name>-setup[.exe] with 96B NYXINST1 trailer (ADR 122).
# Usage: bash tools/package_dist.sh [name]  (default: nyx-main-<short commit>).
# Set CARGO_TARGET_DIR to a folder of your own when several jj workspaces build on one machine.
set -euo pipefail
cd "$(dirname "$0")/.."
name="${1:-nyx-main-$(jj --ignore-working-copy log -r 'main' --no-graph -T 'commit_id.short(8)' 2>/dev/null || echo local)}"
target="${CARGO_TARGET_DIR:-../nyx2-target}"

cargo build --release -p nyx -p nyx-engine -p nyx-mcp -p nyx-script -p nyx-connect -p nyx-install
if [ "${NYX_INSTALLER_NO_GUI:-0}" != "1" ]; then
    cargo build --release -p nyx-install --features gui --bin nyx-setup
fi

exe_ext=""
[ -f "$target/release/nyx.exe" ] && exe_ext=".exe"

out="dist/$name"
rm -rf "$out"
"$target/release/nyx-pack$exe_ext" \
    --target-dir "$target/release" \
    --repo-root "." \
    --dir "$out" \
    --self-extract "dist/$name-setup$exe_ext" \
    --zip "dist/$name.zip"

# Portability verification (ADR 107, QFR-7): nyx3 uses pure-Rust `hdf5-pure` + `zlib-rs`,
# requiring zero external `hdf5.dll` / `libhdf5` runtime dependencies.
checked_hdf5=0
for exe in nyx nyx-mcp nyx-worker nyx-script-worker nyx-connect nyx-uninstall; do
    bin="$out/bin/$exe$exe_ext"
    [ -f "$bin" ] || continue
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
echo "packaged: $out (+ dist/$name-setup$exe_ext, dist/$name.zip)"
