#!/usr/bin/env bash
# Reproducibly generates the 5 clean showcase figures in website/static/figures/
# via headless nyx scripting (`nyx run ... --workspace ...`).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXAMPLES_DIR="${REPO_ROOT}/examples"
OUT_DIR="${REPO_ROOT}/website/static/figures"

NYX_BIN="${NYX_BIN:-}"
if [[ -z "${NYX_BIN}" ]]; then
    if [[ -x "${REPO_ROOT}/target/release/nyx" ]]; then
        NYX_BIN="${REPO_ROOT}/target/release/nyx"
    elif command -v nyx >/dev/null 2>&1; then
        NYX_BIN="$(command -v nyx)"
    else
        echo "error: nyx executable not found; build with 'cargo build --release -p nyx' or set NYX_BIN." >&2
        exit 1
    fi
fi

mkdir -p "${OUT_DIR}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/nyx-site-figures.XXXXXX")"
cleanup() {
    rm -rf "${WORK_DIR}"
}
trap cleanup EXIT

for ex in \
    15_soi_wire.nyx.json \
    18_wgm_cylindrical_pml.nyx.json \
    11_anchor_loss.nyx.json \
    20_cavity_3d_spherical_pml.nyx.json
do
    cp "${EXAMPLES_DIR}/${ex}" "${WORK_DIR}/${ex}"
done

cat > "${WORK_DIR}/render_site_figures.nys" <<'EOF'
# Render the 5 showcase PNG plates for website/static/figures/

nyx.project_open(path="15_soi_wire.nyx.json")
nyx.solve_run()
print("15_soi_wire:", nyx.result_info(physics="ems"))
nyx.figure(physics="ems", kind="gallery", columns=2, max_modes=2, scale=2, path="soi_wire_modes.png")

nyx.project_open(path="18_wgm_cylindrical_pml.nyx.json")
nyx.solve_run()
print("18_wgm_cylindrical_pml:", nyx.result_info(physics="ems"))
nyx.figure(physics="ems", kind="gallery", columns=4, max_modes=4, scale=2, path="wgm_cyl_pml.png")

nyx.project_open(path="11_anchor_loss.nyx.json")
nyx.solve_run()
print("11_anchor_loss:", nyx.result_info(physics="sms"))
nyx.figure(physics="sms", kind="gallery", columns=3, max_modes=6, scale=2, path="anchor_loss_modes.png")

nyx.project_open(path="20_cavity_3d_spherical_pml.nyx.json")
nyx.solve_run()
print("20_cavity_3d_spherical_pml:", nyx.result_info(physics="ems"))
nyx.figure(physics="ems", kind="cutaway3d", scale=2, show_edges="regions", path="cavity_3d_cutaway.png")
nyx.figure(physics="ems", kind="slice3d", plane="z", scale=2, show_edges="regions", path="cavity_3d_slice.png")
EOF

"${NYX_BIN}" run "${WORK_DIR}/render_site_figures.nys" --workspace "${WORK_DIR}" --seconds 120

for plate in \
    soi_wire_modes.png \
    wgm_cyl_pml.png \
    anchor_loss_modes.png \
    cavity_3d_cutaway.png \
    cavity_3d_slice.png
do
    cp "${WORK_DIR}/${plate}" "${OUT_DIR}/${plate}"
done

echo "Generated 5 showcase plates in ${OUT_DIR}"
