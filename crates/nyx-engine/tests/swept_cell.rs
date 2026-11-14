//! A periodic cell made by sweeping a 2D section: the section is meshed, its vertices are renumbered by coordinates so that opposite walls
//! are cut alike, and layers with a `fill` region take that region everywhere (cladding above and below a patterned slab). Floquet pairs
//! then join the walls. Checked on a vacuum cell with a closed form, and by the refusal of a fill that names no region.

use nyx_engine::Engine;
use serde_json::{json, Value};

const C0: f64 = 299_792_458.0;

/// Section x in [0, lx], z in [0, lz], with a rectangular `Core` region inside a `Host` region (both vacuum), periodic along x and z.
fn project(fill: &str, lx: f64, lz: f64, ly: f64, h: f64) -> Value {
    let feats = json!([
        {"id": 1, "name": "cell", "op": {"Rect": {"min": [0.0, 0.0], "max": [lx, lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 2, "name": "core", "op": {"Rect": {"min": [0.3 * lx, 0.3 * lz], "max": [0.6 * lx, 0.7 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 3, "name": "overlay", "op": {"Overlay": {"naming": "Feature", "pins": [], "regions": [
            {"input": 1, "kind": "Material", "name": "Host", "pinned_id": null},
            {"input": 2, "kind": "Material", "name": "Core", "pinned_id": null}]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 4, "name": "pz", "op": {"Periodic": {"sel_a": {"Name": "cell.bottom"}, "sel_b": {"Name": "cell.top"}, "translation": [0.0, 0.0, lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 5, "name": "px", "op": {"Periodic": {"sel_a": {"Name": "cell.left"}, "sel_b": {"Name": "cell.right"}, "translation": [lx, 0.0, 0.0]}}, "origin": "Code", "suppressed": false, "tags": []}
    ]);
    json!({
        "active": "ems", "format": 1, "title": "swept periodic cell",
        "geometry": {"dim": 2, "features": feats, "format": 1, "params": []},
        "mesh": {"h": h, "sweep": [{"thickness": ly / 4.0, "n": 1, "fill": fill}, {"thickness": ly / 2.0, "n": 2}, {"thickness": ly / 4.0, "n": 1, "fill": fill}]},
        "physics": [{"physics": "ems", "format": 1, "order": 2,
            "bcs": [
                {"kind": "periodic", "on": "name(cell.bottom)", "params": {"dest": "name(cell.top)", "phase": 0}},
                {"kind": "periodic", "on": "name(cell.left)", "params": {"dest": "name(cell.right)", "phase": 0}}],
            "materials": {"Host": {"eps_r": 1.0}, "Core": {"eps_r": 1.0}},
            "study": {"kind": "cavity", "params": {"f_search": 5.0e9, "n_modes": 4}}}]
    })
}

#[test]
fn a_swept_periodic_cell_with_filled_cladding_layers_gives_the_vacuum_modes() {
    let (lx, lz, ly) = (0.04, 0.03, 0.025);
    let mut e = Engine::new();
    e.call("project_load", &json!({"json": project("Host", lx, lz, ly, 0.01).to_string()})).unwrap();
    e.call("mesh_build", &json!({})).unwrap();
    // the Core region exists only in the middle layer
    let model = e.model().unwrap().clone();
    let core = model.regions().iter().find(|r| r.1 == "Core").unwrap().0;
    let m = model.solver_mesh();
    let n_core = m.cell_labels.iter().filter(|&&l| l == core).count();
    assert!(n_core > 0, "the unfilled middle layer keeps the core");
    let (lo, hi) = (ly / 4.0, 0.75 * ly);
    for c in 0..m.n_cells() {
        if m.cell_labels[c] == core {
            let y: Vec<f64> = m.cell(c).iter().map(|&v| m.point(v as usize)[1]).collect();
            assert!(
                y.iter().all(|&y| y >= lo - 1e-12 && y <= hi + 1e-12),
                "a core cell outside the middle layer: {y:?}"
            );
        }
    }
    e.call("solve_run", &json!({"physics": "ems"})).unwrap();
    let (res, _) = e.result("ems").unwrap();
    let f: Vec<f64> = res.solutions.iter().filter_map(|s| s.value).collect();
    // Gamma point of the vacuum cell, PEC ends along y: two polarisations at c / (2 ly)
    let exact = C0 / (2.0 * ly);
    let near = f.iter().filter(|&&v| (v / exact - 1.0).abs() < 5e-3).count();
    assert!(near >= 2, "expected a degenerate pair at {exact:.4e} Hz, got {f:?}");
}

#[test]
fn a_fill_region_that_does_not_exist_is_refused_by_name() {
    let mut e = Engine::new();
    e.call("project_load", &json!({"json": project("Nope", 0.04, 0.03, 0.025, 0.01).to_string()})).unwrap();
    let err = e.call("mesh_build", &json!({})).unwrap_err().to_string();
    assert!(err.contains("Nope") && err.contains("Host"), "{err}");
}
