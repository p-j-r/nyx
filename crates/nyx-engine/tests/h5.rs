// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The HDF5 results file (feature `hdf5`; ADR 62, ADR 102): what is written is read back bit for bit (0 ULP), Layout L2 (`Shuffle + LZF`
//! packed 5D) is smaller than both JSON and legacy L1 (`gzip-3`), vanilla `h5py` reads L2 natively at 0 ULP, `nyx-worker` has zero
//! `hdf5.dll` runtime dependency, and the engine commands write and describe the file.
#![cfg(feature = "hdf5")]

use hdf5_pure::{AttrValue, File, FileBuilder};
use nyx_engine::h5io::{FieldStorageLayout, ReadBudgetConfig};
use nyx_engine::Engine;
use nyx_storage::assert_f64_bits_eq;
use serde_json::json;
use std::path::{Path, PathBuf};

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_str()
        .unwrap()
        .to_string()
}

fn tmp(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("nyx-h5-{}-{name}", std::process::id()))
}

fn attr_string(attrs: &std::collections::HashMap<String, AttrValue>, key: &str) -> String {
    let v = &attrs[key];
    v.as_str()
        .map(str::to_string)
        .or_else(|| v.as_strings().and_then(|a| a.first().cloned()))
        .expect("string attribute")
}

#[test]
fn a_results_file_reads_back_exactly_and_is_smaller_than_json() {
    for name in ["15_soi_wire", "04_thermal_bar"] {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example(name) })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let path = tmp(&format!("{name}.h5"));
        let info = e.call("export_results", &json!({ "path": path.to_str().unwrap() })).unwrap();
        let bytes = info["bytes"].as_u64().unwrap();
        let stored = nyx_engine::h5io::read_results(&path).unwrap();
        assert_eq!(stored.project, e.project().to_json());
        let model = e.model().unwrap();
        let mesh = model.solver_mesh();
        assert_eq!(stored.mesh.points, mesh.points);
        assert_eq!(stored.mesh.cells, mesh.cells);
        assert_eq!(stored.mesh.cell_labels, mesh.cell_labels);
        let mut json_bytes = 0;
        for (code, r) in &stored.results {
            let (orig, _) = e.result(code).unwrap();
            assert_eq!(r, orig.as_ref(), "{name}/{code} reads back exactly");
            json_bytes += serde_json::to_vec(orig.as_ref()).unwrap().len();
        }
        assert!(bytes < json_bytes as u64, "{name}: {bytes} bytes of HDF5 against {json_bytes} of JSON");
        let d = e.call("results_file_info", &json!({ "path": path.to_str().unwrap() })).unwrap();
        assert!(!d["results"].as_array().unwrap().is_empty() && d["cells"].as_u64().unwrap() > 0, "{d}");
        assert_eq!(d["cells"].as_u64().unwrap() as usize, mesh.n_cells());
        assert_eq!(d["points"].as_u64().unwrap() as usize, mesh.n_points());
        eprintln!("{name}: HDF5 {bytes} bytes, JSON {json_bytes} bytes");
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(path.with_extension("xdmf"));
    }
}

#[test]
fn a_file_that_nyx_did_not_write_is_refused_and_nothing_is_exported_before_a_solve() {
    let mut e = Engine::new();
    assert!(e.call("export_results", &json!({ "path": tmp("none.h5").to_str().unwrap() })).is_err());
    let other = tmp("other.h5");
    FileBuilder::new().write(&other).unwrap();
    assert!(nyx_engine::h5io::read_results(&other).err().expect("refused").contains("not written by nyx"));
    assert!(nyx_engine::h5io::read_results_summary(&other).expect_err("refused").contains("not written by nyx"));
    let _ = std::fs::remove_file(&other);
    e.call("project_new", &json!({})).unwrap();
}

#[test]
fn every_recorded_run_keeps_its_results_in_an_h5_file() {
    let ws = tmp("ws");
    let _ = std::fs::remove_dir_all(&ws);
    let mut e = Engine::new();
    e.set_workspace(ws.clone(), true);
    e.call("project_open", &json!({ "path": example("06_wr90_waveguide") })).unwrap();
    let r = e.call("solve_run", &json!({})).unwrap();
    assert!(r["run"]["id"].is_string(), "{r}");
    let runs = e.call("store_runs", &json!({})).unwrap();
    let artifacts = runs[0]["artifacts"].as_array().unwrap();
    assert!(artifacts.iter().any(|a| a["path"].as_str().unwrap().ends_with("results.h5")), "{runs}");
    assert!(
        !artifacts.iter().any(|a| a["path"].as_str().unwrap().ends_with(".json")),
        "no JSON copy of the fields: {runs}"
    );
    let _ = std::fs::remove_dir_all(&ws);
}

#[test]
fn the_xdmf_file_names_only_datasets_that_exist_in_the_h5_file() {
    for name in ["04_thermal_bar", "02_cantilever", "03_heated_plate", "15_soi_wire"] {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example(name) })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let path = tmp(&format!("{name}-x.h5"));
        let info = e.call("export_results", &json!({ "path": path.to_str().unwrap() })).unwrap();
        let xdmf = std::fs::read_to_string(info["xdmf"].as_str().unwrap()).unwrap();
        assert!(xdmf.contains("<Xdmf") && xdmf.contains("<Attribute"), "{name}");
        // the file says it is for nyx only and states the cell node order its connectivity is in
        assert!(xdmf.contains("ParaView") && xdmf.contains("CellNodeOrder"), "{name}");
        let file = File::open_streaming(&path).unwrap();
        // no second copy of the nodes and cells: the XDMF names the ones of the results
        assert!(file.group("viz").is_err(), "{name}: there is no /viz group");
        assert!(xdmf.contains("/points") && xdmf.contains("/cells") && !xdmf.contains("/viz/"), "{name}");
        let fname = path.file_name().unwrap().to_str().unwrap();
        let mut refs = 0;
        for part in xdmf.split("<DataItem").skip(1) {
            let body = part.split('>').nth(1).unwrap().split('<').next().unwrap().trim();
            let (f, p) = body.split_once(':').expect("file:path");
            assert_eq!(f, fname);
            let ds = file
                .dataset(p.trim_start_matches('/'))
                .unwrap_or_else(|err| panic!("{name}: {p} does not exist: {err}"));
            let elems: u64 = ds.shape().unwrap().iter().product();
            assert!(elems > 0, "{name}: {p} is empty");
            refs += 1;
        }
        assert!(refs > 6, "{name}: {refs} references");
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(path.with_extension("xdmf"));
    }
}

#[test]
fn a_results_file_opens_as_the_project_with_its_results_and_analyses_work_on_them() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("15_soi_wire") })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let path = tmp("open.h5");
    e.call("export_results", &json!({ "path": path.to_str().unwrap() })).unwrap();
    let want_project = e.project().clone();
    let active = want_project.active.clone().unwrap_or_else(|| want_project.physics[0].physics.clone());
    let (solved, _) = e.result(&active).unwrap();
    let solved = solved.as_ref().clone();

    // a new engine (another program run) opens the file: the same project, the stored mesh, the stored results, none of them stale
    let mut f = Engine::new();
    let v = f.call("open_results", &json!({ "path": path.to_str().unwrap() })).unwrap();
    assert!(v["cells"].as_u64().unwrap() > 0 && !v["results"].as_array().unwrap().is_empty(), "{v}");
    assert_eq!(f.project(), &want_project);
    let (opened, stale) = f.result(&active).unwrap();
    assert!(!stale, "a result of the opened file is current");
    assert_eq!(opened.as_ref(), &solved, "the stored result is exactly the solved one");
    assert!(f.model().is_some());
    // what works on a solve works on the opened file
    assert!(f.call("result_info", &json!({})).is_ok());
    let _ = std::fs::remove_file(&path);
    let _ = std::fs::remove_file(path.with_extension("xdmf"));
}

#[test]
fn a_file_that_cannot_be_opened_leaves_the_project_as_it_was() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("04_thermal_bar") })).unwrap();
    let before = e.project().clone();
    let other = tmp("notnyx.h5");
    FileBuilder::new().write(&other).unwrap();
    assert!(e.call("open_results", &json!({ "path": other.to_str().unwrap() })).is_err());
    assert!(e.call("open_results", &json!({ "path": tmp("missing.h5").to_str().unwrap() })).is_err());
    assert_eq!(e.project(), &before);
    let _ = std::fs::remove_file(&other);
}

#[test]
fn curved_edges_are_stored_in_the_file_and_a_file_without_them_reads_as_straight() {
    // ADR 82: two optional datasets in /mesh; a straight mesh writes none
    let mut opts = femcore::meshing::geomdrive::GeomMeshOpts::new(0.5);
    opts.curved = true;
    let curved = femcore::meshing::geomdrive3::sphere_cavity_with(1.0, &opts).unwrap().mesh;
    assert!(curved.has_curved());
    let straight = femcore::meshing::geomdrive3::sphere_cavity(1.0, 0.5).unwrap().mesh;
    assert!(!straight.has_curved());
    let project = nyx_engine::Project::new("curved edges", 3);
    for (name, mesh) in [("curved", &curved), ("straight", &straight)] {
        let path = tmp(&format!("edge-nodes-{name}.h5"));
        nyx_engine::h5io::write_results(&path, &project, mesh, &[]).unwrap();
        let stored = nyx_engine::h5io::read_results(&path).unwrap();
        assert_eq!(stored.mesh.edge_node_list(), mesh.edge_node_list(), "{name}");
        assert_eq!(stored.mesh.has_curved(), name == "curved");
        let file = File::open_streaming(&path).unwrap();
        let has = file.group("mesh").unwrap().dataset("edge_node_pairs").is_ok();
        assert_eq!(has, name == "curved", "the datasets exist only for a curved mesh");
        let _ = std::fs::remove_file(&path);
    }
}

/// The documented cell node order of layout 4 is the real one: in every cell the edge nodes sit at the midpoints of the edges the
/// `cell_order` text names, and cell fields have one value per cell node in both L2 (`cell_fields_packed`) and L1 (`solution_<k>/cell_fields`).
#[test]
fn the_cell_node_order_of_a_results_file_is_the_documented_one() {
    for name in ["15_soi_wire", "04_thermal_bar"] {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example(name) })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let mesh = e.model().unwrap().solver_mesh();
        let results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)> = e
            .project()
            .physics
            .iter()
            .filter_map(|s| e.result(&s.physics).map(|(r, _)| (s.physics.clone(), r.clone())))
            .collect();
        let refs: Vec<(String, &nyx_physics::SolveResult)> =
            results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();

        let path_l2 = tmp(&format!("{name}-order-l2.h5"));
        e.call("export_results", &json!({ "path": path_l2.to_str().unwrap() })).unwrap();
        let path_l1 = tmp(&format!("{name}-order-l1.h5"));
        nyx_engine::h5io::write_results_with_layout(
            &path_l1,
            e.project(),
            &mesh,
            &refs,
            FieldStorageLayout::L1LegacyGzip3,
        )
        .unwrap();

        for (tag, path) in [("L2", &path_l2), ("L1", &path_l1)] {
            let file = File::open_streaming(path).unwrap();
            assert_eq!(file.root().attrs().unwrap()["version"], AttrValue::I32(4), "{name} {tag}: version is 4i32");
            for code in file.group("results").unwrap().groups().unwrap() {
                let g = file.group(&format!("results/{code}")).unwrap();
                let g_attrs = g.attrs().unwrap();
                let order = attr_string(&g_attrs, "cell_order");
                let (pts, cells) = (g.dataset("points").unwrap(), g.dataset("cells").unwrap());
                let dim = pts.shape().unwrap()[1] as usize;
                let (p, c) = (pts.read_f64().unwrap(), cells.read_u32().unwrap());
                let nloc = cells.shape().unwrap()[1] as usize;
                assert_eq!(order.split(' ').count(), nloc, "{name}/{code}: {order}");
                // the edge nodes: `eij` is the midpoint of vertices i and j
                let words: Vec<&str> = order.split(' ').collect();
                let at = |n: u32, k: usize| p[n as usize * dim + k];
                for row in c.chunks(nloc) {
                    for (a, w) in words.iter().enumerate() {
                        if let Some(rest) = w.strip_prefix('e') {
                            let (i, j) = (rest[0..1].parse::<usize>().unwrap(), rest[1..2].parse::<usize>().unwrap());
                            for k in 0..dim {
                                let mid = 0.5 * (at(row[i], k) + at(row[j], k));
                                assert!(
                                    (at(row[a], k) - mid).abs() < 1e-9 * (1.0 + mid.abs()),
                                    "{name}/{code}: node {w} is not the midpoint of its edge"
                                );
                            }
                        }
                    }
                }
                if tag == "L2" {
                    if let Ok(packed) = g.dataset("cell_fields_packed") {
                        let sh = packed.shape().unwrap();
                        assert_eq!(sh.len(), 5);
                        assert_eq!(
                            (sh[2] * sh[3]) as usize,
                            c.len(),
                            "{name}/{code} L2: n_cells * nloc == cells.len()"
                        );
                        assert_eq!(sh[4], 1, "{name}/{code} L2: real-only solve has R=1");
                    }
                } else {
                    let mut checked = 0usize;
                    for s in g.groups().unwrap().iter().filter(|n| n.starts_with("solution_")) {
                        let Ok(cf) = g.group(&format!("{s}/cell_fields")) else { continue };
                        for f in cf.groups().unwrap() {
                            let re = cf.dataset(&format!("{f}/re")).unwrap().read_f64().unwrap();
                            assert_eq!(re.len(), c.len(), "{name}/{code}/{s}/{f}: one value per cell node");
                            checked += 1;
                        }
                    }
                    assert!(
                        checked > 0
                            || refs
                                .iter()
                                .find(|(k, _)| k == &code)
                                .unwrap()
                                .1
                                .solutions
                                .iter()
                                .all(|s| s.cell_fields.is_empty())
                    );
                }
            }
        }
        let _ = std::fs::remove_file(&path_l2);
        let _ = std::fs::remove_file(path_l2.with_extension("xdmf"));
        let _ = std::fs::remove_file(&path_l1);
    }
}

#[test]
fn layouts_l2_l3_and_l1_round_trip_at_zero_ulp_for_real_and_complex_solvers() {
    // 15_soi_wire: real-only lossless 2D waveguide (im: vec![])
    // 04_thermal_bar: coupled dfs + sms stationary (im: vec![])
    // 09_pml_slab: complex PML leaky waveguide modes (non-empty im!)
    for (name, expect_complex_im) in [("15_soi_wire", false), ("04_thermal_bar", false), ("09_pml_slab", true)] {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example(name) })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let mesh = e.model().unwrap().solver_mesh();
        let results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)> = e
            .project()
            .physics
            .iter()
            .filter_map(|s| e.result(&s.physics).map(|(r, _)| (s.physics.clone(), r.clone())))
            .collect();
        let refs: Vec<(String, &nyx_physics::SolveResult)> =
            results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();

        for layout in [
            FieldStorageLayout::L2ShuffleLzfPacked,
            FieldStorageLayout::L3ContiguousPacked,
            FieldStorageLayout::L1LegacyGzip3,
        ] {
            let path = tmp(&format!("rt-{name}-{}.h5", layout.as_str()));
            nyx_engine::h5io::write_results_with_layout(&path, e.project(), &mesh, &refs, layout).unwrap();
            let stored = nyx_engine::h5io::read_results(&path).unwrap();
            assert_eq!(stored.project, e.project().to_json());
            assert_f64_bits_eq(&format!("{name}/{}/mesh.points", layout.as_str()), &stored.mesh.points, &mesh.points)
                .unwrap();
            assert_eq!(stored.mesh.cells, mesh.cells);
            assert_eq!(stored.mesh.cell_labels, mesh.cell_labels);
            assert_eq!(stored.mesh.facets, mesh.facets);
            assert_eq!(stored.mesh.facet_labels, mesh.facet_labels);

            let mut saw_cf = false;
            for (code, r) in &stored.results {
                let (orig, _) = e.result(code).unwrap();
                assert_eq!(r, orig.as_ref(), "{name}/{code}/{:?} SolveResult equality", layout);
                for (k, (sa, sb)) in r.solutions.iter().zip(&orig.solutions).enumerate() {
                    for (fa, fb) in sa.fields.iter().zip(&sb.fields) {
                        assert_f64_bits_eq(&format!("{name}/{code}/sol_{k}/field_{}", fa.name), &fa.values, &fb.values)
                            .unwrap();
                    }
                    for (ca, cb) in sa.cell_fields.iter().zip(&sb.cell_fields) {
                        saw_cf = true;
                        assert_eq!(
                            ca.im.is_empty(),
                            expect_complex_im.not_expected_empty(),
                            "{name}/{code}/{:?} im emptiness",
                            layout
                        );
                        assert_f64_bits_eq(&format!("{name}/{code}/sol_{k}/cf_{}/re", ca.name), &ca.re, &cb.re)
                            .unwrap();
                        assert_f64_bits_eq(&format!("{name}/{code}/sol_{k}/cf_{}/im", ca.name), &ca.im, &cb.im)
                            .unwrap();
                    }
                }
            }
            assert!(saw_cf, "{name} produced cell fields");

            if layout != FieldStorageLayout::L1LegacyGzip3 {
                let raw_file = File::open_streaming(&path).unwrap();
                for (code, _) in &stored.results {
                    let g = raw_file.group(&format!("results/{code}")).unwrap();
                    let sh = g.dataset("cell_fields_packed").unwrap().shape().unwrap();
                    assert_eq!(sh.len(), 5);
                    assert_eq!(sh[4], if expect_complex_im { 2 } else { 1 }, "{name}/{code}/{:?} R axis", layout);
                }
            }

            // ReadBudgetConfig enforcement check: a 64-byte budget rejects whole-dataset reads cleanly
            let tiny_budget =
                ReadBudgetConfig { max_dataset_bytes: 64, max_group_depth: 24, ..ReadBudgetConfig::default() };
            let err = nyx_engine::h5io::read_results_with_budget(&path, tiny_budget).err().expect("budget enforced");
            assert!(err.contains("exceeds read budget"), "{err}");

            let _ = std::fs::remove_file(&path);
        }
    }
}

trait BoolHelper {
    fn not_expected_empty(self) -> bool;
}
impl BoolHelper for bool {
    fn not_expected_empty(self) -> bool {
        !self
    }
}

#[test]
fn l2_shuffle_lzf_packed_is_strictly_smaller_than_l1_legacy_gzip3_on_soi_wire() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("15_soi_wire") })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let mesh = e.model().unwrap().solver_mesh();
    let results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)> = e
        .project()
        .physics
        .iter()
        .filter_map(|s| e.result(&s.physics).map(|(r, _)| (s.physics.clone(), r.clone())))
        .collect();
    let refs: Vec<(String, &nyx_physics::SolveResult)> = results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();

    let (p_l1, p_l2, p_l3) = (tmp("size-soi-l1.h5"), tmp("size-soi-l2.h5"), tmp("size-soi-l3.h5"));
    nyx_engine::h5io::write_results_with_layout(&p_l1, e.project(), &mesh, &refs, FieldStorageLayout::L1LegacyGzip3)
        .unwrap();
    nyx_engine::h5io::write_results_with_layout(
        &p_l2,
        e.project(),
        &mesh,
        &refs,
        FieldStorageLayout::L2ShuffleLzfPacked,
    )
    .unwrap();
    nyx_engine::h5io::write_results_with_layout(
        &p_l3,
        e.project(),
        &mesh,
        &refs,
        FieldStorageLayout::L3ContiguousPacked,
    )
    .unwrap();

    let (b_l1, b_l2, b_l3) =
        (nyx_engine::h5io::file_bytes(&p_l1), nyx_engine::h5io::file_bytes(&p_l2), nyx_engine::h5io::file_bytes(&p_l3));
    eprintln!("15_soi_wire file sizes: L1={b_l1} B, L2={b_l2} B, L3={b_l3} B");
    assert!(b_l2 < b_l1, "L2 ({b_l2} B) must be strictly smaller than L1 ({b_l1} B) on 15_soi_wire");
    assert!(b_l2 < b_l3, "L2 ({b_l2} B) must be strictly smaller than uncompressed L3 ({b_l3} B)");

    for p in [p_l1, p_l2, p_l3] {
        let _ = std::fs::remove_file(p);
    }
}

#[test]
fn h5py_reads_export_results_l2_output_at_zero_ulp() {
    let py = Path::new(r"C:\ProgramData\anaconda3\python.exe");
    assert!(py.exists(), "expected Anaconda python with h5py at {}", py.display());

    let script = r#"
import sys, json, os, h5py, numpy as np

h5_path, out_dir = sys.argv[1], sys.argv[2]
with h5py.File(h5_path, "r") as f:
    fmt = f.attrs["format"]
    if isinstance(fmt, bytes):
        fmt = fmt.decode("utf-8")
    assert fmt == "nyx-results", f"unexpected format: {fmt}"
    assert int(f.attrs["version"]) == 4, f"unexpected version: {f.attrs['version']}"
    assert f["mesh/points"].dtype == np.float64

    np.asarray(f["mesh/points"][:], dtype="<f8").tofile(os.path.join(out_dir, "mesh_points.bin"))
    g = f["results/ems"]
    fl = g.attrs["field_layout"]
    if isinstance(fl, bytes):
        fl = fl.decode("utf-8")
    assert fl == "L2ShuffleLzfPacked", f"unexpected field_layout: {fl}"

    np.asarray(g["points"][:], dtype="<f8").tofile(os.path.join(out_dir, "ems_points.bin"))
    np.asarray(g["eigenvalues"][:], dtype="<f8").tofile(os.path.join(out_dir, "ems_eigenvalues.bin"))
    packed = np.asarray(g["cell_fields_packed"][:], dtype="<f8")
    packed.tofile(os.path.join(out_dir, "cell_fields_packed.bin"))
    meta = json.loads(bytes(g["cell_fields_meta"][:]).decode("utf-8"))
    with open(os.path.join(out_dir, "meta.json"), "w", encoding="utf-8") as mf:
        json.dump({"packed_shape": list(packed.shape), "meta": meta}, mf)

    f0_names = sorted(g["solution_0/fields"].keys())
    np.asarray(g["solution_0/fields"][f0_names[0]][:], dtype="<f8").tofile(os.path.join(out_dir, "sol0_f0.bin"))
"#;

    for (name, expected_ri) in [("15_soi_wire", 1usize), ("09_pml_slab", 2usize)] {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example(name) })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let h5_path = tmp(&format!("h5py-interop-{name}.h5"));
        let dump_dir = tmp(&format!("h5py-dump-{name}"));
        let _ = std::fs::remove_dir_all(&dump_dir);
        std::fs::create_dir_all(&dump_dir).unwrap();

        e.call("export_results", &json!({ "path": h5_path.to_str().unwrap() })).unwrap();

        let out = std::process::Command::new(py)
            .args(["-c", script, h5_path.to_str().unwrap(), dump_dir.to_str().unwrap()])
            .output()
            .expect("run python h5py check");
        assert!(
            out.status.success(),
            "h5py check failed on {name}:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );

        let read_f64_bin = |fname: &str| -> Vec<f64> {
            let bytes = std::fs::read(dump_dir.join(fname)).unwrap();
            assert_eq!(bytes.len() % 8, 0);
            bytes.as_chunks::<8>().0.iter().map(|&c| f64::from_le_bytes(c)).collect()
        };

        let meta_val: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(dump_dir.join("meta.json")).unwrap()).unwrap();
        assert_eq!(meta_val["packed_shape"][4].as_u64().unwrap() as usize, expected_ri, "{name} R dimension in h5py");

        let mesh = e.model().unwrap().solver_mesh();
        assert_f64_bits_eq(&format!("{name} h5py mesh/points"), &read_f64_bin("mesh_points.bin"), &mesh.points)
            .unwrap();
        let (orig, _) = e.result("ems").unwrap();
        let ord = if orig.order == 2 { 2 } else { 1 };
        assert_f64_bits_eq(
            &format!("{name} h5py results/ems/points"),
            &read_f64_bin("ems_points.bin"),
            &mesh.dof_coordinates(ord),
        )
        .unwrap();

        let mut expected_eigs = Vec::new();
        for s in &orig.solutions {
            let [re, im] = s.eigenvalue.unwrap_or([s.value.unwrap_or(0.0), 0.0]);
            expected_eigs.push(re);
            expected_eigs.push(im);
        }
        assert_f64_bits_eq(
            &format!("{name} h5py results/ems/eigenvalues"),
            &read_f64_bin("ems_eigenvalues.bin"),
            &expected_eigs,
        )
        .unwrap();
        assert_f64_bits_eq(
            &format!("{name} h5py solution_0/fields/00"),
            &read_f64_bin("sol0_f0.bin"),
            &orig.solutions[0].fields[0].values,
        )
        .unwrap();

        let packed = read_f64_bin("cell_fields_packed.bin");
        let n_cells = mesh.n_cells();
        let nloc = orig.solutions[0].cell_fields[0].nloc;
        let field_len = n_cells * nloc;
        let c_comp = orig.solutions[0].cell_fields.len();
        for (k, sol) in orig.solutions.iter().enumerate() {
            for (c, cf) in sol.cell_fields.iter().enumerate() {
                if expected_ri == 1 {
                    let base = (k * c_comp + c) * field_len;
                    let comp_re = &packed[base..base + field_len];
                    assert_f64_bits_eq(
                        &format!("{name} h5py packed mode {k} comp {c} ({}) re", cf.name),
                        comp_re,
                        &cf.re,
                    )
                    .unwrap();
                } else {
                    let base = (k * c_comp + c) * field_len * 2;
                    let pairs = packed[base..base + field_len * 2].as_chunks::<2>().0;
                    let comp_re: Vec<f64> = pairs.iter().map(|p| p[0]).collect();
                    let comp_im: Vec<f64> = pairs.iter().map(|p| p[1]).collect();
                    assert_f64_bits_eq(
                        &format!("{name} h5py packed mode {k} comp {c} ({}) re", cf.name),
                        &comp_re,
                        &cf.re,
                    )
                    .unwrap();
                    assert_f64_bits_eq(
                        &format!("{name} h5py packed mode {k} comp {c} ({}) im", cf.name),
                        &comp_im,
                        &cf.im,
                    )
                    .unwrap();
                }
            }
        }

        let _ = std::fs::remove_file(&h5_path);
        let _ = std::fs::remove_file(h5_path.with_extension("xdmf"));
        let _ = std::fs::remove_dir_all(&dump_dir);
    }
}

/// Parses the PE import directory (`IMAGE_DIRECTORY_ENTRY_IMPORT`) of a 64-bit Windows executable and returns the imported DLL names.
fn pe_imported_dlls(exe_bytes: &[u8]) -> Vec<String> {
    assert_eq!(&exe_bytes[0..2], b"MZ", "DOS MZ magic");
    let pe_off = u32::from_le_bytes(exe_bytes[0x3c..0x40].try_into().unwrap()) as usize;
    assert_eq!(&exe_bytes[pe_off..pe_off + 4], b"PE\0\0", "PE signature");
    let coff_off = pe_off + 4;
    let num_sections = u16::from_le_bytes(exe_bytes[coff_off + 2..coff_off + 4].try_into().unwrap()) as usize;
    let opt_size = u16::from_le_bytes(exe_bytes[coff_off + 16..coff_off + 18].try_into().unwrap()) as usize;
    let opt_off = coff_off + 20;
    let magic = u16::from_le_bytes(exe_bytes[opt_off..opt_off + 2].try_into().unwrap());
    let data_dir_off = match magic {
        0x20b => opt_off + 112, // PE32+
        0x10b => opt_off + 96,  // PE32
        other => panic!("unexpected PE optional header magic 0x{other:x}"),
    };
    let import_rva = u32::from_le_bytes(exe_bytes[data_dir_off + 8..data_dir_off + 12].try_into().unwrap());
    let import_size = u32::from_le_bytes(exe_bytes[data_dir_off + 12..data_dir_off + 16].try_into().unwrap());
    if import_rva == 0 || import_size == 0 {
        return Vec::new();
    }

    let sec_base = opt_off + opt_size;
    let rva_to_off = |rva: u32| -> Option<usize> {
        for i in 0..num_sections {
            let s = sec_base + i * 40;
            let vsize = u32::from_le_bytes(exe_bytes[s + 8..s + 12].try_into().unwrap());
            let va = u32::from_le_bytes(exe_bytes[s + 12..s + 16].try_into().unwrap());
            let raw_size = u32::from_le_bytes(exe_bytes[s + 16..s + 20].try_into().unwrap());
            let raw_ptr = u32::from_le_bytes(exe_bytes[s + 20..s + 24].try_into().unwrap());
            let span = vsize.max(raw_size);
            if rva >= va && rva < va + span {
                return Some((raw_ptr + (rva - va)) as usize);
            }
        }
        None
    };

    let mut dlls = Vec::new();
    let mut desc_off = rva_to_off(import_rva).expect("import table RVA in section");
    while desc_off + 20 <= exe_bytes.len() {
        let name_rva = u32::from_le_bytes(exe_bytes[desc_off + 12..desc_off + 16].try_into().unwrap());
        let first_thunk = u32::from_le_bytes(exe_bytes[desc_off + 16..desc_off + 20].try_into().unwrap());
        if name_rva == 0 && first_thunk == 0 {
            break;
        }
        if let Some(noff) = rva_to_off(name_rva) {
            let end = exe_bytes[noff..].iter().position(|&b| b == 0).unwrap() + noff;
            dlls.push(String::from_utf8_lossy(&exe_bytes[noff..end]).to_ascii_lowercase());
        }
        desc_off += 20;
    }
    dlls
}

#[test]
fn nyx_worker_has_zero_hdf5_dll_dependency_and_runs_without_anaconda_in_path() {
    let worker = PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"));
    let exe_bytes = std::fs::read(&worker).expect("read nyx-worker.exe");
    let dlls = pe_imported_dlls(&exe_bytes);
    eprintln!("nyx-worker imported DLLs: {dlls:?}");
    assert!(!dlls.is_empty(), "PE import parser must see standard Windows OS DLLs");
    for d in &dlls {
        assert!(!d.contains("hdf5"), "nyx-worker must not import any hdf5 DLL, found {d} in {dlls:?}");
        assert!(!d.contains("zlib"), "nyx-worker must not import any C zlib DLL, found {d} in {dlls:?}");
        assert!(!d.contains("mpi"), "nyx-worker must not import any MPI DLL, found {d} in {dlls:?}");
    }

    // Also run nyx-worker in a clean environment with Anaconda / HDF5_DIR stripped from PATH
    let job_dir = tmp("worker-clean-path-job");
    let _ = std::fs::remove_dir_all(&job_dir);
    std::fs::create_dir_all(&job_dir).unwrap();

    let mut ref_eng = Engine::new();
    ref_eng.call("project_open", &json!({ "path": example("04_thermal_bar") })).unwrap();
    ref_eng.call("solve_run", &json!({ "isolate": false })).unwrap();
    let lib = ref_eng.library().unwrap();
    let materials = nyx_engine::worker::material_texts(&lib);
    std::fs::write(job_dir.join("project.nyx.json"), ref_eng.project().to_json()).unwrap();
    std::fs::write(
        job_dir.join("job.json"),
        json!({ "targets": ["dfs", "sms"], "materials": materials, "warm": false }).to_string(),
    )
    .unwrap();

    let mut child = std::process::Command::new(&worker)
        .arg("--job")
        .arg(&job_dir)
        .env("PATH", r"C:\Windows\System32")
        .env_remove("HDF5_DIR")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn nyx-worker without Anaconda on PATH");
    let _keep_stdin = child.stdin.take();
    let out = child.wait_with_output().expect("wait for nyx-worker");
    assert!(
        out.status.success(),
        "nyx-worker failed without Anaconda PATH:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    let stored = nyx_engine::h5io::read_results(&job_dir.join("results.h5")).unwrap();
    for (code, r) in &stored.results {
        let (orig, _) = ref_eng.result(code).unwrap();
        assert_eq!(r.solutions, orig.solutions, "worker clean-PATH {code} solutions match bit-for-bit");
    }

    let _ = std::fs::remove_dir_all(&job_dir);
}

fn write_synthetic_packed_h5(path: &Path, n_modes: usize, n_cells: usize, nloc: usize) {
    let mut fb = FileBuilder::new();
    fb.set_attr("writer", AttrValue::VarLenString("serial".into()));
    fb.set_attr("processes", AttrValue::I64(1));
    fb.set_attr("format", AttrValue::VarLenString("nyx-results".into()));
    fb.set_attr("version", AttrValue::I32(4));
    let proj = nyx_engine::Project::new("synthetic_budget", 3).to_json();
    fb.create_dataset("project").with_u8_data(proj.as_bytes()).with_shape(&[proj.len() as u64]);
    fb.set_attr("fields_saved", AttrValue::I64(1));

    let mut m = fb.create_group("mesh");
    m.create_dataset("points")
        .with_f64_data(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        .with_shape(&[4, 3]);
    m.create_dataset("cells").with_u32_data(&[0, 1, 2, 3]).with_shape(&[1, 4]);
    m.create_dataset("cell_labels").with_i32_data(&[1]).with_shape(&[1]);
    m.create_dataset("facets").with_u32_data(&[0, 1, 2]).with_shape(&[1, 3]);
    m.create_dataset("facet_labels").with_i32_data(&[1]).with_shape(&[1]);
    fb.add_group(m.finish());

    let mut root = fb.create_group("results");
    let mut g = root.create_group("ems");
    g.set_attr("physics", AttrValue::VarLenString("ems".into()));
    g.set_attr("study", AttrValue::VarLenString("cavity".into()));
    g.set_attr("order", AttrValue::I64(2));
    g.set_attr("dim", AttrValue::I64(3));
    g.set_attr("dofs", AttrValue::I64(10));
    g.set_attr("seconds", AttrValue::F64(0.1));
    g.set_attr("n_solutions", AttrValue::I64(n_modes as i64));
    g.set_attr("field_layout", AttrValue::VarLenString("L2ShuffleLzfPacked".into()));
    g.create_dataset("scalars").with_u8_data(b"[]").with_shape(&[2]);
    g.create_dataset("log").with_u8_data(b"[]").with_shape(&[2]);
    g.create_dataset("eigenvalues").with_f64_data(&vec![1.0; n_modes * 2]).with_shape(&[n_modes as u64, 2]);

    let meta = vec![nyx_storage::PackedComponentMeta { name: "Ex".into(), unit: "V/m".into(), nloc, has_im: false }];
    let meta_bytes = serde_json::to_vec(&meta).unwrap();
    g.create_dataset("cell_fields_meta").with_u8_data(&meta_bytes).with_shape(&[meta_bytes.len() as u64]);

    let slice_elems = n_cells * nloc;
    let mut packed = vec![0.0_f64; n_modes * slice_elems];
    for k in 0..n_modes {
        packed[k * slice_elems] = (k + 1) as f64;
        packed[(k + 1) * slice_elems - 1] = -((k + 1) as f64);
    }
    g.create_dataset("cell_fields_packed")
        .with_f64_data(&packed)
        .with_shape(&[n_modes as u64, 1, n_cells as u64, nloc as u64, 1])
        .with_chunks(&[1, 1, n_cells as u64, nloc as u64, 1])
        .with_shuffle()
        .with_lzf();

    for k in 0..n_modes {
        let mut sg = g.create_group(&format!("solution_{k}"));
        sg.set_attr("label", AttrValue::VarLenString(format!("mode {k}")));
        sg.set_attr("eig_re", AttrValue::F64((k + 1) as f64));
        sg.set_attr("eig_im", AttrValue::F64(0.0));
        let mut fg = sg.create_group("fields");
        fg.set_attr("order", AttrValue::VarLenString("[]".into()));
        sg.add_group(fg.finish());
        g.add_group(sg.finish());
    }
    root.add_group(g.finish());
    fb.add_group(root.finish());
    fb.write(path).unwrap();
}

#[test]
fn nyx_h5_max_bytes_child_probe() {
    let Ok(target) = std::env::var("NYX_TEST_H5_BUDGET_PROBE_PATH") else {
        return;
    };
    let expect_ok = std::env::var("NYX_TEST_H5_BUDGET_EXPECT_OK").as_deref() == Ok("1");
    let res = nyx_engine::h5io::read_results(Path::new(&target));
    if expect_ok {
        let stored = res.expect("expected read_results to succeed under NYX_H5_MAX_BYTES");
        assert_eq!(stored.results[0].1.solutions.len(), 4);
    } else {
        let err = res.err().expect("expected read_results to fail under NYX_H5_MAX_BYTES");
        assert!(err.contains("exceeds read budget"), "unexpected error: {err}");
    }
}

#[test]
fn per_mode_packed_budget_allows_multi_mode_above_512_mib_rejects_oversized_slice_and_honors_env_limit() {
    // 1. 4 modes x 140 MB = 560 MB (> 512 MiB = 536,870,912 B total), with per-mode slice = 140 MB (<= 512 MiB): MUST PASS!
    let path_ok = tmp("packed-560mb-4modes.h5");
    write_synthetic_packed_h5(&path_ok, 4, 1_750_000, 10);
    let stored = nyx_engine::h5io::read_results(&path_ok)
        .expect("560 MB total packed result with 140 MB per-mode slices must pass 512 MiB slice / 8 GiB total budget");
    assert_eq!(stored.results.len(), 1);
    let sols = &stored.results[0].1.solutions;
    assert_eq!(sols.len(), 4);
    for (k, sol) in sols.iter().enumerate() {
        assert_eq!(sol.cell_fields[0].re.len(), 17_500_000);
        assert_eq!(sol.cell_fields[0].re[0], (k + 1) as f64);
        assert_eq!(*sol.cell_fields[0].re.last().unwrap(), -((k + 1) as f64));
    }
    drop(stored);

    // 2. NYX_H5_MAX_BYTES environment variable overrides max_total_result_bytes in `read_results`
    let cur_exe = std::env::current_exe().unwrap();
    let out_reject = std::process::Command::new(&cur_exe)
        .args(["--exact", "nyx_h5_max_bytes_child_probe", "--nocapture"])
        .env("NYX_TEST_H5_BUDGET_PROBE_PATH", &path_ok)
        .env("NYX_TEST_H5_BUDGET_EXPECT_OK", "0")
        .env("NYX_H5_MAX_BYTES", "200000000") // 200 MB: slice (140 MB) fits, total (560 MB) exceeds!
        .output()
        .unwrap();
    assert!(
        out_reject.status.success(),
        "child probe reject failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&out_reject.stdout),
        String::from_utf8_lossy(&out_reject.stderr)
    );

    let out_allow = std::process::Command::new(&cur_exe)
        .args(["--exact", "nyx_h5_max_bytes_child_probe", "--nocapture"])
        .env("NYX_TEST_H5_BUDGET_PROBE_PATH", &path_ok)
        .env("NYX_TEST_H5_BUDGET_EXPECT_OK", "1")
        .env("NYX_H5_MAX_BYTES", "700000000") // 700 MB > 560 MB total: allowed!
        .output()
        .unwrap();
    assert!(
        out_allow.status.success(),
        "child probe allow failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&out_allow.stdout),
        String::from_utf8_lossy(&out_allow.stderr)
    );
    let _ = std::fs::remove_file(&path_ok);

    // 3. Corrupted / oversized single-mode slice > 512 MiB (1 mode x 544 MB > 536,870,912 B): MUST FAIL!
    let path_bad_slice = tmp("packed-544mb-1mode.h5");
    write_synthetic_packed_h5(&path_bad_slice, 1, 6_800_000, 10);
    let err =
        nyx_engine::h5io::read_results(&path_bad_slice).err().expect("single-mode slice > 512 MiB must be rejected");
    assert!(err.contains("exceeds read budget"), "expected 'exceeds read budget', got: {err}");
    let _ = std::fs::remove_file(&path_bad_slice);
}

#[test]
fn write_results_with_save_fields_false_omits_mesh_and_fields_and_reads_cleanly() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("15_soi_wire") })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let mesh = e.model().unwrap().solver_mesh();
    let (orig, _) = e.result("ems").unwrap();
    let refs = vec![("ems".to_string(), orig.as_ref())];

    let path = tmp("save-fields-false.h5");
    nyx_engine::h5io::write_results_with_options(
        &path,
        e.project(),
        &mesh,
        &refs,
        FieldStorageLayout::L2ShuffleLzfPacked,
        false,
    )
    .unwrap();

    let bytes = nyx_engine::h5io::file_bytes(&path);
    assert!(bytes < 64 * 1024, "field-less results.h5 must be tiny, got {bytes} B");

    let raw = File::open_streaming(&path).unwrap();
    assert!(raw.group("mesh").is_err(), "/mesh group omitted when save_fields == false");
    let g = raw.group("results/ems").unwrap();
    assert!(g.dataset("points").is_err() && g.dataset("cells").is_err());
    assert!(g.dataset("cell_fields_packed").is_err());
    let sg = g.group("solution_0").unwrap();
    assert!(sg.group("fields").is_err() && sg.group("cell_fields").is_err());

    let stored = nyx_engine::h5io::read_results(&path).unwrap();
    assert!(!stored.fields_saved);
    assert_eq!(stored.mesh.n_cells(), 0);
    assert_eq!(stored.project, e.project().to_json());
    assert_eq!(stored.results.len(), 1);
    let r = &stored.results[0].1;
    assert_eq!(r.scalars, orig.scalars);
    assert_eq!(r.dofs, orig.dofs);
    assert_eq!(r.solutions.len(), orig.solutions.len());
    for (sa, sb) in r.solutions.iter().zip(&orig.solutions) {
        assert_eq!(sa.label, sb.label);
        assert_eq!(sa.value, sb.value);
        assert_eq!(sa.eigenvalue, sb.eigenvalue);
        assert!(sa.fields.is_empty() && sa.cell_fields.is_empty());
    }

    let summary = nyx_engine::h5io::read_results_summary(&path).unwrap();
    assert_eq!(summary.dim, mesh.dim);
    assert_eq!(summary.points, mesh.n_points());
    assert_eq!(summary.cells, mesh.n_cells());
    assert_eq!(summary.results[0].1.scalars, orig.scalars);

    // `open_results` also succeeds cleanly when `!fields_saved` (re-meshing project and restoring scalars/eigenvalues/logs)
    let mut reopened = Engine::new();
    let opened = reopened.call("open_results", &json!({ "path": path.to_str().unwrap() })).unwrap();
    assert_eq!(opened["results"][0]["physics"], "ems");
    let (loaded_res, stale) = reopened.result("ems").unwrap();
    assert!(!stale);
    assert_eq!(loaded_res.scalars, orig.scalars);
    assert_eq!(loaded_res.solutions.len(), orig.solutions.len());

    let _ = std::fs::remove_file(&path);
}

#[test]
fn results_h5_is_byte_deterministic_and_guards_forged_n_solutions() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("04_thermal_bar") })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let mesh = e.model().unwrap().solver_mesh();
    let results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)> = e
        .project()
        .physics
        .iter()
        .filter_map(|s| e.result(&s.physics).map(|(r, _)| (s.physics.clone(), r.clone())))
        .collect();
    let refs: Vec<(String, &nyx_physics::SolveResult)> = results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();

    let p1 = tmp("det-1.h5");
    let p2 = tmp("det-2.h5");
    nyx_engine::h5io::write_results(&p1, e.project(), &mesh, &refs).unwrap();
    std::thread::sleep(std::time::Duration::from_millis(20));
    nyx_engine::h5io::write_results(&p2, e.project(), &mesh, &refs).unwrap();

    let b1 = std::fs::read(&p1).unwrap();
    let b2 = std::fs::read(&p2).unwrap();
    assert_eq!(b1, b2, "results.h5 must be byte-deterministic across repeated writes");
    let raw = File::open_streaming(&p1).unwrap();
    let attrs = raw.root().attrs().unwrap();
    assert!(!attrs.contains_key("created_utc"), "created_utc must not be written into deterministic results.h5");
    let _ = std::fs::remove_file(&p1);
    let _ = std::fs::remove_file(&p2);

    // Forged huge `n_solutions` attribute with `!fields_saved` must be rejected by read budget before Vec::with_capacity OOM
    let forged = tmp("forged-n-sols.h5");
    let mut fb = FileBuilder::new();
    fb.set_attr("writer", AttrValue::VarLenString("serial".into()));
    fb.set_attr("processes", AttrValue::I64(1));
    fb.set_attr("format", AttrValue::VarLenString("nyx-results".into()));
    fb.set_attr("version", AttrValue::I32(4));
    fb.set_attr("fields_saved", AttrValue::I64(0));
    let proj = nyx_engine::Project::new("forged", 2).to_json();
    fb.create_dataset("project").with_u8_data(proj.as_bytes()).with_shape(&[proj.len() as u64]);
    let mut root = fb.create_group("results");
    let mut g = root.create_group("ems");
    g.set_attr("physics", AttrValue::VarLenString("ems".into()));
    g.set_attr("study", AttrValue::VarLenString("modes".into()));
    g.set_attr("order", AttrValue::I64(2));
    g.set_attr("dim", AttrValue::I64(2));
    g.set_attr("dofs", AttrValue::I64(10));
    g.set_attr("seconds", AttrValue::F64(0.1));
    g.set_attr("n_solutions", AttrValue::I64(1_000_000_000));
    g.create_dataset("scalars").with_u8_data(b"[]").with_shape(&[2]);
    g.create_dataset("log").with_u8_data(b"[]").with_shape(&[2]);
    root.add_group(g.finish());
    fb.add_group(root.finish());
    fb.write(&forged).unwrap();

    let err = nyx_engine::h5io::read_results_with_budget(&forged, ReadBudgetConfig::default())
        .err()
        .expect("forged n_solutions=1e9 must be rejected");
    assert!(err.contains("exceeds read budget") || err.contains("n_solutions"), "{err}");
    let _ = std::fs::remove_file(&forged);
}
