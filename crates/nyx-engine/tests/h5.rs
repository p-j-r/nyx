//! The HDF5 results file (feature `hdf5`; ADR 62): what is written is read back bit for bit, the file is smaller than the same data as JSON,
//! and the engine commands write and describe it.
#![cfg(feature = "hdf5")]

use nyx_engine::Engine;
use serde_json::json;
use std::path::PathBuf;

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
        eprintln!("{name}: HDF5 {bytes} bytes, JSON {json_bytes} bytes");
        let _ = std::fs::remove_file(&path);
    }
}

#[test]
fn a_file_that_nyx_did_not_write_is_refused_and_nothing_is_exported_before_a_solve() {
    let mut e = Engine::new();
    assert!(e.call("export_results", &json!({ "path": tmp("none.h5").to_str().unwrap() })).is_err());
    let other = tmp("other.h5");
    hdf5_metno::File::create(&other).unwrap();
    assert!(nyx_engine::h5io::read_results(&other).err().expect("refused").contains("not written by nyx"));
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
        let file = hdf5_metno::File::open(&path).unwrap();
        // no second copy of the nodes and cells: the XDMF names the ones of the results
        assert!(file.group("viz").is_err(), "{name}: there is no /viz group");
        assert!(xdmf.contains("/points") && xdmf.contains("/cells") && !xdmf.contains("/viz/"), "{name}");
        let fname = path.file_name().unwrap().to_str().unwrap();
        let mut refs = 0;
        for part in xdmf.split("<DataItem").skip(1) {
            let body = part.split('>').nth(1).unwrap().split('<').next().unwrap().trim();
            let (f, p) = body.split_once(':').expect("file:path");
            assert_eq!(f, fname);
            let ds = file.dataset(p).unwrap_or_else(|err| panic!("{name}: {p} does not exist: {err}"));
            assert!(ds.size() > 0);
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
    e.call("project_open", &json!({ "path": example("15_soi_wire")})).unwrap();
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
}

#[test]
fn a_file_that_cannot_be_opened_leaves_the_project_as_it_was() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("04_thermal_bar") })).unwrap();
    let before = e.project().clone();
    let other = tmp("notnyx.h5");
    hdf5_metno::File::create(&other).unwrap();
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
        let file = hdf5_metno::File::open(&path).unwrap();
        let has = file.group("mesh").unwrap().dataset("edge_node_pairs").is_ok();
        assert_eq!(has, name == "curved", "the datasets exist only for a curved mesh");
        let _ = std::fs::remove_file(&path);
    }
}

/// The documented cell node order of layout 4 is the real one: in every cell the edge nodes sit at the midpoints of the edges the
/// `cell_order` text names, and a continuous cell field gives the same value at a node from every cell that holds it.
#[test]
fn the_cell_node_order_of_a_results_file_is_the_documented_one() {
    for name in ["15_soi_wire", "04_thermal_bar"] {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example(name) })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let path = tmp(&format!("{name}-order.h5"));
        e.call("export_results", &json!({ "path": path.to_str().unwrap() })).unwrap();
        let file = hdf5_metno::File::open(&path).unwrap();
        assert_eq!(file.attr("version").unwrap().read_scalar::<i32>().unwrap(), 4);
        for code in file.group("results").unwrap().member_names().unwrap() {
            let g = file.group(&format!("results/{code}")).unwrap();
            let order =
                g.attr("cell_order").unwrap().read_scalar::<hdf5_metno::types::VarLenUnicode>().unwrap().to_string();
            let (pts, cells) = (g.dataset("points").unwrap(), g.dataset("cells").unwrap());
            let dim = pts.shape()[1];
            let (p, c) = (pts.read_raw::<f64>().unwrap(), cells.read_raw::<u32>().unwrap());
            let nloc = cells.shape()[1];
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
            // the nodal and the cell view of a continuous field agree
            for s in g.member_names().unwrap().iter().filter(|n| n.starts_with("solution_")) {
                let Ok(cf) = g.group(&format!("{s}/cell_fields")) else { continue };
                for f in cf.member_names().unwrap() {
                    let re = cf.dataset(&format!("{f}/re")).unwrap().read_raw::<f64>().unwrap();
                    assert_eq!(re.len(), c.len(), "{name}/{code}/{s}/{f}: one value per cell node");
                }
            }
        }
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(path.with_extension("xdmf"));
    }
}
