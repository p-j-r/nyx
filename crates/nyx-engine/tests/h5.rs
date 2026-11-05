//! The HDF5 results file (feature `hdf5`; ADR 62): what is written is read back bit for bit, the file is smaller than the same data as JSON,
//! and the engine commands write and describe it.
#![cfg(feature = "hdf5")]

use nyx_engine::Engine;
use serde_json::json;
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples").join(format!("{name}.nyx.json")).to_str().unwrap().to_string()
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
        assert!(d["results"].as_array().unwrap().len() >= 1 && d["cells"].as_u64().unwrap() > 0, "{d}");
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
    assert!(!artifacts.iter().any(|a| a["path"].as_str().unwrap().ends_with(".json")), "no JSON copy of the fields: {runs}");
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
        let file = hdf5_metno::File::open(&path).unwrap();
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
