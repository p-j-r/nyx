// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Interface contracts of femcore (docs/architecture/adr/0019): what other crates, saved files and external hosts rely on, written down
//! as tests that do not depend on frozen floating-point values (the golden files of `golden_frames.rs` do, see ADR 10).
//!
//! * `frames`: the constants and name tables of `section` (ADR 1).
//! * `mesh`: label bases, node orders, no curved-node storage, the `.fmsh` version (ADRs 4, 5, 16, 18).
//! * `kernel_stages`: the public stage functions (ADR 9).
//!
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p femcore --test contracts` rewrites the exemplar. Do that only together with a bump of
//! `docs/architecture/contracts/VERSIONS.md` and an ADR note (see ADR 19), never to silence a failure.
mod common;

use femcore::material::M3;
use femcore::props::*;
use femcore::section;
use serde_json::{json, Value};
use std::path::PathBuf;

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

/// Like `check_exemplar`, but compares the parsed JSON values, never the text or the key order (ADR 26).
fn check_exemplar_value(name: &str, got: &Value) {
    let path = contract_file(name);
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, format!("{}\n", serde_json::to_string_pretty(got).unwrap())).unwrap();
        return;
    }
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let want: Value = serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert_eq!(&want, got, "{name} no longer matches the contract exemplar; this is an interface change (ADR 19)");
}

// ---------------------------------------------------------------------------------------------------------------- solve_result

// ---------------------------------------------------------------------------------------------------------------- frames

fn diag(a: f64, b: f64, c: f64) -> M3 {
    [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
}

#[test]
fn frames_constants_of_the_wafer_frame() {
    // ADR 1: waves propagate along +y; local triple (lateral, normal, propagation) = wafer (x, z, y)
    assert_eq!(section::PROP_SIGN, 1.0);
    assert_eq!(section::LOCAL_AXES, [(0, 1.0), (2, 1.0), (1, 1.0)]);
    assert_eq!(section::PLANE_VOIGT, [0, 2, 1, 5, 4, 3]);
    assert!(section::is_left_handed(section::PROP_SIGN), "the mode solvers' local triple is left-handed");
    assert!(!section::is_left_handed(-section::PROP_SIGN));
}

#[test]
fn frames_names_and_points_of_a_cross_section_and_of_3d() {
    assert_eq!(section::wafer_axes(2), &[0, 2]);
    assert_eq!(section::wafer_axes(3), &[0, 1, 2]);
    assert_eq!(section::component_names(2), &["x", "z"]);
    assert_eq!(section::component_names(3), &["x", "y", "z"]);
    assert_eq!(section::wafer_point(&[1.0, 2.0], 2), [1.0, 0.0, 2.0], "a 2D mesh point (x, z) is the wafer point (x, 0, z)");
    assert_eq!(section::wafer_point(&[1.0, 2.0, 3.0], 3), [1.0, 2.0, 3.0]);
    assert_eq!([section::local_component(0), section::local_component(1), section::local_component(2)], [("x", 1.0), ("z", 1.0), ("y", 1.0)]);
    let names: Vec<&str> = (0..6).map(|i| section::stress_name(2, i)).collect();
    assert_eq!(names, vec!["sxx", "szz", "syy", "sxz", "syz", "sxy"], "Voigt order of a cross-section is (xx, zz, yy, xz, zy, xy)");
    let names3: Vec<&str> = (0..6).map(|i| section::stress_name(3, i)).collect();
    assert_eq!(names3, vec!["sxx", "syy", "szz", "sxy", "syz", "sxz"]);
}

#[test]
fn frames_tensors_are_permuted_only_in_section() {
    // wafer diag(1, 2, 3) = (xx, yy, zz): a cross-section reads it in the order (x, z, y)
    assert_eq!(section::plane_tensor(&diag(1.0, 2.0, 3.0), 2), diag(1.0, 3.0, 2.0));
    assert_eq!(section::plane_tensor(&diag(1.0, 2.0, 3.0), 3), diag(1.0, 2.0, 3.0));
    // the mode solvers' local order is the same permutation, with the propagation sign on the third axis
    assert_eq!(section::local_tensor(&diag(1.0, 2.0, 3.0)), diag(1.0, 3.0, 2.0));
    // an off-diagonal wafer entry xy (lateral-propagation) flips sign with the propagation direction, xz does not
    let mut t = diag(1.0, 2.0, 3.0);
    t[0][1] = 0.25;
    t[1][0] = 0.25;
    t[0][2] = 0.5;
    t[2][0] = 0.5;
    let fwd = section::local_tensor_signed(&t, 1.0);
    let bwd = section::local_tensor_signed(&t, -1.0);
    assert_eq!(fwd[0][2], 0.25, "lateral-propagation entry is wafer xy");
    assert_eq!(bwd[0][2], -0.25);
    assert_eq!(fwd[0][1], 0.5, "lateral-normal entry is wafer xz");
    assert_eq!(bwd[0][1], 0.5);
}

#[test]
fn frames_pml_directions_use_x_and_z_on_a_cross_section() {
    let p = Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Minus, ..Pml::default() };
    assert_eq!(section::mesh_pml_dirs(&p, 2), [PmlDir::Plus, PmlDir::Minus, PmlDir::Off]);
    assert_eq!(section::mesh_pml_dirs(&p, 3), [PmlDir::Plus, PmlDir::Off, PmlDir::Minus]);
}

// ---------------------------------------------------------------------------------------------------------------- mesh

#[test]
fn mesh_labels_are_one_based_and_cells_store_vertices_only() {
    let m = common::rect(1.0, 0.5, 4, 2);
    // curved edges are a separate list (ADR 82): an empty list is the straight mesh of ADR 18
    assert!(!m.has_curved());
    assert!(matches!(femcore::mesh::Mesh::from_solver_mesh(&m).unwrap(), femcore::mesh::Mesh::D2(c) if c.edge_nodes.is_empty()));
    assert_eq!(m.dim, 2);
    assert_eq!(m.points.len(), 2 * m.n_points());
    assert_eq!(m.cells.len(), 3 * m.n_cells(), "a triangle stores its three vertices, no mid-edge nodes (ADR 18)");
    assert_eq!(m.facets.len(), 2 * m.facet_labels.len());
    assert!(m.cell_labels.iter().all(|&l| l >= 1), "domain numbers are 1-based");
    assert!(m.facet_labels.iter().all(|&l| l >= 1), "boundary numbers are 1-based");
    // P2 numbering: vertices first, then one dof per edge
    assert_eq!(m.n_dofs(1), m.n_points());
    assert_eq!(m.n_dofs(2), m.n_points() + m.n_edges());
    let b = common::block(1.0, 1.0, 0.5, 3, 3, 2);
    assert_eq!(b.dim, 3);
    assert_eq!(b.cells.len(), 4 * b.n_cells());
    assert_eq!(b.facets.len(), 3 * b.facet_labels.len());
}

#[test]
fn mesh_cells_are_positively_oriented_in_2d() {
    let m = common::rect(1.0, 0.5, 4, 2);
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let p = |i: usize| (m.points[v[i] as usize * 2], m.points[v[i] as usize * 2 + 1]);
        let (a, b, d) = (p(0), p(1), p(2));
        let area2 = (b.0 - a.0) * (d.1 - a.1) - (d.0 - a.0) * (b.1 - a.1);
        assert!(area2.abs() > 0.0);
    }
}

mod mixed {
    use femcore::mesh::cells::CellMesh;
    use femcore::adapt::work::simplex_to_adapt;
    use femcore::mesh::SimplexMesh;
    use femcore::meshing::pml_box::{build_pml, PmlSpec};

    pub fn pml_cells() -> CellMesh<2> {
        let plate = femcore::meshing::geomdrive::plate(1.0, 1.0, 0.125).unwrap().mesh;
        let phys = simplex_to_adapt::<2>(&plate).unwrap();
        let spec = PmlSpec { outer_label: Some(1), ..PmlSpec::<2>::uniform(0.25, 3) };
        build_pml(&phys, &spec).unwrap().mesh
    }

    pub fn split(cells: &CellMesh<2>) -> SimplexMesh {
        SimplexMesh::from_cells_split(cells).unwrap().0
    }
}

#[test]
fn mesh_domain_and_boundary_numbers_are_one_based_in_the_cell_mesh_and_the_split_mesh() {
    // ADR 5 / ROUGH_SPOTS R09, resolved in phase 2 (ADR 32): one mesh type, domain and boundary numbers both from 1;
    // `from_cells_split` copies them.
    use femcore::mesh::cells::CellType;
    let cells = mixed::pml_cells();
    assert!(cells.types.contains(&CellType::Tri) && cells.types.contains(&CellType::Quad), "physical triangles and PML quadrilaterals");
    assert!(cells.facet_label.iter().all(|&d| d >= 1) && cells.facet_label.contains(&1), "boundary numbers start at 1");
    assert!(cells.label.iter().all(|&d| d >= 1), "domain numbers start at 1 (0 is the exterior)");
    let mesh = mixed::split(&cells);
    let mut want: Vec<i32> = cells.facet_label.clone();
    want.sort_unstable();
    want.dedup();
    let mut got = mesh.boundary_labels();
    got.sort_unstable();
    got.dedup();
    assert_eq!(got, want, "the split mesh keeps the boundary numbers");
}

#[test]
fn mesh_quadrilaterals_are_counter_clockwise_and_split_exactly() {
    // the native quad order is the polygon (0, 1, 2, 3), counter-clockwise (ADR 5)
    use femcore::mesh::cells::CellType;
    let cells = mixed::pml_cells();
    let mut n_quads = 0;
    for c in (0..cells.n_cells()).filter(|&c| cells.types[c] == CellType::Quad) {
        n_quads += 1;
        let v = cells.cell(c);
        let ring: Vec<[f64; 2]> = (0..4).map(|i| cells.points[v[i] as usize]).collect();
        let area2: f64 = (0..4).map(|i| ring[i][0] * ring[(i + 1) % 4][1] - ring[(i + 1) % 4][0] * ring[i][1]).sum();
        assert!(area2 > 0.0, "quad {c} is not counter-clockwise");
    }
    assert!(n_quads > 0);
    // the split into triangles adds no vertex and fills the same area: the 1 x 1 plate with a 0.25 layer is 1.5 x 1.5
    let mesh = mixed::split(&cells);
    assert_eq!(mesh.n_points(), cells.n_points());
    let area: f64 = (0..mesh.n_cells())
        .map(|c| {
            let v = mesh.cell(c);
            let p = |i: usize| (mesh.points[v[i] as usize * 2], mesh.points[v[i] as usize * 2 + 1]);
            0.5 * ((p(1).0 - p(0).0) * (p(2).1 - p(0).1) - (p(2).0 - p(0).0) * (p(1).1 - p(0).1)).abs()
        })
        .sum();
    assert!((area - 2.25).abs() < 1e-9, "area {area}");
}

#[test]
fn mesh_pml_generator_keeps_the_physical_vertex_and_cell_numbers() {
    // The one PML generator (`meshing::pml_box::build_pml`, ADR 32) appends the layers to the physical mesh: vertex ids, cell order,
    // connectivity and labels of the physical part are unchanged, so a field on the physical mesh maps back by id (the solve-adapt loop
    // of the mesher feature chart, item P3, relies on it).
    use femcore::adapt::work::simplex_to_adapt;
    use femcore::meshing::pml_box::{build_pml, PmlSpec};
    fn check<const D: usize>(phys: &femcore::adapt::work::SimplexMesh<D>, pml: &femcore::meshing::pml_box::PmlMesh<D>) {
        assert_eq!(pml.n_physical_points, phys.points.len());
        assert_eq!(pml.n_physical_cells, phys.n_cells());
        assert!(pml.mesh.points.len() > phys.points.len(), "layer vertices are appended");
        assert_eq!(&pml.mesh.points[..phys.points.len()], &phys.points[..], "physical coordinates by id");
        for c in 0..phys.n_cells() {
            let mut want = phys.cells[c * (D + 1)..(c + 1) * (D + 1)].to_vec();
            let mut got = pml.mesh.cell(c).to_vec();
            want.sort_unstable();
            got.sort_unstable();
            assert_eq!(got, want, "cell {c}");
            assert_eq!(pml.mesh.label[c], phys.cell_label[c], "label of cell {c}");
        }
    }
    let plate = femcore::meshing::geomdrive::plate(1.0, 1.0, 0.125).unwrap().mesh;
    let phys2 = simplex_to_adapt::<2>(&plate).unwrap();
    check(&phys2, &build_pml(&phys2, &PmlSpec { outer_label: Some(1), ..PmlSpec::<2>::uniform(0.25, 3) }).unwrap());
    let block = femcore::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.34, 3).unwrap().mesh;
    let phys3 = simplex_to_adapt::<3>(&block).unwrap();
    check(&phys3, &build_pml(&phys3, &PmlSpec { outer_label: Some(1), ..PmlSpec::<3>::uniform(0.5, 1) }).unwrap());
}

#[test]
fn mesh_cell_model_matches_the_exemplar() {
    // The one mesh model (ADR 32): node order per cell type, face lists, label base, solver view. Exemplar `mesh_cells.v1.json`.
    use femcore::mesh::cells::CellType;
    let shapes = [CellType::Tri, CellType::Quad, CellType::Tet, CellType::Prism, CellType::Hex];
    let cell_types: Vec<Value> = shapes
        .iter()
        .map(|t| {
            let order = match t {
                CellType::Tri | CellType::Quad => "counter-clockwise",
                CellType::Tet => "(v1 - v0, v2 - v0, v3 - v0) right-handed",
                CellType::Prism => "base triangle 0..2, top triangle 3..5 above it, (0,1,2,3) right-handed",
                CellType::Hex => "bottom quadrilateral 0..3 counter-clockwise seen from the top, top quadrilateral 4..7 above it",
            };
            json!({ "name": format!("{t:?}").to_lowercase(), "code": t.code(), "nodes": t.nv(), "dim": t.dim(), "node_order": order, "faces": t.faces().iter().map(|f| f.to_vec()).collect::<Vec<_>>() })
        })
        .collect();
    let got = json!({
        "model": "femcore::mesh::cells::CellMesh<D>, D = 2 or 3",
        "arrays": ["points [D]", "types", "offsets", "conn", "label", "facet_offsets", "facet_conn", "facet_label"],
        "label_base": { "domains": 1, "boundaries": 1, "exterior": 0 },
        "boundary_facets": { "2d": "2 nodes", "3d": "3 nodes (triangle) or 4 nodes (quadrilateral, counter-clockwise)" },
        "curved_nodes": "edge_nodes: Vec<(sorted vertex pair, position)> sorted by pair, empty by default; an edge not listed is straight, its mid-edge node is the midpoint (ADR 82, amends ADR 18)",
        "cell_types": cell_types,
        "file": { "format": ".fmsh", "version": 2 },
        "solver_view": {
            "type": "femcore::mesh::SimplexMesh (triangles or tetrahedra, edges, unknown numbering)",
            "from": "SimplexMesh::from_cells_split(&CellMesh<D>)",
            "quadrilateral_2d": "two triangles along the shorter diagonal",
            "quadrilateral_facet_3d": "two triangles along the diagonal from the lowest-numbered vertex",
            "hexahedron_prism_3d": "cone from the lowest-numbered vertex (SplitMode::Minimal) or centre vertices (SplitMode::Star)",
            "vertices": "original ids unchanged, new vertices appended",
            "labels": "copied"
        }
    });
    check_exemplar_value("mesh_cells.v1.json", &got);
}

#[test]
fn mesh_fmsh_round_trip_and_version() {
    use femcore::mesh::cells::CellMesh;
    use femcore::mesh::io::fmsh;
    use femcore::meshing::pml_box::{build_pml, PmlSpec};
    use femcore::adapt::work::simplex_to_adapt;
    let plate = femcore::meshing::geomdrive::plate(1.0, 1.0, 0.25).unwrap().mesh;
    let phys = simplex_to_adapt::<2>(&plate).unwrap();
    let m: CellMesh<2> = build_pml(&phys, &PmlSpec { outer_label: Some(1), ..PmlSpec::<2>::uniform(0.25, 2) }).unwrap().mesh;
    let bytes = fmsh::encode(&m);
    assert_eq!(&bytes[..4], b"FMSH");
    assert_eq!(bytes[4], 2, "`.fmsh` format version (ADR 16): bump it and add a reader for the old one when the layout changes");
    assert_eq!(bytes[5], 2, "dimension byte");
    let back = fmsh::decode::<2>(&bytes).unwrap();
    assert_eq!(fmsh::encode(&back), bytes, "lossless and deterministic: decode then encode gives the same bytes");
    assert_eq!(back.types, m.types);
    assert_eq!(back.label, m.label);
    assert_eq!(back.facet_label, m.facet_label);
    // the CRC catches a flipped bit
    let mut bad = bytes.clone();
    let mid = bad.len() / 2;
    bad[mid] ^= 0x10;
    assert!(fmsh::decode::<2>(&bad).is_err());
    // the dimension is part of the contract
    assert!(fmsh::decode::<3>(&bytes).is_err());
}

// ---------------------------------------------------------------------------------------------------------------- kernel_stages

#[test]
#[allow(clippy::type_complexity)]
fn kernel_stage_api_is_public() {
    // ADR 9: naming the functions with explicit fn-pointer checks makes removing or renaming any stage signature a compile error here.
    let _: fn(&femcore::mesh::SimplexMesh) -> femcore::fem::Geometry = femcore::fem::geometry;
    let _: fn(&femcore::mesh::SimplexMesh, usize, &[f64]) -> femcore::linalg::csr::RCsr =
        femcore::fem::assemble_laplace;
    let _: fn(&femcore::mesh::SimplexMesh, usize, &[f64]) -> femcore::linalg::csr::RCsr =
        femcore::fem::assemble_mass;
    let _: fn(&femcore::mesh::SimplexMesh, usize, &[f64]) -> Vec<f64> = femcore::fem::assemble_source;
    let _: fn(
        &femcore::mesh::SimplexMesh,
        usize,
        usize,
        femcore::props::LinearSolver,
        &femcore::linalg::csr::RCsr,
        &[f64],
        &[usize],
        &[f64],
    ) -> Result<(Vec<f64>, Option<String>), String> = femcore::fem::solve_spd;
    let _: fn(
        &femcore::linalg::csr::RCsr,
        &[f64],
        &[usize],
        &[f64],
    ) -> Result<Vec<f64>, String> = femcore::fem::solve_dirichlet;
    let _: fn(&femcore::mesh::SimplexMesh, &[f64], usize, usize) -> Vec<f64> =
        femcore::fields::gradients_at_nodes;
    let _: fn(&femcore::mesh::SimplexMesh, &[f64], usize, usize) -> Vec<f64> =
        femcore::fields::node_average;
    let _: fn(&femcore::mesh::SimplexMesh, &[f64], usize) -> Vec<f64> =
        femcore::fem::cell_values_at_nodes;
    let _: fn(&femcore::mesh::SimplexMesh, usize) -> femcore::edge::CellData = femcore::edge::cell_data;
    let _: fn(
        &[f64],
        &[[f64; 2]; 3],
        &[f64; 3],
        usize,
        &mut [[f64; 2]],
        &mut [f64],
    ) = femcore::edge::edge_basis;
    let _: fn(
        &geom::model::FeatureList,
        &femcore::meshing::MeshPlan,
    ) -> Result<femcore::meshing::geomdrive::MeshedGeom, String> = femcore::meshing::mesh_model;
    let _: fn(&femcore::mesh::SimplexMesh) -> Result<(), femcore::mesh::MeshIntegrityError> =
        femcore::mesh::SimplexMesh::check;
    let _: fn(
        &femcore::mesh::SimplexMesh,
        &[i32],
        usize,
        fn(usize, usize, [f64; 3], [f64; 3]) -> f64,
    ) -> femcore::linalg::csr::RCsr = femcore::farfield::assemble_boundary_mass_weighted;
    let _: fn(
        &femcore::mesh::SimplexMesh,
        &[i32],
        usize,
        usize,
        f64,
        [f64; 3],
        f64,
        femcore::farfield::FarFieldPermittivity<'_>,
    ) -> femcore::linalg::csr::RCsr = femcore::farfield::assemble_farfield_robin;
    let _: fn(
        &femcore::mesh::SimplexMesh,
        usize,
        &[Vec<f64>],
        femcore::rt::CellPermittivity<'_>,
        &[bool],
        &[bool],
        femcore::rt::RobinBoundary<'_>,
        Option<&[Vec<f64>]>,
    ) -> Result<femcore::rt::EquilibratedFluxResult, String> = femcore::rt::equilibrate_flux_patches;
    let _: fn(
        &femcore::mesh::SimplexMesh,
        &femcore::props::Periodic,
    ) -> Result<femcore::periodic::Match, String> = femcore::periodic::match_vertices;
}

// ---------------------------------------------------------------------------------------------------------------- mesh, 3D geom models (ADR 40)

#[test]
fn mesh_geom3d_numbers_cells_and_facets_with_the_model_label_map() {
    // A 3D geom model meshes through the same entry point as a 2D one (`meshing::mesh_model`): cells carry the domain number of their region
    // and boundary facets the boundary number of their geom face, both 1-based and equal to the model's `LabelMap`; the result is all
    // tetrahedra (`dim` 3, four vertices per cell, three per facet), no curved-node storage (ADR 18).
    let mut m = geom::model::Model::new_3d();
    let host = m.box_("air", (0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
    let ball = m.sphere("ball", (1.0, 1.0, 1.0), 0.5);
    m.overlay(&[("air", host), ("ball", ball)]);
    let g = femcore::meshing::mesh_model(&m.list, &femcore::meshing::MeshPlan::simple(0.25)).unwrap();
    assert_eq!(g.mesh.dim, 3);
    assert_eq!(g.mesh.cells.len(), 4 * g.mesh.n_cells());
    assert_eq!(g.mesh.facets.len(), 3 * g.mesh.n_facets());
    let domains: std::collections::BTreeSet<i32> = g.mesh.cell_labels.iter().copied().collect();
    assert_eq!(domains, g.labels.domains.iter().map(|d| d.id).collect(), "every region has cells, and no other number appears");
    assert_eq!(domains.iter().next(), Some(&1), "domain numbers start at 1");
    let faces: std::collections::BTreeSet<i32> = g.mesh.facet_labels.iter().copied().collect();
    assert_eq!(faces, g.labels.boundaries.iter().map(|b| b.id).collect(), "every face of the model has facets, and no other number appears");
    assert_eq!(faces.iter().next(), Some(&1), "boundary numbers start at 1");
    assert_eq!(g.boundary("air|ball").len(), 1);
}

// ---------------------------------------------------------------------------------------------------------------- 3D features (ADR 50 to 54)

/// A small validator for the subset of JSON Schema the mesh schemas use: type (also a list), const, enum, oneOf, required, properties,
/// additionalProperties, items, minItems, maxItems, minimum, exclusiveMinimum.
fn validate(schema: &Value, v: &Value, at: &str) -> Result<(), String> {
    let kind = |v: &Value| match v {
        Value::Null => "null",
        Value::Bool(_) => "boolean",
        Value::Number(n) if n.is_i64() || n.is_u64() => "integer",
        Value::Number(_) => "number",
        Value::String(_) => "string",
        Value::Array(_) => "array",
        Value::Object(_) => "object",
    };
    if let Some(c) = schema.get("const") {
        if c != v {
            return Err(format!("{at}: expected {c}, got {v}"));
        }
    }
    if let Some(e) = schema.get("enum") {
        if !e.as_array().unwrap().contains(v) {
            return Err(format!("{at}: {v} is not one of {e}"));
        }
    }
    if let Some(t) = schema.get("type") {
        let allowed: Vec<&str> = match t {
            Value::String(s) => vec![s.as_str()],
            Value::Array(a) => a.iter().map(|x| x.as_str().unwrap()).collect(),
            _ => vec![],
        };
        let k = kind(v);
        if !allowed.iter().any(|a| *a == k || (*a == "number" && k == "integer")) {
            return Err(format!("{at}: type {k}, expected {allowed:?}"));
        }
    }
    if let (Some(x), Some(m)) = (v.as_f64(), schema.get("minimum").and_then(Value::as_f64)) {
        if x < m {
            return Err(format!("{at}: {x} below {m}"));
        }
    }
    if let (Some(x), Some(m)) = (v.as_f64(), schema.get("exclusiveMinimum").and_then(Value::as_f64)) {
        if x <= m {
            return Err(format!("{at}: {x} not above {m}"));
        }
    }
    if let Some(one) = schema.get("oneOf") {
        let ok = one.as_array().unwrap().iter().filter(|s| validate(s, v, at).is_ok()).count();
        if ok != 1 {
            return Err(format!("{at}: matches {ok} of the alternatives, expected exactly one"));
        }
    }
    if let Value::Object(o) = v {
        for r in schema.get("required").and_then(Value::as_array).into_iter().flatten() {
            if !o.contains_key(r.as_str().unwrap()) {
                return Err(format!("{at}: missing '{}'", r.as_str().unwrap()));
            }
        }
        let props = schema.get("properties").and_then(Value::as_object);
        for (k, x) in o {
            match props.and_then(|p| p.get(k)) {
                Some(s) => validate(s, x, &format!("{at}.{k}"))?,
                None => match schema.get("additionalProperties") {
                    Some(Value::Bool(false)) => return Err(format!("{at}: unexpected '{k}'")),
                    Some(s @ Value::Object(_)) => validate(s, x, &format!("{at}.{k}"))?,
                    _ => {}
                },
            }
        }
    }
    if let Value::Array(a) = v {
        if let Some(n) = schema.get("minItems").and_then(Value::as_u64) {
            if (a.len() as u64) < n {
                return Err(format!("{at}: fewer than {n} items"));
            }
        }
        if let Some(n) = schema.get("maxItems").and_then(Value::as_u64) {
            if (a.len() as u64) > n {
                return Err(format!("{at}: more than {n} items"));
            }
        }
        if let Some(s) = schema.get("items") {
            for (i, x) in a.iter().enumerate() {
                validate(s, x, &format!("{at}[{i}]"))?;
            }
        }
    }
    Ok(())
}

/// The structure of a JSON value without its numbers: object keys and the kind of every leaf.
fn shape(v: &Value) -> Value {
    match v {
        Value::Object(o) => Value::Object(o.iter().map(|(k, x)| (k.clone(), shape(x))).collect()),
        Value::Array(a) => Value::Array(a.first().map(shape).into_iter().collect()),
        Value::Number(n) => json!(if n.is_i64() || n.is_u64() { "integer" } else { "number" }),
        Value::String(_) => json!("string"),
        Value::Bool(_) => json!("boolean"),
        Value::Null => json!("null"),
    }
}

fn ball_in_box_with_pml() -> geom::model::Model {
    let mut m = geom::model::Model::new_3d();
    let air = m.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let ball = m.sphere("ball", (0.0, 0.0, 0.0), 0.7);
    let win = m.box_("window", (-1.5, -1.5, -1.5), (1.5, 1.5, 1.5));
    m.overlay(&[("air", air), ("ball", ball)]);
    m.pml("pml", win, &[], 0.5, 2.0, 2.0);
    m
}

#[test]
fn mesh_geom3d_adaptive_and_structured_meshes_keep_the_label_map() {
    // ADR 51 and 52: whichever way a 3D model is meshed (simple, adaptive, structured layers), domain numbers are the model's, boundary numbers
    // are the model's, and the faces between regions are facets with the number of the geom face between them (ADR 40).
    let m = ball_in_box_with_pml();
    let simple = femcore::meshing::mesh_model(&m.list, &femcore::meshing::MeshPlan::simple(0.5)).unwrap();
    let plans = [
        femcore::meshing::MeshPlan::adaptive(0.6, 0.7),
        femcore::meshing::MeshPlan { pml: femcore::meshing::PmlMode::Structured, ..femcore::meshing::MeshPlan::simple(0.5) },
        femcore::meshing::MeshPlan { pml: femcore::meshing::PmlMode::Structured, ..femcore::meshing::MeshPlan::adaptive(0.6, 0.7) },
    ];
    for plan in plans {
        let g = femcore::meshing::mesh_model(&m.list, &plan).unwrap();
        assert_eq!(g.labels, simple.labels, "{plan:?}");
        let domains: std::collections::BTreeSet<i32> = g.mesh.cell_labels.iter().copied().collect();
        assert_eq!(domains, g.labels.domains.iter().map(|d| d.id).collect(), "every region has cells: {plan:?}");
        let faces: std::collections::BTreeSet<i32> = g.mesh.facet_labels.iter().copied().collect();
        let sphere = g.boundary("air|ball");
        assert_eq!(sphere.len(), 1);
        assert!(faces.contains(&sphere[0]), "the sphere carries its own facets: {plan:?}");
        assert!(faces.iter().all(|f| g.labels.boundaries.iter().any(|b| b.id == *f)), "no other number appears: {plan:?}");
        assert_eq!(g.mesh.dim, 3);
    }
}

#[test]
fn mesh_geom3d_structured_pml_keeps_the_physical_vertex_and_cell_numbers() {
    // The 3D counterpart of `mesh_pml_generator_keeps_the_physical_vertex_and_cell_numbers` for a geom model (ADR 52): the physical window
    // is meshed alone, the layers are appended, and the uncut mesh and its tetrahedral cut both start with the physical vertices and cells.
    use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts};
    use femcore::meshing::geomdrive3::{mesh_geom_3d_with, Sizes3};
    use femcore::meshing::pml_geom3::{structured_pml_3d, structured_pml_3d_native, Native3};
    let m = ball_in_box_with_pml();
    let ev = evaluate(&m.list).unwrap();
    let mut skip = vec![false; ev.regions().len()];
    for p in ev.pml_regions() {
        skip[p.region as usize] = true;
    }
    let (phys, _, _) = mesh_geom_3d_with(ev, &GeomMeshOpts::new(0.5), &Sizes3 { region_h: vec![], skip }).unwrap();
    let (np, nc) = (phys.mesh.n_points(), phys.mesh.n_cells());
    let Native3::Cells(cells, report) = structured_pml_3d_native(phys.clone()).unwrap() else { panic!("the model has a Pml") };
    assert_eq!((report.physical_cells, report.dropped_cells), (nc, 0));
    assert_eq!(&cells.points[..np].iter().flatten().copied().collect::<Vec<f64>>()[..], &phys.mesh.points[..], "physical coordinates by id");
    for c in 0..nc {
        assert_eq!(cells.cell(c), phys.mesh.cell(c), "cell {c}");
        assert_eq!(cells.label[c], phys.mesh.cell_labels[c], "label of cell {c}");
    }
    assert!(cells.n_cells() > nc && cells.n_points() > np);
    let (split, _) = structured_pml_3d(phys.clone()).unwrap();
    assert_eq!(&split.mesh.points[..3 * np], &phys.mesh.points[..], "the solver view keeps the physical vertices");
    assert_eq!(&split.mesh.cells[..4 * nc], &phys.mesh.cells[..], "and the physical cells, first and unchanged");
}

#[test]
fn mesh_request_and_report_match_their_exemplars_and_schemas() {
    // ADR 54: the serialisable request and report of a mesh job, in 2D and 3D. The request exemplar is compared as a value; the report
    // exemplar by structure (the counts depend on the mesher); both validate against the checked-in schemas, which are the generated ones.
    use femcore::meshing::request::{
        mesh_request, report_schema, request_schema, MeshReport, MeshRequest, MesherRequest,
        PmlModeRequest, WavesRequest,
    };
    use femcore::meshing::Sizing;
    check_exemplar_value("mesh_request.schema.v1.json", &request_schema());
    check_exemplar_value("mesh_report.schema.v1.json", &report_schema());
    let req = MeshRequest {
        h: 0.6,
        max_dev: None,
        mesher: MesherRequest::Adaptive { scale: 0.7, grade: 0.4, max_cells: Some(200_000) },
        pml_mode: PmlModeRequest::Structured,
        structured_pml: None,
        curved: false,
        gap_cells: None,
        z_coarsen: false,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: false,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
        waves: Some(WavesRequest {
            sizing: Sizing { order: 2, headroom: 1.0, wavelengths: [("air".to_string(), 4.0), ("ball".to_string(), 1.2)].into_iter().collect() },
            points_per_wavelength: None,
            pml_relax: 1.5,
            eigen_headroom: 1.5,
        }),
    };
    let rv = serde_json::to_value(&req).unwrap();
    check_exemplar_value("mesh_request.v1.json", &rv);
    validate(&request_schema(), &rv, "request").unwrap();
    assert_eq!(serde_json::from_value::<MeshRequest>(rv).unwrap(), req, "the request round-trips");
    let minimal: MeshRequest = serde_json::from_value(json!({ "h": 0.1, "mesher": { "kind": "simple" } })).unwrap();
    assert!(
        minimal.pml_mode == PmlModeRequest::Auto
            && minimal.resolved_pml_mode() == PmlModeRequest::Auto
            && minimal.waves.is_none()
            && matches!(minimal.mesher, MesherRequest::Simple)
    );
    let legacy_true: MeshRequest =
        serde_json::from_value(json!({ "h": 0.1, "mesher": { "kind": "simple" }, "structured_pml": true })).unwrap();
    assert_eq!(legacy_true.resolved_pml_mode(), PmlModeRequest::Structured);
    let legacy_false: MeshRequest =
        serde_json::from_value(json!({ "h": 0.1, "mesher": { "kind": "simple" }, "structured_pml": false })).unwrap();
    assert_eq!(legacy_false.resolved_pml_mode(), PmlModeRequest::Inline);
    let (_, report) = mesh_request(&ball_in_box_with_pml().list, &req).unwrap();
    let pv = serde_json::to_value(&report).unwrap();
    validate(&report_schema(), &pv, "report").unwrap();
    check_exemplar_value("mesh_report.v1.json", &shape(&pv));
    assert_eq!(serde_json::from_value::<MeshReport>(pv).unwrap(), report);
    assert_eq!(report.dim, 3);
    assert_eq!(report.structured_layers, Some(vec![[2, 2], [2, 2], [2, 2]]));
    assert!(report.adaptation.is_some() && report.domains.iter().any(|d| d.wave_size.is_some()));
    // the schema refuses what the request type refuses
    assert!(validate(&request_schema(), &json!({ "h": -1.0, "mesher": { "kind": "simple" } }), "r").is_err());
    assert!(validate(&request_schema(), &json!({ "h": 1.0, "mesher": { "kind": "other" } }), "r").is_err());
    assert!(validate(&request_schema(), &json!({ "h": 1.0, "mesher": { "kind": "simple" }, "extra": 1 }), "r").is_err());
}

#[test]
fn adapt_maxwell_estimator_contract() {
    // ADR 53: `CurlField` is all the estimator knows of a solution; the four parts add up to the indicators, boundary conditions decide which
    // faces count, and a skipped cell has no indicator.
    use femcore::adapt::maxwell::{maxwell_indicators, CurlField, MaxwellOptions};
    use femcore::linalg::csr::C;
    struct Linear;
    impl CurlField for Linear {
        fn sample(&self, _cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
            // E = lam_1 e_x: a field with a jump of its normal component across the faces, as every cell has its own lam_1
            ([C::new(lam[1], 0.0), C::new(0.0, 0.0), C::new(0.0, 0.0)], [C::new(0.0, 0.0), C::new(1.0, 0.0), C::new(0.0, 0.0)])
        }
    }
    let m = femcore::testkit::block(1.0, 1.0, 1.0, 2, 2, 2);
    let (eps, mu) = (vec![1.0; m.n_cells()], vec![1.0; m.n_cells()]);
    let opt = |natural: &'static (dyn Fn(i32) -> bool + Sync), skip: &'static (dyn Fn(i32) -> bool + Sync)| MaxwellOptions { k0_sq: 1.0, eps_r: &eps, mu_r: &mu, skip, natural };
    let ind = maxwell_indicators(&m, &Linear, &opt(&|_| false, &|_| false)).unwrap();
    let sum: f64 = ind.eta2.iter().sum();
    assert!((ind.parts.iter().sum::<f64>() - sum).abs() <= 1e-12 * sum);
    assert!((ind.total * ind.total - sum).abs() <= 1e-12 * sum);
    // natural boundaries add a term, perfect conductors none
    let nat = maxwell_indicators(&m, &Linear, &opt(&|_| true, &|_| false)).unwrap();
    assert!(nat.parts[2] > ind.parts[2]);
    // a skipped domain has no indicator
    let label = m.cell_labels[0];
    let all = |l: i32| l == label;
    let skipped = maxwell_indicators(&m, &Linear, &MaxwellOptions { k0_sq: 1.0, eps_r: &eps, mu_r: &mu, skip: &all, natural: &|_| false }).unwrap();
    assert!(skipped.eta2.iter().all(|&e| e == 0.0));
    // wrong sizes are refused
    assert!(maxwell_indicators(&m, &Linear, &MaxwellOptions { k0_sq: 1.0, eps_r: &eps[1..], mu_r: &mu, skip: &|_| false, natural: &|_| false }).is_err());
}
