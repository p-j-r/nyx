// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified `bolero` fuzz and regression tests for `nyx_storage::NyxReader` parsing and `ReadBudgetConfig` / `max_solutions` caps.

use bolero::TypeGenerator;
use nyx_storage::{
    fuzz_walk_untrusted, write_nyx_packed, write_nyx_v3, FieldStorageLayout, MeshData, ModeHeader, NyxReader,
    NyxSolutionData, NyxV3Summary, PackedComponentMeta, ReadBudgetConfig,
};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

static TEMP_SEQ: AtomicUsize = AtomicUsize::new(0);

fn scratch_file(tag: &str) -> PathBuf {
    let n = TEMP_SEQ.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("nyx_fuzz_storage_{}_{n}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("create temp dir");
    dir.join(format!("{tag}.nyx"))
}

fn tight_budget() -> ReadBudgetConfig {
    ReadBudgetConfig {
        max_dataset_bytes: 256 * 1024,
        max_total_result_bytes: 512 * 1024,
        max_group_depth: 12,
        max_solutions: 16,
    }
}

fn sample_mesh() -> MeshData {
    MeshData {
        cell_topology: "tetrahedron".into(),
        element_family: "Nedelec".into(),
        element_degree: 1,
        space_dimension: 3,
        cell_labels: vec![1, 1],
        cells_shape: [2, 4],
        cells: vec![0, 1, 2, 3, 1, 2, 3, 4],
        facet_labels: vec![10, 10],
        facets_shape: [2, 3],
        facets: vec![0, 1, 2, 1, 3, 4],
        points_shape: [5, 3],
        points: vec![
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
        ],
        edge_node_pairs: Vec::new(),
        edge_node_points: Vec::new(),
    }
}

fn build_seed_bytes(layout: FieldStorageLayout) -> Vec<u8> {
    let path = scratch_file(layout.as_str());
    let mesh = sample_mesh();
    let headers = vec![
        ModeHeader { mode_index: 0, eig_re: 1.5, eig_im: 0.01, label: "mode_0".into() },
        ModeHeader { mode_index: 1, eig_re: 2.1, eig_im: 0.02, label: "mode_1".into() },
        ModeHeader { mode_index: 2, eig_re: 2.8, eig_im: 0.03, label: "mode_2".into() },
    ];
    let components = vec![PackedComponentMeta {
        name: "Ex".into(),
        unit: "V/m".into(),
        nloc: 4,
        has_im: true,
    }];
    let n_vals = headers.len() * components.len() * mesh.cells_shape[0] * mesh.cells_shape[1] * 2;
    let packed = (0..n_vals).map(|i| (i as f64) * 0.125).collect::<Vec<_>>();
    write_nyx_packed(
        &path,
        b"{\"format\":1,\"physics\":\"ems\"}",
        &mesh,
        "ems",
        "eigen",
        Some(1.55e-6),
        &headers,
        &components,
        &packed,
        layout,
    )
    .expect("seed packed file");
    let bytes = std::fs::read(&path).expect("read seed bytes");
    let _ = std::fs::remove_dir_all(path.parent().expect("parent"));
    bytes
}

fn build_v3_seed_bytes() -> Vec<u8> {
    let path = scratch_file("v3_summary");
    let mesh = sample_mesh();
    let solutions = (0..3usize)
        .map(|k| {
            let mut cell_fields = BTreeMap::new();
            cell_fields.insert(
                "00_Ex".into(),
                ("Ex".into(), "V/m".into(), 4, vec![k as f64; 8], vec![0.5; 8]),
            );
            let mut fields = BTreeMap::new();
            fields.insert("00_V".into(), ("V".into(), "V".into(), "[]".into(), vec![k as f64; 5]));
            let mut viz_fields = BTreeMap::new();
            viz_fields.insert("V".into(), vec![k as f64; 5]);
            NyxSolutionData {
                eig_re: 1.5 + 0.5 * (k as f64),
                eig_im: 0.01 * ((k + 1) as f64),
                label: format!("mode_{k}"),
                cell_fields,
                fields_order_json: "[\"V\"]".into(),
                fields,
                viz_fields,
            }
        })
        .collect();
    let summary = NyxV3Summary {
        format: "nyx-results".into(),
        version: 3,
        writer: "serial".into(),
        processes: 1,
        project_bytes: b"{\"format\":1,\"physics\":\"ems\"}".to_vec(),
        mesh_cell_labels: mesh.cell_labels.clone(),
        mesh_cells_shape: mesh.cells_shape,
        mesh_cells: mesh.cells.clone(),
        mesh_facet_labels: mesh.facet_labels.clone(),
        mesh_facets_shape: mesh.facets_shape,
        mesh_facets: mesh.facets.clone(),
        mesh_points_shape: mesh.points_shape,
        mesh_points: mesh.points.clone(),
        physics: "ems".into(),
        study: "eigen".into(),
        dim: 3,
        dofs: 5,
        n_solutions: 3,
        order: 1,
        seconds: 0.1,
        log_bytes: b"[]".to_vec(),
        scalars_bytes: b"[]".to_vec(),
        viz_topology: "tetrahedron".into(),
        viz_cells_shape: mesh.cells_shape,
        viz_cells: mesh.cells.clone(),
        viz_points_shape: mesh.points_shape,
        viz_points: mesh.points.clone(),
        solutions,
    };
    write_nyx_v3(&summary, &path).expect("seed v3 file");
    let bytes = std::fs::read(&path).expect("read v3 bytes");
    let _ = std::fs::remove_dir_all(path.parent().expect("parent"));
    bytes
}

fn seed_corpus() -> &'static [Vec<u8>; 4] {
    static SEEDS: OnceLock<[Vec<u8>; 4]> = OnceLock::new();
    SEEDS.get_or_init(|| {
        [
            build_v3_seed_bytes(),
            build_seed_bytes(FieldStorageLayout::L1LegacyGzip3),
            build_seed_bytes(FieldStorageLayout::L2ShuffleLzfPacked),
            build_seed_bytes(FieldStorageLayout::L3ContiguousPacked),
        ]
    })
}

fn exercise_reader(bytes: &[u8], budget: ReadBudgetConfig, mode_idx: usize) {
    let _ = fuzz_walk_untrusted(bytes, budget);
    if let Ok(reader) = NyxReader::from_bytes(bytes.to_vec(), budget) {
        let _ = reader.read_mesh();
        let _ = reader.scan_mode_headers("ems");
        let _ = reader.read_mode_cell_fields("ems", mode_idx);
        let _ = reader.read_nyx_v3_summary_for_physics("ems", true);
        let _ = reader.read_nyx_v3_summary_for_physics("sms", false);
    }
}

#[derive(Debug, TypeGenerator)]
struct MutationOp {
    layout_idx: u8,
    offset: u16,
    xor_mask: u8,
    patch_bytes: [u8; 8],
    use_patch: bool,
    truncate_to: u16,
    do_truncate: bool,
    mode_idx: u8,
    max_solutions: u8,
    max_dataset_kb: u8,
    max_depth: u8,
}

#[test]
fn fuzz_nyx_storage_arbitrary_bytes_no_panic() {
    let budget = tight_budget();
    bolero::check!().with_iterations(512).for_each(|bytes: &[u8]| {
        exercise_reader(bytes, budget, 0);
    });
}

#[test]
fn fuzz_nyx_storage_mutated_hdf5_images_and_budgets_no_panic() {
    let corpus = seed_corpus();
    let mut work = Vec::new();
    bolero::check!()
        .with_type::<MutationOp>()
        .with_iterations(512)
        .for_each(|op: &MutationOp| {
            let base = &corpus[(op.layout_idx as usize) % corpus.len()];
            work.clear();
            work.extend_from_slice(base);
            if work.is_empty() {
                return;
            }
            let pos = (op.offset as usize) % work.len();
            if op.use_patch {
                let end = (pos + op.patch_bytes.len()).min(work.len());
                let len = end - pos;
                work[pos..end].copy_from_slice(&op.patch_bytes[..len]);
            } else {
                work[pos] ^= op.xor_mask;
            }
            if op.do_truncate {
                let len = (op.truncate_to as usize) % (work.len() + 1);
                work.truncate(len);
            }
            let budget = ReadBudgetConfig {
                max_dataset_bytes: u64::from(op.max_dataset_kb.max(1)) * 1024,
                max_total_result_bytes: u64::from(op.max_dataset_kb.max(1)) * 2048,
                max_group_depth: usize::from(op.max_depth % 16),
                max_solutions: usize::from(op.max_solutions % 8),
            };
            exercise_reader(&work, budget, usize::from(op.mode_idx % 6));
        });
}

#[test]
fn regression_read_budget_and_max_solutions_strict_enforcement() {
    let full = ReadBudgetConfig::default();

    // Verify v3 summary + max_solutions / dataset / total_result / depth budget enforcement.
    let v3_bytes = &seed_corpus()[0];
    let v3_reader = NyxReader::from_bytes(v3_bytes.clone(), full).expect("open valid v3 seed");
    assert_eq!(v3_reader.read_mesh().expect("v3 mesh").cells_shape, [2, 4]);
    let summary = v3_reader.read_nyx_v3_summary_for_physics("ems", true).expect("v3 summary");
    assert_eq!(summary.n_solutions, 3);
    assert_eq!(summary.solutions.len(), 3);

    let v3_capped = NyxReader::from_bytes(v3_bytes.clone(), ReadBudgetConfig { max_solutions: 2, ..tight_budget() })
        .expect("open v3 header");
    assert!(
        v3_capped.read_nyx_v3_summary_for_physics("ems", false).is_err(),
        "max_solutions = 2 must reject a 3-solution v3 file"
    );

    let v3_tiny_total = NyxReader::from_bytes(
        v3_bytes.clone(),
        ReadBudgetConfig { max_total_result_bytes: 16, ..tight_budget() },
    )
    .expect("open v3 header");
    assert!(
        v3_tiny_total.read_nyx_v3_summary_for_physics("ems", false).is_err(),
        "max_total_result_bytes = 16 must reject v3 summary"
    );

    // Verify packed L1/L2/L3 + budget enforcement.
    for bytes in &seed_corpus()[1..] {
        let reader = NyxReader::from_bytes(bytes.clone(), full).expect("open valid packed seed");
        assert_eq!(reader.read_mesh().expect("mesh").cells_shape, [2, 4]);
        let headers = reader.scan_mode_headers("ems").expect("headers");
        assert_eq!(headers.len(), 3);
        let fields = reader.read_mode_cell_fields("ems", 0).expect("mode 0 fields");
        assert!(!fields.is_empty());

        let capped_sol = ReadBudgetConfig { max_solutions: 2, ..tight_budget() };
        let capped_reader = NyxReader::from_bytes(bytes.clone(), capped_sol).expect("open packed header");
        assert!(
            capped_reader.scan_mode_headers("ems").is_err(),
            "max_solutions = 2 must reject a 3-mode packed file"
        );
    }

    for bytes in seed_corpus() {
        let tiny_ds = ReadBudgetConfig { max_dataset_bytes: 8, ..tight_budget() };
        let tiny_reader = NyxReader::from_bytes(bytes.clone(), tiny_ds).expect("open header");
        assert!(tiny_reader.read_mesh().is_err(), "max_dataset_bytes = 8 must reject mesh");

        let shallow = ReadBudgetConfig { max_group_depth: 0, ..tight_budget() };
        assert!(
            fuzz_walk_untrusted(bytes, shallow).is_err(),
            "max_group_depth = 0 must reject nested groups"
        );
    }
}
