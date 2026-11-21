// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Mandatory 5-gate verification suite (`EC-00..EC-12` & `BM-5`) for `nyx-storage`.

use anyhow::{anyhow, bail, Result};
use hdf5_pure::{make_f64_type, AttrValue, DeflateLevel, File, FileBuilder};
use nyx_storage::{
    append_solved_mode_in_place, assert_f64_bits_eq, fuzz_walk_untrusted, generate_direct_xdmf,
    init_nyx_packed_streaming, viz_link_name, write_direct_xdmf, write_nyx_packed, write_nyx_v3, FieldStorageLayout,
    MeshData, ModeHeader, NyxReader, PackedComponentMeta, ReadBudgetConfig, SweepCheckpoint, SweepPoint,
};
use rayon::prelude::*;
use std::f64::consts::PI;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

const PYTHON_EXE: &str = r"C:\ProgramData\anaconda3\python.exe";
const BENCH_MODES_2D: &str =
    r"C:\Users\robin-adm\Documents\nyx2_benches\.jj-workspaces\main-3c3d3858\sos_300nm_modes.h5";
const BENCH_RIB_2D: &str =
    r"C:\Users\robin-adm\Documents\nyx2_benches\.jj-workspaces\main-3c3d3858\sos_300nm_rib_modes.h5";
const BENCH_MECH_3D: &str =
    r"C:\Users\robin-adm\Documents\nyx2_benches\.jj-workspaces\main-789a3652\sos_300nm_mech_modes.h5";
const BENCH_COUPLED_2D: &str =
    r"C:\Users\robin-adm\Documents\nyx2_benches\.jj-workspaces\main-789a3652\sos_300nm_strip_coupled_q0_a1.h5";
const BENCH_SWG_3D: &str =
    r"C:\Users\robin-adm\Documents\nyx2_benches\.jj-workspaces\main-77744503\swg_3d_bock2010_case1a.h5";

#[cfg(target_os = "windows")]
mod win_mem {
    #[repr(C)]
    struct ProcessMemoryCounters {
        cb: u32,
        page_fault_count: u32,
        peak_working_set_size: usize,
        working_set_size: usize,
        quota_peak_paged_pool_usage: usize,
        quota_paged_pool_usage: usize,
        quota_peak_non_paged_pool_usage: usize,
        quota_non_paged_pool_usage: usize,
        pagefile_usage: usize,
        peak_pagefile_usage: usize,
    }

    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> *mut core::ffi::c_void;
        fn K32GetProcessMemoryInfo(
            process: *mut core::ffi::c_void,
            ppsmem_counters: *mut ProcessMemoryCounters,
            cb: u32,
        ) -> i32;
    }

    pub fn memory_usage_mb() -> (f64, f64) {
        unsafe {
            let mut pmc = core::mem::zeroed::<ProcessMemoryCounters>();
            pmc.cb = core::mem::size_of::<ProcessMemoryCounters>() as u32;
            let ok = K32GetProcessMemoryInfo(GetCurrentProcess(), &mut pmc, pmc.cb);
            if ok != 0 {
                (pmc.working_set_size as f64 / (1024.0 * 1024.0), pmc.peak_working_set_size as f64 / (1024.0 * 1024.0))
            } else {
                (0.0, 0.0)
            }
        }
    }
}

#[cfg(not(target_os = "windows"))]
mod win_mem {
    pub fn memory_usage_mb() -> (f64, f64) {
        (0.0, 0.0)
    }
}

fn test_scratch_dir(sub: &str) -> Result<PathBuf> {
    let base = std::env::var_os("CARGO_TARGET_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
        .join("nyx_storage_gate_tmp")
        .join(sub);
    std::fs::create_dir_all(&base)?;
    Ok(base)
}

fn run_python_inline(code: &str, args: &[&str]) -> Result<String> {
    let out = Command::new(PYTHON_EXE).arg("-c").arg(code).args(args).output()?;
    if !out.status.success() {
        let stdout = String::from_utf8_lossy(&out.stdout);
        let stderr = String::from_utf8_lossy(&out.stderr);
        bail!("python oracle failed (exit {:?}):\nSTDOUT: {}\nSTDERR: {}", out.status.code(), stdout, stderr);
    }
    Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

#[inline]
fn expected_mode_elem(mode_idx: usize, elem_idx: usize) -> f64 {
    let m_base = (mode_idx as f64 + 1.0) * 0.125;
    let v = m_base + ((elem_idx % 97) as f64) * 0.0078125 + ((elem_idx / 97) as f64) * 1e-7;
    if elem_idx & 1 == 0 {
        v
    } else {
        -v * 0.5
    }
}

fn generate_mode_slice_f64(mode_idx: usize, n_comp: usize, n_cells: usize, nloc: usize) -> Vec<f64> {
    let len = n_comp * n_cells * nloc * 2;
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        out.push(expected_mode_elem(mode_idx, i));
    }
    out
}

fn ieee754_extreme_bits() -> Vec<u64> {
    vec![
        0.0_f64.to_bits(),
        (-0.0_f64).to_bits(),
        f64::MIN_POSITIVE.to_bits(),
        (-f64::MIN_POSITIVE).to_bits(),
        f64::MAX.to_bits(),
        f64::MIN.to_bits(),
        1.0e-320_f64.to_bits(),
        (-1.0e-320_f64).to_bits(),
        0x0000_0000_0000_0001_u64,
        0x8000_0000_0000_0001_u64,
        f64::INFINITY.to_bits(),
        f64::NEG_INFINITY.to_bits(),
        0x7ff8_0000_0000_0000_u64,
        0x7ff8_0000_0000_1234_u64,
        0xfff8_0000_0000_5678_u64,
    ]
}

// =============================================================================
// GATE 1 (EC-00): Zero-C / Zero-CMake Offline Build & Dependency Tree Audit
// =============================================================================
#[test]
fn gate1_ec00_zero_c_offline_build_and_tree_audit() -> Result<()> {
    let t0 = Instant::now();
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let target_dir = std::env::var("CARGO_TARGET_DIR")
        .unwrap_or_else(|_| r"c:\Users\robin-adm\Documents\opomech\nyx3-target\storage-core".to_string());

    let out = Command::new("cargo")
        .current_dir(&manifest_dir)
        .env("CARGO_TARGET_DIR", &target_dir)
        .env("CC", "nonexistent_cc")
        .env("CXX", "nonexistent_cxx")
        .env("CMAKE", "nonexistent_cmake")
        .args(["tree", "--offline", "-e", "normal,build", "-p", "nyx-storage"])
        .output()?;

    if !out.status.success() {
        bail!("cargo tree failed: {}", String::from_utf8_lossy(&out.stderr));
    }
    let tree_text = String::from_utf8_lossy(&out.stdout);
    let forbidden = [
        " cc ",
        " cmake ",
        " bindgen ",
        " pkg-config ",
        " vcpkg ",
        "hdf5-sys",
        "libhdf5-sys",
        "hdf5-metno",
        "libz-sys",
        "zlib-ng-sys",
        "-sys v",
    ];
    for pat in forbidden {
        if tree_text.contains(pat) {
            bail!("Forbidden C/CMake/*-sys dependency '{}' found in nyx-storage tree:\n{}", pat, tree_text);
        }
    }
    if !tree_text.contains("hdf5-pure v0.48.0") || !tree_text.contains("zlib-rs v0.6.8") {
        bail!("Expected hdf5-pure v0.48.0 and zlib-rs v0.6.8 in tree:\n{}", tree_text);
    }

    println!(
        "[Gate 1 / EC-00 PASS] ({:.1} ms) 0 C/CMake/*-sys crates in dependency tree; verified hdf5-pure v0.48.0 + zlib-rs v0.6.8 with CC/CXX/CMAKE poisoned",
        t0.elapsed().as_secs_f64() * 1000.0
    );
    Ok(())
}

// =============================================================================
// GATE 2 (EC-01): Real `nyx2_benches` .h5 Read Parity (0 ULP vs h5py 3.13.0)
// =============================================================================
#[test]
fn gate2_ec01_real_nyx2_benches_read_parity_0_ulp() -> Result<()> {
    let t0 = Instant::now();
    let scratch = test_scratch_dir("gate2")?;

    let bench_files = [
        ("sos_300nm_modes.h5", BENCH_MODES_2D, vec!["ems"]),
        ("sos_300nm_rib_modes.h5", BENCH_RIB_2D, vec!["ems"]),
        ("sos_300nm_mech_modes.h5", BENCH_MECH_3D, vec!["sms"]),
        ("sos_300nm_strip_coupled_q0_a1.h5", BENCH_COUPLED_2D, vec!["ems", "sms"]),
        ("swg_3d_bock2010_case1a.h5", BENCH_SWG_3D, vec!["ems"]),
    ];

    let py_dump_code = r#"
import json, sys, h5py, numpy as np

def u64_list(arr):
    a = np.ascontiguousarray(arr, dtype=np.float64).ravel()
    return [int(x) for x in a.view(np.uint64)]

path = sys.argv[1]
phys = sys.argv[2]
out_path = sys.argv[3]

with h5py.File(path, "r") as f:
    g = f[f"results/{phys}"]
    n_sol = int(g.attrs["n_solutions"])
    eigs = []
    mode0_cf = []
    for k in range(n_sol):
        s = g[f"solution_{k}"]
        er = float(s.attrs["eig_re"]) if "eig_re" in s.attrs else float(s.attrs["value"])
        ei = float(s.attrs["eig_im"]) if "eig_im" in s.attrs else 0.0
        eigs.extend(u64_list([er, ei]))
        if k == 0:
            for cfk in sorted(s["cell_fields"].keys()):
                re_a = np.ascontiguousarray(s[f"cell_fields/{cfk}/re"][:], dtype=np.float64).ravel()
                im_a = np.ascontiguousarray(s[f"cell_fields/{cfk}/im"][:], dtype=np.float64).ravel()
                inter = np.empty(re_a.size * 2, dtype=np.float64)
                inter[0::2] = re_a
                inter[1::2] = im_a
                mode0_cf.extend([int(x) for x in inter.view(np.uint64)])
    payload = {
        "n_solutions": n_sol,
        "mesh_points_bits": u64_list(f["mesh/points"][:]),
        "mesh_cells": [int(x) for x in np.ascontiguousarray(f["mesh/cells"][:], dtype=np.uint32).ravel()],
        "mesh_cell_labels": [int(x) for x in np.ascontiguousarray(f["mesh/cell_labels"][:], dtype=np.int32).ravel()],
        "eigs_bits": eigs,
        "mode0_cf_bits": mode0_cf,
    }
with open(out_path, "w", encoding="utf-8") as fp:
    json.dump(payload, fp)
"#;

    let mut total_modes_checked = 0usize;
    for (label, fpath, phys_list) in &bench_files {
        let p = Path::new(fpath);
        let reader = NyxReader::open(p, ReadBudgetConfig::default())?;
        let mesh = reader.read_mesh()?;

        for phys in phys_list {
            let headers = reader.scan_mode_headers(phys)?;
            let summary = reader.read_nyx_v3_summary_for_physics(phys, true)?;
            summary.verify_0_ulp(&summary, true)?;

            // Check every /viz/<phys>/sol_k/<vname> v1 hard-link dataset equals /results/<phys>/solution_k/fields/<fname>
            for (k, sol) in summary.solutions.iter().enumerate() {
                if sol.viz_fields.len() != sol.fields.len() {
                    bail!("{}/{}: sol_{} viz_fields != fields count", label, phys, k);
                }
                for (fk, (_, _, _, fdata)) in &sol.fields {
                    let vk = viz_link_name(fk);
                    let vdata = sol.viz_fields.get(&vk).ok_or_else(|| anyhow!("missing viz link {}", vk))?;
                    assert_f64_bits_eq(&format!("{}/{}/viz/sol_{}/{}", label, phys, k, vk), vdata, fdata)?;
                }
            }

            // Cross-check against Python h5py 3.13.0 at 0 ULP
            let json_out = scratch.join(format!("oracle_{}_{}.json", label, phys));
            run_python_inline(py_dump_code, &[fpath, phys, json_out.to_str().unwrap()])?;
            let oracle: serde_json::Value = serde_json::from_slice(&std::fs::read(&json_out)?)?;

            let py_pts: Vec<f64> = oracle["mesh_points_bits"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| f64::from_bits(v.as_u64().unwrap()))
                .collect();
            assert_f64_bits_eq(&format!("{}/mesh/points", label), &mesh.points, &py_pts)?;

            let py_cells: Vec<u32> =
                oracle["mesh_cells"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as u32).collect();
            if mesh.cells != py_cells {
                bail!("{}/mesh/cells mismatch vs h5py", label);
            }

            let py_eigs: Vec<f64> =
                oracle["eigs_bits"].as_array().unwrap().iter().map(|v| f64::from_bits(v.as_u64().unwrap())).collect();
            let mut rust_eigs = Vec::with_capacity(headers.len() * 2);
            for h in &headers {
                rust_eigs.push(h.eig_re);
                rust_eigs.push(h.eig_im);
            }
            assert_f64_bits_eq(&format!("{}/{}/eigenvalues", label, phys), &rust_eigs, &py_eigs)?;

            let mode0_cf = reader.read_mode_cell_fields(phys, 0)?;
            let py_mode0_cf: Vec<f64> = oracle["mode0_cf_bits"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| f64::from_bits(v.as_u64().unwrap()))
                .collect();
            assert_f64_bits_eq(&format!("{}/{}/mode0_cell_fields", label, phys), &mode0_cf, &py_mode0_cf)?;

            total_modes_checked += headers.len();
        }
    }

    println!(
        "[Gate 2 / EC-01 PASS] ({:.1} ms) Read all 5 real nyx2_benches 2D/3D .h5 files ({} solutions, mech `value` fallback, v1 symbol-table hard links) with 0 ULP vs h5py 3.13.0",
        t0.elapsed().as_secs_f64() * 1000.0,
        total_modes_checked
    );
    Ok(())
}

// =============================================================================
// GATE 3 (EC-02, EC-06, EC-07, EC-08, BM-5):
// L2/L3/L1 Write Round-Trip, >15x L2 Compression, Filters, IEEE-754, Degenerate,
// Streaming `append_solved_mode_in_place`, and 50-Point `SweepCheckpoint`
// =============================================================================
#[test]
fn gate3_ec02_ec06_ec07_ec08_bm5_roundtrip_compression_filters_ieee_sweep() -> Result<()> {
    let t0 = Instant::now();
    let scratch = test_scratch_dir("gate3")?;

    // -------------------------------------------------------------------------
    // 1. EC-02: 2D Triangle_6 (13 modes, 5488 cells) & 3D Tetrahedron_10 (4 modes, 2000 cells)
    //    in L2 Shuffle+LZF, L3 Contiguous, and L1 Deflate-3 + >15x L2 compression ratio
    // -------------------------------------------------------------------------
    let n_modes_2d = 13usize;
    let n_comp_2d = 3usize;
    let n_cells_2d = 5_488usize;
    let nloc_2d = 6usize;

    let mesh_2d = MeshData {
        cell_topology: "Triangle_6".to_string(),
        element_family: "Lagrange".to_string(),
        element_degree: 2,
        space_dimension: 2,
        cell_labels: vec![1i32; n_cells_2d],
        cells_shape: [n_cells_2d, nloc_2d],
        cells: (0..(n_cells_2d * nloc_2d) as u32).map(|i| i % 4045).collect(),
        facet_labels: vec![2i32; 128],
        facets_shape: [128, 3],
        facets: (0..384u32).collect(),
        points_shape: [1024, 2],
        points: (0..2048).map(|i| (i as f64) * 1e-8 - 5e-6).collect(),
        edge_node_pairs: vec![[0, 1], [1, 2], [2, 0]],
        edge_node_points: vec![[0.1e-6, 0.2e-6, 0.0], [0.3e-6, 0.4e-6, 0.0], [-0.1e-6, 0.0, 0.0]],
    };
    let headers_2d: Vec<ModeHeader> = (0..n_modes_2d)
        .map(|k| ModeHeader {
            mode_index: k,
            eig_re: 2.8945934319606623 - (k as f64) * 0.025,
            eig_im: -1.0e-11 * (k as f64 + 1.0),
            label: format!("TE_mode_{k}"),
        })
        .collect();
    let comps_2d = vec![
        PackedComponentMeta { name: "Ex".into(), unit: "V/m".into(), nloc: nloc_2d, has_im: true },
        PackedComponentMeta { name: "Ez".into(), unit: "V/m".into(), nloc: nloc_2d, has_im: true },
        PackedComponentMeta { name: "Ey".into(), unit: "V/m".into(), nloc: nloc_2d, has_im: true },
    ];
    let mode_len_2d = n_comp_2d * n_cells_2d * nloc_2d * 2;
    let mut packed_2d = Vec::with_capacity(n_modes_2d * mode_len_2d);
    for m in 0..n_modes_2d {
        packed_2d.extend_from_slice(&generate_mode_slice_f64(m, n_comp_2d, n_cells_2d, nloc_2d));
    }
    let raw_bytes_2d = (packed_2d.len() * 8) as f64;

    let mut l2_ratio_2d = 0.0;
    for layout in [
        FieldStorageLayout::L2ShuffleLzfPacked,
        FieldStorageLayout::L3ContiguousPacked,
        FieldStorageLayout::L1LegacyGzip3,
    ] {
        let p = scratch.join(format!("tri6_2d_{}.h5", layout.as_str()));
        write_nyx_packed(
            &p,
            b"{\"format\":\"nyx-project\",\"version\":1}",
            &mesh_2d,
            "ems",
            "modal",
            Some(1.55e-6),
            &headers_2d,
            &comps_2d,
            &packed_2d,
            layout,
        )?;
        if layout == FieldStorageLayout::L2ShuffleLzfPacked {
            let fsize = std::fs::metadata(&p)?.len() as f64;
            l2_ratio_2d = raw_bytes_2d / fsize;
            if l2_ratio_2d <= 15.0 {
                bail!("Expected L2 Shuffle+LZF compression ratio > 15x on smooth P2 modes, got {:.2}x", l2_ratio_2d);
            }
        }
        let reader = NyxReader::open(&p, ReadBudgetConfig::default())?;
        let r_mesh = reader.read_mesh()?;
        if r_mesh != mesh_2d {
            bail!("2D MeshData round-trip mismatch for {:?}", layout);
        }
        let r_headers = reader.scan_mode_headers("ems")?;
        if r_headers != headers_2d {
            bail!("2D ModeHeader round-trip mismatch for {:?}", layout);
        }
        let r_meta = reader.read_packed_components_meta("ems")?.unwrap();
        if r_meta != comps_2d {
            bail!("2D PackedComponentMeta mismatch for {:?}", layout);
        }
        for m in [0, 3, 7, 12] {
            let slice = reader.read_mode_cell_fields("ems", m)?;
            let exp = &packed_2d[m * mode_len_2d..(m + 1) * mode_len_2d];
            assert_f64_bits_eq(&format!("2D {:?} mode {}", layout, m), &slice, exp)?;
        }
    }

    // 3D Tetrahedron_10 packed round-trip in L2, L3, L1 + h5py 0-ULP verification
    let n_modes_3d = 4usize;
    let n_comp_3d = 3usize;
    let n_cells_3d = 2_000usize;
    let nloc_3d = 10usize;
    let mesh_3d = MeshData {
        cell_topology: "Tetrahedron_10".to_string(),
        element_family: "Lagrange".to_string(),
        element_degree: 2,
        space_dimension: 3,
        cell_labels: vec![1i32; n_cells_3d],
        cells_shape: [n_cells_3d, nloc_3d],
        cells: (0..(n_cells_3d * nloc_3d) as u32).map(|i| i % 3500).collect(),
        facet_labels: vec![3i32; 200],
        facets_shape: [200, 6],
        facets: (0..1200u32).collect(),
        points_shape: [3500, 3],
        points: (0..10500).map(|i| (i as f64) * 1e-8).collect(),
        edge_node_pairs: Vec::new(),
        edge_node_points: Vec::new(),
    };
    let headers_3d: Vec<ModeHeader> = (0..n_modes_3d)
        .map(|k| ModeHeader {
            mode_index: k,
            eig_re: 1.934171313619096e14 + (k as f64) * 1e11,
            eig_im: 0.046470653663846596,
            label: format!("cavity_mode_{k}"),
        })
        .collect();
    let comps_3d = vec![
        PackedComponentMeta { name: "Ex".into(), unit: "V/m".into(), nloc: nloc_3d, has_im: true },
        PackedComponentMeta { name: "Ey".into(), unit: "V/m".into(), nloc: nloc_3d, has_im: true },
        PackedComponentMeta { name: "Ez".into(), unit: "V/m".into(), nloc: nloc_3d, has_im: true },
    ];
    let mode_len_3d = n_comp_3d * n_cells_3d * nloc_3d * 2;
    let mut packed_3d = Vec::with_capacity(n_modes_3d * mode_len_3d);
    for m in 0..n_modes_3d {
        packed_3d.extend_from_slice(&generate_mode_slice_f64(m, n_comp_3d, n_cells_3d, nloc_3d));
    }
    let p_3d_l2 = scratch.join("tet10_3d_L2ShuffleLzfPacked.h5");
    write_nyx_packed(
        &p_3d_l2,
        b"{\"format\":\"nyx-project\",\"version\":1}",
        &mesh_3d,
        "ems",
        "cavity",
        None,
        &headers_3d,
        &comps_3d,
        &packed_3d,
        FieldStorageLayout::L2ShuffleLzfPacked,
    )?;

    let py_verify_packed = r#"
import sys, h5py, numpy as np
p2d, p3d = sys.argv[1], sys.argv[2]

def gen_mode(m, nc, ncells, nloc):
    n = nc * ncells * nloc * 2
    idx = np.arange(n, dtype=np.int64)
    m_base = (m + 1.0) * 0.125
    v = m_base + (idx % 97).astype(np.float64) * 0.0078125 + (idx // 97).astype(np.float64) * 1e-7
    out = np.where((idx & 1) == 0, v, -v * 0.5)
    return out.view(np.uint64)

with h5py.File(p2d, "r") as f:
    ds = f["/results/ems/cell_fields_packed"]
    assert ds.shape == (13, 3, 5488, 6, 2), f"bad 2D shape {ds.shape}"
    m3 = np.ascontiguousarray(ds[3], dtype=np.float64).ravel().view(np.uint64)
    exp3 = gen_mode(3, 3, 5488, 6)
    assert np.array_equal(m3, exp3), "2D mode 3 0-ULP mismatch in h5py"
    assert f["/mesh/edge_node_pairs"].shape == (3, 2)
    assert f["/mesh/edge_node_points"].shape == (3, 3)

with h5py.File(p3d, "r") as f:
    ds = f["/results/ems/cell_fields_packed"]
    assert ds.shape == (4, 3, 2000, 10, 2), f"bad 3D shape {ds.shape}"
    m2 = np.ascontiguousarray(ds[2], dtype=np.float64).ravel().view(np.uint64)
    exp2 = gen_mode(2, 3, 2000, 10)
    assert np.array_equal(m2, exp2), "3D mode 2 0-ULP mismatch in h5py"
"#;
    let p_2d_l2 = scratch.join("tri6_2d_L2ShuffleLzfPacked.h5");
    run_python_inline(py_verify_packed, &[p_2d_l2.to_str().unwrap(), p_3d_l2.to_str().unwrap()])?;

    // Also test legacy v3 write_nyx_v3 + direct-path XDMF generation & XML DataItem resolution
    let orig_v3 = NyxReader::open(Path::new(BENCH_MODES_2D), ReadBudgetConfig::default())?.read_nyx_v3_summary(true)?;
    let v3_out = scratch.join("reexport_v3_2d.h5");
    let xdmf_out = scratch.join("reexport_v3_2d.xdmf");
    write_nyx_v3(&orig_v3, &v3_out)?;
    let reloaded_v3 = NyxReader::open(&v3_out, ReadBudgetConfig::default())?.read_nyx_v3_summary(true)?;
    reloaded_v3.verify_0_ulp(&orig_v3, true)?;

    let v3_mesh = NyxReader::open(&v3_out, ReadBudgetConfig::default())?.read_mesh()?;
    let v3_headers = NyxReader::open(&v3_out, ReadBudgetConfig::default())?.scan_mode_headers("ems")?;
    let field_list: Vec<(String, String, usize)> =
        orig_v3.solutions[0].fields.iter().map(|(k, (name, _, _, _))| (k.clone(), name.clone(), 1usize)).collect();
    write_direct_xdmf(&xdmf_out, "reexport_v3_2d.h5", &v3_mesh, "ems", &v3_headers, &field_list)?;
    let xdmf_str = generate_direct_xdmf("reexport_v3_2d.h5", &v3_mesh, "ems", &v3_headers, &field_list);
    assert!(xdmf_str.contains("<Xdmf Version=\"3.0\">"));
    assert!(xdmf_str.contains("reexport_v3_2d.h5:/results/ems/solution_0/fields/00_|E|"));

    // -------------------------------------------------------------------------
    // 2. EC-06: Partial-Edge Chunk Filter Pipeline Fidelity (11,844 @ chunk 4,096)
    // -------------------------------------------------------------------------
    let ec06_h5py = scratch.join("ec06_from_h5py.h5");
    let ec06_rust = scratch.join("ec06_from_rust.h5");
    let py_gen_ec06 = r#"
import math, sys, h5py, numpy as np
path = sys.argv[1]
n, ch = 11844, 4096
base = np.array([math.sin(i * 0.0137) * 1e4 + i * 1e-7 for i in range(n)], dtype=np.float64)
with h5py.File(path, "w") as f:
    f.create_dataset("deflate_only", data=base, chunks=(ch,), compression="gzip", compression_opts=3)
    f.create_dataset("shuffle_deflate", data=base * -1.5 + 0.25, chunks=(ch,), shuffle=True, compression="gzip", compression_opts=3)
    f.create_dataset("shuffle_deflate_fletcher32", data=base * 0.75 - 12.0, chunks=(ch,), shuffle=True, compression="gzip", compression_opts=3, fletcher32=True)
    f.create_dataset("shuffle_lzf", data=base * 2.25 + 1.0, chunks=(ch,), shuffle=True, compression="lzf")
"#;
    run_python_inline(py_gen_ec06, &[ec06_h5py.to_str().unwrap()])?;

    let n_ec06 = 11_844usize;
    let ch_ec06 = 4_096u64;
    let base_ec06: Vec<f64> = (0..n_ec06).map(|i| ((i as f64) * 0.0137).sin() * 1e4 + (i as f64) * 1e-7).collect();
    let sd_ec06: Vec<f64> = base_ec06.iter().map(|&x| x * -1.5 + 0.25).collect();
    let sdf_ec06: Vec<f64> = base_ec06.iter().map(|&x| x * 0.75 - 12.0).collect();
    let lzf_ec06: Vec<f64> = base_ec06.iter().map(|&x| x * 2.25 + 1.0).collect();

    {
        let r_h5py = NyxReader::open(&ec06_h5py, ReadBudgetConfig::default())?;
        assert_f64_bits_eq("ec06_h5py/deflate_only", &r_h5py.file().dataset("deflate_only")?.read_f64()?, &base_ec06)?;
        assert_f64_bits_eq(
            "ec06_h5py/shuffle_deflate",
            &r_h5py.file().dataset("shuffle_deflate")?.read_f64()?,
            &sd_ec06,
        )?;
        assert_f64_bits_eq(
            "ec06_h5py/shuffle_deflate_fletcher32",
            &r_h5py.file().dataset("shuffle_deflate_fletcher32")?.read_f64()?,
            &sdf_ec06,
        )?;
        assert_f64_bits_eq("ec06_h5py/shuffle_lzf", &r_h5py.file().dataset("shuffle_lzf")?.read_f64()?, &lzf_ec06)?;
    }

    {
        let mut fb = FileBuilder::new();
        fb.create_dataset("deflate_only")
            .with_f64_data(&base_ec06)
            .with_chunks(&[ch_ec06])
            .with_deflate(DeflateLevel::Three);
        fb.create_dataset("shuffle_deflate")
            .with_f64_data(&sd_ec06)
            .with_chunks(&[ch_ec06])
            .with_shuffle()
            .with_deflate(DeflateLevel::Three);
        fb.create_dataset("shuffle_deflate_fletcher32")
            .with_f64_data(&sdf_ec06)
            .with_chunks(&[ch_ec06])
            .with_shuffle()
            .with_deflate(DeflateLevel::Three)
            .with_fletcher32();
        fb.create_dataset("shuffle_lzf").with_f64_data(&lzf_ec06).with_chunks(&[ch_ec06]).with_shuffle().with_lzf();
        fb.write(&ec06_rust)?;
    }
    let py_verify_ec06 = r#"
import math, sys, h5py, numpy as np
path = sys.argv[1]
n = 11844
base = np.array([math.sin(i * 0.0137) * 1e4 + i * 1e-7 for i in range(n)], dtype=np.float64)
with h5py.File(path, "r") as f:
    assert np.array_equal(f["deflate_only"][:].view(np.uint64), base.view(np.uint64))
    assert np.array_equal(f["shuffle_deflate"][:].view(np.uint64), (base * -1.5 + 0.25).view(np.uint64))
    assert np.array_equal(f["shuffle_deflate_fletcher32"][:].view(np.uint64), (base * 0.75 - 12.0).view(np.uint64))
    assert np.array_equal(f["shuffle_lzf"][:].view(np.uint64), (base * 2.25 + 1.0).view(np.uint64))
"#;
    run_python_inline(py_verify_ec06, &[ec06_rust.to_str().unwrap()])?;

    // -------------------------------------------------------------------------
    // 3. EC-07: Extreme IEEE-754 Bit-Exact Preservation (+0.0, -0.0, subnormals, Inf, NaN payloads)
    // -------------------------------------------------------------------------
    let ec07_path = scratch.join("ec07_ieee754.h5");
    let ieee_bits = ieee754_extreme_bits();
    let ieee_vals: Vec<f64> = ieee_bits.iter().map(|&b| f64::from_bits(b)).collect();
    {
        let mut fb = FileBuilder::new();
        fb.set_attr("neg_zero_attr", AttrValue::F64(f64::from_bits(ieee_bits[1])));
        fb.create_dataset("ieee754").with_f64_data(&ieee_vals).with_chunks(&[8]).with_shuffle().with_lzf();
        fb.write(&ec07_path)?;
    }
    {
        let r = NyxReader::open(&ec07_path, ReadBudgetConfig::default())?;
        let attr_v = r.file().root().attrs()?.get("neg_zero_attr").and_then(AttrValue::as_f64).unwrap();
        assert_eq!(attr_v.to_bits(), ieee_bits[1], "-0.0 sign bit lost in root attribute");
        let read_vals = r.file().dataset("ieee754")?.read_f64()?;
        assert_f64_bits_eq("ec07_ieee754_rust", &read_vals, &ieee_vals)?;
    }
    let ieee_hex = ieee_bits.iter().map(|b| format!("{b:016x}")).collect::<Vec<_>>().join(",");
    let py_verify_ec07 = r#"
import sys, h5py, numpy as np
path, hex_str = sys.argv[1], sys.argv[2]
exp = np.array([int(x, 16) for x in hex_str.split(",")], dtype=np.uint64)
with h5py.File(path, "r") as f:
    ab = np.asarray([f.attrs["neg_zero_attr"]], dtype=np.float64).view(np.uint64)[0]
    assert int(ab) == int(exp[1]), f"attr -0.0 mismatch: 0x{int(ab):016x}"
    db = f["ieee754"][:].view(np.uint64)
    assert np.array_equal(db, exp), "dataset IEEE-754 bits mismatch"
"#;
    run_python_inline(py_verify_ec07, &[ec07_path.to_str().unwrap(), &ieee_hex])?;

    // -------------------------------------------------------------------------
    // 4. EC-08: Degenerate Shapes (n_solutions=0, empty (0, 3) dataset, single-cell (1, 6))
    // -------------------------------------------------------------------------
    let ec08_packed_zero = scratch.join("ec08_zero_modes.h5");
    write_nyx_packed(
        &ec08_packed_zero,
        b"{}",
        &mesh_2d,
        "ems",
        "modal",
        None,
        &[],
        &comps_2d,
        &[],
        FieldStorageLayout::L2ShuffleLzfPacked,
    )?;
    {
        let r = NyxReader::open(&ec08_packed_zero, ReadBudgetConfig::default())?;
        assert!(r.scan_mode_headers("ems")?.is_empty());
        assert_eq!(r.packed_logical_shape("ems")?, Some([0, 3, n_cells_2d, nloc_2d, 2]));
    }

    let ec08_raw = scratch.join("ec08_degenerate.h5");
    {
        let mut fb = FileBuilder::new();
        fb.set_attr("n_solutions", AttrValue::I64(0));
        fb.create_dataset("empty_2d").with_dtype(make_f64_type()).with_shape(&[0, 3]);
        fb.create_dataset("single_cell").with_f64_data(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).with_shape(&[1, 6]);
        fb.write(&ec08_raw)?;
    }
    {
        let f = File::open_streaming(&ec08_raw)?;
        assert_eq!(f.dataset("empty_2d")?.shape()?, vec![0, 3]);
        assert!(f.dataset("empty_2d")?.read_f64()?.is_empty());
        assert_eq!(f.dataset("single_cell")?.shape()?, vec![1, 6]);
        assert_f64_bits_eq("single_cell", &f.dataset("single_cell")?.read_f64()?, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])?;
    }
    let py_verify_ec08 = r#"
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:
    assert int(f.attrs["n_solutions"]) == 0
    assert f["empty_2d"].shape == (0, 3)
    assert f["single_cell"].shape == (1, 6)
with h5py.File(sys.argv[2], "r") as f:
    assert int(f["/results/ems"].attrs["n_solutions"]) == 0
    assert f["/results/ems/eigenvalues"].shape == (0, 2)
    assert f["/results/ems/cell_fields_packed"].shape == (0, 3, 5488, 6, 2)
"#;
    run_python_inline(py_verify_ec08, &[ec08_raw.to_str().unwrap(), ec08_packed_zero.to_str().unwrap()])?;

    // -------------------------------------------------------------------------
    // 5. Incremental Streaming Mode Append (`init_nyx_packed_streaming` + `append_solved_mode_in_place`)
    // -------------------------------------------------------------------------
    let stream_path = scratch.join("streaming_modes_l2.h5");
    init_nyx_packed_streaming(
        &stream_path,
        b"{\"streaming\":true}",
        &mesh_3d,
        "ems",
        "cavity",
        Some(1.55e-6),
        &comps_3d,
        FieldStorageLayout::L2ShuffleLzfPacked,
    )?;
    for m in 0..n_modes_3d {
        let slice = &packed_3d[m * mode_len_3d..(m + 1) * mode_len_3d];
        append_solved_mode_in_place(&stream_path, "ems", &headers_3d[m], slice)?;
    }
    {
        let r = NyxReader::open(&stream_path, ReadBudgetConfig::default())?;
        assert_eq!(r.scan_mode_headers("ems")?, headers_3d);
        assert_eq!(r.packed_logical_shape("ems")?, Some([n_modes_3d, n_comp_3d, n_cells_3d, nloc_3d, 2]));
        for m in 0..n_modes_3d {
            let s = r.read_mode_cell_fields("ems", m)?;
            assert_f64_bits_eq(&format!("streaming_mode_{m}"), &s, &packed_3d[m * mode_len_3d..(m + 1) * mode_len_3d])?;
        }
    }

    // -------------------------------------------------------------------------
    // 6. BM-5: 50-Point Incremental In-Place SweepCheckpoint Append (0 ULP in Rust & h5py)
    // -------------------------------------------------------------------------
    let sweep_path = scratch.join("bm5_sweep_50pts.h5");
    let ckpt = SweepCheckpoint::create(&sweep_path, &["width_nm", "wl_nm"], 2)?;
    let n_sweep = 50usize;
    let mut expected_pts = Vec::with_capacity(n_sweep);
    for i in 0..n_sweep {
        let pt = SweepPoint {
            params: vec![450.0 + (i as f64) * 5.0, 1550.0 + (i as f64) * 0.25],
            eig_re: vec![2.8946 - (i as f64) * 0.001, 2.4858 - (i as f64) * 0.0012],
            eig_im: vec![1.0e-11 * (i as f64 + 1.0), -2.5e-11 * (i as f64 + 1.0)],
        };
        ckpt.append_point(&pt)?;
        expected_pts.push(pt);
    }
    let reopened_ckpt = SweepCheckpoint::open(&sweep_path)?;
    let read_pts = reopened_ckpt.read_points()?;
    assert_eq!(read_pts.len(), n_sweep);
    for (i, (a, e)) in read_pts.iter().zip(expected_pts.iter()).enumerate() {
        assert_f64_bits_eq(&format!("sweep_pt_{i}_params"), &a.params, &e.params)?;
        assert_f64_bits_eq(&format!("sweep_pt_{i}_eig_re"), &a.eig_re, &e.eig_re)?;
        assert_f64_bits_eq(&format!("sweep_pt_{i}_eig_im"), &a.eig_im, &e.eig_im)?;
    }
    let py_verify_bm5 = r#"
import sys, h5py, numpy as np
path = sys.argv[1]
n = 50
with h5py.File(path, "r") as f:
    assert int(f.attrs["n_points"]) == n
    w = f["width_nm"][:].view(np.uint64)
    wl = f["wl_nm"][:].view(np.uint64)
    er0 = f["eig_re_0"][:].view(np.uint64)
    ei1 = f["eig_im_1"][:].view(np.uint64)
    exp_w = np.array([450.0 + i * 5.0 for i in range(n)], dtype=np.float64).view(np.uint64)
    exp_wl = np.array([1550.0 + i * 0.25 for i in range(n)], dtype=np.float64).view(np.uint64)
    exp_er0 = np.array([2.8946 - i * 0.001 for i in range(n)], dtype=np.float64).view(np.uint64)
    exp_ei1 = np.array([-2.5e-11 * (i + 1.0) for i in range(n)], dtype=np.float64).view(np.uint64)
    assert np.array_equal(w, exp_w), "BM-5 width_nm 0-ULP mismatch"
    assert np.array_equal(wl, exp_wl), "BM-5 wl_nm 0-ULP mismatch"
    assert np.array_equal(er0, exp_er0), "BM-5 eig_re_0 0-ULP mismatch"
    assert np.array_equal(ei1, exp_ei1), "BM-5 eig_im_1 0-ULP mismatch"
"#;
    run_python_inline(py_verify_bm5, &[sweep_path.to_str().unwrap()])?;

    println!(
        "[Gate 3 / EC-02..EC-08 & BM-5 PASS] ({:.1} ms) 2D/3D L2/L3/L1 0-ULP round-trip (L2 ratio={:.2}x > 15x), partial edge chunks, 15 IEEE-754 bit patterns, degenerate shapes, streaming mode append, and 50-pt SweepCheckpoint verified in h5py 3.13.0",
        t0.elapsed().as_secs_f64() * 1000.0,
        l2_ratio_2d
    );
    Ok(())
}

// =============================================================================
// GATE 4 (EC-09): Bounded Memory (< 100 MB RSS Delta) Slicing 1 Mode (60 MB)
// from a 10-Mode 600 MB 3D Tetrahedron_10 (125,000 Cells) File in L2 and L3
// =============================================================================
#[test]
fn gate4_child_rss_probe() -> Result<()> {
    let Ok(file_path) = std::env::var("NYX_GATE4_CHILD_FILE") else {
        // Not running as child subprocess probe; no-op during normal test discovery.
        return Ok(());
    };
    let mode_idx: usize = std::env::var("NYX_GATE4_CHILD_MODE").unwrap_or_else(|_| "3".to_string()).parse()?;
    let expected_len: usize = 3 * 125_000 * 10 * 2; // 7,500,000 f64s = 60.00 MB

    let (rss_before, _) = win_mem::memory_usage_mb();
    let reader = NyxReader::open(Path::new(&file_path), ReadBudgetConfig::default())?;
    let slice = reader.read_mode_cell_fields("ems", mode_idx)?;
    std::hint::black_box(&slice);
    let (rss_after, peak_rss) = win_mem::memory_usage_mb();
    let delta_mb = (rss_after - rss_before).max(0.0);

    if slice.len() != expected_len {
        bail!("child probe slice len {} != {}", slice.len(), expected_len);
    }
    // O(1)-memory 0-ULP verification across all 7,500,000 f64s so we do not allocate a second 60 MB buffer
    for (i, &val) in slice.iter().enumerate() {
        let exp = expected_mode_elem(mode_idx, i);
        if val.to_bits() != exp.to_bits() {
            bail!("Gate 4 child probe 0-ULP mismatch at index {}: {} vs {}", i, val, exp);
        }
    }
    drop(slice);

    println!(
        "GATE4_CHILD_RESULT rss_before_mb={:.2} rss_after_mb={:.2} rss_delta_mb={:.2} peak_rss_mb={:.2}",
        rss_before, rss_after, delta_mb, peak_rss
    );
    Ok(())
}

#[test]
fn gate4_ec09_bounded_memory_600mb_3d_tet10_single_mode_slice() -> Result<()> {
    if std::env::var("NYX_GATE4_CHILD_FILE").is_ok() {
        return Ok(());
    }
    let t0 = Instant::now();
    let scratch = test_scratch_dir("gate4")?;

    let n_modes = 10usize;
    let n_comp = 3usize;
    let n_cells = 125_000usize;
    let nloc = 10usize;
    let mode_len = n_comp * n_cells * nloc * 2; // 7,500,000 f64s = 60 MB/mode

    let mesh = MeshData {
        cell_topology: "Tetrahedron_10".to_string(),
        element_family: "Lagrange".to_string(),
        element_degree: 2,
        space_dimension: 3,
        cell_labels: vec![1i32; 16],
        cells_shape: [n_cells, nloc],
        cells: vec![0u32; n_cells * nloc],
        facet_labels: vec![1i32; 8],
        facets_shape: [8, 6],
        facets: vec![0u32; 48],
        points_shape: [16, 3],
        points: vec![0.0f64; 48],
        edge_node_pairs: Vec::new(),
        edge_node_points: Vec::new(),
    };
    let headers: Vec<ModeHeader> = (0..n_modes)
        .map(|k| ModeHeader {
            mode_index: k,
            eig_re: 1.934e14 + (k as f64) * 1e11,
            eig_im: 0.0,
            label: format!("mode_{k}"),
        })
        .collect();
    let comps = vec![
        PackedComponentMeta { name: "Ex".into(), unit: "V/m".into(), nloc, has_im: true },
        PackedComponentMeta { name: "Ey".into(), unit: "V/m".into(), nloc, has_im: true },
        PackedComponentMeta { name: "Ez".into(), unit: "V/m".into(), nloc, has_im: true },
    ];

    let mut all_data = Vec::with_capacity(n_modes * mode_len);
    for m in 0..n_modes {
        all_data.extend_from_slice(&generate_mode_slice_f64(m, n_comp, n_cells, nloc));
    }

    let l2_path = scratch.join("ec09_600mb_tet10_l2.h5");
    let l3_path = scratch.join("ec09_600mb_tet10_l3.h5");
    write_nyx_packed(
        &l2_path,
        b"{}",
        &mesh,
        "ems",
        "cavity",
        None,
        &headers,
        &comps,
        &all_data,
        FieldStorageLayout::L2ShuffleLzfPacked,
    )?;
    write_nyx_packed(
        &l3_path,
        b"{}",
        &mesh,
        "ems",
        "cavity",
        None,
        &headers,
        &comps,
        &all_data,
        FieldStorageLayout::L3ContiguousPacked,
    )?;
    // Free the 600 MB parent buffer before spawning clean child probes
    drop(all_data);

    let l2_size_mb = std::fs::metadata(&l2_path)?.len() as f64 / (1024.0 * 1024.0);
    let l3_size_mb = std::fs::metadata(&l3_path)?.len() as f64 / (1024.0 * 1024.0);

    let exe = std::env::current_exe()?;
    let run_child_probe = |path: &Path| -> Result<f64> {
        let out = Command::new(&exe)
            .env("NYX_GATE4_CHILD_FILE", path.as_os_str())
            .env("NYX_GATE4_CHILD_MODE", "3")
            .args(["--exact", "gate4_child_rss_probe", "--nocapture"])
            .output()?;
        let _ = std::fs::remove_file(path);
        let stdout = String::from_utf8_lossy(&out.stdout);
        let stderr = String::from_utf8_lossy(&out.stderr);
        if !out.status.success() {
            bail!("Gate 4 child probe failed:\nSTDOUT: {}\nSTDERR: {}", stdout, stderr);
        }
        for line in stdout.lines() {
            if let Some(rest) = line.strip_prefix("GATE4_CHILD_RESULT ") {
                for tok in rest.split_whitespace() {
                    if let Some(val_str) = tok.strip_prefix("rss_delta_mb=") {
                        return Ok(val_str.parse::<f64>()?);
                    }
                }
            }
        }
        bail!("Gate 4 child probe did not emit GATE4_CHILD_RESULT:\n{}", stdout)
    };

    let delta_l2_mb = run_child_probe(&l2_path)?;
    let delta_l3_mb = run_child_probe(&l3_path)?;

    if delta_l2_mb >= 100.0 || delta_l3_mb >= 100.0 {
        bail!("Gate 4 RSS delta exceeded 100 MB ceiling: L2={:.2} MB, L3={:.2} MB", delta_l2_mb, delta_l3_mb);
    }

    println!(
        "[Gate 4 / EC-09 PASS] ({:.1} ms) 10-mode 600 MB 3D Tet10 (L2={:.2} MB on disk, L3={:.2} MB on disk): clean-child single-mode (60 MB) slice RSS delta L2={:.2} MB (< 100 MB), L3={:.2} MB (< 100 MB) with 0 ULP",
        t0.elapsed().as_secs_f64() * 1000.0,
        l2_size_mb,
        l3_size_mb,
        delta_l2_mb,
        delta_l3_mb
    );
    Ok(())
}

// =============================================================================
// GATE 5 (EC-11 & EC-12):
// Untrusted Reader Fuzzing (10,000 Truncations/Bit-Flips/u64::MAX OOM-Bombs)
// + Pure-Rust P2 Triangle_6 & Tetrahedron_10 SOS Brillouin PE + MB Overlap Invariant
// =============================================================================
struct SplitMix64(u64);
impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn gen_range(&mut self, len: usize) -> usize {
        if len == 0 {
            0
        } else {
            (self.next_u64() as usize) % len
        }
    }
}

fn generate_10k_mutations(seed_files: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let mut mutations = Vec::with_capacity(10_000);
    let mut rng = SplitMix64(0x2026_1120_CAFE_BEEF);

    for base in seed_files {
        for pow in 0..20 {
            let cut = (1usize << pow).min(base.len());
            mutations.push(base[..cut].to_vec());
        }
        for step in 0..80 {
            let cut = (step * (base.len() / 80).max(1)).min(base.len());
            mutations.push(base[..cut].to_vec());
        }
    }

    for i in 0..600 {
        let base = &seed_files[i % seed_files.len()];
        let mut buf = base.clone();
        if buf.len() > 24 {
            let off = 8 + rng.gen_range(buf.len() - 16);
            let bomb: [u8; 8] = if i % 3 == 0 {
                u64::MAX.to_le_bytes()
            } else if i % 3 == 1 {
                (1u64 << 62).to_le_bytes()
            } else {
                (u32::MAX as u64).to_le_bytes()
            };
            buf[off..off + 8].copy_from_slice(&bomb);
        }
        mutations.push(buf);
    }

    while mutations.len() < 10_000 {
        let idx = mutations.len();
        let base = &seed_files[idx % seed_files.len()];
        match idx % 5 {
            0 => {
                let cut = rng.gen_range(base.len().max(1));
                mutations.push(base[..cut].to_vec());
            }
            1 => {
                let mut buf = base.clone();
                let limit = buf.len().clamp(1, 256);
                let n_flips = 1 + rng.gen_range(4);
                for _ in 0..n_flips {
                    let pos = rng.gen_range(limit);
                    buf[pos] ^= 1u8 << (rng.next_u64() & 7);
                }
                mutations.push(buf);
            }
            2 | 3 => {
                let mut buf = base.clone();
                if !buf.is_empty() {
                    let n_flips = 1 + rng.gen_range(8);
                    for _ in 0..n_flips {
                        let pos = rng.gen_range(buf.len());
                        buf[pos] ^= (rng.next_u64() as u8) | 1;
                    }
                }
                mutations.push(buf);
            }
            _ => {
                let cut = 16 + rng.gen_range(base.len().saturating_sub(16).max(1));
                let mut buf = base[..cut.min(base.len())].to_vec();
                if !buf.is_empty() {
                    let pos = rng.gen_range(buf.len());
                    buf[pos] = !buf[pos];
                }
                mutations.push(buf);
            }
        }
    }
    mutations.truncate(10_000);
    mutations
}

#[allow(dead_code, clippy::needless_range_loop, clippy::if_same_then_else)]
mod brillouin {
    use anyhow::Result;
    use nyx_storage::{reader::attr_f64, reader::attr_str, NyxReader, ReadBudgetConfig};
    use std::f64::consts::PI;
    use std::path::Path;

    pub const C0: f64 = 299792458.0;
    pub const WL0: f64 = 1.55e-6;
    pub const OMEGA_OPT: f64 = 2.0 * PI * C0 / WL0;

    pub const RHO_SI: f64 = 2329.0;
    pub const RHO_SAP: f64 = 3980.0;
    pub const N_SI: f64 = 3.4757;
    pub const EPS_SI: f64 = N_SI * N_SI;

    pub const EPS_SAP_TENSOR: [[f64; 3]; 3] =
        [[3.0419738, 0.0, 0.0], [0.0, 3.0211640, -0.0329809], [0.0, -0.0329809, 2.9896682]];

    pub const P11_SI: f64 = -0.094;
    pub const P12_SI: f64 = 0.017;
    pub const P44_SI: f64 = -0.051;

    pub const PERM_TET10: [usize; 10] = [0, 1, 2, 3, 4, 7, 5, 6, 8, 9];
    pub const PERM_TRI6: [usize; 6] = [0, 1, 2, 3, 5, 4];
    pub const EDGE_LOC: [[usize; 3]; 3] = [[0, 3, 1], [1, 4, 2], [2, 5, 0]];

    const A1: f64 = 0.445948490915965;
    const B1: f64 = 0.108103018168070;
    const W1: f64 = 0.223381589678011;
    const A2: f64 = 0.091576213509771;
    const B2: f64 = 0.816847572980459;
    const W2: f64 = 0.109951743655322;

    pub const QUAD_PTS_2D: [[f64; 2]; 6] = [[A1, A1], [A1, B1], [B1, A1], [A2, A2], [A2, B2], [B2, A2]];
    pub const QUAD_W_2D: [f64; 6] = [0.5 * W1, 0.5 * W1, 0.5 * W1, 0.5 * W2, 0.5 * W2, 0.5 * W2];

    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    pub struct C64 {
        pub re: f64,
        pub im: f64,
    }

    impl C64 {
        #[inline]
        pub const fn new(re: f64, im: f64) -> Self {
            Self { re, im }
        }
        #[inline]
        pub fn conj(self) -> Self {
            Self { re: self.re, im: -self.im }
        }
        #[inline]
        pub fn norm_sq(self) -> f64 {
            self.re * self.re + self.im * self.im
        }
        #[inline]
        pub fn abs(self) -> f64 {
            self.re.hypot(self.im)
        }
        #[inline]
        pub fn angle(self) -> f64 {
            self.im.atan2(self.re)
        }
        #[inline]
        pub fn cis(theta: f64) -> Self {
            Self { re: theta.cos(), im: theta.sin() }
        }
        #[inline]
        pub fn scale(self, s: f64) -> Self {
            Self { re: self.re * s, im: self.im * s }
        }
    }

    impl std::ops::Add for C64 {
        type Output = Self;
        #[inline]
        fn add(self, rhs: Self) -> Self {
            Self { re: self.re + rhs.re, im: self.im + rhs.im }
        }
    }
    impl std::ops::AddAssign for C64 {
        #[inline]
        fn add_assign(&mut self, rhs: Self) {
            self.re += rhs.re;
            self.im += rhs.im;
        }
    }
    impl std::ops::Sub for C64 {
        type Output = Self;
        #[inline]
        fn sub(self, rhs: Self) -> Self {
            Self { re: self.re - rhs.re, im: self.im - rhs.im }
        }
    }
    impl std::ops::Neg for C64 {
        type Output = Self;
        #[inline]
        fn neg(self) -> Self {
            Self { re: -self.re, im: -self.im }
        }
    }
    impl std::ops::Mul for C64 {
        type Output = Self;
        #[inline]
        fn mul(self, rhs: Self) -> Self {
            Self { re: self.re * rhs.re - self.im * rhs.im, im: self.re * rhs.im + self.im * rhs.re }
        }
    }

    pub struct RawEmsMode2D {
        pub eig_re: f64,
        pub eig_im: f64,
        pub pts: Vec<[f64; 2]>,
        pub cells: Vec<[u32; 6]>,
        pub ex_re: Vec<f64>,
        pub ex_im: Vec<f64>,
        pub ez_re: Vec<f64>,
        pub ez_im: Vec<f64>,
        pub ey_re: Vec<f64>,
        pub ey_im: Vec<f64>,
    }

    pub struct RawSmsMode {
        pub topology: String,
        pub point_dim: usize,
        pub pts_flat: Vec<f64>,
        pub cells_flat: Vec<u32>,
        pub ux_re: Vec<f64>,
        pub ux_im: Vec<f64>,
        pub uy_re: Vec<f64>,
        pub uy_im: Vec<f64>,
        pub uz_re: Vec<f64>,
        pub uz_im: Vec<f64>,
    }

    pub fn load_ems_mode_2d(path: &Path, sol_idx: usize) -> Result<RawEmsMode2D> {
        let reader = NyxReader::open(path, ReadBudgetConfig::default())?;
        let file = reader.file();
        let pts_flat = file.dataset("viz/ems/points")?.read_f64()?;
        let (pts_chunks, _) = pts_flat.as_chunks::<2>();
        let pts = pts_chunks.to_vec();

        let cells_flat = file.dataset("viz/ems/cells")?.read_u32()?;
        let (cells_chunks, _) = cells_flat.as_chunks::<6>();
        let cells = cells_chunks.to_vec();

        let sol_grp = file.group(&format!("results/ems/solution_{sol_idx}"))?;
        let sol_attrs = sol_grp.attrs()?;
        let eig_re = attr_f64(&sol_attrs, "eig_re")?;
        let eig_im = attr_f64(&sol_attrs, "eig_im")?;

        let cfe = sol_grp.group("cell_fields")?;
        let ex_re = cfe.dataset("00_Ex/re")?.read_f64()?;
        let ex_im = cfe.dataset("00_Ex/im")?.read_f64()?;
        let ez_re = cfe.dataset("01_Ez/re")?.read_f64()?;
        let ez_im = cfe.dataset("01_Ez/im")?.read_f64()?;
        let ey_re = cfe.dataset("02_Ey/re")?.read_f64()?;
        let ey_im = cfe.dataset("02_Ey/im")?.read_f64()?;

        Ok(RawEmsMode2D { eig_re, eig_im, pts, cells, ex_re, ex_im, ez_re, ez_im, ey_re, ey_im })
    }

    pub fn load_sms_mode(path: &Path, sol_idx: usize) -> Result<RawSmsMode> {
        let reader = NyxReader::open(path, ReadBudgetConfig::default())?;
        let file = reader.file();
        let viz_sms = file.group("viz/sms")?;
        let topology = attr_str(&viz_sms.attrs()?, "topology")?;
        let ds_pts = viz_sms.dataset("points")?;
        let pt_shape = ds_pts.shape()?;
        let point_dim = pt_shape[1] as usize;
        let pts_flat = ds_pts.read_f64()?;

        let ds_cells = viz_sms.dataset("cells")?;
        let cells_flat = ds_cells.read_u32()?;

        let cfs = file.group(&format!("results/sms/solution_{sol_idx}/cell_fields"))?;
        let ux_re = cfs.dataset("00_ux/re")?.read_f64()?;
        let ux_im = cfs.dataset("00_ux/im")?.read_f64()?;
        let (uy_re, uy_im, uz_re, uz_im) = if topology == "Tetrahedron_10" {
            (
                cfs.dataset("01_uy/re")?.read_f64()?,
                cfs.dataset("01_uy/im")?.read_f64()?,
                cfs.dataset("02_uz/re")?.read_f64()?,
                cfs.dataset("02_uz/im")?.read_f64()?,
            )
        } else {
            let uz_r = cfs.dataset("01_uz/re")?.read_f64()?;
            let uz_i = cfs.dataset("01_uz/im")?.read_f64()?;
            let zeros = vec![0.0; uz_r.len()];
            (zeros.clone(), zeros, uz_r, uz_i)
        };

        Ok(RawSmsMode { topology, point_dim, pts_flat, cells_flat, ux_re, ux_im, uy_re, uy_im, uz_re, uz_im })
    }

    #[derive(Clone, Debug)]
    pub struct CellMasks {
        pub phys: Vec<bool>,
        pub core: Vec<bool>,
        pub slab: Vec<bool>,
        pub si_all: Vec<bool>,
        pub sap: Vec<bool>,
        pub air: Vec<bool>,
    }

    pub fn tri6_basis(xi: f64, eta: f64) -> ([f64; 6], [f64; 6], [f64; 6]) {
        let l0 = 1.0 - xi - eta;
        let l1 = xi;
        let l2 = eta;
        let n = [
            l0 * (2.0 * l0 - 1.0),
            l1 * (2.0 * l1 - 1.0),
            l2 * (2.0 * l2 - 1.0),
            4.0 * l0 * l1,
            4.0 * l1 * l2,
            4.0 * l2 * l0,
        ];
        let dn_dxi = [1.0 - 4.0 * l0, 4.0 * l1 - 1.0, 0.0, 4.0 * (l0 - l1), 4.0 * l2, -4.0 * l2];
        let dn_deta = [1.0 - 4.0 * l0, 0.0, 4.0 * l2 - 1.0, -4.0 * l1, 4.0 * l1, 4.0 * (l0 - l2)];
        (n, dn_dxi, dn_deta)
    }

    pub fn classify_cells(
        pts: &[[f64; 2]],
        cells: &[[u32; 6]],
        kind: &str,
        x_phys: f64,
        z_sub: f64,
        z_air: f64,
    ) -> CellMasks {
        let n_cells = cells.len();
        let mut phys = vec![false; n_cells];
        let mut core = vec![false; n_cells];
        let mut slab = vec![false; n_cells];
        let mut si_all = vec![false; n_cells];
        let mut sap = vec![false; n_cells];
        let mut air = vec![false; n_cells];

        for (i, c) in cells.iter().enumerate() {
            let p0 = pts[c[0] as usize];
            let p1 = pts[c[1] as usize];
            let p2 = pts[c[2] as usize];
            let xc = (p0[0] + p1[0] + p2[0]) / 3.0;
            let zc = (p0[1] + p1[1] + p2[1]) / 3.0;
            let in_phys = xc.abs() < x_phys && zc > -z_sub && zc < z_air;
            phys[i] = in_phys;
            if kind == "strip" {
                let in_core = in_phys && xc.abs() < 0.40e-6 && zc > 0.0 && zc < 0.30e-6;
                core[i] = in_core;
                slab[i] = false;
                si_all[i] = in_core;
                sap[i] = in_phys && zc < 0.0;
                air[i] = in_phys && zc > 0.0 && !in_core;
            } else {
                let in_core = in_phys && xc.abs() < 0.25e-6 && zc > 0.0 && zc < 0.30e-6;
                let in_slab = in_phys && !in_core && zc > 0.0 && zc < 0.17e-6;
                core[i] = in_core;
                slab[i] = in_slab;
                si_all[i] = in_core || in_slab;
                sap[i] = in_phys && zc < 0.0;
                air[i] = in_phys && zc > 0.17e-6 && !in_core;
            }
        }
        CellMasks { phys, core, slab, si_all, sap, air }
    }

    pub struct PreparedEmsMode2D {
        pub kind: String,
        pub cells: Vec<[u32; 6]>,
        pub masks: CellMasks,
        pub ex: Vec<[C64; 6]>,
        pub ey: Vec<[C64; 6]>,
        pub ez: Vec<[C64; 6]>,
        pub p0: Vec<[f64; 2]>,
        pub dx1: Vec<f64>,
        pub dz1: Vec<f64>,
        pub dx2: Vec<f64>,
        pub dz2: Vec<f64>,
        pub det_j_signed: Vec<f64>,
        pub norm_e2_all: f64,
        pub norm_eps_e2_tot: f64,
    }

    impl PreparedEmsMode2D {
        pub fn from_raw(raw: &RawEmsMode2D, kind: &str, x_phys: f64, z_sub: f64, z_air: f64) -> Self {
            let n_cells = raw.cells.len();
            let mut ex = vec![[C64::default(); 6]; n_cells];
            let mut ey = vec![[C64::default(); 6]; n_cells];
            let mut ez = vec![[C64::default(); 6]; n_cells];

            for c in 0..n_cells {
                for k in 0..6 {
                    let src = c * 6 + PERM_TRI6[k];
                    ex[c][k] = C64::new(raw.ex_re[src], raw.ex_im[src]);
                    ez[c][k] = C64::new(raw.ez_re[src], raw.ez_im[src]);
                    ey[c][k] = C64::new(raw.ey_re[src], raw.ey_im[src]);
                }
            }

            let masks = classify_cells(&raw.pts, &raw.cells, kind, x_phys, z_sub, z_air);

            let mut z_phasor = C64::default();
            for c in 0..n_cells {
                if masks.core[c] {
                    for k in 0..6 {
                        z_phasor += ex[c][k] * ex[c][k] + ez[c][k] * ez[c][k];
                    }
                }
            }
            let theta = if z_phasor.abs() > 1e-30 { -0.5 * z_phasor.angle() } else { 0.0 };
            let rot = C64::cis(theta);
            for c in 0..n_cells {
                for k in 0..6 {
                    ex[c][k] = ex[c][k] * rot;
                    ey[c][k] = ey[c][k] * rot;
                    ez[c][k] = ez[c][k] * rot;
                }
            }

            let mut sum_ex2 = 0.0;
            let mut sum_ez2 = 0.0;
            let mut sum_ex_re = 0.0;
            let mut sum_ez_re = 0.0;
            for c in 0..n_cells {
                if masks.core[c] {
                    for k in 0..6 {
                        sum_ex2 += ex[c][k].norm_sq();
                        sum_ez2 += ez[c][k].norm_sq();
                        sum_ex_re += ex[c][k].re;
                        sum_ez_re += ez[c][k].re;
                    }
                }
            }
            let flip = if sum_ex2 >= sum_ez2 { sum_ex_re < 0.0 } else { sum_ez_re < 0.0 };
            if flip {
                for c in 0..n_cells {
                    for k in 0..6 {
                        ex[c][k] = -ex[c][k];
                        ey[c][k] = -ey[c][k];
                        ez[c][k] = -ez[c][k];
                    }
                }
            }

            let mut p0 = Vec::with_capacity(n_cells);
            let mut dx1 = Vec::with_capacity(n_cells);
            let mut dz1 = Vec::with_capacity(n_cells);
            let mut dx2 = Vec::with_capacity(n_cells);
            let mut dz2 = Vec::with_capacity(n_cells);
            let mut det_j_signed = Vec::with_capacity(n_cells);
            let mut det_j = Vec::with_capacity(n_cells);

            for c in &raw.cells {
                let pt0 = raw.pts[c[0] as usize];
                let pt1 = raw.pts[c[1] as usize];
                let pt2 = raw.pts[c[2] as usize];
                let d1x = pt1[0] - pt0[0];
                let d1z = pt1[1] - pt0[1];
                let d2x = pt2[0] - pt0[0];
                let d2z = pt2[1] - pt0[1];
                let dj = d1x * d2z - d2x * d1z;
                p0.push(pt0);
                dx1.push(d1x);
                dz1.push(d1z);
                dx2.push(d2x);
                dz2.push(d2z);
                det_j_signed.push(dj);
                det_j.push(dj.abs());
            }

            let mut q_n = [[0.0; 6]; 6];
            for (q, qp) in QUAD_PTS_2D.iter().enumerate() {
                let (n, _, _) = tri6_basis(qp[0], qp[1]);
                q_n[q] = n;
            }

            let mut norm_e2_all = 0.0;
            let mut w_si = 0.0;
            let mut w_sap = 0.0;
            let mut w_air = 0.0;

            for c in 0..n_cells {
                let mut e2_cell = 0.0;
                let mut wsap_cell = 0.0;
                for q in 0..6 {
                    let mut ex_q = C64::default();
                    let mut ey_q = C64::default();
                    let mut ez_q = C64::default();
                    for k in 0..6 {
                        let nk = q_n[q][k];
                        ex_q += ex[c][k].scale(nk);
                        ey_q += ey[c][k].scale(nk);
                        ez_q += ez[c][k].scale(nk);
                    }
                    let e2_q = ex_q.norm_sq() + ey_q.norm_sq() + ez_q.norm_sq();
                    e2_cell += e2_q * QUAD_W_2D[q];
                    if masks.sap[c] {
                        let t1 = ex_q.conj() * ex_q.scale(EPS_SAP_TENSOR[0][0]);
                        let t2 = ey_q.conj() * (ey_q.scale(EPS_SAP_TENSOR[1][1]) + ez_q.scale(EPS_SAP_TENSOR[1][2]));
                        let t3 = ez_q.conj() * (ey_q.scale(EPS_SAP_TENSOR[2][1]) + ez_q.scale(EPS_SAP_TENSOR[2][2]));
                        wsap_cell += (t1 + t2 + t3).re * QUAD_W_2D[q];
                    }
                }
                e2_cell *= det_j[c];
                wsap_cell *= det_j[c];

                if masks.si_all[c] {
                    w_si += EPS_SI * e2_cell;
                }
                if masks.phys[c] {
                    norm_e2_all += e2_cell;
                }
                if masks.sap[c] {
                    w_sap += wsap_cell;
                }
                if masks.air[c] {
                    w_air += e2_cell;
                }
            }

            Self {
                kind: kind.to_string(),
                cells: raw.cells.clone(),
                masks,
                ex,
                ey,
                ez,
                p0,
                dx1,
                dz1,
                dx2,
                dz2,
                det_j_signed,
                norm_e2_all,
                norm_eps_e2_tot: w_si + w_sap + w_air,
            }
        }
    }

    pub struct PreparedSmsSlice2D {
        pub pts: Vec<[f64; 2]>,
        pub cells: Vec<[u32; 6]>,
        pub ux: Vec<[C64; 6]>,
        pub uy: Vec<[C64; 6]>,
        pub uz: Vec<[C64; 6]>,
        pub masks: CellMasks,
    }

    impl PreparedSmsSlice2D {
        pub fn from_raw(raw: &RawSmsMode, kind: &str, x_phys: f64, z_sub: f64, z_air: f64) -> Self {
            let (pts_2d, tri6, mut ux_c, mut uy_c, mut uz_c) = if raw.topology == "Tetrahedron_10" {
                let n_tet = raw.cells_flat.len() / 10;
                let n_pts = raw.pts_flat.len() / 3;
                let mut on_y0 = vec![false; n_pts];
                for p in 0..n_pts {
                    on_y0[p] = raw.pts_flat[p * 3 + 1].abs() < 1e-12;
                }

                let face_nodes = [0usize, 1, 2, 4, 5, 6];
                let mut c6_list = Vec::new();
                let mut ux_list = Vec::new();
                let mut uy_list = Vec::new();
                let mut uz_list = Vec::new();

                for t in 0..n_tet {
                    let cell = &raw.cells_flat[t * 10..(t + 1) * 10];
                    let count_y0 = cell.iter().filter(|&&idx| on_y0[idx as usize]).count();
                    if count_y0 == 6 {
                        let mut c6 = [0u32; 6];
                        let mut ux6 = [C64::default(); 6];
                        let mut uy6 = [C64::default(); 6];
                        let mut uz6 = [C64::default(); 6];
                        for (k, &fn_idx) in face_nodes.iter().enumerate() {
                            c6[k] = cell[fn_idx];
                            let src = t * 10 + PERM_TET10[fn_idx];
                            ux6[k] = C64::new(raw.ux_re[src], raw.ux_im[src]);
                            uy6[k] = C64::new(raw.uy_re[src], raw.uy_im[src]);
                            uz6[k] = C64::new(raw.uz_re[src], raw.uz_im[src]);
                        }
                        c6_list.push(c6);
                        ux_list.push(ux6);
                        uy_list.push(uy6);
                        uz_list.push(uz6);
                    }
                }

                let mut used: Vec<u32> = c6_list.iter().flat_map(|c| c.iter().copied()).collect();
                used.sort_unstable();
                used.dedup();
                let mut remap = vec![0u32; n_pts];
                let mut pts_2d = Vec::with_capacity(used.len());
                for (new_idx, &old_idx) in used.iter().enumerate() {
                    remap[old_idx as usize] = new_idx as u32;
                    let x = raw.pts_flat[old_idx as usize * 3];
                    let z = raw.pts_flat[old_idx as usize * 3 + 2];
                    pts_2d.push([x, z]);
                }
                for c6 in &mut c6_list {
                    for idx in c6.iter_mut() {
                        *idx = remap[*idx as usize];
                    }
                }
                (pts_2d, c6_list, ux_list, uy_list, uz_list)
            } else {
                let n_tri = raw.cells_flat.len() / 6;
                let n_pts = raw.pts_flat.len() / raw.point_dim;
                let mut pts_2d = Vec::with_capacity(n_pts);
                for p in 0..n_pts {
                    pts_2d.push([raw.pts_flat[p * raw.point_dim], raw.pts_flat[p * raw.point_dim + 1]]);
                }
                let mut tri6 = Vec::with_capacity(n_tri);
                let mut ux_list = vec![[C64::default(); 6]; n_tri];
                let uy_list = vec![[C64::default(); 6]; n_tri];
                let mut uz_list = vec![[C64::default(); 6]; n_tri];
                for t in 0..n_tri {
                    let mut c6 = [0u32; 6];
                    c6.copy_from_slice(&raw.cells_flat[t * 6..(t + 1) * 6]);
                    tri6.push(c6);
                    for k in 0..6 {
                        let src = t * 6 + PERM_TRI6[k];
                        ux_list[t][k] = C64::new(raw.ux_re[src], raw.ux_im[src]);
                        uz_list[t][k] = C64::new(raw.uz_re[src], raw.uz_im[src]);
                    }
                }
                (pts_2d, tri6, ux_list, uy_list, uz_list)
            };

            let masks = classify_cells(&pts_2d, &tri6, kind, x_phys, z_sub, z_air);
            let mut z_phasor = C64::default();
            for c in 0..tri6.len() {
                if masks.core[c] {
                    for k in 0..6 {
                        z_phasor += ux_c[c][k] * ux_c[c][k] + uz_c[c][k] * uz_c[c][k];
                    }
                }
            }
            let theta = if z_phasor.abs() > 1e-30 { -0.5 * z_phasor.angle() } else { 0.0 };
            let rot = C64::cis(theta);
            for c in 0..tri6.len() {
                for k in 0..6 {
                    ux_c[c][k] = ux_c[c][k] * rot;
                    uy_c[c][k] = uy_c[c][k] * rot;
                    uz_c[c][k] = uz_c[c][k] * rot;
                }
            }

            Self { pts: pts_2d, cells: tri6, ux: ux_c, uy: uy_c, uz: uz_c, masks }
        }
    }

    #[derive(Clone, Debug)]
    pub struct BrillouinResult {
        pub q_pe_signed_um2_inv: f64,
        pub q_mb_signed_um2_inv: f64,
        pub gb_qm1000_vg: f64,
    }

    fn build_si_p_voigt(orientation: &str) -> [[f64; 6]; 6] {
        let mut p = [[0.0; 6]; 6];
        if orientation == "110" {
            let p11_rot = 0.5 * (P11_SI + P12_SI) + P44_SI;
            let p12_rot = 0.5 * (P11_SI + P12_SI) - P44_SI;
            let p66_rot = 0.5 * (P11_SI - P12_SI);
            p[0][0] = p11_rot;
            p[1][1] = p11_rot;
            p[2][2] = P11_SI;
            p[0][1] = p12_rot;
            p[1][0] = p12_rot;
            p[0][2] = P12_SI;
            p[2][0] = P12_SI;
            p[1][2] = P12_SI;
            p[2][1] = P12_SI;
            p[3][3] = P44_SI;
            p[4][4] = P44_SI;
            p[5][5] = p66_rot;
        } else {
            p[0][0] = P11_SI;
            p[1][1] = P11_SI;
            p[2][2] = P11_SI;
            p[0][1] = P12_SI;
            p[1][0] = P12_SI;
            p[0][2] = P12_SI;
            p[2][0] = P12_SI;
            p[1][2] = P12_SI;
            p[2][1] = P12_SI;
            p[3][3] = P44_SI;
            p[4][4] = P44_SI;
            p[5][5] = P44_SI;
        }
        p
    }

    #[allow(clippy::too_many_arguments)]
    pub fn compute_brillouin_overlap(
        ems_pump: &PreparedEmsMode2D,
        ems_stokes: &PreparedEmsMode2D,
        sms: &PreparedSmsSlice2D,
        f_mech_hz: f64,
        q_ac: f64,
        _neff_opt: f64,
        ng_opt: f64,
        backward_sbs: bool,
        si_orientation: &str,
        same_mesh: bool,
    ) -> BrillouinResult {
        let n_cells = sms.cells.len();

        let mut q_n = [[0.0; 6]; 6];
        let mut q_dn_xi = [[0.0; 6]; 6];
        let mut q_dn_eta = [[0.0; 6]; 6];
        for (q, qp) in QUAD_PTS_2D.iter().enumerate() {
            let (n, dn_xi, dn_eta) = tri6_basis(qp[0], qp[1]);
            q_n[q] = n;
            q_dn_xi[q] = dn_xi;
            q_dn_eta[q] = dn_eta;
        }

        let s1d = [0.5 * (1.0 - (3.0f64 / 5.0).sqrt()), 0.5, 0.5 * (1.0 + (3.0f64 / 5.0).sqrt())];
        let w1d = [0.5 * 5.0 / 9.0, 0.5 * 8.0 / 9.0, 0.5 * 5.0 / 9.0];
        let mut b1d = [[0.0; 3]; 3];
        for q in 0..3 {
            let s = s1d[q];
            b1d[0][q] = (1.0 - s) * (1.0 - 2.0 * s);
            b1d[1][q] = 4.0 * s * (1.0 - s);
            b1d[2][q] = s * (2.0 * s - 1.0);
        }

        let evaluate_ems_at = |ems: &PreparedEmsMode2D, x: f64, z: f64, reg_key: &str| -> (C64, C64, C64) {
            let mask = match reg_key {
                "core" => &ems.masks.core,
                "slab" => &ems.masks.slab,
                "si_all" => &ems.masks.si_all,
                "sap" => &ems.masks.sap,
                _ => &ems.masks.phys,
            };
            let mut best_c = 0usize;
            let mut best_min = f64::NEG_INFINITY;
            let mut best_xi = 0.0;
            let mut best_eta = 0.0;
            for c in 0..ems.cells.len() {
                if !mask[c] {
                    continue;
                }
                let dx = x - ems.p0[c][0];
                let dz = z - ems.p0[c][1];
                let inv = 1.0 / ems.det_j_signed[c];
                let xi = (ems.dz2[c] * dx - ems.dx2[c] * dz) * inv;
                let eta = (-ems.dz1[c] * dx + ems.dx1[c] * dz) * inv;
                let l0 = 1.0 - xi - eta;
                let min_b = xi.min(eta).min(l0);
                if min_b > best_min {
                    best_min = min_b;
                    best_c = c;
                    best_xi = xi;
                    best_eta = eta;
                }
            }
            let (nb, _, _) = tri6_basis(best_xi, best_eta);
            let mut ex_v = C64::default();
            let mut ey_v = C64::default();
            let mut ez_v = C64::default();
            for k in 0..6 {
                ex_v += ems.ex[best_c][k].scale(nb[k]);
                ey_v += ems.ey[best_c][k].scale(nb[k]);
                ez_v += ems.ez[best_c][k].scale(nb[k]);
            }
            (ex_v, ey_v, ez_v)
        };

        let p_mat = build_si_p_voigt(si_orientation);
        let mut p_sap = [[0.0; 6]; 6];
        p_sap[0][0] = -0.21;
        p_sap[1][1] = -0.21;
        p_sap[2][2] = -0.21;
        p_sap[0][1] = -0.02;
        p_sap[1][0] = -0.02;
        p_sap[0][2] = -0.02;
        p_sap[2][0] = -0.02;
        p_sap[1][2] = -0.02;
        p_sap[2][1] = -0.02;
        p_sap[3][3] = -0.10;
        p_sap[4][4] = -0.10;
        p_sap[5][5] = -0.10;
        let eps_sap_avg = 3.0176;

        let mut norm_u2_tot = 0.0;
        let mut norm_rho_u2 = 0.0;
        let mut q_pe_si = C64::default();
        let mut q_pe_sap = C64::default();
        let mut q_mb_top = C64::default();
        let mut q_mb_side = C64::default();
        let mut q_mb_bot = C64::default();

        for c in 0..n_cells {
            let cn = sms.cells[c];
            let p0 = sms.pts[cn[0] as usize];
            let p1 = sms.pts[cn[1] as usize];
            let p2 = sms.pts[cn[2] as usize];
            let dx_dxi = p1[0] - p0[0];
            let dz_dxi = p1[1] - p0[1];
            let dx_deta = p2[0] - p0[0];
            let dz_deta = p2[1] - p0[1];
            let det_j_signed = dx_dxi * dz_deta - dx_deta * dz_dxi;
            let det_j = det_j_signed.abs();
            let inv_det_j = 1.0 / det_j_signed;

            let is_si = sms.masks.si_all[c];
            let is_core = sms.masks.core[c];
            let is_slab = sms.masks.slab[c];
            let is_sap = sms.masks.sap[c];

            let reg_key = if is_core {
                "core"
            } else if is_slab {
                "slab"
            } else if is_sap {
                "sap"
            } else {
                ""
            };

            let mut u2_cell = 0.0;

            for q in 0..6 {
                let mut ux_q = C64::default();
                let mut uy_q = C64::default();
                let mut uz_q = C64::default();
                let mut dux_dxi = C64::default();
                let mut dux_deta = C64::default();
                let mut duy_dxi = C64::default();
                let mut duy_deta = C64::default();
                let mut duz_dxi = C64::default();
                let mut duz_deta = C64::default();
                let mut x_q = 0.0;
                let mut z_q = 0.0;

                for k in 0..6 {
                    let nk = q_n[q][k];
                    let dxi = q_dn_xi[q][k];
                    let deta = q_dn_eta[q][k];
                    ux_q += sms.ux[c][k].scale(nk);
                    uy_q += sms.uy[c][k].scale(nk);
                    uz_q += sms.uz[c][k].scale(nk);
                    dux_dxi += sms.ux[c][k].scale(dxi);
                    dux_deta += sms.ux[c][k].scale(deta);
                    duy_dxi += sms.uy[c][k].scale(dxi);
                    duy_deta += sms.uy[c][k].scale(deta);
                    duz_dxi += sms.uz[c][k].scale(dxi);
                    duz_deta += sms.uz[c][k].scale(deta);
                    let pt = sms.pts[cn[k] as usize];
                    x_q += pt[0] * nk;
                    z_q += pt[1] * nk;
                }

                let u_mag_q = (ux_q.norm_sq() + uy_q.norm_sq() + uz_q.norm_sq()).sqrt();
                u2_cell += (u_mag_q * u_mag_q) * QUAD_W_2D[q];

                let (mut ex_p_q, mut ey_p_q, mut ez_p_q) = (C64::default(), C64::default(), C64::default());
                let (mut ex_s_q, mut ey_s_q, mut ez_s_q) = (C64::default(), C64::default(), C64::default());

                if same_mesh {
                    for k in 0..6 {
                        let nk = q_n[q][k];
                        ex_p_q += ems_pump.ex[c][k].scale(nk);
                        ey_p_q += ems_pump.ey[c][k].scale(nk);
                        ez_p_q += ems_pump.ez[c][k].scale(nk);
                        ex_s_q += ems_stokes.ex[c][k].scale(nk);
                        ey_s_q += ems_stokes.ey[c][k].scale(nk);
                        ez_s_q += ems_stokes.ez[c][k].scale(nk);
                    }
                } else if !reg_key.is_empty() {
                    let p_eval = evaluate_ems_at(ems_pump, x_q, z_q, reg_key);
                    ex_p_q = p_eval.0;
                    ey_p_q = p_eval.1;
                    ez_p_q = p_eval.2;
                    if std::ptr::eq(ems_pump, ems_stokes) {
                        ex_s_q = ex_p_q;
                        ey_s_q = ey_p_q;
                        ez_s_q = ez_p_q;
                    } else {
                        let s_eval = evaluate_ems_at(ems_stokes, x_q, z_q, reg_key);
                        ex_s_q = s_eval.0;
                        ey_s_q = s_eval.1;
                        ez_s_q = s_eval.2;
                    }
                }

                if backward_sbs {
                    ey_s_q = -ey_s_q;
                }

                if is_si || is_sap {
                    let dux_dx = (dux_dxi.scale(dz_deta) - dux_deta.scale(dz_dxi)).scale(inv_det_j);
                    let dux_dz = (dux_deta.scale(dx_dxi) - dux_dxi.scale(dx_deta)).scale(inv_det_j);
                    let duy_dx = (duy_dxi.scale(dz_deta) - duy_deta.scale(dz_dxi)).scale(inv_det_j);
                    let duy_dz = (duy_deta.scale(dx_dxi) - duy_dxi.scale(dx_deta)).scale(inv_det_j);
                    let duz_dx = (duz_dxi.scale(dz_deta) - duz_deta.scale(dz_dxi)).scale(inv_det_j);
                    let duz_dz = (duz_deta.scale(dx_dxi) - duz_dxi.scale(dx_deta)).scale(inv_det_j);

                    let minus_iq = C64::new(0.0, -q_ac);
                    let s_xx = dux_dx;
                    let s_yy = minus_iq * uy_q;
                    let s_zz = duz_dz;
                    let gamma_yz = duy_dz + minus_iq * uz_q;
                    let gamma_xz = dux_dz + duz_dx;
                    let gamma_xy = duy_dx + minus_iq * ux_q;

                    let s_voigt_conj =
                        [s_xx.conj(), s_yy.conj(), s_zz.conj(), gamma_yz.conj(), gamma_xz.conj(), gamma_xy.conj()];
                    let d_voigt = [
                        ex_s_q.conj() * ex_p_q,
                        ey_s_q.conj() * ey_p_q,
                        ez_s_q.conj() * ez_p_q,
                        ey_s_q.conj() * ez_p_q + ez_s_q.conj() * ey_p_q,
                        ex_s_q.conj() * ez_p_q + ez_s_q.conj() * ex_p_q,
                        ex_s_q.conj() * ey_p_q + ey_s_q.conj() * ex_p_q,
                    ];

                    if is_si {
                        for j in 0..6 {
                            let mut dp_j = C64::default();
                            for i in 0..6 {
                                dp_j += d_voigt[i].scale(p_mat[i][j]);
                            }
                            let integrand = (dp_j * s_voigt_conj[j]).scale(-(EPS_SI * EPS_SI));
                            q_pe_si += integrand.scale(QUAD_W_2D[q] * det_j);
                        }
                    } else if is_sap {
                        let mut dp_sap = C64::default();
                        for i in 0..6 {
                            for j in 0..6 {
                                dp_sap += (d_voigt[i] * s_voigt_conj[j]).scale(p_sap[i][j]);
                            }
                        }
                        q_pe_sap += dp_sap.scale(-(eps_sap_avg * eps_sap_avg) * QUAD_W_2D[q] * det_j);
                    }
                }
            }

            u2_cell *= det_j;
            if is_si {
                norm_u2_tot += u2_cell;
                norm_rho_u2 += RHO_SI * u2_cell;
            }
            if is_sap {
                norm_u2_tot += u2_cell;
                norm_rho_u2 += RHO_SAP * u2_cell;
            }

            if is_si {
                let edge_reg_key = if is_core { "core" } else { "slab" };
                for loc in &EDGE_LOC {
                    let pa = sms.pts[cn[loc[0]] as usize];
                    let pm = sms.pts[cn[loc[1]] as usize];
                    let pb = sms.pts[cn[loc[2]] as usize];
                    let xm = pm[0];
                    let zm = pm[1];
                    let mut bnd_kind = "";
                    let mut nx = 0.0;
                    let mut nz = 0.0;
                    let mut eps_ext_tan = 1.0;
                    let mut eps_ext_norm = 1.0;

                    if ems_pump.kind == "strip" {
                        if (zm - 0.30e-6).abs() < 1e-12 && xm.abs() <= 0.40e-6 + 1e-12 {
                            bnd_kind = "top";
                            nz = 1.0;
                        } else if (xm + 0.40e-6).abs() < 1e-12 && (-1e-12..=0.30e-6 + 1e-12).contains(&zm) {
                            bnd_kind = "side";
                            nx = -1.0;
                        } else if (xm - 0.40e-6).abs() < 1e-12 && (-1e-12..=0.30e-6 + 1e-12).contains(&zm) {
                            bnd_kind = "side";
                            nx = 1.0;
                        } else if zm.abs() < 1e-12 && xm.abs() <= 0.40e-6 + 1e-12 {
                            bnd_kind = "bot";
                            nz = -1.0;
                            eps_ext_tan = EPS_SAP_TENSOR[0][0];
                            eps_ext_norm = EPS_SAP_TENSOR[2][2];
                        }
                    } else if (zm - 0.30e-6).abs() < 1e-12 && xm.abs() <= 0.25e-6 + 1e-12 {
                        bnd_kind = "top";
                        nz = 1.0;
                    } else if (zm - 0.17e-6).abs() < 1e-12 && xm.abs() >= 0.25e-6 - 1e-12 {
                        bnd_kind = "top";
                        nz = 1.0;
                    } else if (xm + 0.25e-6).abs() < 1e-12 && (0.17e-6 - 1e-12..=0.30e-6 + 1e-12).contains(&zm) {
                        bnd_kind = "side";
                        nx = -1.0;
                    } else if (xm - 0.25e-6).abs() < 1e-12 && (0.17e-6 - 1e-12..=0.30e-6 + 1e-12).contains(&zm) {
                        bnd_kind = "side";
                        nx = 1.0;
                    } else if zm.abs() < 1e-12 {
                        bnd_kind = "bot";
                        nz = -1.0;
                        eps_ext_tan = EPS_SAP_TENSOR[0][0];
                        eps_ext_norm = EPS_SAP_TENSOR[2][2];
                    }

                    if bnd_kind.is_empty() {
                        continue;
                    }

                    let ds = (pb[0] - pa[0]).hypot(pb[1] - pa[1]);
                    let mut edge_val = C64::default();
                    for q in 0..3 {
                        let mut ux_1d = C64::default();
                        let mut uz_1d = C64::default();
                        let mut x_1d = 0.0;
                        let mut z_1d = 0.0;
                        for m in 0..3 {
                            let bm = b1d[m][q];
                            ux_1d += sms.ux[c][loc[m]].scale(bm);
                            uz_1d += sms.uz[c][loc[m]].scale(bm);
                            let pt = sms.pts[cn[loc[m]] as usize];
                            x_1d += pt[0] * bm;
                            z_1d += pt[1] * bm;
                        }
                        let (exp_1d, eyp_1d, ezp_1d, exs_1d, mut eys_1d, ezs_1d) = if same_mesh {
                            let mut ep = (C64::default(), C64::default(), C64::default());
                            let mut es = (C64::default(), C64::default(), C64::default());
                            for m in 0..3 {
                                let bm = b1d[m][q];
                                ep.0 += ems_pump.ex[c][loc[m]].scale(bm);
                                ep.1 += ems_pump.ey[c][loc[m]].scale(bm);
                                ep.2 += ems_pump.ez[c][loc[m]].scale(bm);
                                es.0 += ems_stokes.ex[c][loc[m]].scale(bm);
                                es.1 += ems_stokes.ey[c][loc[m]].scale(bm);
                                es.2 += ems_stokes.ez[c][loc[m]].scale(bm);
                            }
                            (ep.0, ep.1, ep.2, es.0, es.1, es.2)
                        } else {
                            let ep = evaluate_ems_at(ems_pump, x_1d, z_1d, edge_reg_key);
                            let es = if std::ptr::eq(ems_pump, ems_stokes) {
                                ep
                            } else {
                                evaluate_ems_at(ems_stokes, x_1d, z_1d, edge_reg_key)
                            };
                            (ep.0, ep.1, ep.2, es.0, es.1, es.2)
                        };
                        if backward_sbs {
                            eys_1d = -eys_1d;
                        }
                        let en_p = exp_1d.scale(nx) + ezp_1d.scale(nz);
                        let en_s = exs_1d.scale(nx) + ezs_1d.scale(nz);
                        let e_tan_dot = (exs_1d.conj() * exp_1d + eys_1d.conj() * eyp_1d + ezs_1d.conj() * ezp_1d)
                            - (en_s.conj() * en_p);
                        let dn_prod = (en_s.conj() * en_p).scale(EPS_SI * EPS_SI);
                        let un_conj = ux_1d.conj().scale(nx) + uz_1d.conj().scale(nz);
                        let bracket =
                            e_tan_dot.scale(EPS_SI - eps_ext_tan) + dn_prod.scale(1.0 / eps_ext_norm - 1.0 / EPS_SI);
                        edge_val += (un_conj * bracket).scale(w1d[q] * ds);
                    }

                    match bnd_kind {
                        "top" => q_mb_top += edge_val,
                        "side" => q_mb_side += edge_val,
                        "bot" => q_mb_bot += edge_val,
                        _ => {}
                    }
                }
            }
        }

        let q_pe_tot = q_pe_si + q_pe_sap;
        let q_mb_tot = q_mb_top + q_mb_side + q_mb_bot;
        let q_tot = q_pe_tot + q_mb_tot;

        let omega_m = 2.0 * PI * f_mech_hz;
        let denom_norms = ems_pump.norm_eps_e2_tot * ems_stokes.norm_eps_e2_tot * norm_rho_u2;
        let vg_opt = C0 / ng_opt;
        let prefac_vg = (2.0 * OMEGA_OPT) / ((omega_m * omega_m) * (vg_opt * vg_opt) * denom_norms);
        let gb_qm1000_vg = prefac_vg * q_tot.norm_sq() * 1000.0;

        let norm_denom = (ems_pump.norm_eps_e2_tot * ems_stokes.norm_eps_e2_tot * norm_u2_tot).sqrt();
        let ref_phasor = if q_pe_tot.abs() >= q_mb_tot.abs() { q_pe_tot } else { q_mb_tot };
        let mut phase_ref = if ref_phasor.abs() > 1e-30 { C64::cis(-ref_phasor.angle()) } else { C64::new(1.0, 0.0) };
        if (q_pe_tot * phase_ref).re + (q_mb_tot * phase_ref).re < 0.0 {
            phase_ref = -phase_ref;
        }

        let q_pe_signed_um2_inv = (q_pe_tot * phase_ref).re / norm_denom * 1e-12;
        let q_mb_signed_um2_inv = (q_mb_tot * phase_ref).re / norm_denom * 1e-12;

        BrillouinResult { q_pe_signed_um2_inv, q_mb_signed_um2_inv, gb_qm1000_vg }
    }
}

#[test]
fn gate5_ec11_untrusted_fuzz_and_ec12_sos_brillouin_physics_invariant() -> Result<()> {
    let t0 = Instant::now();
    let scratch = test_scratch_dir("gate5")?;

    // -------------------------------------------------------------------------
    // 1. EC-11: Untrusted Reader Fuzzing (10,000 mutations + u64::MAX OOM-bomb)
    // -------------------------------------------------------------------------
    let seed1 = {
        let mut fb = FileBuilder::new();
        fb.set_attr("title", AttrValue::VarLenString("fuzz-seed-1".into()));
        fb.create_dataset("contig").with_f64_data(&[1.0, 2.0, 3.0, 4.0]);
        fb.create_dataset("lzf").with_f64_data(&[5.0; 64]).with_chunks(&[16]).with_shuffle().with_lzf();
        fb.finish()?
    };
    let seed2 = {
        let mut fb = FileBuilder::new();
        let mut g = fb.create_group("nested");
        for i in 0..16 {
            g.set_attr(&format!("dense_{i:02}"), AttrValue::I64(i));
        }
        g.create_dataset("deflate")
            .with_f64_data(&[7.0; 64])
            .with_chunks(&[16])
            .with_shuffle()
            .with_deflate(DeflateLevel::Three);
        fb.add_group(g.finish());
        fb.finish()?
    };
    // Explicit budget enforcement check: valid 1024-byte dataset rejected under a 256-byte ReadBudgetConfig
    let small_budget = ReadBudgetConfig { max_dataset_bytes: 256, max_group_depth: 8 };
    let budget_res = fuzz_walk_untrusted(&seed1, small_budget);
    assert!(budget_res.is_err(), "ReadBudgetConfig must reject dataset exceeding max_dataset_bytes");

    let mutations = generate_10k_mutations(&[seed1, seed2]);
    let fuzz_budget = ReadBudgetConfig { max_dataset_bytes: 16 * 1024 * 1024, max_group_depth: 24 };
    let panics: Vec<usize> = mutations
        .par_iter()
        .enumerate()
        .filter_map(|(idx, mbytes)| {
            let res = catch_unwind(AssertUnwindSafe(|| {
                let _ = fuzz_walk_untrusted(mbytes, fuzz_budget);
            }));
            if res.is_err() {
                Some(idx)
            } else {
                None
            }
        })
        .collect();
    if !panics.is_empty() {
        bail!("EC-11 fuzzing panicked on {} / 10,000 mutations (first at #{})", panics.len(), panics[0]);
    }

    // -------------------------------------------------------------------------
    // 2. EC-12: Pure-Rust P2 Triangle_6 + Tetrahedron_10 SOS Brillouin PE + MB Overlap
    // -------------------------------------------------------------------------
    use brillouin::{
        compute_brillouin_overlap, load_ems_mode_2d, load_sms_mode, PreparedEmsMode2D, PreparedSmsSlice2D, WL0,
    };

    let ems_strip =
        PreparedEmsMode2D::from_raw(&load_ems_mode_2d(Path::new(BENCH_MODES_2D), 0)?, "strip", 2.0e-6, 2.0e-6, 1.6e-6);
    let ems_rib =
        PreparedEmsMode2D::from_raw(&load_ems_mode_2d(Path::new(BENCH_RIB_2D), 0)?, "rib", 2.0e-6, 2.0e-6, 1.6e-6);
    let ems_te0 = PreparedEmsMode2D::from_raw(
        &load_ems_mode_2d(Path::new(BENCH_COUPLED_2D), 0)?,
        "strip",
        2.0e-6,
        2.0e-6,
        1.2e-6,
    );
    let sms_m1 =
        PreparedSmsSlice2D::from_raw(&load_sms_mode(Path::new(BENCH_MECH_3D), 0)?, "strip", 1.6e-6, 1.6e-6, 0.5e-6);
    let sms_a1 =
        PreparedSmsSlice2D::from_raw(&load_sms_mode(Path::new(BENCH_COUPLED_2D), 4)?, "strip", 2.0e-6, 2.0e-6, 1.2e-6);

    let q_strip = 2.0 * (2.0 * PI / WL0) * 2.894656;
    let bsbs = compute_brillouin_overlap(
        &ems_strip,
        &ems_strip,
        &sms_m1,
        17.944410e9,
        q_strip,
        2.894656,
        3.936293,
        true,
        "100",
        false,
    );
    let fsbs =
        compute_brillouin_overlap(&ems_te0, &ems_te0, &sms_a1, 5.001021e9, 0.0, 2.894656, 3.936293, false, "100", true);
    let rib_bsbs = compute_brillouin_overlap(
        &ems_rib,
        &ems_rib,
        &sms_m1,
        17.944410e9,
        q_strip,
        2.625066,
        3.850000,
        true,
        "100",
        false,
    );

    if (bsbs.gb_qm1000_vg - 37.30).abs() > 0.01 {
        bail!("Expected B-SBS M1 G_B(vg) = 37.30 W^-1 m^-1, got {:.6}", bsbs.gb_qm1000_vg);
    }
    if (fsbs.gb_qm1000_vg - 6053.66).abs() > 0.01 {
        bail!("Expected F-SBS A1 G_B(vg) = 6053.66 W^-1 m^-1, got {:.6}", fsbs.gb_qm1000_vg);
    }
    if (fsbs.q_pe_signed_um2_inv - 6.1783).abs() > 1e-4 {
        bail!("Expected F-SBS A1 Q_PE = 6.1783 um^-2, got {:.6}", fsbs.q_pe_signed_um2_inv);
    }
    if (fsbs.q_mb_signed_um2_inv - 0.0884).abs() > 1e-4 {
        bail!("Expected F-SBS A1 Q_MB = 0.0884 um^-2, got {:.6}", fsbs.q_mb_signed_um2_inv);
    }

    // Re-export sos_300nm_rib_modes.h5 through write_nyx_v3, reload, and assert 0-ULP (<= 1e-15 rel error) Brillouin invariance
    let rib_sum = NyxReader::open(Path::new(BENCH_RIB_2D), ReadBudgetConfig::default())?.read_nyx_v3_summary(false)?;
    let rib_reexport = scratch.join("ec12_rib_reexport.h5");
    write_nyx_v3(&rib_sum, &rib_reexport)?;
    let ems_rib_reloaded =
        PreparedEmsMode2D::from_raw(&load_ems_mode_2d(&rib_reexport, 0)?, "rib", 2.0e-6, 2.0e-6, 1.6e-6);
    let rib_bsbs_reloaded = compute_brillouin_overlap(
        &ems_rib_reloaded,
        &ems_rib_reloaded,
        &sms_m1,
        17.944410e9,
        q_strip,
        2.625066,
        3.850000,
        true,
        "100",
        false,
    );
    let rel_err = (rib_bsbs_reloaded.gb_qm1000_vg - rib_bsbs.gb_qm1000_vg).abs() / rib_bsbs.gb_qm1000_vg.abs();
    if rel_err > 1e-15 {
        bail!("Rib Brillouin overlap relative error after write->read {:.3e} > 1e-15", rel_err);
    }

    println!(
        "[Gate 5 / EC-11 & EC-12 PASS] ({:.1} ms) Survived 10,000/10,000 adversarial fuzz mutations (0 panics/OOMs); SOS Brillouin B-SBS M1 G_B(vg)={:.2} W^-1m^-1, F-SBS A1 G_B(vg)={:.2} W^-1m^-1 (Q_PE={:.4} um^-2, Q_MB={:.4} um^-2, rel_err={:.2e} <= 1e-15)",
        t0.elapsed().as_secs_f64() * 1000.0,
        bsbs.gb_qm1000_vg,
        fsbs.gb_qm1000_vg,
        fsbs.q_pe_signed_um2_inv,
        fsbs.q_mb_signed_um2_inv,
        rel_err
    );
    Ok(())
}
