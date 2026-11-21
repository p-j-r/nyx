// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-Rust HDF5 writers for packed 5D (`L2 Shuffle+LZF`, `L3 Contiguous`, `L1 Deflate-3`),
//! incremental mode-by-mode streaming (`File::open_rw` + `append_staged`), and legacy `nyx-results` v3.

use crate::reader::{attr_i64, attr_str, NyxV3Summary};
use crate::sweep::{atomic_rename, staging_path_for};
use crate::xdmf::viz_link_name;
use crate::{FieldStorageLayout, MeshData, ModeHeader, PackedComponentMeta};
use anyhow::{anyhow, bail, Result};
use hdf5_pure::{
    make_f64_type, make_i32_type, make_u32_type, AttrValue, DeflateLevel, File, FileBuilder, GroupBuilder, MaxExtent,
};
use std::path::Path;

fn write_mesh_group(fb: &mut FileBuilder, mesh: &MeshData) -> Result<()> {
    if mesh.edge_node_pairs.len() != mesh.edge_node_points.len() {
        bail!(
            "MeshData curved edge mismatch: {} pairs vs {} points",
            mesh.edge_node_pairs.len(),
            mesh.edge_node_points.len()
        );
    }

    let expected_cells =
        mesh.cells_shape[0].checked_mul(mesh.cells_shape[1]).ok_or_else(|| anyhow!("mesh.cells_shape overflow"))?;
    if mesh.cells.len() != expected_cells {
        bail!("mesh.cells len {} != shape {:?} product {}", mesh.cells.len(), mesh.cells_shape, expected_cells);
    }

    let expected_facets =
        mesh.facets_shape[0].checked_mul(mesh.facets_shape[1]).ok_or_else(|| anyhow!("mesh.facets_shape overflow"))?;
    if mesh.facets.len() != expected_facets {
        bail!("mesh.facets len {} != shape {:?} product {}", mesh.facets.len(), mesh.facets_shape, expected_facets);
    }

    let expected_points =
        mesh.points_shape[0].checked_mul(mesh.points_shape[1]).ok_or_else(|| anyhow!("mesh.points_shape overflow"))?;
    if mesh.points.len() != expected_points {
        bail!("mesh.points len {} != shape {:?} product {}", mesh.points.len(), mesh.points_shape, expected_points);
    }

    let mut mg = fb.create_group("mesh");
    mg.set_attr("cell_topology", AttrValue::VarLenString(mesh.cell_topology.clone()));
    mg.set_attr("element_family", AttrValue::VarLenString(mesh.element_family.clone()));
    mg.set_attr("element_degree", AttrValue::I64(mesh.element_degree));
    mg.set_attr("space_dimension", AttrValue::I64(mesh.space_dimension));

    put_i32_1d(&mut mg, "cell_labels", &mesh.cell_labels);
    put_u32_2d(&mut mg, "cells", &mesh.cells, mesh.cells_shape);
    put_i32_1d(&mut mg, "facet_labels", &mesh.facet_labels);
    put_u32_2d(&mut mg, "facets", &mesh.facets, mesh.facets_shape);
    put_f64_2d(&mut mg, "points", &mesh.points, mesh.points_shape);

    if !mesh.edge_node_pairs.is_empty() {
        let n_curved = mesh.edge_node_pairs.len();
        let pairs_flat: Vec<u32> = mesh.edge_node_pairs.iter().flat_map(|p| p.iter().copied()).collect();
        let pts_flat: Vec<f64> = mesh.edge_node_points.iter().flat_map(|p| p.iter().copied()).collect();
        put_u32_2d(&mut mg, "edge_node_pairs", &pairs_flat, [n_curved, 2]);
        put_f64_2d(&mut mg, "edge_node_points", &pts_flat, [n_curved, 3]);
    }

    fb.add_group(mg.finish());
    Ok(())
}

fn put_i32_1d(grp: &mut GroupBuilder, name: &str, data: &[i32]) {
    let ds = grp.create_dataset(name);
    if data.is_empty() {
        ds.with_dtype(make_i32_type()).with_shape(&[0]);
    } else {
        ds.with_i32_data(data).with_shape(&[data.len() as u64]);
        if data.len() >= 4096 {
            ds.with_chunks(&[data.len().min(32768) as u64]).with_shuffle().with_lzf();
        }
    }
}

fn put_u32_2d(grp: &mut GroupBuilder, name: &str, data: &[u32], shape: [usize; 2]) {
    let ds = grp.create_dataset(name);
    let r = shape[0] as u64;
    let c = shape[1] as u64;
    if data.is_empty() || r == 0 || c == 0 {
        ds.with_dtype(make_u32_type()).with_shape(&[r, c]);
    } else {
        ds.with_u32_data(data).with_shape(&[r, c]);
        if data.len() >= 4096 {
            let chunk_r = (32768 / shape[1].max(1)).max(1).min(shape[0]) as u64;
            ds.with_chunks(&[chunk_r, c]).with_shuffle().with_lzf();
        }
    }
}

fn put_f64_2d(grp: &mut GroupBuilder, name: &str, data: &[f64], shape: [usize; 2]) {
    let ds = grp.create_dataset(name);
    let r = shape[0] as u64;
    let c = shape[1] as u64;
    if data.is_empty() || r == 0 || c == 0 {
        ds.with_dtype(make_f64_type()).with_shape(&[r, c]);
    } else {
        ds.with_f64_data(data).with_shape(&[r, c]);
        if data.len() >= 4096 {
            let chunk_r = (32768 / shape[1].max(1)).max(1).min(shape[0]) as u64;
            ds.with_chunks(&[chunk_r, c]).with_shuffle().with_lzf();
        }
    }
}

/// Writes a complete `nyx-results` file with packed 5D `cell_fields_packed [M, C, N_c, N_v, 2]`,
/// `eigenvalues [M, 2]`, `cell_fields_meta` JSON, and optional curved mesh edges,
/// atomically replacing `target_path` via `<target>.tmp` -> `std::fs::rename`.
#[allow(clippy::too_many_arguments)]
pub fn write_nyx_packed(
    target_path: &Path,
    project_json: &[u8],
    mesh: &MeshData,
    physics: &str,
    kind: &str,
    wavelength: Option<f64>,
    headers: &[ModeHeader],
    components: &[PackedComponentMeta],
    packed_cell_fields: &[f64],
    layout: FieldStorageLayout,
) -> Result<()> {
    let m = headers.len();
    let c = components.len();
    let n_c = mesh.cells_shape[0];
    let n_v = mesh.cells_shape[1];
    let expected_len = m
        .checked_mul(c)
        .and_then(|x| x.checked_mul(n_c))
        .and_then(|x| x.checked_mul(n_v))
        .and_then(|x| x.checked_mul(2))
        .ok_or_else(|| anyhow!("packed_cell_fields shape overflow"))?;
    if packed_cell_fields.len() != expected_len {
        bail!(
            "packed_cell_fields length {} != expected [M={}, C={}, N_c={}, N_v={}, 2] = {}",
            packed_cell_fields.len(),
            m,
            c,
            n_c,
            n_v,
            expected_len
        );
    }

    let mut fb = FileBuilder::new();
    fb.set_attr("format", AttrValue::VarLenString("nyx-results".into()));
    fb.set_attr("version", AttrValue::I64(4));
    fb.set_attr("writer", AttrValue::VarLenString("nyx".into()));
    fb.set_attr("processes", AttrValue::I64(1));
    fb.create_dataset("project").with_u8_data(project_json);

    write_mesh_group(&mut fb, mesh)?;

    let mut res = fb.create_group("results");
    let mut pg = res.create_group(physics);
    pg.set_attr("physics", AttrValue::VarLenString(physics.to_string()));
    pg.set_attr("n_solutions", AttrValue::I64(m as i64));
    pg.set_attr("kind", AttrValue::VarLenString(kind.to_string()));
    pg.set_attr("field_layout", AttrValue::VarLenString(layout.as_str().to_string()));
    if let Some(wl) = wavelength {
        pg.set_attr("wavelength", AttrValue::F64(wl));
    }

    let labels: Vec<&str> = headers.iter().map(|h| h.label.as_str()).collect();
    let labels_json = serde_json::to_string(&labels)?;
    pg.set_attr("mode_labels_json", AttrValue::VarLenString(labels_json));

    let meta_json = serde_json::to_vec(components)?;
    pg.create_dataset("cell_fields_meta").with_u8_data(&meta_json);

    let mut eig_flat = Vec::with_capacity(m * 2);
    for h in headers {
        eig_flat.push(h.eig_re);
        eig_flat.push(h.eig_im);
    }
    let ds_eig = pg.create_dataset("eigenvalues");
    if m == 0 {
        ds_eig.with_dtype(make_f64_type()).with_shape(&[0, 2]);
    } else {
        ds_eig.with_f64_data(&eig_flat).with_shape(&[m as u64, 2]);
    }

    let shape_5d = [m as u64, c as u64, n_c as u64, n_v as u64, 2];
    let ds_cf = pg.create_dataset("cell_fields_packed");
    if expected_len == 0 {
        ds_cf.with_dtype(make_f64_type()).with_shape(&shape_5d);
    } else {
        ds_cf.with_f64_data(packed_cell_fields).with_shape(&shape_5d);
        let chunk_5d = [1, c as u64, n_c as u64, n_v as u64, 2];
        match layout {
            FieldStorageLayout::L2ShuffleLzfPacked => {
                ds_cf.with_chunks(&chunk_5d).with_shuffle().with_lzf();
            }
            FieldStorageLayout::L1LegacyGzip3 => {
                ds_cf.with_chunks(&chunk_5d).with_shuffle().with_deflate(DeflateLevel::Three);
            }
            FieldStorageLayout::L3ContiguousPacked => {}
        }
    }

    res.add_group(pg.finish());
    fb.add_group(res.finish());

    if let Some(parent) = target_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let tmp_path = staging_path_for(target_path)?;
    fb.write(&tmp_path)?;
    atomic_rename(&tmp_path, target_path)?;
    Ok(())
}

/// Initializes a streaming `nyx-results` file (`n_solutions = 0`) with extensible chunked datasets
/// ready for incremental mode-by-mode appends via [`append_solved_mode_in_place`].
#[allow(clippy::too_many_arguments)]
pub fn init_nyx_packed_streaming(
    target_path: &Path,
    project_json: &[u8],
    mesh: &MeshData,
    physics: &str,
    kind: &str,
    wavelength: Option<f64>,
    components: &[PackedComponentMeta],
    layout: FieldStorageLayout,
) -> Result<()> {
    let c = components.len() as u64;
    let n_c = mesh.cells_shape[0] as u64;
    let n_v = mesh.cells_shape[1] as u64;
    let mode_elems = c
        .checked_mul(n_c)
        .and_then(|x| x.checked_mul(n_v))
        .and_then(|x| x.checked_mul(2))
        .ok_or_else(|| anyhow!("streaming mode_elems overflow"))?;

    let mut fb = FileBuilder::new();
    fb.set_attr("format", AttrValue::VarLenString("nyx-results".into()));
    fb.set_attr("version", AttrValue::I64(4));
    fb.set_attr("writer", AttrValue::VarLenString("nyx".into()));
    fb.set_attr("processes", AttrValue::I64(1));
    fb.create_dataset("project").with_u8_data(project_json);

    write_mesh_group(&mut fb, mesh)?;

    let mut res = fb.create_group("results");
    let mut pg = res.create_group(physics);
    pg.set_attr("physics", AttrValue::VarLenString(physics.to_string()));
    pg.set_attr("n_solutions", AttrValue::I64(0));
    pg.set_attr("kind", AttrValue::VarLenString(kind.to_string()));
    pg.set_attr("field_layout", AttrValue::VarLenString(layout.as_str().to_string()));
    pg.set_attr("mode_labels_json", AttrValue::VarLenString("[]".to_string()));
    if let Some(wl) = wavelength {
        pg.set_attr("wavelength", AttrValue::F64(wl));
    }

    let meta_json = serde_json::to_vec(components)?;
    pg.create_dataset("cell_fields_meta").with_u8_data(&meta_json);

    // `hdf5-pure 0.48.0` supports `append_staged` on rank-1 unlimited chunked datasets;
    // `@logical_shape_tail` records the trailing dimensions so `NyxReader` slices one mode chunk identically.
    pg.create_dataset("eigenvalues")
        .with_dtype(make_f64_type())
        .with_shape(&[0])
        .with_maxshape(&[MaxExtent::Unlimited])
        .with_chunks(&[2])
        .set_attr("logical_shape_tail", AttrValue::I64Array(vec![2]));

    let ds_cf = pg
        .create_dataset("cell_fields_packed")
        .with_dtype(make_f64_type())
        .with_shape(&[0])
        .with_maxshape(&[MaxExtent::Unlimited])
        .with_chunks(&[mode_elems.max(1)])
        .set_attr("logical_shape_tail", AttrValue::I64Array(vec![c as i64, n_c as i64, n_v as i64, 2]));

    match layout {
        FieldStorageLayout::L2ShuffleLzfPacked => {
            ds_cf.with_shuffle().with_lzf();
        }
        FieldStorageLayout::L1LegacyGzip3 => {
            ds_cf.with_shuffle().with_deflate(DeflateLevel::Three);
        }
        FieldStorageLayout::L3ContiguousPacked => {}
    }

    res.add_group(pg.finish());
    fb.add_group(res.finish());

    if let Some(parent) = target_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let tmp_path = staging_path_for(target_path)?;
    fb.write(&tmp_path)?;
    atomic_rename(&tmp_path, target_path)?;
    Ok(())
}

/// Appends one solved mode in place (`eigenvalues`, `cell_fields_packed`, `@n_solutions`, `@mode_labels_json`)
/// via `File::open_rw` + `append_staged` + `commit()`, without re-encoding previously written mode chunks.
pub fn append_solved_mode_in_place(
    path: &Path,
    physics: &str,
    header: &ModeHeader,
    one_mode_slice: &[f64],
) -> Result<()> {
    let phys_path = format!("results/{physics}");
    let eig_path = format!("{phys_path}/eigenvalues");
    let cf_path = format!("{phys_path}/cell_fields_packed");

    let session = File::open_rw(path)?;
    {
        let phys_grp = session.group(&phys_path)?;
        let pattrs = phys_grp.attrs()?;
        let cur_n = attr_i64(&pattrs, "n_solutions")?.max(0) as usize;
        let mut labels: Vec<String> =
            attr_str(&pattrs, "mode_labels_json").ok().and_then(|j| serde_json::from_str(&j).ok()).unwrap_or_default();
        while labels.len() < cur_n {
            labels.push(format!("mode_{}", labels.len()));
        }
        labels.push(header.label.clone());

        let mut ds_cf = session.dataset(&cf_path)?;
        let cf_attrs = ds_cf.attrs()?;
        let tail = cf_attrs
            .get("logical_shape_tail")
            .and_then(AttrValue::to_i64s)
            .ok_or_else(|| anyhow!("{} is not an extensible streaming dataset", cf_path))?;
        let expected_elems = tail
            .iter()
            .try_fold(1usize, |acc, &d| acc.checked_mul(d.max(0) as usize))
            .ok_or_else(|| anyhow!("logical_shape_tail overflow"))?;
        if one_mode_slice.len() != expected_elems {
            bail!(
                "append_solved_mode_in_place slice len {} != expected mode_elems {}",
                one_mode_slice.len(),
                expected_elems
            );
        }

        let mut ds_eig = session.dataset(&eig_path)?;
        ds_eig.append_staged(|b| {
            b.append_f64(&[header.eig_re, header.eig_im]);
        })?;
        if !one_mode_slice.is_empty() {
            ds_cf.append_staged(|b| {
                b.append_f64(one_mode_slice);
            })?;
        }

        phys_grp.set_attr("n_solutions", AttrValue::I64((cur_n + 1) as i64))?;
        phys_grp.set_attr("mode_labels_json", AttrValue::VarLenString(serde_json::to_string(&labels)?))?;
    }
    session.commit()?;
    drop(session);
    Ok(())
}

/// Writes a legacy `nyx-results` v3 file from [`NyxV3Summary`] (`Layout L1` per-solution `Deflate-3`)
/// atomically via `<target>.tmp` -> `std::fs::rename`.
pub fn write_nyx_v3(summary: &NyxV3Summary, target_path: &Path) -> Result<()> {
    let mut fb = FileBuilder::new();
    fb.set_attr("format", AttrValue::VarLenString(summary.format.clone()));
    fb.set_attr("version", AttrValue::I64(summary.version));
    fb.set_attr("writer", AttrValue::VarLenString(summary.writer.clone()));
    fb.set_attr("processes", AttrValue::I64(summary.processes));

    fb.create_dataset("project").with_u8_data(&summary.project_bytes);

    let mut mesh_grp = fb.create_group("mesh");
    mesh_grp.create_dataset("cell_labels").with_i32_data(&summary.mesh_cell_labels);
    mesh_grp
        .create_dataset("cells")
        .with_u32_data(&summary.mesh_cells)
        .with_shape(&[summary.mesh_cells_shape[0] as u64, summary.mesh_cells_shape[1] as u64])
        .with_chunks(&[summary.mesh_cells_shape[0].max(1) as u64, summary.mesh_cells_shape[1].max(1) as u64])
        .with_deflate(DeflateLevel::Three);
    mesh_grp.create_dataset("facet_labels").with_i32_data(&summary.mesh_facet_labels);
    mesh_grp
        .create_dataset("facets")
        .with_u32_data(&summary.mesh_facets)
        .with_shape(&[summary.mesh_facets_shape[0] as u64, summary.mesh_facets_shape[1] as u64]);
    mesh_grp
        .create_dataset("points")
        .with_f64_data(&summary.mesh_points)
        .with_shape(&[summary.mesh_points_shape[0] as u64, summary.mesh_points_shape[1] as u64]);
    fb.add_group(mesh_grp.finish());

    let mut res_grp = fb.create_group("results");
    let mut phys_grp = res_grp.create_group(&summary.physics);
    phys_grp.set_attr("physics", AttrValue::VarLenString(summary.physics.clone()));
    phys_grp.set_attr("study", AttrValue::VarLenString(summary.study.clone()));
    phys_grp.set_attr("dim", AttrValue::I64(summary.dim));
    phys_grp.set_attr("dofs", AttrValue::I64(summary.dofs));
    phys_grp.set_attr("n_solutions", AttrValue::I64(summary.n_solutions));
    phys_grp.set_attr("order", AttrValue::I64(summary.order));
    phys_grp.set_attr("seconds", AttrValue::F64(summary.seconds));
    phys_grp.create_dataset("log").with_u8_data(&summary.log_bytes);
    phys_grp.create_dataset("scalars").with_u8_data(&summary.scalars_bytes);

    for (k, sol) in summary.solutions.iter().enumerate() {
        let mut sol_grp = phys_grp.create_group(&format!("solution_{k}"));
        sol_grp.set_attr("eig_re", AttrValue::F64(sol.eig_re));
        sol_grp.set_attr("eig_im", AttrValue::F64(sol.eig_im));
        sol_grp.set_attr("label", AttrValue::VarLenString(sol.label.clone()));

        let mut cf_grp = sol_grp.create_group("cell_fields");
        for (cf_key, (c_name, c_unit, c_nloc, re, im)) in &sol.cell_fields {
            let mut comp_grp = cf_grp.create_group(cf_key);
            comp_grp.set_attr("name", AttrValue::VarLenString(c_name.clone()));
            comp_grp.set_attr("unit", AttrValue::VarLenString(c_unit.clone()));
            comp_grp.set_attr("nloc", AttrValue::I64(*c_nloc));
            comp_grp
                .create_dataset("re")
                .with_f64_data(re)
                .with_chunks(&[re.len().max(1) as u64])
                .with_deflate(DeflateLevel::Three);
            comp_grp
                .create_dataset("im")
                .with_f64_data(im)
                .with_chunks(&[im.len().max(1) as u64])
                .with_deflate(DeflateLevel::Three);
            cf_grp.add_group(comp_grp.finish());
        }
        sol_grp.add_group(cf_grp.finish());

        let mut f_grp = sol_grp.create_group("fields");
        f_grp.set_attr("order", AttrValue::VarLenString(sol.fields_order_json.clone()));
        for (f_key, (f_name, f_unit, f_comp, data)) in &sol.fields {
            f_grp
                .create_dataset(f_key)
                .with_f64_data(data)
                .set_attr("name", AttrValue::VarLenString(f_name.clone()))
                .set_attr("unit", AttrValue::VarLenString(f_unit.clone()))
                .set_attr("components", AttrValue::VarLenString(f_comp.clone()));
        }
        sol_grp.add_group(f_grp.finish());
        phys_grp.add_group(sol_grp.finish());
    }
    res_grp.add_group(phys_grp.finish());
    fb.add_group(res_grp.finish());

    let mut viz_grp = fb.create_group("viz");
    let mut vphys_grp = viz_grp.create_group(&summary.physics);
    vphys_grp.set_attr("topology", AttrValue::VarLenString(summary.viz_topology.clone()));
    vphys_grp
        .create_dataset("cells")
        .with_u32_data(&summary.viz_cells)
        .with_shape(&[summary.viz_cells_shape[0] as u64, summary.viz_cells_shape[1] as u64])
        .with_chunks(&[summary.viz_cells_shape[0].max(1) as u64, summary.viz_cells_shape[1].max(1) as u64])
        .with_deflate(DeflateLevel::Three);
    vphys_grp
        .create_dataset("points")
        .with_f64_data(&summary.viz_points)
        .with_shape(&[summary.viz_points_shape[0] as u64, summary.viz_points_shape[1] as u64])
        .with_chunks(&[summary.viz_points_shape[0].max(1) as u64, summary.viz_points_shape[1].max(1) as u64])
        .with_deflate(DeflateLevel::Three);

    for (k, sol) in summary.solutions.iter().enumerate() {
        let mut vsol_grp = vphys_grp.create_group(&format!("sol_{k}"));
        for (f_key, (_f_name, _f_unit, _f_comp, data)) in &sol.fields {
            let v_key = viz_link_name(f_key);
            vsol_grp.create_dataset(&v_key).with_f64_data(data);
        }
        vphys_grp.add_group(vsol_grp.finish());
    }
    viz_grp.add_group(vphys_grp.finish());
    fb.add_group(viz_grp.finish());

    if let Some(parent) = target_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let tmp_path = staging_path_for(target_path)?;
    fb.write(&tmp_path)?;
    atomic_rename(&tmp_path, target_path)?;
    Ok(())
}
