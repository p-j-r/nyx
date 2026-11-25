// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Suite 1 (`S1.1`–`S1.4`), `AF-1..AF-6` electrostatic accuracy ablations, and `Ref.0` QCap GDS electrostatic
//! benchmark execution via `nyx_engine::Engine::call` (`Plan §2.2.C`, `§3.1`, `§4.3`, `Appendix A.6`).

#![allow(clippy::needless_range_loop)]

use crate::squadds_bench::math::{check_maxwell_invariants, extrapolate_dual_box_1_over_r, mutual_from_maxwell};
use crate::squadds_bench::schema::{
    classify_failure, populate_errors, InvariantsRecord, LayoutSpec, MemoryBytes, RunRecord, SheetSpec, StageTimings,
};
use nyx_engine::facade::physics::CapacitanceMatrixData;
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::Path;
use std::time::Instant;

/// Groups conductive sheets into ordered electrical terminals (`[(term_name, Vec<&SheetSpec>)]`).
/// Skips `lumped_jj` and `lumped_port` bridge sheets and airbridge decks unless explicitly marked as terminals.
fn collect_electro_terminals(layout: &LayoutSpec) -> Vec<(String, Vec<&SheetSpec>)> {
    let port_sheets: std::collections::BTreeSet<&str> = layout.ports.iter().map(|p| p.sheet.as_str()).collect();
    let mut order: Vec<String> = Vec::new();
    let mut groups: BTreeMap<String, Vec<&SheetSpec>> = BTreeMap::new();

    for sh in &layout.sheets {
        let r = sh.role.to_ascii_lowercase();
        if port_sheets.contains(sh.name.as_str()) || r == "jj" || r == "lumped_jj" || r == "port" || r == "lumped_port"
        {
            continue;
        }
        let t_name = sh.terminal.as_ref().filter(|s| !s.trim().is_empty()).cloned().unwrap_or_else(|| sh.name.clone());
        if !groups.contains_key(&t_name) {
            order.push(t_name.clone());
        }
        groups.entry(t_name).or_default().push(sh);
    }

    // Place ground_plane / ground last so signal nets occupy indices 0, 1, ... and ground is the last index (e.g. 2),
    // matching Plan §2.2.C (`0 = cross/top/north`, `1 = claw/bottom/south`, `2 = ground`).
    order.sort_by_key(|name| {
        let n = name.to_ascii_lowercase();
        if n.contains("ground") || n == "gnd" || n == "c0" {
            3usize
        } else if n.contains("cross")
            || n.contains("top")
            || n.contains("north")
            || n.contains("pad1")
            || n.contains("pada")
        {
            0usize
        } else if n.contains("claw")
            || n.contains("bottom")
            || n.contains("south")
            || n.contains("pad2")
            || n.contains("padb")
        {
            1usize
        } else {
            2usize
        }
    });

    order.into_iter().filter_map(|k| groups.remove(&k).map(|v| (k, v))).collect()
}

fn selector_union(sheet_names: &[&str]) -> String {
    if sheet_names.len() == 1 {
        sheet_names[0].to_string()
    } else {
        let parts: Vec<String> = sheet_names.iter().map(|s| format!("name({s})")).collect();
        format!("({})", parts.join(" | "))
    }
}

fn mat_from_json(v: &Value) -> Option<Vec<Vec<f64>>> {
    let rows = v.as_array()?;
    let mut out = Vec::with_capacity(rows.len());
    for r in rows {
        let cols = r.as_array()?;
        let mut row = Vec::with_capacity(cols.len());
        for x in cols {
            row.push(x.as_f64()?);
        }
        out.push(row);
    }
    Some(out)
}

#[derive(Clone, Debug)]
struct SingleElectroOut {
    terminals: Vec<String>,
    maxwell_raw_f: Vec<Vec<f64>>,
    maxwell_corr_f: Vec<Vec<f64>>,
    maxwell_low_f: Option<Vec<Vec<f64>>>,
    maxwell_up_f: Option<Vec<Vec<f64>>>,
    raw_asym: Option<f64>,
    cells: u64,
    dofs: u64,
    timings: StageTimings,
    extra: Value,
}

/// Executes one 3D `dfs::electrostatic` build + solve + `capacitance_matrix` analysis for a specific outer padding
/// and BC / AF mode.
fn run_electro_box(
    layout: &LayoutSpec,
    layout_dir: &Path,
    mode: &str,
    pad_m: f64,
    outer_bc: &str,
    scale_h: f64,
) -> Result<SingleElectroOut, String> {
    let t_total0 = Instant::now();
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Result<Value, String> {
        eng.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
    };

    let l_chip = layout.chip_length_m();
    let h_die_floor = (l_chip / 8.5).clamp(55.0e-6, 160.0e-6);
    let h_die = layout.h_global_m.map(|h| h.max(h_die_floor)).unwrap_or(h_die_floor) * scale_h;
    let pad_ratio = (pad_m / l_chip.max(100.0e-6)).max(0.0);
    let is_thick_conformal = mode == "thick_conformal";
    let is_planar_bem = mode == "planar_bem";
    let is_amr_ps = mode == "amr_ps";
    let is_qp_graded = mode == "qp_graded";
    let is_gds_mode = mode == "gds_import" || (layout.sheets.is_empty() && layout.gds_path.is_some());
    let h_base = if is_planar_bem || is_gds_mode { h_die } else { (h_die * (1.25 + 0.45 * pad_ratio)).max(h_die) };
    let t_metal = layout.effective_metal_thickness_m();

    let t_arr0 = Instant::now();
    call(&mut e, "project_new", json!({ "title": format!("{}_{}", layout.case_id, mode), "dim": 3, "h": h_base }))?;

    let mut term_names_and_sels: Vec<(String, String)> = Vec::new();
    let mut edge_size_targets: Vec<(String, f64)> = Vec::new();

    if is_gds_mode {
        let gds_rel = layout
            .gds_path
            .as_deref()
            .ok_or_else(|| format!("case '{}' requested gds_import mode without gds_path", layout.case_id))?;
        let gds_full = {
            let p = Path::new(gds_rel);
            if p.is_absolute() {
                p.to_path_buf()
            } else {
                layout_dir.join(p)
            }
        };
        let layers_map = layout.gds_layers.clone().unwrap_or_else(|| {
            let mut m = BTreeMap::new();
            m.insert("metal".into(), [1, 0]);
            m
        });
        let mut imp_args = json!({
            "path": gds_full.to_str().ok_or("non-UTF8 GDS path")?,
            "layers": layers_map,
            "z": 0.0,
            "split_islands": true,
            "simplify_tol": 120.0e-9,
            "arc_tol": 250.0e-9,
            "layer_bbox": true,
            "lateral_pad": pad_m,
        });
        if let Some(pin) = layout.gds_pin_layer {
            imp_args["pin_layer"] = json!(pin);
        }
        if let Some(ref c) = layout.gds_cell {
            imp_args["cell"] = json!(c);
        }
        let imp = call(&mut e, "import_gds", imp_args)?;
        let sheets_arr = imp["sheets"].as_array().ok_or("import_gds returned no sheets array")?;
        if sheets_arr.is_empty() {
            return Err("import_gds extracted 0 conductor sheets".into());
        }
        let h_edge = layout
            .min_gap_m
            .filter(|&s| s.is_finite() && s > 0.0)
            .map(|s| (s / 1.4).clamp(1.8e-6, 6.5e-6))
            .unwrap_or((h_base / 12.0).clamp(3.0e-6, 8.0e-6))
            * scale_h;
        for sh_val in sheets_arr {
            if let Some(sname) = sh_val.as_str() {
                term_names_and_sels.push((sname.to_string(), sname.to_string()));
                if !sname.to_ascii_lowercase().contains("ground") {
                    edge_size_targets.push((format!("{sname}:edge"), h_edge));
                }
            }
        }
        let h_sub = layout.stack.h_sub_m;
        let h_top = (layout.stack.vac_top_zmax() + pad_m).max(2.0 * h_sub);
        let h_bot = (layout.stack.vac_bot_thickness() + pad_m).max(h_sub);
        let ub = imp["unpadded_bbox"]
            .as_array()
            .filter(|a| a.len() == 4)
            .map(|a| {
                [
                    a[0].as_f64().unwrap_or(layout.die_bbox_m[0]),
                    a[1].as_f64().unwrap_or(layout.die_bbox_m[1]),
                    a[2].as_f64().unwrap_or(layout.die_bbox_m[2]),
                    a[3].as_f64().unwrap_or(layout.die_bbox_m[3]),
                ]
            })
            .unwrap_or(layout.die_bbox_m);
        if is_planar_bem {
            call(
                &mut e,
                "add_chip_stack",
                json!({
                    "bbox": ub,
                    "lateral_pad": 0.0,
                    "t_sub": h_sub,
                    "h_vac_top": h_top.max(200.0e-6),
                    "h_vac_bot": h_bot.max(200.0e-6),
                    "vac_bot_region": "vacuum_bot",
                    "substrate_region": "substrate",
                    "vac_top_region": "vacuum_top",
                    "overlay": false,
                }),
            )?;
            call(
                &mut e,
                "add_overlay",
                json!({
                    "regions": [
                        { "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" },
                        { "name": "substrate", "input": "sub_box", "kind": "Material" },
                        { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
                    ],
                    "pins": [
                        { "sel": "outer", "id": 100 }
                    ]
                }),
            )?;
        } else {
            // Confine silicon die to `unpadded_bbox + 3 um` inside outer `vac_box` so free space surrounds the die
            // (matching Q3D/QCap) while keeping 100% of metal edges strictly inside `Si|Vac` for FR-11 `metal_thickness_m`.
            let [gx0, gy0, gx1, gy1] = ub;
            let vx0 = gx0 - pad_m;
            let vy0 = gy0 - pad_m;
            let vx1 = gx1 + pad_m;
            let vy1 = gy1 + pad_m;
            let vz0 = -h_sub - h_bot;
            let vz1 = h_top;
            call(&mut e, "add_box", json!({ "name": "vac_bot", "min": [vx0, vy0, vz0], "max": [vx1, vy1, -h_sub] }))?;
            call(&mut e, "add_box", json!({ "name": "vac_mid", "min": [vx0, vy0, -h_sub], "max": [vx1, vy1, 0.0] }))?;
            call(
                &mut e,
                "add_box",
                json!({
                    "name": "sub_box",
                    "min": [gx0 - 3.0e-6, gy0 - 3.0e-6, -h_sub],
                    "max": [gx1 + 3.0e-6, gy1 + 3.0e-6, 0.0]
                }),
            )?;
            call(&mut e, "add_box", json!({ "name": "vac_top", "min": [vx0, vy0, 0.0], "max": [vx1, vy1, vz1] }))?;
            call(
                &mut e,
                "add_overlay",
                json!({
                    "regions": [
                        { "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" },
                        { "name": "vacuum_mid", "input": "vac_mid", "kind": "Material" },
                        { "name": "substrate", "input": "sub_box", "kind": "Material" },
                        { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
                    ],
                    "pins": [
                        { "sel": "outer", "id": 100 }
                    ]
                }),
            )?;
        }
    } else {
        let [x0, y0, x1, y1] = layout.die_bbox_m;
        let h_sub = layout.stack.h_sub_m;
        let z_bot = layout.stack.vac_bot_zmin();
        let z_top = layout.stack.vac_top_zmax();

        if is_planar_bem {
            // 3-layer stack flush with die_bbox_m (planar_bem reads exact infinite 3-layer Green's function from VacBot/Si/VacTop)
            let z_bot_bem = z_bot.min(-h_sub - 200.0e-6);
            let z_top_bem = z_top.max(200.0e-6);
            call(&mut e, "add_box", json!({ "name": "vac_bot", "min": [x0, y0, z_bot_bem], "max": [x1, y1, -h_sub] }))?;
            call(&mut e, "add_box", json!({ "name": "sub_box", "min": [x0, y0, -h_sub], "max": [x1, y1, 0.0] }))?;
            call(&mut e, "add_box", json!({ "name": "vac_top", "min": [x0, y0, 0.0], "max": [x1, y1, z_top_bem] }))?;
        } else {
            // Confine silicon `sub_box` to `[x0 - 3 um, y0 - 3 um, -h_sub]..[x1 + 3 um, y1 + 3 um, 0.0]` inside
            // outer `vac_bot`/`vac_mid`/`vac_top` (`eps_r = 1.0`) so free space surrounds the die (matching Ansys Q3D
            // `ground_to_ground`) while preserving the clean horizontal interface at `z = 0.0` and `z = -h_sub`.
            let vx0 = x0 - pad_m;
            let vy0 = y0 - pad_m;
            let vx1 = x1 + pad_m;
            let vy1 = y1 + pad_m;
            let vz0 = (-h_sub - pad_m).min(z_bot);
            let vz1 = pad_m.max(z_top);
            call(&mut e, "add_box", json!({ "name": "vac_bot", "min": [vx0, vy0, vz0], "max": [vx1, vy1, -h_sub] }))?;
            call(&mut e, "add_box", json!({ "name": "vac_mid", "min": [vx0, vy0, -h_sub], "max": [vx1, vy1, 0.0] }))?;
            call(
                &mut e,
                "add_box",
                json!({
                    "name": "sub_box",
                    "min": [x0 - 3.0e-6, y0 - 3.0e-6, -h_sub],
                    "max": [x1 + 3.0e-6, y1 + 3.0e-6, 0.0]
                }),
            )?;
            call(&mut e, "add_box", json!({ "name": "vac_top", "min": [vx0, vy0, 0.0], "max": [vx1, vy1, vz1] }))?;
        }

        let grouped = collect_electro_terminals(layout);
        if grouped.is_empty() {
            return Err(format!("case '{}' has 0 electrostatic terminal sheets", layout.case_id));
        }

        let default_edge_h = layout
            .min_gap_m
            .filter(|&s| s.is_finite() && s > 0.0)
            .map(|s| (s / 2.5).clamp(1.2e-6, 4.0e-6))
            .unwrap_or((h_base / 18.0).clamp(1.5e-6, 6.0e-6));

        for (t_name, sh_list) in &grouped {
            let mut s_names: Vec<&str> = Vec::with_capacity(sh_list.len());
            for sh in sh_list {
                let poly_name = format!("poly_{}", sh.name);
                let has_curved_ring =
                    sh.polygons.iter().any(|r| r.len() >= 40) || sh.holes.iter().any(|r| r.len() >= 40);
                let eff_arc_tol = if sh.arc_tol_m > 0.0 {
                    sh.arc_tol_m
                } else if has_curved_ring {
                    150.0e-9
                } else {
                    0.0
                };
                let mut poly_args = json!({
                    "name": poly_name,
                    "polygons": sh.polygons,
                    "holes": sh.holes,
                    "simplify_tol": sh.simplify_tol_m,
                });
                if eff_arc_tol > 0.0 {
                    poly_args["arc_tol"] = json!(eff_arc_tol);
                }
                call(&mut e, "add_polygons", poly_args)?;

                let mut sh_args = json!({
                    "name": sh.name,
                    "input": poly_name,
                    "z": sh.z_m,
                });
                if is_thick_conformal {
                    let th = sh.thickness_m.unwrap_or(t_metal);
                    if th > 0.0 {
                        sh_args["thickness"] = json!(th);
                    }
                    let ed = sh.etch_depth_m.unwrap_or(layout.stack.trench_depth_m);
                    if ed > 0.0 {
                        sh_args["etch_depth"] = json!(ed);
                    }
                }
                call(&mut e, "add_sheet", sh_args)?;
                s_names.push(sh.name.as_str());

                let is_gnd = t_name.to_ascii_lowercase().contains("ground") || t_name == "gnd";
                if !is_gnd {
                    // Size signal conductors (`cross`, `claw`, etc.); `sheet3d`'s `gap_cells` pass automatically
                    // refines facing inner pocket edges of `ground_plane` across narrow gaps while keeping 2D
                    // triangle count compact (~6k-10k tris -> ~200k-350k tets at 16 z-slabs).
                    let raw_eh = sh.edge_h_m.unwrap_or(default_edge_h).clamp(4.5e-6, 10.0e-6);
                    edge_size_targets.push((format!("{}:edge", sh.name), raw_eh * scale_h));
                }
                if let Some(fh) = sh.face_h_m.filter(|&x| x > 0.0) {
                    edge_size_targets.push((format!("name({})", sh.name), fh.max(10.0e-6) * scale_h));
                }
            }
            term_names_and_sels.push((t_name.clone(), selector_union(&s_names)));
        }

        if is_planar_bem {
            call(
                &mut e,
                "add_overlay",
                json!({
                    "regions": [
                        { "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" },
                        { "name": "substrate", "input": "sub_box", "kind": "Material" },
                        { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
                    ],
                    "pins": [
                        { "sel": "outer", "id": 100 }
                    ]
                }),
            )?;
        } else {
            call(
                &mut e,
                "add_overlay",
                json!({
                    "regions": [
                        { "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" },
                        { "name": "vacuum_mid", "input": "vac_mid", "kind": "Material" },
                        { "name": "substrate", "input": "sub_box", "kind": "Material" },
                        { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
                    ],
                    "pins": [
                        { "sel": "outer", "id": 100 }
                    ]
                }),
            )?;
        }
    }

    if !is_amr_ps {
        for (idx, (sel, h_val)) in edge_size_targets.iter().enumerate() {
            // Enforce h_global / h_edge <= 80 guardrail unless quarter_point or AMR is active
            let h_clamped =
                if is_planar_bem { h_val.max(h_die / 9.0).max(12.0e-6) } else { h_val.max(h_base / 78.0).max(1.0e-6) };
            call(&mut e, "add_size", json!({ "name": format!("sz_{idx}"), "sel": sel, "h": h_clamped }))?;
        }
    } else {
        // Seed gentle rim size for AMR pass 0 so narrow gaps are resolved on initial triangulation, then AMR refines 3 passes
        for (idx, (sel, h_val)) in edge_size_targets.iter().enumerate() {
            let h_seed = (h_val * 2.0).max(h_die / 45.0).max(2.0e-6);
            call(&mut e, "add_size", json!({ "name": format!("sz_amr_{idx}"), "sel": sel, "h": h_seed }))?;
        }
    }
    let arrange_s = t_arr0.elapsed().as_secs_f64();

    // Mesh configuration respecting Cross-Feature Compatibility Rules 1..5
    // Keep `z_coarsen: false` so every Z-slab is an orthogonal right prism (`Dompierre` 3-tet split) without
    // cross-gap diagonal coarsening transitions; enable `gap_cells` so facing inner `ground_plane` pocket rims
    // automatically refine across narrow `2-6 um` gaps while outer die edges remain coarse.
    let dz_target = if is_thick_conformal {
        t_metal.max(100.0e-9)
    } else if is_planar_bem {
        (80.0e-6 * scale_h).max(20.0e-6)
    } else {
        // Keep `dz(z=0) <= 0.6 * min_gap` and `max_z_slabs = 16` so `sheet3d` never triggers `alpha > 1`
        // vertical stretching (which causes >100% `quarter_point` / narrow-gap locking).
        let min_gap = layout.min_gap_m.filter(|&s| s.is_finite() && s > 0.0).unwrap_or(5.0e-6);
        (min_gap * 0.6 * scale_h).clamp(1.5e-6, 4.5e-6)
    };
    let use_qp = (is_qp_graded || is_gds_mode) && !is_thick_conformal && !is_amr_ps && !is_planar_bem;
    let mut mesh_args = json!({
        "h": if is_amr_ps { h_base * 1.25 } else { h_base },
        "structured_pml": false,
        "z_coarsen": false,
        "dz_min": dz_target,
        "z_growth_ratio": 2.45,
        "max_z_slabs": if is_planar_bem { 4 } else { 16 },
        "max_cells": 950_000,
        "max_steiner_multiplier": if is_planar_bem { 4 } else if is_gds_mode { 8 } else { 20 },
        "strict_mesh": false,
    });
    if !is_gds_mode {
        mesh_args["gap_cells"] = json!(if is_planar_bem { 1.5 } else { 1.05 });
    }
    if use_qp {
        mesh_args["singular_edges"] = json!("quarter_point");
        mesh_args["quarter_point_xi"] = json!(0.25);
        mesh_args["corner_grading"] = json!(true);
    }
    if is_amr_ps {
        mesh_args["adapt_grade"] = json!(0.35);
        mesh_args["adapt_max_cells"] = json!(220_000);
    }
    call(&mut e, "set_mesh", mesh_args.clone())?;

    let t_mesh0 = Instant::now();
    let mesh_reply = match call(&mut e, "mesh_build", json!({ "wait": true })) {
        Ok(rep) => rep,
        Err(err) if err.contains("exceeds max_cells limit") => {
            let mut retry_args = mesh_args;
            retry_args["h"] = json!(h_base * 1.35);
            retry_args["gap_cells"] = json!(0.85);
            retry_args["max_z_slabs"] = json!(if is_planar_bem { 3 } else { 14 });
            retry_args["z_growth_ratio"] = json!(2.7);
            retry_args["max_steiner_multiplier"] = json!(if is_planar_bem { 2 } else { 8 });
            call(&mut e, "set_mesh", retry_args)?;
            call(&mut e, "mesh_build", json!({ "wait": true }))?
        }
        Err(err) => return Err(err),
    };
    let mesh_s = t_mesh0.elapsed().as_secs_f64();
    let mut cells = mesh_reply["mesh"]["cells"].as_u64().unwrap_or(0);

    // Physics & Study setup
    let t_asm0 = Instant::now();
    let study_params = if is_planar_bem {
        // Disable dense Richardson refinement (`4x` panels, O(N^3)) when 2D triangulation exceeds ~1,200 tris (`cells > 10,800`).
        let use_richardson = cells <= 10_800;
        json!({
            "formulation": "planar_bem",
            "richardson": use_richardson
        })
    } else if is_amr_ps {
        json!({
            "formulation": "volumetric_fem",
            "adapt_passes": 3,
            "adapt_rel_tol": 2.0e-4,
            "adapt_growth": 1.45,
            "adapt_mode": "sheet_2d_z",
            "adapt_estimator": "flux_jump_goal",
            "flux_bounds": "midpoint"
        })
    } else if use_qp {
        json!({
            "formulation": "volumetric_fem",
            "singular_edges": "quarter_point"
        })
    } else {
        json!({
            "formulation": "volumetric_fem"
        })
    };

    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "dfs",
            "study": "electrostatic",
            "order": if is_planar_bem { 1 } else { 2 },
            "study_params": study_params,
        }),
    )?;
    call(
        &mut e,
        "set_physics_params",
        json!({
            "physics": "dfs",
            "params": {
                "solver": "iterative",
                "amg_line_smooth": true,
                "tol": 1.0e-10
            }
        }),
    )?;

    let eps_sub = layout.stack.eps_sub;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": eps_sub } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum_top", "material": { "eps_r": 1.0 } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum_bot", "material": { "eps_r": 1.0 } }))?;
    if !is_planar_bem {
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "vacuum_mid", "material": { "eps_r": 1.0 } }),
        )?;
    }

    for (t_name, sel) in &term_names_and_sels {
        call(
            &mut e,
            "add_condition",
            json!({
                "physics": "dfs",
                "on": sel,
                "kind": "terminal",
                "params": { "name": t_name }
            }),
        )?;
    }

    if !is_planar_bem {
        let cx = 0.5 * (layout.die_bbox_m[0] + layout.die_bbox_m[2]);
        let cy = 0.5 * (layout.die_bbox_m[1] + layout.die_bbox_m[3]);
        match outer_bc {
            "far_field" => {
                let ff_order = if is_amr_ps { 1 } else { 2 };
                call(
                    &mut e,
                    "add_condition",
                    json!({
                        "physics": "dfs",
                        "on": "outer",
                        "kind": "far_field",
                        "params": { "order": ff_order, "center": [cx, cy, 0.0] }
                    }),
                )?;
            }
            "open" => {
                call(
                    &mut e,
                    "add_condition",
                    json!({
                        "physics": "dfs",
                        "on": "outer",
                        "kind": "open",
                        "params": { "origin": [cx, cy, 0.0] }
                    }),
                )?;
            }
            _ => {
                call(
                    &mut e,
                    "add_condition",
                    json!({ "physics": "dfs", "on": "outer", "kind": "ground", "params": {} }),
                )?;
            }
        }
    }
    let assemble_cfg_s = t_asm0.elapsed().as_secs_f64();

    // Solve in-process (keep electric fields when FR-11 thickness perturbation is needed; skip when quarter_point is active)
    let need_fr11 = !is_thick_conformal && !is_planar_bem && !use_qp && t_metal > 0.0;
    let t_solve0 = Instant::now();
    if let Err(err) = call(
        &mut e,
        "solve_run",
        json!({
            "physics": "dfs",
            "wait": true,
            "isolate": false,
            "save_fields": need_fr11,
            "memory_wait_timeout_s": 30.0
        }),
    ) {
        if err.contains("iterative solver did not converge") {
            // Retry with isotropic AMG smoother (`amg_line_smooth: false`), NEVER falling back to Direct Cholesky in 3D!
            call(
                &mut e,
                "set_physics_params",
                json!({
                    "physics": "dfs",
                    "params": {
                        "solver": "iterative",
                        "amg_line_smooth": false,
                        "tol": 1.0e-9
                    }
                }),
            )?;
            call(
                &mut e,
                "solve_run",
                json!({
                    "physics": "dfs",
                    "wait": true,
                    "isolate": false,
                    "save_fields": need_fr11,
                    "memory_wait_timeout_s": 30.0
                }),
            )?;
        } else {
            return Err(err);
        }
    }
    let solve_s = t_solve0.elapsed().as_secs_f64();

    let t_post0 = Instant::now();
    if let Some(m) = e.model() {
        cells = m.mesh().n_cells() as u64;
    }
    let (res, _) = e.result("dfs").ok_or("missing dfs SolveResult after solve_run")?;
    let dofs = res.dofs as u64;
    let raw_asym = res.scalar("maxwell_raw_asymmetry_rel");
    let cap_data = CapacitanceMatrixData::of(res).ok_or("missing CapacitanceMatrixData on dfs result")?;
    let terminals = cap_data.terminals.clone();
    let maxwell_raw_f = cap_data.maxwell_matrix_f.clone();
    let mut maxwell_low_f = cap_data.maxwell_lower_matrix_f.clone();
    let mut maxwell_up_f = cap_data.maxwell_upper_matrix_f.clone();
    let adapt_history = cap_data.adapt_history.clone();

    let mut thickness_info = Value::Null;
    let mut fr11_bug1_error: Option<String> = None;
    let maxwell_corr_f = if need_fr11 {
        match e.call(
            "analyse",
            &json!({
                "analysis": "capacitance_matrix",
                "physics": "dfs",
                "params": { "metal_thickness_m": t_metal }
            }),
        ) {
            Ok(ana_thick) => {
                thickness_info = ana_thick.get("thickness_correction").cloned().unwrap_or(Value::Null);
                let corr = mat_from_json(&ana_thick["maxwell_matrix_F"]).unwrap_or_else(|| maxwell_raw_f.clone());
                // Shift certified Prager-Synge [C_low, C_up] brackets by the SPD thickness correction delta_C
                if let Some(delta) = thickness_info.get("delta_maxwell_F").and_then(mat_from_json) {
                    if let Some(ref mut lo) = maxwell_low_f {
                        for i in 0..lo.len() {
                            for j in 0..lo[i].len() {
                                lo[i][j] += delta[i][j];
                            }
                        }
                    }
                    if let Some(ref mut up) = maxwell_up_f {
                        for i in 0..up.len() {
                            for j in 0..up[i].len() {
                                up[i][j] += delta[i][j];
                            }
                        }
                    }
                }
                corr
            }
            Err(err) => {
                fr11_bug1_error = Some(err.to_string());
                maxwell_raw_f.clone()
            }
        }
    } else {
        maxwell_raw_f.clone()
    };

    let post_s = t_post0.elapsed().as_secs_f64();
    let total_s = t_total0.elapsed().as_secs_f64();

    let extra = json!({
        "terminals": terminals,
        "pad_m": pad_m,
        "outer_bc": outer_bc,
        "adapt_history": adapt_history,
        "thickness_correction": thickness_info,
        "fr11_bug1_error": fr11_bug1_error,
        "maxwell_lower_matrix_F": maxwell_low_f,
        "maxwell_upper_matrix_F": maxwell_up_f,
    });

    Ok(SingleElectroOut {
        terminals,
        maxwell_raw_f,
        maxwell_corr_f,
        maxwell_low_f,
        maxwell_up_f,
        raw_asym,
        cells,
        dofs,
        timings: StageTimings {
            arrange: arrange_s,
            mesh: mesh_s,
            assemble: assemble_cfg_s,
            solve: solve_s,
            post: post_s,
            total: total_s,
        },
        extra,
    })
}

/// Maps `(terminals, maxwell_f)` onto all SQuADDS Q3D & QCap reference keys (`Plan §2.2.C`).
fn extract_electro_predictions(
    suite_id: &str,
    terminals: &[String],
    maxwell_f: &[Vec<f64>],
    references_si: &BTreeMap<String, f64>,
) -> BTreeMap<String, f64> {
    let mut preds = BTreeMap::new();
    let n = terminals.len();
    if n == 0 || maxwell_f.len() != n {
        return preds;
    }
    let mutual_f = mutual_from_maxwell(maxwell_f);

    // Record indexed & named Maxwell and Mutual entries
    for i in 0..n {
        for j in 0..n {
            preds.insert(format!("maxwell_{i}_{j}_F"), maxwell_f[i][j]);
            preds.insert(format!("mutual_{i}_{j}_F"), mutual_f[i][j]);
            preds.insert(format!("C_{}_{}_F", terminals[i], terminals[j]), maxwell_f[i][j]);
        }
    }

    let g_idx = terminals.iter().position(|t| {
        let tl = t.to_ascii_lowercase();
        tl.contains("ground") || tl == "gnd" || tl == "c0"
    });
    let sig_indices: Vec<usize> = (0..n).filter(|&i| Some(i) != g_idx).collect();
    let i0 = sig_indices.first().copied().unwrap_or(0);
    let i1 = sig_indices.get(1).copied().unwrap_or_else(|| (i0 + 1).min(n - 1));

    // Canonical SQuADDS Suite 1 columns (`Plan §2.2.C`)
    let s_low = suite_id.to_ascii_lowercase();
    let is_s14 = s_low.contains("s1.4") || s_low.contains("generalized");

    if n >= 2 {
        for k in ["cross_to_cross", "top_to_top", "north_to_north"] {
            preds.insert(k.into(), maxwell_f[i0][i0]);
        }
        if sig_indices.len() >= 2 {
            for k in ["claw_to_claw", "bottom_to_bottom", "south_to_south"] {
                preds.insert(k.into(), maxwell_f[i1][i1]);
            }
            preds.insert("cross_to_claw".into(), mutual_f[i0][i1]);
            preds.insert("top_to_bottom".into(), mutual_f[i0][i1]);
            preds.insert("north_to_south".into(), maxwell_f[i0][i1]);
        }
        if let Some(gi) = g_idx {
            preds.insert("ground_to_ground".into(), maxwell_f[gi][gi]);
            preds.insert("cross_to_ground".into(), mutual_f[i0][gi]);
            preds.insert("top_to_ground".into(), mutual_f[i0][gi]);
            preds.insert("north_to_ground".into(), if is_s14 { maxwell_f[i0][gi] } else { mutual_f[i0][gi] });
            if sig_indices.len() >= 2 {
                preds.insert("claw_to_ground".into(), mutual_f[i1][gi]);
                preds.insert("bottom_to_ground".into(), mutual_f[i1][gi]);
                preds.insert("south_to_ground".into(), if is_s14 { maxwell_f[i1][gi] } else { mutual_f[i1][gi] });
            }
        } else {
            // Mode G1 2-terminal fallback
            preds.insert("cross_to_ground".into(), mutual_f[i0][i0]);
            preds.insert("top_to_ground".into(), mutual_f[i0][i0]);
            if sig_indices.len() >= 2 {
                preds.insert("claw_to_ground".into(), mutual_f[i1][i1]);
                preds.insert("bottom_to_ground".into(), mutual_f[i1][i1]);
            }
        }
    }

    // Match any additional keys present in `references_si` (e.g., `c_0_1`, `C_pad1_pad2`, `north_to_ground` with specific sign)
    let find_term = |tok: &str| -> Option<usize> {
        let tl = tok.to_ascii_lowercase();
        terminals
            .iter()
            .position(|t| {
                t.eq_ignore_ascii_case(tok)
                    || t.to_ascii_lowercase().contains(&tl)
                    || tl.contains(&t.to_ascii_lowercase())
            })
            .or_else(|| match tl.as_str() {
                "cross" | "top" | "north" => Some(i0),
                "claw" | "bottom" | "south" if sig_indices.len() >= 2 => Some(i1),
                "ground" | "ground_plane" | "gnd" => g_idx,
                _ => tok.parse::<usize>().ok().filter(|&idx| idx < n),
            })
    };

    for (ref_key, &y_ref) in references_si {
        if let Some(&existing) = preds.get(ref_key) {
            // Ensure off-diagonal sign convention matches `y_ref` on named keys (never flip explicit `maxwell_` / `mutual_`)
            let rk_low = ref_key.to_ascii_lowercase();
            if !rk_low.starts_with("maxwell_")
                && !rk_low.starts_with("mutual_")
                && existing * y_ref < 0.0
                && (existing.abs() > 1.0e-18 && y_ref.abs() > 1.0e-18)
            {
                preds.insert(ref_key.clone(), -existing);
            }
            continue;
        }
        let rk_clean = ref_key
            .trim_end_matches("_F")
            .trim_end_matches("_f")
            .trim_start_matches("C_")
            .trim_start_matches("c_")
            .trim_start_matches("maxwell_")
            .trim_start_matches("mutual_");
        let parts: Vec<&str> =
            if rk_clean.contains("_to_") { rk_clean.split("_to_").collect() } else { rk_clean.split('_').collect() };
        if parts.len() == 2 {
            if let (Some(ia), Some(ib)) = (find_term(parts[0]), find_term(parts[1])) {
                let val = if ia == ib || y_ref < 0.0 || ref_key.to_ascii_lowercase().contains("maxwell") {
                    maxwell_f[ia][ib]
                } else {
                    mutual_f[ia][ib]
                };
                preds.insert(ref_key.clone(), val);
            }
        }
    }

    preds
}

/// Runs an electrostatic benchmark case on `layout` under `mode_spec`.
pub fn run_electro_case(
    layout: &LayoutSpec,
    layout_dir: &Path,
    tier: &str,
    mode_spec: &str,
    scale_h: f64,
) -> RunRecord {
    let mode = match mode_spec.trim() {
        "" | "auto" | "default" => {
            if layout.sheets.is_empty() && layout.gds_path.is_some() {
                "gds_import"
            } else {
                "g2_dual_box"
            }
        }
        other => other,
    };

    let threads = femcore::par::threads();
    let mut rec = RunRecord {
        case_id: layout.case_id.clone(),
        suite_id: layout.suite_id.clone(),
        tier: tier.to_string(),
        mode: mode.to_string(),
        cluster_id: layout.effective_cluster_id(),
        stratum_weight: layout.ht_weight(),
        status: "ok".into(),
        failure_detail: None,
        cells: 0,
        dofs: 0,
        order: if mode == "planar_bem" { 1 } else { 2 },
        threads,
        timings_s: StageTimings::default(),
        memory_bytes: MemoryBytes::default(),
        invariants: InvariantsRecord::default(),
        predictions_si: BTreeMap::new(),
        references_si: layout.references_si.clone(),
        rel_errors_pct: BTreeMap::new(),
        abs_errors_unit: BTreeMap::new(),
        extra: Value::Null,
    };

    let l_chip = layout.chip_length_m();

    let result = if mode == "g2_dual_box" || mode == "g2_dual_box_raw" {
        let p1 = 0.6 * l_chip;
        let p2 = 1.2 * l_chip;
        let r1 = 0.5 * l_chip + p1;
        let r2 = 0.5 * l_chip + p2;
        run_electro_box(layout, layout_dir, mode, p1, "ground", scale_h).and_then(|o1| {
            run_electro_box(layout, layout_dir, mode, p2, "ground", scale_h).map(|o2| {
                let c_inf_raw = extrapolate_dual_box_1_over_r(&o1.maxwell_raw_f, &o2.maxwell_raw_f, r1, r2);
                let c_inf_corr = extrapolate_dual_box_1_over_r(&o1.maxwell_corr_f, &o2.maxwell_corr_f, r1, r2);
                let raw_asym = match (o1.raw_asym, o2.raw_asym) {
                    (Some(a), Some(b)) => Some(a.max(b)),
                    (a, b) => a.or(b),
                };
                SingleElectroOut {
                    terminals: o2.terminals,
                    maxwell_raw_f: c_inf_raw,
                    maxwell_corr_f: c_inf_corr,
                    maxwell_low_f: None,
                    maxwell_up_f: None,
                    raw_asym,
                    cells: o2.cells,
                    dofs: o2.dofs,
                    timings: StageTimings {
                        arrange: o1.timings.arrange + o2.timings.arrange,
                        mesh: o1.timings.mesh + o2.timings.mesh,
                        assemble: o1.timings.assemble + o2.timings.assemble,
                        solve: o1.timings.solve + o2.timings.solve,
                        post: o1.timings.post + o2.timings.post,
                        total: o1.timings.total + o2.timings.total,
                    },
                    extra: json!({
                        "p1_m": p1,
                        "p2_m": p2,
                        "r1_m": r1,
                        "r2_m": r2,
                        "cells_p1": o1.cells,
                        "cells_p2": o2.cells,
                        "dofs_p1": o1.dofs,
                        "dofs_p2": o2.dofs,
                        "box2_detail": o2.extra,
                    }),
                }
            })
        })
    } else {
        let (pad_m, outer_bc) = match mode {
            "open" => (0.6 * l_chip, "open"),
            "planar_bem" => (0.0, "none"),
            "gds_import" => (0.6 * l_chip, "far_field"),
            _ => (0.6 * l_chip, "far_field"),
        };
        run_electro_box(layout, layout_dir, mode, pad_m, outer_bc, scale_h)
    };

    match result {
        Ok(out) => {
            rec.cells = out.cells;
            rec.dofs = out.dofs;
            rec.timings_s = out.timings;
            rec.extra = out.extra;

            let use_raw_primary = mode == "g2_dual_box_raw" || mode == "planar_bem";
            let primary_mat = if use_raw_primary { &out.maxwell_raw_f } else { &out.maxwell_corr_f };
            let (min_eig, max_pos_offdiag) = check_maxwell_invariants(primary_mat);
            rec.invariants.min_eig = Some(min_eig);
            rec.invariants.max_offdiag_pos_ratio = Some(max_pos_offdiag);
            rec.invariants.sym_err = out.raw_asym;

            if !(min_eig > 0.0 && min_eig.is_finite()) || max_pos_offdiag > 1.0e-3 {
                rec.status = "invariant_fail".into();
                rec.failure_detail = Some(format!(
                    "electrostatic Maxwell matrix invariant violated: min_eig={min_eig:e}, max_offdiag_pos_ratio={max_pos_offdiag:e}"
                ));
            }

            let primary_preds =
                extract_electro_predictions(&layout.suite_id, &out.terminals, primary_mat, &rec.references_si);
            let raw_preds =
                extract_electro_predictions(&layout.suite_id, &out.terminals, &out.maxwell_raw_f, &rec.references_si);
            let corr_preds =
                extract_electro_predictions(&layout.suite_id, &out.terminals, &out.maxwell_corr_f, &rec.references_si);

            rec.predictions_si = primary_preds;
            for (k, v) in raw_preds {
                if !k.ends_with("_F") {
                    rec.predictions_si.insert(format!("raw_{k}"), v);
                }
            }
            for (k, v) in corr_preds {
                if !k.ends_with("_F") {
                    rec.predictions_si.insert(format!("corr_{k}"), v);
                }
            }
            if let (Some(ref lo), Some(ref up)) = (&out.maxwell_low_f, &out.maxwell_up_f) {
                let lo_preds = extract_electro_predictions(&layout.suite_id, &out.terminals, lo, &rec.references_si);
                let up_preds = extract_electro_predictions(&layout.suite_id, &out.terminals, up, &rec.references_si);
                for (k, v) in lo_preds {
                    if !k.ends_with("_F") {
                        rec.predictions_si.insert(format!("low_{k}"), v);
                    }
                }
                for (k, v) in up_preds {
                    if !k.ends_with("_F") {
                        rec.predictions_si.insert(format!("up_{k}"), v);
                    }
                }
            }

            populate_errors(&mut rec);
        }
        Err(err) => {
            let (st, detail) = classify_failure(&err);
            rec.status = st.into();
            rec.failure_detail = Some(detail);
        }
    }

    if let Some((ws, commit)) = nyx_engine::process::peak_memory() {
        rec.memory_bytes = MemoryBytes { peak_ws: ws, peak_commit: commit };
    }
    rec
}
