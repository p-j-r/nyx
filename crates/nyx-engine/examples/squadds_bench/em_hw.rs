// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Electromagnetic cavity (`ems::cavity`), driven-modal PROM (`ems::driven`), Unified `QubitCavity` 4-way
//! equivalence (`M1`–`M4`), 3D CPW airbridges (`S3.4`), and cryogenic hardware validation (`WM1`, `MUNINN`,
//! `Dissipator`) execution via `nyx_engine::Engine::call` (`Plan §3.2`–`§3.4`, `Appendix A.1`–`A.5`).

#![allow(clippy::needless_range_loop)]

use crate::squadds_bench::electro::run_electro_case;
use crate::squadds_bench::math::{
    check_s_matrix_invariants, fit_damped_lorentzian_pole, fit_driven_capacitance_omega_zero, fit_resonator_s21,
    squadds_compositional_qubit_cavity, squadds_find_kappa_half_wave, synthesize_kubo_feedline_2port, Cpx,
    TransmonKuboModel,
};
use crate::squadds_bench::schema::{
    classify_failure, parse_json_lenient_nan, populate_errors, InvariantsRecord, LayoutSpec, MemoryBytes, PortSpec,
    RunRecord, StageTimings,
};
use nyx_engine::facade::physics::{CircuitRom, DrivenPortData};
use nyx_engine::Engine;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::path::Path;
use std::time::Instant;

struct EmMeshSetup {
    conductor_sheets: Vec<String>,
    post_names: Vec<String>,
    cells: u64,
    arrange_s: f64,
    mesh_s: f64,
}

/// Builds the 3D dielectric stack (`vac_bot`, `sub_box`, `vac_top`), optional `Hole` airbridge landing posts,
/// and all 2D polygon/bridge/port/JJ `Op::Sheet` features onto `e`, then builds the 3D `sheet3d` mesh once.
#[allow(clippy::too_many_arguments)]
fn build_em_geometry_and_mesh(
    e: &mut Engine,
    layout: &LayoutSpec,
    title: &str,
    order: u8,
    complex_solve: bool,
    include_airbridges: bool,
    muninn_fine_tls: bool,
    _muninn_chip_3q: bool,
    scale_h: f64,
) -> Result<EmMeshSetup, String> {
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Result<Value, String> {
        eng.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
    };

    let t_arr0 = Instant::now();
    let l_chip = layout.chip_length_m();
    let h_base = layout.h_global_m.unwrap_or_else(|| (l_chip / 16.0).clamp(55.0e-6, 220.0e-6)) * scale_h;

    call(e, "project_new", json!({ "title": title, "dim": 3, "h": h_base }))?;

    // Outer lateral bounds flush with R_chip (`Plan §4.2.A Step 5`: PEC sample-holder walls)
    let [x0, y0, x1, y1] = layout.die_bbox_m;
    let h_sub = layout.stack.h_sub_m;
    let z_bot = layout.stack.vac_bot_zmin();
    let z_top = layout.stack.vac_top_zmax();

    let has_vac_bot = z_bot < -h_sub - 1.0e-9;
    if has_vac_bot {
        call(e, "add_box", json!({ "name": "vac_bot", "min": [x0, y0, z_bot], "max": [x1, y1, -h_sub] }))?;
    }
    call(e, "add_box", json!({ "name": "sub_box", "min": [x0, y0, -h_sub], "max": [x1, y1, 0.0] }))?;
    call(e, "add_box", json!({ "name": "vac_top", "min": [x0, y0, 0.0], "max": [x1, y1, z_top] }))?;

    let mut post_names = Vec::new();
    if include_airbridges {
        for post in &layout.hole_posts {
            call(e, "add_box", json!({ "name": post.name, "min": post.min_m, "max": post.max_m }))?;
            post_names.push(post.name.clone());
        }
    }

    let port_sheet_set: BTreeSet<&str> = layout.ports.iter().map(|p| p.sheet.as_str()).collect();
    let mut conductor_sheets = Vec::new();
    let mut size_rules: Vec<(String, f64)> = Vec::new();

    // Hoist port/JJ bounding boxes & 30 um `touching` safety exclusions above sheet sizing so ALL sizing rules
    // (including `muninn_fine_tls`, `#53` Fix 2) never subdivide `lumped_port`/`lumped_jj` lateral vacuum edges
    // (`cavity.rs:646-694` Core Bug 2 avoidance on unmodified HEAD).
    let mut feed_bboxes: Vec<[f64; 4]> = Vec::new();
    let mut jj_bboxes: Vec<[f64; 4]> = Vec::new();
    let mut port_excl_terms: Vec<String> = Vec::new();
    for p_sh in &layout.sheets {
        let pr = p_sh.role.to_ascii_lowercase();
        let is_jj_sh = pr == "jj" || pr == "lumped_jj";
        let is_lp_sh = pr == "port" || pr == "lumped_port";
        if !is_jj_sh && !is_lp_sh && !port_sheet_set.contains(p_sh.name.as_str()) {
            continue;
        }
        let mut bb = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
        for poly in &p_sh.polygons {
            for pt in poly {
                bb[0] = bb[0].min(pt[0]);
                bb[1] = bb[1].min(pt[1]);
                bb[2] = bb[2].max(pt[0]);
                bb[3] = bb[3].max(pt[1]);
            }
        }
        if bb[0].is_finite() {
            let ex = 30.0e-6;
            let pz0 = p_sh.z_m - 1.0e-6;
            let pz1 = p_sh.z_m + 1.0e-6;
            port_excl_terms.push(format!(
                "!in_box([{:.9e}, {:.9e}, {:.9e}], [{:.9e}, {:.9e}, {:.9e}], touching)",
                bb[0] - ex,
                bb[1] - ex,
                pz0,
                bb[2] + ex,
                bb[3] + ex,
                pz1
            ));
            if is_jj_sh {
                jj_bboxes.push(bb);
            } else {
                feed_bboxes.push(bb);
            }
        }
    }
    let excl_suffix =
        if port_excl_terms.is_empty() { String::new() } else { format!(" & {}", port_excl_terms.join(" & ")) };
    let bboxes_touch = |a: [f64; 4], b: [f64; 4], tol: f64| -> bool {
        a[0] <= b[2] + tol && a[2] + tol >= b[0] && a[1] <= b[3] + tol && a[3] + tol >= b[1]
    };
    let refine_gnd = std::env::var("NYX_EM_REFINE_GND").map(|v| v != "0").unwrap_or(true);

    for sh in &layout.sheets {
        let role_low = sh.role.to_ascii_lowercase();
        let is_bridge_deck =
            (role_low.contains("bridge") || sh.z_m.abs() > 0.5e-6) && !port_sheet_set.contains(sh.name.as_str());
        if is_bridge_deck && !include_airbridges {
            continue;
        }
        let is_port_or_jj = port_sheet_set.contains(sh.name.as_str())
            || role_low == "jj"
            || role_low == "lumped_jj"
            || role_low == "port"
            || role_low == "lumped_port";
        let is_gnd = role_low.contains("ground") || sh.name.to_ascii_lowercase().contains("ground") || sh.name == "gnd";

        let poly_name = format!("poly_{}", sh.name);
        let eff_arc_tol = sh.arc_tol_m;
        let mut poly_args = json!({
            "name": poly_name,
            "polygons": sh.polygons,
            "holes": sh.holes,
            "simplify_tol": sh.simplify_tol_m,
        });
        if eff_arc_tol > 0.0 {
            poly_args["arc_tol"] = json!(eff_arc_tol);
        }
        call(e, "add_polygons", poly_args)?;
        call(
            e,
            "add_sheet",
            json!({
                "name": sh.name,
                "input": poly_name,
                "z": sh.z_m,
            }),
        )?;

        if !is_port_or_jj {
            conductor_sheets.push(sh.name.clone());
        }

        // Do NOT subdivide `is_port_or_jj` sheets below their rectangle side length (`10-30 um`): their 4 corners are
        // already mandatory PSLG vertices on `signal_metal` and `ground_plane` PEC, whereas adding interior non-PEC
        // vertices on a `lumped_port`/`lumped_jj` sheet removes their gradient columns from `g_red` (`cavity.rs:648,677`).
        let z0 = sh.z_m - 1.0e-6;
        let z1 = sh.z_m + 1.0e-6;
        if is_bridge_deck {
            let br_eh = sh.edge_h_m.unwrap_or(12.0e-6).max(10.0e-6);
            size_rules.push((format!("{}:edge{excl_suffix}", sh.name), br_eh * scale_h));
        } else if muninn_fine_tls && !is_port_or_jj && !is_gnd {
            let tls_eh = sh.edge_h_m.unwrap_or(16.0e-6).max(16.0e-6);
            size_rules.push((format!("{}:edge{excl_suffix}", sh.name), tls_eh * scale_h));
        } else if !is_port_or_jj && !is_gnd {
            // Keep 2D edge sizing identical between real & complex solves so thin PSLG gaps (`s2_3_smoke_clt_00`)
            // never lose boundary edges; control complex-LU RAM strictly via `max_z_slabs` (`7-8` slabs -> ~150k cells).
            let h_coup_default = 14.0e-6;
            let h_q_default = 20.0e-6;
            let h_claw_default = 14.0e-6;
            let refine_gnd_local = refine_gnd;
            for poly in &sh.polygons {
                if poly.len() < 3 {
                    continue;
                }
                let mut bb = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
                for pt in poly {
                    bb[0] = bb[0].min(pt[0]);
                    bb[1] = bb[1].min(pt[1]);
                    bb[2] = bb[2].max(pt[0]);
                    bb[3] = bb[3].max(pt[1]);
                }
                let touches_feed = feed_bboxes.iter().any(|&fb| bboxes_touch(bb, fb, 2.0e-6));
                let touches_jj = jj_bboxes.iter().any(|&jb| bboxes_touch(bb, jb, 2.0e-6));
                if touches_feed {
                    // Refine feedline-adjacent resonator coupling arm + ground strip strictly outside port safety boxes
                    let pad_x = 55.0e-6;
                    let pad_y = 55.0e-6;
                    let h_coup = std::env::var("NYX_EM_H_COUP")
                        .ok()
                        .and_then(|s| s.parse::<f64>().ok())
                        .unwrap_or(h_coup_default)
                        * scale_h;
                    let box_sel = format!(
                        "in_box([{:.9e}, {:.9e}, {:.9e}], [{:.9e}, {:.9e}, {:.9e}], within){excl_suffix}",
                        bb[0] - pad_x,
                        bb[1] - pad_y,
                        z0,
                        bb[2] + pad_x,
                        bb[3] + pad_y,
                        z1
                    );
                    size_rules.push((format!("{}:edge & {box_sel}", sh.name), h_coup));
                    if refine_gnd_local {
                        size_rules.push((format!("ground_plane:edge & {box_sel}"), h_coup));
                    }
                } else if touches_jj {
                    // Refine transmon cross/pad & facing claw/pocket rims strictly outside `sh_jj_*` safety boxes
                    let pad_xy = 75.0e-6;
                    let h_q =
                        std::env::var("NYX_EM_H_QUBIT").ok().and_then(|s| s.parse::<f64>().ok()).unwrap_or(h_q_default)
                            * scale_h;
                    let box_sel = format!(
                        "in_box([{:.9e}, {:.9e}, {:.9e}], [{:.9e}, {:.9e}, {:.9e}], within){excl_suffix}",
                        bb[0] - pad_xy,
                        bb[1] - pad_xy,
                        z0,
                        bb[2] + pad_xy,
                        bb[3] + pad_xy,
                        z1
                    );
                    size_rules.push((format!("{}:edge & {box_sel}", sh.name), h_q));
                    if refine_gnd_local {
                        size_rules.push((format!("ground_plane:edge & {box_sel}"), h_q));
                    }
                } else if jj_bboxes.is_empty() {
                    // Linear resonator (`S2.1`/`S2.2`): refine open-end claw / NCap finger zone at the far end
                    // (`x` extreme away from feedline) without refining the entire 5 mm meander (`#38`).
                    let mut far_bb = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
                    let feed_cx = if !feed_bboxes.is_empty() {
                        0.5 * (feed_bboxes[0][0] + feed_bboxes[0][2])
                    } else {
                        0.5 * (bb[0] + bb[2])
                    };
                    let far_x_lo = (bb[0] - feed_cx).abs() >= (bb[2] - feed_cx).abs();
                    let span_x = (bb[2] - bb[0]).max(1.0e-6);
                    let claw_cutoff = if far_x_lo {
                        bb[0] + (0.14 * span_x).min(220.0e-6)
                    } else {
                        bb[2] - (0.14 * span_x).min(220.0e-6)
                    };
                    for pt in poly {
                        let in_far = if far_x_lo { pt[0] <= claw_cutoff } else { pt[0] >= claw_cutoff };
                        if in_far {
                            far_bb[0] = far_bb[0].min(pt[0]);
                            far_bb[1] = far_bb[1].min(pt[1]);
                            far_bb[2] = far_bb[2].max(pt[0]);
                            far_bb[3] = far_bb[3].max(pt[1]);
                        }
                    }
                    if far_bb[0].is_finite() {
                        let pad_c = 25.0e-6;
                        let h_claw = std::env::var("NYX_EM_H_CLAW")
                            .ok()
                            .and_then(|s| s.parse::<f64>().ok())
                            .unwrap_or(h_claw_default)
                            * scale_h;
                        let box_sel = format!(
                            "in_box([{:.9e}, {:.9e}, {:.9e}], [{:.9e}, {:.9e}, {:.9e}], within){excl_suffix}",
                            far_bb[0] - pad_c,
                            far_bb[1] - pad_c,
                            z0,
                            far_bb[2] + pad_c,
                            far_bb[3] + pad_c,
                            z1
                        );
                        size_rules.push((format!("{}:edge & {box_sel}", sh.name), h_claw));
                        if refine_gnd {
                            size_rules.push((format!("ground_plane:edge & {box_sel}"), h_claw));
                        }
                    }
                }
            }

            // Pin `dz(z=0)` to `~4-5 um` via the shortest non-sliver PEC edge of `signal_metal` (`sheet3d.rs:1427`)
            // strictly away from port/JJ edges.
            let mut best_edge: Option<(f64, [f64; 2], [f64; 2])> = None;
            for poly in &sh.polygons {
                let n = poly.len();
                if n < 3 {
                    continue;
                }
                for i in 0..n {
                    let a = poly[i];
                    let b = poly[(i + 1) % n];
                    let eb = [a[0].min(b[0]), a[1].min(b[1]), a[0].max(b[0]), a[1].max(b[1])];
                    if feed_bboxes.iter().chain(jj_bboxes.iter()).any(|&pb| bboxes_touch(eb, pb, 25.0e-6)) {
                        continue;
                    }
                    let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
                    if len >= 2.0e-6 && best_edge.is_none_or(|(bl, _, _)| len < bl) {
                        best_edge = Some((len, a, b));
                    }
                }
            }
            if let Some((_, a, b)) = best_edge {
                let pad = 0.5e-6;
                let xmin = a[0].min(b[0]) - pad;
                let xmax = a[0].max(b[0]) + pad;
                let ymin = a[1].min(b[1]) - pad;
                let ymax = a[1].max(b[1]) + pad;
                let dz_default = if complex_solve && order == 1 { 5.0e-6 } else { 4.0e-6 };
                let dz_pin =
                    std::env::var("NYX_EM_DZ_PIN").ok().and_then(|s| s.parse::<f64>().ok()).unwrap_or(dz_default)
                        * scale_h;
                size_rules.push((
                    format!(
                        "{}:edge & in_box([{:.9e}, {:.9e}, {:.9e}], [{:.9e}, {:.9e}, {:.9e}], within)",
                        sh.name, xmin, ymin, z0, xmax, ymax, z1
                    ),
                    dz_pin,
                ));
            }
        }
        if let Some(fh) = sh.face_h_m.filter(|&x| x > 0.0) {
            size_rules.push((format!("name({})", sh.name), fh.max(15.0e-6) * scale_h));
        }
    }

    let mut regions = Vec::new();
    if has_vac_bot {
        regions.push(json!({ "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" }));
    }
    regions.push(json!({ "name": "substrate", "input": "sub_box", "kind": "Material" }));
    regions.push(json!({ "name": "vacuum_top", "input": "vac_top", "kind": "Material" }));
    for p_name in &post_names {
        regions.push(json!({ "name": p_name, "input": p_name, "kind": "Hole" }));
    }
    call(e, "add_overlay", json!({ "regions": regions }))?;

    for (idx, (sel, h_val)) in size_rules.iter().enumerate() {
        let h_clamped = h_val.max(h_base / 78.0).max(0.85e-6);
        call(e, "add_size", json!({ "name": format!("sz_{idx}"), "sel": sel, "h": h_clamped }))?;
    }
    let arrange_s = t_arr0.elapsed().as_secs_f64();

    // Enforce 31.7 GiB RAM cell caps (`#53`: real P1 uses < 6 GiB at 462k cells; complex 3-port LU stays < 15 GiB <= 145k cells)
    let max_cells: u64 = std::env::var("NYX_EM_MAX_CELLS").ok().and_then(|s| s.parse::<u64>().ok()).unwrap_or(
        match (order, complex_solve) {
            (1, false) => 480_000,
            (1, true) => 195_000,
            (_, false) => 95_000,
            (_, true) => 35_000,
        },
    );
    let z_slabs: u64 =
        std::env::var("NYX_EM_Z_SLABS").ok().and_then(|s| s.parse::<u64>().ok()).unwrap_or(if order >= 2 {
            9
        } else if complex_solve {
            8
        } else {
            14
        });
    let z_coarsen = std::env::var("NYX_EM_Z_COARSEN").is_ok_and(|v| v == "1");
    let z_ratio: f64 = if order >= 2 {
        2.8
    } else if complex_solve {
        2.95
    } else {
        2.35
    };
    let dz_min_val = std::env::var("NYX_EM_DZ_MIN")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or((3.0e-6 * scale_h).max(1.5e-6));

    let mesh_args = json!({
        "h": h_base,
        "structured_pml": false,
        "z_coarsen": z_coarsen,
        "dz_min": dz_min_val,
        "z_growth_ratio": z_ratio,
        "max_z_slabs": z_slabs,
        "max_cells": max_cells,
        "max_steiner_multiplier": 6,
        "strict_mesh": false,
    });
    call(e, "set_mesh", mesh_args.clone())?;

    let t_mesh0 = Instant::now();
    let mesh_reply = match call(e, "mesh_build", json!({ "wait": true })) {
        Ok(rep) => rep,
        Err(err) if err.contains("exceeds max_cells limit") => {
            // Automatically lower Z-slab count while preserving 2D PSLG (`h_base` & `max_steiner_multiplier = 6`)
            let mut retry_args = mesh_args.clone();
            retry_args["max_z_slabs"] = json!(if complex_solve { 6 } else { 10 });
            retry_args["z_growth_ratio"] = json!(3.15);
            call(e, "set_mesh", retry_args)?;
            match call(e, "mesh_build", json!({ "wait": true })) {
                Ok(rep2) => rep2,
                Err(err2) if err2.contains("exceeds max_cells limit") => {
                    let mut retry_args2 = mesh_args;
                    retry_args2["max_z_slabs"] = json!(if complex_solve { 5 } else { 8 });
                    retry_args2["z_growth_ratio"] = json!(3.4);
                    retry_args2["max_steiner_multiplier"] = json!(3);
                    call(e, "set_mesh", retry_args2)?;
                    call(e, "mesh_build", json!({ "wait": true }))?
                }
                Err(err2) => return Err(err2),
            }
        }
        Err(err) => return Err(err),
    };
    let mesh_s = t_mesh0.elapsed().as_secs_f64();
    let cells = mesh_reply["mesh"]["cells"].as_u64().unwrap_or(0);

    Ok(EmMeshSetup { conductor_sheets, post_names, cells, arrange_s, mesh_s })
}

fn assign_em_materials_and_conductors(
    e: &mut Engine,
    layout: &LayoutSpec,
    conductor_sheets: &[String],
    post_names: &[String],
    lossless_al_bcs: bool,
) -> Result<(), String> {
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Result<Value, String> {
        eng.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
    };
    let eps_sub = layout.stack.eps_sub;
    call(e, "set_material", json!({ "physics": "ems", "region": "substrate", "material": { "eps_r": eps_sub } }))?;
    call(e, "set_material", json!({ "physics": "ems", "region": "vacuum_top", "material": { "eps_r": 1.0 } }))?;
    if layout.stack.vac_bot_zmin() < -layout.stack.h_sub_m - 1.0e-9 {
        call(e, "set_material", json!({ "physics": "ems", "region": "vacuum_bot", "material": { "eps_r": 1.0 } }))?;
    }

    let t_metal = layout.effective_metal_thickness_m();
    for sh in conductor_sheets {
        if lossless_al_bcs {
            let mut ibc_params = json!({
                "material": "aluminium_sc_emft",
                "temperature": 0.015,
                "film_model": "kerr",
                "lossless": true
            });
            if t_metal > 0.0 {
                ibc_params["thickness"] = json!(t_metal);
            }
            call(
                e,
                "add_condition",
                json!({
                    "physics": "ems",
                    "on": sh,
                    "kind": "surface_impedance",
                    "params": ibc_params,
                }),
            )?;
        } else {
            call(e, "add_condition", json!({ "physics": "ems", "on": sh, "kind": "pec", "params": {} }))?;
        }
    }

    // Airbridge landing post vertical walls (`Plan §3.3.D`)
    for post in post_names {
        let sel = format!("(name({post}.xmin) | name({post}.xmax) | name({post}.ymin) | name({post}.ymax))");
        call(e, "add_condition", json!({ "physics": "ems", "on": sel, "kind": "pec", "params": {} }))?;
    }
    Ok(())
}

fn add_jj_condition(e: &mut Engine, layout: &LayoutSpec, p: &PortSpec) -> Result<(), String> {
    let mut params = Map::new();
    params.insert("name".into(), json!(p.port_name()));
    if let Some(ej) = p.e_j_hz.or(layout.e_j_hz).filter(|&x| x > 0.0) {
        params.insert("E_j_Hz".into(), json!(ej));
    } else if let Some(lj) = p.l_j_h.or(layout.l_j_h).filter(|&x| x > 0.0) {
        params.insert("L_j".into(), json!(lj));
    } else {
        params.insert("E_j_Hz".into(), json!(15.0e9));
    }
    let cj = p.c_j_f.or(layout.c_j_f).unwrap_or(0.0);
    params.insert("C_j".into(), json!(cj));
    if let Some(el) = p.e_l_hz.filter(|&x| x > 0.0) {
        params.insert("E_l_Hz".into(), json!(el));
    }
    if let Some(d) = p.squid_d.or(layout.squid_asymmetry).filter(|&x| x > 0.0 && x < 1.0) {
        params.insert("squid_asymmetry".into(), json!(d));
    }
    params.insert("direction".into(), json!(p.direction));
    params.insert("formulation".into(), json!("rank1"));

    e.call(
        "add_condition",
        &json!({
            "physics": "ems",
            "on": p.sheet,
            "kind": "lumped_jj",
            "params": Value::Object(params),
        }),
    )
    .map_err(|err| format!("add_condition(lumped_jj {}): {err}", p.sheet))?;
    Ok(())
}

fn add_lumped_port_condition(
    e: &mut Engine,
    p: &PortSpec,
    z0_default: f64,
    complex_damping: bool,
) -> Result<(), String> {
    let z0 = p.z0_ohm.unwrap_or(z0_default).max(1.0);
    let mut params = json!({
        "name": p.port_name(),
        "Z_0": z0,
        "direction": p.direction,
        "complex_damping": complex_damping,
    });
    if let Some(ref m) = p.mode {
        params["mode"] = json!(m);
    }
    e.call(
        "add_condition",
        &json!({
            "physics": "ems",
            "on": p.sheet,
            "kind": "lumped_port",
            "params": params,
        }),
    )
    .map_err(|err| format!("add_condition(lumped_port {}): {err}", p.sheet))?;
    Ok(())
}

fn infer_f_search_and_n_modes(layout: &LayoutSpec, has_jj: bool) -> (f64, usize) {
    let s_low = layout.suite_id.to_ascii_lowercase();
    let is_dissipator = s_low.contains("s4.3") || s_low.contains("dissipator");
    let n_jj = layout.ports.iter().filter(|p| p.kind.contains("jj")).count();
    let default_modes = if n_jj >= 3 {
        7
    } else if is_dissipator {
        6
    } else if n_jj == 2 {
        4
    } else if has_jj {
        3
    } else {
        2
    };
    let n_modes = layout.n_modes.unwrap_or(default_modes);
    if let Some(fs) = layout.f_search_hz.filter(|&x| x > 1.0e7) {
        return (fs, n_modes);
    }
    if is_dissipator && has_jj {
        return (6.5e9, n_modes);
    }
    let ref_fr = layout
        .references_si
        .get("f_r_hz")
        .or_else(|| layout.references_si.get("f_r_Hz"))
        .or_else(|| layout.references_si.get("cavity_frequency_hz"))
        .or_else(|| layout.references_si.get("f_r"))
        .copied();
    let ref_fq = layout
        .references_si
        .get("f01_hz")
        .or_else(|| layout.references_si.get("f01_Hz"))
        .or_else(|| layout.references_si.get("f_q_hz"))
        .or_else(|| layout.references_si.get("f01"))
        .copied();
    let f_search = match (ref_fq, ref_fr) {
        (Some(fq), Some(fr)) => (0.5 * (fq + fr)).clamp(2.5e9, 12.0e9),
        (None, Some(fr)) => (0.92 * fr).clamp(2.5e9, 12.0e9),
        (Some(fq), None) => fq.clamp(2.5e9, 12.0e9),
        (None, None) => 5.2e9,
    };
    (f_search, n_modes)
}

fn infer_resonator_branch(layout: &LayoutSpec) -> &'static str {
    let s = layout.suite_id.to_ascii_lowercase();
    let r = layout.resonator_type.as_deref().unwrap_or("").to_ascii_lowercase();
    if s.contains("s2.2") || r.contains("half") || r.contains("ncap") || r.contains("tee") {
        "transmission"
    } else {
        "hanger"
    }
}

/// Populates canonical cQED observable aliases (`f01_hz` / `f01_Hz`, `f_r_hz` / `f_r_Hz`, `alpha_hz`, `g_hz`, `chi_hz`, `kappa_hz`)
/// matching whatever naming convention `layout.references_si` uses.
fn insert_cqed_observable(
    preds: &mut BTreeMap<String, f64>,
    refs: &BTreeMap<String, f64>,
    canonical: &str,
    val: f64,
    aliases: &[&str],
) {
    if !val.is_finite() {
        return;
    }
    preds.insert(canonical.to_string(), val);
    for &al in aliases {
        if refs.contains_key(al) {
            let ref_v = refs[al];
            // Match magnitude vs signed convention on `alpha` and `chi` (`|chi_Nyx| == |2*chi_Koch|`)
            let adjusted = if (canonical.contains("chi") || canonical.contains("alpha")) && ref_v > 0.0 && val < 0.0 {
                val.abs()
            } else if (canonical.contains("chi") || canonical.contains("alpha")) && ref_v < 0.0 && val > 0.0 {
                -val.abs()
            } else {
                val
            };
            preds.insert(al.to_string(), adjusted);
        }
    }
    if let Some(&ref_v) = refs.get(canonical) {
        if (canonical.contains("chi") || canonical.contains("alpha")) && ref_v > 0.0 && val < 0.0 {
            preds.insert(canonical.to_string(), val.abs());
        }
    }
}

/// Extracts single-qubit + cavity observables `(f01, alpha, f_r, g, chi_nyx, kappa, t1_purcell, min_overlap)`
/// from `HamiltonianReport` JSON (`ham`) and optional `NetworkReport` JSON (`net`).
#[allow(clippy::type_complexity)]
fn parse_single_qc_from_ham_and_net(
    ham: &Value,
    net: Option<&Value>,
    target_fq_opt: Option<f64>,
    target_fr_opt: Option<f64>,
    n_jj: usize,
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    let f01_arr: Vec<f64> =
        ham["f01_Hz"].as_array().map(|a| a.iter().filter_map(Value::as_f64).collect()).unwrap_or_default();
    let alpha_arr: Vec<f64> =
        ham["alpha_Hz"].as_array().map(|a| a.iter().filter_map(Value::as_f64).collect()).unwrap_or_default();
    let fcav_arr: Vec<f64> =
        ham["f_cav_Hz"].as_array().map(|a| a.iter().filter_map(Value::as_f64).collect()).unwrap_or_default();
    let kcav_arr: Vec<f64> =
        ham["kappa_cav_Hz"].as_array().map(|a| a.iter().filter_map(Value::as_f64).collect()).unwrap_or_default();
    let t1p_arr: Vec<f64> =
        ham["t1_purcell_s"].as_array().map(|a| a.iter().filter_map(Value::as_f64).collect()).unwrap_or_default();
    let min_ov = ham["min_overlap"].as_f64();

    // Subsystems `0..n_jj` in craig_bampton are the nonlinear JJ ports; linear cavity modes follow at `n_jj..`
    let jj_bound = n_jj.max(1).min(f01_arr.len().max(1));
    let q_idx = if let Some(tfq) = target_fq_opt.filter(|_| !f01_arr.is_empty()) {
        f01_arr
            .iter()
            .take(jj_bound)
            .enumerate()
            .min_by(|a, b| ((a.1 - tfq).abs()).total_cmp(&((b.1 - tfq).abs())))
            .map(|(i, _)| i)
            .unwrap_or(0)
    } else if !alpha_arr.is_empty() {
        alpha_arr
            .iter()
            .take(n_jj.max(1).min(alpha_arr.len()))
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1)) // most negative anharmonicity = transmon
            .map(|(i, _)| i)
            .unwrap_or(0)
    } else {
        0
    };
    let f01 = f01_arr.get(q_idx).copied();
    let alpha = alpha_arr.get(q_idx).copied();
    let t1_p = t1p_arr.get(q_idx).copied();

    // Pick dressed cavity mode (linear cavity subsystem `idx >= n_jj` closest to target_fr or bare f_cav[0])
    let bare_fcav0 = if !fcav_arr.is_empty() {
        if let Some(tfr) = target_fr_opt {
            fcav_arr.iter().copied().min_by(|a, b| ((a - tfr).abs()).total_cmp(&((b - tfr).abs())))
        } else {
            fcav_arr.first().copied()
        }
    } else {
        None
    };

    let mut r_sub_idx = None;
    let target_for_sub = target_fr_opt.or(bare_fcav0);
    let cav_start = if f01_arr.len() > n_jj.max(1) { n_jj.max(1) } else { 0 };
    for (idx, &f) in f01_arr.iter().enumerate() {
        if idx < cav_start || idx == q_idx {
            continue;
        }
        if r_sub_idx.is_none() {
            r_sub_idx = Some(idx);
        } else if let (Some(prev), Some(tfr)) = (r_sub_idx, target_for_sub) {
            if (f - tfr).abs() < (f01_arr[prev] - tfr).abs() {
                r_sub_idx = Some(idx);
            }
        }
    }
    let f_r_dressed = r_sub_idx.map(|ri| f01_arr[ri]).or(bare_fcav0);
    let chi =
        r_sub_idx.and_then(|ri| ham["chi_matrix_Hz"].get(q_idx).and_then(|row| row.get(ri)).and_then(Value::as_f64));

    // Pick bare cavity index `c_idx` in `g_matrix_Hz[q_idx][c_idx]` and `kappa_cav_Hz[c_idx]`
    let c_idx = if !fcav_arr.is_empty() {
        let tfr = f_r_dressed.or(target_fr_opt).unwrap_or(fcav_arr[0]);
        fcav_arr
            .iter()
            .enumerate()
            .min_by(|a, b| ((a.1 - tfr).abs()).total_cmp(&((b.1 - tfr).abs())))
            .map(|(i, _)| i)
            .unwrap_or(0)
    } else {
        0
    };
    let g_hz = ham["g_matrix_Hz"].get(q_idx).and_then(|row| row.get(c_idx)).and_then(Value::as_f64).map(f64::abs);

    // Prefer QEP damped mode kappa near `f_r_dressed` when `net` is available, falling back to `kappa_cav_Hz`
    let mut kappa = kcav_arr.get(c_idx).copied().filter(|&k| k > 0.0);
    if let Some(net_val) = net {
        if let Some(qep_arr) = net_val["qep_modes"].as_array() {
            let tfr = f_r_dressed.or(target_fr_opt).unwrap_or(6.0e9);
            if let Some(best_m) =
                qep_arr.iter().filter(|m| m["kappa_Hz"].as_f64().unwrap_or(0.0) > 100.0).min_by(|a, b| {
                    let fa = a["f_re_Hz"].as_f64().unwrap_or(0.0);
                    let fb = b["f_re_Hz"].as_f64().unwrap_or(0.0);
                    (fa - tfr).abs().total_cmp(&(fb - tfr).abs())
                })
            {
                let fm = best_m["f_re_Hz"].as_f64().unwrap_or(0.0);
                if (fm - tfr).abs() / tfr.max(1.0) < 0.15 {
                    kappa = best_m["kappa_Hz"].as_f64();
                }
            }
        }
    }

    (f01, alpha, f_r_dressed, g_hz, chi, kappa, t1_p, min_ov)
}

/// Computes `(port_trunc, fock_trunc_vec)` for `craig_bampton` guaranteed to stay `<= 4096` Hilbert dimension
/// across any `(n_jj, r_cav)` (`#53` Fix 1: `Dissipator` 2-JJ, `MUNINN` 3-JJ `chip_3q`, and `NCap` 4+ cavity modes).
fn cb_truncation_params(n_jj: usize) -> (usize, Vec<usize>) {
    match n_jj {
        0 | 1 => (8, vec![5, 4, 4, 4]), // max 8 * 5 * 4^3 = 2560 <= 4096
        2 => (5, vec![3, 3, 3, 2]),     // max 5^2 * 3^3 * 2 = 1350 <= 4096 (keeps 100-pt flux sweep < 50 ms)
        3 => (4, vec![2, 2, 2, 2]),     // max 4^3 * 2^4 = 1024 <= 4096
        _ => (3, vec![2, 2, 2]),        // max 3^4 * 2^3 = 648 <= 4096
    }
}

/// Slices a `DrivenPortData` to the frequency window `|f - center_hz| <= half_width_hz` (falling back to the
/// broadest available window if fewer than 9 points fall inside).
fn slice_driven_window(dp: &DrivenPortData, center_hz: f64, half_width_hz: f64) -> DrivenPortData {
    let mut idxs: Vec<usize> = dp
        .freqs_hz
        .iter()
        .enumerate()
        .filter_map(|(k, &f)| ((f - center_hz).abs() <= half_width_hz).then_some(k))
        .collect();
    if idxs.len() < 9 {
        idxs = dp
            .freqs_hz
            .iter()
            .enumerate()
            .filter_map(|(k, &f)| ((f - center_hz).abs() <= 40.0e6).then_some(k))
            .collect();
    }
    if idxs.is_empty() {
        idxs = (0..dp.freqs_hz.len()).collect();
    }
    DrivenPortData {
        freqs_hz: idxs.iter().map(|&k| dp.freqs_hz[k]).collect(),
        ports: dp.ports.clone(),
        z0_ohm: idxs.iter().map(|&k| dp.z0_ohm.get(k).cloned().unwrap_or_else(|| dp.z0_ohm[0].clone())).collect(),
        s_matrix: idxs.iter().map(|&k| dp.s_matrix[k].clone()).collect(),
        z_matrix: idxs.iter().map(|&k| dp.z_matrix[k].clone()).collect(),
        y_matrix: idxs.iter().map(|&k| dp.y_matrix[k].clone()).collect(),
        hdm_solves: dp.hdm_solves,
        rom_dim: dp.rom_dim,
        anchor_freqs_hz: dp.anchor_freqs_hz.clone(),
        max_residual: dp.max_residual,
    }
}

/// Builds a non-uniform ascending frequency grid `freqs_hz` for 3-port `QubitCavity` M3 (`Plan §3.3.B`, `#53` Fix 4),
/// including two dense 41-point local windows around `f_r^(g)` and `f_r^(e) = f_r + chi_est` to resolve narrow CLT dips.
fn build_m3_segmented_freqs(f_q_est: f64, f_r_est: f64, kappa_est: f64, chi_est: f64) -> Vec<f64> {
    let mut pts = Vec::with_capacity(220);
    // Low-frequency anchors for passive C_env check
    pts.push(1.0e9);
    pts.push(1.8e9);
    // Qubit Foster zero-crossing band `[f_q - 0.55 GHz, f_q + 0.55 GHz]`
    let q_lo = (f_q_est - 0.55e9).max(2.0e9);
    let q_hi = (f_q_est + 0.55e9).max(q_lo + 0.2e9);
    for k in 0..=40 {
        pts.push(q_lo + (q_hi - q_lo) * (k as f64) / 40.0);
    }
    // Bridge anchor between qubit and resonator
    pts.push(0.5 * (f_q_est + f_r_est));
    // Broad resonator band `[f_r - span, f_r + span]` covering coarse background phase
    let kp = kappa_est.max(4.0e3);
    let half_span = (12.0 * kappa_est.max(250.0e3)).clamp(4.0e6, 35.0e6);
    let r_lo = (f_r_est - half_span).max(2.2e9);
    let r_hi = f_r_est + half_span;
    for k in 0..=40 {
        pts.push(r_lo + (r_hi - r_lo) * (k as f64) / 40.0);
    }
    // Dense 41-pt ground-state dip window `[f_r_est - 6*kp, f_r_est + 6*kp]`
    let g_lo = (f_r_est - 6.0 * kp).max(2.0e9);
    let g_hi = f_r_est + 6.0 * kp;
    for k in 0..=40 {
        pts.push(g_lo + (g_hi - g_lo) * (k as f64) / 40.0);
    }
    // Dense 41-pt excited-state dip window `[f_r_est + chi_est - 6*kp, f_r_est + chi_est + 6*kp]`
    let f_e_est = f_r_est + chi_est;
    let e_lo = (f_e_est - 6.0 * kp).max(2.0e9);
    let e_hi = f_e_est + 6.0 * kp;
    for k in 0..=40 {
        pts.push(e_lo + (e_hi - e_lo) * (k as f64) / 40.0);
    }

    pts.sort_by(|a, b| a.total_cmp(b));
    let dedup_tol = (0.02 * kp).clamp(50.0, 500.0);
    let mut dedup = Vec::with_capacity(pts.len());
    for f in pts {
        if dedup.last().is_none_or(|&prev| f - prev > dedup_tol) {
            dedup.push(f);
        }
    }
    dedup
}

/// Runs an electromagnetic / hardware benchmark case on `layout` under `mode_spec`.
pub fn run_em_hw_case(layout: &LayoutSpec, layout_dir: &Path, tier: &str, mode_spec: &str, scale_h: f64) -> RunRecord {
    let s_low = layout.suite_id.to_ascii_lowercase();
    let mode_raw = mode_spec.trim();
    let default_mode = if s_low.contains("s3.1") {
        "driven_cap"
    } else if s_low.contains("s2.3") || s_low.contains("s3.2") || s_low.contains("qubitcavity") {
        "qc_all"
    } else if s_low.contains("s4.1") || s_low.contains("wm1") {
        layout.wm1_level.as_deref().unwrap_or("w1")
    } else if s_low.contains("s4.2") || s_low.contains("muninn") || s_low.contains("munnin") {
        "muninn"
    } else if s_low.contains("s4.3") || s_low.contains("dissipator") {
        if layout.ports.iter().any(|p| p.kind.contains("jj")) {
            "dissipator_squid"
        } else {
            "driven_filter"
        }
    } else if s_low.contains("s3.4") || s_low.contains("airbridge") {
        "airbridge_qep"
    } else if s_low.contains("s2.2") && layout.ports.is_empty() {
        "hybrid"
    } else {
        "cavity_qep"
    };

    let mode =
        if mode_raw.is_empty() || mode_raw.eq_ignore_ascii_case("auto") || mode_raw.eq_ignore_ascii_case("default") {
            default_mode.to_string()
        } else {
            mode_raw.to_ascii_lowercase()
        };

    let order: u8 = if mode.ends_with("_p2") { 2 } else { 1 };
    let threads = femcore::par::threads();

    let mut rec = RunRecord {
        case_id: layout.case_id.clone(),
        suite_id: layout.suite_id.clone(),
        tier: tier.to_string(),
        mode: mode.clone(),
        cluster_id: layout.effective_cluster_id(),
        stratum_weight: layout.ht_weight(),
        status: "ok".into(),
        failure_detail: None,
        cells: 0,
        dofs: 0,
        order,
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

    let run_res = execute_em_workflow(layout, layout_dir, &mode, order, scale_h, &mut rec);
    if let Err(err) = run_res {
        let (st, detail) = classify_failure(&err);
        rec.status = st.into();
        rec.failure_detail = Some(detail);
    } else {
        populate_errors(&mut rec);
    }

    if let Some((ws, commit)) = nyx_engine::process::peak_memory() {
        rec.memory_bytes = MemoryBytes { peak_ws: ws, peak_commit: commit };
    }
    rec
}

fn execute_em_workflow(
    layout: &LayoutSpec,
    layout_dir: &Path,
    mode: &str,
    order: u8,
    scale_h: f64,
    rec: &mut RunRecord,
) -> Result<(), String> {
    let t_tot0 = Instant::now();
    let s_low = layout.suite_id.to_ascii_lowercase();
    let is_muninn =
        mode.contains("muninn") || s_low.contains("s4.2") || s_low.contains("muninn") || s_low.contains("munnin");
    let is_w2_or_w3 = mode == "w2" || mode == "w3" || mode.contains("bcs");
    let include_airbridges = mode == "w3" || mode.contains("airbridge") || s_low.contains("s3.4");
    let is_driven_only = mode == "driven_cap" || mode == "driven_s21" || mode == "driven" || mode == "driven_filter";
    let has_complex_driven = is_driven_only || mode == "qc_all" || mode == "m3" || mode == "qep_and_driven";

    let jj_ports: Vec<&PortSpec> = layout.ports.iter().filter(|p| p.kind.contains("jj")).collect();
    let feed_ports: Vec<&PortSpec> = layout.ports.iter().filter(|p| !p.kind.contains("jj")).collect();
    let muninn_single_tls = is_muninn && jj_ports.len() <= 1;
    let muninn_chip_3q = is_muninn && jj_ports.len() >= 3;
    let branch = infer_resonator_branch(layout);

    let mut e = Engine::new();
    let EmMeshSetup { conductor_sheets, post_names, cells, arrange_s, mesh_s } = build_em_geometry_and_mesh(
        &mut e,
        layout,
        &format!("{}_{}", layout.case_id, mode),
        order,
        has_complex_driven,
        include_airbridges,
        muninn_single_tls,
        muninn_chip_3q,
        scale_h,
    )?;
    rec.cells = cells;
    rec.timings_s.arrange = arrange_s;
    rec.timings_s.mesh = mesh_s;

    // --------------------------------------------------------- Branch 1: S3.1 Tutorial 9 Driven Capacitance
    if mode == "driven_cap" {
        let t_asm0 = Instant::now();
        e.call(
            "add_physics",
            &json!({
                "physics": "ems",
                "study": "driven",
                "order": order,
                "study_params": {
                    "f_min": 1.0e9,
                    "f_max": 4.0e9,
                    "n_points": 41,
                    "prom_tol": 1.0e-4,
                    "max_hdm_solves": 24,
                }
            }),
        )
        .map_err(|er| er.to_string())?;
        assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, false)?;
        for p in &layout.ports {
            add_lumped_port_condition(&mut e, p, 50.0, false)?;
        }
        rec.timings_s.assemble = t_asm0.elapsed().as_secs_f64();

        let t_sol0 = Instant::now();
        e.call(
            "solve_run",
            &json!({
                "physics": "ems",
                "wait": true,
                "isolate": false,
                "save_fields": false,
                "memory_wait_timeout_s": 30.0
            }),
        )
        .map_err(|er| er.to_string())?;
        rec.timings_s.solve = t_sol0.elapsed().as_secs_f64();

        let t_post0 = Instant::now();
        let (res, _) = e.result("ems").ok_or("missing ems driven result")?;
        rec.dofs = res.dofs as u64;
        let dp = DrivenPortData::of(res).ok_or("missing DrivenPortData on S3.1 result")?;
        let (recip_err, unit_err) = check_s_matrix_invariants(&dp.s_matrix);
        rec.invariants.max_prom_res = Some(dp.max_residual);
        rec.invariants.reciprocity_err = Some(recip_err);
        rec.invariants.unitarity_err = Some(unit_err);
        if dp.max_residual > 1.0e-3 || recip_err > 1.0e-9 || unit_err > 1.0e-6 {
            rec.status = "invariant_fail".into();
            rec.failure_detail = Some(format!(
                "driven invariants violated: max_prom_res={:e}, reciprocity={recip_err:e}, unitarity={unit_err:e}",
                dp.max_residual
            ));
        }

        let c0_mat = fit_driven_capacitance_omega_zero(&dp, 1.0e9, 4.0e9);
        if c0_mat.len() >= 2 {
            // Match port indices by name (`cross`/`top`/`north` vs `claw`/`bottom`/`south`) so alphabetical port sorting never swaps diagonals
            let i0 = dp
                .ports
                .iter()
                .position(|n| {
                    let nl = n.to_ascii_lowercase();
                    nl.contains("cross") || nl.contains("top") || nl.contains("north")
                })
                .unwrap_or(0);
            let i1 = dp
                .ports
                .iter()
                .position(|n| {
                    let nl = n.to_ascii_lowercase();
                    nl.contains("claw") || nl.contains("bottom") || nl.contains("south")
                })
                .filter(|&i| i != i0)
                .unwrap_or(1 - i0.min(1));
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "c_0_0_f",
                c0_mat[i0][i0],
                &["cross_to_cross", "top_to_top", "north_to_north", "C_0_0_F"],
            );
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "c_1_1_f",
                c0_mat[i1][i1],
                &["claw_to_claw", "bottom_to_bottom", "south_to_south", "C_1_1_F"],
            );
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "c_0_1_mut_f",
                (-c0_mat[i0][i1]).abs(),
                &["cross_to_claw", "top_to_bottom", "C_0_1_mut_F"],
            );
            rec.predictions_si.insert("north_to_south".into(), -(-c0_mat[i0][i1]).abs());
        } else if c0_mat.len() == 1 {
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "c_0_0_f",
                c0_mat[0][0],
                &["cross_to_cross", "top_to_top", "C_0_0_F"],
            );
        }
        rec.extra = json!({
            "ports": dp.ports,
            "c_omega_zero_F": c0_mat,
            "hdm_solves": dp.hdm_solves,
            "rom_dim": dp.rom_dim,
        });
        rec.timings_s.post = t_post0.elapsed().as_secs_f64();
        rec.timings_s.total = t_tot0.elapsed().as_secs_f64();
        return Ok(());
    }

    // --------------------------------------------------------- Branch 2: Pure 2-Port Driven Resonator S21 (`driven_s21` / `driven_filter`)
    if mode == "driven_s21" || mode == "driven" || mode == "driven_filter" {
        let t_asm0 = Instant::now();
        let (f_min, f_max, n_pts, max_hdm) = if mode == "driven_filter" {
            (3.0e9, 16.0e9, 131usize, 32usize)
        } else {
            let (fs, _) = infer_f_search_and_n_modes(layout, false);
            let f_ref = layout
                .references_si
                .get("f_r_hz")
                .or_else(|| layout.references_si.get("f_r_Hz"))
                .copied()
                .unwrap_or(fs / 0.92);
            let k_ref = layout
                .references_si
                .get("kappa_hz")
                .or_else(|| layout.references_si.get("kappa_Hz"))
                .copied()
                .unwrap_or(500.0e3);
            let half_bw = (12.0 * k_ref).clamp(3.0e6, 80.0e6);
            ((f_ref - half_bw).max(1.0e9), f_ref + half_bw, 81usize, 24usize)
        };
        e.call(
            "add_physics",
            &json!({
                "physics": "ems",
                "study": "driven",
                "order": order,
                "study_params": {
                    "f_min": f_min,
                    "f_max": f_max,
                    "n_points": n_pts,
                    "prom_tol": 1.0e-4,
                    "max_hdm_solves": max_hdm,
                }
            }),
        )
        .map_err(|er| er.to_string())?;
        assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, is_w2_or_w3)?;
        for p in &layout.ports {
            add_lumped_port_condition(&mut e, p, 50.0, false)?;
        }
        rec.timings_s.assemble = t_asm0.elapsed().as_secs_f64();

        let t_sol0 = Instant::now();
        e.call(
            "solve_run",
            &json!({
                "physics": "ems",
                "wait": true,
                "isolate": false,
                "save_fields": false,
                "memory_wait_timeout_s": 30.0
            }),
        )
        .map_err(|er| er.to_string())?;
        rec.timings_s.solve = t_sol0.elapsed().as_secs_f64();

        let t_post0 = Instant::now();
        let (res, _) = e.result("ems").ok_or("missing ems driven result")?;
        rec.dofs = res.dofs as u64;
        let dp = DrivenPortData::of(res).ok_or("missing DrivenPortData")?;
        let (recip_err, unit_err) = check_s_matrix_invariants(&dp.s_matrix);
        rec.invariants.max_prom_res = Some(dp.max_residual);
        rec.invariants.reciprocity_err = Some(recip_err);
        rec.invariants.unitarity_err = Some(unit_err);

        if mode == "driven_filter" {
            // Fit heavily damped filter pole from Re[Z_ii(f)] or Re[Y_ii(f)] over [6.0, 12.0] GHz (`Plan §3.4.C`)
            let f_idx = dp
                .ports
                .iter()
                .position(|n| {
                    n.to_ascii_lowercase().contains("dissipator") || n.to_ascii_lowercase().contains("filter")
                })
                .unwrap_or(0);
            let z_re: Vec<f64> = dp.z_matrix.iter().map(|zm| zm[f_idx][f_idx][0]).collect();
            if let Some((f_f, kappa_f, q_f)) = fit_damped_lorentzian_pole(&dp.freqs_hz, &z_re, 5.5e9, 13.5e9) {
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "f_f_hz",
                    f_f,
                    &["f_f_Hz", "f_filter_hz"],
                );
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "kappa_f_hz",
                    kappa_f,
                    &["kappa_f_Hz", "kappa_filter_hz"],
                );
                rec.predictions_si.insert("q_f".into(), q_f);
            }
            let net = e.call("analyse", &json!({ "analysis": "network", "physics": "ems" })).unwrap_or(Value::Null);
            rec.extra = json!({
                "hdm_solves": dp.hdm_solves,
                "rom_dim": dp.rom_dim,
                "bbq_modes": net.get("bbq_modes").cloned().unwrap_or(Value::Null),
            });
        } else if dp.ports.len() >= 2 {
            let net = e
                .call(
                    "analyse",
                    &json!({
                        "analysis": "network",
                        "physics": "ems",
                        "params": {
                            "s21_ports": [dp.ports[0], dp.ports[1]],
                            "resonator_branch": branch,
                            "touchstone_format": "RI"
                        }
                    }),
                )
                .map_err(|er| er.to_string())?;
            if let Some(fit) = net.get("s21_fit") {
                if let Some(fr) = fit["f_r_Hz"].as_f64() {
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "f_r_hz",
                        fr,
                        &["f_r_Hz", "cavity_frequency_hz", "f_r"],
                    );
                }
                if let Some(kp) = fit["kappa_Hz"].as_f64() {
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "kappa_hz",
                        kp,
                        &["kappa_Hz", "kappa"],
                    );
                }
            }
            rec.extra = json!({
                "s21_fit": net.get("s21_fit"),
                "hdm_solves": dp.hdm_solves,
                "rom_dim": dp.rom_dim,
            });
        }
        rec.timings_s.post = t_post0.elapsed().as_secs_f64();
        rec.timings_s.total = t_tot0.elapsed().as_secs_f64();
        return Ok(());
    }

    // --------------------------------------------------------- Branch 3: Unified `QubitCavity` 4-Way Equivalence (`M1..M4`)
    if mode == "qc_all" || mode == "m1" || mode == "m2" || mode == "m3" || mode == "m4" {
        let (f_search, n_modes) = infer_f_search_and_n_modes(layout, true);
        let run_m1 = mode == "qc_all" || mode == "m1";
        let run_m2 = mode == "qc_all" || mode == "m2" || mode == "m3"; // M2 gives exact f_r, kappa to center M3's band!
        let run_m3 = mode == "qc_all" || mode == "m3";
        let run_m4 = mode == "qc_all" || mode == "m4";

        let target_fq = layout
            .references_si
            .get("f01_hz")
            .or_else(|| layout.references_si.get("f01_Hz"))
            .or_else(|| layout.references_si.get("f_q_hz"))
            .copied();
        let target_fr = layout.references_si.get("f_r_hz").or_else(|| layout.references_si.get("f_r_Hz")).copied();

        let mut extra_map = Map::new();
        let (mut m1_vals, mut m2_vals) = ((None, None, None, None, None), (None, None, None, None, None, None));
        let mut c_port_jj_env = None::<f64>;

        // M1: Closed feedline (`pec` on feed ports) + `lumped_jj` -> `craig_bampton`
        if run_m1 {
            let t_asm0 = Instant::now();
            e.call(
                "add_physics",
                &json!({
                    "physics": "ems",
                    "study": "cavity",
                    "order": order,
                    "study_params": { "n_modes": n_modes.max(3), "f_search": f_search }
                }),
            )
            .map_err(|er| er.to_string())?;
            assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, false)?;
            for p in &jj_ports {
                add_jj_condition(&mut e, layout, p)?;
            }
            for p in &feed_ports {
                e.call("add_condition", &json!({ "physics": "ems", "on": p.sheet, "kind": "pec", "params": {} }))
                    .map_err(|er| er.to_string())?;
            }
            rec.timings_s.assemble += t_asm0.elapsed().as_secs_f64();

            let t_s0 = Instant::now();
            e.call(
                "solve_run",
                &json!({
                    "physics": "ems",
                    "wait": true,
                    "isolate": false,
                    "save_fields": false,
                    "memory_wait_timeout_s": 30.0
                }),
            )
            .map_err(|er| er.to_string())?;
            rec.timings_s.solve += t_s0.elapsed().as_secs_f64();

            let t_p0 = Instant::now();
            let (res, _) = e.result("ems").ok_or("missing M1 cavity result")?;
            rec.dofs = res.dofs as u64;
            rec.invariants.equipartition_err = res.scalar("equipartition_error");
            if let Some(rom) = CircuitRom::of(res) {
                if !rom.c_port_f.is_empty() {
                    c_port_jj_env = Some(rom.c_port_f[0][0]);
                }
            }
            let (p_trunc, f_trunc) = cb_truncation_params(jj_ports.len());
            let ham = e
                .call(
                    "analyse",
                    &json!({
                        "analysis": "hamiltonian",
                        "physics": "ems",
                        "params": { "basis": "craig_bampton", "port_trunc": p_trunc, "fock_trunc": f_trunc }
                    }),
                )
                .map_err(|er| er.to_string())?;
            let (f01, alpha, fr, g, chi, _kp, _t1, min_ov) =
                parse_single_qc_from_ham_and_net(&ham, None, target_fq, target_fr, jj_ports.len());
            rec.invariants.min_overlap = min_ov;
            m1_vals = (f01, alpha, fr, g, chi);
            if let Some(v) = f01 {
                rec.predictions_si.insert("m1_f01_hz".into(), v);
            }
            if let Some(v) = alpha {
                rec.predictions_si.insert("m1_alpha_hz".into(), v);
            }
            if let Some(v) = fr {
                rec.predictions_si.insert("m1_f_r_hz".into(), v);
            }
            if let Some(v) = g {
                rec.predictions_si.insert("m1_g_hz".into(), v);
            }
            if let Some(v) = chi {
                rec.predictions_si.insert("m1_chi_hz".into(), v);
            }
            extra_map.insert("M1_hamiltonian".into(), ham);
            rec.timings_s.post += t_p0.elapsed().as_secs_f64();
        }

        // M2: Same identical SimplexMesh, passive `lumped_port` (`Z_0: 50`) + `lumped_jj` -> QEP kappa + CB
        if run_m2 {
            let t_asm0 = Instant::now();
            if e.project().spec("ems").is_some() {
                e.call("remove_physics", &json!({ "physics": "ems" })).map_err(|er| er.to_string())?;
            }
            e.call(
                "add_physics",
                &json!({
                    "physics": "ems",
                    "study": "cavity",
                    "order": order,
                    "study_params": { "n_modes": n_modes.max(3), "f_search": f_search }
                }),
            )
            .map_err(|er| er.to_string())?;
            assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, false)?;
            for p in &jj_ports {
                add_jj_condition(&mut e, layout, p)?;
            }
            for p in &feed_ports {
                add_lumped_port_condition(&mut e, p, 50.0, false)?;
            }
            rec.timings_s.assemble += t_asm0.elapsed().as_secs_f64();

            let t_s0 = Instant::now();
            e.call(
                "solve_run",
                &json!({
                    "physics": "ems",
                    "wait": true,
                    "isolate": false,
                    "save_fields": false,
                    "memory_wait_timeout_s": 30.0
                }),
            )
            .map_err(|er| er.to_string())?;
            rec.timings_s.solve += t_s0.elapsed().as_secs_f64();

            let t_p0 = Instant::now();
            let (res, _) = e.result("ems").ok_or("missing M2 cavity result")?;
            rec.dofs = res.dofs as u64;
            if rec.invariants.equipartition_err.is_none() {
                rec.invariants.equipartition_err = res.scalar("equipartition_error");
            }
            if let Some(rom) = CircuitRom::of(res) {
                if !rom.c_port_f.is_empty() && c_port_jj_env.is_none() {
                    c_port_jj_env = Some(rom.c_port_f[0][0]);
                }
            }
            let (p_trunc, f_trunc) = cb_truncation_params(jj_ports.len());
            let ham = e
                .call(
                    "analyse",
                    &json!({
                        "analysis": "hamiltonian",
                        "physics": "ems",
                        "params": { "basis": "craig_bampton", "port_trunc": p_trunc, "fock_trunc": f_trunc }
                    }),
                )
                .map_err(|er| er.to_string())?;
            let net =
                e.call("analyse", &json!({ "analysis": "network", "physics": "ems" })).map_err(|er| er.to_string())?;
            let (f01, alpha, fr, g, chi, kp, _t1, min_ov) =
                parse_single_qc_from_ham_and_net(&ham, Some(&net), target_fq, target_fr, jj_ports.len());
            if rec.invariants.min_overlap.is_none() {
                rec.invariants.min_overlap = min_ov;
            }
            m2_vals = (f01, alpha, fr, g, chi, kp);
            if let Some(v) = f01 {
                rec.predictions_si.insert("m2_f01_hz".into(), v);
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "f01_hz",
                    v,
                    &["f01_Hz", "f_q_hz", "f01"],
                );
            }
            if let Some(v) = alpha {
                rec.predictions_si.insert("m2_alpha_hz".into(), v);
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "alpha_hz",
                    v,
                    &["alpha_Hz", "alpha"],
                );
            }
            if let Some(v) = fr {
                rec.predictions_si.insert("m2_f_r_hz".into(), v);
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "f_r_hz",
                    v,
                    &["f_r_Hz", "cavity_frequency_hz", "f_r"],
                );
            }
            if let Some(v) = g {
                rec.predictions_si.insert("m2_g_hz".into(), v);
                insert_cqed_observable(&mut rec.predictions_si, &rec.references_si, "g_hz", v, &["g_Hz", "g"]);
            }
            if let Some(v) = chi {
                rec.predictions_si.insert("m2_chi_hz".into(), v);
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "chi_hz",
                    v,
                    &["chi_Hz", "chi", "chi_est_hz"],
                );
            }
            if let Some(v) = kp {
                rec.predictions_si.insert("m2_kappa_hz".into(), v);
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "kappa_hz",
                    v,
                    &["kappa_Hz", "kappa"],
                );
            }
            extra_map.insert("M2_hamiltonian".into(), ham);
            extra_map.insert("M2_network_qep".into(), net.get("qep_modes").cloned().unwrap_or(Value::Null));
            extra_map.insert("M2_network_bbq".into(), net.get("bbq_modes").cloned().unwrap_or(Value::Null));
            rec.timings_s.post += t_p0.elapsed().as_secs_f64();
        }

        // M3: Same identical SimplexMesh, 3-port `ems::driven` + Schur BBQ + Two-Pole Kubo S21^(g/e) punchout chi
        if run_m3 && !jj_ports.is_empty() && feed_ports.len() >= 2 {
            let f_q_est = m2_vals.0.or(m1_vals.0).unwrap_or(4.2e9);
            let f_r_est = m2_vals.2.or(m1_vals.2).unwrap_or(6.2e9);
            let kappa_est = m2_vals.5.unwrap_or(500.0e3);
            let chi_est = m2_vals.4.or(m1_vals.4).unwrap_or(-200.0e3);
            let freqs_m3 = build_m3_segmented_freqs(f_q_est, f_r_est, kappa_est, chi_est);

            let t_asm0 = Instant::now();
            if e.project().spec("ems").is_some() {
                e.call("remove_physics", &json!({ "physics": "ems" })).map_err(|er| er.to_string())?;
            }
            e.call(
                "add_physics",
                &json!({
                    "physics": "ems",
                    "study": "driven",
                    "order": order,
                    "study_params": {
                        "freqs_hz": freqs_m3,
                        "prom_tol": 1.0e-4,
                        "max_hdm_solves": 24,
                        "seed_freqs_hz": [f_q_est, f_r_est],
                    }
                }),
            )
            .map_err(|er| er.to_string())?;
            assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, false)?;
            add_lumped_port_condition(&mut e, feed_ports[0], 50.0, false)?;
            add_lumped_port_condition(&mut e, feed_ports[1], 50.0, false)?;
            add_lumped_port_condition(&mut e, jj_ports[0], 50.0, false)?;
            rec.timings_s.assemble += t_asm0.elapsed().as_secs_f64();

            let t_s0 = Instant::now();
            e.call(
                "solve_run",
                &json!({
                    "physics": "ems",
                    "wait": true,
                    "isolate": false,
                    "save_fields": false,
                    "memory_wait_timeout_s": 30.0
                }),
            )
            .map_err(|er| er.to_string())?;
            rec.timings_s.solve += t_s0.elapsed().as_secs_f64();

            let t_p0 = Instant::now();
            let (res_d, _) = e.result("ems").ok_or("missing M3 driven result")?;
            let dp = DrivenPortData::of(res_d).ok_or("missing M3 DrivenPortData")?;
            rec.invariants.max_prom_res = Some(dp.max_residual);
            let (recip_err, unit_err) = check_s_matrix_invariants(&dp.s_matrix);
            rec.invariants.reciprocity_err = Some(recip_err);
            rec.invariants.unitarity_err = Some(unit_err);

            let jj_name = jj_ports[0].port_name();
            let f0_name = feed_ports[0].port_name();
            let f1_name = feed_ports[1].port_name();
            let jj_idx = dp.ports.iter().position(|n| n == &jj_name).unwrap_or(2);
            let i_f0 = dp.ports.iter().position(|n| n == &f0_name).unwrap_or(0);
            let i_f1 = dp.ports.iter().position(|n| n == &f1_name).unwrap_or(1);

            let p_jj = jj_ports[0];
            let ej_hz = p_jj.e_j_hz.or(layout.e_j_hz).unwrap_or_else(|| {
                let lj = p_jj.l_j_h.or(layout.l_j_h).unwrap_or(10.0e-9);
                nyx_engine::analysis::epr::jj_flux_scale() / lj
            });
            let cj_f = p_jj.c_j_f.or(layout.c_j_f).unwrap_or(0.0);

            // Estimate passive C_env from low-frequency Im[Y_env,JJ(omega)] / omega or CircuitRom::c_port_f[0][0]
            let c_env_lowf = {
                let k0 = 0usize;
                let omega0 = 2.0 * PI * dp.freqs_hz[k0];
                let y3 = &dp.y_matrix[k0];
                // Schur terminate feed ports with 50 Ohm (`Y_load = 0.02 S`)
                let y0 = 1.0 / 50.0;
                let a00 = Cpx::from_arr(y3[i_f0][i_f0]) + Cpx::new(y0, 0.0);
                let a01 = Cpx::from_arr(y3[i_f0][i_f1]);
                let a10 = Cpx::from_arr(y3[i_f1][i_f0]);
                let a11 = Cpx::from_arr(y3[i_f1][i_f1]) + Cpx::new(y0, 0.0);
                let det = (a00 * a11 - a01 * a10).inv();
                let v0 = Cpx::from_arr(y3[i_f0][jj_idx]);
                let v1 = Cpx::from_arr(y3[i_f1][jj_idx]);
                let u0 = Cpx::from_arr(y3[jj_idx][i_f0]);
                let u1 = Cpx::from_arr(y3[jj_idx][i_f1]);
                let w0 = (a11 * v0 - a01 * v1) * det;
                let w1 = (a00 * v1 - a10 * v0) * det;
                let y_env = Cpx::from_arr(y3[jj_idx][jj_idx]) - (u0 * w0 + u1 * w1);
                (y_env.im / omega0).max(1.0e-15)
            };
            let c_env = c_port_jj_env.unwrap_or(c_env_lowf);
            let kubo = TransmonKuboModel::build(c_env, cj_f, ej_hz, 0.0);

            // Slice local ground/excited dip sub-bands around `f_r_est +/- 8*kp` and `(f_r_est + chi_est) +/- 8*kp`
            // so narrow CLT dips (`kappa ~ 10 kHz`) are cleanly resolved for Probst-Taubin S21 circle fitting (`#53` Fix 4).
            let kp_win = kappa_est.max(4.0e3);
            let half_win = (8.0 * kp_win).max(40.0e3);
            let sub_dp_g = slice_driven_window(&dp, f_r_est, half_win);
            let sub_dp_e = slice_driven_window(&dp, f_r_est + chi_est, half_win);

            let dp_g = synthesize_kubo_feedline_2port(&sub_dp_g, [i_f0, i_f1], jj_idx, 50.0, &kubo, false);
            let dp_e = synthesize_kubo_feedline_2port(&sub_dp_e, [i_f0, i_f1], jj_idx, 50.0, &kubo, true);
            if let (Ok((fr_g, kp_g, _qc_g)), Ok((fr_e, _kp_e, _qc_e))) =
                (fit_resonator_s21(&dp_g, branch), fit_resonator_s21(&dp_e, branch))
            {
                let chi_m3 = fr_e - fr_g;
                // Recover g^2 from chi_m3 = 2 g^2 alpha * [1/(Delta*(Delta+alpha)) + 1/(Sigma*(Sigma+alpha))]
                // and compute dressed f01 = f01_bare + g^2 / Delta - g^2 / Sigma (Lamb + Bloch-Siegert)
                let delta = kubo.f01_bare_hz - fr_g;
                let sigma = kubo.f01_bare_hz + fr_g;
                let alpha = kubo.alpha_bare_hz;
                let inv_rwa =
                    if (delta * (delta + alpha)).abs() > 1.0e6 { 1.0 / (delta * (delta + alpha)) } else { 0.0 };
                let inv_cr =
                    if (sigma * (sigma + alpha)).abs() > 1.0e6 { 1.0 / (sigma * (sigma + alpha)) } else { 0.0 };
                let denom_chi = 2.0 * alpha * (inv_rwa + inv_cr);
                let g2 = if denom_chi.abs() > 1.0e-30 { (chi_m3 / denom_chi).max(0.0) } else { 0.0 };
                let lamb_bs = (if delta.abs() > 1.0e6 { g2 / delta } else { 0.0 })
                    - (if sigma.abs() > 1.0e6 { g2 / sigma } else { 0.0 });
                let m3_f01 = kubo.f01_bare_hz + lamb_bs;
                let m3_alpha = kubo.alpha_bare_hz;
                rec.predictions_si.insert("m3_f01_hz".into(), m3_f01);
                rec.predictions_si.insert("m3_alpha_hz".into(), m3_alpha);
                rec.predictions_si.insert("m3_f_r_hz".into(), fr_g);
                rec.predictions_si.insert("m3_f_r_e_hz".into(), fr_e);
                rec.predictions_si.insert("m3_g_hz".into(), g2.sqrt());
                rec.predictions_si.insert("m3_kappa_hz".into(), kp_g);
                rec.predictions_si.insert("m3_chi_hz".into(), chi_m3);
                rec.predictions_si.insert("m3_trk_sum_g".into(), kubo.trk_sum_g);
                rec.predictions_si.insert("m3_trk_sum_e".into(), kubo.trk_sum_e);

                if mode == "m3" {
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "f01_hz",
                        m3_f01,
                        &["f01_Hz", "f01"],
                    );
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "alpha_hz",
                        m3_alpha,
                        &["alpha_Hz", "alpha"],
                    );
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "f_r_hz",
                        fr_g,
                        &["f_r_Hz", "f_r"],
                    );
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "kappa_hz",
                        kp_g,
                        &["kappa_Hz", "kappa"],
                    );
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "chi_hz",
                        chi_m3,
                        &["chi_Hz", "chi"],
                    );
                }
            }
            extra_map.insert(
                "M3_kubo".into(),
                json!({
                    "c_env_F": c_env,
                    "trk_sum_g": kubo.trk_sum_g,
                    "trk_sum_e": kubo.trk_sum_e,
                    "hdm_solves": dp.hdm_solves,
                }),
            );
            rec.timings_s.post += t_p0.elapsed().as_secs_f64();
        }

        // M4: SQuADDS Compositional Synthesis (`Plan §3.3.B` & `Appendix A.4`)
        if run_m4 && !jj_ports.is_empty() {
            let p_jj = jj_ports[0];
            let ej_hz = p_jj.e_j_hz.or(layout.e_j_hz).unwrap_or(15.0e9);
            let cj_f = p_jj.c_j_f.or(layout.c_j_f).unwrap_or(0.0);
            let c_xx = layout.references_si.get("cross_to_cross").copied().or(c_port_jj_env).unwrap_or(70.0e-15);
            let c_xc = layout.references_si.get("cross_to_claw").copied().unwrap_or(5.0e-15);
            let c_cg = layout.references_si.get("claw_to_ground").copied().unwrap_or(20.0e-15);
            let fr_bare = m2_vals.2.or(m1_vals.2).or(target_fr).unwrap_or(6.2e9);
            let is_qw = branch == "hanger";
            let (m4_f01, m4_alpha, m4_g, m4_chi) =
                squadds_compositional_qubit_cavity(c_xx, c_xc, c_cg, cj_f, ej_hz, fr_bare, is_qw, 50.0);
            rec.predictions_si.insert("m4_f01_hz".into(), m4_f01);
            rec.predictions_si.insert("m4_alpha_hz".into(), m4_alpha);
            rec.predictions_si.insert("m4_f_r_hz".into(), fr_bare);
            rec.predictions_si.insert("m4_g_hz".into(), m4_g);
            rec.predictions_si.insert("m4_chi_hz".into(), m4_chi);
        }

        // Compute M1 vs M2 vs M3 equivalence relative deltas (%) for Gate G1 (`Plan §7.1.2`)
        let rel_pct = |a: Option<f64>, b: Option<f64>| -> Option<f64> {
            let (va, vb) = (a?, b?);
            let den = 0.5 * (va.abs() + vb.abs());
            (den > 1.0e-30).then(|| 100.0 * (va - vb).abs() / den)
        };
        if let Some(d) =
            rel_pct(rec.predictions_si.get("m1_f01_hz").copied(), rec.predictions_si.get("m2_f01_hz").copied())
        {
            rec.predictions_si.insert("m1_m2_df01_pct".into(), d);
        }
        if let Some(d) =
            rel_pct(rec.predictions_si.get("m1_f01_hz").copied(), rec.predictions_si.get("m3_f01_hz").copied())
        {
            rec.predictions_si.insert("m1_m3_df01_pct".into(), d);
        }
        if let Some(d) =
            rel_pct(rec.predictions_si.get("m1_f_r_hz").copied(), rec.predictions_si.get("m2_f_r_hz").copied())
        {
            rec.predictions_si.insert("m1_m2_dfr_pct".into(), d);
        }
        if let Some(d) =
            rel_pct(rec.predictions_si.get("m1_f_r_hz").copied(), rec.predictions_si.get("m3_f_r_hz").copied())
        {
            rec.predictions_si.insert("m1_m3_dfr_pct".into(), d);
        }
        if let Some(d) =
            rel_pct(rec.predictions_si.get("m1_alpha_hz").copied(), rec.predictions_si.get("m2_alpha_hz").copied())
        {
            rec.predictions_si.insert("m1_m2_dalpha_pct".into(), d);
        }
        if let Some(d) =
            rel_pct(rec.predictions_si.get("m1_alpha_hz").copied(), rec.predictions_si.get("m3_alpha_hz").copied())
        {
            rec.predictions_si.insert("m1_m3_dalpha_pct".into(), d);
        }
        // Guard `m1_m2_dchi_pct` / `m1_m3_dchi_pct` with `Plan §5.1.2` `25 kHz` unified `chi` denominator floor (`#53` Fix 4)
        let rel_pct_floor = |a: Option<f64>, b: Option<f64>, floor: f64| -> Option<f64> {
            let (va, vb) = (a?, b?);
            let den = (0.5 * (va.abs() + vb.abs())).max(floor);
            (den > 1.0e-30).then(|| 100.0 * (va - vb).abs() / den)
        };
        if let Some(d) = rel_pct_floor(
            rec.predictions_si.get("m1_chi_hz").copied(),
            rec.predictions_si.get("m2_chi_hz").copied(),
            25.0e3,
        ) {
            rec.predictions_si.insert("m1_m2_dchi_pct".into(), d);
        }
        if let Some(d) = rel_pct_floor(
            rec.predictions_si.get("m1_chi_hz").copied(),
            rec.predictions_si.get("m3_chi_hz").copied(),
            25.0e3,
        ) {
            rec.predictions_si.insert("m1_m3_dchi_pct".into(), d);
        }

        rec.extra = Value::Object(extra_map);
        rec.timings_s.total = t_tot0.elapsed().as_secs_f64();
        return Ok(());
    }

    // --------------------------------------------------------- Branch 4: Cavity / Hardware (`S2.1`, `S2.2`, `S3.4`, `WM1`, `MUNINN`, `Dissipator`)
    let (f_search, n_modes) = infer_f_search_and_n_modes(layout, !jj_ports.is_empty());
    let t_asm0 = Instant::now();
    e.call(
        "add_physics",
        &json!({
            "physics": "ems",
            "study": "cavity",
            "order": order,
            "study_params": { "n_modes": n_modes, "f_search": f_search }
        }),
    )
    .map_err(|er| er.to_string())?;
    assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, is_w2_or_w3)?;
    for p in &jj_ports {
        add_jj_condition(&mut e, layout, p)?;
    }
    if mode != "hybrid" {
        for p in &feed_ports {
            add_lumped_port_condition(&mut e, p, 50.0, false)?;
        }
    }
    rec.timings_s.assemble = t_asm0.elapsed().as_secs_f64();

    let t_sol0 = Instant::now();
    e.call(
        "solve_run",
        &json!({
            "physics": "ems",
            "wait": true,
            "isolate": false,
            "save_fields": muninn_single_tls,
            "memory_wait_timeout_s": 30.0
        }),
    )
    .map_err(|er| er.to_string())?;
    rec.timings_s.solve = t_sol0.elapsed().as_secs_f64();

    let t_post0 = Instant::now();
    let (res, _) = e.result("ems").ok_or("missing ems cavity result")?;
    rec.dofs = res.dofs as u64;
    rec.invariants.equipartition_err = res.scalar("equipartition_error");

    if let Some(eq_err) = rec.invariants.equipartition_err {
        if eq_err > 1.0e-8 {
            rec.status = "invariant_fail".into();
            rec.failure_detail = Some(format!("lossless CircuitRom equipartition_error = {eq_err:e} > 1e-8"));
        }
    }

    let mut extra_map = Map::new();

    // Check EPR phase noncollinearity on nonlinear systems (`Plan §5.1.3`)
    if !jj_ports.is_empty() {
        if let Ok(epr_out) = e.call("analyse", &json!({ "analysis": "epr", "physics": "ems" })) {
            if let Some(pnc) = epr_out["phase_noncollinearity_max"].as_f64() {
                rec.invariants.phase_noncollinearity_max = Some(pnc);
            }
        }
    }

    let net_opt = if mode != "hybrid" && (!feed_ports.is_empty() || !jj_ports.is_empty()) {
        let mut net_params = Map::new();
        if feed_ports.len() >= 2 {
            net_params.insert("s21_ports".into(), json!([feed_ports[0].port_name(), feed_ports[1].port_name()]));
            net_params.insert("resonator_branch".into(), json!(branch));
        }
        e.call(
            "analyse",
            &json!({
                "analysis": "network",
                "physics": "ems",
                "params": Value::Object(net_params)
            }),
        )
        .ok()
    } else {
        None
    };

    if !jj_ports.is_empty() {
        // Nonlinear circuit: WM1, MUNINN, Dissipator, or coupled S2.2
        let is_dissipator_squid =
            mode.contains("dissipator") || s_low.contains("s4.3") || jj_ports.iter().any(|p| p.squid_d.is_some());
        let sq_idx_opt = jj_ports.iter().position(|p| {
            p.squid_d.is_some()
                || p.port_name().to_ascii_lowercase().contains("squid")
                || p.port_name().to_ascii_lowercase().contains("dissipator")
        });
        let (p_trunc, f_trunc) = cb_truncation_params(jj_ports.len());
        let mut ham_params = json!({
            "basis": "craig_bampton",
            "port_trunc": p_trunc,
            "fock_trunc": f_trunc,
        });
        if is_dissipator_squid {
            let sq_port_name = jj_ports[sq_idx_opt.unwrap_or(0)].port_name();
            let flux_grid: Vec<f64> = (0..100).map(|k| 0.5 * (k as f64) / 99.0).collect();
            ham_params["sweep_port"] = json!(sq_port_name);
            ham_params["squid_flux_sweep"] = json!(flux_grid);
        }

        let ham = e
            .call("analyse", &json!({ "analysis": "hamiltonian", "physics": "ems", "params": ham_params }))
            .map_err(|er| er.to_string())?;

        let target_fq = layout
            .references_si
            .get("f01_hz")
            .or_else(|| layout.references_si.get("f01_Hz"))
            .or_else(|| layout.references_si.get("f_q_hz"))
            .copied();
        let target_fr = layout.references_si.get("f_r_hz").or_else(|| layout.references_si.get("f_r_Hz")).copied();
        let (f01, alpha, fr, g, chi, kp, t1_p, min_ov) =
            parse_single_qc_from_ham_and_net(&ham, net_opt.as_ref(), target_fq, target_fr, jj_ports.len());
        rec.invariants.min_overlap = min_ov;

        if let Some(v) = f01 {
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "f01_hz",
                v,
                &["f01_Hz", "f_q_hz", "f01"],
            );
        }
        if let Some(v) = alpha {
            insert_cqed_observable(&mut rec.predictions_si, &rec.references_si, "alpha_hz", v, &["alpha_Hz", "alpha"]);
        }
        if let Some(v) = fr {
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "f_r_hz",
                v,
                &["f_r_Hz", "cavity_frequency_hz", "f_r"],
            );
        }
        if let Some(v) = g {
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "g_hz",
                v,
                &["g_Hz", "g", "g_qd_hz", "g_cd_hz"],
            );
        }
        if let Some(v) = chi {
            insert_cqed_observable(
                &mut rec.predictions_si,
                &rec.references_si,
                "chi_hz",
                v,
                &["chi_Hz", "chi", "chi_est_hz"],
            );
        }
        if let Some(v) = kp {
            insert_cqed_observable(&mut rec.predictions_si, &rec.references_si, "kappa_hz", v, &["kappa_Hz", "kappa"]);
        }
        if let Some(v) = t1_p {
            rec.predictions_si.insert("t1_purcell_s".into(), v);
        }

        // Multi-qubit coupled MUNINN chip (`3Q`): record per-qubit f01_i, alpha_i, and ZZ chi_matrix_Hz
        if jj_ports.len() >= 2 {
            if let Some(f01_arr) = ham["f01_Hz"].as_array() {
                for (i, v) in f01_arr.iter().enumerate() {
                    if let Some(f) = v.as_f64() {
                        rec.predictions_si.insert(format!("sub_{i}_f01_hz"), f);
                    }
                }
            }
            if let Some(chi_rows) = ham["chi_matrix_Hz"].as_array() {
                for (i, r_val) in chi_rows.iter().enumerate() {
                    if let Some(cols) = r_val.as_array() {
                        for (j, c_val) in cols.iter().enumerate().skip(i + 1) {
                            if let Some(zz) = c_val.as_f64() {
                                rec.predictions_si.insert(format!("zz_{i}_{j}_hz"), zz);
                            }
                        }
                    }
                }
            }
        }

        // DC-SQUID sweep summary & 4-mode Dissipator couplings (`S4.3`, `Plan §3.4.C`)
        if is_dissipator_squid {
            let sq_i = sq_idx_opt.unwrap_or(0);
            if let Some(fd0) = ham["f01_Hz"].get(sq_i).and_then(Value::as_f64) {
                rec.predictions_si.insert("f_dissipator_0_hz".into(), fd0);
            }
            if let Some(ad0) = ham["alpha_Hz"].get(sq_i).and_then(Value::as_f64) {
                rec.predictions_si.insert("alpha_dissipator_hz".into(), ad0);
            }
            if let (Some(fcavs), Some(g_rows)) = (ham["f_cav_Hz"].as_array(), ham["g_matrix_Hz"].as_array()) {
                let fc_vec: Vec<f64> = fcavs.iter().filter_map(Value::as_f64).collect();
                let t_read = fr.or(target_fr).unwrap_or(5.594e9);
                let c_read = fc_vec
                    .iter()
                    .enumerate()
                    .min_by(|a, b| ((a.1 - t_read).abs()).total_cmp(&((b.1 - t_read).abs())))
                    .map(|(i, _)| i);
                let c_filt = fc_vec
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| Some(*i) != c_read)
                    .min_by(|a, b| ((a.1 - 8.6e9).abs()).total_cmp(&((b.1 - 8.6e9).abs())))
                    .map(|(i, _)| i);
                if let (Some(cr), Some(row_sq)) = (c_read, g_rows.get(sq_i).and_then(Value::as_array)) {
                    if let Some(g_cd) = row_sq.get(cr).and_then(Value::as_f64) {
                        rec.predictions_si.insert("g_cd_hz".into(), g_cd.abs());
                    }
                }
                if let (Some(cf), Some(row_sq)) = (c_filt, g_rows.get(sq_i).and_then(Value::as_array)) {
                    if let Some(g_df) = row_sq.get(cf).and_then(Value::as_f64) {
                        rec.predictions_si.insert("g_df_hz".into(), g_df.abs());
                    }
                }
            }
            // Effective qubit-dissipator coupling `g_qd_hz` from cross-Kerr `chi_matrix_Hz[q_i][sq_i]`
            let q_i = if sq_i == 0 && jj_ports.len() >= 2 { 1 } else { 0 };
            if let (Some(chi_qd), Some(fq_v), Some(fd_v)) = (
                ham["chi_matrix_Hz"].get(q_i).and_then(|r| r.get(sq_i)).and_then(Value::as_f64),
                ham["f01_Hz"].get(q_i).and_then(Value::as_f64),
                ham["f01_Hz"].get(sq_i).and_then(Value::as_f64),
            ) {
                let aq = ham["alpha_Hz"].get(q_i).and_then(Value::as_f64).unwrap_or(-175.0e6);
                let ad = ham["alpha_Hz"].get(sq_i).and_then(Value::as_f64).unwrap_or(-350.0e6);
                let delta_qd = (fd_v - fq_v).abs().max(100.0e6);
                let g_qd_est = (0.5 * chi_qd.abs() * delta_qd * (delta_qd + (aq + ad).abs())
                    / (aq + ad).abs().max(1.0e6))
                .max(0.0)
                .sqrt();
                rec.predictions_si.insert("g_qd_hz".into(), g_qd_est);
            }
        }

        if let Some(sw_arr) = ham.get("sweep").and_then(Value::as_array) {
            let sw_ms = ham.get("sweep_elapsed_s").and_then(Value::as_f64).unwrap_or(0.0) * 1.0e3;
            rec.predictions_si.insert("squid_sweep_ms".into(), sw_ms);
            let sq_i = sq_idx_opt.unwrap_or(0);
            let mut f_sq_max = 0.0_f64;
            let mut f_sq_min = f64::INFINITY;
            for pt in sw_arr {
                if let Some(evals) = pt["f01_Hz"].as_array() {
                    if let Some(fv) = evals.get(sq_i).and_then(Value::as_f64) {
                        f_sq_max = f_sq_max.max(fv);
                        if fv > 1.0e8 {
                            f_sq_min = f_sq_min.min(fv);
                        }
                    } else {
                        for fv in evals.iter().filter_map(Value::as_f64) {
                            f_sq_max = f_sq_max.max(fv);
                            if fv > 1.0e8 {
                                f_sq_min = f_sq_min.min(fv);
                            }
                        }
                    }
                }
            }
            if f_sq_max > 0.0 {
                rec.predictions_si.insert("f_squid_max_hz".into(), f_sq_max);
            }
            if f_sq_min.is_finite() {
                rec.predictions_si.insert("f_squid_min_hz".into(), f_sq_min);
            }
        }

        // S4.2 MUNINN 3D TLS participation + Purcell T_1 (`Plan Appendix A.5`, single-qubit windows)
        if muninn_single_tls {
            let (res_ems, _) = e.result("ems").ok_or("missing ems result for MUNINN TLS")?;
            let sol_idx = CircuitRom::of(res_ems)
                .and_then(|rom| {
                    rom.modes
                        .iter()
                        .enumerate()
                        .max_by(|a, b| {
                            let ua: f64 = a.1.u_ind_jj_j.iter().sum();
                            let ub: f64 = b.1.u_ind_jj_j.iter().sum();
                            ua.total_cmp(&ub)
                        })
                        .map(|(idx, _)| idx)
                })
                .unwrap_or(0);

            let mut tls_params = json!({
                "solution": sol_idx,
                "d": layout.effective_metal_thickness_m().max(150.0e-9),
                "t_ma": 3.0e-9, "eps_ma": 9.8,  "tan_delta_ma": 1.5e-3,
                "t_ms": 2.0e-9, "eps_ms": 11.45, "tan_delta_ms": 5.0e-4,
                "t_sa": 2.0e-9, "eps_sa": 3.9,  "tan_delta_sa": 1.5e-3,
            });
            if layout.stack.trench_depth_m > 0.0 {
                tls_params["trench_depth_m"] = json!(layout.stack.trench_depth_m);
            }
            let tls_out = e
                .call("analyse", &json!({ "analysis": "tls_participation", "physics": "ems", "params": tls_params }))
                .map_err(|er| format!("MUNINN tls_participation: {er}"))?;

            if let Some(q_tls) = tls_out["q_tls"].as_f64().filter(|&q| q.is_finite() && q > 0.0) {
                rec.predictions_si.insert("q_tls".into(), q_tls);
                if let Some(fq) = f01 {
                    let gamma_tls = 2.0 * PI * fq / q_tls;
                    let t1_tls = 1.0 / gamma_tls.max(1.0e-30);
                    let gamma_purcell = t1_p.filter(|&tp| tp > 0.0).map(|tp| 1.0 / tp).unwrap_or(0.0);
                    let t1_tot = 1.0 / (gamma_tls + gamma_purcell).max(1.0e-30);
                    rec.predictions_si.insert("t1_tls_s".into(), t1_tls);
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "t1_s",
                        t1_tot,
                        &["t1_S", "T1_s", "t1"],
                    );
                }
            }
            for k in ["p_ma", "p_ms", "p_sa"] {
                if let Some(v) = tls_out[k].as_f64() {
                    rec.predictions_si.insert(k.into(), v);
                }
            }
            extra_map.insert("tls_participation".into(), tls_out);
        }

        extra_map.insert("hamiltonian".into(), ham);
        if let Some(nv) = net_opt {
            extra_map.insert("network_qep".into(), nv.get("qep_modes").cloned().unwrap_or(Value::Null));
        }
    } else {
        // Linear Cavity Resonator (`S2.1` CLT, `S2.2` NCap, `S3.4` Airbridges)
        let target_fr = layout
            .references_si
            .get("f_r_hz")
            .or_else(|| layout.references_si.get("f_r_Hz"))
            .or_else(|| layout.references_si.get("cavity_frequency_hz"))
            .copied();

        if let Some(ref net) = net_opt {
            if let Some(qep_arr) = net["qep_modes"].as_array() {
                let best_m = if let Some(tfr) = target_fr {
                    qep_arr.iter().min_by(|a, b| {
                        let fa = a["f_re_Hz"].as_f64().unwrap_or(0.0);
                        let fb = b["f_re_Hz"].as_f64().unwrap_or(0.0);
                        (fa - tfr).abs().total_cmp(&(fb - tfr).abs())
                    })
                } else {
                    qep_arr.first()
                };
                if let Some(m) = best_m {
                    if let Some(fr) = m["f_re_Hz"].as_f64() {
                        insert_cqed_observable(
                            &mut rec.predictions_si,
                            &rec.references_si,
                            "f_r_hz",
                            fr,
                            &["f_r_Hz", "cavity_frequency_hz", "f_r"],
                        );
                    }
                    if let Some(kp) = m["kappa_Hz"].as_f64() {
                        insert_cqed_observable(
                            &mut rec.predictions_si,
                            &rec.references_si,
                            "kappa_hz",
                            kp,
                            &["kappa_Hz", "kappa"],
                        );
                    }
                    if let Some(qc) = m["q_c"].as_f64() {
                        rec.predictions_si.insert("q_c".into(), qc);
                    }
                }
            }
            extra_map.insert("qep_modes".into(), net.get("qep_modes").cloned().unwrap_or(Value::Null));
            if let Some(fit) = net.get("s21_fit") {
                extra_map.insert("rom_s21_fit".into(), fit.clone());
            }
        } else {
            // Uncoupled PEC cavity (`S2.2` hybrid `f_rough` or closed airbridge cavity)
            let (res_ems, _) = e.result("ems").ok_or("missing uncoupled ems result")?;
            let f_rough = res_ems
                .solutions
                .iter()
                .filter_map(|s| s.value)
                .min_by(|a, b| {
                    if let Some(tfr) = target_fr {
                        (a - tfr).abs().total_cmp(&(b - tfr).abs())
                    } else {
                        a.total_cmp(b)
                    }
                })
                .ok_or("no converged cavity eigenfrequencies")?;
            rec.predictions_si.insert("f_rough_hz".into(), f_rough);

            if mode == "hybrid" || s_low.contains("s2.2") {
                // Combine uncoupled `f_rough` with coupler `(C_tb, C_tg)` via `squadds_find_kappa_half_wave`
                let (mut c_tb, mut c_tg) = (
                    layout.references_si.get("top_to_bottom").copied().unwrap_or(0.0),
                    layout.references_si.get("top_to_ground").copied().unwrap_or(0.0),
                );
                if let Some(ref c_rel) = layout.coupler_layout_path {
                    let c_full = {
                        let p = Path::new(c_rel);
                        if p.is_absolute() {
                            p.to_path_buf()
                        } else {
                            layout_dir.join(p)
                        }
                    };
                    if let Ok(txt) = std::fs::read_to_string(&c_full) {
                        if let Ok(c_lay) = parse_json_lenient_nan::<LayoutSpec>(&txt) {
                            let c_rec = run_electro_case(
                                &c_lay,
                                c_full.parent().unwrap_or(layout_dir),
                                "p2",
                                "far_field",
                                scale_h,
                            );
                            if c_rec.status == "ok" {
                                c_tb = c_rec.predictions_si.get("top_to_bottom").copied().unwrap_or(c_tb);
                                c_tg = c_rec.predictions_si.get("top_to_ground").copied().unwrap_or(c_tg);
                                extra_map.insert("coupler_electro_dofs".into(), json!(c_rec.dofs));
                            }
                        }
                    }
                }
                if c_tb > 0.0 {
                    let is_tee = layout.resonator_type.as_deref().unwrap_or("").to_ascii_lowercase().contains("tee")
                        || layout.case_id.to_ascii_lowercase().contains("tee");
                    let z_load = if is_tee { 25.0 } else { 50.0 };
                    let (fr_hyb, kp_hyb) = squadds_find_kappa_half_wave(f_rough, c_tb, c_tg, 50.0, z_load);
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "f_r_hz",
                        fr_hyb,
                        &["f_r_Hz", "cavity_frequency_hz", "f_r"],
                    );
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "kappa_hz",
                        kp_hyb,
                        &["kappa_Hz", "kappa"],
                    );
                } else {
                    insert_cqed_observable(
                        &mut rec.predictions_si,
                        &rec.references_si,
                        "f_r_hz",
                        f_rough,
                        &["f_r_Hz", "cavity_frequency_hz", "f_r"],
                    );
                }
            } else {
                insert_cqed_observable(
                    &mut rec.predictions_si,
                    &rec.references_si,
                    "f_r_hz",
                    f_rough,
                    &["f_r_Hz", "cavity_frequency_hz", "f_r"],
                );
            }
        }

        // Optional `qep_and_driven` cross-check on identical mesh (`S2.1` 10-case overlap)
        if mode == "qep_and_driven" && feed_ports.len() >= 2 {
            let fr_qep = rec.predictions_si.get("f_r_hz").copied().unwrap_or(f_search / 0.92);
            let kp_qep = rec.predictions_si.get("kappa_hz").copied().unwrap_or(400.0e3).max(20.0e3);
            let span = (8.0 * kp_qep).clamp(1.0e6, 50.0e6);
            e.call("remove_physics", &json!({ "physics": "ems" })).map_err(|er| er.to_string())?;
            e.call(
                "add_physics",
                &json!({
                    "physics": "ems",
                    "study": "driven",
                    "order": order,
                    "study_params": {
                        "f_min": (fr_qep - span).max(1.0e8),
                        "f_max": fr_qep + span,
                        "n_points": 81,
                        "prom_tol": 1.0e-5,
                        "max_hdm_solves": 16,
                    }
                }),
            )
            .map_err(|er| er.to_string())?;
            assign_em_materials_and_conductors(&mut e, layout, &conductor_sheets, &post_names, is_w2_or_w3)?;
            for p in &feed_ports {
                add_lumped_port_condition(&mut e, p, 50.0, false)?;
            }
            let t_sd = Instant::now();
            e.call(
                "solve_run",
                &json!({
                    "physics": "ems",
                    "wait": true,
                    "isolate": false,
                    "save_fields": false,
                    "memory_wait_timeout_s": 30.0
                }),
            )
            .map_err(|er| er.to_string())?;
            rec.timings_s.solve += t_sd.elapsed().as_secs_f64();

            let net_d = e
                .call(
                    "analyse",
                    &json!({
                        "analysis": "network",
                        "physics": "ems",
                        "params": {
                            "s21_ports": [feed_ports[0].port_name(), feed_ports[1].port_name()],
                            "resonator_branch": branch,
                        }
                    }),
                )
                .map_err(|er| er.to_string())?;
            if let Some(fit) = net_d.get("s21_fit") {
                if let (Some(fr_d), Some(kp_d)) = (fit["f_r_Hz"].as_f64(), fit["kappa_Hz"].as_f64()) {
                    rec.predictions_si.insert("driven_f_r_hz".into(), fr_d);
                    rec.predictions_si.insert("driven_kappa_hz".into(), kp_d);
                    rec.predictions_si
                        .insert("qep_vs_driven_dfr_pct".into(), 100.0 * (fr_d - fr_qep).abs() / fr_qep.max(1.0));
                    rec.predictions_si
                        .insert("qep_vs_driven_dkappa_pct".into(), 100.0 * (kp_d - kp_qep).abs() / kp_qep.max(1.0));
                }
                extra_map.insert("driven_s21_fit".into(), fit.clone());
            }
        }
    }

    rec.extra = Value::Object(extra_map);
    rec.timings_s.post = t_post0.elapsed().as_secs_f64();
    rec.timings_s.total = t_tot0.elapsed().as_secs_f64();
    Ok(())
}
