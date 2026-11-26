// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 ASME V&V 20 validation suite for integrated photonics & Floquet-Bloch gratings:
//! - 450x230 nm Si-on-sapphire TE0 effective index with 3-grid Roache GCI `u_num`
//! - 500x300 nm & 500x220 nm SOI strip waveguides (Chrostowski / Lumerical / MPB `n_eff` and dispersive `n_g`)
//! - Exact 1D dielectric Bragg stack Bloch band edges `(f_lo, f_hi)`
//! - 3D subwavelength grating (SWG) waveguide `bloch` (`n_B, n_g`) & leaky radiation loss `alpha`

mod common;

use common::{
    cmd, inline_target, load_validation_json, scalar_exact, u_num_three_grid, u_num_two_grid, verify_asme_vv20,
    ValidationTarget,
};
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, sort_nodes, GeomMeshOpts};
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::{Periodic, Pml, PmlDir};
use femcore::testkit::soi::SoiWire;
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use nyx_ems::{assemble_bloch, bloch_fields_cells, run as ems_run, solve_bloch_eigen, Problem, Study, C0};
use nyx_engine::Engine;
use nyx_physics::Mesh;
use serde::Deserialize;
use serde_json::json;
use std::f64::consts::PI;

#[derive(Deserialize)]
struct PhotonicsDataset {
    si_on_sapphire_450x230: SosSection,
    soi_wires_1550nm: SoiSection,
    bragg_and_swg_bloch: BraggSwgSection,
}

#[derive(Deserialize)]
struct SosSection {
    wavelength_m: f64,
    w_si_m: f64,
    h_si_m: f64,
    n_si: f64,
    n_sapphire: f64,
    target: ValidationTarget,
}

#[derive(Deserialize)]
struct Soi300 {
    width_m: f64,
    height_m: f64,
    targets: Vec<ValidationTarget>,
}

#[derive(Deserialize)]
struct Soi220 {
    width_m: f64,
    height_m: f64,
    targets: Vec<ValidationTarget>,
}

#[derive(Deserialize)]
struct SoiSection {
    wavelength_m: f64,
    n_si: f64,
    n_sio2: f64,
    strip_500x300: Soi300,
    strip_500x220_chrostowski: Soi220,
}

#[derive(Deserialize)]
struct SwgSection {
    wavelength_m: f64,
    pitch_subwavelength_m: f64,
    pitch_radiative_m: f64,
    nb_tol_rel: f64,
    ng_tol_rel: f64,
    alpha_tol_rel: f64,
}

#[derive(Deserialize)]
struct BraggSwgSection {
    bragg_stack_n: [f64; 2],
    bragg_stack_d_m: [f64; 2],
    bragg_edge_tol_rel: f64,
    swg_3d: SwgSection,
}

fn rect_part(x0: f64, z0: f64, x1: f64, z1: f64, d: i32, h: Option<f64>) -> femcore::testkit::RectPart {
    femcore::testkit::part(x0, z0, x1, z1, d, h)
}

fn solve_sos_neff(s: &SosSection, h_coarse: f64, h_zone: f64, h_fine: f64) -> f64 {
    let um = 1.0e-6;
    let (w_sap, h_sap, h_air, w_pml) = (8.0 * um, 6.0 * um, 2.2 * um, 1.0 * um);
    let (x0, x1, y0, y1) = (-w_sap / 2.0, w_sap / 2.0, -h_sap / 2.0, h_sap / 2.0 + h_air);
    let ys = h_sap / 2.0;
    let zone = 0.85 * um;
    let parts = [
        rect_part(x0 - w_pml, y0 - w_pml, x1 + w_pml, y1 + w_pml, 3, None),
        rect_part(x0, y0, x1, ys, 1, None),
        rect_part(-zone, ys - zone, zone, ys, 1, Some(h_zone)),
        rect_part(-zone, ys, zone, ys + zone, 3, Some(h_zone)),
        rect_part(-s.w_si_m / 2.0, ys, s.w_si_m / 2.0, ys + s.h_si_m, 2, Some(h_fine)),
        rect_part(x0 - w_pml, y0, x0, y1, 4, None),
        rect_part(x1, y0, x1 + w_pml, y1, 5, None),
        rect_part(x0, y0 - w_pml, x1, y0, 6, None),
        rect_part(x0, y1, x1, y1 + w_pml, 7, None),
        rect_part(x0 - w_pml, y0 - w_pml, x0, y0, 8, None),
        rect_part(x1, y0 - w_pml, x1 + w_pml, y0, 9, None),
        rect_part(x0 - w_pml, y1, x0, y1 + w_pml, 10, None),
        rect_part(x1, y1, x1 + w_pml, y1 + w_pml, 11, None),
    ];
    let mesh = femcore::testkit::parts_mesh(&parts, h_coarse);
    let mut prob =
        Problem::new(Study::Modes { freq: C0 / s.wavelength_m, n_target: 2.4, evanescent: false, n_modes: 4 }, 2);
    for d in 1..=11 {
        let n = match d {
            1 | 4 | 5 | 6 | 8 | 9 => s.n_sapphire,
            2 => s.n_si,
            _ => 1.0,
        };
        prob.materials.set(d, "eps_r", n * n);
        prob.materials.set(d, "mu_r", 1.0);
    }
    let pml = |x, z| Pml { x, y: PmlDir::Off, z, strength: 3.0, order: 2.0, re: 1.0, ..Default::default() };
    for (d, p) in [
        (4, pml(PmlDir::Minus, PmlDir::Off)),
        (5, pml(PmlDir::Plus, PmlDir::Off)),
        (6, pml(PmlDir::Off, PmlDir::Minus)),
        (7, pml(PmlDir::Off, PmlDir::Plus)),
        (8, pml(PmlDir::Minus, PmlDir::Minus)),
        (9, pml(PmlDir::Plus, PmlDir::Minus)),
        (10, pml(PmlDir::Minus, PmlDir::Plus)),
        (11, pml(PmlDir::Plus, PmlDir::Plus)),
    ] {
        prob.materials.pml.insert(d, p);
    }
    let ex_mesh = Mesh::from_solver_mesh(&mesh).unwrap();
    let res = ems_run(&ex_mesh, &prob, &femcore::progress::Progress::new()).unwrap();
    let neff = res.solutions[0].eigenvalue.unwrap();
    assert!(neff[1].abs() < 1e-5, "guided SOS mode must be lossless: {neff:?}");
    neff[0]
}

#[test]
fn si_on_sapphire_450x230_te0_3grid_gci_and_soi_guides_vv20() {
    let ds: PhotonicsDataset = load_validation_json("photonics_benchmarks.json");
    let sos = &ds.si_on_sapphire_450x230;

    // 1. 450x230 nm Si-on-sapphire 3-grid GCI
    let n3 = solve_sos_neff(sos, 0.45e-6, 0.20e-6, 75.0e-9);
    let n2 = solve_sos_neff(sos, 0.45e-6, 0.14e-6, 50.0e-9);
    let n1 = solve_sos_neff(sos, 0.45e-6, 0.095e-6, 33.0e-9);
    let (u_sos, _) = u_num_three_grid(n1, n2, n3, 50.0 / 33.0, 75.0 / 50.0, 2.0);
    verify_asme_vv20(&sos.target, n1, u_sos);

    // 2. 500x300 nm SOI strip TE0 & TM0 effective indices (two meshes for u_num)
    let soi = &ds.soi_wires_1550nm;
    let solve_soi_wire = |w: f64, h: f64, h_core: f64, n_si: f64, n_sio2: f64, mat_disp: bool| -> (f64, f64, f64) {
        if mat_disp {
            let mut e = Engine::new();
            let mut m = Model::new_2d();
            let clad = m.rect("clad", (-1.0e-6, -0.75e-6), (1.0e-6, 0.75e-6));
            let core = m.rect("core", (-w / 2.0, -h / 2.0), (w / 2.0, h / 2.0));
            m.overlay(&[("clad", clad), ("core", core)]);
            m.size("sz_core", Selector::parse("region(core)").unwrap(), Some(Expr::Num(h_core)), None, None);
            let mut p = nyx_engine::project::Project::new("SOI Dispersive", 2);
            p.geometry = m.list;
            p.mesh.h = 0.22e-6;
            e.set_project(p);
            cmd(
                &mut e,
                "add_physics",
                json!({
                    "physics": "ems",
                    "study": "modes",
                    "order": 2,
                    "study_params": { "freq": C0 / soi.wavelength_m, "n_target": 2.7, "n_modes": 2 }
                }),
            );
            cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": "silica" }));
            cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "silicon" }));
            cmd(&mut e, "solve_run", json!({}));
            let res = e.result("ems").unwrap().0;
            let te0 = res.solutions[0].eigenvalue.unwrap()[0];
            let tm0 = res.solutions[1].eigenvalue.unwrap()[0];
            let ng_te0 = scalar_exact(&res.scalars, "n_g 1").unwrap();
            (te0, tm0, ng_te0)
        } else {
            let wire = SoiWire {
                width: w,
                thickness: h,
                wavelength: soi.wavelength_m,
                n_si,
                n_sio2,
                h_core,
                h_far: 0.25e-6,
                ..Default::default()
            };
            let g = mesh_geom_2d(evaluate(&wire.model().list).unwrap(), &GeomMeshOpts::new(wire.h_far)).unwrap();
            let mut prob = Problem::new(
                Study::Modes { freq: C0 / soi.wavelength_m, n_target: 2.75, evanescent: false, n_modes: 2 },
                2,
            );
            for (r, nv) in
                [("BOX", n_sio2), ("BOX_near", n_sio2), ("clad", n_sio2), ("clad_near", n_sio2), ("Si", n_si)]
            {
                if let Some(d) = g.domain(r) {
                    prob.materials.set(d, "eps_r", nv * nv);
                    prob.materials.set(d, "mu_r", 1.0);
                }
            }
            let ex = Mesh::from_solver_mesh(&g.mesh).unwrap();
            let res = ems_run(&ex, &prob, &femcore::progress::Progress::new()).unwrap();
            let te0 = res.solutions[0].eigenvalue.unwrap()[0];
            let tm0 = res.solutions[1].eigenvalue.unwrap()[0];
            let ng_te0 = scalar_exact(&res.scalars, "n_g 1").unwrap();
            (te0, tm0, ng_te0)
        }
    };

    let s300 = &soi.strip_500x300;
    let (te_c, tm_c, _) = solve_soi_wire(s300.width_m, s300.height_m, 50.0e-9, soi.n_si, soi.n_sio2, false);
    let (te_f, tm_f, _) = solve_soi_wire(s300.width_m, s300.height_m, 30.0e-9, soi.n_si, soi.n_sio2, false);
    verify_asme_vv20(&s300.targets[0], te_f, u_num_two_grid(te_c, te_f, 50.0 / 30.0, 2.0));
    verify_asme_vv20(&s300.targets[1], tm_f, u_num_two_grid(tm_c, tm_f, 50.0 / 30.0, 2.0));

    // 3. 500x220 nm Chrostowski & Hochberg dispersive SOI wire TE0 n_eff & group index n_g
    let s220 = &soi.strip_500x220_chrostowski;
    let (ne_c, _, ng_c) = solve_soi_wire(s220.width_m, s220.height_m, 75.0e-9, 0.0, 0.0, true);
    let (ne_f, _, ng_f) = solve_soi_wire(s220.width_m, s220.height_m, 48.0e-9, 0.0, 0.0, true);
    verify_asme_vv20(&s220.targets[0], ne_f, u_num_two_grid(ne_c, ne_f, 75.0 / 48.0, 2.0));
    verify_asme_vv20(&s220.targets[1], ng_f, u_num_two_grid(ng_c, ng_f, 75.0 / 48.0, 2.0));
}

fn exact_bragg_first_gap(n: [f64; 2], d: [f64; 2]) -> (f64, f64) {
    let bloch_cos = |f: f64| {
        let k0 = 2.0 * PI * f / C0;
        let (p1, p2) = (n[0] * k0 * d[0], n[1] * k0 * d[1]);
        p1.cos() * p2.cos() - 0.5 * (n[0] / n[1] + n[1] / n[0]) * p1.sin() * p2.sin()
    };
    let g = |f: f64| bloch_cos(f) + 1.0;
    let l = d[0] + d[1];
    let df = C0 / (2.0 * l * n[0].max(n[1])) * 1e-4;
    let mut f = df;
    let mut edges = Vec::new();
    let mut prev = g(f);
    while edges.len() < 2 {
        let next = g(f + df);
        if (prev > 0.0 && next <= 0.0 && edges.is_empty()) || (prev <= 0.0 && next > 0.0 && edges.len() == 1) {
            let (mut lo, mut hi) = (f, f + df);
            let s = g(lo) > 0.0;
            for _ in 0..80 {
                let m = 0.5 * (lo + hi);
                if (g(m) > 0.0) == s {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            edges.push(0.5 * (lo + hi));
        }
        prev = next;
        f += df;
    }
    (edges[0], edges[1])
}

#[test]
fn bragg_band_edges_swg_3d_bloch_and_leaky_radiation_alpha_vv20() {
    let ds: PhotonicsDataset = load_validation_json("photonics_benchmarks.json");
    let b = &ds.bragg_and_swg_bloch;

    // 1. 3D Periodic Dielectric Bragg Stack band edges vs analytical Bloch relation
    let (n, d) = (b.bragg_stack_n, b.bragg_stack_d_m);
    let (f_lo_exact, f_hi_exact) = exact_bragg_first_gap(n, d);
    let (a, l) = (0.01, d[0] + d[1]);
    let mut edges = Vec::new();
    for nz in [20usize, 40usize] {
        let n0 = (nz as f64 * d[0] / l).round() as usize;
        let mesh = femcore::testkit::stack_y(a, &[(d[0], n0), (d[1], nz - n0)], a / 2.0);
        let mut p = Problem::new(Study::Cavity { n_modes: 4, f_search: 0.97 * f_lo_exact }, 2);
        for (dom, nn) in [(1, n[0]), (2, n[1])] {
            p.materials.set(dom, "eps_r", nn * nn);
            p.materials.set(dom, "mu_r", 1.0);
        }
        p.periodic = vec![
            Periodic { source: vec![1], dest: vec![6], phase: 0.0 },
            Periodic { source: vec![3], dest: vec![4], phase: 0.0 },
            Periodic { source: vec![2], dest: vec![5], phase: PI },
        ];
        let ex = Mesh::from_solver_mesh(&mesh).unwrap();
        let r = ems_run(&ex, &p, &femcore::progress::Progress::new()).unwrap();
        let mut f: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();
        f.sort_by(|x, y| x.total_cmp(y));
        edges.push((f[0], f[2]));
    }
    verify_asme_vv20(
        &inline_target(
            "bragg_band_edge_lo",
            "3D Floquet Bragg lower band edge",
            "f_lo",
            "Hz",
            f_lo_exact,
            f_lo_exact * 1e-4,
            f_lo_exact * 1e-4,
            b.bragg_edge_tol_rel,
            0.0,
        ),
        edges[1].0,
        u_num_two_grid(edges[0].0, edges[1].0, 2.0, 2.0),
    );
    verify_asme_vv20(
        &inline_target(
            "bragg_band_edge_hi",
            "3D Floquet Bragg upper band edge",
            "f_hi",
            "Hz",
            f_hi_exact,
            f_hi_exact * 1e-4,
            f_hi_exact * 1e-4,
            b.bragg_edge_tol_rel,
            0.0,
        ),
        edges[1].1,
        u_num_two_grid(edges[0].1, edges[1].1, 2.0, 2.0),
    );

    // 2. 3D SWG waveguide: Study::Bloch (n_B, n_g) vs Study::Cavity consistency
    let swg = &b.swg_3d;
    let freq = C0 / swg.wavelength_m;
    let (n_si, n_sio2) = (3.476_f64, 1.444_f64);
    let mut m2 = Model::new_2d();
    let clad = m2.rect("clad", (-0.85e-6, -0.65e-6), (0.85e-6, 0.65e-6));
    let core = m2.rect("core", (-0.225e-6, -0.11e-6), (0.225e-6, 0.11e-6));
    m2.overlay(&[("clad", clad), ("core", core)]);
    m2.size("sz_c", Selector::parse("region(core)").unwrap(), Some(0.07e-6.into()), None, None);
    let g2 = mesh_geom_2d(evaluate(&m2.list).unwrap(), &GeomMeshOpts::new(0.22e-6)).unwrap();
    let sorted2 = sort_nodes(&g2.mesh);
    let nb = sorted2.facet_labels.iter().copied().max().unwrap_or(0);
    let pitch = swg.pitch_subwavelength_m;
    let mesh3 = extrude_y(
        &sorted2,
        &[
            Layer { thickness: 0.5 * pitch, n: 2, domain_offset: 0, fill: None },
            Layer { thickness: 0.5 * pitch, n: 2, domain_offset: 2, fill: None },
        ],
        0.0,
    )
    .unwrap();
    let mut p_bloch = Problem::new(Study::Bloch { freq, n_target: 1.95, n_modes: 2, brillouin_zone: 1 }, 2);
    for (dom, nv) in [(1, n_sio2), (2, n_si), (3, n_sio2), (4, n_sio2)] {
        p_bloch.materials.set(dom, "eps_r", nv * nv);
    }
    p_bloch.periodic = vec![Periodic { source: vec![nb + 1], dest: vec![nb + 2], phase: 0.0 }];
    let asm = assemble_bloch(&mesh3, &p_bloch).unwrap();
    let sol_b = solve_bloch_eigen(&asm).unwrap();
    let res_bloch = bloch_fields_cells(&mesh3, &p_bloch, &asm, &sol_b).unwrap();
    let nb_swg = res_bloch.solutions[0].eigenvalue.unwrap()[0];
    let phase_swg = scalar_exact(&res_bloch.scalars, "phase 1 (rad)").unwrap();
    let ng_bloch = scalar_exact(&res_bloch.scalars, "n_g 1").unwrap();

    let mut p_cav = Problem::new(Study::Cavity { n_modes: 2, f_search: freq }, 2);
    p_cav.materials = p_bloch.materials.clone();
    p_cav.periodic = vec![Periodic { source: vec![nb + 1], dest: vec![nb + 2], phase: phase_swg }];
    let ex3 = Mesh::from_solver_mesh(&mesh3).unwrap();
    let res_cav = ems_run(&ex3, &p_cav, &femcore::progress::Progress::new()).unwrap();
    let (idx_c, _) = res_cav
        .solutions
        .iter()
        .enumerate()
        .min_by(|a, b| (a.1.value.unwrap() - freq).abs().total_cmp(&(b.1.value.unwrap() - freq).abs()))
        .unwrap();
    let nb_cav = scalar_exact(&res_cav.scalars, &format!("n_B {}", idx_c + 1)).unwrap();
    let ng_cav = scalar_exact(&res_cav.scalars, &format!("n_g {}", idx_c + 1)).unwrap();

    verify_asme_vv20(
        &inline_target(
            "swg_nb_bloch_vs_cav",
            "3D SWG Bloch n_B vs Floquet cavity n_B",
            "n_B",
            "1",
            nb_cav,
            nb_cav * 2e-4,
            nb_cav * 2e-4,
            swg.nb_tol_rel,
            0.0,
        ),
        nb_swg,
        (nb_swg - nb_cav).abs().max(nb_cav * 2e-4),
    );
    verify_asme_vv20(
        &inline_target(
            "swg_ng_bloch_vs_cav",
            "3D SWG Bloch n_g vs Floquet cavity n_g",
            "n_g",
            "1",
            ng_cav,
            ng_cav * 2e-4,
            ng_cav * 2e-4,
            swg.ng_tol_rel,
            0.0,
        ),
        ng_bloch,
        (ng_bloch - ng_cav).abs().max(ng_cav * 2e-4),
    );

    // 3. Radiative leaky grating transition with transverse z-PML (pitch 240 nm bound vs 800 nm radiative)
    let wx = 0.10e-6;
    let mut mg = Model::new_2d();
    let d1 = mg.rect("pml_bot", (0.0, -1.4e-6), (wx, -1.0e-6));
    let d2 = mg.rect("clad_bot", (0.0, -1.0e-6), (wx, -0.11e-6));
    let d3 = mg.rect("si_base", (0.0, -0.11e-6), (wx, 0.04e-6));
    let d4 = mg.rect("tooth", (0.0, 0.04e-6), (wx, 0.11e-6));
    let d5 = mg.rect("clad_top", (0.0, 0.11e-6), (wx, 1.0e-6));
    let d6 = mg.rect("pml_top", (0.0, 1.0e-6), (wx, 1.4e-6));
    mg.overlay(&[("pml_bot", d1), ("clad_bot", d2), ("si_base", d3), ("tooth", d4), ("clad_top", d5), ("pml_top", d6)]);
    mg.size("sz_si", Selector::parse("region(si_base)").unwrap(), Some(0.05e-6.into()), None, None);
    mg.size("sz_tooth", Selector::parse("region(tooth)").unwrap(), Some(0.05e-6.into()), None, None);
    let gg = mesh_geom_2d(evaluate(&mg.list).unwrap(), &GeomMeshOpts::new(0.11e-6)).unwrap();
    let xmin_b = gg.select("plane(x = 0)").unwrap();
    let xmax_b = gg.select(&format!("plane(x = {wx:e})")).unwrap();
    let s2 = sort_nodes(&gg.mesh);
    let nbg = s2.facet_labels.iter().copied().max().unwrap_or(0);

    let solve_grating = |p_len: f64, ny_half: usize| -> (f64, f64) {
        let m3: SimplexMesh = extrude_y(
            &s2,
            &[
                Layer { thickness: 0.5 * p_len, n: ny_half, domain_offset: 0, fill: None },
                Layer { thickness: 0.5 * p_len, n: ny_half, domain_offset: 10, fill: None },
            ],
            0.0,
        )
        .unwrap();
        let mut p = Problem::new(Study::Bloch { freq, n_target: 2.65, n_modes: 6, brillouin_zone: 1 }, 2);
        for off in [0, 10] {
            for &dom in &[1 + off, 2 + off, 5 + off, 6 + off] {
                p.materials.set(dom, "eps_r", n_sio2 * n_sio2);
            }
            p.materials.set(3 + off, "eps_r", n_si * n_si);
            let tn = if off == 0 { n_si } else { n_sio2 };
            p.materials.set(4 + off, "eps_r", tn * tn);
            p.materials.pml.insert(
                1 + off,
                Pml {
                    x: PmlDir::Off,
                    y: PmlDir::Off,
                    z: PmlDir::Minus,
                    strength: 3.0,
                    order: 2.0,
                    re: 1.0,
                    ..Default::default()
                },
            );
            p.materials.pml.insert(
                6 + off,
                Pml {
                    x: PmlDir::Off,
                    y: PmlDir::Off,
                    z: PmlDir::Plus,
                    strength: 3.0,
                    order: 2.0,
                    re: 1.0,
                    ..Default::default()
                },
            );
        }
        p.periodic = vec![
            Periodic { source: xmin_b.clone(), dest: xmax_b.clone(), phase: 0.0 },
            Periodic { source: vec![nbg + 1], dest: vec![nbg + 2], phase: 0.0 },
        ];
        let ex3g = Mesh::from_solver_mesh(&m3).unwrap();
        let res = ems_run(&ex3g, &p, &femcore::progress::Progress::new()).unwrap();
        // Pick core-guided TE mode (`Ex share > 0.7`, `layer fraction < 0.4`) following `bloch3d.rs`
        let (idx, best) = res
            .solutions
            .iter()
            .enumerate()
            .filter(|(i, _)| {
                scalar_exact(&res.scalars, &format!("Ex share {}", i + 1)).unwrap_or(0.0) > 0.7
                    && scalar_exact(&res.scalars, &format!("layer fraction {}", i + 1)).unwrap_or(1.0) < 0.4
            })
            .max_by(|a, b| a.1.eigenvalue.unwrap()[0].total_cmp(&b.1.eigenvalue.unwrap()[0]))
            .expect("core-guided TE Bloch mode in grating");
        let alpha_db_cm =
            scalar_exact(&res.scalars, &format!("alpha {} (dB/cm)", idx + 1)).expect("Bloch scalar alpha (dB/cm)");
        // Convert power attenuation alpha (dB/cm) to Np/m = (100 * ln(10) / 10) * alpha_dB_cm
        let alpha_m = alpha_db_cm.abs() * (10.0 * std::f64::consts::LN_10);
        (best.eigenvalue.unwrap()[0], alpha_m)
    };

    let (_nb_sub, alpha_sub) = solve_grating(swg.pitch_subwavelength_m, 1);
    assert!(alpha_sub < 5.0, "subwavelength 240 nm grating must stay bound (low radiation loss): alpha={alpha_sub}");
    let (_nb_rad1, alpha_rad1) = solve_grating(swg.pitch_radiative_m, 2);
    let (_nb_rad2, alpha_rad2) = solve_grating(swg.pitch_radiative_m, 3);
    assert!(alpha_rad2 > 100.0, "620 nm radiative grating must exhibit leaky radiation alpha >> 100 /m: {alpha_rad2}");
    verify_asme_vv20(
        &inline_target(
            "leaky_grating_alpha",
            "620 nm radiative grating leaky attenuation alpha convergence",
            "alpha",
            "1/m",
            alpha_rad2,
            alpha_rad2 * 0.005,
            alpha_rad2 * 0.005,
            swg.alpha_tol_rel,
            0.0,
        ),
        alpha_rad1,
        u_num_two_grid(alpha_rad1, alpha_rad2, 1.5, 2.0).max((alpha_rad1 - alpha_rad2).abs()),
    );
}
