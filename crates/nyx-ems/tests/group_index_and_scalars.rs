// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for FR-3 (single-solve dispersive Poynting group index `v_g`/`n_g`, 3D mode scalars,
//! all-PMC / Gamma-point null-space pin, `pmc` across all studies, and midplane `symmetry` BCs):
//! 1. Single-solve `n_g` vs 5-point central finite difference (`< 5e-4` rel) on dispersive `silicon` + `silica`
//!    for 2D `modes`, 3D `bloch`, and 3D `cavity`.
//! 2. Component shares (`Ex` + `Ey` + `Ez` = 1) and domain shares (`sum` = 1) + all-PMC & `phase = 0` cavity null-space pin.
//! 3. `pmc` BC and `symmetry` (`pec_x`, `pmc_x`, `pec_x_pmc_z`, `pmc_x_pec_z`) on 1/2 and 1/4 meshes matching
//!    full-mesh `TE0` and `TM0` (`< 1e-4` rel) in 2D `modes`, 3D `bloch`, and 3D `cavity`.

mod common;

use femcore::material::Library;
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, sort_nodes, GeomMeshOpts};
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::Periodic;
use geom::model::Model;
use geom::select::Selector;
use nyx_ems::{Ems, Problem, Study, Symmetry, C0};
use nyx_physics::{BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, Scalar, SolveResult};
use serde_json::json;
use std::f64::consts::PI;
use std::path::Path;
use std::sync::Arc;

fn scalar_num(r: &SolveResult, name: &str) -> Option<f64> {
    r.scalars.iter().find_map(|(k, v)| match (k == name, v) {
        (true, Scalar::Num(x)) => Some(*x),
        _ => None,
    })
}

fn load_repo_library() -> Arc<Library> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Arc::new(Library::load_dir(&dir).expect("load repo materials library"))
}

struct WireResolver {
    y0: i32,
    y1: i32,
}

impl LabelResolver for WireResolver {
    fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
        match r {
            "clad" => Ok(vec![1]),
            "core" => Ok(vec![2]),
            _ => Err(Error::invalid(format!("no domain {r}"))),
        }
    }
    fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
        match s {
            "ymin" => Ok(vec![self.y0]),
            "ymax" => Ok(vec![self.y1]),
            _ => Err(Error::invalid(format!("no boundary {s}"))),
        }
    }
}

#[test]
fn dispersive_poynting_group_index_matches_5pt_fd_in_2d_modes_and_3d_bloch_and_cavity() {
    let lib = load_repo_library();

    // Build SOI wire cross-section: 450 x 220 nm silicon core in 1.7 x 1.3 um silica cladding
    let mut m = Model::new_2d();
    let clad = m.rect("clad", (-0.85e-6, -0.65e-6), (0.85e-6, 0.65e-6));
    let core = m.rect("core", (-0.225e-6, -0.11e-6), (0.225e-6, 0.11e-6));
    m.overlay(&[("clad", clad), ("core", core)]);
    m.size("sz_core", Selector::parse("region(core)").unwrap(), Some(0.06e-6.into()), None, None);
    let g2 = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.20e-6)).unwrap();
    let mesh2 = sort_nodes(&g2.mesh);
    let ex2 = Mesh::from_solver_mesh(&mesh2).unwrap();

    let pitch = 0.20e-6;
    let nb = mesh2.facet_labels.iter().copied().max().unwrap_or(0);
    let mesh3 = extrude_y(&mesh2, &[Layer { thickness: pitch, n: 2, domain_offset: 0, fill: None }], 0.0).unwrap();
    let ex3 = Mesh::from_solver_mesh(&mesh3).unwrap();
    let resolver = WireResolver { y0: nb + 1, y1: nb + 2 };

    // Evaluate at f0 = 196.0 THz (lambda0 = 1.52955 um, strictly inside silicon.toml's smooth (1.50, 1.55) um PCHIP cubic
    // interval so the 5-point FD stencil does not straddle a C^1 table knot!)
    let f0 = 196.0e12;
    let df = 0.15e12;

    // 1. 2D `modes` dispersive single-solve n_g vs 5-point central FD of beta(f) = (2 pi f / c0) * n_eff(f)
    let solve_2d_at = |f: f64| -> SolveResult {
        let mut spec = ProblemSpec::new("ems", 2, "modes");
        spec.study.params = json!({ "freq": f, "n_target": 2.5, "n_modes": 2 });
        spec.materials.insert("core".into(), MaterialAssign::Named("silicon".into()));
        spec.materials.insert("clad".into(), MaterialAssign::Named("silica".into()));
        let mut ctx = Ctx::new(&resolver);
        ctx.library = Some(lib.clone());
        Ems.solve(&ex2, &spec, &mut ctx).unwrap()
    };

    let r2_m2 = solve_2d_at(f0 - 2.0 * df);
    let r2_m1 = solve_2d_at(f0 - df);
    let r2_0 = solve_2d_at(f0);
    let r2_p1 = solve_2d_at(f0 + df);
    let r2_p2 = solve_2d_at(f0 + 2.0 * df);

    for mode_idx in 0..2 {
        let beta_of = |r: &SolveResult, f: f64| 2.0 * PI * f / C0 * r.solutions[mode_idx].eigenvalue.unwrap()[0];
        let dbeta_df = (beta_of(&r2_m2, f0 - 2.0 * df) - 8.0 * beta_of(&r2_m1, f0 - df)
            + 8.0 * beta_of(&r2_p1, f0 + df)
            - beta_of(&r2_p2, f0 + 2.0 * df))
            / (12.0 * df);
        let ng_fd = (C0 / (2.0 * PI)) * dbeta_df;
        let ng_single = scalar_num(&r2_0, &format!("n_g {}", mode_idx + 1)).unwrap();
        let rel = ((ng_single - ng_fd) / ng_fd).abs();
        assert!(
            rel < 5e-4,
            "2D dispersive mode {}: single-solve n_g = {ng_single:.8} vs 5-pt FD {ng_fd:.8} (rel = {rel:.2e})",
            mode_idx + 1
        );
    }

    // 2. 3D `bloch` dispersive single-solve n_g vs 5-point central FD of k_y(f)
    let solve_3d_bloch_at = |f: f64| -> SolveResult {
        let mut spec = ProblemSpec::new("ems", 2, "bloch");
        spec.study.params = json!({ "freq": f, "n_target": 2.5, "n_modes": 2 });
        spec.materials.insert("core".into(), MaterialAssign::Named("silicon".into()));
        spec.materials.insert("clad".into(), MaterialAssign::Named("silica".into()));
        spec.bcs.push(BcSpec { on: "ymin".into(), kind: "periodic".into(), params: json!({ "dest": "ymax" }) });
        let mut ctx = Ctx::new(&resolver);
        ctx.library = Some(lib.clone());
        Ems.solve(&ex3, &spec, &mut ctx).unwrap()
    };

    let b_m2 = solve_3d_bloch_at(f0 - 2.0 * df);
    let b_m1 = solve_3d_bloch_at(f0 - df);
    let b_0 = solve_3d_bloch_at(f0);
    let b_p1 = solve_3d_bloch_at(f0 + df);
    let b_p2 = solve_3d_bloch_at(f0 + 2.0 * df);

    let ky_of = |r: &SolveResult| scalar_num(r, "k_y re 1 (rad/m)").unwrap();
    let dky_df = (ky_of(&b_m2) - 8.0 * ky_of(&b_m1) + 8.0 * ky_of(&b_p1) - ky_of(&b_p2)) / (12.0 * df);
    let ng_bloch_fd = (C0 / (2.0 * PI)) * dky_df;
    let ng_bloch = scalar_num(&b_0, "n_g 1").unwrap();
    let rel_bloch = ((ng_bloch - ng_bloch_fd) / ng_bloch_fd).abs();
    assert!(
        rel_bloch < 5e-4,
        "3D Bloch dispersive TE0: single-solve n_g = {ng_bloch:.8} vs 5-pt FD {ng_bloch_fd:.8} (rel = {rel_bloch:.2e})"
    );

    // 3. 3D `cavity` dispersive single-solve n_g at self-consistent (f_search = f0, phase = k_y(f0) * pitch)
    let phase_0 = scalar_num(&b_0, "phase 1 (rad)").unwrap();
    let mut spec_cav = ProblemSpec::new("ems", 2, "cavity");
    spec_cav.study.params = json!({ "f_search": f0, "n_modes": 2 });
    spec_cav.materials.insert("core".into(), MaterialAssign::Named("silicon".into()));
    spec_cav.materials.insert("clad".into(), MaterialAssign::Named("silica".into()));
    spec_cav.bcs.push(BcSpec {
        on: "ymin".into(),
        kind: "periodic".into(),
        params: json!({ "dest": "ymax", "phase": phase_0 }),
    });
    let mut ctx_cav = Ctx::new(&resolver);
    ctx_cav.library = Some(lib);
    let r_cav = Ems.solve(&ex3, &spec_cav, &mut ctx_cav).unwrap();
    let idx_cav = r_cav
        .solutions
        .iter()
        .enumerate()
        .min_by(|a, b| (a.1.value.unwrap() - f0).abs().total_cmp(&(b.1.value.unwrap() - f0).abs()))
        .unwrap()
        .0;
    let ng_cav = scalar_num(&r_cav, &format!("n_g {}", idx_cav + 1)).unwrap();
    let rel_cav = ((ng_cav - ng_bloch_fd) / ng_bloch_fd).abs();
    assert!(
        rel_cav < 5e-4,
        "3D Cavity dispersive TE0: single-solve n_g = {ng_cav:.8} vs 5-pt FD {ng_bloch_fd:.8} (rel = {rel_cav:.2e})"
    );

    // 4. Verify component shares and domain shares sum to 1.0 in 2D modes, 3D bloch, and 3D cavity
    for (tag, res) in [("2D modes", &r2_0), ("3D bloch", &b_0), ("3D cavity", &r_cav)] {
        for m_idx in 1..=res.solutions.len() {
            let ex_s = scalar_num(res, &format!("Ex share {m_idx}")).unwrap();
            let ey_s = scalar_num(res, &format!("Ey share {m_idx}")).unwrap();
            let ez_s = scalar_num(res, &format!("Ez share {m_idx}")).unwrap();
            assert!(
                ((ex_s + ey_s + ez_s) - 1.0).abs() < 1e-12,
                "{tag} mode {m_idx}: component shares sum to {}",
                ex_s + ey_s + ez_s
            );
            let d1_s = scalar_num(res, &format!("share in domain 1 {m_idx}")).unwrap();
            let d2_s = scalar_num(res, &format!("share in domain 2 {m_idx}")).unwrap();
            assert!(((d1_s + d2_s) - 1.0).abs() < 1e-12, "{tag} mode {m_idx}: domain shares sum to {}", d1_s + d2_s);
        }
    }
}

#[test]
fn all_pmc_box_and_zero_phase_periodic_cavity_null_space_pin() {
    // 1. All-PMC 3D box (pec_vertex empty, no periodic pairs): fundamental mode of PMC box [a, b, d]
    //    has same eigenfrequency spectrum c/2 sqrt((m/a)^2 + (n/b)^2 + (p/d)^2) as PEC box!
    let (a, b, d) = (0.04, 0.03, 0.02);
    let mesh = common::block(a, b, d, 5, 4, 3);
    let f_110 = 0.5 * C0 * ((1.0 / a).powi(2) + (1.0 / b).powi(2)).sqrt();
    let mut p_all_pmc = Problem::new(Study::Cavity { n_modes: 3, f_search: 0.92 * f_110 }, 2);
    p_all_pmc.pmc = vec![1, 2, 3, 4, 5, 6];
    let res = common::solve(&mesh, &p_all_pmc).expect("all-PMC cavity solves cleanly with pinned scalar node");
    let got0 = res.solutions[0].value.unwrap();
    assert!(((got0 - f_110) / f_110).abs() < 2e-3, "all-PMC fundamental {got0} vs {f_110}");

    // 2. Gamma-point (phase = 0) 3D periodic cell with PMC transverse walls (pec_vertex empty + periodic phase = 0):
    //    note on S^1_y x [0, a] x [0, d] with no PEC walls, dy is a static cohomology mode at 0 Hz; search near f_100 = c0 / (2 a).
    let f_100 = 0.5 * C0 / a;
    let mut p_gamma = Problem::new(Study::Cavity { n_modes: 3, f_search: 1.02 * f_100 }, 2);
    p_gamma.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase: 0.0 }];
    p_gamma.pmc = vec![1, 3, 4, 6];
    let res_g = common::solve(&mesh, &p_gamma).expect("Gamma-point PMC waveguide cell solves cleanly");
    assert!(
        res_g.solutions.iter().any(|s| ((s.value.unwrap() - f_100) / f_100).abs() < 2e-3),
        "expected f_100 = {f_100:.4e} among {:?}",
        res_g.solutions.iter().map(|s| s.value.unwrap()).collect::<Vec<_>>()
    );
}

#[test]
fn pmc_and_midplane_symmetry_half_and_quarter_meshes_match_full_mesh_te0_and_tm0() {
    let lambda0 = 1.55e-6;
    let freq = C0 / lambda0;
    let (n_si, n_sio2) = (3.476_f64, 1.444_f64);
    let (w_box, h_box) = (1.6e-6, 1.2e-6);
    let (w_si, h_si) = (0.46e-6, 0.22e-6);
    let (h_clad, h_core) = (0.11e-6, 0.030e-6);

    let make_section = |x0: f64, z0: f64| -> (SimplexMesh, Vec<i32>, Vec<i32>) {
        let mut m = Model::new_2d();
        let clad = m.rect("clad", (x0, z0), (0.5 * w_box, 0.5 * h_box));
        let core = m.rect("core", (x0.max(-0.5 * w_si), z0.max(-0.5 * h_si)), (0.5 * w_si, 0.5 * h_si));
        m.overlay(&[("clad", clad), ("core", core)]);
        m.size("sz_core", Selector::parse("region(core)").unwrap(), Some(h_core.into()), None, None);
        let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h_clad)).unwrap();
        let x0_labels = if x0.abs() < 1e-12 { g.select("plane(x = 0)").unwrap() } else { vec![] };
        let z0_labels = if z0.abs() < 1e-12 { g.select("plane(z = 0)").unwrap() } else { vec![] };
        (sort_nodes(&g.mesh), x0_labels, z0_labels)
    };

    let (full2, _, _) = make_section(-0.5 * w_box, -0.5 * h_box);
    let (half2_x, x0_half_labels, _) = make_section(0.0, -0.5 * h_box);
    let (quart2, x0_q_labels, z0_q_labels) = make_section(0.0, 0.0);

    let base_2d = |sym: Symmetry| -> Problem {
        let mut p = Problem::new(Study::Modes { freq, n_target: 2.5, evanescent: false, n_modes: 2 }, 2);
        p.materials.set(1, "eps_r", n_sio2 * n_sio2);
        p.materials.set(2, "eps_r", n_si * n_si);
        p.symmetry = sym;
        p
    };

    // Full 2D reference: mode 0 = TE0, mode 1 = TM0
    let r_full2 = common::solve(&full2, &base_2d(Symmetry::None)).unwrap();
    let neff_te0 = r_full2.solutions[0].eigenvalue.unwrap()[0];
    let neff_tm0 = r_full2.solutions[1].eigenvalue.unwrap()[0];
    let ng_te0 = scalar_num(&r_full2, "n_g 1").unwrap();
    let ng_tm0 = scalar_num(&r_full2, "n_g 2").unwrap();

    // Half mesh (x >= 0): pec_x -> TE0, pmc_x -> TM0, plus explicit `pmc` BC on Modes
    let r_half_te = common::solve(&half2_x, &base_2d(Symmetry::PecX)).unwrap();
    let r_half_tm = common::solve(&half2_x, &base_2d(Symmetry::PmcX)).unwrap();
    assert!(((r_half_te.solutions[0].eigenvalue.unwrap()[0] - neff_te0) / neff_te0).abs() < 1e-4);
    assert!(((r_half_tm.solutions[0].eigenvalue.unwrap()[0] - neff_tm0) / neff_tm0).abs() < 1e-4);

    let mut p_explicit_pmc = base_2d(Symmetry::None);
    p_explicit_pmc.pmc = x0_half_labels;
    let r_half_pmc_bc = common::solve(&half2_x, &p_explicit_pmc).unwrap();
    assert!(((r_half_pmc_bc.solutions[0].eigenvalue.unwrap()[0] - neff_tm0) / neff_tm0).abs() < 1e-4);

    // Quarter mesh (x >= 0, z >= 0): TE0 = pec_x_pmc_z, TM0 = pmc_x_pec_z
    let mut q_te = base_2d(Symmetry::PecXPmcZ);
    q_te.study = Study::Modes { freq, n_target: 2.5, evanescent: false, n_modes: 1 };
    let r_q_te = common::solve(&quart2, &q_te).unwrap();
    assert!(
        ((r_q_te.solutions[0].eigenvalue.unwrap()[0] - neff_te0) / neff_te0).abs() < 1e-4,
        "2D 1/4 TE0 {} vs full {neff_te0}",
        r_q_te.solutions[0].eigenvalue.unwrap()[0]
    );
    assert!(((scalar_num(&r_q_te, "n_g 1").unwrap() - ng_te0) / ng_te0).abs() < 5e-4);

    let mut q_tm = base_2d(Symmetry::PmcXPecZ);
    q_tm.study = Study::Modes { freq, n_target: 2.0, evanescent: false, n_modes: 1 };
    let r_q_tm = common::solve(&quart2, &q_tm).unwrap();
    assert!(
        ((r_q_tm.solutions[0].eigenvalue.unwrap()[0] - neff_tm0) / neff_tm0).abs() < 1e-4,
        "2D 1/4 TM0 {} vs full {neff_tm0}",
        r_q_tm.solutions[0].eigenvalue.unwrap()[0]
    );
    assert!(((scalar_num(&r_q_tm, "n_g 1").unwrap() - ng_tm0) / ng_tm0).abs() < 5e-4);

    // Also test explicit `pmc` BC via ProblemSpec envelope on `modes` (previously forbidden!)
    struct QuartTable {
        x0: Vec<i32>,
        z0: Vec<i32>,
    }
    impl LabelResolver for QuartTable {
        fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
            match r {
                "clad" => Ok(vec![1]),
                "core" => Ok(vec![2]),
                _ => Err(Error::invalid(format!("no domain {r}"))),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            match s {
                "x0" => Ok(self.x0.clone()),
                "z0" => Ok(self.z0.clone()),
                _ => Err(Error::invalid(format!("no boundary {s}"))),
            }
        }
    }
    let q_table = QuartTable { x0: x0_q_labels, z0: z0_q_labels };
    let ex_q2 = Mesh::from_solver_mesh(&quart2).unwrap();
    let mut spec_pmc_modes = ProblemSpec::new("ems", 2, "modes");
    spec_pmc_modes.study.params = json!({ "freq": freq, "n_target": 2.5, "n_modes": 1 });
    spec_pmc_modes.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), n_sio2 * n_sio2)].into()));
    spec_pmc_modes.materials.insert("core".into(), MaterialAssign::Props([("eps_r".into(), n_si * n_si)].into()));
    spec_pmc_modes.bcs.push(BcSpec { on: "z0".into(), kind: "pmc".into(), params: json!({}) });
    let r_env_pmc = Ems.solve(&ex_q2, &spec_pmc_modes, &mut Ctx::new(&q_table)).unwrap();
    assert!(((r_env_pmc.solutions[0].eigenvalue.unwrap()[0] - neff_te0) / neff_te0).abs() < 1e-4);

    // 3D `bloch` & `cavity` quarter-mesh symmetry vs full-mesh reference (pitch = 60 nm keeps phase ~ 0.60 rad per P2 cell)
    let pitch = 0.06e-6;
    let nb_q = quart2.facet_labels.iter().copied().max().unwrap_or(0);
    let quart3 = extrude_y(&quart2, &[Layer { thickness: pitch, n: 1, domain_offset: 0, fill: None }], 0.0).unwrap();

    for (sym, want_neff, target) in [(Symmetry::PecXPmcZ, neff_te0, 2.5), (Symmetry::PmcXPecZ, neff_tm0, 1.9)] {
        let mut pb = Problem::new(Study::Bloch { freq, n_target: target, n_modes: 1, brillouin_zone: 1 }, 2);
        pb.materials.set(1, "eps_r", n_sio2 * n_sio2);
        pb.materials.set(2, "eps_r", n_si * n_si);
        pb.periodic = vec![Periodic { source: vec![nb_q + 1], dest: vec![nb_q + 2], phase: 0.0 }];
        pb.symmetry = sym;
        let rb = common::solve(&quart3, &pb).unwrap();
        let got_nb = rb.solutions[0].eigenvalue.unwrap()[0];
        assert!(
            ((got_nb - want_neff) / want_neff).abs() < 1e-4,
            "3D Bloch 1/4 {:?}: {got_nb:.6} vs full {want_neff:.6}",
            sym
        );

        let phase = want_neff * (2.0 * PI / lambda0) * pitch;
        let mut pc = Problem::new(Study::Cavity { n_modes: 1, f_search: freq }, 2);
        pc.materials = pb.materials.clone();
        pc.periodic = vec![Periodic { source: vec![nb_q + 1], dest: vec![nb_q + 2], phase }];
        pc.symmetry = sym;
        let rc = common::solve(&quart3, &pc).unwrap();
        let got_f = rc.solutions[0].value.unwrap();
        assert!(((got_f - freq) / freq).abs() < 1e-4, "3D Cavity 1/4 {:?}: {got_f:.6e} vs {freq:.6e}", sym);
    }
}
