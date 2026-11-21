// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative verification of Milestone 3.2 (`surface_impedance` KICQ IBC + `pec` sheet),
//! Milestone 3.3 (`lumped_jj` rank-1 operator), Milestone 4.1 (`Hx, Hy, Hz`, `Kx, Ky, Kz`,
//! `U_elec = 1 J` normalization), Milestone 4.3 (passive Craig-Bampton `CircuitRom`), and
//! Acceptance Benchmark Test 3 (Park et al., npj Quantum Inf. 2026, Fig. 1d-e).

use femcore::material::{bcs_surface_impedance, FilmModel, Library, Superconducting, EV_TO_J};
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_ems::ibc::MU0;
use nyx_ems::{rom_eigenfrequencies, Ems};
use nyx_physics::{BcSpec, CircuitRom, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, TlsMeta};
use serde_json::json;
use std::collections::BTreeMap;
use std::sync::Arc;

struct GeomResolver<'a>(&'a MeshedGeom);

impl LabelResolver for GeomResolver<'_> {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.domains(region);
        if v.is_empty() {
            Err(Error::invalid(format!("no region '{region}'")))
        } else {
            Ok(v)
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.boundary(selector);
        if !v.is_empty() {
            return Ok(v);
        }
        self.0.select(selector).map_err(Error::invalid)
    }
}

fn park_nb_spec() -> Superconducting {
    Superconducting {
        source: Vec::new(),
        tc: 9.2,
        gap_0: Some(1.395e-3 * EV_TO_J),
        lambda_l: Some(33.3e-9),
        xi_0: Some(39.0e-9),
        mean_free_path: None,
        sigma_n: Some(5.5e7),
        ..Default::default()
    }
}

fn park_nb_library() -> Arc<Library> {
    let toml_text = r#"
format = 2
name = "Nb_Park2026"
description = "100 nm Nb thin film (Park et al. 2026)"

[superconducting]
tc = "9.2 K"
gap_0 = "2.23503640443e-22 J"
lambda_L = "33.3e-9 m"
xi_0 = "39.0e-9 m"
sigma_n = "5.5e7 S/m"
"#;
    let mut texts = BTreeMap::new();
    texts.insert("nb_park2026".to_string(), toml_text.to_string());
    Arc::new(Library::default().with_texts(&texts).expect("valid Nb TOML"))
}

/// Builds a 3D lambda/2 CPW resonator on 500 um Si (`eps_r = 11.45`) using zero-thickness interior 2D sheets
/// (`Model::sheet`) at `z = 0`.
fn build_cpw_3d_mesh(
    w_trace: f64,
    s_gap: f64,
    w_box: f64,
    l_res: f64,
    h_sub: f64,
    h_air: f64,
    h_mesh: f64,
) -> MeshedGeom {
    let mut m = Model::new_3d();
    let si = m.box_("si_box", (-w_box / 2.0, 0.0, -h_sub), (w_box / 2.0, l_res, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, 0.0, 0.0), (w_box / 2.0, l_res, h_air));

    let trace_2d = m.rect("trace_2d", (-w_trace / 2.0, 0.0), (w_trace / 2.0, l_res));
    m.sheet("cpw_trace", trace_2d, 0.0);

    let gnd_l_2d = m.rect("gnd_l_2d", (-w_box / 2.0, 0.0), (-w_trace / 2.0 - s_gap, l_res));
    m.sheet("cpw_gnd_l", gnd_l_2d, 0.0);

    let gnd_r_2d = m.rect("gnd_r_2d", (w_trace / 2.0 + s_gap, 0.0), (w_box / 2.0, l_res));
    m.sheet("cpw_gnd_r", gnd_r_2d, 0.0);

    m.overlay(&[("si", si), ("air", air)]);
    mesh_model(&m.list, &MeshPlan::simple(h_mesh)).expect("conformal 3D sheet-stack CPW mesh")
}

#[test]
fn park_2026_kicq_3d_cpw_resonator_shift_and_bcs_temperature_sweep() {
    let nb = park_nb_spec();
    let lib = park_nb_library();
    let f_ref = 7.064e9;
    let d_film = 100e-9;
    let eval_0 = bcs_surface_impedance(&nb, f_ref, 0.01, d_film, FilmModel::Kerr).unwrap();
    eprintln!(
        "Nb (100 nm, Kerr, 10 mK): L_k_sq = {:.6e} H/sq, R_s = {:.6e} Ohm/sq, lambda_eff = {:.2} nm",
        eval_0.l_k_sq,
        eval_0.r_s,
        eval_0.lambda_eff * 1e9
    );

    // 3D lambda/2 CPW on 500 um Si (eps_r = 11.45) with zero-thickness interior sheets at z = 0
    let h_sub = 500e-6;
    let h_air = 500e-6;
    let w_trace = 49.3e-6;
    let s_gap = 28.95e-6;
    let w_box = 297.0e-6;
    let l_res = 8.080e-3;
    let h_mesh = 1.2e-3;

    let g = build_cpw_3d_mesh(w_trace, s_gap, w_box, l_res, h_sub, h_air, h_mesh);
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);
    eprintln!(
        "3D CPW sheet mesh: {} cells, {} points, {} facets",
        g.mesh.n_cells(),
        g.mesh.n_points(),
        g.mesh.n_facets()
    );

    // 1. Solve with `pec` on the CPW center trace and ground planes
    let mut spec_pec = ProblemSpec::new("ems", 1, "cavity");
    spec_pec.study.params = json!({ "n_modes": 1, "f_search": 7.0e9 });
    spec_pec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec_pec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
        spec_pec.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
    }
    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib.clone());
    let res_pec = Ems.solve(&exchange, &spec_pec, &mut ctx).unwrap();
    let f_pec = res_pec.solutions[0].value.unwrap();
    eprintln!("f_r(PEC) = {:.6} GHz (target 7.064 GHz)", f_pec / 1e9);

    // 2. Solve with `surface_impedance` at T = 10 mK via the material library ("nb_park2026", "kerr")
    let solve_at_temp = |temp_k: f64| -> nyx_physics::SolveResult {
        let mut spec_sc = ProblemSpec::new("ems", 1, "cavity");
        spec_sc.study.params = json!({ "n_modes": 1, "f_search": 6.9e9 });
        spec_sc.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        spec_sc.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sheet in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
            spec_sc.bcs.push(BcSpec {
                on: sheet.into(),
                kind: "surface_impedance".into(),
                params: json!({
                    "material": "nb_park2026",
                    "thickness": d_film,
                    "temperature": temp_k,
                    "f_ref": f_ref,
                    "film_model": "kerr"
                }),
            });
        }
        let mut ctx = Ctx::new(&resolver);
        ctx.library = Some(lib.clone());
        Ems.solve(&exchange, &spec_sc, &mut ctx).unwrap()
    };

    let res_0 = solve_at_temp(0.01);
    let f_0 = res_0.solutions[0].value.unwrap();
    let shift_0_mhz = (f_0 - f_pec) / 1e6;
    eprintln!(
        "f_r(10 mK) = {:.6} GHz, downward shift = {:.2} MHz (target ~ -154 MHz to ~ 6.910 GHz)",
        f_0 / 1e9,
        shift_0_mhz
    );

    assert!(
        (f_pec - 7.064e9).abs() / 7.064e9 < 0.005,
        "PEC resonant frequency {:.6} GHz outside 7.064 GHz +/- 0.5%",
        f_pec / 1e9
    );
    assert!(
        (f_0 - 6.910e9).abs() / 6.910e9 < 0.005,
        "10 mK superconducting resonant frequency {:.6} GHz outside 6.910 GHz +/- 0.5%",
        f_0 / 1e9
    );
    assert!(
        (shift_0_mhz - (-154.0)).abs() < 10.0,
        "0 K kinetic inductance shift {:.2} MHz should be near -154 MHz",
        shift_0_mhz
    );

    // Verify TlsMeta and CircuitRom emission on the surface_impedance solve
    let tls = TlsMeta::of(&res_0).expect("TlsMeta emitted on surface_impedance cavity solve");
    assert_eq!(tls.domain_eps_r.get(&g.domain("si").unwrap()), Some(&11.45));
    assert!(!tls.sheet_labels.is_empty(), "interior sheet labels recorded in TlsMeta");
    let rom_0 = CircuitRom::of(&res_0).expect("CircuitRom emitted on surface_impedance cavity solve");
    assert_eq!(rom_0.modes.len(), res_0.solutions.len());
    assert!((rom_0.modes[0].u_elec_j - 1.0).abs() < 1e-12);
    assert!(rom_0.modes[0].u_kin_ibc_j > 0.02, "kinetic energy fraction > 2%: {:?}", rom_0.modes[0]);
    let eq_err_0 = (rom_0.modes[0].u_elec_j - (rom_0.modes[0].u_mag_vol_j + rom_0.modes[0].u_kin_ibc_j)).abs();
    assert!(eq_err_0 < 1e-12, "equipartition error at 10 mK: {eq_err_0:e}");

    // 3. Temperature sweep over T / T_c in [0.01, 0.85] against analytical BCS kinetic inductance shift
    // Since K_IBC(T) = (1 / L_eff(T)) K_geom, the analytical BCS curve anchored at (f_pec, f_0) is
    // f_BCS(T) = f_pec / sqrt(1 + gamma_geom * L_eff(T)) where gamma_geom = ((f_pec / f_0)^2 - 1) / L_eff(T_0).
    let gamma_geom = ((f_pec / f_0).powi(2) - 1.0) / eval_0.l_eff;
    let t_over_tc_list = [0.01, 0.25, 0.50, 0.65, 0.75, 0.85];
    let mut max_bcs_residual_mhz = 0.0_f64;

    for &t_ratio in &t_over_tc_list {
        let temp_k = t_ratio * nb.tc;
        let bcs = bcs_surface_impedance(&nb, f_ref, temp_k, d_film, FilmModel::Kerr).unwrap();
        let f_bcs_analytic = f_pec / (1.0 + gamma_geom * bcs.l_eff).sqrt();
        let res_t = solve_at_temp(temp_k);
        let f_fem = res_t.solutions[0].value.unwrap();
        let diff_mhz = (f_fem - f_bcs_analytic).abs() / 1e6;
        max_bcs_residual_mhz = max_bcs_residual_mhz.max(diff_mhz);
        eprintln!(
            "  T/Tc = {:.2} (T = {:.3} K): L_k_sq = {:.4e} H/sq, f_FEM = {:.6} GHz, f_BCS = {:.6} GHz, |diff| = {:.4} MHz",
            t_ratio,
            temp_k,
            bcs.l_k_sq,
            f_fem / 1e9,
            f_bcs_analytic / 1e9,
            diff_mhz
        );
        assert!(
            diff_mhz < 5.0,
            "T/Tc = {t_ratio}: FEM frequency {:.6} GHz differs from analytical BCS {:.6} GHz by {:.3} MHz (>= 5 MHz)",
            f_fem / 1e9,
            f_bcs_analytic / 1e9,
            diff_mhz
        );
    }
    eprintln!(
        "Max BCS temperature-sweep residual across T/Tc in [0.01, 0.85]: {:.4} MHz (< 5 MHz)",
        max_bcs_residual_mhz
    );
}

#[test]
fn lumped_jj_rank1_equipartition_physical_fields_and_circuit_rom() {
    // 3D planar transmon / lumped-JJ resonator inside a 3D cavity:
    // two superconducting pads (`cpw_pad1`, `cpw_pad2`) bridged by a `lumped_jj` sheet (`Q1_JJ`) at z = 0
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.2e-3, 1.6e-3, 300e-6, 300e-6);
    let si = m.box_("si_box", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));

    // Pad 1: y in [-500 um, -40 um], Pad 2: y in [+40 um, +500 um], JJ bridge: y in [-40 um, +40 um], x in [-40 um, +40 um]
    let pad1_2d = m.rect("pad1_2d", (-160e-6, -500e-6), (160e-6, -40e-6));
    m.sheet("pad1", pad1_2d, 0.0);

    let pad2_2d = m.rect("pad2_2d", (-160e-6, 40e-6), (160e-6, 500e-6));
    m.sheet("pad2", pad2_2d, 0.0);

    let jj_2d = m.rect("jj_2d", (-40e-6, -40e-6), (40e-6, 40e-6));
    m.sheet("jj_sheet", jj_2d, 0.0);

    m.overlay(&[("si", si), ("air", air)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(220e-6)).expect("3D transmon sheet mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let mut spec = ProblemSpec::new("ems", 2, "cavity");
    spec.study.params = json!({ "n_modes": 1, "f_search": 5.5e9 });
    spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    spec.bcs.push(BcSpec {
        on: "pad1".into(),
        kind: "surface_impedance".into(),
        params: json!({ "L_k_sq": 1.0e-13, "R_s": 0.0 }),
    });
    spec.bcs.push(BcSpec {
        on: "pad2".into(),
        kind: "surface_impedance".into(),
        params: json!({ "L_k_sq": 1.0e-13, "R_s": 0.0 }),
    });
    spec.bcs.push(BcSpec {
        on: "jj_sheet".into(),
        kind: "lumped_jj".into(),
        params: json!({
            "name": "Q1_JJ",
            "E_j_Hz": 18.0e9,
            "C_j": 35.0e-15,
            "E_l_Hz": 2.0e9,
            "phi_ext_rad": 0.25,
            "n_g": 0.1,
            "direction": [0.0, 1.0, 0.0],
            "formulation": "rank1"
        }),
    });

    let mut ctx = Ctx::new(&resolver);
    let res = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
    assert!(!res.solutions.is_empty());

    // 1. Check equipartition_error < 1e-12
    let eq_scalar = res.scalar("equipartition_error").expect("equipartition_error scalar recorded");
    eprintln!("Lumped JJ rank-1 max equipartition_error = {eq_scalar:.3e}");
    assert!(eq_scalar < 1e-12, "equipartition_error {eq_scalar:e} >= 1e-12");

    let rom = CircuitRom::of(&res).expect("CircuitRom present in SolveResult::scalars");
    assert_eq!(rom.ports.len(), 1);
    assert_eq!(rom.ports[0].name, "Q1_JJ");
    assert!(rom.ports[0].l_j_h.unwrap() > 0.0);
    assert!((rom.ports[0].e_j_hz.unwrap() - 18.0e9).abs() < 1.0);
    assert!((rom.ports[0].e_l_hz.unwrap() - 2.0e9).abs() < 1.0);
    assert!(
        rom.c_port_f.len() == 1 && rom.c_port_f[0][0] > 1.0e-15,
        "3D pad-to-pad C_port = {:.3e} F",
        rom.c_port_f[0][0]
    );
    eprintln!(
        "CircuitRom: r_static = {}, r_dyn = {}, C_port = {:.4} fF",
        rom.r_static,
        rom.r_dyn,
        rom.c_port_f[0][0] * 1e15
    );

    for (m_idx, row) in rom.modes.iter().enumerate() {
        let tot_mag = row.u_mag_vol_j + row.u_kin_ibc_j + row.u_ind_jj_j.iter().sum::<f64>();
        let err = (row.u_elec_j - tot_mag).abs() / row.u_elec_j;
        eprintln!(
            "  mode {}: f = {:.6} GHz, U_elec = {:.6} J, U_mag_vol = {:.6} J, U_kin_ibc = {:.6} J, U_ind_jj = {:.6} J, equip_err = {:.3e}",
            m_idx + 1,
            row.freq_hz / 1e9,
            row.u_elec_j,
            row.u_mag_vol_j,
            row.u_kin_ibc_j,
            row.u_ind_jj_j[0],
            err
        );
        assert!((row.u_elec_j - 1.0).abs() < 1e-12);
        assert!(err < 1e-12, "mode {} equipartition error {err:e} >= 1e-12", m_idx + 1);
    }
    assert!(rom.modes[0].u_ind_jj_j[0] > 0.1, "fundamental transmon mode has significant JJ inductive participation");

    // 2. Check physical CellFields (`Ex, Ey, Ez` in "V/m", `Hx, Hy, Hz` in "A/m", `Kx, Ky, Kz` in "A/m")
    // and verify that integrating (1/4) mu_0 |H|^2 over the mesh matches `u_mag_vol_j`!
    let sol0 = &res.solutions[0];
    for name in ["Ex", "Ey", "Ez"] {
        let cf = sol0.cell_fields.iter().find(|f| f.name == name).expect(name);
        assert_eq!(cf.unit, "V/m");
    }
    for name in ["Hx", "Hy", "Hz", "Kx", "Ky", "Kz"] {
        let cf = sol0.cell_fields.iter().find(|f| f.name == name).expect(name);
        assert_eq!(cf.unit, "A/m");
    }
    let kx_cf = sol0.cell_fields.iter().find(|f| f.name == "Ky").unwrap();
    let max_ky = kx_cf.re.iter().zip(&kx_cf.im).map(|(&r, &im)| r.hypot(im)).fold(0.0_f64, f64::max);
    assert!(max_ky > 0.0, "sheet surface current Ky is non-zero on cells adjacent to sheets");

    // Verify Hx, Hy, Hz physical accuracy: integrating (1/4) mu_0 |H|^2 over the mesh matches u_mag_vol_j
    let hx_cf = sol0.cell_fields.iter().find(|f| f.name == "Hx").unwrap();
    let hy_cf = sol0.cell_fields.iter().find(|f| f.name == "Hy").unwrap();
    let hz_cf = sol0.cell_fields.iter().find(|f| f.name == "Hz").unwrap();
    let geo = femcore::fem::geometry(&g.mesh);
    let quad = femcore::fem::quad_rule(3, 4);
    let nloc = 10;
    let mut u_mag_integrated = 0.0;
    for ci in 0..g.mesh.n_cells() {
        for q in 0..quad.w.len() {
            let (sh, _) = femcore::fem::shape(&quad.lam[q], 2);
            let (mut hx_re, mut hx_im) = (0.0, 0.0);
            let (mut hy_re, mut hy_im) = (0.0, 0.0);
            let (mut hz_re, mut hz_im) = (0.0, 0.0);
            for (a, &na) in sh.iter().enumerate().take(nloc) {
                let idx = ci * nloc + a;
                hx_re += na * hx_cf.re[idx];
                hx_im += na * hx_cf.im[idx];
                hy_re += na * hy_cf.re[idx];
                hy_im += na * hy_cf.im[idx];
                hz_re += na * hz_cf.re[idx];
                hz_im += na * hz_cf.im[idx];
            }
            let h2 = hx_re * hx_re + hx_im * hx_im + hy_re * hy_re + hy_im * hy_im + hz_re * hz_re + hz_im * hz_im;
            u_mag_integrated += 0.25 * MU0 * quad.w[q] * geo.vol[ci] * h2;
        }
    }
    let rel_h_err = (u_mag_integrated - rom.modes[0].u_mag_vol_j).abs() / rom.modes[0].u_mag_vol_j.max(1e-30);
    eprintln!(
        "Integrated (1/4) mu_0 |H|^2 = {:.6e} J vs exact u_mag_vol_j = {:.6e} J (rel diff {:.3e})",
        u_mag_integrated, rom.modes[0].u_mag_vol_j, rel_h_err
    );
    assert!(rel_h_err < 1e-9, "H field energy integral disagrees with u_mag_vol_j: {rel_h_err:e}");

    // 3. Verify CircuitRom eigenvalues match the direct 3D FEM eigenvalues!
    let rom_freqs = rom_eigenfrequencies(&rom).unwrap();
    eprintln!("ROM eigenfrequencies (GHz): {:?}", rom_freqs.iter().map(|f| f / 1e9).collect::<Vec<_>>());
    for row in &rom.modes {
        let nearest_rel =
            rom_freqs.iter().map(|&rf| (rf - row.freq_hz).abs() / row.freq_hz).fold(f64::INFINITY, f64::min);
        assert!(
            nearest_rel < 1e-9,
            "FEM mode {:.6} GHz not reproduced by CircuitRom eigenvalues {:?} (nearest rel diff {nearest_rel:e})",
            row.freq_hz / 1e9,
            rom_freqs
        );
    }
}
