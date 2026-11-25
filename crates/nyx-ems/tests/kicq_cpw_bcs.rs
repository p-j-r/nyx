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
use nyx_ems::{rom_complex_eigenfrequencies, rom_eigenfrequencies, Ems};
use nyx_physics::{
    BcSpec, CircuitRom, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, TlsMeta,
};
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

#[test]
fn fr12_cavity_lumped_port_qep_vs_rayleigh_vs_driven_and_f_search_invariance() {
    // 3D coupled feedline + transmon/resonator chip with lumped_jj + 50-Ohm lumped_port(s)
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.4e-3, 1.4e-3, 300e-6, 300e-6);
    let si = m.box_("si_box", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));
    let p1_2d = m.rect("p1_2d", (-520e-6, -700e-6), (-340e-6, -460e-6));
    m.sheet("port1_sheet", p1_2d, 0.0);
    let feed_2d = m.rect("feed_2d", (-520e-6, -460e-6), (-340e-6, 460e-6));
    m.sheet("feed_trace", feed_2d, 0.0);
    let p2_2d = m.rect("p2_2d", (-520e-6, 460e-6), (-340e-6, 700e-6));
    m.sheet("port2_sheet", p2_2d, 0.0);
    let pad1_2d = m.rect("pad1_2d", (-140e-6, -500e-6), (420e-6, -60e-6));
    m.sheet("res_pad1", pad1_2d, 0.0);
    let pad2_2d = m.rect("pad2_2d", (-140e-6, 60e-6), (420e-6, 500e-6));
    m.sheet("res_pad2", pad2_2d, 0.0);
    let jj_2d = m.rect("jj_2d", (40e-6, -60e-6), (240e-6, 60e-6));
    m.sheet("jj_sheet", jj_2d, 0.0);
    m.overlay(&[("si", si), ("air", air)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(240e-6)).expect("coupled resonator mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let e_j_hz = 20.0e9;
    let c_j_f = 35.0e-15;

    let build_cav_spec = |f_search: f64, complex_damping: bool| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "cavity");
        s.study.params = json!({ "n_modes": 1, "f_search": f_search });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
            s.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
        }
        s.bcs.push(BcSpec {
            on: "jj_sheet".into(),
            kind: "lumped_jj".into(),
            params: json!({
                "name": "Q1_JJ",
                "E_j_Hz": e_j_hz,
                "C_j": c_j_f,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }),
        });
        s.bcs.push(BcSpec {
            on: "port1_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({
                "name": "port1",
                "Z_0": 50.0,
                "direction": [0.0, 1.0, 0.0],
                "complex_damping": complex_damping
            }),
        });
        s.bcs.push(BcSpec {
            on: "port2_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({
                "name": "port2",
                "Z_0": 50.0,
                "direction": [0.0, -1.0, 0.0],
                "complex_damping": complex_damping
            }),
        });
        s
    };

    // 1. Real cavity path (`complex_damping: false`) + CircuitRom complex QEP (including f_search +/- 25% invariance)
    let res_real = Ems.solve(&exchange, &build_cav_spec(3.2e9, false), &mut Ctx::new(&resolver)).unwrap();
    let rom_real = CircuitRom::of(&res_real).expect("CircuitRom emitted");
    let ev_real = res_real.solutions[0].eigenvalue.unwrap();
    assert_eq!(ev_real[1], 0.0, "complex_damping: false must take the purely real eigensolve path");
    assert_eq!(rom_real.ports.len(), 3);
    assert_eq!(rom_real.ports[0].kind, "lumped_jj");
    assert_eq!(rom_real.ports[0].z0_ohm, None);
    assert_eq!(rom_real.ports[1].kind, "lumped_port");
    assert_eq!(rom_real.ports[1].z0_ohm, Some(50.0));
    assert_eq!(rom_real.ports[2].kind, "lumped_port");
    assert_eq!(rom_real.ports[2].z0_ohm, Some(50.0));
    assert_eq!(rom_real.c_port_f.len(), 3);
    assert!(!rom_real.c_r.is_empty(), "CircuitRom::c_r must contain lumped_port damping");

    let qep_freqs = rom_complex_eigenfrequencies(&rom_real).expect("rom_complex_eigenfrequencies");
    let f_qep = *qep_freqs
        .iter()
        .min_by(|a, b| (a.re - ev_real[0]).abs().total_cmp(&(b.re - ev_real[0]).abs()))
        .expect("QEP mode near cavity resonance");
    let q_qep = f_qep.re / (2.0 * f_qep.im.abs());

    // 2. Complex Rayleigh path (`complex_damping: true`) at f0
    let f0 = ev_real[0];
    let res_ray_0 = Ems.solve(&exchange, &build_cav_spec(f0, true), &mut Ctx::new(&resolver)).unwrap();
    let ev_ray_0 = res_ray_0.solutions[0].eigenvalue.unwrap();
    let q_ray_0 = ev_ray_0[0] / (2.0 * ev_ray_0[1].abs());

    // Verify equipartition_error < 1e-12 on damped complex cavity (FR-22)
    let eq_damped = res_ray_0.scalar("equipartition_error").unwrap();
    assert!(eq_damped < 1e-12, "complex damped cavity equipartition_error = {eq_damped:e}");

    // 3. Driven S21 complex-circle Lorentzian linewidth comparison (exact Fano-free |S21(f) - S21_bg|^2)
    let half_bw = 25.0 * (f0 / q_qep);
    let mut s_drv = ProblemSpec::new("ems", 1, "driven");
    s_drv.study.params = json!({
        "f_min": f0 - half_bw,
        "f_max": f0 + half_bw,
        "n_points": 1001,
        "prom_tol": 1e-6,
        "strict_prom_tol": true
    });
    s_drv.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    s_drv.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
        s_drv.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
    }
    s_drv.bcs.push(BcSpec {
        on: "jj_sheet".into(),
        kind: "lumped_jj".into(),
        params: json!({
            "name": "Q1_JJ",
            "E_j_Hz": e_j_hz,
            "C_j": c_j_f,
            "direction": [0.0, 1.0, 0.0],
            "formulation": "rank1"
        }),
    });
    s_drv.bcs.push(BcSpec {
        on: "port1_sheet".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
    });
    s_drv.bcs.push(BcSpec {
        on: "port2_sheet".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }),
    });
    let res_drv = Ems.solve(&exchange, &s_drv, &mut Ctx::new(&resolver)).unwrap();
    let d_drv = DrivenPortData::of(&res_drv).unwrap();
    // |S21(f) - S21_bg|^2 is an exact symmetric Lorentzian regardless of Fano phase rotation
    let n_pts = d_drv.s_matrix.len();
    let s21_bg = [
        0.5 * (d_drv.s_matrix[0][1][0][0] + d_drv.s_matrix[n_pts - 1][1][0][0]),
        0.5 * (d_drv.s_matrix[0][1][0][1] + d_drv.s_matrix[n_pts - 1][1][0][1]),
    ];
    let lor: Vec<f64> = d_drv
        .s_matrix
        .iter()
        .map(|s| {
            let dr = s[1][0][0] - s21_bg[0];
            let di = s[1][0][1] - s21_bg[1];
            dr * dr + di * di
        })
        .collect();
    let (k_pk, &lor_max) = lor.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap();
    let half = 0.5 * lor_max;
    let mut f_lo_hm = None;
    for k in (0..k_pk).rev() {
        if lor[k] <= half && lor[k + 1] >= half {
            let t = (half - lor[k]) / (lor[k + 1] - lor[k]);
            f_lo_hm = Some(d_drv.freqs_hz[k] + t * (d_drv.freqs_hz[k + 1] - d_drv.freqs_hz[k]));
            break;
        }
    }
    let mut f_hi_hm = None;
    for k in k_pk..(lor.len() - 1) {
        if lor[k] >= half && lor[k + 1] <= half {
            let t = (lor[k] - half) / (lor[k] - lor[k + 1]);
            f_hi_hm = Some(d_drv.freqs_hz[k] + t * (d_drv.freqs_hz[k + 1] - d_drv.freqs_hz[k]));
            break;
        }
    }
    let q_driven = d_drv.freqs_hz[k_pk] / (f_hi_hm.unwrap() - f_lo_hm.unwrap());

    eprintln!(
        "FR-12 coupled feedline+transmon: f_QEP = {:.6} GHz (Q={:.1}), f_Ray = {:.6} GHz (Q={:.1}), f_drv = {:.6} GHz (Q={:.1})",
        f_qep.re / 1e9,
        q_qep,
        ev_ray_0[0] / 1e9,
        q_ray_0,
        d_drv.freqs_hz[k_pk] / 1e9,
        q_driven
    );

    assert!((ev_ray_0[0] - f_qep.re).abs() / f_qep.re < 1e-2, "Rayleigh vs QEP freq > 1%");
    assert!((q_ray_0 - q_qep).abs() / q_qep < 1e-2, "Rayleigh vs QEP Q > 1%");
    assert!((d_drv.freqs_hz[k_pk] - f_qep.re).abs() / f_qep.re < 1e-2, "Driven vs QEP freq > 1%");
    assert!((q_driven - q_qep).abs() / q_qep < 1e-2, "Driven vs QEP Q > 1%: {q_driven:.1} vs {q_qep:.1}");

    // 4. f_search +/- 25% invariance of mode-exact Rayleigh damping quotient (`ruiz_unscale_and_rayleigh_damped`)
    // on both (a) resistive surface_impedance CPW cavity and (b) direct-shunt lumped_port(complex_damping: true) cavity
    let g_cpw = build_cpw_3d_mesh(49.3e-6, 28.95e-6, 297.0e-6, 8.080e-3, 500e-6, 500e-6, 1.2e-3);
    let ex_cpw = Mesh::from_solver_mesh(&g_cpw.mesh).unwrap();
    let resv_cpw = GeomResolver(&g_cpw);
    let solve_cpw_ray = |f_s: f64| -> (f64, f64) {
        let mut s = ProblemSpec::new("ems", 1, "cavity");
        s.study.params = json!({ "n_modes": 1, "f_search": f_s });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sheet in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
            s.bcs.push(BcSpec {
                on: sheet.into(),
                kind: "surface_impedance".into(),
                params: json!({ "L_k_sq": 1.0e-13, "R_s": 5.0e-4, "f_ref": 7.0e9 }),
            });
        }
        let r = Ems.solve(&ex_cpw, &s, &mut Ctx::new(&resv_cpw)).unwrap();
        let ev = r.solutions[0].eigenvalue.unwrap();
        (ev[0], ev[0] / (2.0 * ev[1].abs()))
    };
    let (f_cpw_0, q_cpw_0) = solve_cpw_ray(6.9e9);
    let (_, q_cpw_lo) = solve_cpw_ray(0.75 * f_cpw_0);
    let (_, q_cpw_hi) = solve_cpw_ray(1.25 * f_cpw_0);
    eprintln!(
        "FR-12 Rayleigh f_search +/-25% invariance (IBC CPW): Q(f0)={:.2}, Q(0.75*f0)={:.2}, Q(1.25*f0)={:.2}",
        q_cpw_0, q_cpw_lo, q_cpw_hi
    );
    assert!((q_cpw_lo - q_cpw_0).abs() / q_cpw_0 < 1e-3, "IBC Q(0.75*f0) drifted: {q_cpw_lo} vs {q_cpw_0}");
    assert!((q_cpw_hi - q_cpw_0).abs() / q_cpw_0 < 1e-3, "IBC Q(1.25*f0) drifted: {q_cpw_hi} vs {q_cpw_0}");

    // Direct-shunt lumped_port on box cavity face `y0` (high Z_0 = 50 kOhm -> Q ~ 150) at f_search +/- 25%
    let box_mesh = femcore::testkit::block(20e-3, 20e-3, 10e-3, 4, 4, 2);
    let box_ex = Mesh::from_solver_mesh(&box_mesh).unwrap();
    struct BoxResv;
    impl LabelResolver for BoxResv {
        fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
            match region {
                "air" => Ok(vec![1]),
                _ => Err(Error::invalid(region)),
            }
        }
        fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
            match selector {
                "y0" => Ok(vec![2]),
                _ => Err(Error::invalid(selector)),
            }
        }
    }
    let solve_box_lp = |f_s: f64| -> (f64, f64) {
        let mut s = ProblemSpec::new("ems", 2, "cavity");
        s.study.params = json!({ "n_modes": 1, "f_search": f_s });
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        s.bcs.push(BcSpec {
            on: "y0".into(),
            kind: "lumped_port".into(),
            params: json!({
                "name": "shunt_p",
                "Z_0": 50_000.0,
                "direction": [0.0, 0.0, 1.0],
                "complex_damping": true
            }),
        });
        let r = Ems.solve(&box_ex, &s, &mut Ctx::new(&BoxResv)).unwrap();
        let ev = r.solutions[0].eigenvalue.unwrap();
        (ev[0], ev[0] / (2.0 * ev[1].abs()))
    };
    let (f_box_0, q_box_0) = solve_box_lp(8.38e9);
    let (_, q_box_lo) = solve_box_lp(0.75 * f_box_0);
    let (_, q_box_hi) = solve_box_lp(1.25 * f_box_0);
    eprintln!(
        "FR-12 Rayleigh f_search +/-25% invariance (lumped_port box): Q(f0)={:.2}, Q(0.75*f0)={:.2}, Q(1.25*f0)={:.2}",
        q_box_0, q_box_lo, q_box_hi
    );
    assert!((q_box_lo - q_box_0).abs() / q_box_0 < 1e-3, "LP box Q(0.75*f0) drifted: {q_box_lo} vs {q_box_0}");
    assert!((q_box_hi - q_box_0).abs() / q_box_0 < 1e-3, "LP box Q(1.25*f0) drifted: {q_box_hi} vs {q_box_0}");

    // 4. FR-12 Pt 3: CircuitRom PML exclusion (add a PML cap above air and verify ROM builds cleanly from physical cells)
    let mut m_pml = Model::new_3d();
    let si_p = m_pml.box_("si_box", (-400e-6, -500e-6, -200e-6), (400e-6, 500e-6, 0.0));
    let air_p = m_pml.box_("air_box", (-400e-6, -500e-6, 0.0), (400e-6, 500e-6, 200e-6));
    let pml_p = m_pml.box_("pml_box", (-400e-6, -500e-6, 200e-6), (400e-6, 500e-6, 350e-6));
    let pad1_p = m_pml.rect("pad1_p", (-120e-6, -350e-6), (120e-6, -40e-6));
    m_pml.sheet("pad1", pad1_p, 0.0);
    let pad2_p = m_pml.rect("pad2_p", (-120e-6, 40e-6), (120e-6, 350e-6));
    m_pml.sheet("pad2", pad2_p, 0.0);
    let jj_p = m_pml.rect("jj_p", (-40e-6, -40e-6), (40e-6, 40e-6));
    m_pml.sheet("jj_sheet", jj_p, 0.0);
    m_pml.overlay(&[("si", si_p), ("air", air_p), ("pml_cap", pml_p)]);
    let g_pml = mesh_model(&m_pml.list, &MeshPlan::simple(220e-6)).expect("PML ROM mesh");
    let ex_pml = Mesh::from_solver_mesh(&g_pml.mesh).unwrap();
    let resv_pml = GeomResolver(&g_pml);

    let mut s_pml = ProblemSpec::new("ems", 1, "cavity");
    s_pml.study.params = json!({ "n_modes": 1, "f_search": 5.5e9 });
    s_pml.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    s_pml.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    s_pml.materials.insert("pml_cap".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    s_pml.bcs.push(BcSpec { on: "pml_cap".into(), kind: "pml".into(), params: json!({ "z": "plus" }) });
    s_pml.bcs.push(BcSpec { on: "pad1".into(), kind: "pec".into(), params: json!({}) });
    s_pml.bcs.push(BcSpec { on: "pad2".into(), kind: "pec".into(), params: json!({}) });
    s_pml.bcs.push(BcSpec {
        on: "jj_sheet".into(),
        kind: "lumped_jj".into(),
        params: json!({ "name": "Q1", "E_j_Hz": 18.0e9, "C_j": 35.0e-15, "direction": [0.0, 1.0, 0.0] }),
    });
    let res_pml = Ems.solve(&ex_pml, &s_pml, &mut Ctx::new(&resv_pml)).unwrap();
    let rom_pml = CircuitRom::of(&res_pml).expect("CircuitRom built with PML present");
    let pml_rom_freqs = rom_eigenfrequencies(&rom_pml).expect("ROM eigenfrequencies with PML excluded");
    assert!(!pml_rom_freqs.is_empty());
    let f_pml_fem = res_pml.solutions[0].value.unwrap();
    let nearest_pml = pml_rom_freqs.iter().map(|&rf| (rf - f_pml_fem).abs() / f_pml_fem).fold(f64::INFINITY, f64::min);
    assert!(nearest_pml < 0.02, "PML-excluded ROM eigenvalue differs from FEM by {nearest_pml:e}");
}

#[test]
fn fr18_fr21_fr22_fr23_lossless_al_sc_squid_asymmetry_q50_and_tls_meta() {
    // Load the real material library from the workspace root (`../../materials`)
    let mat_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    let lib = Arc::new(Library::load_dir(&mat_dir).expect("load workspace materials library"));

    // Build 3D transmon mesh (`pad1`, `pad2`, `jj_sheet`)
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.2e-3, 1.6e-3, 300e-6, 300e-6);
    let si = m.box_("si_box", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));
    let pad1_2d = m.rect("pad1_2d", (-160e-6, -500e-6), (160e-6, -40e-6));
    m.sheet("pad1", pad1_2d, 0.0);
    let pad2_2d = m.rect("pad2_2d", (-160e-6, 40e-6), (160e-6, 500e-6));
    m.sheet("pad2", pad2_2d, 0.0);
    let jj_2d = m.rect("jj_2d", (-40e-6, -40e-6), (40e-6, 40e-6));
    m.sheet("jj_sheet", jj_2d, 0.0);
    m.overlay(&[("si", si), ("air", air)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(220e-6)).expect("3D transmon mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    // 1. FR-21 (lossless: true on aluminium_sc_emft omitting thickness & film_model) + FR-18 (squid_asymmetry) + FR-23 (TlsMeta)
    let mut spec = ProblemSpec::new("ems", 1, "cavity");
    spec.study.params = json!({ "n_modes": 1, "f_search": 5.5e9 });
    spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for pad in ["pad1", "pad2"] {
        spec.bcs.push(BcSpec {
            on: pad.into(),
            kind: "surface_impedance".into(),
            // Deliberately omit `thickness` and `film_model` to test fallback to nominal_thickness (150 nm) and film_model ("kerr")
            params: json!({
                "material": "aluminium_sc_emft",
                "temperature": 0.01,
                "lossless": true
            }),
        });
    }
    spec.bcs.push(BcSpec {
        on: "jj_sheet".into(),
        kind: "lumped_jj".into(),
        params: json!({
            "name": "Q1_SQUID",
            "E_j_Hz": 18.0e9,
            "C_j": 35.0e-15,
            "phi_ext_rad": 0.4,
            "squid_asymmetry": 0.25,
            "direction": [0.0, 1.0, 0.0]
        }),
    });

    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib.clone());
    let res_ll = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
    let ev_ll = res_ll.solutions[0].eigenvalue.unwrap();
    assert_eq!(ev_ll[1], 0.0, "lossless: true must take the real eigensolve path (Im f == 0)");
    let eq_ll = res_ll.scalar("equipartition_error").unwrap();
    assert!(eq_ll < 1e-11, "lossless aluminium_sc_emft equipartition_error = {eq_ll:e}");

    let rom = CircuitRom::of(&res_ll).expect("CircuitRom emitted");
    assert!(rom.modes[0].u_kin_ibc_j > 0.0, "kinetic inductance from aluminium_sc_emft must be active");
    assert_eq!(rom.ports[0].squid_asymmetry, Some(0.25), "FR-18 squid_asymmetry copied into RomPort");

    // FR-23 checks: jj_labels excluded from sheet_labels, sheet_thickness_m populated (150 nm), interfaces.ma populated
    let tls = TlsMeta::of(&res_ll).expect("TlsMeta emitted");
    let jj_lbls = g.boundary("jj_sheet");
    for &jl in &jj_lbls {
        assert!(!tls.sheet_labels.contains(&jl), "FR-23: jj_sheet label {jl} must NOT be in TlsMeta::sheet_labels");
    }
    let pad1_lbls = g.boundary("pad1");
    for &pl in &pad1_lbls {
        assert!(tls.sheet_labels.contains(&pl), "pad1 label {pl} must be in TlsMeta::sheet_labels");
        let thick = tls.sheet_thickness_m.get(&pl).copied().expect("nominal_thickness populated in sheet_thickness_m");
        assert!((thick - 150e-9).abs() < 1e-12, "expected 150 nm thickness, got {thick:e}");
    }
    let ifaces = tls.interfaces.as_ref().expect("FR-23: TlsMeta::interfaces populated from aluminium_sc_emft");
    let ma = ifaces.ma.as_ref().expect("aluminium_sc_emft defines [interfaces.ma]");
    assert!((ma.thickness_m - 3.0e-9).abs() < 1e-15);
    assert!((ma.eps_r - 9.8).abs() < 1e-12);

    // FR-18 validation: squid_asymmetry = 1.0 must be rejected with Error::Invalid
    let mut spec_bad_d = spec.clone();
    spec_bad_d.bcs.last_mut().unwrap().params["squid_asymmetry"] = json!(1.0);
    let mut ctx_bad = Ctx::new(&resolver);
    ctx_bad.library = Some(lib.clone());
    let err_d = Ems.solve(&exchange, &spec_bad_d, &mut ctx_bad).unwrap_err();
    assert!(matches!(err_d, Error::Invalid(_)), "squid_asymmetry=1.0 must return Error::Invalid, got {err_d:?}");

    // 2. FR-22: strong damping (Q ~ 50) via parallel lumped_port damping or lossy CPW still achieves equipartition_error < 1e-12 with Re(k0^2)
    let g_cpw = build_cpw_3d_mesh(49.3e-6, 28.95e-6, 297.0e-6, 8.080e-3, 500e-6, 500e-6, 1.2e-3);
    let ex_cpw = Mesh::from_solver_mesh(&g_cpw.mesh).unwrap();
    let resv_cpw = GeomResolver(&g_cpw);
    let mut spec_q50 = ProblemSpec::new("ems", 1, "cavity");
    spec_q50.study.params = json!({ "n_modes": 1, "f_search": 6.9e9 });
    spec_q50.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec_q50.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
        spec_q50.bcs.push(BcSpec {
            on: sheet.into(),
            kind: "surface_impedance".into(),
            // Large R_s on the CPW resonator tuned to give Q ~ 50
            params: json!({ "L_k_sq": 1.0e-13, "R_s": 2.0e-3 }),
        });
    }
    let res_q50 = Ems.solve(&ex_cpw, &spec_q50, &mut Ctx::new(&resv_cpw)).unwrap();
    let q_val = res_q50.scalar("Q1").expect("lossy mode has finite Q1");
    let eq_q50 = res_q50.scalar("equipartition_error").expect("equipartition_error emitted");
    eprintln!("FR-22 strongly damped CPW cavity: Q = {:.2}, equipartition_error = {:.3e}", q_val, eq_q50);
    assert!(q_val > 10.0 && q_val < 150.0, "expected low-Q (~50) regime, got Q = {q_val:.2}");
    assert!(eq_q50 < 1e-12, "FR-22 equipartition_error at Q={q_val:.1} was {eq_q50:e} (>= 1e-12)");
}

#[test]
fn fr26_refined_rank1_lumped_port_and_floating_jj_order1_and_order2() {
    // Both a grounded lumped_port AND a floating rank-1 lumped_jj sheet are refined so they have
    // interior and lateral-vacuum vertices; both order = 1 and order = 2 must converge cleanly
    // with 0 spurious near-zero modes and pass the post-assembly G_red^T b_red = 0 guard.
    let mut m = Model::new_3d();
    let si = m.box_("si_box", (-500e-6, -900e-6, -240e-6), (500e-6, 900e-6, 0.0));
    let air = m.box_("air_box", (-500e-6, -900e-6, 0.0), (500e-6, 900e-6, 240e-6));

    let feed_r = m.rect("feed_r", (-360e-6, -650e-6), (-220e-6, 650e-6));
    m.sheet("feed", feed_r, 0.0);
    let port_r = m.rect("port_r", (-360e-6, -900e-6), (-220e-6, -650e-6));
    m.sheet("port_in", port_r, 0.0);

    let pad_a = m.rect("pad_a_r", (60e-6, -480e-6), (340e-6, -70e-6));
    m.sheet("pad_a", pad_a, 0.0);
    let pad_b = m.rect("pad_b_r", (60e-6, 70e-6), (340e-6, 480e-6));
    m.sheet("pad_b", pad_b, 0.0);
    let jj_r = m.rect("jj_r", (130e-6, -70e-6), (270e-6, 70e-6));
    m.sheet("jj_sheet", jj_r, 0.0);

    m.overlay(&[("si", si), ("air", air)]);
    m.size("sz_port", geom::select::Selector::parse("port_in:edge").unwrap(), Some(50e-6.into()), None, None);
    m.size("sz_jj", geom::select::Selector::parse("jj_sheet:edge").unwrap(), Some(45e-6.into()), None, None);
    let plan = MeshPlan::simple(300e-6);
    let g = mesh_model(&m.list, &plan).expect("FR-26 refined port+JJ mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    for order in [1u8, 2u8] {
        let mut spec = ProblemSpec::new("ems", order, "cavity");
        spec.study.params = json!({ "n_modes": 2, "f_search": 6.0e9 });
        spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for s in ["feed", "pad_a", "pad_b"] {
            spec.bcs.push(BcSpec { on: s.into(), kind: "pec".into(), params: json!({}) });
        }
        spec.bcs.push(BcSpec {
            on: "port_in".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
        });
        spec.bcs.push(BcSpec {
            on: "jj_sheet".into(),
            kind: "lumped_jj".into(),
            params: json!({
                "name": "Q1_JJ",
                "E_j_Hz": 16.0e9,
                "C_j": 30.0e-15,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }),
        });
        let res = Ems.solve(&exchange, &spec, &mut Ctx::new(&resolver)).expect("FR-26 solve must converge");
        assert!(!res.solutions.is_empty());
        for sol in &res.solutions {
            let f_re = sol.eigenvalue.expect("eigenvalue present")[0];
            assert!(f_re > 1.0e8, "FR-26 order {order}: spurious near-zero mode {f_re:e} Hz");
        }
        let eq = res.scalar("equipartition_error").expect("equipartition_error emitted");
        assert!(eq < 1e-9, "FR-26 order {order}: equipartition_error = {eq:e} >= 1e-9");
    }
}

#[test]
fn fr29_assemble_maxwell_bitwise_deterministic_across_chunks() {
    // Build a > 3000-cell 3D mesh spanning multiple FIXED_CELL_CHUNK=2048 blocks and verify
    // assemble_maxwell produces bit-identical K and M values across repeated runs.
    let mesh = femcore::testkit::block(0.03, 0.02, 0.015, 9, 9, 9);
    assert!(mesh.n_cells() > 3000, "expected > 3000 cells, got {}", mesh.n_cells());
    let ned = nyx_ems::cavity::Ned3::new(&mesh, 1, (0..mesh.n_points() as u32).collect());
    let mut e = [[femcore::linalg::csr::c(0.0, 0.0); 3]; 3];
    for (a, row) in e.iter_mut().enumerate() {
        row[a] = femcore::linalg::csr::c(11.45, -1e-4);
    }
    let eps = vec![e; mesh.n_cells()];
    let mut mu_e = [[femcore::linalg::csr::c(0.0, 0.0); 3]; 3];
    for (a, row) in mu_e.iter_mut().enumerate() {
        row[a] = femcore::linalg::csr::c(1.0, 0.0);
    }
    let mu_inv = vec![mu_e; mesh.n_cells()];
    let (k1, m1) = nyx_ems::cavity::assemble_maxwell(&mesh, &ned, &eps, &mu_inv);
    let (k2, m2) = nyx_ems::cavity::assemble_maxwell(&mesh, &ned, &eps, &mu_inv);
    assert_eq!(k1.indptr(), k2.indptr());
    assert_eq!(k1.indices(), k2.indices());
    assert_eq!(
        k1.values().iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect::<Vec<_>>(),
        k2.values().iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect::<Vec<_>>()
    );
    assert_eq!(
        m1.values().iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect::<Vec<_>>(),
        m2.values().iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect::<Vec<_>>()
    );
}

#[test]
fn fr30_perturbative_ibc_parity_with_direct_and_coarse_p1_stability() {
    // Part 1: Parity on a fine P2 3D cavity wall (single-sided exterior face) with small kinetic inductance & R_s
    let boxm = femcore::testkit::block(0.024, 0.018, 0.012, 6, 6, 3);
    let box_ex = Mesh::from_solver_mesh(&boxm).unwrap();
    struct BoxRes;
    impl LabelResolver for BoxRes {
        fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
            match r {
                "air" => Ok(vec![1]),
                _ => Err(Error::invalid(format!("no region {r}"))),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            match s {
                "y0" => Ok(vec![2]),
                _ => Err(Error::invalid(format!("no boundary {s}"))),
            }
        }
    }
    let bres = BoxRes;
    let run_box = |mode: Option<&str>| -> nyx_physics::SolveResult {
        let mut s = ProblemSpec::new("ems", 2, "cavity");
        s.study.params = json!({ "n_modes": 1, "f_search": 10.0e9 });
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        if let Some(m) = mode {
            s.bcs.push(BcSpec {
                on: "y0".into(),
                kind: "surface_impedance".into(),
                params: json!({ "L_k_sq": 1.0e-10, "R_s": 5.0e-3, "f_ref": 10.417e9, "ibc_mode": m }),
            });
        }
        Ems.solve(&box_ex, &s, &mut Ctx::new(&bres)).unwrap()
    };
    let r_pec = run_box(None);
    let r_dir = run_box(Some("direct"));
    let r_pert = run_box(Some("perturbative"));
    let f_pec = r_pec.solutions[0].value.unwrap();
    let f_dir = r_dir.solutions[0].value.unwrap();
    let f_pert = r_pert.solutions[0].value.unwrap();
    let df_dir = f_dir - f_pec;
    let df_pert = f_pert - f_pec;
    let rel_shift_err = (df_pert - df_dir).abs() / df_dir.abs().max(1e-30);
    let q_dir = r_dir.scalar("Q1").expect("direct Q1");
    let q_pert = r_pert.scalar("Q1").expect("perturbative Q1");
    let rel_q_err = (q_pert - q_dir).abs() / q_dir.max(1e-30);
    eprintln!(
        "FR-30 P2 exterior cavity wall: df_dir = {:.4e} Hz, df_pert = {:.4e} Hz (rel err = {:.3}%), Q_dir = {:.3e}, Q_pert = {:.3e} (rel err = {:.3}%)",
        df_dir,
        df_pert,
        rel_shift_err * 100.0,
        q_dir,
        q_pert,
        rel_q_err * 100.0
    );
    assert!(df_dir < 0.0 && df_pert < 0.0, "kinetic inductance must lower resonant frequency");
    assert!(
        rel_shift_err < 0.05,
        "FR-30 exterior wall perturbative vs direct df err = {:.3}% (>= 5%)",
        rel_shift_err * 100.0
    );
    assert!(rel_q_err < 0.05, "FR-30 exterior wall perturbative vs direct Q err = {:.3}% (>= 5%)", rel_q_err * 100.0);

    // Check CircuitRom matches perturbed eigenvalue & Q_cond
    let rom_pert = CircuitRom::of(&r_pert).expect("CircuitRom on perturbative IBC");
    let rom_f = rom_eigenfrequencies(&rom_pert).unwrap();
    let rom_fc = rom_complex_eigenfrequencies(&rom_pert).unwrap();
    assert!((rom_f[0] - f_pert).abs() / f_pert < 1e-5, "CircuitRom f_re {} vs {}", rom_f[0], f_pert);
    let q_rom = rom_fc[0].re / (2.0 * rom_fc[0].im.abs());
    assert!((q_rom - q_pert).abs() / q_pert < 1e-4, "CircuitRom Q {} vs {}", q_rom, q_pert);

    // Part 2: Parity on a fine-P2 interior 2-sided microstrip/parallel-plate resonator (h << w, tests exact [n x H] jump!)
    let mut mm = Model::new_3d();
    let si = mm.box_("si", (-1200e-6, 0.0, -300e-6), (1200e-6, 3600e-6, 0.0));
    let air = mm.box_("air", (-1200e-6, 0.0, 0.0), (1200e-6, 3600e-6, 300e-6));
    let strip = mm.rect("strip_r", (-1200e-6, 0.0), (1200e-6, 3600e-6));
    mm.sheet("strip", strip, 0.0);
    mm.overlay(&[("si", si), ("air", air)]);
    let mut ms_plan = MeshPlan::simple(700e-6);
    ms_plan.opts.dz_min = Some(100e-6);
    ms_plan.opts.max_z_slabs = Some(6);
    let g_ms = mesh_model(&mm.list, &ms_plan).expect("fine P2 2-sided interior sheet mesh");
    let ex_ms = Mesh::from_solver_mesh(&g_ms.mesh).unwrap();
    let res_ms = GeomResolver(&g_ms);
    let run_ms = |mode: Option<&str>| -> nyx_physics::SolveResult {
        let mut s = ProblemSpec::new("ems", 2, "cavity");
        s.study.params = json!({ "n_modes": 1, "f_search": 12.5e9 });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        match mode {
            None => s.bcs.push(BcSpec { on: "strip".into(), kind: "pec".into(), params: json!({}) }),
            Some(m) => s.bcs.push(BcSpec {
                on: "strip".into(),
                kind: "surface_impedance".into(),
                params: json!({ "L_k_sq": 1.0e-12, "R_s": 2.0e-4, "f_ref": 12.5e9, "ibc_mode": m }),
            }),
        }
        Ems.solve(&ex_ms, &s, &mut Ctx::new(&res_ms)).unwrap()
    };
    let ms_pec = run_ms(None).solutions[0].value.unwrap();
    let ms_dir_res = run_ms(Some("direct"));
    let ms_pert_res = run_ms(Some("perturbative"));
    let ms_dir = ms_dir_res.solutions[0].value.unwrap();
    let ms_pert = ms_pert_res.solutions[0].value.unwrap();
    let ms_df_dir = ms_dir - ms_pec;
    let ms_df_pert = ms_pert - ms_pec;
    let ms_rel = (ms_df_pert - ms_df_dir).abs() / ms_df_dir.abs().max(1e-30);
    eprintln!(
        "FR-30 P2 interior 2-sided sheet: f_pec = {:.6} GHz, df_dir = {:.4e} Hz, df_pert = {:.4e} Hz (rel err = {:.2}%)",
        ms_pec / 1e9,
        ms_df_dir,
        ms_df_pert,
        ms_rel * 100.0
    );
    assert!(ms_df_dir < 0.0 && ms_df_pert < 0.0);
    assert!(ms_rel < 0.05, "FR-30 interior 2-sided sheet perturbative vs direct err = {:.2}% (>= 5%)", ms_rel * 100.0);

    // Part 3: Coarse-P1 narrow strip stability (`auto` -> `perturbative`) + `direct` coarse warning check
    let g_coarse = build_cpw_3d_mesh(20e-6, 12e-6, 400e-6, 6.0e-3, 400e-6, 400e-6, 800e-6);
    let ex_c = Mesh::from_solver_mesh(&g_coarse.mesh).unwrap();
    let res_c = GeomResolver(&g_coarse);
    let run_coarse = |mode: Option<&str>| -> nyx_physics::SolveResult {
        let mut s = ProblemSpec::new("ems", 1, "cavity");
        s.study.params = json!({ "n_modes": 1, "f_search": 9.5e9 });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sh in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
            match mode {
                None => s.bcs.push(BcSpec { on: sh.into(), kind: "pec".into(), params: json!({}) }),
                Some(m) => s.bcs.push(BcSpec {
                    on: sh.into(),
                    kind: "surface_impedance".into(),
                    params: json!({ "L_k_sq": 1.0e-13, "R_s": 0.0, "ibc_mode": m }),
                }),
            }
        }
        Ems.solve(&ex_c, &s, &mut Ctx::new(&res_c)).unwrap()
    };
    let c_pec = run_coarse(None).solutions[0].value.unwrap();
    let c_dir_res = run_coarse(Some("direct"));
    let c_auto_res = run_coarse(Some("auto"));
    let c_auto = c_auto_res.solutions[0].value.unwrap();
    assert!(
        c_dir_res.log.iter().any(|l| l.contains("coarse strip resolution")),
        "FR-30 coarse P1 direct mode must emit coarse strip resolution warning; log = {:?}",
        c_dir_res.log
    );
    assert!(
        c_auto_res.log.iter().any(|l| l.contains("perturbative mode active")),
        "FR-30 coarse P1 auto mode must select perturbative; log = {:?}",
        c_auto_res.log
    );
    let auto_shift_pct = (c_auto - c_pec).abs() / c_pec * 100.0;
    assert!(
        c_auto < c_pec && auto_shift_pct < 2.0,
        "FR-30 coarse P1 auto/perturbative shift must be small and physical (< 2%), got {:.3}%",
        auto_shift_pct
    );
}
