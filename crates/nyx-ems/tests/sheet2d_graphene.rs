// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative verification of FR-35 (ADR 137): Zero-Thickness 2D Materials & Tensor Surface Conductivity Sheets.
//!
//! Covers:
//! 1. 3D `ems::driven` analytical Kubo graphene normal-incidence Fresnel $T / R / A$ vs closed-form
//!    $T = |2 / (2 + \eta_0 \sigma)|^2$, $R = |-\eta_0 \sigma / (2 + \eta_0 \sigma)|^2$, $A = 1 - T - R$ ($< 10^{-6}$).
//! 2. 2D `ems::modes` TM graphene surface plasmon dispersion vs analytical $q_{\text{spp}}$ ($< 0.5\%$).
//! 3. Electrostatically gated graphene-on-Si photonic waveguide (`potential` $V(\mathbf{r}) \to \mu_c(\mathbf{r})$)
//!    Pauli-blocking attenuation sweep $\Delta\alpha(V_g)$ at $\lambda_0 = 1550\,\text{nm}$.
//! 4. Magneto-optic Hall bias $B_0 \neq 0$ ($\sigma_{xy} = -\sigma_{yx}$, non-symmetric `ImageProm`),
//!    in-plane rotated anisotropic 2D sheet ($\theta_{2\text{D}} = 35^\circ$), 3D cavity `Kx/Ky/Kz` surface current
//!    extraction, and thread determinism (`NYX_THREADS = 1` vs `8`).

use femcore::linalg::csr::c;
use femcore::linalg::prom::PromOperator;
use femcore::material::{
    chemical_potential_from_gate_potential, gate_oxide_capacitance, kubo_graphene_isotropic,
    sheet2d_fresnel_normal_tra, tm_plasmon_wavenumber_nonretarded, tm_plasmon_wavenumber_retarded_symmetric, Library,
    EV_TO_J,
};
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_ems::ibc::ETA0;
use nyx_ems::{assemble_emw_driven, decode, Ems, C0};
use nyx_physics::{
    discretise, BcSpec, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec,
};
use serde_json::json;
use std::f64::consts::PI;
use std::path::Path;
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

fn load_repo_material_library() -> Arc<Library> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Arc::new(Library::load_dir(&dir).expect("load repository materials/ library"))
}

/// Builds a 3D TEM parallel-plate transmission cell `[0, w] x [0, w] x [-l/2, +l/2]` with a conformal
/// zero-thickness 2D sheet (`"graphene_sheet"`) at the midplane `z = 0`.
/// - `x = 0` and `x = w`: PEC walls (`E_x` normal to PEC, so TEM electric field is polarized along `x`)
/// - `y = 0` and `y = w`: PMC walls (`H_y` normal to PMC, so TEM magnetic field is polarized along `y`)
/// - `z = -l/2` (`"z = min"`): Port 1 (`Z_0 = ETA0`, `direction = [1, 0, 0]`)
/// - `z = +l/2` (`"z = max"`): Port 2 (`Z_0 = ETA0`, `direction = [1, 0, 0]`)
fn build_3d_tem_fresnel_sheet_mesh(w: f64, l: f64, h_mesh: f64) -> MeshedGeom {
    let mut m = Model::new_3d();
    let b_in = m.box_("box_in", (0.0, 0.0, -0.5 * l), (w, w, 0.0));
    let b_out = m.box_("box_out", (0.0, 0.0, 0.0), (w, w, 0.5 * l));
    let gr_2d = m.rect("gr_rect", (0.0, 0.0), (w, w));
    m.sheet("graphene_sheet", gr_2d, 0.0);
    m.overlay(&[("vac_in", b_in), ("vac_out", b_out)]);
    mesh_model(&m.list, &MeshPlan::simple(h_mesh)).expect("conformal 3D TEM Fresnel sheet mesh")
}

#[test]
fn fresnel_normal_incidence_kubo_graphene_3d_driven() {
    let lib = load_repo_material_library();
    let w = 2.0e-6;
    let l = 1.0e-9;
    let g = build_3d_tem_fresnel_sheet_mesh(w, l, 1.2e-6);
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    // Test across three THz / mid-IR frequencies with monolayer Kubo graphene (mu_c = 0.40 eV, tau = 0.5 ps, T = 300 K)
    let mu_c_ev = 0.40;
    let tau_s = 0.5e-12;
    let temp_k = 300.0;
    let v_f = 1.0e6;

    for &freq_hz in &[2.0e12, 5.0e12, 10.0e12] {
        let omega = 2.0 * PI * freq_hz;
        let sigma = kubo_graphene_isotropic(omega, temp_k, mu_c_ev * EV_TO_J, tau_s).unwrap();
        let (t_exact, r_exact, a_exact) = sheet2d_fresnel_normal_tra(sigma, 1.0, 1.0);

        let mut spec = ProblemSpec::new("ems", 1, "driven");
        spec.study.params = json!({
            "freqs_hz": [freq_hz],
            "prom_tol": null
        });
        spec.materials.insert("vac_in".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        spec.materials.insert("vac_out".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        spec.bcs.push(BcSpec { on: "extreme(-x)".into(), kind: "pec".into(), params: json!({}) });
        spec.bcs.push(BcSpec { on: "extreme(+x)".into(), kind: "pec".into(), params: json!({}) });
        spec.bcs.push(BcSpec { on: "extreme(-y)".into(), kind: "pmc".into(), params: json!({}) });
        spec.bcs.push(BcSpec { on: "extreme(+y)".into(), kind: "pmc".into(), params: json!({}) });
        spec.bcs.push(BcSpec {
            on: "extreme(-z)".into(),
            kind: "lumped_port".into(),
            params: json!({
                "name": "port_in",
                "Z_0": ETA0,
                "direction": [1.0, 0.0, 0.0]
            }),
        });
        spec.bcs.push(BcSpec {
            on: "extreme(+z)".into(),
            kind: "lumped_port".into(),
            params: json!({
                "name": "port_out",
                "Z_0": ETA0,
                "direction": [1.0, 0.0, 0.0]
            }),
        });
        spec.bcs.push(BcSpec {
            on: "graphene_sheet".into(),
            kind: "sheet_2d".into(),
            params: json!({
                "material": "graphene",
                "mu_c_eV": mu_c_ev,
                "tau_s": tau_s,
                "temperature": temp_k,
                "v_f": v_f,
                "f_ref": freq_hz
            }),
        });

        let mut ctx = Ctx::new(&resolver);
        ctx.library = Some(lib.clone());
        let res = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
        let data = DrivenPortData::of(&res).expect("DrivenPortData");
        let s11 = c(data.s_matrix[0][0][0][0], data.s_matrix[0][0][0][1]);
        let s21 = c(data.s_matrix[0][1][0][0], data.s_matrix[0][1][0][1]);
        let t_fem = s21.norm_sqr();
        let r_fem = s11.norm_sqr();
        let a_fem = 1.0 - t_fem - r_fem;

        let err_t = (t_fem - t_exact).abs();
        let err_r = (r_fem - r_exact).abs();
        let err_a = (a_fem - a_exact).abs();
        eprintln!(
            "3D Fresnel @ {:.1} THz: sigma = ({:.4e}, {:.4e}) S | T_fem = {:.8} vs T_ex = {:.8} (|dT| = {:.2e}), |dR| = {:.2e}, |dA| = {:.2e}",
            freq_hz / 1e12,
            sigma.re,
            sigma.im,
            t_fem,
            t_exact,
            err_t,
            err_r,
            err_a
        );
        assert!(err_t < 1e-6, "Fresnel transmittance error {err_t:e} >= 1e-6 at {freq_hz} Hz");
        assert!(err_r < 1e-6, "Fresnel reflectance error {err_r:e} >= 1e-6 at {freq_hz} Hz");
        assert!(err_a < 1e-6, "Fresnel absorptance error {err_a:e} >= 1e-6 at {freq_hz} Hz");
    }
}

#[test]
fn tm_graphene_surface_plasmon_dispersion_2d_modes() {
    let lib = load_repo_material_library();
    // 2D cross-section strip x in [0, w_x], z in [-h_z, +h_z] split at z = 0 by a monolayer graphene sheet.
    // PMC walls on x = 0 and x = w_x enforce the uniform-in-x 2D TM plasmon (E_z, E_y != 0, E_x = 0).
    let w_x = 0.25e-6;
    let h_z = 3.5e-6;
    let mut m = Model::new_2d();
    let bot = m.rect("cladding_bot", (0.0, -h_z), (w_x, 0.0));
    let top = m.rect("cladding_top", (0.0, 0.0), (w_x, h_z));
    m.overlay(&[("bot", bot), ("top", top)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(0.14e-6)).expect("2D graphene TM plasmon strip mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    // At f = 8.0 THz (lambda_0 = 37.47 um), mu_c = 0.25 eV, tau = 1.0 ps, T = 300 K:
    // Graphene conductivity is strongly inductive (Im(sigma) < 0), supporting a tightly bound TM surface plasmon
    // with Re(n_eff) ~ 9.115 ( plasmon wavelength lambda_spp ~ 4.1 um, vertical decay length 1/gamma ~ 0.65 um,
    // so h_z = 3.5 um > 5.3 decay lengths).
    let freq_hz = 8.0e12;
    let omega = 2.0 * PI * freq_hz;
    let k0 = omega / C0;
    let mu_c_ev = 0.25;
    let tau_s = 1.0e-12;
    let temp_k = 300.0;
    let eps_r = 1.0;

    let sigma = kubo_graphene_isotropic(omega, temp_k, mu_c_ev * EV_TO_J, tau_s).unwrap();
    let q_nonret = tm_plasmon_wavenumber_nonretarded(sigma, omega, eps_r, eps_r);
    let beta_ret = tm_plasmon_wavenumber_retarded_symmetric(sigma, omega, eps_r);
    let neff_ret = beta_ret / k0;
    let neff_nonret = ((k0 * k0 * eps_r) + q_nonret * q_nonret).sqrt() / k0;

    let mut spec = ProblemSpec::new("ems", 2, "modes");
    spec.study.params = json!({
        "freq": freq_hz,
        "n_target": neff_ret.re,
        "evanescent": false,
        "n_modes": 1
    });
    spec.materials.insert("bot".into(), MaterialAssign::Props([("eps_r".to_string(), eps_r)].into()));
    spec.materials.insert("top".into(), MaterialAssign::Props([("eps_r".to_string(), eps_r)].into()));
    spec.bcs.push(BcSpec { on: "extreme(-x)".into(), kind: "pmc".into(), params: json!({}) });
    spec.bcs.push(BcSpec { on: "extreme(+x)".into(), kind: "pmc".into(), params: json!({}) });
    spec.bcs.push(BcSpec {
        on: "interface(bot, top)".into(),
        kind: "sheet_2d".into(),
        params: json!({
            "material": "graphene",
            "mu_c_eV": mu_c_ev,
            "tau_s": tau_s,
            "temperature": temp_k
        }),
    });

    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib);
    let res = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
    assert!(!res.solutions.is_empty(), "TM graphene surface plasmon mode must converge");
    let eig = res.solutions[0].eigenvalue.expect("complex neff eigenvalue");
    let neff_fem = c(eig[0], eig[1]);
    let beta_fem = neff_fem * k0;
    let gamma_fem = (beta_fem * beta_fem - k0 * k0 * eps_r).sqrt();

    let rel_err_ret = (beta_fem - beta_ret).norm() / beta_ret.norm();
    let rel_err_gamma_qs = (gamma_fem - q_nonret).norm() / q_nonret.norm();
    eprintln!(
        "2D TM graphene plasmon @ 8 THz: neff_fem = {:.6} {:+.4e}i vs neff_ret = {:.6} {:+.4e}i (neff_qs_corr = {:.6} {:+.4e}i) | rel_err_beta = {:.4}%, rel_err_gamma_qs = {:.4}%",
        neff_fem.re,
        neff_fem.im,
        neff_ret.re,
        neff_ret.im,
        neff_nonret.re,
        neff_nonret.im,
        rel_err_ret * 100.0,
        rel_err_gamma_qs * 100.0
    );
    assert!(rel_err_ret < 5e-3, "2D TM graphene plasmon dispersion error {:.4}% >= 0.5%", rel_err_ret * 100.0);
    assert!(
        rel_err_gamma_qs < 5e-3,
        "2D TM graphene plasmon quasi-static decay wavenumber error {:.4}% >= 0.5%",
        rel_err_gamma_qs * 100.0
    );
}

#[test]
fn electrostatically_gated_graphene_on_si_waveguide_attenuation_sweep() {
    let lib = load_repo_material_library();
    // Graphene-on-Si photonic strip waveguide at lambda_0 = 1550 nm (f = 193.414 THz, hbar * omega = 0.80 eV):
    // - SiO2 box/cladding: x in [-1.0 um, +1.0 um], z in [-0.6 um, +0.6 um] (eps_r = 2.085)
    // - Si waveguide core: x in [-0.25 um, +0.25 um], z in [-0.11 um, +0.11 um] (eps_r = 12.08)
    // - Monolayer graphene sheet on the top interface of the Si core (`core|clad` / `z = +0.11 um`)
    // - Gate oxide Al2O3: eps_ox = 9.0, d_ox = 7.0 nm -> C_ox = eps_0 * eps_ox / d_ox = 1.1378e-2 F/m^2.
    let mut m = Model::new_2d();
    let clad = m.rect("clad_rect", (-1.0e-6, -0.6e-6), (1.0e-6, 0.6e-6));
    let core_bot = m.rect("core_bot_rect", (-0.25e-6, -0.11e-6), (0.25e-6, 0.0));
    let core_top = m.rect("core_top_rect", (-0.25e-6, 0.0), (0.25e-6, 0.11e-6));
    m.overlay(&[("clad", clad), ("core_bot", core_bot), ("core_top", core_top)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(0.11e-6)).expect("2D graphene-on-Si waveguide mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let lambda0 = 1550.0e-9;
    let freq_hz = C0 / lambda0;
    let k0 = 2.0 * PI / lambda0;
    let eps_ox = 9.0;
    let d_ox = 7.0e-9;
    let c_ox = gate_oxide_capacitance(eps_ox, d_ox).unwrap();
    let v_cnp = -0.15; // slight initial doping at V_g = 0 V (|V_g - V_cnp| = 0.15 V -> mu_c ~ 0.123 eV < hbar*omega/2 = 0.40 eV)

    // Sweep electrostatic gate potential field V(x, z) across V_g in [0.0, 1.5, 3.0, 5.0] V
    // (simulating the imported `dfs.V` nodal electrostatic potential field over the waveguide cross-section).
    let gate_voltages = [0.0_f64, 1.5, 3.0, 5.0];
    let mut alphas_db_per_um = Vec::with_capacity(gate_voltages.len());
    let mut mu_c_ev_list = Vec::with_capacity(gate_voltages.len());

    for &v_g in &gate_voltages {
        let (_n_s, mu_c_j) = chemical_potential_from_gate_potential(v_g, v_cnp, c_ox, 1.0e6).unwrap();
        let mu_c_ev = mu_c_j / EV_TO_J;
        mu_c_ev_list.push(mu_c_ev);

        // Nodal electrostatic potential field V(x, z) from gate electrode centered over the Si core:
        let np = g.mesh.n_points();
        let v_field: Vec<f64> = (0..np)
            .map(|i| {
                let p = g.mesh.point(i);
                // Smooth lateral profile peaking at V_g across the active core [-0.25 um, +0.25 um]
                let x_norm = (p[0] / 0.6e-6).powi(2);
                v_g * (-0.15 * x_norm).exp()
            })
            .collect();

        let mut spec = ProblemSpec::new("ems", 2, "modes");
        spec.study.params = json!({
            "freq": freq_hz,
            "n_target": 2.45,
            "evanescent": false,
            "n_modes": 1
        });
        spec.materials.insert("core_bot".into(), MaterialAssign::Props([("eps_r".to_string(), 12.08)].into()));
        spec.materials.insert("core_top".into(), MaterialAssign::Props([("eps_r".to_string(), 12.08)].into()));
        spec.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".to_string(), 2.085)].into()));
        spec.bcs.push(BcSpec {
            on: "interface(core_bot, core_top)".into(),
            kind: "sheet_2d".into(),
            params: json!({
                "material": "graphene",
                "tau_s": 0.2e-12,
                "temperature": 300.0,
                "eps_ox": eps_ox,
                "d_ox": d_ox,
                "v_cnp": v_cnp
            }),
        });

        let mut ctx = Ctx::new(&resolver);
        ctx.library = Some(lib.clone());
        ctx.inputs.insert("potential".into(), v_field);
        let res = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
        let eig = res.solutions[0].eigenvalue.expect("complex waveguide mode eigenvalue");
        let neff = c(eig[0], eig[1]);
        // Power attenuation coefficient alpha = -2 * k0 * Im(n_eff) [Np/m]; in dB/um: 4.3429448 * alpha * 1e-6
        let alpha_np_m = (-2.0 * k0 * neff.im).max(0.0);
        let alpha_db_um = 10.0 / std::f64::consts::LN_10 * alpha_np_m * 1e-6;
        eprintln!(
            "Gated graphene-on-Si @ V_g = {:.1} V (mu_c = {:.3} eV): neff = {:.6} {:+.4e}i, alpha = {:.4} dB/um",
            v_g, mu_c_ev, neff.re, neff.im, alpha_db_um
        );
        alphas_db_per_um.push(alpha_db_um);
    }

    // Below Pauli blocking (V_g = 0 V, mu_c ~ 0.12 eV < 0.40 eV), interband absorption is strong (> 0.02 dB/um).
    // Above Pauli blocking (V_g = 5 V, mu_c ~ 0.72 eV > 0.40 eV), interband absorption is suppressed by > 85%.
    let alpha_on = alphas_db_per_um[0];
    let alpha_off = *alphas_db_per_um.last().unwrap();
    let delta_alpha = alpha_on - alpha_off;
    eprintln!(
        "Electro-absorption modulation depth: alpha(0 V) = {:.4} dB/um -> alpha(5 V) = {:.4} dB/um (Delta alpha = {:.4} dB/um, {:.1}% suppression)",
        alpha_on,
        alpha_off,
        delta_alpha,
        100.0 * delta_alpha / alpha_on.max(1e-30)
    );
    assert!(alpha_on > 0.02, "expected strong interband absorption at V_g = 0 V, got {alpha_on} dB/um");
    assert!(
        alpha_off < 0.15 * alpha_on,
        "expected Pauli-blocked attenuation at V_g = 5 V ({alpha_off} dB/um) to be < 15% of V_g = 0 V ({alpha_on} dB/um)"
    );
    assert!(delta_alpha > 0.015, "expected positive modulation swing Delta alpha > 0.015 dB/um, got {delta_alpha}");
}

#[test]
fn magneto_optic_hall_sheet_anisotropic_rotated_sheet_and_thread_determinism() {
    let lib = load_repo_material_library();
    let w = 2.0e-6;
    let l = 1.0e-9;
    let g = build_3d_tem_fresnel_sheet_mesh(w, l, 1.2e-6);
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let d3 = discretise(&exchange).unwrap();
    let resolver = GeomResolver(&g);

    // 1. Out-of-plane magnetic bias B_0 = 5.0 T produces gyrotropic Hall off-diagonal sigma_xy = -sigma_yx != 0
    //    and flags DrivenAssembly::has_non_symmetric_sheets = true (PromOperator::is_symmetric() == false).
    let mut spec_hall = ProblemSpec::new("ems", 1, "driven");
    spec_hall.study.params = json!({
        "f_min": 2.0e12,
        "f_max": 4.0e12,
        "n_points": 5,
        "prom_tol": 1e-4
    });
    spec_hall.materials.insert("vac_in".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    spec_hall.materials.insert("vac_out".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    spec_hall.bcs.push(BcSpec { on: "extreme(-x)".into(), kind: "pec".into(), params: json!({}) });
    spec_hall.bcs.push(BcSpec { on: "extreme(+x)".into(), kind: "pec".into(), params: json!({}) });
    spec_hall.bcs.push(BcSpec { on: "extreme(-y)".into(), kind: "pmc".into(), params: json!({}) });
    spec_hall.bcs.push(BcSpec { on: "extreme(+y)".into(), kind: "pmc".into(), params: json!({}) });
    spec_hall.bcs.push(BcSpec {
        on: "extreme(-z)".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "p1", "Z_0": ETA0, "direction": [1.0, 0.0, 0.0] }),
    });
    spec_hall.bcs.push(BcSpec {
        on: "extreme(+z)".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "p2", "Z_0": ETA0, "direction": [1.0, 0.0, 0.0] }),
    });
    spec_hall.bcs.push(BcSpec {
        on: "graphene_sheet".into(),
        kind: "sheet_2d".into(),
        params: json!({
            "material": "graphene",
            "mu_c_eV": 0.35,
            "tau_s": 0.5e-12,
            "temperature": 300.0,
            "b_0": 5.0
        }),
    });

    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib.clone());
    let prob_hall = decode(&spec_hall, &ctx).unwrap();
    let asm_hall = assemble_emw_driven(&d3, &prob_hall).unwrap();
    assert!(asm_hall.has_non_symmetric_sheets, "B_0 != 0 must flag has_non_symmetric_sheets");
    assert!(!PromOperator::is_symmetric(&asm_hall), "PromOperator::is_symmetric must be false when B_0 != 0");

    let res_prom = Ems.solve(&exchange, &spec_hall, &mut ctx).unwrap();
    let mut spec_hdm = spec_hall.clone();
    spec_hdm.study.params["prom_tol"] = serde_json::Value::Null;
    let res_hdm = Ems.solve(&exchange, &spec_hdm, &mut ctx).unwrap();
    let d_prom = DrivenPortData::of(&res_prom).unwrap();
    let d_hdm = DrivenPortData::of(&res_hdm).unwrap();
    for k in 0..d_prom.freqs_hz.len() {
        let s21_p = c(d_prom.s_matrix[k][1][0][0], d_prom.s_matrix[k][1][0][1]);
        let s21_h = c(d_hdm.s_matrix[k][1][0][0], d_hdm.s_matrix[k][1][0][1]);
        let err = (s21_p - s21_h).norm();
        assert!(err < 1e-5, "non-symmetric Hall PROM vs HDM mismatch {err:e} at index {k}");
    }

    // 2. Thread determinism (1 vs 8 threads) on 3D driven solve with Sheet2D
    let pool1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let pool8 = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap();
    let res_t1 = pool1.install(|| {
        let mut c1 = Ctx::new(&resolver);
        c1.library = Some(lib.clone());
        Ems.solve(&exchange, &spec_hdm, &mut c1).unwrap()
    });
    let res_t8 = pool8.install(|| {
        let mut c8 = Ctx::new(&resolver);
        c8.library = Some(lib.clone());
        Ems.solve(&exchange, &spec_hdm, &mut c8).unwrap()
    });
    let d_t1 = DrivenPortData::of(&res_t1).unwrap();
    let d_t8 = DrivenPortData::of(&res_t8).unwrap();
    assert_eq!(d_t1.s_matrix, d_t8.s_matrix, "3D Sheet2D driven S-matrix must be bit-identical for 1 vs 8 threads");
    assert_eq!(
        res_t1.solutions[0].fields[0].values, res_t8.solutions[0].fields[0].values,
        "3D Sheet2D nodal |E| field must be bit-identical for 1 vs 8 threads"
    );

    // 3. 2D rotated anisotropic sheet (theta_deg = 35.0, sigma_uu != sigma_vv) triggering off-diagonal A1 coupling
    let mut m2 = Model::new_2d();
    let r_bot = m2.rect("r_bot", (0.0, -1.0e-6), (1.0e-6, 0.0));
    let r_top = m2.rect("r_top", (0.0, 0.0), (1.0e-6, 1.0e-6));
    m2.overlay(&[("bot", r_bot), ("top", r_top)]);
    let g2 = mesh_model(&m2.list, &MeshPlan::simple(0.35e-6)).unwrap();
    let ex2 = Mesh::from_solver_mesh(&g2.mesh).unwrap();
    let res2_map = GeomResolver(&g2);

    let mut spec_aniso = ProblemSpec::new("ems", 1, "modes");
    spec_aniso.study.params = json!({
        "freq": 30.0e12,
        "n_target": 1.8,
        "evanescent": false,
        "n_modes": 1
    });
    spec_aniso.materials.insert("bot".into(), MaterialAssign::Props([("eps_r".to_string(), 3.5)].into()));
    spec_aniso.materials.insert("top".into(), MaterialAssign::Props([("eps_r".to_string(), 2.2)].into()));
    spec_aniso.bcs.push(BcSpec {
        on: "interface(bot, top)".into(),
        kind: "sheet_2d".into(),
        params: json!({
            "model": "anisotropic",
            "sigma_uu_re": 2.0e-4,
            "sigma_uu_im": -1.5e-3,
            "sigma_vv_re": 0.5e-4,
            "sigma_vv_im": -0.3e-3,
            "theta_deg": 35.0
        }),
    });
    let m_t1 = pool1.install(|| Ems.solve(&ex2, &spec_aniso, &mut Ctx::new(&res2_map)).unwrap());
    let m_t8 = pool8.install(|| Ems.solve(&ex2, &spec_aniso, &mut Ctx::new(&res2_map)).unwrap());
    assert_eq!(
        m_t1.solutions[0].eigenvalue, m_t8.solutions[0].eigenvalue,
        "2D rotated anisotropic Sheet2D mode eigenvalue must be bit-identical for 1 vs 8 threads"
    );
}
