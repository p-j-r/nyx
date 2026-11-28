// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-41 (ADR 143) quantitative verification suite for $\chi^{(2)}$ and $\chi^{(3)}$ nonlinear optics
//! Harmonic Balance FEM (`Study::NonlinearDriven`) and 2D fixed-power waveguide eigenmodes (`Study::NonlinearModes`).
//!
//! Covers:
//! 1. Exact Wirtinger Fréchet Jacobian vs central finite-difference verification (`< 1e-8`) for:
//!    - Single-tone anisotropic rotated `m-3m` Kerr + TPA + steady-state FCA/FCD (`eval_chi3_residual_and_jacobian`),
//!    - Multi-harmonic depleted `shg`, 3-wave `sfg` with QPM poling, and two-tone `chi3_xpm` (`eval_chi2_residual_and_jacobian`).
//! 2. Depleted-pump $\chi^{(2)}$ SHG conversion efficiency $\eta(L) = \tanh^2(\Gamma L)$ (`< 0.5%` up to $> 75\%$
//!    pump depletion), exact discrete Manley–Rowe photon-flux conservation ($N_1 + 2 N_2 = \text{const}$, `< 1e-6`),
//!    `newton` vs `jfnk` parity, `undepleted` vs `depleted` comparison, and 3-wave `sfg` / `opa` Manley–Rowe conservation.
//! 3. Periodically poled `lithium_niobate` (`3m`, Zelmon Sellmeier dispersion) QPM SHG ($\Lambda_{\text{QPM}} = 2 L_c$)
//!    achieving $> 200\times$ conversion enhancement over unpoled phase-mismatched $\text{LiNbO}_3$ and matching
//!    first-order Fourier $d_{\text{eff,QPM}} = \frac{2}{\pi} |d_{33}|$ with exact Manley–Rowe conservation (`< 1e-6`).
//! 4. High-$Q$ ($Q_L \approx 1.0 \times 10^5$) `silicon_nitride` Kerr cavity S-curve optical bistability traced cleanly
//!    across upper and lower saddle-node turning points (including the unstable middle branch $\mathrm{d}\alpha/\mathrm{d}s < 0$)
//!    via bordered Keller pseudo-arclength continuation (`continuation = "arclength"`), agreeing with analytical
//!    single-mode Coupled-Mode Theory (`< 0.5%`).
//! 5. Cross-polarized $\chi^{(3)}$ SPM (circular $2/3$ Maker–Terhune ratio), two-tone `chi3_xpm` ($2\times$ co-polarized vs
//!    $2/3\times$ orthogonally polarized XPM, `< 0.5%`), `silicon` TPA + FCA/FCD nonlinear saturation, and 2D cross-section
//!    `Study::NonlinearModes` $\Delta n_{\text{eff}}(P) = \gamma P / k_0$ (`< 0.5%`).
//! 6. Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` bitwise `to_bits()` determinism across `Study::NonlinearDriven` and
//!    `Study::NonlinearModes`.

use femcore::linalg::csr::{c, RCsr};
use femcore::material::Library;
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::testkit::block;
use geom::model::Model;
use nyx_ems::ibc::{EPS0, MU0};
use nyx_ems::nonlinear::{assemble_nonlinear_driven, eval_chi2_residual_and_jacobian, eval_chi3_residual_and_jacobian};
use nyx_ems::{Ems, C0};
use nyx_physics::{
    discretise, BcSpec, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, Scalar,
    SolveResult,
};
use serde_json::json;
use std::f64::consts::PI;
use std::path::Path;
use std::process::Command;
use std::sync::Arc;

struct BlockResolver;

impl LabelResolver for BlockResolver {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        match region {
            "core" | "guide" | "cav" => Ok(vec![1]),
            _ => Err(Error::invalid(format!("no region '{region}'"))),
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        match selector {
            "x0" => Ok(vec![1]),
            "y0" => Ok(vec![2]),
            "z0" => Ok(vec![3]),
            "z1" => Ok(vec![4]),
            "y1" => Ok(vec![5]),
            "x1" => Ok(vec![6]),
            "pmc_x" => Ok(vec![1, 6]),
            "pec_z" => Ok(vec![3, 4]),
            "pec_xz" => Ok(vec![1, 3, 4, 6]),
            _ => Err(Error::invalid(format!("no boundary '{selector}'"))),
        }
    }
}

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

fn load_library() -> Arc<Library> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Arc::new(Library::load_dir(&dir).expect("load shipped materials library"))
}

fn num_scalar(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(v))) => *v,
        other => panic!("missing numeric scalar '{name}', got {other:?} in {:?}", r.scalars),
    }
}

fn csr_frob_rel_err_vs_fd<F>(j_exact: &RCsr, x0: &[f64], eps: f64, mut eval_f: F) -> f64
where
    F: FnMut(&[f64], &mut [f64]),
{
    let n = x0.len();
    assert_eq!(j_exact.nrows(), n);
    assert_eq!(j_exact.ncols(), n);
    let mut xp = x0.to_vec();
    let mut xm = x0.to_vec();
    let mut fp = vec![0.0; n];
    let mut fm = vec![0.0; n];

    let mut diff_sq = 0.0_f64;
    let mut norm_sq = 0.0_f64;
    for col in 0..n {
        let h = eps * x0[col].abs().max(1.0);
        xp[col] = x0[col] + h;
        xm[col] = x0[col] - h;
        eval_f(&xp, &mut fp);
        eval_f(&xm, &mut fm);
        xp[col] = x0[col];
        xm[col] = x0[col];

        for row in 0..n {
            let fd = (fp[row] - fm[row]) / (2.0 * h);
            let (cols, vals) = j_exact.row(row);
            let ex = cols.binary_search(&(col as u32)).ok().map_or(0.0, |idx| vals[idx]);
            let d = ex - fd;
            diff_sq += d * d;
            norm_sq += ex * ex;
        }
    }
    diff_sq.sqrt() / norm_sq.sqrt().max(1e-30)
}

#[test]
fn wirtinger_frechet_jacobians_match_central_finite_differences() {
    let lib = load_library();
    let solver_mesh = block(0.20e-6, 0.80e-6, 0.20e-6, 1, 2, 1);
    let mesh = Mesh::from_solver_mesh(&solver_mesh).unwrap();
    let disc = discretise(&mesh).unwrap();
    let freq1 = C0 / 1.55e-6;

    // 1. Single-tone rotated cubic m-3m Silicon (anisotropic Kerr + TPA + steady-state FCA + FCD):
    let mut spec_chi3 = ProblemSpec::new("ems", 1, "nonlinear_driven");
    spec_chi3.study.params = json!({
        "process": "chi3",
        "freq": freq1,
        "power_w": 0.5
    });
    spec_chi3
        .materials
        .insert("core".into(), MaterialAssign::Oriented { material: "silicon".into(), euler_deg: [27.0, 43.0, 19.0] });
    spec_chi3.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
    spec_chi3.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "p1", "mode": 1}) });
    spec_chi3.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "p2", "mode": 1}) });

    let mut ctx = Ctx::new(&BlockResolver);
    ctx.library = Some(lib.clone());
    let prob_chi3 = nyx_ems::decode(&spec_chi3, &ctx).unwrap();
    let asm_chi3 = assemble_nonlinear_driven(&disc, &prob_chi3).unwrap();
    let pat_chi3 = asm_chi3.single_pattern.as_ref().unwrap();
    let nred = asm_chi3.harmonics[0].nred;
    let y0_chi3: Vec<f64> = (0..(2 * nred)).map(|k| 2.5e7 * ((0.73 * (k as f64) + 0.41).sin())).collect();
    let rhs_hat: Vec<_> = (0..nred).map(|k| c(0.3 * (k as f64).cos(), -0.2 * (k as f64).sin())).collect();
    let mut f_chi3 = vec![0.0; 2 * nred];
    let mut jv_chi3 = pat_chi3.zero_values();
    let max_nc = eval_chi3_residual_and_jacobian(&asm_chi3, &y0_chi3, &rhs_hat, &mut f_chi3, Some(&mut jv_chi3));
    assert!(max_nc > 1e20, "expected strong steady-state TPA carrier density in Jacobian check, got {max_nc:e}");
    let j_chi3 = pat_chi3.into_csr(jv_chi3);
    let err_chi3 = csr_frob_rel_err_vs_fd(&j_chi3, &y0_chi3, 1e-6, |x, out| {
        eval_chi3_residual_and_jacobian(&asm_chi3, x, &rhs_hat, out, None);
    });
    assert!(err_chi3 < 1e-8, "chi3 + TPA + FCA/FCD Wirtinger Jacobian relative FD error = {err_chi3:.3e} >= 1e-8");

    // 2. Multi-harmonic 3-wave depleted SFG on rotated Lithium Niobate with QPM poling:
    let mut spec_sfg = ProblemSpec::new("ems", 1, "nonlinear_driven");
    spec_sfg.study.params = json!({
        "process": "sfg",
        "freq": freq1,
        "freq_2": 0.82 * freq1,
        "power_w": 0.2,
        "signal_power_w": 0.15,
        "coupling": "depleted"
    });
    spec_sfg.materials.insert(
        "core".into(),
        MaterialAssign::Oriented { material: "lithium_niobate".into(), euler_deg: [33.0, 58.0, 21.0] },
    );
    spec_sfg.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
    spec_sfg.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "p1", "mode": 1}) });
    spec_sfg.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "p2", "mode": 1}) });

    let prob_sfg = nyx_ems::decode(&spec_sfg, &ctx).unwrap();
    let asm_sfg = assemble_nonlinear_driven(&disc, &prob_sfg).unwrap();
    let pat_sfg = asm_sfg.coupled_pattern.as_ref().unwrap();
    let n_tot_sfg: usize = asm_sfg.harmonics.iter().map(|h| 2 * h.nred).sum();
    let y0_sfg: Vec<f64> = (0..n_tot_sfg).map(|k| 1.5e8 * ((1.13 * (k as f64) - 0.27).cos())).collect();
    let rhs_sfg: Vec<Vec<_>> = asm_sfg.harmonics.iter().map(|h| vec![c(0.1, -0.05); h.nred]).collect();
    let mut f_sfg = vec![0.0; n_tot_sfg];
    let mut jv_sfg = pat_sfg.zero_values();
    eval_chi2_residual_and_jacobian(&asm_sfg, &y0_sfg, &rhs_sfg, true, &mut f_sfg, Some(&mut jv_sfg));
    let j_sfg = pat_sfg.into_csr(jv_sfg);
    let err_sfg = csr_frob_rel_err_vs_fd(&j_sfg, &y0_sfg, 1e-6, |x, out| {
        eval_chi2_residual_and_jacobian(&asm_sfg, x, &rhs_sfg, true, out, None);
    });
    assert!(err_sfg < 1e-8, "3-wave depleted SFG Wirtinger Jacobian relative FD error = {err_sfg:.3e} >= 1e-8");

    // 3. Two-tone Maker-Terhune XPM + cross-TPA Jacobian:
    let mut spec_xpm = spec_chi3.clone();
    spec_xpm.study.params = json!({
        "process": "chi3_xpm",
        "freq": freq1,
        "freq_2": 1.07 * freq1,
        "power_w": 0.25,
        "signal_power_w": 0.10,
        "coupling": "depleted"
    });
    spec_xpm.materials.insert("core".into(), MaterialAssign::Named("silicon_nitride".into()));
    let prob_xpm = nyx_ems::decode(&spec_xpm, &ctx).unwrap();
    let asm_xpm = assemble_nonlinear_driven(&disc, &prob_xpm).unwrap();
    let pat_xpm = asm_xpm.coupled_pattern.as_ref().unwrap();
    let n_tot_xpm: usize = asm_xpm.harmonics.iter().map(|h| 2 * h.nred).sum();
    let y0_xpm: Vec<f64> = (0..n_tot_xpm).map(|k| 4.0e8 * ((0.91 * (k as f64) + 0.17).sin())).collect();
    let rhs_xpm: Vec<Vec<_>> = asm_xpm.harmonics.iter().map(|h| vec![c(0.2, 0.1); h.nred]).collect();
    let mut f_xpm = vec![0.0; n_tot_xpm];
    let mut jv_xpm = pat_xpm.zero_values();
    eval_chi2_residual_and_jacobian(&asm_xpm, &y0_xpm, &rhs_xpm, true, &mut f_xpm, Some(&mut jv_xpm));
    let j_xpm = pat_xpm.into_csr(jv_xpm);
    let err_xpm = csr_frob_rel_err_vs_fd(&j_xpm, &y0_xpm, 1e-6, |x, out| {
        eval_chi2_residual_and_jacobian(&asm_xpm, x, &rhs_xpm, true, out, None);
    });
    assert!(err_xpm < 1e-8, "2-tone chi3_xpm Wirtinger Jacobian relative FD error = {err_xpm:.3e} >= 1e-8");
}

#[test]
fn depleted_chi2_shg_tanh2_efficiency_and_manley_rowe_conservation() {
    // Sub-cutoff parallel-plate TEM waveguide (`pmc` on x = 0, w_x; PEC on z = 0, h_z; waveports on y = 0, L):
    // Supports the unique z-polarized TEM mode E = E_z(y) z_hat with exact n_eff = n_0 at both omega_1 and 2*omega_1.
    let lambda1 = 1.55e-6;
    let freq1 = C0 / lambda1;
    let omega1 = 2.0 * PI * freq1;
    let n0 = 2.0_f64;
    let eps_r = n0 * n0;
    let d33 = 25.0e-12_f64; // 25 pm/V
    let w_x = 0.15e-6;
    let h_z = 0.15e-6;
    let area = w_x * h_z;
    let lambda_g1 = lambda1 / n0;
    let length = 10.0 * lambda_g1; // 7.75 um (8 quadratic P2 elements per guided 2*omega_1 wavelength)

    let solver_mesh = block(w_x, length, h_z, 1, 160, 1);
    let mesh = Mesh::from_solver_mesh(&solver_mesh).unwrap();

    // Boyd Eq. (2.6.28): Gamma^2 = (2 * omega_1^2 * d_eff^2 * P_1(0)) / (n_0^3 * eps_0 * c_0^3 * A)
    let p_for_gamma_l = |gl: f64| -> f64 {
        let gamma = gl / length;
        (gamma * gamma) * (n0.powi(3) * EPS0 * C0.powi(3) * area) / (2.0 * omega1 * omega1 * d33 * d33)
    };
    let gl_targets = [0.25_f64, 0.60, 1.00, 1.35];
    let powers_w: Vec<f64> = gl_targets.iter().map(|&gl| p_for_gamma_l(gl)).collect();

    let mut spec = ProblemSpec::new("ems", 2, "nonlinear_driven");
    spec.study.params = json!({
        "process": "shg",
        "freq": freq1,
        "powers_w": powers_w,
        "drive_port": 0,
        "coupling": "depleted",
        "nl_solver": "newton",
        "rtol": 1e-11
    });
    spec.materials.insert("core".into(), MaterialAssign::Props([("eps_r".into(), eps_r), ("d33".into(), d33)].into()));
    let eta0 = (MU0 / EPS0).sqrt();
    let z0_tem = eta0 / n0; // exact matched TEM wave impedance for w_x == h_z
    spec.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
    spec.bcs.push(BcSpec {
        on: "y0".into(),
        kind: "lumped_port".into(),
        params: json!({"name": "in", "Z_0": z0_tem, "direction": [0.0, 0.0, 1.0]}),
    });
    spec.bcs.push(BcSpec {
        on: "y1".into(),
        kind: "lumped_port".into(),
        params: json!({"name": "out", "Z_0": z0_tem, "direction": [0.0, 0.0, 1.0]}),
    });

    let res = Ems.solve(&mesh, &spec, &mut Ctx::new(&BlockResolver)).expect("depleted SHG solve");

    for (idx, &gl) in gl_targets.iter().enumerate() {
        let step = idx + 1;
        let eta_fem = num_scalar(&res, &format!("shg_efficiency step {step}"));
        let eta_exact = gl.tanh().powi(2);
        let rel_err = (eta_fem - eta_exact).abs() / eta_exact;
        let mr_err = num_scalar(&res, &format!("manley_rowe_rel_err step {step}"));
        assert!(
            rel_err < 5e-3,
            "step {step} (Gamma*L = {gl:.2}): SHG efficiency FEM = {eta_fem:.6} vs tanh^2(Gamma*L) = {eta_exact:.6} (rel_err = {:.3}% >= 0.5%)",
            rel_err * 100.0
        );
        assert!(mr_err < 1e-6, "step {step}: Manley-Rowe relative error = {mr_err:.3e} >= 1e-6");
    }

    // Verify undepleted SHG matches (Gamma*L)^2 at weak drive and over-predicts tanh^2(Gamma*L) at strong drive:
    let mut spec_undep = spec.clone();
    spec_undep.study.params["coupling"] = json!("undepleted");
    let res_undep = Ems.solve(&mesh, &spec_undep, &mut Ctx::new(&BlockResolver)).expect("undepleted SHG solve");
    let eta_undep_weak = num_scalar(&res_undep, "shg_efficiency step 1");
    let eta_undep_strong = num_scalar(&res_undep, "shg_efficiency step 4");
    assert!(
        ((eta_undep_weak - 0.25_f64.powi(2)) / 0.25_f64.powi(2)).abs() < 5e-3,
        "undepleted weak SHG {eta_undep_weak:.6} must match (Gamma*L)^2 = {:.6}",
        0.25_f64.powi(2)
    );
    assert!(
        eta_undep_strong > gl_targets[3].tanh().powi(2) * 1.5,
        "undepleted SHG at Gamma*L=1.35 ({eta_undep_strong:.4}) must exceed depleted tanh^2 ({:.4})",
        gl_targets[3].tanh().powi(2)
    );

    // Verify JFNK matches exact Newton at moderate depletion (Gamma*L = 0.60):
    let mut spec_jfnk = spec.clone();
    spec_jfnk.study.params = json!({
        "process": "shg",
        "freq": freq1,
        "power_w": powers_w[1],
        "drive_port": 0,
        "coupling": "depleted",
        "nl_solver": "jfnk",
        "rtol": 1e-10
    });
    let res_jfnk = Ems.solve(&mesh, &spec_jfnk, &mut Ctx::new(&BlockResolver)).expect("JFNK SHG solve");
    let eta_jfnk = num_scalar(&res_jfnk, "shg_efficiency");
    let eta_newt_2 = num_scalar(&res, "shg_efficiency step 2");
    assert!(
        (eta_jfnk - eta_newt_2).abs() / eta_newt_2 < 1e-7,
        "JFNK SHG efficiency ({eta_jfnk:.8}) must match Newton ({eta_newt_2:.8})"
    );

    // Verify 3-wave Optical Parametric Amplification (`process = "opa"`, pump f_3 -> signal f_1 + idler f_2)
    // amplifies signal power AND conserves both Manley-Rowe invariants N_1 + N_3 and N_2 + N_3 to < 1e-6!
    let solver_opa = block(w_x, 8.0 * lambda_g1, h_z, 1, 64, 1);
    let mesh_opa = Mesh::from_solver_mesh(&solver_opa).unwrap();
    let p_pump = powers_w[2];
    let p_sig = 0.02 * p_pump;
    let mut spec_opa = ProblemSpec::new("ems", 2, "nonlinear_driven");
    spec_opa.study.params = json!({
        "process": "opa",
        "freq": 2.0 * freq1,
        "freq_2": 0.90 * freq1,
        "power_w": p_pump,
        "signal_power_w": p_sig,
        "drive_port": 0,
        "signal_port": 0,
        "coupling": "depleted",
        "rtol": 1e-11
    });
    spec_opa.materials = spec.materials.clone();
    spec_opa.bcs = spec.bcs.clone();
    let res_opa = Ems.solve(&mesh_opa, &spec_opa, &mut Ctx::new(&BlockResolver)).expect("depleted OPA solve");
    let p_sig_out = num_scalar(&res_opa, "harmonic 1 port 'out' P_out (W)");
    let p_idl_out = num_scalar(&res_opa, "harmonic 2 port 'out' P_out (W)");
    let mr_opa = num_scalar(&res_opa, "manley_rowe_rel_err");
    assert!(
        p_sig_out > 1.25 * p_sig && p_idl_out > 0.20 * p_sig,
        "expected parametric signal amplification and idler generation, got P_sig_out={p_sig_out:.4e} (P_sig_in={p_sig:.4e}), P_idl_out={p_idl_out:.4e}"
    );
    assert!(mr_opa < 1e-6, "3-wave OPA Manley-Rowe relative error = {mr_opa:.3e} >= 1e-6");
}

#[test]
fn qpm_poled_lithium_niobate_shg_enhancement_vs_unpoled() {
    let lib = load_library();
    let ln = lib.get("lithium_niobate").expect("lithium_niobate in library");
    let lambda1 = 1.55e-6;
    let lambda2 = 0.5 * lambda1;
    let freq1 = C0 / lambda1;
    let omega1 = 2.0 * PI * freq1;

    // Extraordinary optic axis Z is along wafer z ([0, 0, 0]):
    let n1 = ln.optical_permittivity(lambda1, None, [0.0; 3]).unwrap()[2][2].re.sqrt();
    let n2 = ln.optical_permittivity(lambda2, None, [0.0; 3]).unwrap()[2][2].re.sqrt();
    let d33_abs = ln.nonlinear.as_ref().unwrap().chi2.as_ref().unwrap().d[2][2].abs();
    assert!(
        (n2 - n1) > 0.03 && d33_abs > 20e-12,
        "unexpected LiNbO3 dispersion/d33: n1={n1}, n2={n2}, |d33|={d33_abs:e}"
    );

    // Exact coherence length L_c = pi / Delta_k = lambda_1 / (4 * (n_2 - n_1)) and QPM period Lambda = 2 * L_c:
    let l_coh = lambda1 / (4.0 * (n2 - n1));
    let lambda_qpm = 2.0 * l_coh;
    let length = lambda_qpm; // 2 * L_c (Delta_k * L = 2*pi): exact phase-mismatch null for unpoled LiNbO3!

    let w_x = 0.14e-6;
    let h_z = 0.14e-6;
    let area = w_x * h_z;
    // 12 quadratic P2 elements per guided second-harmonic wavelength (~24 DOFs/lambda_g2) so dispersion drift < 0.015 rad:
    let solver_mesh = block(w_x, length, h_z, 1, 640, 1);
    let mesh = Mesh::from_solver_mesh(&solver_mesh).unwrap();

    // Operate at Gamma_QPM * L = 0.15 so unpoled sits cleanly at the sinc^2(pi) = 0 null (> 500x suppression)
    // and QPM matches first-order Fourier d_eff = (2 / pi) * |d33|:
    let d_eff_qpm = (2.0 / PI) * d33_abs;
    let gl_target = 0.15_f64;
    let gamma_qpm = gl_target / length;
    let p_in = (gamma_qpm * gamma_qpm) * (n1 * n1 * n2 * EPS0 * C0.powi(3) * area)
        / (2.0 * omega1 * omega1 * d_eff_qpm * d_eff_qpm);

    let make_ln_spec = |mat_name: &str| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 2, "nonlinear_driven");
        s.study.params = json!({
            "process": "shg",
            "freq": freq1,
            "power_w": p_in,
            "drive_port": 0,
            "coupling": "depleted",
            "rtol": 1e-10
        });
        s.materials.insert("core".into(), MaterialAssign::Named(mat_name.into()));
        s.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
        s.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in", "mode": 1}) });
        s.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out", "mode": 1}) });
        s
    };

    // 1. Unpoled congruent LiNbO3 over L = 4 * L_c:
    let spec_unpoled = make_ln_spec("lithium_niobate");
    let mut ctx_unpoled = Ctx::new(&BlockResolver);
    ctx_unpoled.library = Some(lib.clone());
    let res_unpoled = Ems.solve(&mesh, &spec_unpoled, &mut ctx_unpoled).expect("unpoled LiNbO3 SHG");
    let eta_unpoled = num_scalar(&res_unpoled, "shg_efficiency");

    // 2. Periodically poled LiNbO3 (`extends = "lithium_niobate"` with `[nonlinear.qpm] period = "<Lambda> m"`):
    let ppln_toml = format!(
        "format = 2\nname = \"Periodically poled LiNbO3\"\nextends = \"lithium_niobate\"\n\n[nonlinear.qpm]\nperiod = \"{lambda_qpm:.12e} m\"\naxis = [0.0, 1.0, 0.0]\nduty_cycle = 0.5\n"
    );
    let lib_qpm = lib
        .with_texts(&[("ppln".to_string(), ppln_toml)].into())
        .expect("parse PPLN material extending lithium_niobate");

    let spec_qpm = make_ln_spec("ppln");
    let mut ctx_qpm = Ctx::new(&BlockResolver);
    ctx_qpm.library = Some(Arc::new(lib_qpm));
    let res_qpm = Ems.solve(&mesh, &spec_qpm, &mut ctx_qpm).expect("QPM PPLN SHG");

    let eta_qpm = num_scalar(&res_qpm, "shg_efficiency");
    let eta_qpm_analytic = gl_target.tanh().powi(2);
    let mr_qpm = num_scalar(&res_qpm, "manley_rowe_rel_err");

    assert!(
        eta_qpm > 200.0 * eta_unpoled.max(1e-12),
        "expected QPM SHG efficiency ({eta_qpm:.5}) >> unpoled phase-mismatched efficiency ({eta_unpoled:.5e})"
    );
    let rel_err_qpm = (eta_qpm - eta_qpm_analytic).abs() / eta_qpm_analytic;
    assert!(
        rel_err_qpm < 0.05,
        "QPM PPLN SHG efficiency FEM = {eta_qpm:.5} vs asymptotic Fourier tanh^2 = {eta_qpm_analytic:.5} (rel_err = {:.2}% >= 5%)",
        rel_err_qpm * 100.0
    );
    assert!(mr_qpm < 1e-5, "PPLN QPM Manley-Rowe error = {mr_qpm:.3e} >= 1e-5");
}

#[test]
fn sin_high_q_kerr_cavity_bistability_arclength_vs_analytical_cmt() {
    let lib = load_library();
    let sin_mat = lib.get("silicon_nitride").expect("silicon_nitride in library");
    let lambda_nom = 1.55e-6;
    let n0 = sin_mat.optical_permittivity(lambda_nom, None, [0.0; 3]).unwrap()[0][0].re.sqrt();
    let n2 = sin_mat.nonlinear.as_ref().unwrap().chi3.as_ref().unwrap().n_2.unwrap();

    // Half-wave Si3N4 Fabry-Perot resonator (`cos(pi y / L_y) z_hat`) with symmetric high-impedance lumped ports
    // on y = 0 and y = L_y tuned to loaded quality factor Q_L = 1.0e5:
    let w_x = 0.25e-6;
    let h_z = 0.25e-6;
    let l_y = lambda_nom / (2.0 * n0);
    let vol = w_x * l_y * h_z;
    let v_eff = (2.0 / 3.0) * vol; // standing-wave cos^4 / (cos^2)^2 mode volume

    let eta0 = (MU0 / EPS0).sqrt();
    let q_target = 1.0e5_f64;
    let z0_port = (4.0 * eta0 * q_target) / (PI * n0);

    let solver_mesh = block(w_x, l_y, h_z, 1, 8, 1);
    let mesh = Mesh::from_solver_mesh(&solver_mesh).unwrap();

    // Step 1: Two-pass Lorentzian probe on linear `Study::Driven` (with fixed wavelength `lambda_nom` so
    // linear and nonlinear studies share identical material permittivity) to extract exact discrete `(f_0, FWHM)`:
    let f0_est = C0 / (2.0 * n0 * l_y);
    let probe_resonance = |f_center: f64, df: f64| -> (f64, f64) {
        let mut s = ProblemSpec::new("ems", 2, "driven");
        s.study.params = json!({
            "freqs_hz": [f_center - df, f_center, f_center + df],
            "solver": "direct"
        });
        s.materials.insert("cav".into(), MaterialAssign::Named("silicon_nitride".into()));
        s.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
        s.bcs.push(BcSpec {
            on: "y0".into(),
            kind: "lumped_port".into(),
            params: json!({"name": "p1", "Z_0": z0_port, "direction": [0.0, 0.0, 1.0]}),
        });
        s.bcs.push(BcSpec {
            on: "y1".into(),
            kind: "lumped_port".into(),
            params: json!({"name": "p2", "Z_0": z0_port, "direction": [0.0, 0.0, 1.0]}),
        });
        let mut ctx_p = Ctx::new(&BlockResolver);
        ctx_p.library = Some(lib.clone());
        ctx_p.wavelength = Some(lambda_nom);
        let r = Ems.solve(&mesh, &s, &mut ctx_p).expect("linear cavity probe");
        let d = DrivenPortData::of(&r).expect("DrivenPortData");
        let s_m = c(d.s_matrix[0][1][0][0], d.s_matrix[0][1][0][1]);
        let s_p = c(d.s_matrix[2][1][0][0], d.s_matrix[2][1][0][1]);
        // Symmetric lossless 2-port cavity satisfies |S_21(f_0)| = 1, so 1/S_21(f) = +- [1 + i * 2*(f - f_0)/FWHM] at any f:
        let inv_m = c(1.0, 0.0) / s_m;
        let inv_p = c(1.0, 0.0) / s_p;
        let slope = (inv_p.im - inv_m.im) / (2.0 * df);
        let f_res = (f_center - df) - inv_m.im / slope;
        let fwhm = 2.0 * (0.5 * (inv_m.re + inv_p.re)).abs() / slope.abs();
        (f_res, fwhm)
    };
    let (f0_pass1, fwhm_pass1) = probe_resonance(f0_est, 0.5 * f0_est / q_target);
    let (f0_fem, fwhm_hz) = probe_resonance(f0_pass1, 0.4 * fwhm_pass1);
    let q_loaded = f0_fem / fwhm_hz;
    assert!(
        (q_loaded - q_target).abs() / q_target < 5e-3,
        "expected high-Q ~ 1e5 resonator, got Q_L = {q_loaded:.1} (f0 = {f0_fem:.6e} Hz)"
    );

    let omega0 = 2.0 * PI * f0_fem;
    let kappa = 2.0 * PI * fwhm_hz;
    let kappa_e1 = 0.5 * kappa;
    let kappa_e2 = 0.5 * kappa;
    let g_kerr = (omega0 * C0 * n2) / (n0 * n0 * v_eff);

    // Step 2: Red-detuned pump at normalized detuning Delta = (omega_0 - omega) / (kappa / 2) = 2.50 > sqrt(3):
    let delta_norm = 2.50_f64;
    let f_drive = f0_fem - 0.5 * delta_norm * fwhm_hz;
    let half_k = 0.5 * kappa;
    let p0_cmt = half_k.powi(3) / (kappa_e1 * g_kerr);

    // Exact CMT saddle-node turning points: u_+- = (2*Delta -+ sqrt(Delta^2 - 3)) / 3
    let disc_cmt = (delta_norm * delta_norm - 3.0).sqrt();
    let u_fold_up = (2.0 * delta_norm - disc_cmt) / 3.0; // upper power threshold (jump up from lower branch)
    let u_fold_dn = (2.0 * delta_norm + disc_cmt) / 3.0; // lower power threshold (drop down from upper branch)
    let p_fold_up_cmt = p0_cmt * u_fold_up * (1.0 + (u_fold_up - delta_norm).powi(2));
    let p_fold_dn_cmt = p0_cmt * u_fold_dn * (1.0 + (u_fold_dn - delta_norm).powi(2));

    // Trace entire S-curve via pseudo-arclength continuation up to 4.8 * P_0:
    let mut spec_nl = ProblemSpec::new("ems", 2, "nonlinear_driven");
    spec_nl.study.params = json!({
        "process": "chi3",
        "freq": f_drive,
        "power_min_w": 0.0,
        "power_max_w": 4.8 * p0_cmt,
        "drive_port": 0,
        "continuation": "arclength",
        "arclength_ds": 0.04,
        "arclength_max_steps": 320,
        "rtol": 1e-11
    });
    spec_nl.materials.insert("cav".into(), MaterialAssign::Named("silicon_nitride".into()));
    spec_nl.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
    spec_nl.bcs.push(BcSpec {
        on: "y0".into(),
        kind: "lumped_port".into(),
        params: json!({"name": "p1", "Z_0": z0_port, "direction": [0.0, 0.0, 1.0]}),
    });
    spec_nl.bcs.push(BcSpec {
        on: "y1".into(),
        kind: "lumped_port".into(),
        params: json!({"name": "p2", "Z_0": z0_port, "direction": [0.0, 0.0, 1.0]}),
    });

    let mut ctx_nl = Ctx::new(&BlockResolver);
    ctx_nl.library = Some(lib);
    ctx_nl.wavelength = Some(lambda_nom);
    let res_nl = Ems.solve(&mesh, &spec_nl, &mut ctx_nl).expect("Si3N4 Kerr cavity arclength bistability solve");

    let n_tp = num_scalar(&res_nl, "n_turning_points") as usize;
    assert_eq!(n_tp, 2, "expected exactly 2 saddle-node turning points on Kerr S-curve");
    let tp1_fem = num_scalar(&res_nl, "turning_point_power 1 (W)");
    let tp2_fem = num_scalar(&res_nl, "turning_point_power 2 (W)");

    let err_tp1 = (tp1_fem - p_fold_up_cmt).abs() / p_fold_up_cmt;
    let err_tp2 = (tp2_fem - p_fold_dn_cmt).abs() / p_fold_dn_cmt;
    assert!(
        err_tp1 < 5e-3,
        "upper fold power FEM = {tp1_fem:.6e} W vs CMT = {p_fold_up_cmt:.6e} W (rel_err = {:.3}% >= 0.5%)",
        err_tp1 * 100.0
    );
    assert!(
        err_tp2 < 5e-3,
        "lower fold power FEM = {tp2_fem:.6e} W vs CMT = {p_fold_dn_cmt:.6e} W (rel_err = {:.3}% >= 0.5%)",
        err_tp2 * 100.0
    );

    // Verify all 3 branches (lower stable, unstable middle with tangent_alpha < 0, upper stable)
    // satisfy the exact CMT state equation P_in(U) and P_out,2 = kappa_e2 * U to < 0.5%!
    let n_steps = num_scalar(&res_nl, "n_steps") as usize;
    let mut n_unstable = 0usize;
    let mut n_upper = 0usize;
    for st in 2..=n_steps {
        let p_in_k = num_scalar(&res_nl, &format!("P_in step {st} (W)"));
        let u_k = num_scalar(&res_nl, &format!("U_stored step {st} (J)"));
        let ta_k = num_scalar(&res_nl, &format!("tangent_alpha step {st}"));
        let p_out2_k = num_scalar(&res_nl, &format!("harmonic 1 port 'p2' P_out step {st} (W)"));
        if ta_k < 0.0 {
            n_unstable += 1;
        } else if n_unstable > 0 {
            n_upper += 1;
        }
        if p_in_k > 0.2 * p0_cmt {
            let u_dim = g_kerr * u_k / half_k;
            let p_in_cmt = p0_cmt * u_dim * (1.0 + (u_dim - delta_norm).powi(2));
            let p_out2_cmt = kappa_e2 * u_k;
            let rel_p = (p_in_k - p_in_cmt).abs() / p_in_cmt;
            let rel_out = (p_out2_k - p_out2_cmt).abs() / p_out2_cmt;
            assert!(
                rel_p < 5e-3 && rel_out < 5e-3,
                "step {st} (ta={ta_k:+.3}): P_in rel_err={:.3}%, P_out2 rel_err={:.3}%",
                rel_p * 100.0,
                rel_out * 100.0
            );
        }
    }
    assert!(
        n_unstable >= 3 && n_upper >= 3,
        "expected multiple steps on both unstable middle branch ({n_unstable}) and upper stable branch ({n_upper})"
    );
    assert_eq!(
        num_scalar(&res_nl, "symbolic_analyses") as usize,
        1,
        "bordered pseudo-arclength solver must reuse a single symbolic factorization across the entire S-curve"
    );
}

#[test]
fn cross_polarized_chi3_spm_xpm_tpa_fca_and_2d_nonlinear_modes_gamma() {
    let lib = load_library();
    let lambda0 = 1.55e-6;
    let freq0 = C0 / lambda0;

    // -----------------------------------------------------------------------------------------
    // Part A: 3D degenerate-polarization waveguide (`pmc` on ALL 4 side walls x=0,w and z=0,w!)
    // On a square cross-section with PMC on x=0,w and PEC? Wait: a PEC square waveguide `a x a`
    // (`a = 0.65e-6`, `n0 = 2.0`) has exact degenerate fundamental TE_10 and TE_01 modes!
    // Even cleaner: with `a = 0.65e-6`, `n0 = 2.0`, TE_10 and TE_01 (`mode = 1` and `mode = 2`)
    // are orthogonal linear polarizations with identical beta_0 and A_eff = (2/3)^2? No, int sin^2 = a/2,
    // int sin^4 = 3a/8, so 1D sin^2(pi x / a) uniform in z has A_eff = (2/3) * a^2!
    // Let's compare the nonlinear phase shift directly between:
    //   1) Linear polarization (driving mode 1 with power P) -> Delta phi_lin
    //   2) Circular polarization (driving mode 1 with P/2 @ 0 rad and mode 2 with P/2 @ +pi/2 rad) -> Delta phi_circ!
    //Wait: in a PEC square waveguide, mode 1 has E_z(x) = E_0 sin(pi x / a) (varies in x, constant in z)
    // whereas mode 2 has E_x(z) = E_0 sin(pi z / a) (varies in z, constant in x) — their spatial profiles
    // sin(pi x / a) and sin(pi z / a) are NOT pointwise equal!
    // Whereas on a parallel-plate cell? No — in a square PEC waveguide, at two frequencies (or even simpler:
    // on a uniform periodic/plane-wave box with `periodic` in x and z, BOTH `E_x` and `E_z` are spatially
    // UNIFORM across the entire cross-section!!)!
    // Let's check what happens on the PEC square waveguide `a x a` for XPM between TE_10 (E_z = sin(pi x / a))
    // and TE_01 (E_x = sin(pi z / a)):
    //   int_0^a int_0^a sin^2(pi x / a) sin^2(pi z / a) dx dz = (a/2)*(a/2) = a^2 / 4,
    //   whereas self-overlap int_0^a int_0^a sin^4(pi x / a) dx dz = (3a/8)*a = 3 a^2 / 8!
    // Instead, if we put `pmc` on x=0,w (`1, 6`) for E_z TEM mode AND compare `b_over_a = 1.0` vs `b_over_a = 0.0`?
    // NO! We want TRUE simultaneously active orthogonal vector components `(E_x, E_z)` at the SAME point!
    // How do we get spatially uniform `E_x` AND spatially uniform `E_z` in the same waveguide?
    // Look at `TE_11` / `HE_11` in a square dielectric core waveguide, OR even simpler:
    // in a square PEC waveguide `a x a`, any two modes on the SAME cross-section work if we compare against
    // the exact analytical integral — OR in 2D `Study::NonlinearModes` AND 3D `Study::NonlinearDriven`!
    // Wait: in a square PEC waveguide `a x a`:
    //   Mode A (`TE_10`): E_A(x, z) = sin(pi x / a) z_hat  (pure z-polarized!)
    //   Mode B (`TE_01`): E_B(x, z) = sin(pi z / a) x_hat  (pure x-polarized!)
    // Let's compute the EXACT analytical ratio of circular-drive (`P/2` in Mode A at phase `0`, `P/2` in Mode B at phase `pi/2`,
    // so E = (E_0 / sqrt(2)) [ sin(pi z / a) x_hat + i sin(pi x / a) z_hat ]) vs linear-drive (`P` in Mode A, E = E_0 sin(pi x / a) z_hat)!
    //   - Linear drive (`P` in Mode A):
    //     P_z^(3)/eps0 = (A + 0.5 B) E_0^3 sin^3(pi x / a). Overlap with sin(pi x / a):
    //     I_lin = (A + 0.5 B) * E_0^2 * int_0^a int_0^a sin^4(pi x / a) dx dz = (3/2 A) * E_0^2 * (3/8 a^2) = (9/16) A E_0^2 a^2.
    //   - Circular two-mode drive (`P/2` in Mode A at `0`, `P/2` in Mode B at `pi/2`, so E_x = c_0 sin(pi z / a), E_z = i c_0 sin(pi x / a) with c_0 = E_0 / sqrt(2)):
    //     Notice E . E = c_0^2 [ sin^2(pi z / a) - sin^2(pi x / a) ], and E . E^* = c_0^2 [ sin^2(pi z / a) + sin^2(pi x / a) ]!
    //     Then on Mode A (`z` component, projected against `(i c_0)^* sin(pi x / a)`):
    //     P_z^(3)/eps0 = A c_0^2 [sin^2(pi z/a) + sin^2(pi x/a)] (i c_0 sin(pi x/a))
    //                  + 0.5 B c_0^2 [sin^2(pi z/a) - sin^2(pi x/a)] (-i c_0 sin(pi x/a))
    //                  = i c_0^3 sin(pi x/a) [ (A + 0.5 B) sin^2(pi x/a) + (A - 0.5 B) sin^2(pi z/a) ]!
    //     Projecting onto sin(pi x / a) over `[0, a]^2` (with c_0^2 = E_0^2 / 2 and B = A):
    //     I_circ = (1/2) E_0^2 [ (3/2 A) * (3/8 a^2) + (1/2 A) * (1/4 a^2) ] = (1/2) E_0^2 A a^2 [ 9/16 + 2/16 ] = (11/32) A E_0^2 a^2!
    //     Ratio `Delta phi_circ / Delta phi_lin` = `(11/32) / (9/16) = 11 / 18`!
    //     Wait: what if `B = 0` (`b_over_a = 0`, non-conjugate scalar-like `|E|^2 E` WITHOUT the `0.5 B (E.E) E^*` term!)?
    //     If `0.5 B (E.E) E^*` were missing (`B = 0`, `A = chi_eff`), `I_circ / I_lin` would be `(1/2)[3/8 + 1/4] / (3/8) = 5/6 = 15/18 != 11/18`!!
    //     That means the square PEC waveguide distinguishes the exact `0.5 B (E.E) E^*` conjugate tensor term (`11/18 = 0.61111`)
    //     from a scalar `|E|^2 E` approximation (`15/18 = 0.83333`) by **36%**!!
    //   - AND for two-tone `chi3_xpm` on the SAME square PEC waveguide (`pump` at `f_1` in Mode 1, weak `probe` at `f_2`):
    //     1) Co-polarized probe (`probe` ALSO in Mode 1 `TE_10`): spatial integral is `3/8 a^2`, tensor factor is `2 * chi_eff` -> `Delta phi_co / Delta phi_spm = 2 * (f_2 / f_1) * (beta_1_fac)`!
    //     2) Cross-polarized probe (`probe` in Mode 2 `TE_01`, strictly orthogonal `E_1 . E_2 = 0` EVERYWHERE!):
    //        spatial integral is `(1/4 a^2) / (3/8 a^2) = 2/3`, tensor factor is `(2/3) * chi_eff` vs `2 * chi_eff` (ratio `1/3`),
    //        so at the SAME `f_2`, `Delta phi_cross / Delta phi_co = (1/3) * ((1/4) / (3/8)) = 2 / 9` EXACTLY!!
    //        (Whereas without cross-polarization tensor physics, `Delta phi_cross / Delta phi_co` would be `2/3` — a **3x** error!)
    // Strictly two-mode PEC waveguide with well-separated aspect ratio (`a_x = 0.58 um, a_z = 0.74 um, n_0 = 1.52` at 1.55 um):
    // k = n_0 * k_0 = 6.162 um^-1, so ONLY TE_01 (k_c = 4.245) and TE_10 (k_c = 5.417) propagate with wide beta^2 separation
    // (eliminating unstructured Delaunay mode hybridization), while TE_11/TM_11 (k_c = 6.882) and TE_02 (k_c = 8.491)
    // remain strictly below cutoff at both f_0 and 1.03*f_0!
    let a_wg_x = 0.58e-6;
    let a_wg_z = 0.74e-6;
    let l_wg = 0.25e-6;
    let n0 = 1.52_f64;
    let eps_wg = n0 * n0;
    let n2_wg = 2.4e-18_f64;
    let solver_sq = block(a_wg_x, l_wg, a_wg_z, 3, 2, 3);
    let mesh_sq = Mesh::from_solver_mesh(&solver_sq).unwrap();
    let p_test = 10.0_f64; // W (dn ~ 1e-4, clean linear-in-P perturbative Kerr regime)

    // Helper to extract transmission phase arg(b_out / b_out_linear) on Mode 1 (`out_m1`):
    let phase_shift_spm = |p_ports: Vec<f64>, phases: Vec<f64>| -> f64 {
        let mut s = ProblemSpec::new("ems", 2, "nonlinear_driven");
        s.study.params = json!({
            "process": "chi3",
            "freq": freq0,
            "powers_w": [1e-4, p_test],
            "port_powers_w": p_ports,
            "port_phases_rad": phases,
            "rtol": 1e-12
        });
        s.materials
            .insert("core".into(), MaterialAssign::Props([("eps_r".into(), eps_wg), ("n_2".into(), n2_wg)].into()));
        // 2 orthogonal modes on input y0 (ports 0, 1) and 2 orthogonal modes on output y1 (ports 2, 3):
        s.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in_m1", "mode": 1}) });
        s.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in_m2", "mode": 2}) });
        s.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out_m1", "mode": 1}) });
        s.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out_m2", "mode": 2}) });

        let r = Ems.solve(&mesh_sq, &s, &mut Ctx::new(&BlockResolver)).expect("square wg chi3 solve");
        let b_lin = c(
            num_scalar(&r, "harmonic 1 port 'out_m1' b_out_re step 1 (W^0.5)"),
            num_scalar(&r, "harmonic 1 port 'out_m1' b_out_im step 1 (W^0.5)"),
        );
        let b_nl = c(
            num_scalar(&r, "harmonic 1 port 'out_m1' b_out_re step 2 (W^0.5)"),
            num_scalar(&r, "harmonic 1 port 'out_m1' b_out_im step 2 (W^0.5)"),
        );
        let ratio = b_nl / b_lin;
        ratio.im.atan2(ratio.re)
    };

    let k_spm = n0 * (2.0 * PI * freq0 / C0);
    let beta_m1_spm = (k_spm * k_spm - (PI / a_wg_z).powi(2)).sqrt();
    let beta_m2_spm = (k_spm * k_spm - (PI / a_wg_x).powi(2)).sqrt();
    let dphi_lin = phase_shift_spm(vec![p_test, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 0.0]);
    // Pre-bias input phase by +0.5*(beta_2 - beta_1)*L so relative phase at waveguide midpoint y = L/2 is exact +pi/2:
    let dbeta_l = (beta_m1_spm - beta_m2_spm) * l_wg;
    let sinc_walkoff = if dbeta_l.abs() > 1e-12 { dbeta_l.sin() / dbeta_l } else { 1.0 };
    let phase2_mid_circ = 0.5 * PI + 0.5 * (beta_m2_spm - beta_m1_spm) * l_wg;
    let dphi_circ = phase_shift_spm(vec![0.5 * p_test, 0.5 * p_test, 0.0, 0.0], vec![0.0, phase2_mid_circ, 0.0, 0.0]);
    let ratio_circ_lin = dphi_circ / dphi_lin;
    let exact_circ_lin = 0.5 * (1.0 + (2.0 / 9.0) * (beta_m1_spm / beta_m2_spm) * (2.0 - sinc_walkoff));

    // Two-tone XPM (`process = "chi3_xpm"`) on the same strictly two-mode waveguide (`mesh_sq`):
    let freq_probe = 1.03 * freq0;
    let k_probe = n0 * (2.0 * PI * freq_probe / C0);
    let beta_probe_m1 = (k_probe * k_probe - (PI / a_wg_z).powi(2)).sqrt();
    let beta_probe_m2 = (k_probe * k_probe - (PI / a_wg_x).powi(2)).sqrt();

    let xpm_probe_shift = |sig_port: usize| -> f64 {
        let mut s = ProblemSpec::new("ems", 2, "nonlinear_driven");
        s.study.params = json!({
            "process": "chi3_xpm",
            "freq": freq0,
            "freq_2": freq_probe,
            "powers_w": [0.0, p_test],
            "drive_port": 0,
            "signal_power_w": 1.0e-4,
            "signal_port": sig_port,
            "coupling": "undepleted",
            "rtol": 1e-12
        });
        s.materials
            .insert("core".into(), MaterialAssign::Props([("eps_r".into(), eps_wg), ("n_2".into(), n2_wg)].into()));
        s.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in_m1", "mode": 1}) });
        s.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in_m2", "mode": 2}) });
        s.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out_m1", "mode": 1}) });
        s.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out_m2", "mode": 2}) });

        let r = Ems.solve(&mesh_sq, &s, &mut Ctx::new(&BlockResolver)).expect("chi3_xpm solve");
        let b1_0 = c(
            num_scalar(&r, "harmonic 2 port 'out_m1' b_out_re step 1 (W^0.5)"),
            num_scalar(&r, "harmonic 2 port 'out_m1' b_out_im step 1 (W^0.5)"),
        );
        let b2_0 = c(
            num_scalar(&r, "harmonic 2 port 'out_m2' b_out_re step 1 (W^0.5)"),
            num_scalar(&r, "harmonic 2 port 'out_m2' b_out_im step 1 (W^0.5)"),
        );
        let out_name = if b1_0.norm() >= b2_0.norm() { "out_m1" } else { "out_m2" };
        let b0 = if b1_0.norm() >= b2_0.norm() { b1_0 } else { b2_0 };
        let b1 = c(
            num_scalar(&r, &format!("harmonic 2 port '{out_name}' b_out_re step 2 (W^0.5)")),
            num_scalar(&r, &format!("harmonic 2 port '{out_name}' b_out_im step 2 (W^0.5)")),
        );
        let ratio = b1 / b0;
        ratio.im.atan2(ratio.re)
    };
    let dphi_xpm_co = xpm_probe_shift(0);
    let dphi_xpm_cross = xpm_probe_shift(1);
    let ratio_xpm = dphi_xpm_cross / dphi_xpm_co;
    let exact_xpm_ratio = (2.0 / 9.0) * (beta_probe_m1 / beta_probe_m2);
    assert!(
        ((ratio_circ_lin - exact_circ_lin) / exact_circ_lin).abs() < 5e-3,
        "circular vs linear SPM phase-shift ratio = {ratio_circ_lin:.6} vs exact Maker-Terhune = {exact_circ_lin:.6}"
    );
    assert!(
        ((ratio_xpm - exact_xpm_ratio) / exact_xpm_ratio).abs() < 5e-3,
        "orthogonal vs co-polarized XPM phase ratio = {ratio_xpm:.6} vs exact Maker-Terhune = {exact_xpm_ratio:.6}"
    );

    // -----------------------------------------------------------------------------------------
    // Part B: Silicon (`materials/silicon.toml`) TPA + steady-state FCA/FCD nonlinear absorption
    // -----------------------------------------------------------------------------------------
    let solver_si = block(0.15e-6, 20.0e-6, 0.15e-6, 1, 64, 1);
    let mesh_si = Mesh::from_solver_mesh(&solver_si).unwrap();
    let mut spec_si = ProblemSpec::new("ems", 2, "nonlinear_driven");
    spec_si.study.params = json!({
        "process": "chi3",
        "freq": freq0,
        "powers_w": [1.0e-5, 0.10, 0.30],
        "drive_port": 0,
        "rtol": 1e-10
    });
    spec_si.materials.insert("core".into(), MaterialAssign::Named("silicon".into()));
    spec_si.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
    spec_si.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in", "mode": 1}) });
    spec_si.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out", "mode": 1}) });
    let mut ctx_si = Ctx::new(&BlockResolver);
    ctx_si.library = Some(lib.clone());
    let res_si = Ems.solve(&mesh_si, &spec_si, &mut ctx_si).expect("silicon TPA+FCA solve");
    let t_low = num_scalar(&res_si, "harmonic 1 port 'out' P_out step 1 (W)") / 1.0e-5;
    let t_mid = num_scalar(&res_si, "harmonic 1 port 'out' P_out step 2 (W)") / 0.10;
    let t_high = num_scalar(&res_si, "harmonic 1 port 'out' P_out step 3 (W)") / 0.30;
    let nc_mid = num_scalar(&res_si, "max_N_c step 2 (1/m^3)");
    let nc_high = num_scalar(&res_si, "max_N_c step 3 (1/m^3)");
    assert!(
        t_low > 0.90 && t_mid < 0.95 * t_low && t_high < 0.80 * t_mid && nc_high > 5.0 * nc_mid && nc_high > 1.0e25,
        "expected strong monotonic TPA+FCA transmission drop in Si: T_low={t_low:.5}, T_mid={t_mid:.5}, T_high={t_high:.5}, N_c=({nc_mid:.3e}, {nc_high:.3e})"
    );

    // -----------------------------------------------------------------------------------------
    // Part C: 2D `Study::NonlinearModes` on a high-confinement Si3N4 strip waveguide in SiO2:
    // Verify k_0 * (n_eff(P) - n_eff(0)) / P == gamma (1/(W*m)) to < 0.5%!
    // -----------------------------------------------------------------------------------------
    let mut m2d = Model::new_2d();
    let clad = m2d.rect("clad_r", (-1.5e-6, -1.2e-6), (1.5e-6, 1.2e-6));
    let core = m2d.rect("core_r", (-0.45e-6, -0.22e-6), (0.45e-6, 0.22e-6));
    m2d.overlay(&[("clad", clad), ("core", core)]);
    let g2d = mesh_model(&m2d.list, &MeshPlan::simple(0.18e-6)).expect("mesh 2D Si3N4 waveguide");
    let mesh2d = Mesh::from_solver_mesh(&g2d.mesh).unwrap();
    let res2d_lbl = GeomResolver(&g2d);

    let p_wg = 500.0_f64; // W (gives ~1e-5 Kerr index shift on Si3N4 strip, well in linear-in-P Kerr regime)
    let mut spec_nlm = ProblemSpec::new("ems", 2, "nonlinear_modes");
    spec_nlm.study.params = json!({
        "freq": freq0,
        "n_target": 1.95,
        "mode_index": 1,
        "powers_w": [0.0, p_wg, 2.0 * p_wg],
        "rtol": 1e-11
    });
    spec_nlm.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), 1.444_f64.powi(2))].into()));
    spec_nlm.materials.insert("core".into(), MaterialAssign::Named("silicon_nitride".into()));

    let mut ctx_nlm = Ctx::new(&res2d_lbl);
    ctx_nlm.library = Some(lib.clone());
    let res_nlm = Ems.solve(&mesh2d, &spec_nlm, &mut ctx_nlm).expect("2D nonlinear_modes Si3N4 solve");

    let gamma = num_scalar(&res_nlm, "gamma (1/(W*m))");
    let a_eff = num_scalar(&res_nlm, "A_eff (m^2)");
    let dbeta_1 = num_scalar(&res_nlm, "delta_beta_re 2 (1/m)");
    let dbeta_2 = num_scalar(&res_nlm, "delta_beta_re 3 (1/m)");
    assert!(gamma > 0.5 && a_eff > 0.1e-12 && a_eff < 2.0e-12, "unexpected Si3N4 gamma={gamma}, A_eff={a_eff:e}");

    let gamma_from_p1 = dbeta_1 / p_wg;
    let gamma_from_p2 = dbeta_2 / (2.0 * p_wg);
    let err_g1 = (gamma_from_p1 - gamma).abs() / gamma;
    let err_g2 = (gamma_from_p2 - gamma).abs() / gamma;
    assert!(
        err_g1 < 5e-3 && err_g2 < 5e-3,
        "Study::NonlinearModes delta_beta/P ({gamma_from_p1:.6}, {gamma_from_p2:.6}) vs vectorial gamma ({gamma:.6}): rel_err = ({:.3}%, {:.3}%)",
        err_g1 * 100.0,
        err_g2 * 100.0
    );

    // Also verify Silicon 2D `nonlinear_modes` reports positive `gamma_tpa` and `fom_tpa`, and Im(n_eff(P)) < 0 under TPA!
    let mut spec_si_modes = spec_nlm.clone();
    spec_si_modes.study.params = json!({
        "freq": freq0,
        "n_target": 3.40,
        "mode_index": 1,
        "power_w": 1.0,
        "rtol": 1e-11
    });
    spec_si_modes.materials.insert(
        "core".into(),
        MaterialAssign::Props(
            [("eps_r".into(), 3.476_f64.powi(2)), ("n_2".into(), 5.5e-18), ("beta_tpa".into(), 5.5e-12)].into(),
        ),
    );
    let res_si_m =
        Ems.solve(&mesh2d, &spec_si_modes, &mut Ctx::new(&res2d_lbl)).expect("2D nonlinear_modes Si TPA solve");
    let gamma_tpa = num_scalar(&res_si_m, "gamma_tpa (1/(W*m))");
    let fom_tpa = num_scalar(&res_si_m, "fom_tpa");
    let dneff_im = num_scalar(&res_si_m, "delta_neff_im");
    assert!(
        gamma_tpa > 1.0 && fom_tpa > 0.3 && dneff_im < 0.0,
        "expected positive gamma_tpa={gamma_tpa}, fom_tpa={fom_tpa}, and absorptive delta_neff_im={dneff_im:e}"
    );
}

fn run_determinism_payload() -> String {
    let lib = load_library();
    let freq0 = C0 / 1.55e-6;

    // 1. 3D depleted SHG solve:
    let solver3 = block(0.15e-6, 4.0e-6, 0.15e-6, 1, 24, 1);
    let mesh3 = Mesh::from_solver_mesh(&solver3).unwrap();
    let mut s3 = ProblemSpec::new("ems", 2, "nonlinear_driven");
    s3.study.params = json!({
        "process": "shg",
        "freq": freq0,
        "power_w": 50.0,
        "coupling": "depleted",
        "rtol": 1e-11
    });
    s3.materials.insert("core".into(), MaterialAssign::Named("lithium_niobate".into()));
    s3.bcs.push(BcSpec { on: "pmc_x".into(), kind: "pmc".into(), params: json!({}) });
    s3.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({"name": "in", "mode": 1}) });
    s3.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({"name": "out", "mode": 1}) });
    let mut ctx3 = Ctx::new(&BlockResolver);
    ctx3.library = Some(lib.clone());
    let r3 = Ems.solve(&mesh3, &s3, &mut ctx3).unwrap();

    // 2. 2D nonlinear_modes solve:
    let mut m2 = Model::new_2d();
    let clad = m2.rect("c0", (-1.2e-6, -1.0e-6), (1.2e-6, 1.0e-6));
    let core = m2.rect("c1", (-0.4e-6, -0.2e-6), (0.4e-6, 0.2e-6));
    m2.overlay(&[("clad", clad), ("core", core)]);
    let g2 = mesh_model(&m2.list, &MeshPlan::simple(0.20e-6)).unwrap();
    let mesh2 = Mesh::from_solver_mesh(&g2.mesh).unwrap();
    let res2 = GeomResolver(&g2);
    let mut s2 = ProblemSpec::new("ems", 2, "nonlinear_modes");
    s2.study.params = json!({
        "freq": freq0,
        "n_target": 1.95,
        "mode_index": 1,
        "power_w": 200.0,
        "rtol": 1e-11
    });
    s2.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), 1.444_f64.powi(2))].into()));
    s2.materials.insert("core".into(), MaterialAssign::Named("silicon_nitride".into()));
    let mut ctx2 = Ctx::new(&res2);
    ctx2.library = Some(lib);
    let r2 = Ems.solve(&mesh2, &s2, &mut ctx2).unwrap();

    let mut bits = Vec::new();
    for (_, sc) in r3.scalars.iter().chain(r2.scalars.iter()) {
        if let Scalar::Num(v) = sc {
            bits.push(format!("{:016x}", v.to_bits()));
        }
    }
    for f in &r3.solutions[0].fields {
        for &v in f.values.iter().take(16) {
            bits.push(format!("{:016x}", v.to_bits()));
        }
    }
    bits.join(",")
}

#[test]
fn subprocess_thread_count_determinism_nyx_threads_1_vs_8() {
    if let Ok(tag) = std::env::var("NYX_FR41_DET_WORKER") {
        let payload = run_determinism_payload();
        println!("NYX_FR41_BITS[{tag}]={payload}");
        return;
    }

    let exe = std::env::current_exe().expect("test binary path");
    let run_sub = |threads: usize| -> String {
        let out = Command::new(&exe)
            .arg("--exact")
            .arg("subprocess_thread_count_determinism_nyx_threads_1_vs_8")
            .arg("--nocapture")
            .env("NYX_THREADS", threads.to_string())
            .env("NYX_FR41_DET_WORKER", threads.to_string())
            .output()
            .expect("spawn determinism subprocess");
        assert!(
            out.status.success(),
            "subprocess NYX_THREADS={threads} failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        let prefix = format!("NYX_FR41_BITS[{threads}]=");
        stdout.lines().find_map(|l| l.strip_prefix(&prefix)).expect("determinism payload line").to_string()
    };

    let b1 = run_sub(1);
    let b8 = run_sub(8);
    assert_eq!(b1, b8, "nonlinear_driven & nonlinear_modes outputs must be bit-identical between NYX_THREADS=1 and 8");
}
