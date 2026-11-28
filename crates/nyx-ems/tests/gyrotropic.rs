// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of FR-34 (ADR 136): Full 3x3 non-symmetric and Hermitian gyrotropic electromagnetic
//! tensors across material specification, Bunge Euler axial-vector rotation, 2D P2 waveguide Faraday
//! rotation, Ce:YIG/Si non-reciprocal phase shift (`Delta beta_NR`), biorthogonal left eigenvectors,
//! 3D driven magneto-optic isolators/circulators (`|S_21| != |S_12|`, `||S^H S - I||_inf < 1e-6`),
//! general-LU solver dispatch (`lblt` rejection), and multi-thread determinism (`NYX_THREADS=1` vs `8`).

#![allow(clippy::needless_range_loop)]

mod common;

use femcore::linalg::csr::c;
use femcore::material::{
    axial_from_skew, rotate_axial_vector, rotate_tensor, rotation, skew_from_axial, Anisotropy, Gyrotropy,
    MaterialSpec, Tensor,
};
use femcore::mesh::SimplexMesh;
use femcore::optical::{cell_materials, is_non_symmetric_cells};
use nyx_ems::{
    assemble_modes_complex_with_bcs, biorthogonal_poynting_overlap, solve_mode_pairs_complex, solve_modes_complex, Ems,
    Problem, Study, C0,
};
use nyx_physics::{
    BcSpec, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, PropValue,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::f64::consts::PI;

fn run_on_pool<R: Send>(n_threads: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap().install(f)
}

/// Builds a 1D-layered 2D strip mesh on `[0, wx] x [z_breaks[0], z_breaks.last()]` where each vertical
/// layer `k` (`z in [z_breaks[k], z_breaks[k+1]]`) has domain tag `(k + 1)` and `nz_per_layer[k]`
/// subdivisions.
/// Boundary labels:
/// - `1`: bottom (`z = z_min`)
/// - `2`: top (`z = z_max`)
/// - `3`: left (`x = 0`)
/// - `4`: right (`x = wx`)
fn layered_strip_2d(wx: f64, nx: usize, z_breaks: &[f64], nz_per_layer: &[usize]) -> SimplexMesh {
    assert_eq!(z_breaks.len(), nz_per_layer.len() + 1);
    let mut z_coords = Vec::new();
    let mut z_layer_dom = Vec::new();
    for (layer_idx, &nz) in nz_per_layer.iter().enumerate() {
        let z0 = z_breaks[layer_idx];
        let z1 = z_breaks[layer_idx + 1];
        for j in 0..nz {
            z_coords.push(z0 + (z1 - z0) * (j as f64) / (nz as f64));
            z_layer_dom.push((layer_idx + 1) as i32);
        }
    }
    z_coords.push(*z_breaks.last().unwrap());
    let nz_tot = z_layer_dom.len();
    let vid = |ix: usize, iz: usize| -> u32 { (iz * (nx + 1) + ix) as u32 };

    let mut points = Vec::with_capacity(2 * (nx + 1) * (nz_tot + 1));
    for &z in &z_coords {
        for ix in 0..=nx {
            let x = wx * (ix as f64) / (nx as f64);
            points.push(x);
            points.push(z);
        }
    }

    let mut cells = Vec::with_capacity(6 * nx * nz_tot);
    let mut cell_labels = Vec::with_capacity(2 * nx * nz_tot);
    for iz in 0..nz_tot {
        let dom = z_layer_dom[iz];
        for ix in 0..nx {
            let v00 = vid(ix, iz);
            let v10 = vid(ix + 1, iz);
            let v01 = vid(ix, iz + 1);
            let v11 = vid(ix + 1, iz + 1);
            cells.extend_from_slice(&[v00, v10, v11]);
            cell_labels.push(dom);
            cells.extend_from_slice(&[v00, v11, v01]);
            cell_labels.push(dom);
        }
    }

    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for ix in 0..nx {
        facets.extend_from_slice(&[vid(ix, 0), vid(ix + 1, 0)]);
        facet_labels.push(1); // bottom
        facets.extend_from_slice(&[vid(ix, nz_tot), vid(ix + 1, nz_tot)]);
        facet_labels.push(2); // top
    }
    for iz in 0..nz_tot {
        facets.extend_from_slice(&[vid(0, iz), vid(0, iz + 1)]);
        facet_labels.push(3); // left
        facets.extend_from_slice(&[vid(nx, iz), vid(nx, iz + 1)]);
        facet_labels.push(4); // right
    }

    SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels)
}

#[test]
fn bunge_euler_axial_vector_rotation_and_toml_gyrotropy_models() {
    // 1. Axial vector rotation equivalence: R * [g]_x * R^T == [det(R) * R * g]_x
    let g_local = [0.018, -0.007, 0.031];
    let skew_local = skew_from_axial(&g_local);
    assert_eq!(axial_from_skew(&skew_local), g_local);

    for euler in [[25.0, 40.0, -15.0], [90.0, 90.0, 0.0], [135.0, 62.5, 210.0]] {
        let r = rotation(euler);
        let g_rot = rotate_axial_vector(&g_local, &r);
        let skew_from_rot_g = skew_from_axial(&g_rot);
        let skew_rot_mat = rotate_tensor(&skew_local, &r);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (skew_from_rot_g[a][b] - skew_rot_mat[a][b]).abs() < 1e-14,
                    "axial skew rotation mismatch at ({a},{b}) for euler={euler:?}: {} vs {}",
                    skew_from_rot_g[a][b],
                    skew_rot_mat[a][b]
                );
            }
        }
        // Under an improper reflection R_improper = -R (det = -1), an axial vector stays invariant
        // compared to a polar vector because det(-I) * (-I) = +I.
        let mut r_imp = r;
        for row in &mut r_imp {
            for v in row {
                *v = -*v;
            }
        }
        let g_imp = rotate_axial_vector(&g_local, &r_imp);
        for k in 0..3 {
            assert!((g_imp[k] - g_rot[k]).abs() < 1e-14);
        }
    }

    // 2. TOML `[optical.gyrotropy]` constant, Lorentz, and tabulated parsing + round-trip
    let toml_lorentz = r#"
format = 2
name = "CeYIG_Lorentz"

[optical]
wavelength_unit = "um"
range = [1.1, 2.0]
model = "constant"
n = 2.22

[optical.gyrotropy]
model = "lorentz"
frequency_unit = "eV"
axis = [1.0, 0.0, 0.0]
g_inf = 0.001
plasma_frequency = "0.25 eV"
oscillators = [
    [0.18, 1.05, 0.04]
]
"#;
    let spec_lor = MaterialSpec::from_toml(toml_lorentz).expect("valid Lorentz gyrotropy TOML");
    assert!(spec_lor.is_non_symmetric());
    let opt = spec_lor.optical.as_ref().unwrap();
    assert!(matches!(opt.gyrotropy, Some(Gyrotropy::Lorentz { .. })));
    let eps = opt.permittivity_global(1.55e-6, None, [0.0, 0.0, 0.0]).unwrap();
    // With axis = +x, eps_yz = -i g_x, eps_zy = +i g_x:
    assert!((eps[0][0].re - 2.22 * 2.22).abs() < 1e-12);
    assert!((eps[1][2].im + eps[2][1].im).abs() < 1e-14);
    assert!(eps[2][1].im.abs() > 1e-3, "expected non-zero gyrotropic skew entry eps_zy, got {}", eps[2][1].im);
    // Round-trip via `to_toml`
    let dumped = spec_lor.to_toml();
    let spec_rt = MaterialSpec::from_toml(&dumped).expect("round-trip Lorentz gyrotropy TOML");
    let eps2 = spec_rt.optical.as_ref().unwrap().permittivity_global(1.55e-6, None, [0.0, 0.0, 0.0]).unwrap();
    assert!((eps[1][2] - eps2[1][2]).norm() < 1e-12);

    let toml_tab = r#"
format = 2
name = "CeYIG_Tab"

[optical]
wavelength_unit = "um"
range = [1.3, 1.7]
model = "constant"
n = 2.22

[optical.gyrotropy]
model = "tabulated"
data = [
    [1.30, 0.0, 0.015, 0.0],
    [1.55, 0.0, 0.009, 0.0],
    [1.70, 0.0, 0.006, 0.0],
]
"#;
    let spec_tab = MaterialSpec::from_toml(toml_tab).expect("valid Tabulated gyrotropy TOML");
    assert!(spec_tab.is_non_symmetric());
    let g_155 = spec_tab.optical.as_ref().unwrap().gyrotropy_global(1.55e-6, [0.0, 0.0, 0.0]).unwrap().unwrap();
    assert!((g_155[1].re - 0.009).abs() < 1e-14 && g_155[1].im == 0.0);
}

#[test]
fn faraday_rotation_p2_matches_analytical_closed_form() {
    let lam0 = 1.55e-6;
    let freq = C0 / lam0;
    let k0 = 2.0 * PI / lam0;
    let eps_r = 4.84; // n0 = 2.20

    // Part A: All-orders exact closed-form Voigt dispersion in a magnetized slab (`[0, wx] x [0, hz]`,
    // PEC on z=0, hz (labels 1, 2); PMC on x=0, wx (labels 3, 4)) with transverse gyrotropy `g = (g_x, 0, 0)`.
    // Exact analytical dispersion for TM_1: beta_exact = sqrt(k0^2 * (eps_r^2 - g_x^2) / eps_r - (pi / hz)^2).
    let (wx, hz) = (0.2e-6, 0.8e-6);
    let slab_mesh = layered_strip_2d(wx, 2, &[0.0, hz], &[16]);
    let gx_voigt = 0.35;
    let eps_voigt = (eps_r * eps_r - gx_voigt * gx_voigt) / eps_r;
    let beta_voigt_exact = (k0 * k0 * eps_voigt - (PI / hz).powi(2)).sqrt();
    let neff_voigt_exact = beta_voigt_exact / k0;

    let mut prob_voigt =
        Problem::new(Study::Modes { freq, n_target: neff_voigt_exact, evanescent: false, n_modes: 2 }, 2);
    prob_voigt.materials.set(1, "eps_r", eps_r);
    prob_voigt.materials.set_gyrotropy(1, [gx_voigt, 0.0, 0.0]);
    prob_voigt.pmc = vec![3, 4];
    let (eps_v, mu_v) = cell_materials(&slab_mesh, &prob_voigt.materials).unwrap();
    assert!(is_non_symmetric_cells(&eps_v, &mu_v));

    let res_voigt = common::solve(&slab_mesh, &prob_voigt).unwrap();
    let neff_voigt_fem = res_voigt
        .solutions
        .iter()
        .map(|s| common::neff(s).0)
        .min_by(|a, b| ((a - neff_voigt_exact).abs()).partial_cmp(&(b - neff_voigt_exact).abs()).unwrap())
        .unwrap();
    assert!(
        (neff_voigt_fem - neff_voigt_exact).abs() < 1e-5,
        "Voigt P2 neff={neff_voigt_fem:.8} vs exact={neff_voigt_exact:.8}"
    );

    // Part B: Longitudinal Faraday rotation (`g = (0, g_y, 0)`) in a square PEC waveguide `[0, a] x [0, a]`
    // with P2 Nedelec elements.
    // The degenerate TE_10 / TE_01 doublet splits into RCP and LCP hybrid modes with exact analytical
    // Faraday rotation angle over length L:
    //   theta_F = 1/2 * (beta_+ - beta_-) * L = (4 * k0^2 * g_y * L) / (pi^2 * beta_0) + O(g_y^3),
    // where beta_0 = sqrt(k0^2 * eps_r - (pi / a)^2).
    let a = 1.2e-6;
    let gy = 5.0e-3;
    let length = 50.0e-6;
    let beta0 = (k0 * k0 * eps_r - (PI / a).powi(2)).sqrt();
    let n0_te = beta0 / k0;
    let theta_f_exact = (4.0 * k0 * k0 * gy * length) / (PI * PI * beta0);

    let sq_mesh = common::rect(a, a, 14, 14);
    let mut prob_faraday = Problem::new(Study::Modes { freq, n_target: n0_te, evanescent: false, n_modes: 2 }, 2);
    prob_faraday.materials.aniso.insert(
        1,
        Anisotropy {
            permittivity: Some(Tensor::Diagonal { x: eps_r, y: eps_r, z: eps_r }),
            gyrotropy: Some([0.0, gy, 0.0]),
            ..Default::default()
        },
    );
    let (eps_f, mu_f) = cell_materials(&sq_mesh, &prob_faraday.materials).unwrap();
    assert!(is_non_symmetric_cells(&eps_f, &mu_f));
    let modes_f = solve_modes_complex(&sq_mesh, freq, &eps_f, &mu_f, c(n0_te, 0.0), 2, 2).unwrap();
    assert!(modes_f.is_non_symmetric);
    assert_eq!(modes_f.neff.len(), 2);

    let n_plus = modes_f.neff[0].re.max(modes_f.neff[1].re);
    let n_minus = modes_f.neff[0].re.min(modes_f.neff[1].re);
    let theta_f_fem = 0.5 * k0 * (n_plus - n_minus) * length;
    let err_rad = (theta_f_fem - theta_f_exact).abs();
    assert!(
        err_rad < 1e-4,
        "Faraday rotation P2 error {err_rad:.3e} rad exceeds 1e-4 rad (FEM={theta_f_fem:.6} rad, exact={theta_f_exact:.6} rad)"
    );
}

/// Analytical 3-layer gyrotropic TM slab waveguide dispersion root finder for `beta`.
/// Substrate (`z in [-d_sub, 0]`, permittivity `eps_sub`, PEC at `-d_sub`),
/// Magneto-optic core (`z in [0, h]`, permittivity `eps_c`, transverse gyrotropy `g_x`),
/// Cladding (`z in [h, h + d_clad]`, permittivity `eps_clad`, PEC at `h + d_clad`).
#[allow(clippy::too_many_arguments)]
fn analytical_3layer_tm_gyrotropic_beta(
    k0: f64,
    eps_sub: f64,
    d_sub: f64,
    eps_c: f64,
    h: f64,
    gx: f64,
    eps_clad: f64,
    d_clad: f64,
) -> f64 {
    let eps_v = (eps_c * eps_c - gx * gx) / eps_c;
    let f = |beta: f64| -> f64 {
        let kappa = (k0 * k0 * eps_v - beta * beta).max(1e-30).sqrt();
        let g_sub = (beta * beta - k0 * k0 * eps_sub).max(1e-30).sqrt();
        let g_clad = (beta * beta - k0 * k0 * eps_clad).max(1e-30).sqrt();
        let phi_sub = g_sub * (g_sub * d_sub).tanh() / eps_sub;
        let phi_clad = g_clad * (g_clad * d_clad).tanh() / eps_clad;
        let num = kappa * eps_v * (phi_sub + phi_clad);
        let den = kappa * kappa - (eps_v * phi_sub + beta * gx / eps_c) * (eps_v * phi_clad - beta * gx / eps_c);
        (kappa * h).sin() * den - (kappa * h).cos() * num
    };

    // Fundamental TM_0 root lies in `(k0 * sqrt(eps_sub), k0 * sqrt(eps_v))` with `kappa * h in (0, pi)`.
    let b_min = k0 * eps_sub.sqrt() * 1.001;
    let b_max = (k0 * k0 * eps_v - (0.1 / h).powi(2)).sqrt();
    let n_scan = 400;
    let mut bracket = None;
    let mut prev_b = b_max;
    let mut prev_f = f(prev_b);
    for i in 1..=n_scan {
        let b = b_max - (b_max - b_min) * (i as f64) / (n_scan as f64);
        let fb = f(b);
        if prev_f * fb <= 0.0 {
            bracket = Some((b, prev_b));
            break;
        }
        prev_b = b;
        prev_f = fb;
    }
    let (mut lo, mut hi) = bracket.expect("fundamental TM_0 root bracketed");
    let mut flo = f(lo);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let fmid = f(mid);
        if flo * fmid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = fmid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn ce_yig_on_si_waveguide_non_reciprocal_phase_shift_and_biorthogonality() {
    // Asymmetric 3-layer planar slab waveguide at lambda = 1.55 um:
    // - Substrate (domain 1): SiO2, n_sub = 1.444, d_sub = 1.20 um
    // - Core (domain 2): Ce:YIG/Si hybrid slab, n_core = 2.22, h = 0.45 um, transverse gyrotropy g_x = 0.025
    // - Cladding (domain 3): Air, n_clad = 1.0, d_clad = 1.20 um
    // With PEC on top/bottom (labels 1, 2) and PMC on left/right (labels 3, 4),
    // the TM mode (E_z, E_y, H_x) is strictly 1D-exact in z.
    let lam0 = 1.55e-6;
    let freq = C0 / lam0;
    let k0 = 2.0 * PI / lam0;
    let eps_sub = 1.444_f64.powi(2);
    let eps_core = 2.22_f64.powi(2);
    let eps_clad = 1.0;
    let (d_sub, h_core, d_clad) = (1.20e-6, 0.45e-6, 1.20e-6);
    let gx = 0.025;

    let beta_plus_exact =
        analytical_3layer_tm_gyrotropic_beta(k0, eps_sub, d_sub, eps_core, h_core, gx, eps_clad, d_clad);
    let beta_minus_exact =
        analytical_3layer_tm_gyrotropic_beta(k0, eps_sub, d_sub, eps_core, h_core, -gx, eps_clad, d_clad);
    let delta_beta_exact = beta_plus_exact - beta_minus_exact;
    assert!(delta_beta_exact.abs() > 1e3, "expected significant non-reciprocal shift, got {delta_beta_exact}");

    let n_tm_guess = 0.5 * (beta_plus_exact + beta_minus_exact) / k0;
    let mesh = layered_strip_2d(0.15e-6, 1, &[-d_sub, 0.0, h_core, h_core + d_clad], &[18, 16, 18]);

    let build_prob = |gx_val: f64, n_target: f64| -> Problem {
        let mut p = Problem::new(Study::Modes { freq, n_target, evanescent: false, n_modes: 2 }, 2);
        p.materials.set(1, "eps_r", eps_sub);
        p.materials.set(2, "eps_r", eps_core);
        p.materials.set_gyrotropy(2, [gx_val, 0.0, 0.0]);
        p.materials.set(3, "eps_r", eps_clad);
        p.pmc = vec![3, 4];
        p
    };

    // 1. Solve forward (+beta, n_target > 0) and backward (-beta, n_target < 0) on the SAME problem (+gx)!
    let p_fwd = build_prob(gx, n_tm_guess);
    let (eps, mu) = cell_materials(&mesh, &p_fwd.materials).unwrap();
    let asm = assemble_modes_complex_with_bcs(&mesh, freq, &eps, &mu, 2, &[], &[], &[3, 4]).unwrap();
    assert!(asm.is_non_symmetric);
    let modes_fwd = solve_mode_pairs_complex(&asm, c(n_tm_guess, 0.0), 2).unwrap();
    let modes_bwd = solve_mode_pairs_complex(&asm, c(-n_tm_guess, 0.0), 2).unwrap();
    assert!(modes_fwd.is_non_symmetric && modes_bwd.is_non_symmetric);

    let neff_fwd = modes_fwd
        .neff
        .iter()
        .map(|z| z.re)
        .min_by(|a, b| ((a * k0 - beta_plus_exact).abs()).partial_cmp(&((b * k0 - beta_plus_exact).abs())).unwrap())
        .unwrap();
    let neff_bwd_abs = modes_bwd
        .neff
        .iter()
        .map(|z| -z.re)
        .min_by(|a, b| ((a * k0 - beta_minus_exact).abs()).partial_cmp(&((b * k0 - beta_minus_exact).abs())).unwrap())
        .unwrap();

    // Also verify magnetization reversal (-gx with n_target > 0) gives the identical backward propagation constant
    let p_rev = build_prob(-gx, n_tm_guess);
    let (eps_rev, mu_rev) = cell_materials(&mesh, &p_rev.materials).unwrap();
    let asm_rev = assemble_modes_complex_with_bcs(&mesh, freq, &eps_rev, &mu_rev, 2, &[], &[], &[3, 4]).unwrap();
    let modes_rev = solve_mode_pairs_complex(&asm_rev, c(n_tm_guess, 0.0), 2).unwrap();
    let neff_rev = modes_rev
        .neff
        .iter()
        .map(|z| z.re)
        .min_by(|a, b| ((a * k0 - beta_minus_exact).abs()).partial_cmp(&((b * k0 - beta_minus_exact).abs())).unwrap())
        .unwrap();
    assert!(
        (neff_bwd_abs - neff_rev).abs() < 1e-10,
        "backward root (-beta) {neff_bwd_abs:.10} must match reversed magnetization (+beta, -gx) {neff_rev:.10}"
    );

    let delta_beta_fem = (neff_fwd - neff_bwd_abs) * k0;
    let rel_err = ((delta_beta_fem - delta_beta_exact) / delta_beta_exact).abs();
    assert!(
        rel_err < 5e-3,
        "Ce:YIG/Si Delta beta_NR relative error {:.4}% exceeds 0.5% (FEM={delta_beta_fem:.4} rad/m, exact={delta_beta_exact:.4} rad/m)",
        rel_err * 100.0
    );

    // 2. Verify true left eigenvectors & biorthogonal normalization + cross-section Poynting overlap
    for &norm in &modes_fwd.biorth_norms {
        assert!((norm - c(1.0, 0.0)).norm() < 1e-10, "expected unit biorthogonal norm y^H B x = 1, got {norm:?}");
    }
    // Check that left eigenvector is NOT trivially conj(right eigenvector) when gyrotropy is active
    let diff_conj: f64 =
        modes_fwd.et[0].iter().zip(&modes_fwd.et_left[0]).map(|(r, l)| (*l - r.conj()).norm_sqr()).sum::<f64>().sqrt();
    assert!(diff_conj > 1e-6, "non-symmetric mode must not use y = conj(x) shortcut");

    let poynting_00 = biorthogonal_poynting_overlap(&mesh, freq, &eps, &mu, 2, &modes_fwd, 0, 0);
    assert!(poynting_00.norm() > 0.0 && poynting_00.re.is_finite());

    // 3. Thread determinism (`NYX_THREADS=1` vs `8`)
    let neff_t1 = run_on_pool(1, || solve_mode_pairs_complex(&asm, c(n_tm_guess, 0.0), 2).unwrap().neff);
    let neff_t8 = run_on_pool(8, || solve_mode_pairs_complex(&asm, c(n_tm_guess, 0.0), 2).unwrap().neff);
    assert_eq!(neff_t1, neff_t8, "2D gyrotropic waveguide modes must be bit-identical across thread counts");
}

struct CirculatorResolver;

impl LabelResolver for CirculatorResolver {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        match region {
            "ferrite" => Ok(vec![1]),
            _ => Err(Error::invalid(format!("unknown domain '{region}'"))),
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        match selector {
            "port1" => Ok(vec![1]),         // x = 0
            "port2" => Ok(vec![2]),         // y = 0
            "port3" => Ok(vec![6]),         // x = lx
            "pec_plates" => Ok(vec![3, 4]), // z = 0 and z = lz ground planes
            "pmc_3port" => Ok(vec![5]),     // y = ly open magnetic wall
            "pmc_2port" => Ok(vec![5, 6]),  // y = ly and x = lx open magnetic walls
            _ => Err(Error::invalid(format!("unknown boundary '{selector}'"))),
        }
    }
}

fn ferrite_props(eps_r: f64, eps_r_im: f64, gz: f64) -> MaterialAssign {
    let mut m = BTreeMap::new();
    m.insert("eps_r".to_string(), PropValue::Scalar(eps_r));
    if eps_r_im != 0.0 {
        m.insert("eps_r_im".to_string(), PropValue::Scalar(eps_r_im));
    }
    m.insert("gyrotropy".to_string(), PropValue::Diag([0.45 * gz, -0.35 * gz, gz]));
    m.insert("mu_gyrotropy".to_string(), PropValue::Diag([0.0, 0.0, 0.35 * gz]));
    MaterialAssign::PropsTensor(m)
}

#[test]
fn driven_magneto_optic_isolator_and_circulator_non_reciprocity_and_unitarity() {
    // 3D parallel-plate stripline junction `[0, lx] x [0, ly] x [0, lz]` with PEC ground planes at
    // `z = 0` (label 3) and `z = lz` (label 4), loaded with a gyrotropic / gyromagnetic medium!
    // - Port 1 sits on face `x = 0` (label 1), bridging `z = 0` to `z = lz` (`direction = [0, 0, 1]`)
    // - Port 2 sits on face `y = 0` (label 2), bridging `z = 0` to `z = lz` (`direction = [0, 0, 1]`)
    // - Port 3 sits on face `x = lx` (label 6), bridging `z = 0` to `z = lz` (`direction = [0, 0, 1]`)
    let (lx, ly, lz) = (0.020, 0.020, 0.008);
    let mesh3 = common::block(lx, ly, lz, 3, 3, 2);
    let ex_mesh = Mesh::from_solver_mesh(&mesh3).unwrap();
    let freq = 6.5e9;

    // 1. Verify `solver = "lblt"` is explicitly rejected when the tensor is non-symmetric
    let mut spec_lblt = ProblemSpec::new("ems", 1, "driven");
    spec_lblt.study.params = json!({
        "f_min": freq,
        "f_max": freq,
        "n_points": 1,
        "solver": "lblt",
    });
    spec_lblt.materials.insert("ferrite".into(), ferrite_props(4.5, 0.0, 1.2));
    spec_lblt.bcs = vec![
        BcSpec { on: "pec_plates".into(), kind: "pec".into(), params: json!({}) },
        BcSpec { on: "pmc_2port".into(), kind: "pmc".into(), params: json!({}) },
        BcSpec {
            on: "port1".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "P1", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
        },
        BcSpec {
            on: "port2".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "P2", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
        },
    ];
    let err_lblt = Ems.solve(&ex_mesh, &spec_lblt, &mut Ctx::new(&CirculatorResolver)).unwrap_err();
    assert!(
        err_lblt.to_string().contains("non-symmetric") || err_lblt.to_string().contains("lblt"),
        "expected diagnostic rejecting lblt on gyrotropic problem, got: {err_lblt}"
    );

    // 2. Lossless 3-port magneto-optic circulator / junction:
    //    Simultaneously exhibits amplitude non-reciprocity `|S_21| != |S_12|` AND lossless unitarity
    //    `||S^H S - I||_inf < 1e-6`!
    let build_3port_spec = |prom_tol: serde_json::Value| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({
            "f_min": freq,
            "f_max": freq * 1.02,
            "n_points": 2,
            "solver": "direct",
            "prom_tol": prom_tol,
        });
        s.materials.insert("ferrite".into(), ferrite_props(4.5, 0.0, 1.2));
        s.bcs = vec![
            BcSpec { on: "pec_plates".into(), kind: "pec".into(), params: json!({}) },
            BcSpec { on: "pmc_3port".into(), kind: "pmc".into(), params: json!({}) },
            BcSpec {
                on: "port1".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "P1", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
            },
            BcSpec {
                on: "port2".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "P2", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
            },
            BcSpec {
                on: "port3".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "P3", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
            },
        ];
        s
    };

    let res_direct =
        Ems.solve(&ex_mesh, &build_3port_spec(serde_json::Value::Null), &mut Ctx::new(&CirculatorResolver)).unwrap();
    let dp = DrivenPortData::of(&res_direct).expect("driven port data");
    assert_eq!(dp.ports.len(), 3);
    let s = &dp.s_matrix[0];
    let s21 = c(s[1][0][0], s[1][0][1]);
    let s12 = c(s[0][1][0], s[0][1][1]);
    assert!(
        (s21.norm() - s12.norm()).abs() > 1e-2,
        "expected non-reciprocal amplitude |S_21|={} != |S_12|={}",
        s21.norm(),
        s12.norm()
    );

    // Check lossless unitarity `||S^H S - I||_inf < 1e-6`
    let n_p = 3;
    let mut max_unitarity_err = 0.0_f64;
    for i in 0..n_p {
        let mut row_err = 0.0_f64;
        for j in 0..n_p {
            let mut sum = c(0.0, 0.0);
            for k in 0..n_p {
                let ski = c(s[k][i][0], s[k][i][1]);
                let skj = c(s[k][j][0], s[k][j][1]);
                sum += ski.conj() * skj;
            }
            let target = if i == j { c(1.0, 0.0) } else { c(0.0, 0.0) };
            row_err += (sum - target).norm();
        }
        max_unitarity_err = max_unitarity_err.max(row_err);
    }
    assert!(
        max_unitarity_err < 1e-6,
        "lossless gyrotropic S-matrix unitarity error ||S^H S - I||_inf = {max_unitarity_err:.3e} >= 1e-6"
    );

    // Verify PROM sweep (`PromOperator::is_symmetric() == false`) matches direct HDM solve on all non-reciprocal S_ij
    let res_prom = Ems.solve(&ex_mesh, &build_3port_spec(json!(1e-7)), &mut Ctx::new(&CirculatorResolver)).unwrap();
    let dp_prom = DrivenPortData::of(&res_prom).unwrap();
    for fi in 0..2 {
        for i in 0..n_p {
            for j in 0..n_p {
                let sd = c(dp.s_matrix[fi][i][j][0], dp.s_matrix[fi][i][j][1]);
                let sp = c(dp_prom.s_matrix[fi][i][j][0], dp_prom.s_matrix[fi][i][j][1]);
                assert!(
                    (sd - sp).norm() < 1e-5,
                    "PROM vs HDM mismatch at freq {fi} S_{}{}: {sp:?} vs {sd:?}",
                    i + 1,
                    j + 1
                );
            }
        }
    }

    // 3. 2-port magneto-optic isolator (with resonant damping `eps_r_im` + gyrotropy) showing 2-port `|S_21| != |S_12|`
    //    and lossless 2-port unitarity `||S^H S - I||_inf < 1e-6` when damping is zero.
    let build_2port = |loss_im: f64| -> DrivenPortData {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({
            "f_min": freq,
            "f_max": freq,
            "n_points": 1,
            "solver": "auto",
            "prom_tol": serde_json::Value::Null,
        });
        s.materials.insert("ferrite".into(), ferrite_props(4.5, loss_im, 1.2));
        s.bcs = vec![
            BcSpec { on: "pec_plates".into(), kind: "pec".into(), params: json!({}) },
            BcSpec { on: "pmc_2port".into(), kind: "pmc".into(), params: json!({}) },
            BcSpec {
                on: "port1".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "P1", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
            },
            BcSpec {
                on: "port2".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "P2", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }),
            },
        ];
        let res = Ems.solve(&ex_mesh, &s, &mut Ctx::new(&CirculatorResolver)).unwrap();
        DrivenPortData::of(&res).unwrap()
    };

    let dp_2port_lossless = build_2port(0.0);
    let s2_0 = &dp_2port_lossless.s_matrix[0];
    let s21_0 = c(s2_0[1][0][0], s2_0[1][0][1]);
    let s12_0 = c(s2_0[0][1][0], s2_0[0][1][1]);
    assert!((s21_0 - s12_0).norm() > 1e-2, "lossless 2-port gyrator must have S_21 != S_12");
    let mut err_2p_unit = 0.0_f64;
    for i in 0..2 {
        let mut r_err = 0.0_f64;
        for j in 0..2 {
            let mut sum = c(0.0, 0.0);
            for k in 0..2 {
                sum += c(s2_0[k][i][0], s2_0[k][i][1]).conj() * c(s2_0[k][j][0], s2_0[k][j][1]);
            }
            let target = if i == j { c(1.0, 0.0) } else { c(0.0, 0.0) };
            r_err += (sum - target).norm();
        }
        err_2p_unit = err_2p_unit.max(r_err);
    }
    assert!(err_2p_unit < 1e-6, "lossless 2-port unitarity ||S^H S - I||_inf = {err_2p_unit:.3e} >= 1e-6");

    let dp_2port_isolator = build_2port(-0.35);
    let s2_iso = &dp_2port_isolator.s_matrix[0];
    let s21_iso = c(s2_iso[1][0][0], s2_iso[1][0][1]);
    let s12_iso = c(s2_iso[0][1][0], s2_iso[0][1][1]);
    assert!(
        (s21_iso.norm() - s12_iso.norm()).abs() > 1e-3,
        "2-port damped magneto-optic isolator must have |S_21|={} != |S_12|={}",
        s21_iso.norm(),
        s12_iso.norm()
    );

    // 4. Thread determinism (`NYX_THREADS=1` vs `8`) on 3D driven non-symmetric solve
    let s_t1 = run_on_pool(1, || build_2port(0.0).s_matrix);
    let s_t8 = run_on_pool(8, || build_2port(0.0).s_matrix);
    assert_eq!(s_t1, s_t8, "3D driven gyrotropic S-matrix must be bit-identical across thread counts");
}
