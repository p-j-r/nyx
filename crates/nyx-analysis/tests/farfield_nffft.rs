// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-36 (ADR 138) validation of the Stratton-Chu near-to-far-field transformation and of the vectorial beam overlap on
//! synthetic `SolveResult`s whose per-cell P2 fields are sampled from closed-form solutions: a Hertzian dipole (directivity
//! 1.5, radiated power), a uniformly illuminated focusing aperture (Airy PSF, Strehl ratio), Gaussian beams (self and
//! waist-mismatch overlap), the LP01 fiber mode against its Gaussian approximation, and thread-count determinism.

use nyx_analysis::farfield::{
    sample_surface, CVec3, CellSampler, Equivalence, FieldSelection, Medium, Radiator, SurfaceSpec, C0, C64, EPS0, MU0,
};
use nyx_analysis::{
    beam_overlap, beam_plane_samples, solve_lp_mode, Analysis, Beam, BeamOverlapAnalysis, BeamShape, FarFieldAnalysis,
};
use nyx_physics::{CellField, Mesh, Solution, SolveResult};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Kuhn tetrahedra (6 per hex, conforming) of the hexes `(i, j, k)` of an `n` grid over `[lo, hi]` with `keep(i, j, k)`;
/// exterior faces get the label `label(centroid)` (or are left out for `None`).
fn kuhn_mesh(
    lo: [f64; 3],
    hi: [f64; 3],
    n: [usize; 3],
    keep: impl Fn(usize, usize, usize) -> bool,
    label: impl Fn([f64; 3]) -> Option<i32>,
) -> Mesh {
    let idx = |i: usize, j: usize, k: usize| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    let coord = |id: usize| {
        let i = id % (n[0] + 1);
        let j = (id / (n[0] + 1)) % (n[1] + 1);
        let k = id / ((n[0] + 1) * (n[1] + 1));
        [
            lo[0] + (hi[0] - lo[0]) * i as f64 / n[0] as f64,
            lo[1] + (hi[1] - lo[1]) * j as f64 / n[1] as f64,
            lo[2] + (hi[2] - lo[2]) * k as f64 / n[2] as f64,
        ]
    };
    let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let mut map: HashMap<usize, u32> = HashMap::new();
    let mut pts: Vec<f64> = Vec::new();
    let mut cells: Vec<u32> = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                if !keep(i, j, k) {
                    continue;
                }
                for perm in perms {
                    let mut v = [0usize; 3];
                    let mut tet = vec![idx(i, j, k)];
                    for ax in perm {
                        v[ax] = 1;
                        tet.push(idx(i + v[0], j + v[1], k + v[2]));
                    }
                    for g in tet {
                        let next = map.len() as u32;
                        let id = *map.entry(g).or_insert_with(|| {
                            pts.extend_from_slice(&coord(g));
                            next
                        });
                        cells.push(id);
                    }
                }
            }
        }
    }
    let mut faces: HashMap<[u32; 3], u32> = HashMap::new();
    for c in cells.chunks(4) {
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&q| q != skip).map(|q| c[q]).collect();
            f.sort_unstable();
            *faces.entry([f[0], f[1], f[2]]).or_insert(0) += 1;
        }
    }
    let mut ext: Vec<[u32; 3]> = faces.into_iter().filter(|(_, n)| *n == 1).map(|(f, _)| f).collect();
    ext.sort_unstable();
    let (mut facets, mut flabels) = (Vec::new(), Vec::new());
    for f in ext {
        let p = |v: u32| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]];
        let c = [0, 1, 2].map(|q| (p(f[0])[q] + p(f[1])[q] + p(f[2])[q]) / 3.0);
        if let Some(l) = label(c) {
            facets.extend_from_slice(&f);
            flabels.push(l);
        }
    }
    let nc = cells.len() / 4;
    Mesh::new(3, pts, cells, vec![1; nc], facets, flabels)
}

/// A driven-like result whose P2 per-cell fields `E{prefix}*`, `H{prefix}*` are `f` at the cell nodes.
fn result_from(mesh: &Mesh, freq: f64, prefix: &str, f: impl Fn([f64; 3]) -> (CVec3, CVec3)) -> SolveResult {
    let disc = nyx_physics::discretise(mesh).unwrap();
    let lam = femcore::fields::local_bary(3, 2);
    let mut re = vec![Vec::new(); 6];
    let mut im = vec![Vec::new(); 6];
    for c in 0..disc.n_cells() {
        let v = disc.cell(c);
        let p: Vec<[f64; 3]> = v
            .iter()
            .map(|&i| [disc.point(i as usize)[0], disc.point(i as usize)[1], disc.point(i as usize)[2]])
            .collect();
        for l in &lam {
            let x = [0, 1, 2].map(|q| (0..4).map(|a| l[a] * p[a][q]).sum::<f64>());
            let (e, h) = f(x);
            for q in 0..3 {
                re[q].push(e[q].re);
                im[q].push(e[q].im);
                re[3 + q].push(h[q].re);
                im[3 + q].push(h[q].im);
            }
        }
    }
    let names = ["E", "E", "E", "H", "H", "H"];
    let comps = ["x", "y", "z", "x", "y", "z"];
    let cell_fields = (0..6)
        .map(|q| CellField {
            name: format!("{}{prefix}{}", names[q], comps[q]),
            unit: if q < 3 { "V/m".into() } else { "A/m".into() },
            nloc: 10,
            re: re[q].clone(),
            im: im[q].clone(),
        })
        .collect();
    SolveResult {
        physics: "ems".into(),
        study: "driven".into(),
        order: 2,
        dim: 3,
        solutions: vec![Solution {
            label: "f".into(),
            fields: vec![],
            value: Some(freq),
            eigenvalue: None,
            cell_fields,
        }],
        ..Default::default()
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Exact fields of an oscillating electric dipole `p` (C m) in vacuum, `exp(+j w t)`, `exp(-j k r)`.
fn dipole(x: [f64; 3], p: [f64; 3], k: f64) -> (CVec3, CVec3) {
    let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    if r < 1e-12 {
        return ([C64::ZERO; 3], [C64::ZERO; 3]);
    }
    let n = [x[0] / r, x[1] / r, x[2] / r];
    let np = cross(n, p);
    let npn = cross(np, n);
    let ndp = n[0] * p[0] + n[1] * p[1] + n[2] * p[2];
    let q = [0, 1, 2].map(|i| 3.0 * n[i] * ndp - p[i]);
    let ph = C64::cis(-k * r);
    let hfac = ph * C64::new(1.0, -1.0 / (k * r)) * (C0 * k * k / (4.0 * PI * r));
    let far = ph * (k * k / r);
    let near = ph * C64::new(1.0 / (r * r * r), k / (r * r));
    let e = [0, 1, 2].map(|i| (far * npn[i] + near * q[i]) * (1.0 / (4.0 * PI * EPS0)));
    let h = [0, 1, 2].map(|i| hfac * np[i]);
    (e, h)
}

fn num(v: &Value, key: &str) -> f64 {
    v.get(key).and_then(Value::as_f64).unwrap_or_else(|| panic!("missing {key} in {v}"))
}

fn shell_mesh(n: usize, half: f64, shell: usize) -> Mesh {
    let keep = move |i: usize, j: usize, k: usize| {
        let m = i.min(j).min(k);
        let x = i.max(j).max(k);
        m < shell || x >= n - shell
    };
    let label = move |c: [f64; 3]| {
        let outer = c.iter().any(|v| (v.abs() - half).abs() < 1e-9 * half);
        Some(if outer { 1 } else { 2 })
    };
    kuhn_mesh([-half; 3], [half; 3], [n; 3], keep, label)
}

#[test]
fn hertzian_dipole_directivity_power_and_near_field() {
    let freq = C0; // lambda = 1 m
    let k = 2.0 * PI;
    let p = [0.0, 0.0, 1e-9];
    let (n, half) = (16usize, 0.6);
    let h = 2.0 * half / n as f64;
    let mesh = shell_mesh(n, half, 2);
    let result = result_from(&mesh, freq, "", |x| dipole(x, p, k));
    let z0 = MU0 * C0;
    let p_exact = C0 * C0 * z0 * k.powi(4) * p[2] * p[2] / (12.0 * PI);

    // (a) the labelled outer boundary facets
    let params = json!({
        "surface": { "kind": "labels", "labels": [1] }, "n_theta": 37, "n_phi": 24,
        "kspace": { "n": 11 }, "stratified": { "z0": 0.0, "eps2": 1.0, "side": "below" }
    });
    let out = FarFieldAnalysis.run_with_params(&mesh, &result, &params).unwrap();
    let (dmax, pff, pfl) = (num(&out, "D_max"), num(&out, "P_rad_farfield_w"), num(&out, "P_rad_flux_w"));
    eprintln!(
        "dipole labels: D_max = {dmax:.6} at theta {} | P_ff = {pff:.6e} P_flux = {pfl:.6e} P_exact = {p_exact:.6e} ({:.3e}, {:.3e})",
        out["D_max_theta_deg"],
        pff / p_exact - 1.0,
        pfl / p_exact - 1.0
    );
    assert!((dmax - 1.5).abs() < 1e-3, "D_max {dmax}");
    assert!((num(&out, "D_max_theta_deg") - 90.0).abs() < 1e-9);
    assert!((pff / p_exact - 1.0).abs() < 1e-2, "P_ff {pff} vs {p_exact}");
    assert!((pfl / p_exact - 1.0).abs() < 1e-2, "P_flux {pfl} vs {p_exact}");
    // stratified projection into an identical half-space carries half the dipole power
    let ptx = num(&out["stratified"], "P_transmitted_w");
    eprintln!("dipole stratified (eps2 = 1): P_tx / (P_ff / 2) - 1 = {:.3e}", ptx / (0.5 * pff) - 1.0);
    assert!((ptx / (0.5 * pff) - 1.0).abs() < 1e-2);
    assert_eq!(out["kspace"]["U_w_per_sr"].as_array().map(Vec::len), Some(11));

    // (b) an axis-aligned box on interior mesh faces (two cells per face)
    let b = half - h;
    let params = json!({ "surface": { "kind": "box", "lo": [-b, -b, -b], "hi": [b, b, b], "sampling": "faces" }, "n_theta": 37, "n_phi": 24 });
    let outb = FarFieldAnalysis.run_with_params(&mesh, &result, &params).unwrap();
    let cov = num(&outb["surface"], "coverage");
    let (dmax_b, pff_b) = (num(&outb, "D_max"), num(&outb, "P_rad_farfield_w"));
    eprintln!("dipole box: coverage {cov:.12} D_max = {dmax_b:.6} P_ff/P_exact - 1 = {:.3e}", pff_b / p_exact - 1.0);
    assert!((cov - 1.0).abs() < 1e-9);
    assert!((dmax_b - 1.5).abs() < 1e-3);
    assert!((pff_b / p_exact - 1.0).abs() < 1e-2);

    // (c) the near-zone dyadic kernel reproduces the exact field outside the surface
    let sel = FieldSelection::with_prefix(0, "");
    let sampler = CellSampler::new(&mesh, &result, &sel).unwrap();
    let surf = sample_surface(&sampler, &SurfaceSpec::Labels { labels: vec![1], inside_domains: vec![] }).unwrap();
    let medium = Medium::new(1.0, 1.0, freq).unwrap();
    let rad = Radiator::new(&surf.samples, &medium, Equivalence::Love, None, 4);
    for x in [[1.5, 0.4, 0.9], [0.2, -2.3, -0.7]] {
        let (e, hh) = rad.near(x);
        let (ex, hx) = dipole(x, p, k);
        let err = |a: &CVec3, b: &CVec3| {
            let d: f64 = (0..3).map(|i| (a[i] - b[i]).norm_sqr()).sum();
            let s: f64 = (0..3).map(|i| b[i].norm_sqr()).sum();
            (d / s).sqrt()
        };
        let (ee, eh) = (err(&e, &ex), err(&hh, &hx));
        eprintln!("dipole near field at {x:?}: rel err E {ee:.3e} H {eh:.3e}");
        assert!(ee < 1e-2 && eh < 1e-2);
    }
}

#[test]
fn focusing_aperture_psf_matches_airy() {
    let freq = C0; // lambda = 1 m
    let k = 2.0 * PI;
    let eta0 = MU0 * C0;
    let (na, f) = (0.2f64, 15.0f64);
    let a = f * na.asin().tan();
    let half = 3.3;
    let n = 53usize;
    let hz = 2.0 * half / n as f64;
    let mesh = kuhn_mesh([-half, -half, -hz], [half, half, hz], [n, n, 2], |_, _, _| true, |_| Some(1));
    let result = result_from(&mesh, freq, "", |x| {
        let rho = x[0].hypot(x[1]);
        if rho > a {
            return ([C64::ZERO; 3], [C64::ZERO; 3]);
        }
        let d = [-x[0], -x[1], f - x[2]];
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let s = [d[0] / l, d[1] / l, d[2] / l];
        let ev = [1.0 - s[0] * s[0], -s[0] * s[1], -s[0] * s[2]];
        let hv = cross(s, ev);
        let ph = C64::cis(k * l);
        ([0, 1, 2].map(|i| ph * ev[i]), [0, 1, 2].map(|i| ph * (hv[i] / eta0)))
    });
    let r0 = 0.61 / na;
    let params = json!({
        "surface": { "kind": "plane", "axis": "z", "position": 0.0, "lo": [-half, -half], "hi": [half, half], "sampling": "faces" },
        "n_theta": 5, "n_phi": 4, "power_quadrature": { "n_theta": 16, "n_phi": 16 },
        "focal": { "axis": "z", "position": f, "u": [-2.0 * r0, 2.0 * r0, 121], "na": na }
    });
    let out = FarFieldAnalysis.run_with_params(&mesh, &result, &params).unwrap();
    let fo = &out["focal"];
    let (strehl, z_meas, z_airy) = (num(fo, "strehl"), num(fo, "measured_first_zero_m"), num(fo, "airy_first_zero_m"));
    eprintln!(
        "PSF: equivalence {} coverage {} strehl = {strehl:.4}, first zero {z_meas:.4} vs Airy {z_airy:.4} ({:.3e}), peak at u = {}",
        out["equivalence"],
        out["surface"]["coverage"],
        z_meas / z_airy - 1.0,
        fo["peak_u_m"]
    );
    assert!((strehl - 1.0).abs() < 0.05, "strehl {strehl}");
    assert!((z_meas / z_airy - 1.0).abs() < 0.03, "first zero {z_meas} vs {z_airy}");
    assert!(num(fo, "peak_u_m").abs() < 0.11);
}

fn gaussian_slab(w1: f64, l: f64, n: usize) -> (Mesh, SolveResult, Medium) {
    let freq = C0;
    let medium = Medium::new(1.0, 1.0, freq).unwrap();
    let hz = 2.0 * l / n as f64;
    let mesh = kuhn_mesh([-l, -l, -hz], [l, l, hz], [n, n, 2], |_, _, _| true, |_| Some(1));
    let beam = Beam::new(BeamShape::Gaussian { waist_m: w1 }, [0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let result = result_from(&mesh, freq, "", |x| beam.fields(&medium, x));
    (mesh, result, medium)
}

#[test]
fn gaussian_self_and_mismatch_overlap() {
    let (w1, w2, l) = (2.0, 3.0, 9.0);
    let (mesh, result, _) = gaussian_slab(w1, l, 48);
    let params = json!({
        "plane": { "axis": "z", "position": 0.0, "lo": [-l, -l], "hi": [l, l], "n": [120, 120] },
        "beams": [ { "kind": "gaussian", "waist_m": w1 }, { "kind": "gaussian", "waist_m": w2 }, { "kind": "thin_lens", "na": 0.3 } ]
    });
    let out = BeamOverlapAnalysis.run_with_params(&mesh, &result, &params).unwrap();
    let beams = out["beams"].as_array().unwrap();
    let (e_self, e_mis, e_lens) = (num(&beams[0], "eta"), num(&beams[1], "eta"), num(&beams[2], "eta"));
    let want = (2.0 * w1 * w2 / (w1 * w1 + w2 * w2)).powi(2);
    eprintln!(
        "gaussian overlap: self {e_self:.6} (E-only {:.6}), mismatch {e_mis:.6} vs {want:.6} ({:.3e}), thin lens NA 0.3: {e_lens:.4}",
        num(&beams[0], "eta_e"),
        e_mis - want
    );
    assert!(e_self > 0.999, "self overlap {e_self}");
    assert!((e_mis - want).abs() < 2e-3, "mismatch {e_mis} vs {want}");
    assert!((num(&beams[1], "eta_e") - want).abs() < 2e-3);
    assert!(e_lens > 0.0 && e_lens < 1.0);
}

#[test]
fn lp01_matches_its_gaussian_approximation() {
    let (a, nco, ncl, lam) = (4.1e-6, 1.45, 1.444, 1.55e-6);
    let mode = solve_lp_mode(0, 1, a, nco, ncl, lam).unwrap();
    // Marcuse's mode-field fit w/a = 0.65 + 1.619 V^-1.5 + 2.879 V^-6
    let w = a * (0.65 + 1.619 * mode.v.powf(-1.5) + 2.879 * mode.v.powi(-6));
    let medium = Medium::new(ncl * ncl, 1.0, C0 / lam).unwrap();
    let lp = Beam::new(
        BeamShape::Lp { mode: mode.clone(), sin_orientation: false },
        [0.0; 3],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
    )
    .unwrap();
    let g = Beam::new(BeamShape::Gaussian { waist_m: w }, [0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let s = 6.0 * a;
    let samples = beam_plane_samples(&lp, &medium, 2, 0.0, [-s, -s], [s, s], [160, 160]);
    let r = beam_overlap(&samples, &g, &medium, true, 4);
    eprintln!(
        "LP01: V = {:.4} u = {:.6} w = {:.6} b = {:.6} n_eff = {:.8}, Marcuse w = {w:.4e}, eta = {:.6}",
        mode.v, mode.u, mode.w, mode.b, mode.n_eff, r.eta
    );
    assert!(r.eta > 0.99, "LP01 vs Gaussian {}", r.eta);
    // the root satisfies the characteristic equation of the scalar LP modes
    let lhs = mode.u * nyx_analysis::bessel_j(1, mode.u) / nyx_analysis::bessel_j(0, mode.u);
    let rhs = mode.w * nyx_analysis::bessel_k(1, mode.w) / nyx_analysis::bessel_k(0, mode.w);
    assert!((lhs - rhs).abs() < 1e-9 * rhs.abs(), "{lhs} vs {rhs}");
}

#[test]
fn results_are_bitwise_independent_of_thread_count() {
    let freq = C0;
    let k = 2.0 * PI;
    let p = [0.3e-9, 0.0, 1e-9];
    let mesh = shell_mesh(8, 0.6, 1);
    let result = result_from(&mesh, freq, "s", |x| dipole(x, p, k));
    let base = json!({
        "field_prefix": "s", "surface": { "kind": "labels", "labels": [1] }, "n_theta": 19, "n_phi": 12,
        "power_quadrature": { "n_theta": 24, "n_phi": 48 }, "input_power_w": 1.0, "incident_amplitude_v_per_m": 1.0,
        "kspace": { "n": 9 }, "focal": { "axis": "z", "position": 2.0, "u": [-1.0, 1.0, 9], "v": [-1.0, 1.0, 9] }
    });
    let run = |threads: usize| {
        let mut p = base.clone();
        p["threads"] = json!(threads);
        FarFieldAnalysis.run_with_params(&mesh, &result, &p).unwrap()
    };
    let (a, b) = (run(1), run(8));
    assert_eq!(a, b, "farfield differs between 1 and 8 threads");
    let (mesh2, result2, _) = gaussian_slab(2.0, 6.0, 16);
    let ov = |threads: usize| {
        let p = json!({ "threads": threads, "plane": { "n": [40, 40] }, "beams": [ { "kind": "gaussian", "waist_m": 2.5 }, { "kind": "thin_lens", "na": 0.2 } ] });
        BeamOverlapAnalysis.run_with_params(&mesh2, &result2, &p).unwrap()
    };
    assert_eq!(ov(1), ov(8), "overlap differs between 1 and 8 threads");
    eprintln!(
        "determinism: farfield D_max {} and overlap eta {} identical for 1 and 8 threads",
        a["D_max"],
        ov(1)["eta"]
    );
}
