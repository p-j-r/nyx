// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Floquet diffraction ports (`diffraction_port`, FR-36, ADR 138) on the planar faces of a transversely Floquet-periodic unit cell.
//!
//! On a port face with outward normal $\hat{\mathbf{n}}$, area $A$ and rectangular lattice $(\Lambda_1, \Lambda_2)$ the tangential field
//! is expanded in orthonormal Rayleigh-Floquet harmonics
//! $$\mathbf{e}_{mn}^{\sigma}(\mathbf{r}_t) = \hat{\mathbf{u}}_{mn}^{\sigma}\,\frac{e^{-i\mathbf{k}_{t,mn}\cdot(\mathbf{r}-\mathbf{r}_0)}}{\sqrt{A}},\qquad
//! \mathbf{k}_{t,mn} = \mathbf{k}_{t,0} + 2\pi\Big(\tfrac{m}{\Lambda_1}, \tfrac{n}{\Lambda_2}\Big),$$
//! with $\hat{\mathbf{u}}^{\text{TM}} = \mathbf{k}_t/|\mathbf{k}_t|$ and $\hat{\mathbf{u}}^{\text{TE}} = \hat{\mathbf{a}}_n\times\hat{\mathbf{u}}^{\text{TM}}$
//! ($\hat{\mathbf{a}}_n$ the positive coordinate axis normal to the face, shared by both ports of a cell). The Bloch vector
//! $\mathbf{k}_{t,0}$ follows from the Floquet phases of the periodic pairs ($\mathbf{k}_{t,0}\cdot\mathbf{s} = \varphi$ for a pair shifted
//! by $\mathbf{s}$). For an outgoing harmonic ($e^{-i\beta\zeta}$ along $\hat{\mathbf{n}}$, $\beta = \sqrt{k^2-|\mathbf{k}_t|^2}$, $\operatorname{Im}\beta\le 0$)
//! the exact Dirichlet-to-Neumann map is $\hat{\mathbf{n}}\times\mu_r^{-1}\nabla\times\mathbf{E} = \gamma\,\mathbf{E}_t$ with
//! $\gamma^{\text{TE}} = i\beta/\mu_r$ and $\gamma^{\text{TM}} = i k^2/(\beta\mu_r)$; the port adds the rank-one terms
//! $\gamma_h(\omega)\,\mathbf{c}_h\bar{\mathbf{c}}_h^{T}$ ($\mathbf{c}_h = \int_\Gamma \mathbf{e}_h\cdot\boldsymbol{\psi}\,dS$) to the curl-curl operator,
//! and an incident harmonic of amplitude $a$ adds $2\gamma_h a\,\mathbf{c}_h$ to the right-hand side. The modal admittance is
//! $Y_h = \gamma_h/(i\omega\mu_0)$; propagating harmonics are the driven-port channels of the S-matrix (unit-power waves,
//! $|S_{pq}|^2$ = diffraction efficiency), evanescent harmonics only enter the DtN operator.

use crate::cavity::{basis3_keyed, Ned3};
use crate::ibc::{sorted3, MU0};
use crate::modes::C0;
use crate::problem::DiffractionPortBc;
use femcore::fem::{build_face_to_cell_3d, geometry};
use femcore::linalg::csr::{c, C};
use nyx_physics::Discretisation;
use std::collections::BTreeSet;
use std::f64::consts::PI;

/// One Rayleigh-Floquet harmonic of a diffraction port.
#[derive(Clone, Debug)]
pub struct FloquetHarmonic {
    pub m: i32,
    pub n: i32,
    /// `true` for TE (`E_t` perpendicular to `k_t`), `false` for TM.
    pub te: bool,
    /// Transverse wavevector in the port's `(axis_1, axis_2)` coordinates (1/m).
    pub kt: [f64; 2],
    /// Unit tangential polarization vector in 3D.
    pub u_hat: [f64; 3],
    /// Drive functional `c_h = int e_h . psi dS` on full Nedelec DOFs (sparse `(dof, value)`, sorted by dof).
    pub c_full: Vec<(usize, C)>,
}

/// One assembled diffraction port.
#[derive(Clone, Debug)]
pub struct AssembledFloquetPort {
    pub bc: DiffractionPortBc,
    /// Outward unit normal.
    pub n_out: [f64; 3],
    /// Coordinate axis normal to the face and the two transverse axes.
    pub axis_n: usize,
    pub axes_t: [usize; 2],
    /// Rectangular lattice constants along `axes_t` (m) and port area (m^2).
    pub lattice: [f64; 2],
    pub area: f64,
    /// Phase origin of the harmonics (minimum corner of the face).
    pub origin: [f64; 3],
    /// Real relative permittivity and permeability of the (homogeneous) port medium.
    pub eps_r: f64,
    pub mu_r: f64,
    /// Bloch transverse wavevector from the periodic phases (1/m).
    pub kt0: [f64; 2],
    pub harmonics: Vec<FloquetHarmonic>,
}

impl AssembledFloquetPort {
    /// Port-medium wavenumber $k(\omega) = \omega\sqrt{\varepsilon_r\mu_r}/c_0$.
    #[inline]
    pub fn k(&self, omega: f64) -> f64 {
        omega * (self.eps_r * self.mu_r).sqrt() / C0
    }

    /// $(\beta, \beta', \beta'')$ of harmonic `h` at `omega` ($\operatorname{Im}\beta \le 0$; guarded away from the Wood anomaly).
    pub fn beta(&self, h: usize, omega: f64) -> (C, C, C) {
        let hm = &self.harmonics[h];
        let n = (self.eps_r * self.mu_r).sqrt();
        let k = omega * n / C0;
        let dk = n / C0;
        let kt2 = hm.kt[0] * hm.kt[0] + hm.kt[1] * hm.kt[1];
        let b2 = k * k - kt2;
        let floor = 1e-9 * k.max(1e-300);
        let beta = if b2 >= 0.0 { c(b2.sqrt().max(floor), 0.0) } else { c(0.0, -(-b2).sqrt().max(floor)) };
        let dbeta = c(k * dk, 0.0) / beta;
        let d2beta = (c(dk * dk, 0.0) - dbeta * dbeta) / beta;
        (beta, dbeta, d2beta)
    }

    /// Whether harmonic `h` propagates at `omega`.
    #[inline]
    pub fn propagating(&self, h: usize, omega: f64) -> bool {
        let hm = &self.harmonics[h];
        let k = self.k(omega);
        k * k - (hm.kt[0] * hm.kt[0] + hm.kt[1] * hm.kt[1]) > 1e-10 * k * k
    }

    /// DtN coefficient $(\gamma_h, \gamma_h', \gamma_h'')$ at `omega`.
    pub fn gamma(&self, h: usize, omega: f64) -> (C, C, C) {
        let (b, b1, b2) = self.beta(h, omega);
        let i = c(0.0, 1.0);
        let inv_mu = 1.0 / self.mu_r;
        if self.harmonics[h].te {
            (i * b * inv_mu, i * b1 * inv_mu, i * b2 * inv_mu)
        } else {
            // gamma = i (n^2 / (c^2 mu)) w^2 / beta
            let s = self.eps_r * self.mu_r / (C0 * C0) * inv_mu;
            let w = c(omega, 0.0);
            let r0 = w * w / b;
            let r1 = w * 2.0 / b - w * w * b1 / (b * b);
            let r2 =
                c(2.0, 0.0) / b - w * b1 * 4.0 / (b * b) + w * w * b1 * b1 * 2.0 / (b * b * b) - w * w * b2 / (b * b);
            (i * r0 * s, i * r1 * s, i * r2 * s)
        }
    }

    /// Modal admittance $Y_h = \gamma_h / (i\omega\mu_0)$ (S).
    #[inline]
    pub fn admittance(&self, h: usize, omega: f64) -> C {
        self.gamma(h, omega).0 / c(0.0, omega * MU0)
    }

    /// Power-wave excitation multiplier $g = 2\gamma_h\sqrt{2/\operatorname{Re}Y_h}$ and its derivatives (zero when evanescent).
    pub fn excitation(&self, h: usize, omega: f64) -> (f64, C, C, C) {
        if !self.propagating(h, omega) {
            return (f64::INFINITY, c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        }
        let (b, b1, b2) = self.beta(h, omega);
        let (b, b1, b2) = (b.re, b1.re, b2.re);
        let (hv, h1, h2) = if self.harmonics[h].te {
            let s = 2.0 * MU0 / self.mu_r;
            (s * omega * b, s * (b + omega * b1), s * (2.0 * b1 + omega * b2))
        } else {
            let s = 2.0 * MU0 * self.eps_r * self.mu_r / (C0 * C0 * self.mu_r);
            let w = omega;
            let q0 = w * w * w / b;
            let q1 = 3.0 * w * w / b - w * w * w * b1 / (b * b);
            let q2 = 6.0 * w / b - 6.0 * w * w * b1 / (b * b) + 2.0 * w * w * w * b1 * b1 / (b * b * b)
                - w * w * w * b2 / (b * b);
            (s * q0, s * q1, s * q2)
        };
        let g = hv.max(1e-300);
        let sq = g.sqrt();
        let dsq = 0.5 * h1 / sq;
        let d2sq = 0.5 * h2 / sq - 0.25 * h1 * h1 / (g * sq);
        let y = self.admittance(h, omega).re.max(1e-300);
        // z used by the S formula: S = a sqrt(Re Y / 2) - delta  <=>  z = 2 / Re Y
        (2.0 / y, c(0.0, 2.0 * sq), c(0.0, 2.0 * dsq), c(0.0, 2.0 * d2sq))
    }

    /// Name of the channel of harmonic `h` (`port:(m,n)TE`).
    pub fn channel_name(&self, h: usize) -> String {
        let hm = &self.harmonics[h];
        format!("{}:({},{}){}", self.bc.name, hm.m, hm.n, if hm.te { "TE" } else { "TM" })
    }
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Degree-4 Dunavant rule on 4 congruent sub-triangles (24 points): barycentric points and weights summing to 1.
fn face_rule() -> Vec<([f64; 3], f64)> {
    let base = femcore::fem::quad_rule(2, 4);
    let subs: [[[f64; 3]; 3]; 4] = [
        [[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]],
        [[0.5, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.5, 0.5]],
        [[0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [0.0, 0.0, 1.0]],
        [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]],
    ];
    let mut out = Vec::with_capacity(4 * base.w.len());
    for s in &subs {
        for q in 0..base.w.len() {
            let l = &base.lam[q];
            let mut b = [0.0; 3];
            for k in 0..3 {
                b[k] = l[0] * s[0][k] + l[1] * s[1][k] + l[2] * s[2][k];
            }
            out.push((b, 0.25 * base.w[q]));
        }
    }
    out
}

/// Assembles every `diffraction_port` of `ports` on a straight-faced simplex mesh.
///
/// `periodic_shifts` are the `(shift vector, phase)` pairs of the Floquet-periodic boundaries (`k_t0 . shift = phase`);
/// `eps`, `mur` are the real cell permittivity tensors and permeabilities; `f_max` sets the automatic order truncation.
pub fn assemble_floquet_ports_3d(
    mesh: &Discretisation,
    ned: &Ned3,
    ports: &[DiffractionPortBc],
    pec_labels: &[i32],
    periodic_shifts: &[([f64; 3], f64)],
    eps: &[[[f64; 3]; 3]],
    eps_im: &[[[f64; 3]; 3]],
    mur: &[f64],
    f_max: f64,
) -> Result<Vec<AssembledFloquetPort>, String> {
    if ports.is_empty() {
        return Ok(Vec::new());
    }
    if periodic_shifts.is_empty() {
        return Err(
            "diffraction_port requires transverse Floquet-periodic boundary pairs ('periodic') on the unit cell".into(),
        );
    }
    let geo = geometry(mesh);
    let face_to_cell = build_face_to_cell_3d(mesh);
    let rule = face_rule();
    let n_dofs = ned.n_dofs();
    let mut out = Vec::with_capacity(ports.len());
    for bc in ports {
        // ---- port triangles
        let mut tris = Vec::new();
        let mut seen = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            let l = mesh.facet_labels[f];
            if !bc.labels.contains(&l) || pec_labels.contains(&l) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            if seen.insert(key) {
                tris.push(key);
            }
        }
        if tris.is_empty() {
            return Err(format!("diffraction_port '{}' selects no boundary triangles", bc.name));
        }
        // ---- plane, normal axis, bbox, medium
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        let mut n_acc = [0.0; 3];
        let mut area = 0.0;
        let mut eps_vals = Vec::with_capacity(tris.len());
        let mut cells = Vec::with_capacity(tris.len());
        for &t in &tris {
            let &(ci, _) = face_to_cell
                .get(&t)
                .ok_or_else(|| format!("diffraction_port '{}': facet {t:?} has no adjacent tetrahedron", bc.name))?;
            if mesh.cell_is_curved(ci) {
                return Err(format!("diffraction_port '{}' does not support curved adjacent cells", bc.name));
            }
            let p: Vec<[f64; 3]> = t
                .iter()
                .map(|&v| {
                    let q = mesh.point(v as usize);
                    [q[0], q[1], q[2]]
                })
                .collect();
            for q in &p {
                for k in 0..3 {
                    lo[k] = lo[k].min(q[k]);
                    hi[k] = hi[k].max(q[k]);
                }
            }
            let cv = mesh.cell(ci);
            let opp = cv.iter().find(|v| !t.contains(v)).copied().unwrap_or(cv[0]);
            let po = mesh.point(opp as usize);
            let mut nr = cross3(sub3(p[1], p[0]), sub3(p[2], p[0]));
            if dot3(nr, sub3(p[0], [po[0], po[1], po[2]])) < 0.0 {
                nr = [-nr[0], -nr[1], -nr[2]];
            }
            let a2 = dot3(nr, nr).sqrt();
            area += 0.5 * a2;
            for k in 0..3 {
                n_acc[k] += 0.5 * nr[k];
            }
            eps_vals.push((eps[ci], eps_im[ci], mur[ci]));
            cells.push(ci);
        }
        let nl = dot3(n_acc, n_acc).sqrt();
        if nl <= 0.0 {
            return Err(format!("diffraction_port '{}' has a zero net normal", bc.name));
        }
        let n_out = [n_acc[0] / nl, n_acc[1] / nl, n_acc[2] / nl];
        let axis_n = (0..3).max_by(|&a, &b| n_out[a].abs().total_cmp(&n_out[b].abs())).unwrap_or(2);
        if n_out[axis_n].abs() < 1.0 - 1e-9
            || (hi[axis_n] - lo[axis_n]).abs() > 1e-9 * (hi[0] - lo[0]).abs().max(hi[1] - lo[1]).max(hi[2] - lo[2])
        {
            return Err(format!("diffraction_port '{}' must be a planar face normal to a coordinate axis", bc.name));
        }
        let axes_t = match axis_n {
            0 => [1, 2],
            1 => [2, 0],
            _ => [0, 1],
        };
        let lattice = [hi[axes_t[0]] - lo[axes_t[0]], hi[axes_t[1]] - lo[axes_t[1]]];
        if !(lattice[0] > 0.0 && lattice[1] > 0.0) {
            return Err(format!("diffraction_port '{}' has a degenerate lattice {lattice:?}", bc.name));
        }
        if (area - lattice[0] * lattice[1]).abs() > 1e-6 * lattice[0] * lattice[1] {
            return Err(format!(
                "diffraction_port '{}' must cover the full rectangular unit-cell face (area {area:.6e} vs lattice {:.6e})",
                bc.name,
                lattice[0] * lattice[1]
            ));
        }
        let (e0, ei0, m0) = eps_vals[0];
        let eps_r = e0[axes_t[0]][axes_t[0]];
        for (e, ei, m) in &eps_vals {
            let iso = (0..3).all(|a| {
                (0..3).all(|b| (e[a][b] - if a == b { eps_r } else { 0.0 }).abs() <= 1e-9 * eps_r.abs().max(1.0))
            });
            if !iso
                || (m - m0).abs() > 1e-9
                || ei.iter().flatten().any(|v| v.abs() > 1e-12)
                || ei0.iter().flatten().any(|v| v.abs() > 1e-12)
            {
                return Err(format!(
                    "diffraction_port '{}' requires a homogeneous, isotropic, lossless medium adjacent to the port face",
                    bc.name
                ));
            }
        }
        if !(eps_r > 0.0 && m0 > 0.0) {
            return Err(format!("diffraction_port '{}' needs a positive port-medium eps_r and mu_r", bc.name));
        }
        // ---- Bloch vector from the periodic phases (least squares on k_t . s_t = phase)
        let (mut a11, mut a12, mut a22, mut r1, mut r2) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (s, ph) in periodic_shifts {
            let (s1, s2) = (s[axes_t[0]], s[axes_t[1]]);
            if s[axis_n].abs() > 1e-9 * (s1.abs() + s2.abs()) {
                continue;
            }
            a11 += s1 * s1;
            a12 += s1 * s2;
            a22 += s2 * s2;
            r1 += s1 * ph;
            r2 += s2 * ph;
        }
        let det = a11 * a22 - a12 * a12;
        let kt0 = if det.abs() > 1e-30 * (a11 * a22).max(1e-300) {
            [(a22 * r1 - a12 * r2) / det, (a11 * r2 - a12 * r1) / det]
        } else if a11 > 0.0 {
            [r1 / a11, 0.0]
        } else if a22 > 0.0 {
            [0.0, r2 / a22]
        } else {
            return Err(format!("diffraction_port '{}': no periodic pair is transverse to the port face", bc.name));
        };
        // ---- harmonic set
        let k_max = 2.0 * PI * f_max * (eps_r * m0).sqrt() / C0;
        let auto = |lat: f64, kt: f64| -> i32 { (((k_max + kt.abs()) * lat / (2.0 * PI)).floor() as i32 + 1).max(1) };
        let m_max = bc.m_max.unwrap_or_else(|| auto(lattice[0], kt0[0])).max(0);
        let n_max = bc.n_max.unwrap_or_else(|| auto(lattice[1], kt0[1])).max(0);
        if (2 * m_max + 1) * (2 * n_max + 1) > 2500 {
            return Err(format!("diffraction_port '{}': too many Floquet harmonics ({m_max}, {n_max})", bc.name));
        }
        let a_n = {
            let mut v = [0.0; 3];
            v[axis_n] = 1.0;
            v
        };
        let e1 = {
            let mut v = [0.0; 3];
            v[axes_t[0]] = 1.0;
            v
        };
        let mut origin = lo;
        origin[axis_n] = lo[axis_n];
        let inv_sqrt_a = 1.0 / area.sqrt();
        // Tangential traces of the port triangles at the face rule: (dofs, weights, positions, psi)
        let mut traces = Vec::with_capacity(tris.len());
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        for (&t, &ci) in tris.iter().zip(&cells) {
            let v = mesh.cell(ci);
            let mut lf = [0usize; 3];
            for (k, &tv) in t.iter().enumerate() {
                lf[k] = v.iter().position(|&cv| cv == tv).ok_or("diffraction_port: triangle vertex not in its cell")?;
            }
            let bnd = ned.boundary_dofs(mesh, t);
            let cd = &ned.cell_dofs[ci];
            let slots: Vec<usize> = bnd
                .iter()
                .map(|gd| {
                    cd.iter().position(|d| d == gd).ok_or("diffraction_port: boundary dof not in its cell".to_string())
                })
                .collect::<Result<_, _>>()?;
            let gid = ned.key_pairs(mesh, ci);
            let mut gs = [[0.0; 3]; 4];
            for b in 0..4 {
                for k in 0..3 {
                    gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                }
            }
            let p: Vec<[f64; 3]> = t
                .iter()
                .map(|&vv| {
                    let q = mesh.point(vv as usize);
                    [q[0], q[1], q[2]]
                })
                .collect();
            let tri_area = 0.5
                * dot3(cross3(sub3(p[1], p[0]), sub3(p[2], p[0])), cross3(sub3(p[1], p[0]), sub3(p[2], p[0]))).sqrt();
            let mut pts = Vec::with_capacity(rule.len());
            for (b3, wq) in &rule {
                let mut lam4 = [0.0; 4];
                for k in 0..3 {
                    lam4[lf[k]] = b3[k];
                }
                basis3_keyed(ned.order, &lam4, &gs, &gid, &mut w, &mut cu);
                let x = [
                    b3[0] * p[0][0] + b3[1] * p[1][0] + b3[2] * p[2][0],
                    b3[0] * p[0][1] + b3[1] * p[1][1] + b3[2] * p[2][1],
                    b3[0] * p[0][2] + b3[1] * p[1][2] + b3[2] * p[2][2],
                ];
                let psi: Vec<[f64; 3]> = slots
                    .iter()
                    .map(|&s| {
                        let ws = w[s];
                        let wn = ws[axis_n];
                        let mut r = ws;
                        r[axis_n] -= wn;
                        r
                    })
                    .collect();
                pts.push((wq * tri_area, x, psi));
            }
            traces.push((bnd, pts));
        }
        let mut harmonics = Vec::new();
        for m in -m_max..=m_max {
            for n in -n_max..=n_max {
                let kt = [kt0[0] + 2.0 * PI * m as f64 / lattice[0], kt0[1] + 2.0 * PI * n as f64 / lattice[1]];
                let ktn = (kt[0] * kt[0] + kt[1] * kt[1]).sqrt();
                let tm_hat = if ktn > 1e-12 * k_max.max(1e-300) {
                    let mut v = [0.0; 3];
                    v[axes_t[0]] = kt[0] / ktn;
                    v[axes_t[1]] = kt[1] / ktn;
                    v
                } else {
                    e1
                };
                let te_hat = cross3(a_n, tm_hat);
                let mut kt3 = [0.0; 3];
                kt3[axes_t[0]] = kt[0];
                kt3[axes_t[1]] = kt[1];
                for (te, u_hat) in [(true, te_hat), (false, tm_hat)] {
                    let mut acc: std::collections::BTreeMap<usize, C> = std::collections::BTreeMap::new();
                    for (bnd, pts) in &traces {
                        for (wq, x, psi) in pts {
                            let ph = -dot3(kt3, sub3(*x, origin));
                            let e = c(ph.cos(), ph.sin()) * (inv_sqrt_a * wq);
                            for (s, &d) in bnd.iter().enumerate() {
                                let val = dot3(psi[s], u_hat);
                                if val != 0.0 {
                                    *acc.entry(d).or_insert(c(0.0, 0.0)) += e * val;
                                }
                            }
                        }
                    }
                    let c_full: Vec<(usize, C)> = acc.into_iter().filter(|(d, _)| *d < n_dofs).collect();
                    harmonics.push(FloquetHarmonic { m, n, te, kt, u_hat, c_full });
                }
            }
        }
        out.push(AssembledFloquetPort {
            bc: bc.clone(),
            n_out,
            axis_n,
            axes_t,
            lattice,
            area,
            origin,
            eps_r,
            mu_r: m0,
            kt0,
            harmonics,
        });
    }
    Ok(out)
}
