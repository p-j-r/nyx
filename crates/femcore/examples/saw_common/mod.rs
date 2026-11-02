//! Shared helpers of the SAW examples (`saw_s0`, `saw_s1`): silicon (001) film on r-plane sapphire in the wafer frame
//! (x lateral, y propagation, z normal; the cross-section is the x-z plane), bulk velocities of the sapphire cut,
//! materials with the crystal rotations of `chips`, absorbing layers, energy fractions of the complex modes.
#![allow(dead_code)]
use femcore::chips::crystal_rotations;
use femcore::csparse::C;
use femcore::material::*;
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use femcore::solid_kz::{assemble_stiffness_kz, solve_solid_kz_full, KzModes};
use std::time::Instant;

pub const UM: f64 = 1e-6;
pub const T_SI: f64 = 300e-9;
pub const RHO_SI: f64 = 2329.0;
pub const RHO_SAPPHIRE: f64 = 3980.0;

pub fn envf(k: &str, d: f64) -> f64 {
    std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(d)
}

pub fn sapphire_stiffness() -> Stiffness {
    Stiffness::Trigonal { c11: 497.7e9, c12: 162.8e9, c13: 116.0e9, c14: -21.9e9, c33: 500.9e9, c44: 146.8e9 }
}

pub fn silicon_stiffness() -> Stiffness {
    Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }
}

/// Elastic materials of the wafer-frame problem: `sapphire` and `silicon` are lists of domain numbers.
pub fn set_materials(spec: &mut ProblemSpec, psi_deg: f64, domains: &[i32], sapphire: &[i32], silicon: &[i32]) {
    spec.materials.frame = Frame::Wafer;
    let (r_si, r_sap) = crystal_rotations(psi_deg);
    for &d in domains {
        spec.materials.set(d, "E", 100e9);
        spec.materials.set(d, "nu", 0.25);
        spec.materials.set(d, "rho", 1.0);
    }
    for &d in sapphire {
        spec.materials.aniso.insert(d, Anisotropy { stiffness: Some(sapphire_stiffness()), euler_deg: euler_from_rotation(&r_sap), ..Default::default() });
        spec.materials.set(d, "rho", RHO_SAPPHIRE);
    }
    for &d in silicon {
        spec.materials.aniso.insert(d, Anisotropy { stiffness: Some(silicon_stiffness()), euler_deg: euler_from_rotation(&r_si), ..Default::default() });
        spec.materials.set(d, "rho", RHO_SI);
    }
}

pub fn pml(x: PmlDir, z: PmlDir) -> Pml {
    Pml { x, y: PmlDir::Off, z, strength: envf("PML_STRENGTH", 3.0), order: 2.0, re: 1.0 }
}

/// Labels of the boundary facets on the outside of the mesh (used by one cell) whose two end points satisfy `on`.
pub fn outer_labels(mesh: &SimplexMesh, on: &dyn Fn([f64; 2]) -> bool) -> Vec<i32> {
    let ext = mesh.exterior_facets();
    let mut labels = vec![];
    for f in 0..mesh.n_facets() {
        if !ext[f] {
            continue;
        }
        let (a, b) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
        if on([a[0], a[1]]) && on([b[0], b[1]]) {
            labels.push(mesh.facet_labels[f]);
        }
    }
    labels.sort();
    labels.dedup();
    labels
}

/// Bulk waves of the sapphire cut in the sagittal plane (y, z) of the wafer, propagation direction `n = (0, cos t, sin t)`.
pub struct Bulk {
    /// velocities (m/s) along the surface direction (t = 0): shear polarised along x (SH), quasi-shear sagittal, quasi-longitudinal
    pub v_sh: f64,
    pub v_sv: f64,
    pub v_l: f64,
    /// thresholds: the smallest surface phase velocity omega/q at which a bulk wave of that kind exists, min over t of v(t)/cos t
    pub thr_sh: f64,
    pub thr_sv: f64,
    pub thr_l: f64,
    /// largest |Gamma_xy|, |Gamma_xz| relative to the diagonal, and largest stiffness coupling between even and odd x components (GPa)
    pub mirror_gamma: f64,
    pub mirror_c: f64,
}

fn vidx(i: usize, j: usize) -> usize {
    if i == j {
        i
    } else {
        match (i.min(j), i.max(j)) {
            (0, 1) => 3,
            (1, 2) => 4,
            _ => 5,
        }
    }
}

fn gamma(c6: &M6, rho: f64, n: [f64; 3]) -> [[f64; 3]; 3] {
    let mut g = [[0.0; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            for j in 0..3 {
                for l in 0..3 {
                    g[i][k] += c6[vidx(i, j)][vidx(k, l)] * n[j] * n[l] / rho;
                }
            }
        }
    }
    g
}

pub fn sapphire_bulk(psi_deg: f64) -> Bulk {
    let (_, r_sap) = crystal_rotations(psi_deg);
    let c6 = stiffness_global(&sapphire_stiffness(), euler_from_rotation(&r_sap));
    let mut mirror_c = 0.0f64;
    for i in [0usize, 1, 2, 4] {
        for j in [3usize, 5] {
            mirror_c = mirror_c.max(c6[i][j].abs() * 1e-9);
        }
    }
    let (mut v0, mut thr) = ([0.0; 3], [f64::MAX; 3]);
    let mut mirror_gamma = 0.0f64;
    for step in -850..=850 {
        let t = (step as f64 * 0.1).to_radians();
        let g = gamma(&c6, RHO_SAPPHIRE, [0.0, t.cos(), t.sin()]);
        mirror_gamma = mirror_gamma.max(g[0][1].abs().max(g[0][2].abs()) / g[0][0]);
        let sh = g[0][0].sqrt();
        let (a, b, d) = (g[1][1], g[2][2], g[1][2]);
        let (mid, rad) = (0.5 * (a + b), (0.25 * (a - b) * (a - b) + d * d).sqrt());
        let v = [sh, (mid - rad).sqrt(), (mid + rad).sqrt()];
        for k in 0..3 {
            thr[k] = thr[k].min(v[k] / t.cos());
        }
        if step == 0 {
            v0 = v;
        }
    }
    Bulk { v_sh: v0[0], v_sv: v0[1], v_l: v0[2], thr_sh: thr[0], thr_sv: thr[1], thr_l: thr[2], mirror_gamma, mirror_c }
}

/// Strain energy `1/2 Re(u^H K u)` of every mode in the cells where `pred` holds (the stiffness without the layer stretch).
pub fn elastic_energy_where(mesh: &SimplexMesh, modes: &KzModes, order: usize, pred: &dyn Fn(usize) -> bool) -> Vec<f64> {
    let zero = [[0.0; 6]; 6];
    let c6: Vec<M6> = (0..mesh.n_cells()).map(|c| if pred(c) { modes.c6[c] } else { zero }).collect();
    let k = assemble_stiffness_kz(mesh, &c6, order, modes.kz, None);
    modes
        .vectors
        .iter()
        .map(|u| {
            let ku = k.matvec(u);
            0.5 * u.iter().zip(&ku).map(|(a, b)| (a.conj() * *b).re).sum::<f64>()
        })
        .collect()
}

pub fn centroid_y(mesh: &SimplexMesh, c: usize) -> f64 {
    mesh.cell(c).iter().map(|&v| mesh.point(v as usize)[1]).sum::<f64>() / 3.0
}

pub struct Solved {
    pub modes: KzModes,
    pub dofs: usize,
    pub seconds: f64,
}

/// Eigenfrequencies near `f_search` at the wave number `q` along +y (the out-of-plane wave number of the wafer frame; the solver gets -q of its z = -y).
pub fn solve_q(mesh: &SimplexMesh, spec: &ProblemSpec, q: f64, n_modes: usize, f_search: f64) -> Result<Solved, String> {
    let mut s = spec.clone();
    s.study = Study::Eigenfrequency { n_modes, f_search };
    s.out_of_plane = Some(q);
    let inner = femcore::frame_adapter::to_model_frame(&s, 2);
    let t = Instant::now();
    let (r, modes) = solve_solid_kz_full(mesh, &inner, q)?;
    Ok(Solved { modes, dofs: r.dofs, seconds: t.elapsed().as_secs_f64() })
}

pub fn q_factor(f: C) -> f64 {
    if f.im.abs() > 0.0 {
        f.re / (2.0 * f.im.abs())
    } else {
        f64::INFINITY
    }
}

/// Displacement magnitude at the mesh vertices and the three components (wafer names: x lateral, y propagation, z normal;
/// the model frame stores x, z, -y).
pub fn vertex_fields(mesh: &SimplexMesh, modes: &KzModes, order: usize, j: usize) -> Vec<[f64; 4]> {
    let n1 = mesh.n_dofs(order);
    let u = &modes.vectors[j];
    (0..mesh.n_points())
        .map(|i| {
            let (ux, uz, uy) = (u[i].norm(), u[n1 + i].norm(), u[2 * n1 + i].norm());
            [(ux * ux + uy * uy + uz * uz).sqrt(), ux, uy, uz]
        })
        .collect()
}
