//! Linear elastic solid mechanics (linear elastic material) on P1/P2
//! triangles (plane stress or plane strain) and tetrahedra: static solve with boundary loads,
//! prescribed displacements and thermal expansion, eigenfrequencies, stress and von Mises.
//! Materials may be anisotropic (a 6x6 stiffness per domain, see `material`).
//!
//! Everything is in the wafer axes (`section`): a 2D mesh is the x-z cross-section, y the axis out of it. Dof layout: all displacements
//! along the first mesh axis, then along the second (then the third): (x, z) on a cross-section, (x, y, z) in 3D. Voigt order of stresses
//! and strains with engineering shear strains: (xx, yy, zz, xy, yz, xz) in 3D; on a cross-section the Voigt order of the axes (x, z, y),
//! (xx, zz, yy, xz, zy, xy), that `section::plane_stiffness` gives the stiffness in.
//!
//! 2D: plane strain uses the (xx, zz, xz) block of the stiffness with eps_yy = 0; plane stress
//! condenses eps_yy out (sigma_yy = 0). Out-of-plane shear strains are taken as zero. The
//! thermal strain is `alpha (T - Tref)` (a tensor for anisotropic expansion); in plane strain
//! its out-of-plane part loads the in-plane stress through the stiffness, which for an isotropic
//! material is the (1 + nu) alpha dT in-plane free expansion that plane strain gives.

use crate::linalg::driver::EigOpts;
use crate::linalg::driver::eigs_sym;
use crate::fem::*;
use crate::linalg::csr::{c, C};
use crate::fields::{gradients_at_nodes, node_average, CellFields};
use crate::material::{M3, M6};
use crate::mesh::SimplexMesh;
use crate::problem::*;
use crate::dfs::check_labels;
use crate::linalg::csr::RCsr as Csr;
use rayon::prelude::*;
use std::time::Instant;

/// Indices of the in-plane components (xx, zz, xz) in the solver's Voigt order on a cross-section.
const INPLANE: [usize; 3] = [0, 1, 3];

/// Constitutive matrix, row-major `ns x ns` (ns = 3 in 2D, 6 in 3D).
pub(crate) fn d_matrix(c: &M6, dim: usize, plane: Plane) -> Vec<f64> {
    if dim == 3 {
        return c.iter().flat_map(|r| r.iter().copied()).collect();
    }
    let mut d = vec![0.0; 9];
    for (a, &i) in INPLANE.iter().enumerate() {
        for (b, &j) in INPLANE.iter().enumerate() {
            d[a * 3 + b] = match plane {
                Plane::Strain => c[i][j],
                Plane::Stress => c[i][j] - c[i][2] * c[2][j] / c[2][2],
            };
        }
    }
    d
}

/// Thermal strain vector (solver Voigt order, engineering shear) for expansion tensor `a` and
/// temperature rise `dt`.
fn thermal_strain(a: &M3, dt: f64) -> [f64; 6] {
    [a[0][0] * dt, a[1][1] * dt, a[2][2] * dt, 2.0 * a[0][1] * dt, 2.0 * a[1][2] * dt, 2.0 * a[0][2] * dt]
}

/// The stress (in the components that load the model: all six in 3D, xx, yy, xy in 2D) that the
/// thermal strain `eth` would produce if the body were fully constrained.
fn thermal_stress(c: &M6, dim: usize, plane: Plane, eth: &[f64; 6]) -> Vec<f64> {
    if dim == 3 {
        return (0..6).map(|i| (0..6).map(|j| c[i][j] * eth[j]).sum()).collect();
    }
    match plane {
        Plane::Strain => INPLANE.iter().map(|&i| (0..6).map(|j| c[i][j] * eth[j]).sum()).collect(),
        Plane::Stress => {
            let d = d_matrix(c, 2, Plane::Stress);
            (0..3).map(|a| (0..3).map(|b| d[a * 3 + b] * eth[INPLANE[b]]).sum()).collect()
        }
    }
}

/// Strain-displacement matrix `B[s * (d*nloc) + comp*nloc + a]` from shape gradients `grad[a*d + k]`.
pub(crate) fn b_matrix(grad: &[f64], nloc: usize, d: usize, b: &mut [f64]) {
    let ns = if d == 2 { 3 } else { 6 };
    b[..ns * d * nloc].iter_mut().for_each(|x| *x = 0.0);
    let nd = d * nloc;
    for k in 0..d {
        for a in 0..nloc {
            b[k * nd + k * nloc + a] = grad[a * d + k];
        }
    }
    let shear: &[(usize, usize)] = if d == 2 { &[(0, 1)] } else { &[(0, 1), (1, 2), (0, 2)] };
    for (r, &(p, q)) in shear.iter().enumerate() {
        let row = (d + r) * nd;
        for a in 0..nloc {
            b[row + p * nloc + a] = grad[a * d + q];
            b[row + q * nloc + a] = grad[a * d + p];
        }
    }
}

pub(crate) fn dof_layout(mesh: &SimplexMesh, order: usize, c: usize, scratch: &mut Vec<usize>, full: &mut Vec<usize>) {
    let n1 = mesh.n_dofs(order);
    mesh.cell_dofs(c, order, scratch);
    full.clear();
    for k in 0..mesh.dim {
        full.extend(scratch.iter().map(|&d| d + k * n1));
    }
}

/// Stiffness matrix with a (possibly anisotropic) stiffness per cell.
pub fn assemble_solid(mesh: &SimplexMesh, c6: &[M6], order: usize, plane: Plane, thickness: f64) -> Csr {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let n = d * mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let (mut scratch, mut full) = (Vec::new(), Vec::new());
        let mut grad = vec![0.0; nloc * d];
        let mut b = vec![0.0; ns * nd];
        let mut db = vec![0.0; ns * nd];
        let mut ke = vec![0.0; nd * nd];
        for c in range {
            let dm = d_matrix(&c6[c], d, plane);
            ke.iter_mut().for_each(|x| *x = 0.0);
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                b_matrix(&grad, nloc, d, &mut b);
                for s in 0..ns {
                    for j in 0..nd {
                        let mut acc = 0.0;
                        for t in 0..ns {
                            acc += dm[s * ns + t] * b[t * nd + j];
                        }
                        db[s * nd + j] = acc;
                    }
                }
                for i in 0..nd {
                    for j in 0..nd {
                        let mut acc = 0.0;
                        for s in 0..ns {
                            acc += b[s * nd + i] * db[s * nd + j];
                        }
                        ke[i * nd + j] += re.w[q] * acc;
                    }
                }
            }
            let sc = geo.vol[c] * if d == 2 { thickness } else { 1.0 };
            ke.iter_mut().for_each(|x| *x *= sc);
            dof_layout(mesh, order, c, &mut scratch, &mut full);
            coo.push_block(&full, &full, &ke);
        }
    })
}

/// Load vector of the thermal strain `alpha (T - Tref)`; `t` is the temperature at the scalar
/// dofs of the mesh and `alpha` the expansion tensor of every cell.
#[allow(clippy::too_many_arguments)]
pub fn thermal_load(mesh: &SimplexMesh, c6: &[M6], alpha: &[M3], t_field: &[f64], tref: f64, order: usize, plane: Plane, thickness: f64) -> Vec<f64> {
    let d = mesh.dim;
    let re = RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let n1 = mesh.n_dofs(order);
    let nc = mesh.n_cells();
    // cells in blocks on several threads, each with its own load vector, added up at the end
    let nt = if nc < 4000 { 1 } else { crate::par::threads_up_to(8) };
    let chunk = nc.div_ceil(nt).max(1);
    let (re, geo) = (&re, &geo);
    let parts: Vec<Vec<f64>> = (0..nt)
        .into_par_iter()
        .map(|t| {
            let mut f = vec![0.0; d * n1];
            let (mut scratch, mut full) = (Vec::new(), Vec::new());
            let mut grad = vec![0.0; nloc * d];
            let mut b = vec![0.0; ns * nd];
            for c in t * chunk..((t + 1) * chunk).min(nc) {
                mesh.cell_dofs(c, order, &mut scratch);
                let mut fe = vec![0.0; nd];
                for q in 0..re.w.len() {
                    let tq: f64 = (0..nloc).map(|a| re.n[q][a] * t_field[scratch[a]]).sum();
                    let sig = thermal_stress(&c6[c], d, plane, &thermal_strain(&alpha[c], tq - tref));
                    shape_gradients(re, geo, c, q, &mut grad);
                    b_matrix(&grad, nloc, d, &mut b);
                    for i in 0..nd {
                        let mut acc = 0.0;
                        for s in 0..ns {
                            acc += b[s * nd + i] * sig[s];
                        }
                        fe[i] += re.w[q] * acc;
                    }
                }
                let sc = geo.vol[c] * if d == 2 { thickness } else { 1.0 };
                dof_layout(mesh, order, c, &mut scratch, &mut full);
                for (i, &dof) in full.iter().enumerate() {
                    f[dof] += fe[i] * sc;
                }
            }
            f
        })
        .collect();
    let mut parts = parts.into_iter();
    let mut f = parts.next().unwrap_or_else(|| vec![0.0; d * n1]);
    for p in parts {
        f.iter_mut().zip(&p).for_each(|(a, b)| *a += b);
    }
    f
}

/// Cauchy stress (solver Voigt order, six components, also in 2D) at every local node of every
/// cell from the displacement gradients `grad[p*d*d + i*d + k] = d u_i / d x_k`. `eth` holds the
/// thermal strain vector per cell node, or is empty.
pub fn stress_at_nodes(grad: &[f64], nc: usize, nloc: usize, d: usize, c6: &[M6], plane: Plane, eth: &[[f64; 6]]) -> Vec<[f64; 6]> {
    let mut out = vec![[0.0; 6]; nc * nloc];
    for c in 0..nc {
        let dm = if d == 2 && plane == Plane::Stress { Some(d_matrix(&c6[c], 2, Plane::Stress)) } else { None };
        for a in 0..nloc {
            let p = c * nloc + a;
            let g = &grad[p * d * d..(p + 1) * d * d];
            let mut e = [0.0; 6];
            for i in 0..d {
                e[i] = g[i * d + i];
            }
            e[3] = g[1] + g[d]; // gamma_xy = du/dy + dv/dx
            if d == 3 {
                e[4] = g[3 + 2] + g[3 * 2 + 1]; // gamma_yz
                e[5] = g[2] + g[3 * 2]; // gamma_xz
            }
            let th = if eth.is_empty() { [0.0; 6] } else { eth[p] };
            let mut m = [0.0; 6];
            for i in 0..6 {
                m[i] = e[i] - th[i];
            }
            match &dm {
                None => {
                    for i in 0..6 {
                        out[p][i] = (0..6).map(|j| c6[c][i][j] * m[j]).sum();
                    }
                }
                Some(dm) => {
                    // plane stress: only the in-plane components carry stress
                    for (ai, &i) in INPLANE.iter().enumerate() {
                        out[p][i] = (0..3).map(|bj| dm[ai * 3 + bj] * m[INPLANE[bj]]).sum();
                    }
                }
            }
        }
    }
    out
}

/// von Mises stress of a six-component stress.
pub fn von_mises(s: &[f64; 6]) -> f64 {
    let (sx, sy, sz, sxy, syz, sxz) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2)) + 3.0 * (sxy * sxy + syz * syz + sxz * sxz)).sqrt()
}

pub struct Setup {
    pub d: usize,
    pub n1: usize,
    pub c6: Vec<M6>,
    pub fixed: Vec<usize>,
    pub vals: Vec<f64>,
}

pub fn setup(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<Setup, String> {
    let d = mesh.dim;
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let c6: Vec<M6> = spec.materials.cell_stiffness_skipping(mesh, &spec.inactive)?.iter().map(|c| crate::section::plane_stiffness(c, d)).collect();
    let n1 = mesh.n_dofs(order);
    let mut seen = std::collections::BTreeMap::new();
    for dof in untouched_dofs(mesh, order, d, &spec.inactive) {
        seen.insert(dof, 0.0);
    }
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, values } = &bc.kind {
            let bd = boundary_dofs(mesh, &bc.labels, order);
            for (k, &w) in crate::section::wafer_axes(d).iter().enumerate() {
                if fixed[w] {
                    for &dof in &bd {
                        seen.entry(dof + k * n1).or_insert(values[w]);
                    }
                }
            }
        }
    }
    let fixed: Vec<usize> = seen.keys().copied().collect();
    let vals: Vec<f64> = seen.values().copied().collect();
    Ok(Setup { d, n1, c6, fixed, vals })
}

pub struct Thermal {
    pub alpha: Vec<M3>,
    pub t: Vec<f64>,
    pub tref: f64,
}

/// Stage 3 of the real solid kernels: plot fields (`|u|`, `von Mises`, `u{x,y,z}`, stresses under wafer names), the
/// displacement per dof, and the same components per cell without averaging.
fn fields_of(mesh: &SimplexMesh, u: &[f64], order: usize, st: &Setup, plane: Plane, thermal: Option<&Thermal>, with_stress: bool) -> (Vec<Field>, Vec<f64>, CellFields) {
    let d = st.d;
    let n1 = st.n1;
    let nloc = mesh.nloc(order);
    let nc = mesh.n_cells();
    let mut fields = vec![];
    let mut cells = CellFields::new(mesh, order);
    let unorm: Vec<f64> = (0..n1).map(|i| (0..d).map(|k| u[k * n1 + i] * u[k * n1 + i]).sum::<f64>().sqrt()).collect();
    fields.push(Field { name: "|u|".into(), unit: "m".into(), values: unorm });
    if with_stress {
        let mut tl = Instant::now();
        let grad = gradients_at_nodes(mesh, u, d, order);
        lap(&mut tl, "  gradients");
        let mut eth = vec![];
        if let Some(th) = thermal {
            let mut dofs = Vec::new();
            for c in 0..nc {
                mesh.cell_dofs(c, order, &mut dofs);
                for a in 0..nloc {
                    eth.push(thermal_strain(&th.alpha[c], th.t[dofs[a]] - th.tref));
                }
            }
        }
        let s = stress_at_nodes(&grad, nc, nloc, d, &st.c6, plane, &eth);
        lap(&mut tl, "  stress at nodes");
        let vm: Vec<f64> = s.iter().map(von_mises).collect();
        fields.push(Field { name: "von Mises".into(), unit: "Pa".into(), values: node_average(mesh, &vm, 1, order) });
        cells.push_real("von Mises", "Pa", vm);
        for (k, ax) in crate::section::component_names(d).iter().enumerate() {
            fields.push(Field { name: format!("u{ax}"), unit: "m".into(), values: u[k * n1..(k + 1) * n1].to_vec() });
            cells.push_real(format!("u{ax}"), "m", crate::dfs::cell_values_at_nodes(mesh, &u[k * n1..(k + 1) * n1], order));
        }
        let comps: &[usize] = if d == 2 { &[0, 1, 3, 2] } else { &[0, 1, 2, 3, 4, 5] };
        for &idx in comps {
            let v: Vec<f64> = s.iter().map(|x| x[idx]).collect();
            fields.push(Field { name: crate::section::stress_name(d, idx).into(), unit: "Pa".into(), values: node_average(mesh, &v, 1, order) });
            cells.push_real(crate::section::stress_name(d, idx), "Pa", v);
        }
    } else {
        for (k, ax) in crate::section::component_names(d).iter().enumerate() {
            fields.push(Field { name: format!("u{ax}"), unit: "m".into(), values: u[k * n1..(k + 1) * n1].to_vec() });
            cells.push_real(format!("u{ax}"), "m", crate::dfs::cell_values_at_nodes(mesh, &u[k * n1..(k + 1) * n1], order));
        }
    }
    let mut disp = vec![0.0; n1 * d];
    for i in 0..n1 {
        for k in 0..d {
            disp[i * d + k] = u[k * n1 + i];
        }
    }
    (fields, disp, cells)
}

/// Stage 3 of the complex-frequency solid kernels (PML, Floquet, low-reflecting): the mode shape `u` over all dofs
/// (component blocks) with its largest component made real and of unit modulus, as the plot fields `|u|`, `Re u{x,y,z}`,
/// the displacement of the real parts, and the complex components per cell.
pub(crate) fn complex_mode_fields(mesh: &SimplexMesh, order: usize, u: &[C]) -> (Vec<Field>, Vec<f64>, CellFields) {
    let (d, n1) = (mesh.dim, mesh.n_dofs(order));
    let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let scale = 1.0 / big.norm().max(1e-300);
    let mag: Vec<f64> = (0..n1).map(|a| (0..d).map(|k| u[k * n1 + a].norm_sqr()).sum::<f64>().sqrt() * scale).collect();
    let mut fields = vec![Field { name: "|u|".into(), unit: "a.u.".into(), values: mag }];
    let mut disp = vec![0.0; n1 * d];
    let mut cells = CellFields::new(mesh, order);
    let nloc = mesh.nloc(order);
    for (comp, ax) in crate::section::component_names(d).iter().enumerate() {
        let v: Vec<f64> = (0..n1).map(|a| (u[comp * n1 + a] * rot * scale).re).collect();
        for a in 0..n1 {
            disp[a * d + comp] = v[a];
        }
        fields.push(Field { name: format!("Re u{ax}"), unit: "a.u.".into(), values: v });
        let (mut re_c, mut im_c) = (Vec::with_capacity(mesh.n_cells() * nloc), Vec::with_capacity(mesh.n_cells() * nloc));
        let mut dofs = Vec::new();
        for cell in 0..mesh.n_cells() {
            mesh.cell_dofs(cell, order, &mut dofs);
            for &dof in dofs.iter().take(nloc) {
                let z = u[comp * n1 + dof] * rot * scale;
                re_c.push(z.re);
                im_c.push(z.im);
            }
        }
        cells.push_complex(format!("u{ax}"), "a.u.", re_c, im_c);
    }
    (fields, disp, cells)
}

/// Stage 2 of the real solid eigenfrequency kernel: the `nm` eigenpairs of `K u = w^2 M u` on the free dofs nearest the
/// search frequency (shift-invert Lanczos; for the lowest modes of a large problem LOBPCG with a multigrid
/// preconditioner, with a factorisation as the fallback), and a note on which route was taken.
#[allow(clippy::too_many_arguments)]
pub fn solve_real_modes(mesh: &SimplexMesh, spec: &ProblemSpec, st: &Setup, k: &Csr, m: &Csr, free: &[usize], nm: usize, f_search: f64) -> Result<(crate::linalg::driver::EigPairs, Option<String>), String> {
    let (d, order) = (st.d, spec.order);
    let kff = k.submatrix(&free, &free);
    let mff = m.submatrix(&free, &free);
    let mut eig_note: Option<String> = None;
    let sigma = if f_search > 0.0 { (2.0 * std::f64::consts::PI * f_search).powi(2) } else { -1.0 };
    // for the lowest modes (shift below the spectrum) K - sigma M is positive definite and large
    // 3D problems are solved with multigrid-preconditioned CG instead of a sparse factorisation
    let big = if mesh.dim == 3 { free.len() > crate::fem::LOBPCG_MIN_FREE } else { free.len() > crate::fem::AMG_THRESHOLD };
    let want_amg = match spec.solver {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => sigma <= 0.0,
        LinearSolver::Auto => sigma <= 0.0 && big,
    };
    let pairs = if want_amg {
        // no factorisation: block preconditioned CG (LOBPCG) with one multigrid cycle as preconditioner
        let shifted_full = k.add_scaled(&m, -sigma);
        let sys = crate::fem::AmgSystem::new(mesh, order, d, &shifted_full, &st.fixed)?;
        let factor = sys.convergence_factor();
        // multigrid that converges slowly (very thin or stretched parts) loses to a factorisation
        // unless the factorisation would not fit
        let use_amg = spec.solver == LinearSolver::Iterative || factor <= crate::fem::LOBPCG_MAX_FACTOR || free.len() > crate::fem::FACTORISE_MAX_FREE;
        if use_amg {
            // in Auto mode a stalled iteration gives up early and the factorisation takes over
            let guard = spec.solver != LinearSolver::Iterative;
            let res = crate::linalg::driver::lobpcg_guarded(&kff, &mff, nm, nm.clamp(4, 8), &|r: &[f64]| sys.precondition_free(r), 1e-8, 400, guard);
            match res {
                Ok(r) if r.converged || spec.solver == LinearSolver::Iterative => {
                    eig_note = Some(format!("LOBPCG with AMG preconditioner: {} levels {:?}, convergence factor {:.2}, setup {:.1} s{}", sys.n_levels(), sys.sizes(), factor, sys.setup_seconds, if r.converged { "" } else { " (not fully converged)" }));
                    r
                }
                Err(e) if e.starts_with(crate::linalg::driver::STALLED) => {
                    eig_note = Some(format!("{e}; used a sparse factorisation instead"));
                    eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?
                }
                _ => eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?,
            }
        } else {
            eig_note = Some(format!("multigrid converges slowly here (factor {factor:.2}), used a sparse factorisation"));
            eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?
        }
    } else {
        eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?
    };
    Ok((pairs, eig_note))
}

/// Stage 1 of the static solid kernel: the load vector (boundary tractions in wafer axes, thermal expansion) next to the
/// stiffness matrix. A static problem needs a fixed boundary.
pub struct StaticAssembly {
    pub k: Csr,
    pub f: Vec<f64>,
    pub thermal: Option<Thermal>,
}

pub fn assemble_static(mesh: &SimplexMesh, spec: &ProblemSpec, state: &SolveState, st: &Setup, k: Csr) -> Result<StaticAssembly, String> {
    if st.fixed.is_empty() {
        return Err("a static solid problem needs a fixed boundary (otherwise it floats freely)".into());
    }
    let (d, n1, order, thk) = (st.d, st.n1, spec.order, spec.thickness);
    let n = d * n1;
    let mut f = vec![0.0; n];
    for bc in &spec.bcs {
        if let BcKind::Load(t) = &bc.kind {
            for (comp, &w) in crate::section::wafer_axes(d).iter().enumerate() {
                if t[w] != 0.0 {
                    let tc = t[w] * if d == 2 { thk } else { 1.0 };
                    add_boundary_load(mesh, &bc.labels, &|_| tc, order, &mut f, d, comp);
                }
            }
        }
    }
    let mut thermal: Option<Thermal> = None;
    if let Some(th) = &spec.thermal {
        let alpha: Vec<M3> = spec.materials.cell_tensor(mesh, TensorKind::Expansion, Some(0.0))?.iter().map(|a| crate::section::plane_tensor(a, d)).collect();
        let t = match &th.source {
            ThermalSource::Uniform(tv) => vec![*tv; n1],
            ThermalSource::FromHeat => match &state.heat_t {
                Some(t) if t.len() == n1 => t.clone(),
                _ => return Err("no heat-transfer result on this mesh and element order; solve Heat transfer first".into()),
            },
        };
        let fth = thermal_load(mesh, &st.c6, &alpha, &t, th.t_ref, order, spec.plane, thk);
        f.iter_mut().zip(&fth).for_each(|(a, b)| *a += b);
        thermal = Some(Thermal { alpha, t, tref: th.t_ref });
    }
    Ok(StaticAssembly { k, f, thermal })
}

/// Stage 2: the displacement (component blocks) by the linear solver `spec.solver` names, and its note.
pub fn solve_static(mesh: &SimplexMesh, spec: &ProblemSpec, st: &Setup, asm: &StaticAssembly) -> Result<(Vec<f64>, Option<String>), String> {
    solve_spd(mesh, spec.order, st.d, spec.solver, &asm.k, &asm.f, &st.fixed, &st.vals)
}

/// Records the time since `t` under the label `what` in the stage log (and on stderr when NYX_PROFILE is set),
/// then restarts `t`. See `progress::lap` for how labels are treated.
pub(crate) fn lap(t: &mut Instant, what: &str) {
    crate::progress::lap(what, t.elapsed().as_secs_f64());
    *t = Instant::now();
}

pub fn solve_solid(mesh: &SimplexMesh, spec: &ProblemSpec, state: &SolveState) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let mut tl = Instant::now();
    let order = spec.order;
    if spec.materials.has_pml() {
        return match spec.study {
            Study::Eigenfrequency { .. } => crate::sms_pml::solve_solid_pml(mesh, spec),
            _ => Err("perfectly matched layers are used with the Eigenfrequency study (they give complex frequencies and Q)".into()),
        };
    }
    let st = setup(mesh, spec)?;
    lap(&mut tl, "setup (materials, constraints)");
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let thk = spec.thickness;
    let k = assemble_solid(mesh, &st.c6, order, spec.plane, thk);
    lap(&mut tl, "stiffness assembly");
    match &spec.study {
        Study::Static => {
            let asm = assemble_static(mesh, spec, state, &st, k)?;
            lap(&mut tl, "loads");
            let (u, note) = solve_static(mesh, spec, &st, &asm)?;
            lap(&mut tl, "linear solve");
            let (fields, disp, cells) = fields_of(mesh, &u, order, &st, spec.plane, asm.thermal.as_ref(), true);
            lap(&mut tl, "fields (gradients, stress, averaging)");
            let umax = fields[0].values.iter().cloned().fold(0.0, f64::max);
            let vmax = fields.iter().find(|f| f.name == "von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
            Ok(SolveResult {
                order,
                dim: d,
                solutions: vec![Solution { label: "static".into(), fields, cell_fields: Some(cells), disp: Some(disp), ..Default::default() }],
                scalars: vec![("max |u| (m)".into(), Scalar::Num(umax)), ("max von Mises (Pa)".into(), Scalar::Num(vmax))],
                log: note.into_iter().collect(),
                seconds: t0.elapsed().as_secs_f64(),
                dofs: n,
                ..Default::default()
            })
        }
        Study::Eigenfrequency { n_modes, f_search } => {
            let mut rho = spec.materials.cell_values(mesh, "rho", None)?;
            for (c, r) in rho.iter_mut().enumerate() {
                if spec.inactive.contains(&mesh.cell_labels[c]) {
                    *r = 0.0;
                }
            }
            let rho_t: Vec<f64> = rho.iter().map(|r| r * if d == 2 { thk } else { 1.0 }).collect();
            let m1 = assemble_mass(mesh, order, &rho_t);
            let blocks: Vec<&Csr> = (0..d).map(|_| &m1).collect();
            let m = Csr::block_diag(&blocks);
            if !spec.periodic.is_empty() {
                return crate::sms_floquet::solve_floquet(mesh, spec, &st, &k, &m, *n_modes, *f_search, t0);
            }
            let lrb_labels: Vec<i32> = spec.bcs.iter().filter(|b| matches!(b.kind, BcKind::LowReflecting)).flat_map(|b| b.labels.iter().copied()).collect();
            if !lrb_labels.is_empty() {
                return crate::sms_lrb::solve_lrb(mesh, spec, &st, &k, &m, &rho, &lrb_labels, *n_modes, *f_search, t0);
            }
            let free = free_dofs(n, &st.fixed);
            let nm = *n_modes;
            if nm < 1 || nm + 1 >= free.len() {
                return Err("number of modes must be between 1 and the number of free dofs".into());
            }
            let (pairs, eig_note) = solve_real_modes(mesh, spec, &st, &k, &m, &free, nm, *f_search)?;
            let mut order_idx: Vec<usize> = (0..pairs.values.len()).collect();
            order_idx.sort_by(|&a, &b| pairs.values[a].partial_cmp(&pairs.values[b]).unwrap());
            let mut solutions = vec![];
            let mut scalars = vec![];
            for (j, &i) in order_idx.iter().enumerate() {
                let freq = pairs.values[i].max(0.0).sqrt() / (2.0 * std::f64::consts::PI);
                let mut u = vec![0.0; n];
                for (a, &dof) in free.iter().enumerate() {
                    u[dof] = pairs.vectors[i][a];
                }
                let mx = u.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                if mx > 0.0 {
                    u.iter_mut().for_each(|x| *x /= mx);
                }
                let (fields, disp, cells) = fields_of(mesh, &u, order, &st, spec.plane, None, false);
                solutions.push(Solution { label: format!("mode {}: {:.6} Hz", j + 1, freq), fields, cell_fields: Some(cells), disp: Some(disp), value: Some(freq), ..Default::default() });
                scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
            }
            let mut log = vec![format!("Krylov dimension {}{}", pairs.dim, if pairs.converged { "" } else { " (not fully converged)" })];
            log.extend(eig_note);
            if st.fixed.is_empty() {
                log.push("no fixed boundary: the first six (3D) or three (2D) modes are rigid-body motions at ~0 Hz".into());
            }
            Ok(SolveResult { order, dim: d, solutions, scalars, log, seconds: t0.elapsed().as_secs_f64(), dofs: n, ..Default::default() })
        }
        _ => Err("solid mechanics supports Static and Eigenfrequency studies".into()),
    }
}

/// Element-wise stress at the centroid of every cell, one row per cell:
/// `[sxx, sxz, szx, szz, syy, von Mises]` in 2D (the in-plane axes are x and z), `[s (3x3 row-major), von Mises]` in 3D.
pub fn centroid_stress(mesh: &SimplexMesh, u: &[f64], c6: &[M6], plane: Plane, order: usize) -> Vec<Vec<f64>> {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let geo = geometry(mesh);
    let lam: Vec<f64> = vec![1.0 / (d as f64 + 1.0); d + 1];
    let (_, dn) = shape(&lam, order);
    let nloc = mesh.nloc(order);
    let mut out = vec![];
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        let base = c * (d + 1) * d;
        let mut grad = vec![0.0; d * d];
        for i in 0..d {
            for k in 0..d {
                let mut s = 0.0;
                for a in 0..nloc {
                    let mut ga = 0.0;
                    for b in 0..=d {
                        ga += dn[a * (d + 1) + b] * geo.g[base + b * d + k];
                    }
                    s += u[i * n1 + dofs[a]] * ga;
                }
                grad[i * d + k] = s;
            }
        }
        let s = stress_at_nodes(&grad, 1, 1, d, &c6[c..c + 1], plane, &[])[0];
        let vm = von_mises(&s);
        out.push(if d == 2 {
            vec![s[0], s[3], s[3], s[1], s[2], vm]
        } else {
            vec![s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2], vm]
        });
    }
    out
}
