//! Elastic modes of a 2D cross-section with an out-of-plane wave number (solid mechanics
//! with an out-of-plane wave number on a 2D component): all three displacement components
//! `u(x, y) exp(-i kz z)` live on the cross-section and the strains contain `-i kz` terms,
//!
//! ```text
//!   e_xx = du_x/dx   e_yy = du_y/dy   e_zz = -i kz u_z
//!   g_xy = du_x/dy + du_y/dx   g_yz = du_z/dy - i kz u_y   g_xz = du_z/dx - i kz u_x
//! ```
//!
//! The stiffness matrix is complex Hermitian, the mass matrix real. With `kz` equal to the phase
//! constant of a travelling wave this gives the dispersion of acoustic modes of waveguides, as
//! needed for Brillouin scattering. Materials use the full 6x6 stiffness (anisotropic and
//! rotated crystals included); no plane assumption is made.

use crate::csparse::{c, CCoo, CCsr, C};
use crate::fem::*;
use crate::material::M6;
use crate::mesh::SimplexMesh;
use crate::problem::*;
use crate::scalar::check_labels;
use crate::sparse::Csr;
use std::time::Instant;

fn ccsr_block_diag(blocks: &[&CCsr]) -> CCsr {
    let n: usize = blocks.iter().map(|b| b.nrows).sum();
    let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
    let mut off = 0;
    for b in blocks {
        for r in 0..b.nrows {
            for k in b.indptr[r]..b.indptr[r + 1] {
                indices.push(b.indices[k] + off);
                data.push(b.data[k]);
            }
            indptr.push(indices.len());
        }
        off += b.ncols;
    }
    CCsr { nrows: n, ncols: off, indptr, indices, data }
}

fn to_complex(a: &Csr) -> CCsr {
    CCsr { nrows: a.nrows, ncols: a.ncols, indptr: a.indptr.clone(), indices: a.indices.clone(), data: a.data.iter().map(|&v| c(v, 0.0)).collect() }
}

/// Stiffness matrix over `3 * n_dofs(order)` unknowns (component blocks): complex Hermitian
/// without absorbing layers. With perfectly matched layers (`stretch`: the complex coordinate
/// stretches of every cell) the in-plane derivatives are divided by the stretch, the cell weight
/// is the determinant of the stretch, and the test strain uses the transpose of the strain operator
/// with the conjugate wave-number terms, so the matrix is complex symmetric in the layers.
// The B-matrix tables below are written `row * nd + component * nloc + a` for every entry, including the ones with a
// zero factor, so that all lines read alike.
#[allow(clippy::erasing_op)]
pub fn assemble_stiffness_kz(mesh: &SimplexMesh, c6: &[M6], order: usize, kz: f64, stretch: Option<&[[C; 3]]>) -> CCsr {
    let re = RefElement::new(2, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let beta = c(0.0, -kz);
    let mut coo = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    // B = B0 + beta B1, six strain rows over the 3 * nloc local dofs (comp * nloc + a)
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    for ci in 0..mesh.n_cells() {
        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let dm = &c6[ci];
        let (sx, sy) = stretch.map_or((c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1]));
        let det = sx * sy;
        for q in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, q, &mut grad);
            b0.iter_mut().for_each(|x| *x = c(0.0, 0.0));
            b1.iter_mut().for_each(|x| *x = 0.0);
            for a in 0..nloc {
                let (gx, gy, nv) = (c(grad[a * 2], 0.0) / sx, c(grad[a * 2 + 1], 0.0) / sy, re.n[q][a]);
                // rows: 0 xx, 1 yy, 2 zz, 3 xy, 4 yz, 5 xz
                b0[0 * nd + 0 * nloc + a] = gx;
                b0[nd + nloc + a] = gy;
                b1[2 * nd + 2 * nloc + a] = nv;
                b0[3 * nd + 0 * nloc + a] = gy;
                b0[3 * nd + nloc + a] = gx;
                b0[4 * nd + 2 * nloc + a] = gy;
                b1[4 * nd + nloc + a] = nv;
                b0[5 * nd + 2 * nloc + a] = gx;
                b1[5 * nd + 0 * nloc + a] = nv;
            }
            // D B
            for s in 0..6 {
                for j in 0..nd {
                    let mut acc = c(0.0, 0.0);
                    for t in 0..6 {
                        acc += (b0[t * nd + j] + beta * b1[t * nd + j]) * dm[s][t];
                    }
                    db[s * nd + j] = acc;
                }
            }
            // K += w B^H (D B)
            for i in 0..nd {
                for j in 0..nd {
                    let mut acc = c(0.0, 0.0);
                    for s in 0..6 {
                        let bh = b0[s * nd + i] + beta.conj() * b1[s * nd + i];
                        acc += bh * db[s * nd + j];
                    }
                    ke[i * nd + j] += acc * re.w[q];
                }
            }
        }
        let vol = c(geo.vol[ci], 0.0) * det;
        ke.iter_mut().for_each(|x| *x *= vol);
        mesh.cell_dofs(ci, order, &mut scratch);
        full.clear();
        for k in 0..3 {
            full.extend(scratch.iter().map(|&d| d + k * n1));
        }
        coo.push_block(&full, &full, &ke);
    }
    coo.to_csr()
}

/// The complex eigenvectors behind a `solve_solid_kz` result (for coupling integrals such as the
/// Brillouin gain): displacement over all `3 * n_dofs(order)` dofs (component blocks), the complex
/// frequencies in Hz, in the order of the solutions.
pub struct KzModes {
    pub kz: f64,
    pub freqs: Vec<C>,
    pub vectors: Vec<Vec<C>>,
    pub c6: Vec<M6>,
    pub rho: Vec<f64>,
    pub converged: bool,
}

/// Eigenfrequencies of the cross-section at the out-of-plane wave number `kz` (rad/m).
pub fn solve_solid_kz(mesh: &SimplexMesh, spec: &ProblemSpec, kz: f64) -> Result<SolveResult, String> {
    solve_solid_kz_full(mesh, spec, kz).map(|r| r.0)
}

/// Like `solve_solid_kz`, and also returns the complex mode shapes.
pub fn solve_solid_kz_full(mesh: &SimplexMesh, spec: &ProblemSpec, kz: f64) -> Result<(SolveResult, KzModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Kz(kz))
}

/// What the third displacement component and the third dimension are: the out-of-plane direction
/// with wave number `kz`, or the azimuth of an axisymmetric body with azimuthal order `m`
/// (mesh coordinates x = r, y = z).
#[derive(Clone, Copy, Debug)]
pub enum CrossSection {
    Kz(f64),
    Axisymmetric(f64),
}

/// Elastic modes of a body of revolution with azimuthal dependence `exp(-i m theta)` on the (r, z) half
/// plane (mesh x = r, y = z): the displacement components are (u_r, u_z, u_theta) in this order.
pub fn solve_solid_axisym(mesh: &SimplexMesh, spec: &ProblemSpec, m: f64) -> Result<(SolveResult, KzModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Axisymmetric(m))
}

pub fn solve_cross_section(mesh: &SimplexMesh, spec: &ProblemSpec, kind: CrossSection) -> Result<(SolveResult, KzModes), String> {
    let kz = match kind {
        CrossSection::Kz(k) => k,
        CrossSection::Axisymmetric(m) => m,
    };
    let t0 = Instant::now();
    if mesh.dim != 2 {
        return Err("an out-of-plane wave number applies to a 2D cross-section".into());
    }
    let Study::Eigenfrequency { n_modes, f_search } = spec.study.clone() else {
        return Err("with an out-of-plane wave number the solid mechanics study must be Eigenfrequency".into());
    };
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let mut c6 = spec.materials.cell_stiffness_skipping(mesh, &spec.inactive)?;
    if let CrossSection::Axisymmetric(_) = kind {
        // a crystal's axes stay fixed while the local frame turns around the ring: the stiffness of the axisymmetric
        // model is the average over rotations about the ring axis (the y direction of the (r, z) plane); exact for isotropic materials
        for m in c6.iter_mut() {
            *m = crate::solid_axisym::average_about_axis(m, 1, 16);
        }
    }
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    // fixed boundaries (zero displacement components)
    let mut fixed_set: std::collections::BTreeSet<usize> = untouched_dofs(mesh, order, 3, &spec.inactive).into_iter().collect();
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, .. } = &bc.kind {
            let bd = boundary_dofs(mesh, &bc.labels, order);
            for k in 0..3 {
                if fixed[k] {
                    for &dof in &bd {
                        fixed_set.insert(dof + k * n1);
                    }
                }
            }
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let stretch = if spec.materials.has_pml() { Some(crate::emw_pml::cell_stretches_model_frame(mesh, &spec.materials.pml)) } else { None };
    let k = match kind {
        CrossSection::Kz(_) => assemble_stiffness_kz(mesh, &c6, order, kz, stretch.as_deref()),
        CrossSection::Axisymmetric(m) => crate::solid_axisym::assemble_stiffness_axisym(mesh, &c6, order, m, stretch.as_deref())?,
    };
    let mut rho = spec.materials.cell_values(mesh, "rho", None)?;
    for (c, r) in rho.iter_mut().enumerate() {
        if spec.inactive.contains(&mesh.cell_labels[c]) {
            *r = 0.0;
        }
    }
    // mass: weighted by the determinant of the stretch in the layers (complex there)
    let mc = match &stretch {
        _ if matches!(kind, CrossSection::Axisymmetric(_)) => crate::solid_axisym::assemble_mass_axisym(mesh, &rho, order, stretch.as_deref()),
        None => to_complex(&Csr::block_diag(&[&assemble_mass(mesh, order, &rho), &assemble_mass(mesh, order, &rho), &assemble_mass(mesh, order, &rho)])),
        Some(st) => {
            let w: Vec<C> = (0..mesh.n_cells()).map(|i| st[i][0] * st[i][1]).collect();
            let re_part: Vec<f64> = rho.iter().zip(&w).map(|(r, z)| r * z.re).collect();
            let im_part: Vec<f64> = rho.iter().zip(&w).map(|(r, z)| r * z.im).collect();
            let (mr, mi) = (assemble_mass(mesh, order, &re_part), assemble_mass(mesh, order, &im_part));
            let m1 = |a: &Csr, b: &Csr| -> CCsr { to_complex(a).add_scaled(&to_complex(b), c(0.0, 1.0)) };
            let one = m1(&mr, &mi);
            let blocks: Vec<&CCsr> = vec![&one, &one, &one];
            ccsr_block_diag(&blocks)
        }
    };
    let free = free_dofs(n, &fixed);
    if n_modes < 1 || n_modes + 1 >= free.len() {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }
    let kf = k.submatrix(&free, &free);
    let mf = mc.submatrix(&free, &free);
    // the shift: the requested frequency, or just below the spectrum for the lowest modes
    let sigma = if f_search > 0.0 { (2.0 * std::f64::consts::PI * f_search).powi(2) } else { -1.0 };
    let lrb_labels: Vec<i32> = spec.bcs.iter().filter(|b| matches!(b.kind, BcKind::LowReflecting)).flat_map(|b| b.labels.iter().copied()).collect();
    // (frequencies, displacement vectors over the free dofs, converged, Krylov dimension)
    type Found = (Vec<(C, usize)>, Vec<Vec<C>>, bool, usize);
    let (mut modes, vectors, converged, kdim): Found = if lrb_labels.is_empty() {
        let pairs = crate::csparse::eigs_complex(&kf, &mf, c(sigma, 0.0), n_modes, 1e-9, 400)?;
        let modes = pairs
            .values
            .iter()
            .enumerate()
            .filter(|(i, _)| pairs.residuals.get(*i).is_none_or(|r| *r < 1e-6))
            .map(|(i, l)| {
                let mut f = l.sqrt() / (2.0 * std::f64::consts::PI);
                if f.re < 0.0 {
                    f = -f;
                }
                (f, i)
            })
            .collect();
        (modes, pairs.vectors.clone(), pairs.converged, pairs.dim)
    } else {
        if f_search <= 0.0 {
            return Err("low-reflecting boundaries need a search frequency (the problem is quadratic in the frequency)".into());
        }
        let cfull = crate::solid_lrb::damping_matrix(mesh, &c6, &rho, &lrb_labels, order, 3);
        let cf = cfull.submatrix(&free, &free);
        let (omegas, vecs, conv, dim) = crate::solid_lrb::quadratic_modes(&kf, &mf, &cf, 2.0 * std::f64::consts::PI * f_search, n_modes)?;
        let modes = omegas.iter().enumerate().map(|(i, w)| (*w / (2.0 * std::f64::consts::PI), i)).collect();
        (modes, vecs, conv, dim)
    };
    let pairs = crate::csparse::CEigPairs { values: vec![], vectors, residuals: vec![], dim: kdim, converged };
    let has_loss = stretch.is_some() || !lrb_labels.is_empty();
    modes.sort_by(|a, b| a.0.re.partial_cmp(&b.0.re).unwrap());
    let mut solutions = vec![];
    let mut scalars = vec![];
    let mut out_freqs: Vec<C> = vec![];
    let mut out_vectors: Vec<Vec<C>> = vec![];
    for (j, &(fz, i)) in modes.iter().enumerate() {
        let freq = fz.re;
        let mut u = vec![C::new(0.0, 0.0); n];
        for (a, &dof) in free.iter().enumerate() {
            u[dof] = pairs.vectors[i][a];
        }
        out_freqs.push(fz);
        out_vectors.push(u.clone());
        // rotate the phase so the largest component is real, scale to unit maximum modulus
        let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);
        let comp = |k: usize| -> Vec<C> { (0..n1).map(|a| u[k * n1 + a] * rot * scale).collect() };
        let (ux, uy, uz) = (comp(0), comp(1), comp(2));
        let mag: Vec<f64> = (0..n1).map(|a| (ux[a].norm_sqr() + uy[a].norm_sqr() + uz[a].norm_sqr()).sqrt()).collect();
        let re = |v: &[C]| -> Vec<f64> { v.iter().map(|z| z.re).collect() };
        let fields = vec![
            Field { name: "|u|".into(), unit: "a.u.".into(), values: mag },
            Field { name: "ux".into(), unit: "a.u.".into(), values: re(&ux) },
            Field { name: "uy".into(), unit: "a.u.".into(), values: re(&uy) },
            Field { name: "uz".into(), unit: "a.u.".into(), values: re(&uz) },
            Field { name: "|uz|".into(), unit: "a.u.".into(), values: uz.iter().map(|z| z.norm()).collect() },
        ];
        let mut disp = vec![0.0; n1 * 2];
        for a in 0..n1 {
            disp[a * 2] = ux[a].re;
            disp[a * 2 + 1] = uy[a].re;
        }
        if has_loss && fz.im.abs() > 0.0 {
            let q = fz.re / (2.0 * fz.im.abs());
            solutions.push(Solution { label: format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, fz.re, fz.im, q), fields, disp: Some(disp), value: Some(freq), ..Default::default() });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", fz.re, fz.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        } else {
            solutions.push(Solution { label: format!("mode {}: {:.6} Hz", j + 1, freq), fields, disp: Some(disp), value: Some(freq), ..Default::default() });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
        }
    }
    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let what = match kind {
        CrossSection::Kz(_) => format!("out-of-plane wave number kz = {kz:.6e} rad/m"),
        CrossSection::Axisymmetric(_) => format!("axisymmetric, azimuthal order m = {kz}"),
    };
    let mut log = vec![format!("{what}; complex Hermitian eigenproblem, Krylov dimension {}{}", pairs.dim, if pairs.converged { "" } else { " (not every pair converged; those are omitted)" })];
    if fixed.is_empty() {
        log.push("no fixed boundary: the lowest modes may be rigid-body-like".into());
    }
    let result = SolveResult { order, dim: 2, solutions, scalars, log, seconds: t0.elapsed().as_secs_f64(), dofs: n, ..Default::default() };
    Ok((result, KzModes { kz, freqs: out_freqs, vectors: out_vectors, c6, rho, converged: pairs.converged }))
}
