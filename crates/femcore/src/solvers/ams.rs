//! **Auxiliary-space Maxwell solver** (AMS, Hiptmair–Xu) for lowest-order `H(curl)` matrices
//! (port item P1), for the curl-curl-plus-mass operator `A = K + βM` (also complex symmetric).
//!
//! The null space of `K` (the gradients) defeats ordinary multigrid. AMS removes it with two
//! auxiliary nodal spaces, using the discrete gradient `G` and the vertex coordinates:
//!
//! * `G`: edges × nodes; `K G = 0`. The correction `G · AMG(Gᵀ A G)⁻¹ · Gᵀ` handles the gradient
//!   (kernel) part.
//! * `Π = [Π_x Π_y Π_z]`: the Nédélec interpolation of nodal vector fields,
//!   `(Π_c)_{e,j} = ½ |G_{e,j}| · (t_e)_c` with `t_e = Σ_j G_{e,j} x_j` the edge vector. The correction
//!   `Σ_c Π_c · AMG(Π_cᵀ A Π_c)⁻¹ · Π_cᵀ` handles the smooth (non-gradient) part. As in the simplest
//!   hypre variant the vector Laplacian is approximated by its three scalar diagonal blocks.
//!
//! The cycle is the symmetric multiplicative one (`01210` in hypre's notation): forward
//! Gauss–Seidel, `G` correction, `Π` correction (the three components additive), `G` correction,
//! backward Gauss–Seidel. The
//! resulting preconditioner is symmetric (`Mᵀ = M`), so it drives CG (SPD `A`) and COCG (complex
//! symmetric `A`) alike.
//!
//! Limits: lowest-order elements only; `G` must be the *interior* gradient consistent with the
//! boundary conditions of `A` (same constrained edges, constrained nodes removed); no
//! rotated/anisotropic coefficient handling; Gauss–Seidel is sequential.
//!
//! Written from the Hiptmair–Xu paper (2007) and Kolev–Vassilevski (2009), not from hypre's code.

use crate::solvers::amg::{Amg, AmgOptions, gauss_seidel};
use crate::solvers::csr::Csr;
use crate::solvers::direct::SolverError;
use crate::solvers::faer_glue::FaerScalar;
use crate::solvers::linop::Preconditioner;
use std::cell::RefCell;

#[derive(Debug, Clone)]
pub struct AmsOptions {
    /// Options of the AMG hierarchies on the auxiliary nodal operators.
    pub amg: AmgOptions,
    /// Gauss–Seidel sweeps of the edge smoother (each direction).
    pub smoother_sweeps: usize,
}

impl Default for AmsOptions {
    fn default() -> Self {
        Self { amg: AmgOptions::default(), smoother_sweeps: 1 }
    }
}

struct Scratch<T> {
    x: Vec<T>,
    res: Vec<T>,
    tmp_edge: Vec<T>,
    tmp_node: Vec<T>,
    corr_node: Vec<T>,
}

pub struct Ams<T: FaerScalar> {
    a: Csr<T>,
    inv_diag: Vec<T>,
    g: Csr<T>,
    gt: Csr<T>,
    amg_g: Amg<T>,
    pi: [Csr<T>; 3],
    pit: [Csr<T>; 3],
    amg_pi: [Amg<T>; 3],
    sweeps: usize,
    scratch: RefCell<Scratch<T>>,
}

impl<T: FaerScalar> Ams<T> {
    /// `a`: the edge matrix (SPD, or complex symmetric); `g`: discrete gradient (edges × interior
    /// nodes); `edge_vectors[e]`: the vector `x_b − x_a` from the first to the second end node of
    /// edge `e` (the orientation `g`'s row uses: `−1` at `a`, `+1` at `b`).
    ///
    /// The interpolation `Π_c` takes `½ (t_e)_c` at every *interior* end node of `e`; an end node
    /// removed by the boundary condition carries the value zero, so it contributes nothing.
    pub fn new(a: &Csr<T>, g: &Csr<f64>, edge_vectors: &[[f64; 3]], opts: &AmsOptions) -> Result<Self, SolverError> {
        let ne = a.nrows();
        assert_eq!(a.ncols(), ne);
        assert_eq!(g.nrows(), ne);
        assert_eq!(edge_vectors.len(), ne);
        let nn = g.ncols();
        let gt_real = g.transpose();
        let gm: Csr<T> = g.map(T::from_f64);
        let gtm: Csr<T> = gt_real.map(T::from_f64);
        let a_g = gtm.matmul(&a.matmul(&gm));
        let amg_g = Amg::new(&a_g, opts.amg.clone())?;
        let mut pis = Vec::new();
        let mut pits = Vec::new();
        let mut amgs = Vec::new();
        for c in 0..3 {
            let mut trip = Vec::with_capacity(g.nnz());
            for (e, t) in edge_vectors.iter().enumerate() {
                for &j in g.row(e).0 {
                    trip.push((e, j as usize, 0.5 * t[c]));
                }
            }
            let pic: Csr<f64> = Csr::from_triplets(ne, nn, &trip).unwrap();
            let pm: Csr<T> = pic.map(T::from_f64);
            let pmt: Csr<T> = pic.transpose().map(T::from_f64);
            let a_pi = pmt.matmul(&a.matmul(&pm));
            amgs.push(Amg::new(&a_pi, opts.amg.clone())?);
            pis.push(pm);
            pits.push(pmt);
        }
        let inv_diag: Vec<T> =
            a.diagonal().into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect();
        let z = |n: usize| vec![T::ZERO; n];
        let to3 = |mut v: Vec<Csr<T>>| -> [Csr<T>; 3] {
            let c = v.pop().unwrap();
            let b = v.pop().unwrap();
            let a = v.pop().unwrap();
            [a, b, c]
        };
        let to3a = |mut v: Vec<Amg<T>>| -> [Amg<T>; 3] {
            let c = v.pop().unwrap();
            let b = v.pop().unwrap();
            let a = v.pop().unwrap();
            [a, b, c]
        };
        Ok(Self {
            a: a.clone(),
            inv_diag,
            g: gm,
            gt: gtm,
            amg_g,
            pi: to3(pis),
            pit: to3(pits),
            amg_pi: to3a(amgs),
            sweeps: opts.smoother_sweeps,
            scratch: RefCell::new(Scratch { x: z(ne), res: z(ne), tmp_edge: z(ne), tmp_node: z(nn), corr_node: z(nn) }),
        })
    }

    /// `x += Σ_i Q_i · B_i⁻¹ · Q_iᵀ (b − A x)`: the corrections are *additive* (all computed from the
    /// same residual), which for the three `Π_c` is the block-diagonal vector-Laplacian solve and keeps
    /// the whole cycle symmetric.
    fn correct(&self, parts: &[(&Csr<T>, &Csr<T>, &Amg<T>)], b: &[T], sc: &mut Scratch<T>) {
        self.a.spmv(&sc.x, &mut sc.res);
        for (r, bb) in sc.res.iter_mut().zip(b) {
            *r = *bb - *r;
        }
        for &(q, qt, amg) in parts {
            qt.spmv(&sc.res, &mut sc.tmp_node);
            amg.apply(&sc.tmp_node, &mut sc.corr_node);
            q.spmv(&sc.corr_node, &mut sc.tmp_edge);
            for (x, t) in sc.x.iter_mut().zip(&sc.tmp_edge) {
                *x += *t;
            }
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for Ams<T> {
    fn apply(&self, b: &[T], z: &mut [T]) {
        let mut guard = self.scratch.borrow_mut();
        let sc = &mut *guard;
        sc.x.fill(T::ZERO);
        for _ in 0..self.sweeps {
            gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, true);
        }
        let grad = [(&self.g, &self.gt, &self.amg_g)];
        let pis = [0, 1, 2].map(|c| (&self.pi[c], &self.pit[c], &self.amg_pi[c]));
        self.correct(&grad, b, sc);
        self.correct(&pis, b, sc);
        self.correct(&grad, b, sc);
        for _ in 0..self.sweeps {
            gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, false);
        }
        z.copy_from_slice(&sc.x);
    }
}
