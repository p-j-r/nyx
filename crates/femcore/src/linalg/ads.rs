//! **Auxiliary-space divergence solver** (ADS, Hiptmair–Xu 2007; Kolev–Vassilevski 2012 for the practical variant)
//! for lowest-order `H(div)` (Raviart–Thomas, face element) matrices `A = Dᵀ W D + β M_f` (port item P1, the
//! `H(div)` companion of [`crate::linalg::ams`]).
//!
//! The null space of the grad-div term is the set of curls, so, as for AMS, multigrid on the face matrix alone
//! fails and the cure is the *regular decomposition* `RT = Π(H¹)³ + curl(Nédélec) + smooth`: three auxiliary
//! spaces in addition to the smoother.
//!
//! * `C`: faces × edges, the discrete curl (`D C = 0`). The correction `C · AMS(Cᵀ A C + ε M_e)⁻¹ · Cᵀ` removes the
//!   curl (kernel) part; the edge matrix `Cᵀ A C = β CᵀM_fC` is the *curl-curl* matrix, which is singular on gradients,
//!   so a mass term `ε M_e` (edge mass matrix) regularises it, and the edge problem is solved with the existing
//!   [`Ams`] preconditioner (itself a three-space method with the discrete gradient `G`).
//! * `Π = [Π_x Π_y Π_z]`, faces × nodes: the Raviart–Thomas interpolation of nodal vector fields. The correction
//!   `Σ_c Π_c · AMG(Π_cᵀ A Π_c)⁻¹ · Π_cᵀ` handles the smooth part (vector Laplacian approximated by its three scalar
//!   diagonal blocks, as in AMS).
//! * Smoother: symmetric Gauss–Seidel on `A`.
//!
//! The cycle is the symmetric multiplicative one `S, C, Π, C, Sᵀ` (smoother, curl correction, nodal correction, curl
//! correction, transposed smoother), so the preconditioner is symmetric and drives CG on SPD `A` (or COCG on
//! complex symmetric `A`).
//!
//! Limits: lowest-order elements, tetrahedral/any mesh as long as `C`, `G`, `Π` are consistent with each other and
//! with the boundary conditions (all of them restricted to interior dofs, as the generator in [`crate::linalg::maxwell`]
//! does); no rotated coefficients. Written from the two papers, not from hypre's code.

use crate::linalg::amg::{Amg, AmgOptions, gauss_seidel};
use crate::linalg::ams::{Ams, AmsOptions};
use crate::linalg::csr::Csr;
use crate::linalg::direct::SolverError;
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::linop::Preconditioner;
use std::cell::RefCell;

#[derive(Debug, Clone)]
pub struct AdsOptions {
    /// Options of the AMG hierarchies on the nodal operators `Π_cᵀ A Π_c`.
    pub amg: AmgOptions,
    /// Options of the edge-space preconditioner.
    pub ams: AmsOptions,
    /// Gauss–Seidel sweeps of the face smoother (each direction).
    pub smoother_sweeps: usize,
    /// Mass-term weight of the regularised curl-curl matrix, relative to `mean diag(CᵀAC) / mean diag(M_e)`; the
    /// default 0.05 puts the gradient part of the edge problem a little below the top of the curl spectrum.
    pub edge_shift: f64,
}

impl Default for AdsOptions {
    fn default() -> Self {
        Self { amg: AmgOptions::default(), ams: AmsOptions::default(), smoother_sweeps: 1, edge_shift: 0.05 }
    }
}

struct Scratch<T> {
    x: Vec<T>,
    res: Vec<T>,
    tmp_face: Vec<T>,
    tmp_edge: Vec<T>,
    corr_edge: Vec<T>,
    tmp_node: Vec<T>,
    corr_node: Vec<T>,
}

pub struct Ads<T: FaerScalar> {
    a: Csr<T>,
    inv_diag: Vec<T>,
    c: Csr<T>,
    ct: Csr<T>,
    ams: Ams<T>,
    pi: [Csr<T>; 3],
    pit: [Csr<T>; 3],
    amg_pi: [Amg<T>; 3],
    sweeps: usize,
    scratch: RefCell<Scratch<T>>,
}

impl<T: FaerScalar> Ads<T> {
    /// * `a`: the face matrix;
    /// * `curl`: faces × edges (`D C = 0`); `gradient`, `edge_vectors`: the edge-space data of [`Ams::new`]
    ///   (edges × nodes and the edge vectors, same interior edges as the columns of `curl`);
    /// * `edge_mass`: the edge-element mass matrix (regularisation of `Cᵀ A C`);
    /// * `pi[c]`: faces × nodes, the face-flux interpolation of the nodal field `u e_c`.
    pub fn new(
        a: &Csr<T>,
        curl: &Csr<f64>,
        gradient: &Csr<f64>,
        edge_vectors: &[[f64; 3]],
        edge_mass: &Csr<f64>,
        pi: &[Csr<f64>; 3],
        opts: &AdsOptions,
    ) -> Result<Self, SolverError> {
        let nf = a.nrows();
        assert_eq!(a.ncols(), nf);
        assert_eq!(curl.nrows(), nf);
        let ne = curl.ncols();
        assert_eq!(edge_mass.nrows(), ne);
        let cm: Csr<T> = curl.map(T::from_f64);
        let ctm: Csr<T> = curl.transpose().map(T::from_f64);
        let a_edge0 = ctm.matmul(&a.matmul(&cm));
        let em: Csr<T> = edge_mass.map(T::from_f64);
        let mean = |v: Vec<T>| v.iter().map(|x| x.modulus()).sum::<f64>() / v.len().max(1) as f64;
        let ratio = mean(a_edge0.diagonal()) / mean(em.diagonal()).max(f64::MIN_POSITIVE);
        let a_edge = Csr::axpby(T::ONE, &a_edge0, T::from_f64(opts.edge_shift * ratio), &em);
        let ams = Ams::new(&a_edge, gradient, edge_vectors, &opts.ams)?;
        let mut pis = Vec::new();
        let mut pits = Vec::new();
        let mut amgs = Vec::new();
        for pc in pi.iter() {
            let pm: Csr<T> = pc.map(T::from_f64);
            let pmt: Csr<T> = pc.transpose().map(T::from_f64);
            let a_pi = pmt.matmul(&a.matmul(&pm));
            amgs.push(Amg::new(&a_pi, opts.amg.clone())?);
            pis.push(pm);
            pits.push(pmt);
        }
        let inv_diag: Vec<T> =
            a.diagonal().into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect();
        let nn = pi[0].ncols();
        let z = |n: usize| vec![T::ZERO; n];
        fn arr3<X>(mut v: Vec<X>) -> [X; 3] {
            let c = v.pop().unwrap();
            let b = v.pop().unwrap();
            let a = v.pop().unwrap();
            [a, b, c]
        }
        Ok(Self {
            a: a.clone(),
            inv_diag,
            c: cm,
            ct: ctm,
            ams,
            pi: arr3(pis),
            pit: arr3(pits),
            amg_pi: arr3(amgs),
            sweeps: opts.smoother_sweeps,
            scratch: RefCell::new(Scratch {
                x: z(nf),
                res: z(nf),
                tmp_face: z(nf),
                tmp_edge: z(ne),
                corr_edge: z(ne),
                tmp_node: z(nn),
                corr_node: z(nn),
            }),
        })
    }

    fn residual(&self, b: &[T], sc: &mut Scratch<T>) {
        self.a.spmv(&sc.x, &mut sc.res);
        for (r, bb) in sc.res.iter_mut().zip(b) {
            *r = *bb - *r;
        }
    }

    /// `x += C · AMS⁻¹ · Cᵀ (b − A x)`.
    fn curl_correction(&self, b: &[T], sc: &mut Scratch<T>) {
        self.residual(b, sc);
        self.ct.spmv(&sc.res, &mut sc.tmp_edge);
        self.ams.apply(&sc.tmp_edge, &mut sc.corr_edge);
        self.c.spmv(&sc.corr_edge, &mut sc.tmp_face);
        for (x, t) in sc.x.iter_mut().zip(&sc.tmp_face) {
            *x += *t;
        }
    }

    /// `x += Σ_c Π_c AMG_c⁻¹ Π_cᵀ (b − A x)` (the three components additive).
    fn nodal_correction(&self, b: &[T], sc: &mut Scratch<T>) {
        self.residual(b, sc);
        for c in 0..3 {
            self.pit[c].spmv(&sc.res, &mut sc.tmp_node);
            self.amg_pi[c].apply(&sc.tmp_node, &mut sc.corr_node);
            self.pi[c].spmv(&sc.corr_node, &mut sc.tmp_face);
            for (x, t) in sc.x.iter_mut().zip(&sc.tmp_face) {
                *x += *t;
            }
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for Ads<T> {
    fn apply(&self, b: &[T], z: &mut [T]) {
        let mut guard = self.scratch.borrow_mut();
        let sc = &mut *guard;
        sc.x.fill(T::ZERO);
        for _ in 0..self.sweeps {
            gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, true);
        }
        self.curl_correction(b, sc);
        self.nodal_correction(b, sc);
        self.curl_correction(b, sc);
        for _ in 0..self.sweeps {
            gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, false);
        }
        z.copy_from_slice(&sc.x);
    }
}
