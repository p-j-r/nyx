//! Overlapping **additive and restricted additive Schwarz** preconditioners (port item P4; the role of
//! PETSc `PCASM`) with exact local solves by sparse factorisation and an optional Nicolaides
//! coarse space (one constant per subdomain), for real and complex matrices.
//!
//! For subdomains `Ω_i` (index sets `I_i`, grown by `overlap` layers of the matrix graph) and
//! local matrices `A_i = A[I_i, I_i]`:
//!
//! * [`SchwarzKind::Basic`]: `M⁻¹ r = Σ_i R_iᵀ A_i⁻¹ R_i r` — symmetric for symmetric `A`, so usable in
//!   CG / COCG.
//! * [`SchwarzKind::Restricted`] (RAS, PETSc's default): the prolongation drops the overlap-only dofs,
//!   `M⁻¹ r = Σ_i R̃_iᵀ A_i⁻¹ R_i r` — non-symmetric but usually fewer iterations; use GMRES/FGMRES.
//! * Two-level: `+ Z A₀⁻¹ Zᵀ r` with `A₀ = ZᵀAZ` and `Z` one column per subdomain holding an algebraic
//!   *partition of unity* `w_i(g) = d_i(g)/Σ_j d_j(g)` (`d_i` = graph distance of dof `g` from the boundary of the
//!   overlapping subdomain `i`), i.e. smoothed Nicolaides constants. (Plain indicators of the
//!   non-overlapping parts, tried first, are discontinuous and gave no benefit: 28/41/54 iterations for
//!   4/16/64 subdomains.) The coarse correction needs overlap ≥ 1 to be useful.
//!
//! The subdomain solves are independent. `workers > 1` runs them on scoped threads; the local
//! corrections are added in subdomain order, so the result does not depend on the number of
//! workers (bit for bit). Partitions are given as a dof → part map: [`contiguous_partition`] (what PETSc
//! uses for a sequential `PCASM` with `-pc_asm_blocks`) or [`graph_partition`] (breadth-first
//! region growing on the matrix graph, better shaped subdomains for unstructured numbering).
//! This is groundwork for MPI (item D1): each subdomain is what one rank would own.

use crate::solvers::csr::Csr;
use crate::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use crate::solvers::faer_glue::FaerScalar;
use crate::solvers::linop::Preconditioner;
use crate::solvers::scalar::Scalar;
use std::cell::RefCell;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SchwarzKind {
    Basic,
    Restricted,
}

#[derive(Debug, Clone)]
pub struct SchwarzOptions {
    pub overlap: usize,
    pub kind: SchwarzKind,
    /// Add the Nicolaides coarse correction.
    pub coarse: bool,
    /// Threads for the subdomain solves (1 = sequential).
    pub workers: usize,
    /// Local matrices are SPD (Cholesky first) or general.
    pub definiteness: Definiteness,
}

impl Default for SchwarzOptions {
    fn default() -> Self {
        Self {
            overlap: 1,
            kind: SchwarzKind::Restricted,
            coarse: false,
            workers: 1,
            definiteness: Definiteness::General,
        }
    }
}

/// `part[i] = i·parts/n`: contiguous index ranges of (almost) equal size.
pub fn contiguous_partition(n: usize, parts: usize) -> Vec<usize> {
    (0..n).map(|i| i * parts / n).collect()
}

/// Breadth-first region growing on the (symmetrised) graph of `a`: seeds are spread by repeated
/// farthest-node selection and each region grows by one layer in turn until all nodes are taken.
pub fn graph_partition<T: Scalar>(a: &Csr<T>, parts: usize) -> Vec<usize> {
    let n = a.nrows();
    let adj = symmetric_adjacency(a);
    const NONE: usize = usize::MAX;
    let mut part = vec![NONE; n];
    let mut frontier: Vec<Vec<usize>> = Vec::new();
    // seeds: node 0, then repeatedly the node farthest (graph distance) from all seeds
    let mut dist = vec![usize::MAX; n];
    let mut seeds = vec![0usize];
    let bfs = |seed: usize, dist: &mut Vec<usize>| {
        let mut q = std::collections::VecDeque::new();
        let mut d2 = vec![usize::MAX; n];
        d2[seed] = 0;
        q.push_back(seed);
        while let Some(u) = q.pop_front() {
            for &v in &adj[u] {
                if d2[v as usize] == usize::MAX {
                    d2[v as usize] = d2[u] + 1;
                    q.push_back(v as usize);
                }
            }
        }
        for (a, b) in dist.iter_mut().zip(&d2) {
            *a = (*a).min(*b);
        }
    };
    bfs(0, &mut dist);
    while seeds.len() < parts.min(n) {
        let far = (0..n).filter(|&i| dist[i] != usize::MAX).max_by_key(|&i| dist[i]).unwrap_or(seeds.len());
        seeds.push(far);
        bfs(far, &mut dist);
    }
    for (p, &s) in seeds.iter().enumerate() {
        part[s] = p;
        frontier.push(vec![s]);
    }
    let mut remaining = n - seeds.len();
    while remaining > 0 {
        let mut progressed = false;
        for (p, fr) in frontier.iter_mut().enumerate() {
            let mut next = Vec::new();
            for &u in fr.iter() {
                for &v in &adj[u] {
                    if part[v as usize] == NONE {
                        part[v as usize] = p;
                        next.push(v as usize);
                        remaining -= 1;
                        progressed = true;
                    }
                }
            }
            *fr = next;
        }
        if !progressed {
            // disconnected leftovers: give them to part 0
            for p in part.iter_mut().filter(|p| **p == NONE) {
                *p = 0;
                remaining -= 1;
            }
        }
    }
    part
}

fn symmetric_adjacency<T: Scalar>(a: &Csr<T>) -> Vec<Vec<u32>> {
    let n = a.nrows();
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for i in 0..n {
        for &j in a.row(i).0 {
            if j as usize != i {
                adj[i].push(j);
                adj[j as usize].push(i as u32);
            }
        }
    }
    for l in &mut adj {
        l.sort_unstable();
        l.dedup();
    }
    adj
}

struct Domain<T: FaerScalar> {
    /// Global dofs of the subdomain (sorted).
    idx: Vec<u32>,
    /// Graph distance (≥ 1) of each dof from the subdomain boundary (dofs adjacent to a dof outside it have 1).
    depth: Vec<u32>,
    /// `true` for dofs owned by this subdomain (its non-overlapping part).
    own: Vec<bool>,
    fac: ShiftedFactorization<T>,
}

struct Coarse<T: FaerScalar> {
    /// Coarse basis `Z` (`n × nparts`: smooth partition-of-unity weights) and its transpose.
    z: Csr<T>,
    zt: Csr<T>,
    fac: ShiftedFactorization<T>,
}

pub struct Schwarz<T: FaerScalar> {
    n: usize,
    domains: Vec<Domain<T>>,
    coarse: Option<Coarse<T>>,
    kind: SchwarzKind,
    workers: usize,
    /// Scratch for the per-domain corrections (reused across applications).
    scratch: RefCell<Vec<Vec<T>>>,
}

impl<T: FaerScalar> Schwarz<T> {
    /// `part[i]` assigns dof `i` to a subdomain (`0..nparts`).
    pub fn new(a: &Csr<T>, part: &[usize], opts: &SchwarzOptions) -> Result<Self, SolverError> {
        let n = a.nrows();
        if a.ncols() != n || part.len() != n {
            return Err(SolverError::InvalidArgument("matrix must be square and part must have one entry per dof"));
        }
        let nparts = part.iter().copied().max().map_or(0, |m| m + 1);
        let adj = symmetric_adjacency(a);
        let mut domains = Vec::with_capacity(nparts);
        let mut local_of = vec![u32::MAX; n];
        for p in 0..nparts {
            let own_set: Vec<u32> = (0..n as u32).filter(|&i| part[i as usize] == p).collect();
            if own_set.is_empty() {
                continue;
            }
            // grow by `overlap` layers
            let mut in_dom = vec![false; n];
            let mut layer: Vec<u32> = own_set.clone();
            for &i in &own_set {
                in_dom[i as usize] = true;
            }
            let mut all = own_set.clone();
            for _ in 0..opts.overlap {
                let mut next = Vec::new();
                for &u in &layer {
                    for &v in &adj[u as usize] {
                        if !in_dom[v as usize] {
                            in_dom[v as usize] = true;
                            next.push(v);
                            all.push(v);
                        }
                    }
                }
                layer = next;
            }
            all.sort_unstable();
            for (l, &g) in all.iter().enumerate() {
                local_of[g as usize] = l as u32;
            }
            // local matrix A[I, I]
            let mut trip = Vec::new();
            for (li, &g) in all.iter().enumerate() {
                let (cs, vs) = a.row(g as usize);
                for (&c, &v) in cs.iter().zip(vs) {
                    let lc = local_of[c as usize];
                    if lc != u32::MAX {
                        trip.push((li, lc as usize, v));
                    }
                }
            }
            let local = Csr::from_triplets(all.len(), all.len(), &trip)
                .map_err(|_| SolverError::InvalidArgument("local matrix assembly failed"))?;
            let mut fac = ShiftedFactorization::new(&local, None, opts.definiteness)?;
            fac.set_shift(T::ZERO)?;
            let own: Vec<bool> = all.iter().map(|&g| part[g as usize] == p).collect();
            // depth from the subdomain boundary by multi-source BFS inside the subdomain
            let mut depth = vec![0u32; all.len()];
            let mut queue = std::collections::VecDeque::new();
            for (li, &g) in all.iter().enumerate() {
                if adj[g as usize].iter().any(|&v| local_of[v as usize] == u32::MAX) {
                    depth[li] = 1;
                    queue.push_back(li);
                }
            }
            while let Some(li) = queue.pop_front() {
                for &v in &adj[all[li] as usize] {
                    let lv = local_of[v as usize];
                    if lv != u32::MAX && depth[lv as usize] == 0 {
                        depth[lv as usize] = depth[li] + 1;
                        queue.push_back(lv as usize);
                    }
                }
            }
            // a subdomain without any boundary (the whole domain): constant weight
            for d in depth.iter_mut().filter(|d| **d == 0) {
                *d = 1_000_000;
            }
            for &g in &all {
                local_of[g as usize] = u32::MAX;
            }
            domains.push(Domain { idx: all, depth, own, fac });
        }
        let coarse = if opts.coarse {
            // Z[g, i] = d_i(g) / Σ_j d_j(g) over the subdomains containing g
            let mut total = vec![0.0f64; n];
            for d in &domains {
                for (&g, &w) in d.idx.iter().zip(&d.depth) {
                    total[g as usize] += w as f64;
                }
            }
            let mut trip = Vec::new();
            for (k, d) in domains.iter().enumerate() {
                for (&g, &w) in d.idx.iter().zip(&d.depth) {
                    trip.push((g as usize, k, T::from_f64(w as f64 / total[g as usize])));
                }
            }
            let z = Csr::from_triplets(n, domains.len(), &trip)
                .map_err(|_| SolverError::InvalidArgument("coarse basis assembly failed"))?;
            let zt = z.transpose();
            let a0 = zt.matmul(&a.matmul(&z));
            let mut fac = ShiftedFactorization::new(&a0, None, Definiteness::General)?;
            fac.set_shift(T::ZERO)?;
            Some(Coarse { z, zt, fac })
        } else {
            None
        };
        let scratch = RefCell::new(domains.iter().map(|d| vec![T::ZERO; d.idx.len()]).collect());
        Ok(Self { n, domains, coarse, kind: opts.kind, workers: opts.workers.max(1), scratch })
    }

    pub fn subdomains(&self) -> usize {
        self.domains.len()
    }

    /// Sizes of the (overlapping) subdomains.
    pub fn sizes(&self) -> Vec<usize> {
        self.domains.iter().map(|d| d.idx.len()).collect()
    }

    fn solve_domain(d: &Domain<T>, r: &[T], out: &mut [T]) {
        for (o, &g) in out.iter_mut().zip(&d.idx) {
            *o = r[g as usize];
        }
        d.fac.solve_block_in_place(out, 1);
    }
}

impl<T: FaerScalar + Send + Sync> Schwarz<T> {
    fn local_solves(&self, r: &[T], bufs: &mut [Vec<T>]) {
        let w = self.workers.min(self.domains.len()).max(1);
        if w == 1 {
            for (d, b) in self.domains.iter().zip(bufs.iter_mut()) {
                Self::solve_domain(d, r, b);
            }
        } else {
            let chunk = self.domains.len().div_ceil(w);
            std::thread::scope(|sc| {
                for (ds, bs) in self.domains.chunks(chunk).zip(bufs.chunks_mut(chunk)) {
                    sc.spawn(move || {
                        for (d, b) in ds.iter().zip(bs.iter_mut()) {
                            Self::solve_domain(d, r, b);
                        }
                    });
                }
            });
        }
    }
}

impl<T: FaerScalar + Send + Sync> Preconditioner<T> for Schwarz<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        assert_eq!(r.len(), self.n);
        z.fill(T::ZERO);
        let mut bufs = self.scratch.borrow_mut();
        self.local_solves(r, &mut bufs);
        // add the corrections in subdomain order (independent of the number of workers)
        for (d, b) in self.domains.iter().zip(bufs.iter()) {
            for ((&g, &v), &own) in d.idx.iter().zip(b).zip(&d.own) {
                if self.kind == SchwarzKind::Basic || own {
                    z[g as usize] += v;
                }
            }
        }
        if let Some(c) = &self.coarse {
            let mut rc = vec![T::ZERO; c.zt.nrows()];
            c.zt.spmv(r, &mut rc);
            c.fac.solve_block_in_place(&mut rc, 1);
            let mut corr = vec![T::ZERO; self.n];
            c.z.spmv(&rc, &mut corr);
            for (zi, ci) in z.iter_mut().zip(&corr) {
                *zi += *ci;
            }
        }
    }
}
