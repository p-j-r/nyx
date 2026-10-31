//! Incomplete factorisations as preconditioners (port item P3): **ILU(k)**, **ILUT** and
//! **incomplete Cholesky** (IC(k); Hermitian `LLᴴ` or complex-symmetric `LLᵀ`), for real and complex
//! matrices. Sequential and deterministic; natural ordering (apply a permutation to the matrix
//! first if you want RCM/AMD).
//!
//! * [`Ilu::ilu_k`]: fill pattern by *levels of fill* (`level(i,m) = level(i,j) + level(j,m) + 1 ≤ k`,
//!   `k = 0` keeps the pattern of `A`), numeric factorisation in IKJ order. Same algorithm as PETSc's
//!   `PCILU` with `-pc_factor_levels k` in natural ordering, so iteration counts can be compared.
//! * [`Ilu::ilut`]: Saad's dual-threshold ILUT: drop below `droptol·‖a_i‖₂` and keep at most `lfil` entries in
//!   each of the `L` and `U` parts of every row.
//! * [`Ic`]: incomplete Cholesky on the lower triangle with the same level-of-fill pattern; if a pivot is not
//!   safely positive (Hermitian) or zero (complex symmetric) the factorisation is restarted on
//!   `A + α·diag(A)` with `α` doubled from `1e-3` (Manteuffel shift); the shift used is reported.
//!
//! A tiny `U` pivot in ILU is replaced by `±1e-10·‖a_i‖` (reported in `perturbed_pivots`) instead of failing.

use crate::solvers::csr::Csr;
use crate::solvers::linop::Preconditioner;
use crate::solvers::scalar::Scalar;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

const INF: u32 = u32::MAX;

/// Row-wise sparse triangular storage used by both factorisations: strictly-lower part, diagonal, strictly-upper part.
#[derive(Debug, Clone)]
pub struct Ilu<T: Scalar> {
    n: usize,
    l_ptr: Vec<usize>,
    l_idx: Vec<u32>,
    l_val: Vec<T>,
    diag: Vec<T>,
    u_ptr: Vec<usize>,
    u_idx: Vec<u32>,
    u_val: Vec<T>,
    pub perturbed_pivots: usize,
}

/// Fill pattern: for each row the sorted columns (`< i` lower, `> i` upper) with their levels.
struct Pattern {
    l: Vec<Vec<(u32, u32)>>,
    u: Vec<Vec<(u32, u32)>>,
}

fn symbolic_levels<T: Scalar>(a: &Csr<T>, k: u32) -> Pattern {
    let n = a.nrows();
    let mut l: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    let mut u: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    let mut lev = vec![INF; n];
    for i in 0..n {
        let mut inrow: Vec<u32> = Vec::new();
        let mut heap: BinaryHeap<Reverse<u32>> = BinaryHeap::new();
        for &c in a.row(i).0 {
            if lev[c as usize] == INF {
                lev[c as usize] = 0;
                inrow.push(c);
                if (c as usize) < i {
                    heap.push(Reverse(c));
                }
            }
        }
        if lev[i] == INF {
            lev[i] = 0; // the diagonal is always part of the pattern
            inrow.push(i as u32);
        }
        while let Some(Reverse(j)) = heap.pop() {
            let lij = lev[j as usize];
            for &(m, ljm) in &u[j as usize] {
                let nl = lij.saturating_add(ljm).saturating_add(1);
                if nl <= k {
                    if lev[m as usize] == INF {
                        lev[m as usize] = nl;
                        inrow.push(m);
                        if (m as usize) < i {
                            heap.push(Reverse(m));
                        }
                    } else if nl < lev[m as usize] {
                        lev[m as usize] = nl;
                    }
                }
            }
        }
        inrow.sort_unstable();
        for &c in &inrow {
            let entry = (c, lev[c as usize]);
            match (c as usize).cmp(&i) {
                std::cmp::Ordering::Less => l[i].push(entry),
                std::cmp::Ordering::Greater => u[i].push(entry),
                std::cmp::Ordering::Equal => {}
            }
            lev[c as usize] = INF;
        }
    }
    Pattern { l, u }
}

fn row_norm<T: Scalar>(a: &Csr<T>, i: usize) -> f64 {
    a.row(i).1.iter().map(|v| v.modulus_sq()).sum::<f64>().sqrt()
}

impl<T: Scalar> Ilu<T> {
    /// ILU(`k`): `k = 0` is ILU(0).
    pub fn ilu_k(a: &Csr<T>, k: u32) -> Self {
        assert_eq!(a.nrows(), a.ncols(), "square matrix required");
        let n = a.nrows();
        let pat = symbolic_levels(a, k);
        let mut me = Self::empty(n);
        let mut w = vec![T::ZERO; n];
        let mut stamp = vec![usize::MAX; n];
        for i in 0..n {
            let (cs, vs) = a.row(i);
            for &c in &pat.l[i] {
                stamp[c.0 as usize] = i;
                w[c.0 as usize] = T::ZERO;
            }
            for &c in &pat.u[i] {
                stamp[c.0 as usize] = i;
                w[c.0 as usize] = T::ZERO;
            }
            stamp[i] = i;
            w[i] = T::ZERO;
            for (&c, &v) in cs.iter().zip(vs) {
                w[c as usize] += v;
            }
            for &(j, _) in &pat.l[i] {
                let j = j as usize;
                let lij = w[j] / me.diag[j];
                w[j] = lij;
                for p in me.u_ptr[j]..me.u_ptr[j + 1] {
                    let m = me.u_idx[p] as usize;
                    if stamp[m] == i {
                        let t = lij * me.u_val[p];
                        w[m] -= t;
                    }
                }
            }
            let scale = row_norm(a, i).max(f64::MIN_POSITIVE);
            me.push_row(i, &pat.l[i], &pat.u[i], &w, scale);
        }
        me
    }

    /// ILUT with drop tolerance `droptol` (relative to the 2-norm of the row of `A`) and at most `lfil` kept
    /// entries in the `L` part and in the `U` part of each row.
    pub fn ilut(a: &Csr<T>, droptol: f64, lfil: usize) -> Self {
        assert_eq!(a.nrows(), a.ncols(), "square matrix required");
        let n = a.nrows();
        let mut me = Self::empty(n);
        let mut w = vec![T::ZERO; n];
        let mut present = vec![false; n];
        for i in 0..n {
            let (cs, vs) = a.row(i);
            let scale = row_norm(a, i).max(f64::MIN_POSITIVE);
            let tau = droptol * scale;
            let mut touched: Vec<u32> = Vec::new();
            let mut heap: BinaryHeap<Reverse<u32>> = BinaryHeap::new();
            for (&c, &v) in cs.iter().zip(vs) {
                let cu = c as usize;
                if !present[cu] {
                    present[cu] = true;
                    touched.push(c);
                    if cu < i {
                        heap.push(Reverse(c));
                    }
                }
                w[cu] += v;
            }
            if !present[i] {
                present[i] = true;
                touched.push(i as u32);
            }
            while let Some(Reverse(j)) = heap.pop() {
                let ju = j as usize;
                let lij = w[ju] / me.diag[ju];
                if lij.modulus() < tau {
                    w[ju] = T::ZERO;
                    continue;
                }
                w[ju] = lij;
                for p in me.u_ptr[ju]..me.u_ptr[ju + 1] {
                    let m = me.u_idx[p] as usize;
                    if !present[m] {
                        present[m] = true;
                        touched.push(m as u32);
                        if m < i {
                            heap.push(Reverse(m as u32));
                        }
                    }
                    let t = lij * me.u_val[p];
                    w[m] -= t;
                }
            }
            let mut lo: Vec<(u32, T)> = Vec::new();
            let mut up: Vec<(u32, T)> = Vec::new();
            for &c in &touched {
                let cu = c as usize;
                let v = w[cu];
                if cu < i && v.modulus() > 0.0 {
                    lo.push((c, v));
                } else if cu > i && v.modulus() >= tau {
                    up.push((c, v));
                }
            }
            let keep = |v: &mut Vec<(u32, T)>| {
                if v.len() > lfil {
                    v.sort_by(|a, b| b.1.modulus().total_cmp(&a.1.modulus()));
                    v.truncate(lfil);
                }
                v.sort_by_key(|e| e.0);
            };
            keep(&mut lo);
            keep(&mut up);
            let lpat: Vec<(u32, u32)> = lo.iter().map(|e| (e.0, 0)).collect();
            let upat: Vec<(u32, u32)> = up.iter().map(|e| (e.0, 0)).collect();
            // reload kept values into w (others zero) and store
            for &c in &touched {
                let cu = c as usize;
                if cu != i {
                    w[cu] = T::ZERO;
                }
            }
            for (c, v) in lo.iter().chain(&up) {
                w[*c as usize] = *v;
            }
            me.push_row(i, &lpat, &upat, &w, scale);
            for &c in &touched {
                present[c as usize] = false;
                w[c as usize] = T::ZERO;
            }
        }
        me
    }

    fn empty(n: usize) -> Self {
        Self {
            n,
            l_ptr: vec![0],
            l_idx: Vec::new(),
            l_val: Vec::new(),
            diag: Vec::with_capacity(n),
            u_ptr: vec![0],
            u_idx: Vec::new(),
            u_val: Vec::new(),
            perturbed_pivots: 0,
        }
    }

    fn push_row(&mut self, i: usize, lpat: &[(u32, u32)], upat: &[(u32, u32)], w: &[T], scale: f64) {
        for &(c, _) in lpat {
            self.l_idx.push(c);
            self.l_val.push(w[c as usize]);
        }
        self.l_ptr.push(self.l_idx.len());
        let mut d = w[i];
        let floor = 1e-10 * scale;
        if d.modulus() < floor {
            // keep the sign/phase if there is one
            d = if d.modulus() > 0.0 { d.scale(floor / d.modulus()) } else { T::from_f64(floor) };
            self.perturbed_pivots += 1;
        }
        self.diag.push(d);
        for &(c, _) in upat {
            self.u_idx.push(c);
            self.u_val.push(w[c as usize]);
        }
        self.u_ptr.push(self.u_idx.len());
    }

    pub fn nnz(&self) -> usize {
        self.l_val.len() + self.u_val.len() + self.n
    }

    /// `z ← (LU)⁻¹ r`.
    pub fn solve(&self, r: &[T], z: &mut [T]) {
        let n = self.n;
        // L y = r (unit lower)
        for i in 0..n {
            let mut s = r[i];
            for p in self.l_ptr[i]..self.l_ptr[i + 1] {
                s -= self.l_val[p] * z[self.l_idx[p] as usize];
            }
            z[i] = s;
        }
        // U z = y
        for i in (0..n).rev() {
            let mut s = z[i];
            for p in self.u_ptr[i]..self.u_ptr[i + 1] {
                s -= self.u_val[p] * z[self.u_idx[p] as usize];
            }
            z[i] = s / self.diag[i];
        }
    }
}

impl<T: Scalar> Preconditioner<T> for Ilu<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        self.solve(r, z)
    }
}

// ---------------------------------------------------------------------------------------------
// Incomplete Cholesky
// ---------------------------------------------------------------------------------------------

/// Incomplete Cholesky `A ≈ L Lᴴ` (Hermitian) or `L Lᵀ` (complex symmetric / real symmetric).
#[derive(Debug, Clone)]
pub struct Ic<T: Scalar> {
    n: usize,
    hermitian: bool,
    ptr: Vec<usize>,
    idx: Vec<u32>,
    val: Vec<T>,
    diag: Vec<T>,
    /// Diagonal shift `α` of the successful attempt (`0` if none was needed).
    pub shift: f64,
    pub attempts: usize,
}

impl<T: Scalar> Ic<T> {
    /// IC(`k`) of the lower triangle of `a` (assumed Hermitian if `hermitian`, else symmetric).
    pub fn ic_k(a: &Csr<T>, k: u32, hermitian: bool) -> Self {
        assert_eq!(a.nrows(), a.ncols(), "square matrix required");
        let pat = symbolic_levels(a, k);
        let mut alpha = 0.0;
        let mut attempts = 0;
        loop {
            attempts += 1;
            if let Some(mut f) = Self::try_factor(a, &pat, hermitian, alpha) {
                f.shift = alpha;
                f.attempts = attempts;
                return f;
            }
            alpha = if alpha == 0.0 { 1e-3 } else { alpha * 2.0 };
            assert!(alpha < 1e6, "incomplete Cholesky failed even with a large diagonal shift");
        }
    }

    fn try_factor(a: &Csr<T>, pat: &Pattern, hermitian: bool, alpha: f64) -> Option<Self> {
        let n = a.nrows();
        let mut me = Self {
            n,
            hermitian,
            ptr: vec![0],
            idx: Vec::new(),
            val: Vec::new(),
            diag: Vec::with_capacity(n),
            shift: 0.0,
            attempts: 0,
        };
        let mut w = vec![T::ZERO; n];
        let mut stamp = vec![usize::MAX; n];
        for i in 0..n {
            for &(c, _) in &pat.l[i] {
                stamp[c as usize] = i;
                w[c as usize] = T::ZERO;
            }
            let mut aii = T::ZERO;
            let (cs, vs) = a.row(i);
            for (&c, &v) in cs.iter().zip(vs) {
                let cu = c as usize;
                if cu < i && stamp[cu] == i {
                    w[cu] += v;
                } else if cu == i {
                    aii += v;
                }
            }
            aii += aii.scale(alpha);
            // l_ij = (a_ij − Σ_{m<j} l_im conj(l_jm)) / l_jj, ascending j
            let mut sum_sq = 0.0f64;
            let mut sum_bil = T::ZERO;
            for &(j, _) in &pat.l[i] {
                let j = j as usize;
                let mut s = w[j];
                for p in me.ptr[j]..me.ptr[j + 1] {
                    let m = me.idx[p] as usize;
                    if stamp[m] == i {
                        let ljm = if hermitian { me.val[p].conj() } else { me.val[p] };
                        s -= w[m] * ljm;
                    }
                }
                let lij = s / me.diag[j];
                w[j] = lij;
                if hermitian {
                    sum_sq += lij.modulus_sq();
                } else {
                    sum_bil += lij * lij;
                }
            }
            let d2 = if hermitian { T::from_f64(aii.re() - sum_sq) } else { aii - sum_bil };
            let dii = if hermitian {
                let r = d2.re();
                if r <= 1e-14 * aii.re().abs().max(f64::MIN_POSITIVE) {
                    return None;
                }
                T::from_f64(r.sqrt())
            } else {
                if d2.modulus() <= 1e-14 * aii.modulus().max(f64::MIN_POSITIVE) {
                    return None;
                }
                d2.sqrt_principal()
            };
            for &(c, _) in &pat.l[i] {
                me.idx.push(c);
                me.val.push(w[c as usize]);
            }
            me.ptr.push(me.idx.len());
            me.diag.push(dii);
        }
        Some(me)
    }

    pub fn nnz(&self) -> usize {
        self.val.len() + self.n
    }

    /// `z ← (L Lᴴ)⁻¹ r` (or `(L Lᵀ)⁻¹ r`).
    pub fn solve(&self, r: &[T], z: &mut [T]) {
        let n = self.n;
        // L y = r
        for i in 0..n {
            let mut s = r[i];
            for p in self.ptr[i]..self.ptr[i + 1] {
                s -= self.val[p] * z[self.idx[p] as usize];
            }
            z[i] = s / self.diag[i];
        }
        // Lᴴ z = y : column sweep (the rows of L are the columns of Lᴴ)
        for i in (0..n).rev() {
            z[i] = z[i] / if self.hermitian { self.diag[i].conj() } else { self.diag[i] };
            let zi = z[i];
            for p in self.ptr[i]..self.ptr[i + 1] {
                let j = self.idx[p] as usize;
                let l = if self.hermitian { self.val[p].conj() } else { self.val[p] };
                z[j] -= l * zi;
            }
        }
    }
}

impl<T: Scalar> Preconditioner<T> for Ic<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        self.solve(r, z)
    }
}
