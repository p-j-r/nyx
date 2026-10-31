//! Reverse Cuthill-McKee ordering and a banded LU solver for the harness.
//!
//! The LU does not pivot. That is adequate for the PML Helmholtz systems here (complex
//! coefficients keep the pivots away from zero), and the solver refuses to continue on a tiny
//! pivot; callers should still check the residual, which [`BandMatrix::solve`] returns.

use super::cplx::{C64, ZERO};

/// Reverse Cuthill-McKee permutation: `perm[new] = old`.
pub fn rcm(adj: &[Vec<u32>]) -> Vec<u32> {
    let n = adj.len();
    let deg: Vec<usize> = adj.iter().map(|a| a.len()).collect();
    let mut visited = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut starts: Vec<usize> = (0..n).collect();
    starts.sort_by_key(|&i| (deg[i], i));
    for &s in &starts {
        if visited[s] {
            continue;
        }
        // pseudo-peripheral start: repeat BFS from the farthest minimum-degree node
        let mut root = s;
        for _ in 0..3 {
            let (last, _) = bfs_last(adj, &deg, root);
            if last == root {
                break;
            }
            root = last;
        }
        let mut queue = std::collections::VecDeque::new();
        visited[root] = true;
        queue.push_back(root);
        while let Some(u) = queue.pop_front() {
            order.push(u as u32);
            let mut nb: Vec<u32> = adj[u].iter().cloned().filter(|&v| !visited[v as usize]).collect();
            nb.sort_by_key(|&v| (deg[v as usize], v));
            for v in nb {
                if !visited[v as usize] {
                    visited[v as usize] = true;
                    queue.push_back(v as usize);
                }
            }
        }
    }
    order.reverse();
    order
}

fn bfs_last(adj: &[Vec<u32>], deg: &[usize], root: usize) -> (usize, usize) {
    let mut seen = std::collections::HashSet::new();
    let mut level = vec![root];
    seen.insert(root);
    let mut last_level = level.clone();
    let mut depth = 0;
    while !level.is_empty() {
        last_level = level.clone();
        let mut next = Vec::new();
        for &u in &level {
            for &v in &adj[u] {
                if seen.insert(v as usize) {
                    next.push(v as usize);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        level = next;
        depth += 1;
    }
    let best = last_level.into_iter().min_by_key(|&i| (deg[i], i)).unwrap();
    (best, depth)
}

/// Square band matrix with half-bandwidth `bw`.
pub struct BandMatrix {
    n: usize,
    bw: usize,
    /// row `i` holds columns `i - bw ..= i + bw` at `data[i * (2 bw + 1) + (j + bw - i)]`
    data: Vec<C64>,
}

impl BandMatrix {
    pub fn new(n: usize, bw: usize) -> Self {
        BandMatrix { n, bw, data: vec![ZERO; n * (2 * bw + 1)] }
    }

    pub fn bandwidth(&self) -> usize {
        self.bw
    }

    fn idx(&self, i: usize, j: usize) -> usize {
        debug_assert!(j + self.bw >= i && j <= i + self.bw, "entry ({i},{j}) outside band {}", self.bw);
        i * (2 * self.bw + 1) + (j + self.bw - i)
    }

    pub fn add(&mut self, i: usize, j: usize, v: C64) {
        let k = self.idx(i, j);
        self.data[k] += v;
    }

    pub fn get(&self, i: usize, j: usize) -> C64 {
        if j + self.bw < i || j > i + self.bw { ZERO } else { self.data[self.idx(i, j)] }
    }

    /// `y = A x`.
    pub fn mul_vec(&self, x: &[C64]) -> Vec<C64> {
        let mut y = vec![ZERO; self.n];
        for i in 0..self.n {
            let (lo, hi) = (i.saturating_sub(self.bw), (i + self.bw).min(self.n - 1));
            let mut s = ZERO;
            for j in lo..=hi {
                s += self.data[self.idx(i, j)] * x[j];
            }
            y[i] = s;
        }
        y
    }

    /// `self += alpha * other` (same size and bandwidth).
    pub fn axpy(&mut self, alpha: C64, other: &BandMatrix) {
        assert!(self.n == other.n && self.bw == other.bw);
        for (x, y) in self.data.iter_mut().zip(&other.data) {
            *x += alpha * *y;
        }
    }

    /// LU factorisation (no pivoting) that can be reused for several right-hand sides.
    pub fn factor(&self) -> Result<BandLu, String> {
        let mut a = BandMatrix { n: self.n, bw: self.bw, data: self.data.clone() };
        let (n, bw) = (a.n, a.bw);
        for k in 0..n {
            let piv = a.data[a.idx(k, k)];
            let scale =
                (k.saturating_sub(bw)..=(k + bw).min(n - 1)).map(|j| a.data[a.idx(k, j)].abs()).fold(0.0, f64::max);
            if piv.abs() <= 1e-14 * scale || piv.abs() == 0.0 {
                return Err(format!("tiny pivot at row {k}"));
            }
            let jmax = (k + bw).min(n - 1);
            for i in k + 1..=(k + bw).min(n - 1) {
                let aik = a.data[a.idx(i, k)];
                if aik.re == 0.0 && aik.im == 0.0 {
                    continue;
                }
                let f = aik / piv;
                // store the multiplier in the eliminated slot
                let ik = a.idx(i, k);
                a.data[ik] = f;
                for j in k + 1..=jmax {
                    let t = a.data[a.idx(k, j)] * f;
                    let ij = a.idx(i, j);
                    a.data[ij] -= t;
                }
            }
        }
        Ok(BandLu { lu: a })
    }

    /// Solve `A x = b`; returns `(x, relative residual)`.
    pub fn solve(&self, b: &[C64]) -> Result<(Vec<C64>, f64), String> {
        let x = self.factor()?.solve(b);
        let r = self.mul_vec(&x);
        let num: f64 = r.iter().zip(b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
        let den: f64 = b.iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt().max(1e-300);
        Ok((x, num / den))
    }
}

/// Factorisation produced by [`BandMatrix::factor`].
pub struct BandLu {
    lu: BandMatrix,
}

impl BandLu {
    pub fn solve(&self, b: &[C64]) -> Vec<C64> {
        let a = &self.lu;
        let (n, bw) = (a.n, a.bw);
        let mut x = b.to_vec();
        for k in 0..n {
            for i in k + 1..=(k + bw).min(n - 1) {
                let f = a.data[a.idx(i, k)];
                if f.re != 0.0 || f.im != 0.0 {
                    x[i] = x[i] - x[k] * f;
                }
            }
        }
        for k in (0..n).rev() {
            let mut s = x[k];
            for j in k + 1..=(k + bw).min(n - 1) {
                s -= a.data[a.idx(k, j)] * x[j];
            }
            x[k] = s / a.data[a.idx(k, k)];
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rcm_is_a_permutation_and_reduces_bandwidth_of_a_grid() {
        let (nx, ny) = (30usize, 5usize);
        // number the long direction slowly so natural bandwidth is large
        let id = |i: usize, j: usize| (j * nx + i) as u32;
        let mut adj = vec![Vec::new(); nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                for (di, dj) in [(1i32, 0i32), (0, 1)] {
                    let (a, b) = (i as i32 + di, j as i32 + dj);
                    if (a as usize) < nx && (b as usize) < ny {
                        adj[id(i, j) as usize].push(id(a as usize, b as usize));
                        adj[id(a as usize, b as usize) as usize].push(id(i, j));
                    }
                }
            }
        }
        let perm = rcm(&adj);
        let mut seen = vec![false; nx * ny];
        for &p in &perm {
            assert!(!seen[p as usize]);
            seen[p as usize] = true;
        }
        let mut inv = vec![0usize; perm.len()];
        for (new, &old) in perm.iter().enumerate() {
            inv[old as usize] = new;
        }
        let bw = |m: &dyn Fn(usize) -> usize| -> usize {
            (0..nx * ny)
                .flat_map(|u| adj[u].iter().map(move |&v| (u, v as usize)))
                .map(|(u, v)| m(u).abs_diff(m(v)))
                .max()
                .unwrap()
        };
        let natural = bw(&|u| u);
        let reordered = bw(&|u| inv[u]);
        assert!(reordered <= ny + 1 + 1, "{reordered} vs natural {natural}");
        assert!(reordered < natural);
    }

    #[test]
    fn band_solver_matches_known_solution() {
        let n = 200;
        let bw = 3;
        let mut a = BandMatrix::new(n, bw);
        let mut s = 12345u64;
        let mut rnd = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
        };
        for i in 0..n {
            for j in i.saturating_sub(bw)..=(i + bw).min(n - 1) {
                let mut v = C64::new(rnd(), rnd());
                if i == j {
                    v += C64::new(4.0, 1.0);
                }
                a.add(i, j, v);
            }
        }
        let xs: Vec<C64> = (0..n).map(|i| C64::new((i as f64).sin(), (i as f64 * 0.3).cos())).collect();
        let b = a.mul_vec(&xs);
        let (x, res) = a.solve(&b).unwrap();
        assert!(res < 1e-12);
        for (p, q) in x.iter().zip(&xs) {
            assert!((*p - *q).abs() < 1e-10);
        }
    }
}
