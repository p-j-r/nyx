//! Metric field of a tetrahedral background mesh: per-vertex metrics interpolated log-Euclideanly inside the tetrahedra (the 3D
//! counterpart of `BackgroundMetric2`), used to carry the metric of a start mesh to the vertices that adaptation creates.

use super::field::MetricField;
use super::metric::{exp_sym, grade, log_spd, symmetrize, Sym};
use crate::meshers::mesh::SimplexMesh;

pub struct BackgroundMetric3 {
    pts: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
    logs: Vec<Sym<3>>,
    lo: [f64; 3],
    cell: f64,
    n: [usize; 3],
    buckets: Vec<Vec<u32>>,
}

impl BackgroundMetric3 {
    pub fn new(mesh: &SimplexMesh<3>, metrics: &[Sym<3>]) -> Self {
        assert_eq!(mesh.n_points(), metrics.len());
        let bb = mesh.bbox();
        let nt = mesh.n_cells().max(1);
        let vol = ((bb[0][1] - bb[0][0]) * (bb[1][1] - bb[1][0]) * (bb[2][1] - bb[2][0])).max(1e-300);
        // about 2 tetrahedra per bucket
        let cell = (2.0 * vol / nt as f64).cbrt().max(1e-12);
        let n: [usize; 3] = std::array::from_fn(|k| (((bb[k][1] - bb[k][0]) / cell).ceil() as usize).max(1));
        let tets: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
        let mut buckets = vec![Vec::new(); n[0] * n[1] * n[2]];
        for (t, tt) in tets.iter().enumerate() {
            let mut lo = [0usize; 3];
            let mut hi = [0usize; 3];
            for k in 0..3 {
                let a = tt.iter().map(|&v| mesh.points[v as usize][k]).fold(f64::INFINITY, f64::min);
                let b = tt.iter().map(|&v| mesh.points[v as usize][k]).fold(f64::NEG_INFINITY, f64::max);
                lo[k] = (((a - bb[k][0]) / cell) as usize).min(n[k] - 1);
                hi[k] = (((b - bb[k][0]) / cell) as usize).min(n[k] - 1);
            }
            for iz in lo[2]..=hi[2] {
                for iy in lo[1]..=hi[1] {
                    for ix in lo[0]..=hi[0] {
                        buckets[(iz * n[1] + iy) * n[0] + ix].push(t as u32);
                    }
                }
            }
        }
        BackgroundMetric3 { pts: mesh.points.clone(), tets, logs: metrics.iter().map(log_spd).collect(), lo: [bb[0][0], bb[1][0], bb[2][0]], cell, n, buckets }
    }

    /// Like [`BackgroundMetric3::new`] but first grades the vertex metrics so that the requested size grows by at most `g` per unit
    /// distance. Adapting to an ungraded field with sharp features is unreliable; this is the recommended constructor.
    pub fn new_graded(mesh: &SimplexMesh<3>, metrics: &[Sym<3>], g: f64) -> Self {
        let mut m = metrics.to_vec();
        grade(&mut m, &mesh.points, &mesh.edges(), g, 100);
        Self::new(mesh, &m)
    }

    fn bary(&self, t: usize, x: &[f64; 3]) -> [f64; 4] {
        let tt = self.tets[t];
        let p: [[f64; 3]; 4] = std::array::from_fn(|i| self.pts[tt[i] as usize]);
        let col = |k: usize| [p[k][0] - p[0][0], p[k][1] - p[0][1], p[k][2] - p[0][2]];
        let (u, v, w) = (col(1), col(2), col(3));
        let r = [x[0] - p[0][0], x[1] - p[0][1], x[2] - p[0][2]];
        let det3 = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        let d = det3(u, v, w);
        let (l1, l2, l3) = (det3(r, v, w) / d, det3(u, r, w) / d, det3(u, v, r) / d);
        [1.0 - l1 - l2 - l3, l1, l2, l3]
    }

    fn scan(&self, idx: usize, x: &[f64; 3], best: &mut Option<(f64, usize, [f64; 4])>) {
        for &t in &self.buckets[idx] {
            let b = self.bary(t as usize, x);
            let m = b.iter().cloned().fold(f64::INFINITY, f64::min);
            if best.as_ref().is_none_or(|(bm, bt, _)| m > *bm || (m == *bm && (t as usize) < *bt)) {
                *best = Some((m, t as usize, b));
            }
        }
    }
}

impl MetricField<3> for BackgroundMetric3 {
    fn at(&self, x: &[f64; 3]) -> Sym<3> {
        let i: [usize; 3] = std::array::from_fn(|k| (((x[k] - self.lo[k]) / self.cell).max(0.0) as usize).min(self.n[k] - 1));
        let mut best: Option<(f64, usize, [f64; 4])> = None;
        self.scan((i[2] * self.n[1] + i[1]) * self.n[0] + i[0], x, &mut best);
        if best.as_ref().is_none_or(|b| b.0 < -1e-12) {
            for jz in i[2].saturating_sub(1)..=(i[2] + 1).min(self.n[2] - 1) {
                for jy in i[1].saturating_sub(1)..=(i[1] + 1).min(self.n[1] - 1) {
                    for jx in i[0].saturating_sub(1)..=(i[0] + 1).min(self.n[0] - 1) {
                        self.scan((jz * self.n[1] + jy) * self.n[0] + jx, x, &mut best);
                    }
                }
            }
        }
        let (_, t, b) = best.expect("background mesh is empty");
        let b = if b.iter().any(|&w| w < 0.0) {
            // outside the mesh: clamp to the closest point of the barycentric range
            let c = b.map(|w| w.max(0.0));
            let s: f64 = c.iter().sum();
            c.map(|w| w / s)
        } else {
            b
        };
        let tt = self.tets[t];
        let mut acc = [[0.0; 3]; 3];
        for k in 0..4 {
            let l = &self.logs[tt[k] as usize];
            for i in 0..3 {
                for j in 0..3 {
                    acc[i][j] += b[k] * l[i][j];
                }
            }
        }
        exp_sym(&symmetrize(&acc))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshers::adapt::metric::{isotropic, sizes};
    use crate::meshers::fixtures::box_tet;

    #[test]
    fn background_reproduces_log_linear_fields_and_clamps_outside() {
        let m = box_tet([[0.0, 1.0], [0.0, 2.0], [0.0, 1.0]], [4, 8, 4]);
        // sizes h = exp(-1 - x): log-linear, so log-Euclidean interpolation is exact
        let h = |x: &[f64; 3]| (-1.0 - x[0]).exp();
        let ms: Vec<_> = m.points.iter().map(|p| isotropic::<3>(h(p))).collect();
        let bg = BackgroundMetric3::new(&m, &ms);
        for p in [[0.123, 0.456, 0.3], [0.999, 1.999, 0.99], [0.0, 0.0, 0.0], [0.5, 1.0, 0.5], [0.77, 0.31, 0.12]] {
            let (hh, _) = sizes(&bg.at(&p));
            assert!((0..3).all(|k| (hh[k] - h(&p)).abs() < 1e-12), "{p:?}");
        }
        assert!(sizes(&bg.at(&[3.0, 5.0, -2.0])).0[0].is_finite());
    }

    #[test]
    fn graded_constructor_limits_size_growth() {
        let m = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [10, 10, 10]);
        let ms: Vec<_> = m.points.iter().map(|p| isotropic::<3>(if p[0] < 0.05 { 0.005 } else { 0.5 })).collect();
        let bg = BackgroundMetric3::new_graded(&m, &ms, 0.5);
        let h = |x: f64| sizes(&bg.at(&[x, 0.5, 0.5])).0[0];
        assert!((h(0.0) - 0.005).abs() < 1e-9);
        assert!(h(0.35) <= 0.005 + 0.5 * 0.35 + 1e-9 && h(0.35) > h(0.1));
    }
}
