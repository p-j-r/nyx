//! Metric fields: anything that can return the requested metric at a point.
//!
//! `BackgroundMetric2` interpolates a vertex-based metric defined on a (typically coarser or
//! previous) triangle mesh: log-Euclidean interpolation inside the containing triangle, found
//! through a uniform bucket grid. This is how a metric computed on one mesh drives the adaptation
//! of the next one.

use super::metric::{Sym, exp_sym, isotropic, log_spd, symmetrize};
use crate::meshers::mesh::SimplexMesh;

/// A metric defined everywhere in the domain.
pub trait MetricField<const D: usize>: Sync {
    fn at(&self, x: &[f64; D]) -> Sym<D>;
}

impl<const D: usize, F: Fn(&[f64; D]) -> Sym<D> + Sync> MetricField<D> for F {
    fn at(&self, x: &[f64; D]) -> Sym<D> {
        self(x)
    }
}

/// Uniform isotropic size `h` everywhere.
pub struct UniformSize(pub f64);

impl<const D: usize> MetricField<D> for UniformSize {
    fn at(&self, _: &[f64; D]) -> Sym<D> {
        isotropic::<D>(self.0)
    }
}

/// Metric interpolated from per-vertex values on a triangle mesh.
pub struct BackgroundMetric2 {
    pts: Vec<[f64; 2]>,
    tris: Vec<[u32; 3]>,
    logs: Vec<Sym<2>>,
    lo: [f64; 2],
    cell: f64,
    nx: usize,
    ny: usize,
    buckets: Vec<Vec<u32>>,
}

impl BackgroundMetric2 {
    pub fn new(mesh: &SimplexMesh<2>, metrics: &[Sym<2>]) -> Self {
        assert_eq!(mesh.n_points(), metrics.len());
        let bb = mesh.bbox();
        let nt = mesh.n_cells().max(1);
        let area = ((bb[0][1] - bb[0][0]) * (bb[1][1] - bb[1][0])).max(1e-300);
        // about 2 triangles per bucket
        let cell = (2.0 * area / nt as f64).sqrt().max(1e-12);
        let nx = (((bb[0][1] - bb[0][0]) / cell).ceil() as usize).max(1);
        let ny = (((bb[1][1] - bb[1][0]) / cell).ceil() as usize).max(1);
        let tris: Vec<[u32; 3]> =
            (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2]]).collect();
        let mut buckets = vec![Vec::new(); nx * ny];
        for (t, tr) in tris.iter().enumerate() {
            let p: Vec<[f64; 2]> = tr.iter().map(|&v| mesh.points[v as usize]).collect();
            let (x0, x1) = (
                p.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min),
                p.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max),
            );
            let (y0, y1) = (
                p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min),
                p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max),
            );
            let ix0 = (((x0 - bb[0][0]) / cell) as usize).min(nx - 1);
            let ix1 = (((x1 - bb[0][0]) / cell) as usize).min(nx - 1);
            let iy0 = (((y0 - bb[1][0]) / cell) as usize).min(ny - 1);
            let iy1 = (((y1 - bb[1][0]) / cell) as usize).min(ny - 1);
            for iy in iy0..=iy1 {
                for ix in ix0..=ix1 {
                    buckets[iy * nx + ix].push(t as u32);
                }
            }
        }
        BackgroundMetric2 {
            pts: mesh.points.clone(),
            tris,
            logs: metrics.iter().map(log_spd).collect(),
            lo: [bb[0][0], bb[1][0]],
            cell,
            nx,
            ny,
            buckets,
        }
    }

    /// Like [`BackgroundMetric2::new`] but first grades the vertex metrics so that the requested
    /// size grows by at most `g` per unit distance (see `metric::grade`). Adapting to an ungraded
    /// field with sharp features is unreliable; this is the recommended constructor.
    pub fn new_graded(mesh: &SimplexMesh<2>, metrics: &[Sym<2>], g: f64) -> Self {
        let mut m = metrics.to_vec();
        super::metric::grade(&mut m, &mesh.points, &mesh.edges(), g, 100);
        Self::new(mesh, &m)
    }

    fn scan(&self, ix: usize, iy: usize, x: &[f64; 2], best: &mut Option<(f64, usize, [f64; 3])>) {
        for &t in &self.buckets[iy * self.nx + ix] {
            let b = self.bary(t as usize, x);
            let m = b[0].min(b[1]).min(b[2]);
            if best.as_ref().is_none_or(|(bm, bt, _)| m > *bm || (m == *bm && (t as usize) < *bt)) {
                *best = Some((m, t as usize, b));
            }
        }
    }

    fn bary(&self, t: usize, x: &[f64; 2]) -> [f64; 3] {
        let tr = self.tris[t];
        let (a, b, c) = (self.pts[tr[0] as usize], self.pts[tr[1] as usize], self.pts[tr[2] as usize]);
        let det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        let l1 = ((x[0] - a[0]) * (c[1] - a[1]) - (x[1] - a[1]) * (c[0] - a[0])) / det;
        let l2 = ((b[0] - a[0]) * (x[1] - a[1]) - (b[1] - a[1]) * (x[0] - a[0])) / det;
        [1.0 - l1 - l2, l1, l2]
    }
}

impl MetricField<2> for BackgroundMetric2 {
    fn at(&self, x: &[f64; 2]) -> Sym<2> {
        let ix = (((x[0] - self.lo[0]) / self.cell).max(0.0) as usize).min(self.nx - 1);
        let iy = (((x[1] - self.lo[1]) / self.cell).max(0.0) as usize).min(self.ny - 1);
        // best triangle = the one whose smallest barycentric coordinate is largest (>= 0 inside)
        let mut best: Option<(f64, usize, [f64; 3])> = None;
        self.scan(ix, iy, x, &mut best);
        if best.as_ref().is_none_or(|b| b.0 < -1e-12) {
            for jy in iy.saturating_sub(1)..=(iy + 1).min(self.ny - 1) {
                for jx in ix.saturating_sub(1)..=(ix + 1).min(self.nx - 1) {
                    self.scan(jx, jy, x, &mut best);
                }
            }
        }
        let (_, t, b) = best.expect("background mesh is empty");
        let b = if b.iter().any(|&w| w < 0.0) {
            // outside the mesh: clamp to the closest point of the triangle's barycentric range
            let c = [b[0].max(0.0), b[1].max(0.0), b[2].max(0.0)];
            let s: f64 = c.iter().sum();
            [c[0] / s, c[1] / s, c[2] / s]
        } else {
            b
        };
        let tr = self.tris[t];
        let mut acc = [[0.0; 2]; 2];
        for k in 0..3 {
            let l = &self.logs[tr[k] as usize];
            for i in 0..2 {
                for j in 0..2 {
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
    use crate::meshers::adapt::metric::sizes;
    use crate::meshers::fixtures::box_tri;

    #[test]
    fn graded_constructor_limits_size_growth() {
        let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 20, 20);
        let ms: Vec<_> = m.points.iter().map(|p| isotropic::<2>(if p[0] < 0.05 { 0.005 } else { 0.5 })).collect();
        let bg = BackgroundMetric2::new_graded(&m, &ms, 0.5);
        let h = |x: f64| sizes(&bg.at(&[x, 0.5])).0[0];
        assert!((h(0.0) - 0.005).abs() < 1e-9);
        // 0.3 away from the fine region the size can be at most 0.005 + 0.5 * 0.3 (plus interpolation slack)
        assert!(h(0.35) <= 0.005 + 0.5 * 0.35 + 1e-9, "{}", h(0.35));
        assert!(h(0.35) > h(0.1));
    }

    #[test]
    fn background_reproduces_vertex_values_and_log_linear_fields() {
        let m = box_tri([[0.0, 1.0], [0.0, 2.0]], 10, 20);
        // sizes h(x, y) = exp(-1 - x): log-linear, so log-Euclidean interpolation is exact
        let h = |x: &[f64; 2]| (-1.0 - x[0]).exp();
        let ms: Vec<_> = m.points.iter().map(|p| isotropic::<2>(h(p))).collect();
        let bg = BackgroundMetric2::new(&m, &ms);
        for p in [[0.123, 0.456], [0.999, 1.999], [0.0, 0.0], [0.5, 1.0], [0.77, 0.31]] {
            let (hh, _) = sizes(&bg.at(&p));
            assert!((hh[0] - h(&p)).abs() < 1e-12 && (hh[1] - h(&p)).abs() < 1e-12, "{p:?}");
        }
        // outside: clamped, finite
        let (hh, _) = sizes(&bg.at(&[3.0, 5.0]));
        assert!(hh[0].is_finite());
    }
}
