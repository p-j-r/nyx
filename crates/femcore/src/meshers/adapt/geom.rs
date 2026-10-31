//! Geometry oracle for 2D: where may a vertex of a boundary or interface chain go?
//!
//! By default a chain is the polyline of the input mesh, so straight sides and interfaces are
//! preserved exactly (sliding vertices are projected onto the polyline). For smooth boundaries a
//! chain can be given an analytic projector, for example onto a circle or onto the zero set of a
//! signed-distance function, so that new vertices land on the true curve and the discrete domain
//! converges to the exact one under refinement.

use super::topo::{Sig, Work2};

type Projector = std::sync::Arc<dyn Fn([f64; 2]) -> [f64; 2] + Send + Sync>;

pub trait Geometry2: Sync {
    /// Closest admissible point of chain `chain` to `p`.
    fn project(&self, p: [f64; 2], chain: u32) -> [f64; 2];
}

/// Polyline chains taken from the input mesh, with optional analytic overrides per boundary label.
pub struct PolylineGeometry {
    segs: Vec<Vec<([f64; 2], [f64; 2])>>,
    over: Vec<Option<Projector>>,
    sigs: Vec<Sig>,
}

impl PolylineGeometry {
    pub fn from_work(w: &Work2) -> Self {
        let mut segs = vec![Vec::new(); w.chains.len()];
        let mut es: Vec<_> = w.chain_of.iter().collect();
        es.sort_unstable();
        for (e, &c) in es {
            segs[c as usize].push((w.pts[e[0] as usize], w.pts[e[1] as usize]));
        }
        PolylineGeometry {
            segs,
            over: w.chains.iter().map(|_| None).collect(),
            sigs: w.chains.iter().map(|c| c.sig).collect(),
        }
    }

    /// Use `f` to project every chain of boundary label `label`.
    pub fn with_projector(mut self, label: i32, f: impl Fn([f64; 2]) -> [f64; 2] + Send + Sync + 'static) -> Self {
        let f: Projector = std::sync::Arc::new(f);
        for c in 0..self.sigs.len() {
            if self.sigs[c] == Sig::Bnd(label) {
                self.over[c] = Some(f.clone());
            }
        }
        self
    }
}

impl Geometry2 for PolylineGeometry {
    fn project(&self, p: [f64; 2], chain: u32) -> [f64; 2] {
        if let Some(f) = &self.over[chain as usize] {
            return f(p);
        }
        let mut best = (f64::INFINITY, p);
        for (a, b) in &self.segs[chain as usize] {
            let d = [b[0] - a[0], b[1] - a[1]];
            let l2 = d[0] * d[0] + d[1] * d[1];
            let t = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / l2).clamp(0.0, 1.0);
            let q = [a[0] + t * d[0], a[1] + t * d[1]];
            let r2 = (q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2);
            if r2 < best.0 {
                best = (r2, q);
            }
        }
        best.1
    }
}

/// Projector onto the circle of the given centre and radius.
pub fn circle_projector(centre: [f64; 2], radius: f64) -> impl Fn([f64; 2]) -> [f64; 2] + Send + Sync + 'static {
    move |p| {
        let d = [p[0] - centre[0], p[1] - centre[1]];
        let r = d[0].hypot(d[1]);
        [centre[0] + d[0] / r * radius, centre[1] + d[1] / r * radius]
    }
}

/// Projector onto the zero set of a smooth signed-distance-like function (Newton steps along the
/// numerical gradient).
pub fn sdf_projector(
    sdf: impl Fn(&[f64; 2]) -> f64 + Send + Sync + 'static,
) -> impl Fn([f64; 2]) -> [f64; 2] + Send + Sync + 'static {
    move |mut p| {
        for _ in 0..30 {
            let h = 1e-6;
            let f0 = sdf(&p);
            let gx = (sdf(&[p[0] + h, p[1]]) - sdf(&[p[0] - h, p[1]])) / (2.0 * h);
            let gy = (sdf(&[p[0], p[1] + h]) - sdf(&[p[0], p[1] - h])) / (2.0 * h);
            let g2 = gx * gx + gy * gy;
            if f0.abs() < 1e-14 || g2 == 0.0 {
                break;
            }
            p = [p[0] - f0 * gx / g2, p[1] - f0 * gy / g2];
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshers::adapt::field::UniformSize;
    use crate::meshers::adapt::topo::TopoOptions;
    use crate::meshers::fixtures::{box_tri, side_label};

    #[test]
    fn polyline_projects_onto_sides() {
        let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 4, 4);
        let w = Work2::from_mesh(&m, &UniformSize(0.2), &TopoOptions::default()).unwrap();
        let g = PolylineGeometry::from_work(&w);
        let c = w.chains.iter().position(|c| c.sig == Sig::Bnd(side_label(0, true))).unwrap() as u32;
        let q = g.project([1.3, 0.37], c);
        assert!((q[0] - 1.0).abs() < 1e-15 && (q[1] - 0.37).abs() < 1e-15);
        let q = g.project([1.3, 1.9], c);
        assert!((q[0] - 1.0).abs() < 1e-15 && (q[1] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn circle_and_sdf_projectors_agree() {
        let c = circle_projector([0.5, 0.5], 2.0);
        let s = sdf_projector(|p| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - 2.0);
        for p in [[1.0, 1.0], [-3.0, 0.2], [0.5, 3.1]] {
            let (a, b) = (c(p), s(p));
            assert!((a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9);
            assert!((((a[0] - 0.5).powi(2) + (a[1] - 0.5).powi(2)).sqrt() - 2.0).abs() < 1e-12);
        }
    }
}
