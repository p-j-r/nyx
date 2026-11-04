//! Geometry oracle for 3D: where may a vertex of a surface patch or a crease curve go?
//!
//! By default a patch is the set of triangles of the input mesh and a curve the polyline of its crease edges, so flat faces and straight
//! creases are preserved exactly (a sliding vertex is projected onto the closest point of the patch or curve). A patch can be given an
//! analytic projector, for example onto a sphere or onto the zero set of a signed-distance function, so that new vertices land on the true
//! surface and the discrete domain converges to the exact one under refinement. Closest-point queries run over all triangles of the
//! patch, linear in their number: fine for the sizes of the tests and of one adaptation cycle of a chip; a bounding-volume tree is the
//! next step when a patch has millions of triangles.

use super::topo::Sig;
use super::topo3::{Face, Work3};

type Projector = std::sync::Arc<dyn Fn([f64; 3]) -> [f64; 3] + Send + Sync>;

pub trait Geometry3: Sync {
    /// Closest admissible point of surface patch `patch` to `p`.
    fn project_surface(&self, p: [f64; 3], patch: u32) -> [f64; 3];
    /// Closest admissible point of crease curve `curve` to `p`.
    fn project_crease(&self, p: [f64; 3], curve: u32) -> [f64; 3];
}

/// Triangles and polylines taken from the input mesh, with optional analytic overrides per boundary label.
pub struct FacetGeometry {
    tris: Vec<Vec<[[f64; 3]; 3]>>,
    segs: Vec<Vec<([f64; 3], [f64; 3])>>,
    over: Vec<Option<Projector>>,
    sigs: Vec<Sig>,
}

impl FacetGeometry {
    pub fn from_work(w: &Work3) -> Self {
        let mut tris = vec![Vec::new(); w.patches.len()];
        let mut fs: Vec<(&Face, &u32)> = w.patch_of.iter().collect();
        fs.sort_unstable();
        for (f, &p) in fs {
            tris[p as usize].push([w.pts[f[0] as usize], w.pts[f[1] as usize], w.pts[f[2] as usize]]);
        }
        let mut segs = vec![Vec::new(); w.curves.len()];
        let mut es: Vec<_> = w.curve_of.iter().collect();
        es.sort_unstable();
        for (e, &c) in es {
            segs[c as usize].push((w.pts[e[0] as usize], w.pts[e[1] as usize]));
        }
        FacetGeometry { tris, segs, over: w.patches.iter().map(|_| None).collect(), sigs: w.patches.iter().map(|p| p.sig).collect() }
    }

    /// Use `f` to project every patch of boundary label `label`.
    pub fn with_projector(mut self, label: i32, f: impl Fn([f64; 3]) -> [f64; 3] + Send + Sync + 'static) -> Self {
        let f: Projector = std::sync::Arc::new(f);
        for p in 0..self.sigs.len() {
            if self.sigs[p] == Sig::Bnd(label) {
                self.over[p] = Some(f.clone());
            }
        }
        self
    }
}

impl Geometry3 for FacetGeometry {
    fn project_surface(&self, p: [f64; 3], patch: u32) -> [f64; 3] {
        if let Some(f) = &self.over[patch as usize] {
            return f(p);
        }
        let mut best = (f64::INFINITY, p);
        for t in &self.tris[patch as usize] {
            let q = closest_on_triangle(p, t);
            let r2 = d2(q, p);
            if r2 < best.0 {
                best = (r2, q);
            }
        }
        best.1
    }

    fn project_crease(&self, p: [f64; 3], curve: u32) -> [f64; 3] {
        let mut best = (f64::INFINITY, p);
        for (a, b) in &self.segs[curve as usize] {
            let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let l2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let t = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1] + (p[2] - a[2]) * d[2]) / l2).clamp(0.0, 1.0);
            let q = [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]];
            let r2 = d2(q, p);
            if r2 < best.0 {
                best = (r2, q);
            }
        }
        best.1
    }
}

fn d2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Closest point of triangle `t` to `p` (the case analysis of Ericson, *Real-Time Collision Detection*, 5.1.5).
pub fn closest_on_triangle(p: [f64; 3], t: &[[f64; 3]; 3]) -> [f64; 3] {
    let [a, b, c] = *t;
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2_) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2_ <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2_;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return [a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]];
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2_ - d1 * d6;
    if vb <= 0.0 && d2_ >= 0.0 && d6 <= 0.0 {
        let w = d2_ / (d2_ - d6);
        return [a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]];
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return [b[0] + w * (c[0] - b[0]), b[1] + w * (c[1] - b[1]), b[2] + w * (c[2] - b[2])];
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    [a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w]
}

/// Projector onto the sphere of the given centre and radius.
pub fn sphere_projector(centre: [f64; 3], radius: f64) -> impl Fn([f64; 3]) -> [f64; 3] + Send + Sync + 'static {
    move |p| {
        let d = sub(p, centre);
        let r = dot(d, d).sqrt();
        [centre[0] + d[0] / r * radius, centre[1] + d[1] / r * radius, centre[2] + d[2] / r * radius]
    }
}

/// Projector onto the zero set of a smooth signed-distance-like function (Newton steps along the numerical gradient).
pub fn sdf_projector3(sdf: impl Fn(&[f64; 3]) -> f64 + Send + Sync + 'static) -> impl Fn([f64; 3]) -> [f64; 3] + Send + Sync + 'static {
    move |mut p| {
        for _ in 0..30 {
            let h = 1e-6;
            let f0 = sdf(&p);
            let mut g = [0.0; 3];
            for k in 0..3 {
                let (mut hi, mut lo) = (p, p);
                hi[k] += h;
                lo[k] -= h;
                g[k] = (sdf(&hi) - sdf(&lo)) / (2.0 * h);
            }
            let g2 = dot(g, g);
            if f0.abs() < 1e-14 || g2 == 0.0 {
                break;
            }
            p = [p[0] - f0 * g[0] / g2, p[1] - f0 * g[1] / g2, p[2] - f0 * g[2] / g2];
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapt::field::UniformSize;
    use crate::adapt::topo3::{TopoOptions3, VClass3};
    use crate::mesh::fixtures::box_tet;

    #[test]
    fn triangle_closest_point_covers_face_edge_and_vertex_regions() {
        let t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let q = closest_on_triangle([0.2, 0.2, 5.0], &t);
        assert!(d2(q, [0.2, 0.2, 0.0]) < 1e-28);
        assert_eq!(closest_on_triangle([0.5, -1.0, 0.3], &t), [0.5, 0.0, 0.0]);
        assert_eq!(closest_on_triangle([-1.0, -1.0, 1.0], &t), [0.0, 0.0, 0.0]);
        let q = closest_on_triangle([1.0, 1.0, 0.0], &t);
        assert!((q[0] - 0.5).abs() < 1e-15 && (q[1] - 0.5).abs() < 1e-15);
    }

    #[test]
    fn patches_and_creases_of_a_box_project_onto_the_box() {
        let m = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [3, 3, 3]);
        let w = Work3::from_mesh(&m, &UniformSize(0.3), &TopoOptions3::default()).unwrap();
        let g = FacetGeometry::from_work(&w);
        // a face vertex and a crease vertex of the box
        let fv = (0..w.pts.len()).find(|&v| matches!(w.vclass[v], VClass3::Surface(_)) && w.pts[v][0] == 0.0).unwrap();
        let VClass3::Surface(p) = w.vclass[fv] else { unreachable!() };
        let q = g.project_surface([-0.4, w.pts[fv][1], w.pts[fv][2]], p);
        assert!(q[0].abs() < 1e-15);
        let cv = (0..w.pts.len()).find(|&v| matches!(w.vclass[v], VClass3::Crease(_)) && w.pts[v][0] == 0.0 && w.pts[v][1] == 0.0).unwrap();
        let VClass3::Crease(c) = w.vclass[cv] else { unreachable!() };
        let q = g.project_crease([0.3, 0.2, w.pts[cv][2]], c);
        assert!(q[0].abs() < 1e-15 && q[1].abs() < 1e-15);
    }

    #[test]
    fn sphere_and_sdf_projectors_agree() {
        let s = sphere_projector([0.5, 0.0, -1.0], 2.0);
        let f = sdf_projector3(|p| ((p[0] - 0.5).powi(2) + p[1].powi(2) + (p[2] + 1.0).powi(2)).sqrt() - 2.0);
        for p in [[1.0, 1.0, 1.0], [-3.0, 0.2, 0.0], [0.5, 3.1, -2.0]] {
            let (a, b) = (s(p), f(p));
            assert!(d2(a, b) < 1e-16);
        }
    }
}
