//! Test fixtures: structured box meshes, labelled like the output of a CAD mesher, and analytic
//! scalar fields with exact gradients and Hessians.

use crate::adapt::metric::Sym;
use crate::mesh::native::SimplexMesh;

/// Label of the boundary facet on the minimum / maximum side of axis `a`: `1 + 2a + side`.
pub fn side_label(axis: usize, max_side: bool) -> i32 {
    1 + 2 * axis as i32 + max_side as i32
}

/// Tensor-product triangulation of the rectangle spanned by the coordinate lists `xs`, `ys`
/// (strictly increasing). Every cell is cut along the diagonal through its smallest vertex id, so
/// the result is symmetric under the usual index conventions. Cell label 1; boundary facet labels
/// from [`side_label`].
pub fn rect_tri(xs: &[f64], ys: &[f64]) -> SimplexMesh<2> {
    let (nx, ny) = (xs.len(), ys.len());
    let mut m = SimplexMesh::<2>::default();
    for j in 0..ny {
        for i in 0..nx {
            m.points.push([xs[i], ys[j]]);
        }
    }
    let id = |i: usize, j: usize| (j * nx + i) as u32;
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let q = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)];
            let [t0, t1] = crate::mesh::cells::split_quad(&q);
            m.push_cell(&t0, 1);
            m.push_cell(&t1, 1);
        }
    }
    for i in 0..nx - 1 {
        m.push_facet(&[id(i, 0), id(i + 1, 0)], side_label(1, false));
        m.push_facet(&[id(i, ny - 1), id(i + 1, ny - 1)], side_label(1, true));
    }
    for j in 0..ny - 1 {
        m.push_facet(&[id(0, j), id(0, j + 1)], side_label(0, false));
        m.push_facet(&[id(nx - 1, j), id(nx - 1, j + 1)], side_label(0, true));
    }
    m
}

/// Uniform `nx x ny` cells over `[x0, x1] x [y0, y1]`.
pub fn box_tri(b: [[f64; 2]; 2], nx: usize, ny: usize) -> SimplexMesh<2> {
    let lin = |lo: f64, hi: f64, n: usize| (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).collect::<Vec<_>>();
    rect_tri(&lin(b[0][0], b[0][1], nx), &lin(b[1][0], b[1][1], ny))
}

/// Tensor-product tetrahedral mesh (Kuhn / Freudenthal split of every brick along its main
/// diagonal, which is conforming across bricks without any extra rule). Cell label 1; boundary
/// facet labels from [`side_label`].
pub fn brick_tet(xs: &[f64], ys: &[f64], zs: &[f64]) -> SimplexMesh<3> {
    let (nx, ny, nz) = (xs.len(), ys.len(), zs.len());
    let mut m = SimplexMesh::<3>::default();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                m.points.push([xs[i], ys[j], zs[k]]);
            }
        }
    }
    let id = |i: usize, j: usize, k: usize| ((k * ny + j) * nx + i) as u32;
    // the six monotone paths from corner (0,0,0) to (1,1,1)
    const PERMS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    for k in 0..nz - 1 {
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                for perm in PERMS {
                    let mut c = [0usize; 3];
                    let mut vs = [id(i, j, k); 4];
                    for (s, &ax) in perm.iter().enumerate() {
                        c[ax] = 1;
                        vs[s + 1] = id(i + c[0], j + c[1], k + c[2]);
                    }
                    m.push_cell(&vs, 1);
                }
            }
        }
    }
    m.orient_positive();
    // boundary facets: each boundary brick face is a quad -> 2 triangles via the split rule
    let mut quad = |q: [u32; 4], lab: i32| {
        let [t0, t1] = crate::mesh::cells::split_quad(&q);
        m.push_facet(&t0, lab);
        m.push_facet(&t1, lab);
    };
    let mut pending: Vec<([u32; 4], i32)> = Vec::new();
    for k in 0..nz - 1 {
        for j in 0..ny - 1 {
            pending.push(([id(0, j, k), id(0, j + 1, k), id(0, j + 1, k + 1), id(0, j, k + 1)], side_label(0, false)));
            pending.push((
                [id(nx - 1, j, k), id(nx - 1, j + 1, k), id(nx - 1, j + 1, k + 1), id(nx - 1, j, k + 1)],
                side_label(0, true),
            ));
        }
    }
    for k in 0..nz - 1 {
        for i in 0..nx - 1 {
            pending.push(([id(i, 0, k), id(i + 1, 0, k), id(i + 1, 0, k + 1), id(i, 0, k + 1)], side_label(1, false)));
            pending.push((
                [id(i, ny - 1, k), id(i + 1, ny - 1, k), id(i + 1, ny - 1, k + 1), id(i, ny - 1, k + 1)],
                side_label(1, true),
            ));
        }
    }
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            pending.push(([id(i, j, 0), id(i + 1, j, 0), id(i + 1, j + 1, 0), id(i, j + 1, 0)], side_label(2, false)));
            pending.push((
                [id(i, j, nz - 1), id(i + 1, j, nz - 1), id(i + 1, j + 1, nz - 1), id(i, j + 1, nz - 1)],
                side_label(2, true),
            ));
        }
    }
    for (q, lab) in pending {
        quad(q, lab);
    }
    m
}

/// Uniform `nx x ny x nz` bricks over the box `b`.
pub fn box_tet(b: [[f64; 2]; 3], n: [usize; 3]) -> SimplexMesh<3> {
    let lin = |a: usize| (0..=n[a]).map(|i| b[a][0] + (b[a][1] - b[a][0]) * i as f64 / n[a] as f64).collect::<Vec<_>>();
    brick_tet(&lin(0), &lin(1), &lin(2))
}

/// A scalar field with exact derivatives, for metric and interpolation-error tests.
pub trait Field<const D: usize> {
    fn value(&self, x: &[f64; D]) -> f64;
    fn gradient(&self, x: &[f64; D]) -> [f64; D];
    fn hessian(&self, x: &[f64; D]) -> Sym<D>;
}

/// `exp(-sum_a ((x_a - c_a) / w_a)^2)`: anisotropic Gaussian bump.
pub struct Gaussian<const D: usize> {
    pub centre: [f64; D],
    pub width: [f64; D],
}

impl<const D: usize> Field<D> for Gaussian<D> {
    fn value(&self, x: &[f64; D]) -> f64 {
        let mut s = 0.0;
        for a in 0..D {
            let t = (x[a] - self.centre[a]) / self.width[a];
            s += t * t;
        }
        (-s).exp()
    }
    fn gradient(&self, x: &[f64; D]) -> [f64; D] {
        let f = self.value(x);
        std::array::from_fn(|a| -2.0 * (x[a] - self.centre[a]) / (self.width[a] * self.width[a]) * f)
    }
    fn hessian(&self, x: &[f64; D]) -> Sym<D> {
        let f = self.value(x);
        let mut h = [[0.0; D]; D];
        for a in 0..D {
            let wa2 = self.width[a] * self.width[a];
            let da = (x[a] - self.centre[a]) / wa2;
            for b in 0..D {
                let wb2 = self.width[b] * self.width[b];
                let db = (x[b] - self.centre[b]) / wb2;
                h[a][b] = 4.0 * da * db * f;
            }
            h[a][a] -= 2.0 / wa2 * f;
        }
        h
    }
}

/// `tanh((x_axis - at) / delta)`: a boundary or interior layer normal to one axis.
pub struct Layer<const D: usize> {
    pub axis: usize,
    pub at: f64,
    pub delta: f64,
}

impl<const D: usize> Field<D> for Layer<D> {
    fn value(&self, x: &[f64; D]) -> f64 {
        ((x[self.axis] - self.at) / self.delta).tanh()
    }
    fn gradient(&self, x: &[f64; D]) -> [f64; D] {
        let t = self.value(x);
        let mut g = [0.0; D];
        g[self.axis] = (1.0 - t * t) / self.delta;
        g
    }
    fn hessian(&self, x: &[f64; D]) -> Sym<D> {
        let t = self.value(x);
        let mut h = [[0.0; D]; D];
        h[self.axis][self.axis] = -2.0 * t * (1.0 - t * t) / (self.delta * self.delta);
        h
    }
}

/// `tanh((n . x - c) / delta)` with unit normal `n`: a layer in an arbitrary direction.
pub struct ObliqueLayer<const D: usize> {
    pub normal: [f64; D],
    pub offset: f64,
    pub delta: f64,
}

impl<const D: usize> ObliqueLayer<D> {
    fn s(&self, x: &[f64; D]) -> f64 {
        let n: f64 = self.normal.iter().map(|v| v * v).sum::<f64>().sqrt();
        let d: f64 = (0..D).map(|a| self.normal[a] * x[a]).sum::<f64>() / n;
        (d - self.offset) / self.delta
    }
    fn unit(&self) -> [f64; D] {
        let n: f64 = self.normal.iter().map(|v| v * v).sum::<f64>().sqrt();
        std::array::from_fn(|a| self.normal[a] / n)
    }
}

impl<const D: usize> Field<D> for ObliqueLayer<D> {
    fn value(&self, x: &[f64; D]) -> f64 {
        self.s(x).tanh()
    }
    fn gradient(&self, x: &[f64; D]) -> [f64; D] {
        let t = self.value(x);
        let n = self.unit();
        std::array::from_fn(|a| (1.0 - t * t) / self.delta * n[a])
    }
    fn hessian(&self, x: &[f64; D]) -> Sym<D> {
        let t = self.value(x);
        let n = self.unit();
        let c = -2.0 * t * (1.0 - t * t) / (self.delta * self.delta);
        std::array::from_fn(|a| std::array::from_fn(|b| c * n[a] * n[b]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_mesh_is_valid_and_labelled() {
        let m = box_tri([[0.0, 2.0], [-1.0, 1.0]], 7, 5);
        m.check().unwrap();
        assert_eq!(m.n_cells(), 70);
        assert!((m.total_measure() - 4.0).abs() < 1e-14);
        assert_eq!(m.n_facets(), 2 * 7 + 2 * 5);
        assert_eq!(m.boundary_keys().len(), m.n_facets());
    }

    #[test]
    fn brick_mesh_is_valid_and_conforming() {
        let m = box_tet([[0.0, 1.0], [0.0, 2.0], [0.0, 3.0]], [3, 2, 2]);
        m.check().unwrap();
        assert_eq!(m.n_cells(), 3 * 2 * 2 * 6);
        assert!((m.total_measure() - 6.0).abs() < 1e-13);
        assert_eq!(m.boundary_keys().len(), m.n_facets());
        assert_eq!(m.n_facets(), 2 * 2 * (3 * 2 + 3 * 2 + 2 * 2));
    }

    #[test]
    fn analytic_derivatives_match_finite_differences() {
        let g = Gaussian::<2> { centre: [0.3, 0.6], width: [0.2, 0.05] };
        let l = Layer::<3> { axis: 1, at: 0.4, delta: 0.1 };
        let x2 = [0.35, 0.58];
        let x3 = [0.1, 0.45, 0.9];
        let h = 1e-5;
        for a in 0..2 {
            let (mut p, mut q) = (x2, x2);
            p[a] += h;
            q[a] -= h;
            let fd = (g.value(&p) - g.value(&q)) / (2.0 * h);
            assert!((fd - g.gradient(&x2)[a]).abs() < 1e-6);
            let fdh = (g.gradient(&p)[a] - g.gradient(&q)[a]) / (2.0 * h);
            assert!((fdh - g.hessian(&x2)[a][a]).abs() < 1e-4);
        }
        let (mut p, mut q) = (x3, x3);
        p[1] += h;
        q[1] -= h;
        assert!(((l.gradient(&p)[1] - l.gradient(&q)[1]) / (2.0 * h) - l.hessian(&x3)[1][1]).abs() < 1e-4);
    }
}
