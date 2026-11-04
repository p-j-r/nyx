//! Element-quality and edge-length statistics, Euclidean or measured in a metric field.
//!
//! Simplex quality is the mean ratio (1 for the regular simplex, 0 for a degenerate one):
//! `4 sqrt(3) A / sum l^2` in 2D and `12 (3 V)^(2/3) / sum l^2` in 3D. In a metric the lengths and
//! the measure are taken in the element-average metric. Quads, prisms and hexahedra are rated by
//! their smallest corner scaled Jacobian. These are computed straight from the coordinates so they
//! can serve as an independent check of the code that produced the mesh.

use crate::adapt::metric::{Sym, det, identity, length};
use crate::mesh::cells::{CellMesh, CellType};
use crate::adapt::work::{SimplexMesh, simplex_signed_measure};

/// Mean-ratio quality of a simplex with vertices `p` (`D + 1` points), optionally in the metric
/// obtained by averaging the vertex metrics. Negative for inverted simplices.
pub fn simplex_quality<const D: usize>(p: &[[f64; D]], metrics: Option<&[Sym<D>]>) -> f64 {
    let avg: Sym<D> = match metrics {
        None => identity::<D>(),
        Some(ms) => {
            let mut a = [[0.0; D]; D];
            for m in ms {
                for i in 0..D {
                    for j in 0..D {
                        a[i][j] += m[i][j] / ms.len() as f64;
                    }
                }
            }
            a
        }
    };
    let mut sum_l2 = 0.0;
    for i in 0..=D {
        for j in i + 1..=D {
            let mut e = [0.0; D];
            for a in 0..D {
                e[a] = p[j][a] - p[i][a];
            }
            let l = length(&avg, &e);
            sum_l2 += l * l;
        }
    }
    let meas = simplex_signed_measure::<D>(p) * det(&avg).sqrt();
    match D {
        2 => 4.0 * 3f64.sqrt() * meas / sum_l2,
        3 => {
            let s = 3.0 * meas;
            12.0 * s.signum() * s.abs().powf(2.0 / 3.0) / sum_l2
        }
        _ => panic!("only D = 2 and D = 3 are supported"),
    }
}

/// Summary of a set of scalar samples in `[0, 1]` (quality) with a 10-bin histogram.
#[derive(Clone, Debug, Default)]
pub struct QualityStats {
    pub n: usize,
    pub min: f64,
    pub mean: f64,
    /// Counts of samples in `[0, 0.1), [0.1, 0.2), ... [0.9, 1]`; negative samples go to bin 0.
    pub histogram: [usize; 10],
}

impl QualityStats {
    fn from_samples(it: impl Iterator<Item = f64>) -> Self {
        let mut s = QualityStats { min: f64::INFINITY, ..Default::default() };
        let mut sum = 0.0;
        for q in it {
            s.n += 1;
            s.min = s.min.min(q);
            sum += q;
            let b = ((q.max(0.0) * 10.0) as usize).min(9);
            s.histogram[b] += 1;
        }
        s.mean = if s.n > 0 { sum / s.n as f64 } else { 0.0 };
        s
    }
}

/// Mean-ratio statistics of a simplex mesh, with an optional vertex metric field.
pub fn mesh_quality<const D: usize>(mesh: &SimplexMesh<D>, metrics: Option<&[Sym<D>]>) -> QualityStats {
    QualityStats::from_samples((0..mesh.n_cells()).map(|c| {
        let vs = mesh.cell(c);
        let mut p = [[0.0; D]; 4];
        for (i, &v) in vs.iter().enumerate() {
            p[i] = mesh.points[v as usize];
        }
        match metrics {
            None => simplex_quality::<D>(&p[..D + 1], None),
            Some(ms) => {
                let mut mm = [[[0.0; D]; D]; 4];
                for (i, &v) in vs.iter().enumerate() {
                    mm[i] = ms[v as usize];
                }
                simplex_quality::<D>(&p[..D + 1], Some(&mm[..D + 1]))
            }
        }
    }))
}

/// Statistics of edge lengths measured in a metric (or Euclidean when `None`).
#[derive(Clone, Debug, Default)]
pub struct EdgeStats {
    pub n: usize,
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    /// Fraction of edges with length in `[1/sqrt(2), sqrt(2)]`, the usual "conformity to the
    /// metric" figure of merit for adapted meshes.
    pub in_band: f64,
    /// Counts in `[0, 0.3), [0.3, 0.6), [0.6, 0.7071), [0.7071, 0.9), [0.9, 1.3), [1.3, 1.4142),
    /// [1.4142, 2), [2, 5), [5, inf)`.
    pub histogram: [usize; 9],
}

/// Edge-length statistics of `mesh`; lengths use [`crate::adapt::metric::edge_length`].
pub fn edge_stats<const D: usize>(mesh: &SimplexMesh<D>, metrics: Option<&[Sym<D>]>) -> EdgeStats {
    use crate::adapt::metric::edge_length;
    let id = identity::<D>();
    let mut s = EdgeStats { min: f64::INFINITY, ..Default::default() };
    let bounds = [0.3, 0.6, std::f64::consts::FRAC_1_SQRT_2, 0.9, 1.3, std::f64::consts::SQRT_2, 2.0, 5.0];
    let (mut sum, mut band) = (0.0, 0usize);
    for e in mesh.edges() {
        let (a, b) = (e[0] as usize, e[1] as usize);
        let mut v = [0.0; D];
        for k in 0..D {
            v[k] = mesh.points[b][k] - mesh.points[a][k];
        }
        let l = match metrics {
            Some(ms) => edge_length(&ms[a], &ms[b], &v),
            None => length(&id, &v),
        };
        s.n += 1;
        s.min = s.min.min(l);
        s.max = s.max.max(l);
        sum += l;
        if (std::f64::consts::FRAC_1_SQRT_2..=std::f64::consts::SQRT_2).contains(&l) {
            band += 1;
        }
        s.histogram[bounds.iter().position(|&x| l < x).unwrap_or(8)] += 1;
    }
    if s.n > 0 {
        s.mean = sum / s.n as f64;
        s.in_band = band as f64 / s.n as f64;
    }
    s
}

/// Neighbour triples (right-handed on the reference shape) at every corner of a prism or hex.
fn corner_triples(t: CellType) -> &'static [[u8; 3]] {
    use std::sync::OnceLock;
    static HEX: OnceLock<Vec<[u8; 3]>> = OnceLock::new();
    static PRISM: OnceLock<Vec<[u8; 3]>> = OnceLock::new();
    let build = |nbrs: &[[u8; 3]], pts: &[[f64; 3]]| -> Vec<[u8; 3]> {
        nbrs.iter()
            .enumerate()
            .map(|(i, n)| {
                let v = |k: u8| {
                    let q = pts[k as usize];
                    [q[0] - pts[i][0], q[1] - pts[i][1], q[2] - pts[i][2]]
                };
                let (a, b, c) = (v(n[0]), v(n[1]), v(n[2]));
                let d = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0]);
                if d > 0.0 { *n } else { [n[0], n[2], n[1]] }
            })
            .collect()
    };
    match t {
        CellType::Hex => HEX.get_or_init(|| {
            let pts = [
                [0., 0., 0.],
                [1., 0., 0.],
                [1., 1., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [1., 0., 1.],
                [1., 1., 1.],
                [0., 1., 1.],
            ];
            build(&[[1, 3, 4], [0, 2, 5], [1, 3, 6], [0, 2, 7], [0, 5, 7], [1, 4, 6], [2, 5, 7], [3, 4, 6]], &pts)
        }),
        CellType::Prism => PRISM.get_or_init(|| {
            let pts = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]];
            build(&[[1, 2, 3], [0, 2, 4], [0, 1, 5], [0, 4, 5], [1, 3, 5], [2, 3, 4]], &pts)
        }),
        _ => &[],
    }
}

/// Smallest corner scaled Jacobian of a mixed cell (1 for a perfect square, cube or right prism
/// corner; at most 1; negative when inverted). Triangles and tetrahedra return the mean ratio.
pub fn cell_scaled_jacobian_min<const D: usize>(m: &CellMesh<D>, c: usize) -> f64 {
    let vs = m.cell(c);
    let pt = |i: usize| m.points[vs[i] as usize];
    match m.types[c] {
        CellType::Tri => simplex_quality::<D>(&[pt(0), pt(1), pt(2)], None),
        CellType::Tet => simplex_quality::<D>(&[pt(0), pt(1), pt(2), pt(3)], None),
        CellType::Quad => {
            let mut q = f64::INFINITY;
            for i in 0..4 {
                let (a, b, c0) = (pt((i + 3) % 4), pt(i), pt((i + 1) % 4));
                let u = [a[0] - b[0], a[1] - b[1]];
                let v = [c0[0] - b[0], c0[1] - b[1]];
                // interior angle at b runs from (next - b) to (prev - b) counter-clockwise
                let cross = v[0] * u[1] - v[1] * u[0];
                q = q.min(cross / ((u[0].hypot(u[1])) * (v[0].hypot(v[1]))));
            }
            q
        }
        CellType::Prism | CellType::Hex => {
            let mut q = f64::INFINITY;
            for (i, n) in corner_triples(m.types[c]).iter().enumerate() {
                let e = |k: u8| {
                    let p = pt(k as usize);
                    let o = pt(i);
                    [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
                };
                let (a, b, cc) = (e(n[0]), e(n[1]), e(n[2]));
                let d = a[0] * (b[1] * cc[2] - b[2] * cc[1]) - a[1] * (b[0] * cc[2] - b[2] * cc[0])
                    + a[2] * (b[0] * cc[1] - b[1] * cc[0]);
                let nrm = |x: [f64; 3]| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
                q = q.min(d / (nrm(a) * nrm(b) * nrm(cc)));
            }
            q
        }
    }
}

/// Quality statistics over all cells of a mixed mesh (scaled Jacobian / mean ratio).
pub fn cellmesh_quality<const D: usize>(m: &CellMesh<D>) -> QualityStats {
    QualityStats::from_samples((0..m.n_cells()).map(|c| cell_scaled_jacobian_min(m, c)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_simplices_have_quality_one() {
        let h = 3f64.sqrt() / 2.0;
        let t = [[0.0, 0.0], [1.0, 0.0], [0.5, h]];
        assert!((simplex_quality::<2>(&t, None) - 1.0).abs() < 1e-14);
        let tet = [[1.0, 1.0, 1.0], [-1.0, 1.0, -1.0], [1.0, -1.0, -1.0], [-1.0, -1.0, 1.0]];
        assert!((simplex_quality::<3>(&tet, None) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn quality_is_scale_invariant_and_negative_when_inverted() {
        let t = [[0.0, 0.0], [3.0, 0.0], [1.0, 2.0]];
        let s = [[0.0, 0.0], [30.0, 0.0], [10.0, 20.0]];
        assert!((simplex_quality::<2>(&t, None) - simplex_quality::<2>(&s, None)).abs() < 1e-14);
        let inv = [t[0], t[2], t[1]];
        assert!(simplex_quality::<2>(&inv, None) < 0.0);
    }

    #[test]
    fn metric_quality_of_stretched_triangle_is_one_in_its_own_metric() {
        // right isosceles-like triangle stretched 10x in x; metric undoes the stretch
        let h = 3f64.sqrt() / 2.0;
        let t = [[0.0, 0.0], [10.0, 0.0], [5.0, h]];
        let m = crate::adapt::metric::from_sizes(&[10.0, 1.0], &identity::<2>());
        let ms = [m, m, m];
        assert!((simplex_quality::<2>(&t, Some(&ms)) - 1.0).abs() < 1e-13);
        assert!(simplex_quality::<2>(&t, None) < 0.5);
    }

    #[test]
    fn jacobians_of_perfect_shapes() {
        let mut m = CellMesh::<3>::new();
        m.points = vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 1.],
            [1., 1., 1.],
            [0., 1., 1.],
        ];
        m.push_cell(CellType::Hex, &[0, 1, 2, 3, 4, 5, 6, 7], 1);
        assert!((cell_scaled_jacobian_min(&m, 0) - 1.0).abs() < 1e-14);
        let mut p = CellMesh::<3>::new();
        p.points = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]];
        p.push_cell(CellType::Prism, &[0, 1, 2, 3, 4, 5], 1);
        let q = cell_scaled_jacobian_min(&p, 0);
        assert!(q > 0.7 && q <= 1.0, "{q}");
        let mut q2 = CellMesh::<2>::new();
        q2.points = vec![[0., 0.], [2., 0.], [2., 1.], [0., 1.]];
        q2.push_cell(CellType::Quad, &[0, 1, 2, 3], 1);
        assert!((cell_scaled_jacobian_min(&q2, 0) - 1.0).abs() < 1e-14);
        // a flipped quad is negative
        let mut q3 = CellMesh::<2>::new();
        q3.points = q2.points.clone();
        q3.push_cell(CellType::Quad, &[0, 3, 2, 1], 1);
        assert!(cell_scaled_jacobian_min(&q3, 0) < 0.0);
    }
}
