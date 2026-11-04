//! The edge-element basis of the 2D mode solvers and of everything that evaluates their fields: Whitney and degree-2 Nedelec functions on a
//! triangle (`edge_basis`), the per-cell numbering and signs of the edge functions (`cell_data`) and the barycentric gradients of a triangle
//! (`tri_g`). Split out of `ems` when the electromagnetic kernels moved into `nyx-ems` (ADR 45): `optical` (the field sampler of the
//! polarization analysis) and the kernels both use them, so they stay below both in `femcore`.

use crate::fem::Geometry;
use crate::mesh::SimplexMesh;

const PAIRS: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

fn cross(u: [f64; 2], v: [f64; 2]) -> f64 {
    u[0] * v[1] - u[1] * v[0]
}

/// Nedelec basis at barycentric point `l` of one cell with barycentric gradients `g[i]`:
/// `w[s]` and `curl[s]` for the local functions (3 Whitney signed to the global edge
/// direction, then for order 2 the gradients of the three edge bubbles and two interior
/// functions).
pub fn edge_basis(l: &[f64], g: &[[f64; 2]; 3], sign: &[f64; 3], order: usize, w: &mut [[f64; 2]], curl: &mut [f64]) {
    for (s, &(a, b)) in PAIRS.iter().enumerate() {
        for k in 0..2 {
            w[s][k] = (l[a] * g[b][k] - l[b] * g[a][k]) * sign[s];
        }
        curl[s] = 2.0 * cross(g[a], g[b]) * sign[s];
    }
    if order == 2 {
        for (s, &(a, b)) in PAIRS.iter().enumerate() {
            for k in 0..2 {
                w[3 + s][k] = l[a] * g[b][k] + l[b] * g[a][k];
            }
            curl[3 + s] = 0.0;
        }
        let wab = |a: usize, b: usize| [l[a] * g[b][0] - l[b] * g[a][0], l[a] * g[b][1] - l[b] * g[a][1]];
        let cab = |a: usize, b: usize| 2.0 * cross(g[a], g[b]);
        let w01 = wab(0, 1);
        let w12 = wab(1, 2);
        for k in 0..2 {
            w[6][k] = l[2] * w01[k];
            w[7][k] = l[0] * w12[k];
        }
        curl[6] = cross(g[2], w01) + l[2] * cab(0, 1);
        curl[7] = cross(g[0], w12) + l[0] * cab(1, 2);
    }
}

pub struct CellData {
    pub sign: Vec<[f64; 3]>,
    pub gdofs: Vec<Vec<usize>>,
    pub net: usize,
}

pub fn cell_data(mesh: &SimplexMesh, order: usize) -> CellData {
    let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
    let mut sign = vec![[1.0; 3]; nc];
    let mut gdofs = Vec::with_capacity(nc);
    for c in 0..nc {
        let v = mesh.cell(c);
        let eid = mesh.cell_edge_ids(c);
        for (s, &(a, b)) in PAIRS.iter().enumerate() {
            sign[c][s] = if v[a] < v[b] { 1.0 } else { -1.0 };
        }
        let mut g: Vec<usize> = eid.iter().map(|&e| e as usize).collect();
        if order == 2 {
            g.extend(eid.iter().map(|&e| ne + e as usize));
            g.push(2 * ne + 2 * c);
            g.push(2 * ne + 2 * c + 1);
        }
        gdofs.push(g);
    }
    let net = if order == 1 { ne } else { 2 * ne + 2 * nc };
    CellData { sign, gdofs, net }
}

pub fn tri_g(geo: &Geometry, c: usize) -> [[f64; 2]; 3] {
    let b = c * 6;
    [[geo.g[b], geo.g[b + 1]], [geo.g[b + 2], geo.g[b + 3]], [geo.g[b + 4], geo.g[b + 5]]]
}

