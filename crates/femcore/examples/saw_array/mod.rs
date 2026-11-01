//! Array helpers shared by the SAW examples that need arrays of ridges (copied from `saw_s3.rs`): ridge arrays on sapphire as meshes,
//! the problem set-up, per-vertex mass weights, and a point locator for P1 interpolation of the single-ridge mode shape.
#![allow(dead_code)]
use crate::saw_common::*;
use femcore::csparse::C;
use femcore::mesh::SimplexMesh;
use femcore::mesher::*;
use femcore::problem::*;
use femcore::solid_kz::KzModes;

pub const CORE_MIN: f64 = 0.5;
pub const PML_MAX: f64 = 0.2;
pub const PML_DOMAINS: [i32; 5] = [4, 5, 6, 8, 9];

pub struct Array {
    pub lam: f64,
    /// centres of the N ridges (full array), left to right
    pub xc: Vec<f64>,
    pub xe: f64,
    pub xo: f64,
    pub d: f64,
    pub pt: f64,
}

pub fn centres(lam: f64, gaps: &[f64]) -> Vec<f64> {
    let mut x = vec![0.0];
    for g in gaps {
        x.push(x.last().unwrap() + lam + g);
    }
    let mid = 0.5 * x.last().unwrap();
    x.iter().map(|v| v - mid).collect()
}

pub fn build_array(lam: f64, gaps: &[f64], half: bool, scale: f64, side: f64, depth_lam: f64) -> (SimplexMesh, Array) {
    let xc = centres(lam, gaps);
    let xe = xc.last().unwrap() + 0.5 * lam;
    let (xin, pt, d) = (xe + side * lam, lam, depth_lam * lam);
    let xo = xin + pt;
    let (h_far, h_surf, h_si, h_edge) = (lam / 8.0 * scale, lam / 16.0 * scale, 75e-9 * scale, 75e-9 * scale);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: dom, h };
    let xlo = if half { 0.0 } else { -xin };
    let mut parts = vec![rect(xlo, -d, xin, 0.0, 1, None), rect(xlo, -0.6 * lam, xin, 0.0, 1, Some(h_surf))];
    for &x in xc.iter().filter(|&&x| !half || x > -1e-12) {
        let (a, b) = ((x - 0.5 * lam - 0.25 * UM).max(xlo), x + 0.5 * lam + 0.25 * UM);
        parts.push(rect(a, -0.25 * UM, b, 0.0, 1, Some(h_edge)));
    }
    parts.push(rect(xin, -d, xo, 0.0, 5, None));
    parts.push(rect(xlo, -d - pt, xin, -d, 6, None));
    parts.push(rect(xin, -d - pt, xo, -d, 9, None));
    if !half {
        parts.push(rect(-xo, -d, -xin, 0.0, 4, None));
        parts.push(rect(-xo, -d - pt, -xin, -d, 8, None));
    }
    for &x in xc.iter().filter(|&&x| !half || x > -1e-12) {
        parts.push(rect((x - 0.5 * lam).max(xlo), 0.0, x + 0.5 * lam, T_SI, 2, Some(h_si)));
    }
    let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
    (mesh, Array { lam, xc, xe, xo, d, pt })
}

pub fn base_spec(psi: f64, half: bool) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    set_materials(&mut s, psi, &[1, 2, 4, 5, 6, 8, 9], &[1, 4, 5, 6, 8, 9], &[2]);
    let mut dirs = vec![(5, PmlDir::Plus, PmlDir::Off), (6, PmlDir::Off, PmlDir::Minus), (9, PmlDir::Plus, PmlDir::Minus)];
    if !half {
        dirs.extend([(4, PmlDir::Minus, PmlDir::Off), (8, PmlDir::Minus, PmlDir::Minus)]);
    }
    for (dom, px, py) in dirs {
        s.materials.pml.insert(dom, pml(px, py));
    }
    s
}

pub fn cell_x(mesh: &SimplexMesh) -> Vec<f64> {
    (0..mesh.n_cells()).map(|c| mesh.cell(c).iter().map(|&v| mesh.point(v as usize)[0]).sum::<f64>() / 3.0).collect()
}

/// Lumped mass weight of every vertex (rho x area / 3 of the cells around it).
pub fn vertex_weights(mesh: &SimplexMesh) -> Vec<f64> {
    let mut w = vec![0.0; mesh.n_points()];
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let (a, b, d) = (mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize));
        let area = 0.5 * ((b[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (b[1] - a[1])).abs();
        let rho = if mesh.cell_labels[c] == 2 { RHO_SI } else { RHO_SAPPHIRE };
        for &i in v {
            w[i as usize] += rho * area / 3.0;
        }
    }
    w
}

/// P1 interpolation of vertex data on a mesh, with a grid hash for point location.
pub struct Locator<'a> {
    mesh: &'a SimplexMesh,
    h: f64,
    bins: std::collections::HashMap<(i64, i64), Vec<u32>>,
}

impl<'a> Locator<'a> {
    pub fn new(mesh: &'a SimplexMesh) -> Self {
        let h = 0.5 * UM;
        let mut bins: std::collections::HashMap<(i64, i64), Vec<u32>> = Default::default();
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let p: Vec<&[f64]> = v.iter().map(|&i| mesh.point(i as usize)).collect();
            let (x0, x1) = (p.iter().map(|q| q[0]).fold(f64::MAX, f64::min), p.iter().map(|q| q[0]).fold(f64::MIN, f64::max));
            let (y0, y1) = (p.iter().map(|q| q[1]).fold(f64::MAX, f64::min), p.iter().map(|q| q[1]).fold(f64::MIN, f64::max));
            for i in (x0 / h).floor() as i64..=(x1 / h).floor() as i64 {
                for j in (y0 / h).floor() as i64..=(y1 / h).floor() as i64 {
                    bins.entry((i, j)).or_default().push(c as u32);
                }
            }
        }
        Locator { mesh, h, bins }
    }

    /// (cell, barycentric coordinates) of the point (x, y) if it lies in the mesh.
    pub fn locate(&self, x: f64, y: f64) -> Option<(usize, [f64; 3])> {
        let key = ((x / self.h).floor() as i64, (y / self.h).floor() as i64);
        for &c in self.bins.get(&key)? {
            let v = self.mesh.cell(c as usize);
            let (a, b, d) = (self.mesh.point(v[0] as usize), self.mesh.point(v[1] as usize), self.mesh.point(v[2] as usize));
            let det = (b[1] - d[1]) * (a[0] - d[0]) + (d[0] - b[0]) * (a[1] - d[1]);
            let l0 = ((b[1] - d[1]) * (x - d[0]) + (d[0] - b[0]) * (y - d[1])) / det;
            let l1 = ((d[1] - a[1]) * (x - d[0]) + (a[0] - d[0]) * (y - d[1])) / det;
            let l2 = 1.0 - l0 - l1;
            let eps = -1e-9;
            if l0 >= eps && l1 >= eps && l2 >= eps {
                return Some((c as usize, [l0, l1, l2]));
            }
        }
        None
    }

    pub fn eval(&self, data: &[[C; 3]], x: f64, y: f64) -> [C; 3] {
        match self.locate(x, y) {
            Some((c, l)) => {
                let v = self.mesh.cell(c);
                let mut out = [C::new(0.0, 0.0); 3];
                for k in 0..3 {
                    for a in 0..3 {
                        out[k] += data[v[a] as usize][k] * l[a];
                    }
                }
                out
            }
            None => [C::new(0.0, 0.0); 3],
        }
    }
}

pub fn vertex_vectors(mesh: &SimplexMesh, m: &KzModes, j: usize) -> Vec<[C; 3]> {
    let n1 = mesh.n_dofs(2);
    let u = &m.vectors[j];
    (0..mesh.n_points()).map(|i| [u[i], u[n1 + i], u[2 * n1 + i]]).collect()
}

