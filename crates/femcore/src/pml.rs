//! Perfectly matched layers by complex coordinate stretching, shared by the electromagnetic kernels (`nyx-ems`) and the solid kernels: the
//! stretch factor of a point inside a layer (`stretch`, from a `props::Pml` and the bounding box of its domain) and the stretch of every
//! cell of a 2D or 3D mesh (`cell_stretches`). Split out of `ems_pml` when `ems` moved into its crate (ADR 45).

use crate::linalg::csr::{c, C};
use crate::mesh::SimplexMesh;
use crate::props::{Pml, PmlDir};

/// Stretch factors along the two axes of the mesh (`dirs`, the layer directions of each) of a cell centre for a PML specification and the
/// bounding box of its domain.
pub fn stretch(p: &Pml, dirs: [PmlDir; 2], centre: [f64; 2], lo: [f64; 2], hi: [f64; 2]) -> (C, C) {
    let one = |dir: PmlDir, x: f64, lo: f64, hi: f64| -> C {
        let t = match dir {
            PmlDir::Off => return c(1.0, 0.0),
            PmlDir::Plus => (x - lo) / (hi - lo).max(1e-300),
            PmlDir::Minus => (hi - x) / (hi - lo).max(1e-300),
        };
        let t = t.clamp(0.0, 1.0);
        let prof = if p.order == 0.0 { 1.0 } else { t.powf(p.order) };
        c(p.re, -p.strength * prof)
    };
    (one(dirs[0], centre[0], lo[0], hi[0]), one(dirs[1], centre[1], lo[1], hi[1]))
}

/// Complex stretch factors `(sx, sy, sz)` of every cell of a 2D or 3D mesh for the PML domains of
/// `pml` (cells of other domains get 1). The profile is evaluated at the cell centre against the
/// bounding box of its domain.
pub fn cell_stretches(mesh: &SimplexMesh, pml: &std::collections::BTreeMap<i32, Pml>) -> Vec<[C; 3]> {
    cell_stretches_by(mesh, pml, |p, d| crate::section::mesh_pml_dirs(p, d))
}

fn cell_stretches_by(mesh: &SimplexMesh, pml: &std::collections::BTreeMap<i32, Pml>, dirs: impl Fn(&Pml, usize) -> [PmlDir; 3]) -> Vec<[C; 3]> {
    let d = mesh.dim;
    let mut boxes: std::collections::BTreeMap<i32, ([f64; 3], [f64; 3])> = Default::default();
    for ci in 0..mesh.n_cells() {
        let l = mesh.cell_labels[ci];
        if pml.get(&l).is_some_and(|p| p.x != PmlDir::Off || p.y != PmlDir::Off || p.z != PmlDir::Off) {
            let e = boxes.entry(l).or_insert(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]));
            for &v in mesh.cell(ci) {
                let pt = mesh.point(v as usize);
                for k in 0..d {
                    e.0[k] = e.0[k].min(pt[k]);
                    e.1[k] = e.1[k].max(pt[k]);
                }
            }
        }
    }
    (0..mesh.n_cells())
        .map(|ci| {
            let mut s = [c(1.0, 0.0); 3];
            if let (Some(p), Some((lo, hi))) = (pml.get(&mesh.cell_labels[ci]), boxes.get(&mesh.cell_labels[ci])) {
                let mut cen = [0.0; 3];
                for &v in mesh.cell(ci) {
                    for k in 0..d {
                        cen[k] += mesh.point(v as usize)[k] / (d as f64 + 1.0);
                    }
                }
                for (k, dir) in dirs(p, d).into_iter().enumerate().take(d) {
                    s[k] = stretch_one(p, dir, cen[k], lo[k], hi[k]);
                }
            }
            s
        })
        .collect()
}

fn stretch_one(p: &Pml, dir: PmlDir, x: f64, lo: f64, hi: f64) -> C {
    let t = match dir {
        PmlDir::Off => return c(1.0, 0.0),
        PmlDir::Plus => (x - lo) / (hi - lo).max(1e-300),
        PmlDir::Minus => (hi - x) / (hi - lo).max(1e-300),
    };
    let t = t.clamp(0.0, 1.0);
    let prof = if p.order == 0.0 { 1.0 } else { t.powf(p.order) };
    c(p.re, -p.strength * prof)
}
