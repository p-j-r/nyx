//! Shared helpers for the tests of the physics crates and of femcore: small meshes made from geom models with the boundary numbers the
//! tests have always used, and `maxdiff`. Public so that `nyx-dfs`, `nyx-sms` and `nyx-ems` use the same fixtures.
#![allow(dead_code)]

pub fn maxdiff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

use crate::mesh::SimplexMesh;

fn relabel(mesh: &mut SimplexMesh, map: &[(&[i32], i32)]) {
    let old = mesh.facet_labels.clone();
    for (f, l) in mesh.facet_labels.iter_mut().enumerate() {
        *l = map.iter().find(|(from, _)| from.contains(&old[f])).map_or(old[f], |(_, to)| *to);
    }
}

/// A rectangle `[0, lx] x [0, ly]` made from a geom model (`geomdrive::plate`) with element size of the finer of `lx / nx` and `ly / ny`. Its mesh
/// columns are the plane the 2D kernels see; the boundary numbers are the ones these tests have always used: 1 x = 0, 2 y = 0, 3 y = ly,
/// 4 x = lx (the second axis is z of the wafer's x-z plane).
pub fn rect(lx: f64, ly: f64, nx: usize, ny: usize) -> SimplexMesh {
    let p = crate::meshing::geomdrive::plate(lx, ly, (lx / nx as f64).min(ly / ny as f64)).unwrap();
    let mut m = p.mesh;
    relabel(&mut m, &[(&p.xmin, 1), (&p.zmin, 2), (&p.zmax, 3), (&p.xmax, 4)]);
    m
}

/// A box `[0, lx] x [0, ly] x [0, lz]`: the x-z rectangle of `geomdrive::block` swept along y through `ny` layers. Boundary numbers:
/// 1 x = 0, 2 y = 0, 3 z = 0, 4 z = lz, 5 y = ly, 6 x = lx.
pub fn block(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> SimplexMesh {
    let b = crate::meshing::geomdrive::block(lx, ly, lz, (lx / nx as f64).min(lz / nz as f64), ny).unwrap();
    let mut m = b.mesh;
    relabel(&mut m, &[(&b.xmin, 1), (&b.ymin, 2), (&b.zmin, 3), (&b.zmax, 4), (&b.ymax, 5), (&b.xmax, 6)]);
    m
}

fn geom_plate(m: &mut geom::model::Model, parts: &[(&str, (f64, f64), (f64, f64))]) {
    let ids: Vec<(&str, geom::model::FeatureId)> = parts.iter().map(|(n, lo, hi)| (*n, m.rect(n, *lo, *hi))).collect();
    m.overlay(&ids);
}

/// The rectangle `[0, lx] x [0, ly]` in two domains, 1 for `x < xc` and 2 for `x > xc`, meeting along the line `x = xc` (a model of two
/// rectangles, so the interface is a mesh edge). Boundary numbers as `rect`.
pub fn rect_split_x(lx: f64, ly: f64, xc: f64, h: f64) -> SimplexMesh {
    use crate::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
    let mut m = geom::model::Model::new_2d();
    geom_plate(&mut m, &[("a", (0.0, 0.0), (xc, ly)), ("b", (xc, 0.0), (lx, ly))]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h)).unwrap();
    assert_eq!((g.domain("a"), g.domain("b")), (Some(1), Some(2)));
    let (l, b, t, r) = (g.select("plane(x = 0)").unwrap(), g.select("plane(z = 0)").unwrap(), g.select(&format!("plane(z = {ly:e})")).unwrap(), g.select(&format!("plane(x = {lx:e})")).unwrap());
    let mut mesh = g.mesh;
    relabel(&mut mesh, &[(&l, 1), (&b, 2), (&t, 3), (&r, 4)]);
    mesh
}

/// The box `[0, lx] x [0, ly] x [0, lz]` in two domains, 1 for `z < zc` and 2 for `z > zc`: the x-z section of two stacked rectangles swept
/// along y through `ny` layers. Boundary numbers as `block`.
pub fn block_split_z(lx: f64, ly: f64, lz: f64, zc: f64, h: f64, ny: usize) -> SimplexMesh {
    use crate::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
    use crate::meshing::mesher::{extrude_y, Layer};
    let mut m = geom::model::Model::new_2d();
    geom_plate(&mut m, &[("lower", (0.0, 0.0), (lx, zc)), ("upper", (0.0, zc), (lx, lz))]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h)).unwrap();
    assert_eq!((g.domain("lower"), g.domain("upper")), (Some(1), Some(2)));
    let (xmin, xmax, zmin, zmax) = (g.select("plane(x = 0)").unwrap(), g.select(&format!("plane(x = {lx:e})")).unwrap(), g.select("plane(z = 0)").unwrap(), g.select(&format!("plane(z = {lz:e})")).unwrap());
    let nb = g.mesh.facet_labels.iter().copied().max().unwrap();
    let mut mesh = extrude_y(&crate::meshing::geomdrive::sort_nodes(&g.mesh), &[Layer { thickness: ly, n: ny, domain_offset: 0 }], 0.0).unwrap();
    relabel(&mut mesh, &[(&xmin, 1), (&[nb + 1], 2), (&zmin, 3), (&zmax, 4), (&[nb + 2], 5), (&xmax, 6)]);
    mesh
}

/// A stack along y of `a x a` sections (x-z): layer `k` is `layers[k] = (thickness, number of element layers)` and domain `k + 1`.
/// Boundary numbers as `block`; the side walls can be paired by Floquet conditions (1 with 6, 3 with 4) and so can the ends (2 with 5).
pub fn stack_y(a: f64, layers: &[(f64, usize)], h: f64) -> SimplexMesh {
    use crate::meshing::mesher::{extrude_y, Layer};
    let p = crate::meshing::geomdrive::plate(a, a, h).unwrap();
    let nb = p.mesh.facet_labels.iter().copied().max().unwrap();
    let ls: Vec<Layer> = layers.iter().enumerate().map(|(k, (t, n))| Layer { thickness: *t, n: *n, domain_offset: k as i32 }).collect();
    let mut mesh = extrude_y(&crate::meshing::geomdrive::sort_nodes(&p.mesh), &ls, 0.0).unwrap();
    relabel(&mut mesh, &[(&p.xmin, 1), (&[nb + 1], 2), (&p.zmin, 3), (&p.zmax, 4), (&[nb + 2], 5), (&p.xmax, 6)]);
    mesh
}

/// One axis-aligned rectangle of a section: `(x0, z0, x1, z1)`, the domain number it gets, and an element size of its own (else the
/// global one).
pub struct RectPart {
    pub lo: (f64, f64),
    pub hi: (f64, f64),
    pub domain: i32,
    pub h: Option<f64>,
}

pub fn part(x0: f64, z0: f64, x1: f64, z1: f64, domain: i32, h: Option<f64>) -> RectPart {
    RectPart { lo: (x0, z0), hi: (x1, z1), domain, h }
}

/// The rectangles laid over each other in order (later ones win), as a geom model meshed with element size `h`; the cells carry the
/// domain numbers asked for. The boundary numbers are those of the model's edges.
pub fn parts_mesh(parts: &[RectPart], h: f64) -> SimplexMesh {
    use crate::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
    let mut m = geom::model::Model::new_2d();
    let names: Vec<String> = (0..parts.len()).map(|k| format!("part{k}")).collect();
    let ids: Vec<(&str, geom::model::FeatureId)> = parts.iter().zip(&names).map(|(p, n)| (n.as_str(), m.rect(n, p.lo, p.hi))).collect();
    m.overlay(&ids);
    for (p, n) in parts.iter().zip(&names) {
        if let Some(hp) = p.h {
            m.size(&format!("size_{n}"), geom::select::Selector::parse(&format!("region({n})")).unwrap(), Some(hp.into()), None, None);
        }
    }
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h)).unwrap();
    let mut mesh = g.mesh.clone();
    let domain_of: std::collections::BTreeMap<i32, i32> = parts.iter().zip(&names).filter_map(|(p, n)| g.domain(n).map(|d| (d, p.domain))).collect();
    for l in mesh.cell_labels.iter_mut() {
        *l = domain_of[l];
    }
    mesh
}
