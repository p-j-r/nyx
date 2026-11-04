//! Meshes with quadrilaterals, hexahedra and prisms are split into simplices; the result must be
//! conforming, fill the same volume and solve exactly for linear fields.
mod common;
use femcore::mesh::{SimplexMesh, SplitMode};
use femcore::problem::*;
use femcore::sms::solve_solid;
use femcore::mesh::cells::{CellMesh as NativeCells, CellType};
use std::collections::HashMap;

// The fixtures below are written in tensor-product vertex order (i + 2 j + 4 k), the way structured generators number the corners of
// a cell. `native` turns them into the counter-clockwise order of `mesh::cells::CellMesh`.
#[derive(Clone, Debug, Default)]
struct ElementType {
    name: String,
    nep: usize,
    elem: Vec<u32>,
    dom: Vec<i32>,
}

impl ElementType {
    fn len(&self) -> usize {
        self.elem.len().checked_div(self.nep).unwrap_or(0)
    }
    fn row(&self, i: usize) -> &[u32] {
        &self.elem[i * self.nep..(i + 1) * self.nep]
    }
}

#[derive(Clone, Debug, Default)]
struct CellMesh {
    sdim: usize,
    points: Vec<f64>,
    types: Vec<ElementType>,
}

impl CellMesh {
    fn n_points(&self) -> usize {
        self.points.len().checked_div(self.sdim).unwrap_or(0)
    }
    fn ty(&self, name: &str) -> Option<&ElementType> {
        self.types.iter().find(|t| t.name == name)
    }
}

fn native<const D: usize>(nm: &CellMesh) -> NativeCells<D> {
    assert_eq!(nm.sdim, D);
    let mut out = NativeCells::<D>::new();
    out.points = nm.points.chunks(D).map(|p| std::array::from_fn(|k| p[k])).collect();
    for t in &nm.types {
        for i in 0..t.len() {
            let r = t.row(i);
            let lab = t.dom[i];
            match (D, t.name.as_str()) {
                (2, "tri") => out.push_cell(CellType::Tri, r, lab),
                (2, "quad") => out.push_cell(CellType::Quad, &[r[0], r[1], r[3], r[2]], lab),
                (2, "edg") | (3, "tri") => out.push_facet(r, lab),
                (3, "quad") => out.push_facet(&[r[0], r[1], r[3], r[2]], lab),
                (3, "tet") => out.push_cell(CellType::Tet, r, lab),
                (3, "hex") => out.push_cell(CellType::Hex, &[r[0], r[1], r[3], r[2], r[4], r[5], r[7], r[6]], lab),
                (3, "prism") => out.push_cell(CellType::Prism, r, lab),
                _ => panic!("unexpected element type {}", t.name),
            }
        }
    }
    out
}

fn grid_points(nx: usize, ny: usize, nz: usize, lx: f64, ly: f64, lz: f64) -> Vec<f64> {
    let mut p = vec![];
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                p.extend_from_slice(&[lx * i as f64 / nx as f64, ly * j as f64 / ny as f64, lz * k as f64 / nz as f64]);
            }
        }
    }
    p
}

fn et(name: &str, nep: usize, elem: Vec<u32>, dom: Vec<i32>) -> ElementType {
    ElementType { name: name.into(), nep, elem, dom }
}

/// Hexahedral block in the standard vertex order (index = i + 2j + 4k); the six faces are the
/// boundary quads 0..5 (x=0, x=L, y=0, y=L, z=0, z=L).
fn hex_block(nx: usize, ny: usize, nz: usize, lx: f64, ly: f64, lz: f64) -> CellMesh {
    let vid = |i: usize, j: usize, k: usize| ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
    let mut hex = vec![];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                for (dk, dj, di) in [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)] {
                    hex.push(vid(i + di, j + dj, k + dk));
                }
            }
        }
    }
    let nh = nx * ny * nz;
    let (mut quad, mut qdom) = (vec![], vec![]);
    let mut q = |a: u32, b: u32, c: u32, d: u32, lab: i32| {
        quad.extend_from_slice(&[a, b, c, d]);
        qdom.push(lab);
    };
    for k in 0..nz {
        for j in 0..ny {
            q(vid(0, j, k), vid(0, j + 1, k), vid(0, j, k + 1), vid(0, j + 1, k + 1), 1);
            q(vid(nx, j, k), vid(nx, j + 1, k), vid(nx, j, k + 1), vid(nx, j + 1, k + 1), 2);
        }
        for i in 0..nx {
            q(vid(i, 0, k), vid(i + 1, 0, k), vid(i, 0, k + 1), vid(i + 1, 0, k + 1), 3);
            q(vid(i, ny, k), vid(i + 1, ny, k), vid(i, ny, k + 1), vid(i + 1, ny, k + 1), 4);
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            q(vid(i, j, 0), vid(i + 1, j, 0), vid(i, j + 1, 0), vid(i + 1, j + 1, 0), 5);
            q(vid(i, j, nz), vid(i + 1, j, nz), vid(i, j + 1, nz), vid(i + 1, j + 1, nz), 6);
        }
    }
    CellMesh { sdim: 3, points: grid_points(nx, ny, nz, lx, ly, lz), types: vec![et("quad", 4, quad, qdom), et("hex", 8, hex, vec![1; nh])] }
}

/// A block of prisms: each (i, j) square is cut into two triangles that are extruded in z.
fn prism_block(nx: usize, ny: usize, nz: usize, lx: f64, ly: f64, lz: f64) -> CellMesh {
    let vid = |i: usize, j: usize, k: usize| ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
    let (mut prism, mut pd) = (vec![], vec![]);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                for tri in [[(i, j), (i + 1, j), (i + 1, j + 1)], [(i, j), (i + 1, j + 1), (i, j + 1)]] {
                    for kk in [k, k + 1] {
                        for &(a, b) in &tri {
                            prism.push(vid(a, b, kk));
                        }
                    }
                    pd.push(1);
                }
            }
        }
    }
    // boundary: triangles on z = 0 and z = L (labels 4, 5), quads on the four sides (labels 0..3)
    let (mut tri, mut td) = (vec![], vec![]);
    for j in 0..ny {
        for i in 0..nx {
            for t in [[(i, j), (i + 1, j), (i + 1, j + 1)], [(i, j), (i + 1, j + 1), (i, j + 1)]] {
                tri.extend(t.iter().map(|&(a, b)| vid(a, b, 0)));
                td.push(5);
                tri.extend(t.iter().map(|&(a, b)| vid(a, b, nz)));
                td.push(6);
            }
        }
    }
    let (mut quad, mut qd) = (vec![], vec![]);
    for k in 0..nz {
        for j in 0..ny {
            quad.extend_from_slice(&[vid(0, j, k), vid(0, j + 1, k), vid(0, j, k + 1), vid(0, j + 1, k + 1)]);
            qd.push(1);
            quad.extend_from_slice(&[vid(nx, j, k), vid(nx, j + 1, k), vid(nx, j, k + 1), vid(nx, j + 1, k + 1)]);
            qd.push(2);
        }
        for i in 0..nx {
            quad.extend_from_slice(&[vid(i, 0, k), vid(i + 1, 0, k), vid(i, 0, k + 1), vid(i + 1, 0, k + 1)]);
            qd.push(3);
            quad.extend_from_slice(&[vid(i, ny, k), vid(i + 1, ny, k), vid(i, ny, k + 1), vid(i + 1, ny, k + 1)]);
            qd.push(4);
        }
    }
    let np = pd.len();
    CellMesh { sdim: 3, points: grid_points(nx, ny, nz, lx, ly, lz), types: vec![et("tri", 3, tri, td), et("quad", 4, quad, qd), et("prism", 6, prism, vec![1; np])] }
}

fn volume(m: &SimplexMesh) -> f64 {
    let g = femcore::fem::geometry(m);
    g.vol.iter().sum()
}

/// Every face of a cell is shared by one or two cells; faces used once are exactly the boundary facets.
fn assert_conforming(m: &SimplexMesh) {
    let mut faces: HashMap<[u32; 3], usize> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            f.sort_unstable();
            *faces.entry([f[0], f[1], f[2]]).or_insert(0) += 1;
        }
    }
    assert!(faces.values().all(|&n| n <= 2), "a face is shared by more than two cells");
    let mut boundary: Vec<[u32; 3]> = faces.iter().filter(|(_, &n)| n == 1).map(|(f, _)| *f).collect();
    let mut given: Vec<[u32; 3]> = (0..m.n_facets())
        .map(|f| {
            let mut v = m.facet(f).to_vec();
            v.sort_unstable();
            [v[0], v[1], v[2]]
        })
        .collect();
    boundary.sort();
    given.sort();
    assert_eq!(boundary, given, "the boundary of the split cells differs from the boundary facets");
}

#[test]
fn quads_become_triangles_and_solve_exactly() {
    let (nx, ny) = (6usize, 3usize);
    let vid = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    let mut points = vec![];
    for j in 0..=ny {
        for i in 0..=nx {
            points.extend_from_slice(&[i as f64 / nx as f64, 0.2 * j as f64 / ny as f64]);
        }
    }
    let (mut quad, mut edg, mut ed) = (vec![], vec![], vec![]);
    for j in 0..ny {
        for i in 0..nx {
            quad.extend_from_slice(&[vid(i, j), vid(i + 1, j), vid(i, j + 1), vid(i + 1, j + 1)]);
        }
        edg.extend_from_slice(&[vid(0, j), vid(0, j + 1)]);
        ed.push(1);
        edg.extend_from_slice(&[vid(nx, j), vid(nx, j + 1)]);
        ed.push(4);
    }
    for i in 0..nx {
        edg.extend_from_slice(&[vid(i, 0), vid(i + 1, 0)]);
        ed.push(2);
        edg.extend_from_slice(&[vid(i, ny), vid(i + 1, ny)]);
        ed.push(3);
    }
    let nq = nx * ny;
    let nm = CellMesh { sdim: 2, points, types: vec![et("quad", 4, quad, vec![1; nq]), et("edg", 2, edg, ed)] };
    let (m, rep) = SimplexMesh::from_cells_split(&native::<2>(&nm)).unwrap();
    assert_eq!((rep.quads, m.n_cells()), (nq, 2 * nq));
    assert!((volume(&m) - 0.2).abs() < 1e-14);
    assert_eq!(m.boundary_labels(), vec![1, 2, 3, 4]);
    // uniaxial plane stress: ux = sigma x / E, exact in P2
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.plane = Plane::Stress;
    spec.materials.default.insert("E".into(), 1e9);
    spec.materials.default.insert("nu".into(), 0.0);
    spec.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![4], kind: BcKind::Load([1e6, 0.0, 0.0]) },
    ];
    let r = solve_solid(&m, &spec, &SolveState::default()).unwrap();
    let ux = &r.solutions[0].field("ux").unwrap().values;
    let x = m.dof_coordinates(2);
    let exact: Vec<f64> = (0..ux.len()).map(|i| 1e6 * x[2 * i] / 1e9).collect();
    assert!(common::maxdiff(ux, &exact) < 1e-12);
}

fn uniaxial_3d(m: &SimplexMesh, volume_expected: f64) {
    assert!((volume(m) - volume_expected).abs() < 1e-12 * volume_expected);
    assert_conforming(m);
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.materials.default.insert("E".into(), 1e9);
    spec.materials.default.insert("nu".into(), 0.0);
    spec.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![3], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![5], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        Bc { on: None, labels: vec![2], kind: BcKind::Load([1e6, 0.0, 0.0]) },
    ];
    let r = solve_solid(m, &spec, &SolveState::default()).unwrap();
    let ux = &r.solutions[0].field("ux").unwrap().values;
    let x = m.dof_coordinates(2);
    let exact: Vec<f64> = (0..ux.len()).map(|i| 1e6 * x[3 * i] / 1e9).collect();
    assert!(common::maxdiff(ux, &exact) < 1e-12, "{}", common::maxdiff(ux, &exact));
}

#[test]
fn hexahedra_are_split_conformingly() {
    let nm = hex_block(3, 2, 2, 1.0, 0.4, 0.3);
    // star split: vertices at face and cell centres, 24 tetrahedra per hexahedron
    let (m, rep) = SimplexMesh::from_cells_split_with(&native::<3>(&nm), SplitMode::Star).unwrap();
    assert_eq!((rep.hexes, rep.boundary_quads, m.n_cells()), (12, 2 * (2 * 2 + 3 * 2 + 3 * 2), 12 * 24));
    uniaxial_3d(&m, 1.0 * 0.4 * 0.3);
    assert_eq!(m.boundary_labels(), vec![1, 2, 3, 4, 5, 6]);
    // minimal split: no new vertices, 6 tetrahedra per hexahedron, 2 triangles per boundary quad
    let (m, rep) = SimplexMesh::from_cells_split_with(&native::<3>(&nm), SplitMode::Minimal).unwrap();
    assert_eq!((rep.hexes, m.n_cells(), m.n_facets(), m.n_points()), (12, 12 * 6, 2 * 32, nm.n_points()));
    uniaxial_3d(&m, 1.0 * 0.4 * 0.3);
    assert_eq!(m.boundary_labels(), vec![1, 2, 3, 4, 5, 6]);
}

#[test]
fn prisms_are_split_conformingly() {
    let nm = prism_block(2, 2, 2, 1.0, 0.4, 0.3);
    let (m, rep) = SimplexMesh::from_cells_split_with(&native::<3>(&nm), SplitMode::Star).unwrap();
    assert_eq!((rep.prisms, m.n_cells()), (16, 16 * 14));
    uniaxial_3d(&m, 1.0 * 0.4 * 0.3);
    let (m, rep) = SimplexMesh::from_cells_split_with(&native::<3>(&nm), SplitMode::Minimal).unwrap();
    assert_eq!((rep.prisms, m.n_cells(), m.n_points()), (16, 16 * 3, nm.n_points()));
    uniaxial_3d(&m, 1.0 * 0.4 * 0.3);
}

#[test]
fn minimal_split_conforms_for_any_vertex_numbering() {
    // the diagonal rule uses only global vertex numbers, so a scrambled numbering of the same
    // hexahedral block must still give a conforming split that fills the same volume
    let nm = hex_block(3, 3, 2, 1.0, 1.0, 0.6);
    let n = nm.n_points() as u32;
    for seed in [1u64, 7, 99, 12345] {
        let mut perm: Vec<u32> = (0..n).collect();
        let mut x = seed;
        for i in (1..perm.len()).rev() {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            perm.swap(i, (x >> 33) as usize % (i + 1));
        }
        let mut points = vec![0.0; nm.points.len()];
        for v in 0..n as usize {
            points[perm[v] as usize * 3..perm[v] as usize * 3 + 3].copy_from_slice(&nm.points[v * 3..v * 3 + 3]);
        }
        let types = nm.types.iter().map(|t| et(&t.name, t.nep, t.elem.iter().map(|&v| perm[v as usize]).collect(), t.dom.clone())).collect();
        let sc = CellMesh { sdim: 3, points, types };
        let (m, _) = SimplexMesh::from_cells_split_with(&native::<3>(&sc), SplitMode::Minimal).unwrap();
        assert_eq!(m.n_cells(), 18 * 6);
        assert!((volume(&m) - 0.6).abs() < 1e-12, "seed {seed}");
        assert_conforming(&m);
    }
}

#[test]
fn hexahedra_and_prisms_sharing_quad_faces_conform() {
    // a hex block on the left (x in [0, 1]) and a prism block on the right (x in [1, 2])
    let h = hex_block(2, 2, 2, 1.0, 0.4, 0.3);
    let p = prism_block(2, 2, 2, 1.0, 0.4, 0.3);
    let nh = h.n_points() as u32;
    let mut points = h.points.clone();
    // prism points shifted by 1 in x, sharing the plane x = 1 (the last column of the hex grid)
    // with the first column of the prism grid
    let (nx, ny, nz) = (2usize, 2usize, 2usize);
    let shared = |j: usize, k: usize| ((k * (ny + 1) + j) * (nx + 1) + nx) as u32; // hex vertex at x = 1
    let mut map: HashMap<u32, u32> = HashMap::new();
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                let pid = ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
                if i == 0 {
                    map.insert(pid, shared(j, k));
                } else {
                    map.insert(pid, nh + (points.len() / 3 - nh as usize) as u32);
                    points.extend_from_slice(&p.points[pid as usize * 3..pid as usize * 3 + 3]);
                    let l = points.len();
                    points[l - 3] += 1.0;
                }
            }
        }
    }
    let mut types = vec![];
    // hex cells and their outer boundary quads except x = 1, which is interior now
    let hq = h.ty("quad").unwrap();
    let (mut quad, mut qd) = (vec![], vec![]);
    for i in 0..hq.len() {
        if hq.dom[i] != 2 {
            quad.extend_from_slice(hq.row(i));
            qd.push(hq.dom[i]);
        }
    }
    // prism side/boundary quads and triangles (excluding the x = 0 side, interior now)
    let pq = p.ty("quad").unwrap();
    for i in 0..pq.len() {
        if pq.dom[i] != 1 {
            quad.extend(pq.row(i).iter().map(|v| map[v]));
            qd.push(pq.dom[i]);
        }
    }
    let pt = p.ty("tri").unwrap();
    let tri: Vec<u32> = pt.elem.iter().map(|v| map[v]).collect();
    types.push(et("quad", 4, quad, qd));
    types.push(et("tri", 3, tri, pt.dom.clone()));
    types.push(et("hex", 8, h.ty("hex").unwrap().elem.clone(), h.ty("hex").unwrap().dom.clone()));
    let pp = p.ty("prism").unwrap();
    types.push(et("prism", 6, pp.elem.iter().map(|v| map[v]).collect(), vec![2; pp.len()]));
    let nm = CellMesh { sdim: 3, points, types };
    let (m, rep) = SimplexMesh::from_cells_split(&native::<3>(&nm)).unwrap();
    assert_eq!((rep.hexes, rep.prisms), (8, 16));
    assert!((volume(&m) - 2.0 * 0.4 * 0.3).abs() < 1e-12);
    assert_conforming(&m);
    assert_eq!(m.domain_labels(), vec![1, 2]);
}
