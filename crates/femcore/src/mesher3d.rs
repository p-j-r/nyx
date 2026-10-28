//! Tetrahedral meshes of 3D solids built from primitives (box, sphere, cylinder, extruded 2D
//! outline), painted in order like the 2D mesher: later bodies override earlier ones and a body
//! with domain 0 or less is a hole.
//!
//! The method is isosurface stuffing (Labelle and Shewchuk 2007, with the vertex-warping and
//! marching-tetrahedra cuts of Quartet): the solid is filled with the tetrahedra of a body-centred
//! cubic lattice, which have good shape; a lattice vertex that lies within `alpha` of an edge
//! length of the surface is moved onto it, and the few tetrahedra that the surface still crosses
//! are cut along it. It is fast because every step is local and independent: the vertices of the
//! lattice are classified in parallel, tetrahedra are generated in parallel blocks, only
//! the thin layer next to the surface needs any bookkeeping, and cut points are identified by the
//! pair of lattice vertices of their edge so no search structure is needed.
//!
//! The lattice is aligned with the bounding box of the solid, so faces of boxes that coincide with
//! the bounding box are reproduced exactly; other surfaces are represented piecewise linearly with a
//! geometric error of order h^2. Sharp edges and corners of the solid (where two or three faces of the
//! bodies meet) are fitted: a lattice vertex near such a crease is moved onto it.
//!
//! Material interfaces are fitted too. Besides the outer surface of the solid, the surface of every
//! positive body is a cutting surface when there are several bodies: each tetrahedron that crosses one is
//! cut along it (both sides are kept), so every final cell lies in one body and takes its domain number.
//! Cut points are identified by a hash of the two end points and the surface, so all tetrahedra
//! agree on shared faces without any shared state.
//!
//! Boundary numbers: bodies in order, each numbering its faces. Box: 1 x=min, 2 x=max, 3 y=min,
//! 4 y=max, 5 z=min, 6 z=max (the numbering of `SimplexMesh::block`). Cylinder: 1 side, 2 base,
//! 3 top. Sphere: 1. Extrusion: the edges of the outline (as in the 2D mesher), then bottom, top.

use crate::mesh::SimplexMesh;
use crate::mesher::{dist_point_segment, Shape};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Solid {
    Box { min: [f64; 3], max: [f64; 3] },
    Sphere { c: [f64; 3], r: f64 },
    /// circular cylinder along axis 0, 1 or 2 from the disc centred at `base` to `base + len` along the axis
    Cylinder { base: [f64; 3], axis: usize, r: f64, len: f64 },
    /// a 2D outline in the x-y plane extruded from z0 to z1
    Extrusion { profile: Shape, z0: f64, z1: f64 },
}

fn len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Signed distance of a 2D outline (negative inside) and the number of its edges.
fn sdf2(shape: &Shape, p: [f64; 2]) -> f64 {
    match shape {
        Shape::Rect { x0, y0, x1, y1 } => {
            let (cx, cy) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
            let (hx, hy) = (0.5 * (x1 - x0).abs(), 0.5 * (y1 - y0).abs());
            let (qx, qy) = ((p[0] - cx).abs() - hx, (p[1] - cy).abs() - hy);
            (qx.max(0.0).powi(2) + qy.max(0.0).powi(2)).sqrt() + qx.max(qy).min(0.0)
        }
        Shape::Circle { cx, cy, r } => ((p[0] - cx).powi(2) + (p[1] - cy).powi(2)).sqrt() - r,
        Shape::Polygon(v) => {
            let n = v.len();
            let d = (0..n).map(|i| dist_point_segment(p, v[i], v[(i + 1) % n])).fold(f64::INFINITY, f64::min);
            if shape.contains(p) { -d } else { d }
        }
    }
}

fn outline_edges(shape: &Shape) -> i32 {
    match shape {
        Shape::Rect { .. } => 4,
        Shape::Circle { .. } => 1,
        Shape::Polygon(v) => v.len() as i32,
    }
}

/// Which edge (1-based) of the outline is nearest to `p`, with the 2D numbering of the mesher.
fn outline_edge_at(shape: &Shape, p: [f64; 2]) -> i32 {
    match shape {
        Shape::Rect { x0, y0, x1, y1 } => {
            let (xa, xb, ya, yb) = (x0.min(*x1), x0.max(*x1), y0.min(*y1), y0.max(*y1));
            let d = [(p[0] - xa).abs(), (p[1] - ya).abs(), (p[1] - yb).abs(), (p[0] - xb).abs()];
            (0..4).min_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap()).unwrap() as i32 + 1
        }
        Shape::Circle { .. } => 1,
        Shape::Polygon(v) => {
            let n = v.len();
            (0..n).min_by(|&a, &b| dist_point_segment(p, v[a], v[(a + 1) % n]).partial_cmp(&dist_point_segment(p, v[b], v[(b + 1) % n])).unwrap()).unwrap_or(0) as i32 + 1
        }
    }
}

impl Solid {
    pub fn name(&self) -> &'static str {
        match self {
            Solid::Box { .. } => "Box",
            Solid::Sphere { .. } => "Sphere",
            Solid::Cylinder { .. } => "Cylinder",
            Solid::Extrusion { .. } => "Extrusion",
        }
    }

    /// Signed distance to the surface, negative inside (exact for the primitives, a bound for the rest).
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        match self {
            Solid::Box { min, max } => {
                let mut q = [0.0; 3];
                for k in 0..3 {
                    q[k] = (p[k] - 0.5 * (min[k] + max[k])).abs() - 0.5 * (max[k] - min[k]).abs();
                }
                len3([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]) + q[0].max(q[1]).max(q[2]).min(0.0)
            }
            Solid::Sphere { c, r } => len3([p[0] - c[0], p[1] - c[1], p[2] - c[2]]) - r,
            Solid::Cylinder { base, axis, r, len } => {
                let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
                let radial = ((p[a1] - base[a1]).powi(2) + (p[a2] - base[a2]).powi(2)).sqrt() - r;
                let axial = (p[*axis] - base[*axis] - 0.5 * len).abs() - 0.5 * len.abs();
                radial.max(0.0).hypot(axial.max(0.0)) + radial.max(axial).min(0.0)
            }
            Solid::Extrusion { profile, z0, z1 } => {
                let d = sdf2(profile, [p[0], p[1]]);
                let dz = (p[2] - 0.5 * (z0 + z1)).abs() - 0.5 * (z1 - z0).abs();
                d.max(0.0).hypot(dz.max(0.0)) + d.max(dz).min(0.0)
            }
        }
    }

    pub fn n_faces(&self) -> i32 {
        match self {
            Solid::Box { .. } => 6,
            Solid::Sphere { .. } => 1,
            Solid::Cylinder { .. } => 3,
            Solid::Extrusion { profile, .. } => outline_edges(profile) + 2,
        }
    }

    /// The face (1-based) of this solid nearest to `p`.
    pub fn face_at(&self, p: [f64; 3]) -> i32 {
        match self {
            Solid::Box { min, max } => {
                let mut best = (f64::INFINITY, 1);
                for k in 0..3 {
                    for (side, v) in [(0, min[k].min(max[k])), (1, min[k].max(max[k]))] {
                        let d = (p[k] - v).abs();
                        if d < best.0 {
                            best = (d, 2 * k as i32 + side + 1);
                        }
                    }
                }
                best.1
            }
            Solid::Sphere { .. } => 1,
            Solid::Cylinder { base, axis, r, len } => {
                let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
                let radial = (((p[a1] - base[a1]).powi(2) + (p[a2] - base[a2]).powi(2)).sqrt() - r).abs();
                let (d0, d1) = ((p[*axis] - base[*axis]).abs(), (p[*axis] - base[*axis] - len).abs());
                if radial <= d0 && radial <= d1 { 1 } else if d0 <= d1 { 2 } else { 3 }
            }
            Solid::Extrusion { profile, z0, z1 } => {
                let ne = outline_edges(profile);
                let dside = sdf2(profile, [p[0], p[1]]).abs();
                let (d0, d1) = ((p[2] - z0).abs(), (p[2] - z1).abs());
                if dside <= d0 && dside <= d1 { outline_edge_at(profile, [p[0], p[1]]) } else if d0 <= d1 { ne + 1 } else { ne + 2 }
            }
        }
    }

    pub fn bbox(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            Solid::Box { min, max } => ([min[0].min(max[0]), min[1].min(max[1]), min[2].min(max[2])], [min[0].max(max[0]), min[1].max(max[1]), min[2].max(max[2])]),
            Solid::Sphere { c, r } => ([c[0] - r, c[1] - r, c[2] - r], [c[0] + r, c[1] + r, c[2] + r]),
            Solid::Cylinder { base, axis, r, len } => {
                let (mut lo, mut hi) = (*base, *base);
                for k in 0..3 {
                    if k == *axis {
                        lo[k] = base[k].min(base[k] + len);
                        hi[k] = base[k].max(base[k] + len);
                    } else {
                        lo[k] = base[k] - r;
                        hi[k] = base[k] + r;
                    }
                }
                (lo, hi)
            }
            Solid::Extrusion { profile, z0, z1 } => {
                let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
                let mut add = |p: [f64; 2]| {
                    for k in 0..2 {
                        lo[k] = lo[k].min(p[k]);
                        hi[k] = hi[k].max(p[k]);
                    }
                };
                match profile {
                    Shape::Rect { x0, y0, x1, y1 } => {
                        add([*x0, *y0]);
                        add([*x1, *y1]);
                    }
                    Shape::Circle { cx, cy, r } => {
                        add([cx - r, cy - r]);
                        add([cx + r, cy + r]);
                    }
                    Shape::Polygon(v) => v.iter().for_each(|p| add(*p)),
                }
                ([lo[0], lo[1], z0.min(*z1)], [hi[0], hi[1], z0.max(*z1)])
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Body3d {
    pub solid: Solid,
    /// domain number of the cells inside this body (0 or negative: a hole, no cells)
    pub domain: i32,
    /// target element size inside this body when finer than the global one (local refinement)
    #[serde(default)]
    pub h: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Spec3d {
    pub bodies: Vec<Body3d>,
    /// target element (lattice) size
    pub h: f64,
    /// a lattice vertex within this fraction of an edge length of the surface is moved onto it
    /// (0.3 gives dihedral angles of at least about 10 degrees on smooth surfaces)
    #[serde(default = "default_alpha")]
    pub alpha: f64,
    /// fit sharp edges and corners of the bodies (vertices near a crease move onto it)
    #[serde(default = "default_true")]
    pub snap_features: bool,
    /// a body thinner than this many elements along an axis gets a finer lattice along that axis
    /// (at most 4 times finer than `h`), so that it holds at least this many layers; 0 turns it off
    #[serde(default = "default_layers")]
    pub min_layers: f64,
}

fn default_layers() -> f64 {
    2.0
}

fn default_true() -> bool {
    true
}

fn default_alpha() -> f64 {
    0.3
}

impl Default for Spec3d {
    fn default() -> Self {
        Spec3d { bodies: vec![], h: 0.1, alpha: default_alpha(), snap_features: true, min_layers: default_layers() }
    }
}

/// A report of how the mesh was made.
#[derive(Clone, Debug, Default)]
pub struct Report3d {
    pub lattice_vertices: usize,
    pub warped: usize,
    pub cut_points: usize,
    pub cells: usize,
    /// edges split by local refinement
    pub refined: usize,
}

struct Field<'a> {
    bodies: &'a [Body3d],
}

impl Field<'_> {
    /// Negative inside the painted solid, positive outside, zero on its surface.
    fn f(&self, p: [f64; 3]) -> f64 {
        let mut f = f64::INFINITY;
        for b in self.bodies {
            let s = b.solid.sdf(p);
            if b.domain > 0 {
                f = f.min(s);
            } else {
                f = f.max(-s);
            }
        }
        f
    }

    /// Domain number at `p`: the last body that contains it (0 for a hole or nothing).
    fn domain(&self, p: [f64; 3]) -> i32 {
        for b in self.bodies.iter().rev() {
            if b.solid.sdf(p) < 0.0 {
                return b.domain.max(0);
            }
        }
        0
    }

    /// Domain of a cell centroid; a centroid just outside the (piecewise linear) solid takes the
    /// domain of the nearest positive body.
    fn domain_near(&self, p: [f64; 3]) -> i32 {
        let d = self.domain(p);
        if d > 0 {
            return d;
        }
        self.bodies.iter().filter(|b| b.domain > 0).min_by(|a, b| a.solid.sdf(p).partial_cmp(&b.solid.sdf(p)).unwrap()).map_or(1, |b| b.domain)
    }

    /// Boundary number of the surface point `p`: the face of the body whose surface is nearest.
    fn boundary_label(&self, p: [f64; 3]) -> i32 {
        let mut offset = 0;
        let mut best = (f64::INFINITY, 1);
        for b in self.bodies {
            let d = b.solid.sdf(p).abs();
            if d < best.0 {
                best = (d, offset + b.solid.face_at(p));
            }
            offset += b.solid.n_faces();
        }
        best.1
    }
}

const TAG: u64 = 1 << 63;

/// A vertex of a piece during cutting: its id (lattice number, or a tagged hash for a cut point),
/// its position, and whether it lies on the outer surface of the solid.
#[derive(Clone, Copy, Debug)]
struct V {
    id: u64,
    p: [f64; 3],
    surf: bool,
    /// class (-1 inside, 0 on, +1 outside) with respect to every cutting surface
    cl: [i8; MAX_SURFACES],
}

/// Largest number of cutting surfaces: the outer surface of the solid and the faces of the bodies (interfaces).
const MAX_SURFACES: usize = 64;

fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

/// Identity of the cut point on the segment between vertices `a` and `b` for cutting surface `s`.
fn cut_id(a: u64, b: u64, s: usize) -> u64 {
    let (lo, hi) = (a.min(b), a.max(b));
    TAG | (mix64(lo ^ mix64(hi ^ mix64(s as u64 + 1))) >> 1)
}

/// Splits a prism (triangles `v[0..3]` and `v[3..6]`, vertical edges v[i]-v[i+3]) into three
/// tetrahedra by coning from its lowest-numbered vertex; quadrilateral faces are cut along the
/// diagonal from their lowest-numbered vertex so neighbours agree. Ids give the global order.
fn split_prism(v: [V; 6], out: &mut Vec<[V; 4]>) {
    const QUADS: [[usize; 4]; 3] = [[0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]];
    let m = (0..6).min_by_key(|&i| v[i].id).unwrap();
    for q in QUADS {
        if q.contains(&m) {
            continue;
        }
        let k = (0..4).min_by_key(|&i| v[q[i]].id).unwrap();
        let c = |j: usize| v[q[(k + j) % 4]];
        out.push([v[m], c(0), c(1), c(2)]);
        out.push([v[m], c(0), c(2), c(3)]);
    }
    for tri in [[0usize, 1, 2], [3, 4, 5]] {
        if !tri.contains(&m) {
            out.push([v[m], v[tri[0]], v[tri[1]], v[tri[2]]]);
        }
    }
}

/// The two triangles of a quadrilateral (cyclic order) along the diagonal from its lowest id.
fn quad_tris(q: [V; 4], out: &mut Vec<[V; 3]>) {
    let k = (0..4).min_by_key(|&i| q[i].id).unwrap();
    let c = |j: usize| q[(k + j) % 4];
    out.push([c(0), c(1), c(2)]);
    out.push([c(0), c(2), c(3)]);
}

struct Lattice {
    n: [usize; 3],
    origin: [f64; 3],
    h: [f64; 3],
    n_corners: usize,
}

impl Lattice {
    fn corner(&self, i: usize, j: usize, k: usize) -> u32 {
        ((k * (self.n[1] + 1) + j) * (self.n[0] + 1) + i) as u32
    }
    fn center(&self, i: usize, j: usize, k: usize) -> u32 {
        (self.n_corners + (k * self.n[1] + j) * self.n[0] + i) as u32
    }
    fn n_vertices(&self) -> usize {
        self.n_corners + self.n[0] * self.n[1] * self.n[2]
    }
    fn pos(&self, v: u32) -> [f64; 3] {
        let v = v as usize;
        if v < self.n_corners {
            let (nx, ny) = (self.n[0] + 1, self.n[1] + 1);
            let (i, j, k) = (v % nx, (v / nx) % ny, v / (nx * ny));
            [self.origin[0] + i as f64 * self.h[0], self.origin[1] + j as f64 * self.h[1], self.origin[2] + k as f64 * self.h[2]]
        } else {
            let c = v - self.n_corners;
            let (i, j, k) = (c % self.n[0], (c / self.n[0]) % self.n[1], c / (self.n[0] * self.n[1]));
            [self.origin[0] + (i as f64 + 0.5) * self.h[0], self.origin[1] + (j as f64 + 0.5) * self.h[1], self.origin[2] + (k as f64 + 0.5) * self.h[2]]
        }
    }
    /// The (up to 14) lattice neighbours of a vertex.
    fn neighbours(&self, v: u32, out: &mut Vec<u32>) {
        out.clear();
        let vi = v as usize;
        let n = self.n;
        if vi < self.n_corners {
            let (nx, ny) = (n[0] + 1, n[1] + 1);
            let (i, j, k) = ((vi % nx) as i64, ((vi / nx) % ny) as i64, (vi / (nx * ny)) as i64);
            for (di, dj, dk) in [(-1i64, 0i64, 0i64), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)] {
                let (a, b, c) = (i + di, j + dj, k + dk);
                if a >= 0 && b >= 0 && c >= 0 && a <= n[0] as i64 && b <= n[1] as i64 && c <= n[2] as i64 {
                    out.push(self.corner(a as usize, b as usize, c as usize));
                }
            }
            for (di, dj, dk) in [(-1i64, -1i64, -1i64), (0, -1, -1), (-1, 0, -1), (0, 0, -1), (-1, -1, 0), (0, -1, 0), (-1, 0, 0), (0, 0, 0)] {
                let (a, b, c) = (i + di, j + dj, k + dk);
                if a >= 0 && b >= 0 && c >= 0 && a < n[0] as i64 && b < n[1] as i64 && c < n[2] as i64 {
                    out.push(self.center(a as usize, b as usize, c as usize));
                }
            }
        } else {
            let c = vi - self.n_corners;
            let (i, j, k) = ((c % n[0]) as i64, ((c / n[0]) % n[1]) as i64, (c / (n[0] * n[1])) as i64);
            for (di, dj, dk) in [(-1i64, 0i64, 0i64), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)] {
                let (a, b, cc) = (i + di, j + dj, k + dk);
                if a >= 0 && b >= 0 && cc >= 0 && a < n[0] as i64 && b < n[1] as i64 && cc < n[2] as i64 {
                    out.push(self.center(a as usize, b as usize, cc as usize));
                }
            }
            for (di, dj, dk) in [(0i64, 0i64, 0i64), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)] {
                out.push(self.corner((i + di) as usize, (j + dj) as usize, (k + dk) as usize));
            }
        }
    }
}

fn threads() -> usize {
    crate::par::threads()
}

/// Runs `f(range)` on blocks of `0..n` in parallel and returns the results in order.
fn par_blocks<T: Send, F: Fn(std::ops::Range<usize>) -> T + Sync>(n: usize, min_block: usize, f: F) -> Vec<T> {
    let nt = threads().min(n.div_ceil(min_block.max(1))).max(1);
    let chunk = n.div_ceil(nt).max(1);
    use rayon::prelude::*;
    (0..n.div_ceil(chunk)).into_par_iter().map(|t| f(t * chunk..((t + 1) * chunk).min(n))).collect()
}

/// A face of a body extended to an infinite surface (a plane, a cylinder or a sphere): the building
/// block for finding where faces meet. `value` is negative inside the body.
#[derive(Clone, Debug)]
pub enum Surf {
    Plane { n: [f64; 3], d: f64 },
    Cyl { c: [f64; 3], axis: usize, r: f64 },
    Sph { c: [f64; 3], r: f64 },
}

impl Surf {
    fn value(&self, p: [f64; 3]) -> f64 {
        match self {
            Surf::Plane { n, d } => n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - d,
            Surf::Cyl { c, axis, r } => {
                let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
                ((p[a1] - c[a1]).powi(2) + (p[a2] - c[a2]).powi(2)).sqrt() - r
            }
            Surf::Sph { c, r } => len3([p[0] - c[0], p[1] - c[1], p[2] - c[2]]) - r,
        }
    }
    fn grad(&self, p: [f64; 3]) -> [f64; 3] {
        match self {
            Surf::Plane { n, .. } => *n,
            Surf::Cyl { c, axis, .. } => {
                let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
                let (x, y) = (p[a1] - c[a1], p[a2] - c[a2]);
                let l = (x * x + y * y).sqrt().max(1e-300);
                let mut g = [0.0; 3];
                g[a1] = x / l;
                g[a2] = y / l;
                g
            }
            Surf::Sph { c, .. } => {
                let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
                let l = len3(d).max(1e-300);
                [d[0] / l, d[1] / l, d[2] / l]
            }
        }
    }
    /// Projection of `p` onto the surface (one Newton step along the gradient).
    fn project(&self, p: [f64; 3]) -> [f64; 3] {
        let (v, g) = (self.value(p), self.grad(p));
        let l2 = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).max(1e-300);
        [p[0] - v * g[0] / l2, p[1] - v * g[1] / l2, p[2] - v * g[2] / l2]
    }
}

/// A face of a body as an infinite surface together with the box that contains the face itself
/// (the surface only matters where the body has this face).
#[derive(Clone, Debug)]
pub struct FaceSurf {
    pub surf: Surf,
    pub lo: [f64; 3],
    pub hi: [f64; 3],
}

impl FaceSurf {
    fn contains(&self, p: [f64; 3], margin: f64) -> bool {
        (0..3).all(|k| p[k] >= self.lo[k] - margin && p[k] <= self.hi[k] + margin)
    }
    fn same_plane(&self, o: &FaceSurf) -> bool {
        match (&self.surf, &o.surf) {
            (Surf::Plane { n, d }, Surf::Plane { n: n2, d: d2 }) => {
                let same = (0..3).all(|k| (n[k] - n2[k]).abs() < 1e-12) && (d - d2).abs() < 1e-12 * (1.0 + d.abs());
                let flipped = (0..3).all(|k| (n[k] + n2[k]).abs() < 1e-12) && (d + d2).abs() < 1e-12 * (1.0 + d.abs());
                same || flipped
            }
            _ => false,
        }
    }
}

impl Solid {
    /// The faces with the box that contains each of them.
    pub fn faces_with_extent(&self) -> Vec<FaceSurf> {
        let (blo, bhi) = self.bbox();
        match self {
            Solid::Box { min, max } => {
                let mut v = vec![];
                for (k, surf) in self.face_surfaces().into_iter().enumerate() {
                    let axis = k / 2;
                    let (mut lo, mut hi) = (blo, bhi);
                    let at = if k % 2 == 0 { min[axis].min(max[axis]) } else { min[axis].max(max[axis]) };
                    lo[axis] = at;
                    hi[axis] = at;
                    v.push(FaceSurf { surf, lo, hi });
                }
                v
            }
            Solid::Cylinder { base, axis, len, .. } => {
                let mut v = vec![];
                for (k, surf) in self.face_surfaces().into_iter().enumerate() {
                    let (mut lo, mut hi) = (blo, bhi);
                    if k == 1 {
                        lo[*axis] = base[*axis];
                        hi[*axis] = base[*axis];
                    } else if k == 2 {
                        lo[*axis] = base[*axis] + len;
                        hi[*axis] = base[*axis] + len;
                    }
                    v.push(FaceSurf { surf, lo, hi });
                }
                v
            }
            _ => self.face_surfaces().into_iter().map(|surf| FaceSurf { surf, lo: blo, hi: bhi }).collect(),
        }
    }

    /// The faces of the body as infinite surfaces (for fitting sharp edges and corners).
    pub fn face_surfaces(&self) -> Vec<Surf> {
        let plane = |n: [f64; 3], d: f64| Surf::Plane { n, d };
        match self {
            Solid::Box { min, max } => {
                let mut v = vec![];
                for k in 0..3 {
                    let mut n = [0.0; 3];
                    n[k] = 1.0;
                    let (lo, hi) = (min[k].min(max[k]), min[k].max(max[k]));
                    v.push(plane([-n[0], -n[1], -n[2]], -lo));
                    v.push(plane(n, hi));
                }
                v
            }
            Solid::Sphere { c, r } => vec![Surf::Sph { c: *c, r: *r }],
            Solid::Cylinder { base, axis, r, len } => {
                let mut n = [0.0; 3];
                n[*axis] = if *len >= 0.0 { 1.0 } else { -1.0 };
                let (b0, b1) = (base[*axis], base[*axis] + len);
                vec![
                    Surf::Cyl { c: *base, axis: *axis, r: *r },
                    plane([-n[0], -n[1], -n[2]], -(n[*axis] * b0)),
                    plane(n, n[*axis] * b1),
                ]
            }
            Solid::Extrusion { profile, z0, z1 } => {
                let mut v = vec![plane([0.0, 0.0, -1.0], -z0.min(*z1)), plane([0.0, 0.0, 1.0], z0.max(*z1))];
                match profile {
                    Shape::Rect { x0, y0, x1, y1 } => {
                        let (xa, xb, ya, yb) = (x0.min(*x1), x0.max(*x1), y0.min(*y1), y0.max(*y1));
                        v.push(plane([-1.0, 0.0, 0.0], -xa));
                        v.push(plane([1.0, 0.0, 0.0], xb));
                        v.push(plane([0.0, -1.0, 0.0], -ya));
                        v.push(plane([0.0, 1.0, 0.0], yb));
                    }
                    Shape::Circle { cx, cy, r } => v.push(Surf::Cyl { c: [*cx, *cy, 0.0], axis: 2, r: *r }),
                    Shape::Polygon(pts) => {
                        let n = pts.len();
                        let area: f64 = (0..n).map(|i| pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1]).sum::<f64>();
                        let sign = if area >= 0.0 { 1.0 } else { -1.0 };
                        for i in 0..n {
                            let (a, b) = (pts[i], pts[(i + 1) % n]);
                            let (ex, ey) = (b[0] - a[0], b[1] - a[1]);
                            let l = (ex * ex + ey * ey).sqrt().max(1e-300);
                            // outward normal of a counter-clockwise polygon is (ey, -ex)
                            let nn = [sign * ey / l, -sign * ex / l, 0.0];
                            v.push(plane(nn, nn[0] * a[0] + nn[1] * a[1]));
                        }
                    }
                }
                v
            }
        }
    }
}

/// The cutting surfaces: 0 is the boundary of the whole solid; the others are the surfaces of the
/// positive bodies when there are several (the interfaces between materials).
struct Surfaces<'a> {
    field: Field<'a>,
    faces: Vec<FaceSurf>,
    /// probe distance for deciding that a point on a contact face between bodies is interior
    delta: f64,
}

impl Surfaces<'_> {
    fn n(&self) -> usize {
        1 + self.faces.len()
    }
    fn value(&self, s: usize, p: [f64; 3]) -> f64 {
        if s == 0 { self.field.f(p) } else { self.faces[s - 1].surf.value(p) }
    }
    fn class(&self, s: usize, p: [f64; 3], tol: f64) -> i8 {
        let f = self.value(s, p);
        if f < -tol {
            -1
        } else if f > tol {
            1
        } else if s == 0 && self.interior_contact(p) {
            // two bodies touch here: the face between them is not a boundary of the solid
            -1
        } else {
            0
        }
    }
    /// True if a point where the solid's field vanishes is nevertheless interior: moving a little along
    /// every axis in both directions stays inside the solid (the point lies on a contact face between bodies).
    fn interior_contact(&self, p: [f64; 3]) -> bool {
        for k in 0..3 {
            for sgn in [-1.0, 1.0] {
                let mut q = p;
                q[k] += sgn * self.delta;
                if self.field.f(q) >= 0.0 {
                    return false;
                }
            }
        }
        true
    }
}

/// Point on the segment `p0`-`p1` where the sign of surface `s` changes (bisection).
fn bisect_surface(sf: &Surfaces, s: usize, p0: [f64; 3], p1: [f64; 3]) -> (f64, [f64; 3]) {
    let neg0 = sf.value(s, p0) < 0.0;
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    let at = |t: f64| [p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1]), p0[2] + t * (p1[2] - p0[2])];
    for _ in 0..44 {
        let mid = 0.5 * (lo + hi);
        if (sf.value(s, at(mid)) < 0.0) == neg0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, at(t))
}

/// The point where several faces meet nearest to `v` (alternating projections onto the faces), if it is
/// a point of the solid's surface within `reach` of `v`.
fn snap_to_crease(faces: &[&FaceSurf], field: &Field, v: [f64; 3], reach: f64, tol_surface: f64) -> Option<[f64; 3]> {
    let mut p = v;
    for _ in 0..12 {
        for f in faces {
            p = f.surf.project(p);
        }
    }
    if faces.iter().any(|f| f.surf.value(p).abs() > 1e-9 * reach.max(1e-300) * 100.0 || !f.contains(p, 1e-6 * reach)) {
        return None;
    }
    if len3([p[0] - v[0], p[1] - v[1], p[2] - v[2]]) > 1.6 * reach || field.f(p).abs() > tol_surface {
        return None;
    }
    Some(p)
}

/// Moves `p` onto the zero set of `f` with Newton steps along the numerical gradient.
fn newton(f: &dyn Fn([f64; 3]) -> f64, mut p: [f64; 3], scale: f64) -> [f64; 3] {
    let e = 1e-6 * scale;
    for _ in 0..6 {
        let v = f(p);
        if v.abs() < 1e-12 * scale {
            break;
        }
        let mut g = [0.0; 3];
        for k in 0..3 {
            let (mut a, mut b) = (p, p);
            a[k] -= e;
            b[k] += e;
            g[k] = (f(b) - f(a)) / (2.0 * e);
        }
        let g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        if g2 < 1e-300 {
            break;
        }
        for k in 0..3 {
            p[k] -= v * g[k] / g2;
        }
    }
    p
}

/// Meshes the solid with tetrahedra.
pub fn mesh_3d(spec: &Spec3d) -> Result<(SimplexMesh, Report3d), String> {
    if spec.bodies.is_empty() {
        return Err("add at least one body".into());
    }
    if spec.h.is_nan() || spec.h <= 0.0 || !spec.h.is_finite() {
        return Err("the element size must be positive".into());
    }
    if !(0.05..=0.45).contains(&spec.alpha) {
        return Err("alpha must be between 0.05 and 0.45".into());
    }
    let positive: Vec<usize> = (0..spec.bodies.len()).filter(|&i| spec.bodies[i].domain > 0).collect();
    if positive.is_empty() {
        return Err("at least one body needs a positive domain number (domain 0 or less makes a hole)".into());
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &i in &positive {
        let (l, h) = spec.bodies[i].solid.bbox();
        for k in 0..3 {
            lo[k] = lo[k].min(l[k]);
            hi[k] = hi[k].max(h[k]);
        }
    }
    // lattice spacing per axis: a thin body gets a finer lattice across its thickness
    let mut h_axis = [spec.h; 3];
    if spec.min_layers > 0.0 {
        for &i in &positive {
            let (l, h) = spec.bodies[i].solid.bbox();
            for k in 0..3 {
                let t = h[k] - l[k];
                if t > 0.0 && t < spec.min_layers * spec.h {
                    h_axis[k] = h_axis[k].min(t / spec.min_layers).max(0.25 * spec.h);
                }
            }
        }
    }
    let mut n = [0usize; 3];
    let mut hh = [0.0; 3];
    for k in 0..3 {
        let ext = hi[k] - lo[k];
        if ext.is_nan() || ext <= 0.0 || !ext.is_finite() {
            return Err("a body has no extent in some direction".into());
        }
        let cells = (ext / h_axis[k]).round().max(1.0) as usize;
        hh[k] = ext / cells as f64;
        n[k] = cells + 2; // one padding cell on each side keeps the lattice boundary outside the solid
    }
    let total_cells = n[0] * n[1] * n[2];
    if total_cells * 12 > 60_000_000 {
        return Err(format!("that would be about {} million tetrahedra; choose a larger element size (limit 60 million)", total_cells * 12 / 1_000_000));
    }
    let lat = Lattice { n, origin: [lo[0] - hh[0], lo[1] - hh[1], lo[2] - hh[2]], h: hh, n_corners: (n[0] + 1) * (n[1] + 1) * (n[2] + 1) };
    let nv = lat.n_vertices();
    let field = Field { bodies: &spec.bodies };
    // the faces of the positive bodies are cutting surfaces when there are several of them
    let mut faces: Vec<FaceSurf> = vec![];
    let distinct_domains = positive.iter().map(|&i| spec.bodies[i].domain).collect::<std::collections::BTreeSet<_>>().len();
    if distinct_domains >= 2 {
        for &i in &positive {
            faces.extend(spec.bodies[i].solid.faces_with_extent());
        }
    }
    if faces.len() + 1 > MAX_SURFACES {
        return Err(format!("too many body faces for fitted interfaces ({}; the limit is {}): use fewer or simpler bodies", faces.len(), MAX_SURFACES - 1));
    }
    let sf = Surfaces { field, faces, delta: 1e-6 * hh[0].min(hh[1]).min(hh[2]) };
    let ns = sf.n();
    let hmin = hh[0].min(hh[1]).min(hh[2]);
    let tol = 1e-9 * hmin;
    let alpha = spec.alpha;
    let mut tl = std::time::Instant::now();
    let mut lap = |what: &str| crate::solid::lap(&mut tl, what);

    // 1. classify the lattice vertices for every cutting surface: -1 inside, +1 outside, 0 on it
    let mut cls: Vec<Vec<i8>> = Vec::with_capacity(ns);
    for s in 0..ns {
        let blocks = par_blocks(nv, 20_000, |r| r.map(|v| sf.class(s, lat.pos(v as u32), tol)).collect::<Vec<i8>>());
        let mut c = Vec::with_capacity(nv);
        for b in blocks {
            c.extend(b);
        }
        cls.push(c);
    }
    lap("  mesh3d: classify vertices");

    // faces of all bodies, for fitting sharp edges and corners
    let all_faces: Vec<FaceSurf> = if spec.snap_features { spec.bodies.iter().flat_map(|b| b.solid.faces_with_extent()).collect() } else { vec![] };

    // 2. warp vertices that are within alpha of an edge length of a surface onto it; vertices near a
    // crease or corner of the solid go to it
    let warps = par_blocks(nv, 20_000, |r| {
        let mut out: Vec<(u32, [f64; 3], usize)> = vec![];
        let mut nb = Vec::with_capacity(14);
        for v in r {
            let pv = lat.pos(v as u32);
            // (1) a crease or corner of the solid within reach
            if !all_faces.is_empty() && (0..ns).all(|s| s == 0 || cls[s][v] != 0) {
                let reach = alpha * hmin;
                // coplanar faces (of touching bodies) count once, with their extents joined
                let mut near: Vec<FaceSurf> = vec![];
                for f in all_faces.iter().filter(|f| f.surf.value(pv).abs() < reach) {
                    if let Some(g) = near.iter_mut().find(|g| g.same_plane(f)) {
                        for k in 0..3 {
                            g.lo[k] = g.lo[k].min(f.lo[k]);
                            g.hi[k] = g.hi[k].max(f.hi[k]);
                        }
                    } else {
                        near.push(f.clone());
                    }
                }
                if near.len() >= 2 {
                    {
                        // at most three, the nearest
                        let mut sorted = near;
                        sorted.sort_by(|a, b| a.surf.value(pv).abs().partial_cmp(&b.surf.value(pv).abs()).unwrap());
                        let pick: Vec<&FaceSurf> = sorted.iter().take(3).collect();
                        let mut found = None;
                        for take in [pick.len(), 2] {
                            if take < 2 || take > pick.len() {
                                continue;
                            }
                            if let Some(q) = snap_to_crease(&pick[..take], &sf.field, pv, reach, 1e-6 * hmin) {
                                found = Some(q);
                                break;
                            }
                        }
                        if let Some(q) = found {
                            // a vertex that already sits on the crease stays
                            if len3([q[0] - pv[0], q[1] - pv[1], q[2] - pv[2]]) > 1e-9 * hmin {
                                out.push((v as u32, q, 0));
                            }
                            continue;
                        }
                    }
                }
            }
            // a vertex that is already on some surface stays where it is
            if (0..ns).any(|s| cls[s][v] == 0) {
                continue;
            }
            lat.neighbours(v as u32, &mut nb);
            // (2) nearest cut point on an edge to an opposite-sign vertex, over all surfaces, and
            // (preferred, for the solid's own surface) the point where an edge that runs inside the
            // surface (to a vertex exactly on it) leaves the solid: moving there keeps a vertex on a flat face
            let mut best: Option<(f64, [f64; 3], usize)> = None;
            let mut in_surface: Option<(f64, [f64; 3], usize)> = None;
            for &w in &nb {
                let pw = lat.pos(w);
                let elen = len3([pw[0] - pv[0], pw[1] - pv[1], pw[2] - pv[2]]);
                for s in 0..ns {
                    let (cv, cw) = (cls[s][v], cls[s][w as usize]);
                    if cw == -cv {
                        let pa = [pv[0] + alpha * (pw[0] - pv[0]), pv[1] + alpha * (pw[1] - pv[1]), pv[2] + alpha * (pw[2] - pv[2])];
                        let sa = sf.value(s, pa);
                        if (sa < 0.0) != (cv < 0) || sa.abs() <= tol {
                            let (t, c) = bisect_surface(&sf, s, pv, pw);
                            if t < alpha {
                                let d = t * elen;
                                if best.is_none_or(|(bd, _, _)| d < bd) {
                                    best = Some((d, c, s));
                                }
                            }
                        }
                    } else if s == 0 && cw == 0 && cv > 0 {
                        let pa = [pv[0] + alpha * (pw[0] - pv[0]), pv[1] + alpha * (pw[1] - pv[1]), pv[2] + alpha * (pw[2] - pv[2])];
                        if sf.value(0, pa) <= tol {
                            // leaving the surface region along the edge: bisect on f <= tol
                            let (mut lo_t, mut hi_t) = (0.0f64, 1.0f64);
                            let at = |t: f64| [pv[0] + t * (pw[0] - pv[0]), pv[1] + t * (pw[1] - pv[1]), pv[2] + t * (pw[2] - pv[2])];
                            for _ in 0..44 {
                                let mid = 0.5 * (lo_t + hi_t);
                                if sf.value(0, at(mid)) <= tol {
                                    hi_t = mid;
                                } else {
                                    lo_t = mid;
                                }
                            }
                            let t = 0.5 * (lo_t + hi_t);
                            if t < alpha {
                                let d = t * elen;
                                if in_surface.is_none_or(|(bd, _, _)| d < bd) {
                                    in_surface = Some((d, at(t), 0));
                                }
                            }
                        }
                    }
                }
            }
            if let Some((_, c, s)) = in_surface.or(best) {
                out.push((v as u32, c, s));
            }
        }
        out
    });
    let mut warped: HashMap<u32, [f64; 3]> = HashMap::new();
    let mut warped_to: Vec<(u32, usize)> = vec![];
    for w in warps {
        for (v, c, st) in w {
            warped.insert(v, c);
            warped_to.push((v, st));
        }
    }
    // classes at the new positions: zero for the surface a vertex was moved onto, by position for the others
    for &(v, st) in &warped_to {
        let p = warped[&v];
        for s in 0..ns {
            cls[s][v as usize] = if s == st { 0 } else { sf.class(s, p, tol) };
        }
    }
    // a warp that puts all four vertices of a lattice tetrahedron on one surface would leave a flat
    // tetrahedron: undo the largest of the warps in it, and repeat
    for _round in 0..6 {
        let flat = par_blocks(n[0] * n[1] * n[2], 4_000, |range| {
            let mut undo: Vec<u32> = vec![];
            for c in range {
                let (i, j, k) = (c % n[0], (c / n[0]) % n[1], c / (n[0] * n[1]));
                let a = lat.center(i, j, k);
                for axis in 0..3 {
                    let (b, corners): (u32, [u32; 4]) = match axis {
                        0 if i + 1 < n[0] => (lat.center(i + 1, j, k), [lat.corner(i + 1, j, k), lat.corner(i + 1, j + 1, k), lat.corner(i + 1, j + 1, k + 1), lat.corner(i + 1, j, k + 1)]),
                        1 if j + 1 < n[1] => (lat.center(i, j + 1, k), [lat.corner(i, j + 1, k), lat.corner(i + 1, j + 1, k), lat.corner(i + 1, j + 1, k + 1), lat.corner(i, j + 1, k + 1)]),
                        2 if k + 1 < n[2] => (lat.center(i, j, k + 1), [lat.corner(i, j, k + 1), lat.corner(i + 1, j, k + 1), lat.corner(i + 1, j + 1, k + 1), lat.corner(i, j + 1, k + 1)]),
                        _ => continue,
                    };
                    for e in 0..4 {
                        let v = [a, b, corners[e], corners[(e + 1) % 4]];
                        let pp: Vec<[f64; 3]> = v.iter().map(|&w| warped.get(&w).copied().unwrap_or_else(|| lat.pos(w))).collect();
                        let (e1, e2, e3) = ([pp[1][0] - pp[0][0], pp[1][1] - pp[0][1], pp[1][2] - pp[0][2]], [pp[2][0] - pp[0][0], pp[2][1] - pp[0][1], pp[2][2] - pp[0][2]], [pp[3][0] - pp[0][0], pp[3][1] - pp[0][1], pp[3][2] - pp[0][2]]);
                        let vol = (e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])).abs() / 6.0;
                        if vol < 2e-3 * hmin * hmin * hmin {
                            // the warped vertex that moved furthest
                            let mut best: Option<(f64, u32)> = None;
                            for &w in &v {
                                if let Some(p) = warped.get(&w) {
                                    let q = lat.pos(w);
                                    let d = len3([p[0] - q[0], p[1] - q[1], p[2] - q[2]]);
                                    if best.is_none_or(|(bd, _)| d > bd) {
                                        best = Some((d, w));
                                    }
                                }
                            }
                            if let Some((_, w)) = best {
                                undo.push(w);
                            }
                        }
                    }
                }
            }
            undo
        });
        let mut undo: Vec<u32> = flat.into_iter().flatten().collect();
        undo.sort_unstable();
        undo.dedup();
        if undo.is_empty() {
            break;
        }
        for w in undo {
            warped.remove(&w);
            let p = lat.pos(w);
            for s in 0..ns {
                cls[s][w as usize] = sf.class(s, p, tol);
            }
        }
    }
    let n_warped = warped.len();
    lap("  mesh3d: warp");
    let pos_of = |v: u32| -> [f64; 3] { warped.get(&v).copied().unwrap_or_else(|| lat.pos(v)) };

    // 3. generate tetrahedra: whole ones in the interior, cut pieces next to the surfaces
    struct Part {
        inner: Vec<[u32; 4]>,
        surf: Vec<[u64; 4]>,
        /// cut points made here: id, position, on the outer surface
        cuts: Vec<(u64, [f64; 3], bool)>,
    }
    let make_v = |v: u32| -> V {
        let mut cl = [0i8; MAX_SURFACES];
        for s in 0..ns {
            cl[s] = cls[s][v as usize];
        }
        V { id: v as u64, p: pos_of(v), surf: cls[0][v as usize] == 0, cl }
    };
    let parts = par_blocks(n[0] * n[1] * n[2], 4_000, |range| {
        let mut part = Part { inner: vec![], surf: vec![], cuts: vec![] };
        let mut pieces: Vec<[V; 4]> = vec![];
        let mut next: Vec<[V; 4]> = vec![];
        for c in range {
            let (i, j, k) = (c % n[0], (c / n[0]) % n[1], c / (n[0] * n[1]));
            let a = lat.center(i, j, k);
            for axis in 0..3 {
                let (b, corners): (u32, [u32; 4]) = match axis {
                    0 if i + 1 < n[0] => (lat.center(i + 1, j, k), [lat.corner(i + 1, j, k), lat.corner(i + 1, j + 1, k), lat.corner(i + 1, j + 1, k + 1), lat.corner(i + 1, j, k + 1)]),
                    1 if j + 1 < n[1] => (lat.center(i, j + 1, k), [lat.corner(i, j + 1, k), lat.corner(i + 1, j + 1, k), lat.corner(i + 1, j + 1, k + 1), lat.corner(i, j + 1, k + 1)]),
                    2 if k + 1 < n[2] => (lat.center(i, j, k + 1), [lat.corner(i, j, k + 1), lat.corner(i + 1, j, k + 1), lat.corner(i + 1, j + 1, k + 1), lat.corner(i, j + 1, k + 1)]),
                    _ => continue,
                };
                for e in 0..4 {
                    let v = [a, b, corners[e], corners[(e + 1) % 4]];
                    // the solid's surface: nothing inside, nothing to cut, or cut
                    let c0 = [cls[0][v[0] as usize], cls[0][v[1] as usize], cls[0][v[2] as usize], cls[0][v[3] as usize]];
                    if c0.iter().all(|&x| x >= 0) {
                        continue;
                    }
                    let straddles0 = c0.iter().any(|&x| x > 0);
                    // the interfaces between materials: does any surface cross this tetrahedron?
                    let mut straddles = straddles0;
                    for s in 1..ns {
                        let (mut ni, mut no) = (0, 0);
                        for &w in &v {
                            match cls[s][w as usize] {
                                -1 => ni += 1,
                                1 => no += 1,
                                _ => {}
                            }
                        }
                        if ni > 0 && no > 0 {
                            straddles = true;
                        }
                    }
                    if !straddles {
                        part.inner.push(v);
                        continue;
                    }
                    pieces.clear();
                    pieces.push([make_v(v[0]), make_v(v[1]), make_v(v[2]), make_v(v[3])]);
                    for s in 0..ns {
                        next.clear();
                        for t in pieces.iter() {
                            let sides: &[i8] = if s == 0 { &[-1] } else { &[-1, 1] };
                            for &side in sides {
                                clip_piece(&sf, s, side, t, tol, &mut part.cuts, &mut next);
                            }
                        }
                        std::mem::swap(&mut pieces, &mut next);
                        if pieces.is_empty() {
                            break;
                        }
                    }
                    for t in &pieces {
                        if t.iter().all(|w| w.id & TAG == 0) {
                            part.inner.push([t[0].id as u32, t[1].id as u32, t[2].id as u32, t[3].id as u32]);
                        } else {
                            part.surf.push([t[0].id, t[1].id, t[2].id, t[3].id]);
                        }
                    }
                }
            }
        }
        part
    });
    lap("  mesh3d: tetrahedra");

    // 4. number the vertices: used lattice vertices, then the cut points
    let mut used = vec![false; nv];
    let mut n_inner = 0;
    let mut cut_list: Vec<(u64, [f64; 3], bool)> = vec![];
    for p in &parts {
        n_inner += p.inner.len();
        for t in &p.inner {
            for &v in t {
                used[v as usize] = true;
            }
        }
        for t in &p.surf {
            for &k in t {
                if k & TAG == 0 {
                    used[k as usize] = true;
                }
            }
        }
        cut_list.extend(p.cuts.iter().copied());
    }
    cut_list.sort_unstable_by_key(|c| c.0);
    cut_list.dedup_by_key(|c| c.0);
    let n_cuts = cut_list.len();
    let mut remap = vec![u32::MAX; nv];
    let mut points: Vec<f64> = Vec::new();
    let mut vflag: Vec<bool> = vec![];
    let mut next_id = 0u32;
    for v in 0..nv {
        if used[v] {
            remap[v] = next_id;
            next_id += 1;
            points.extend_from_slice(&pos_of(v as u32));
            vflag.push(cls[0][v] == 0);
        }
    }
    let n_lattice_used = next_id as usize;
    for (_, p, f) in &cut_list {
        points.extend_from_slice(p);
        vflag.push(*f);
    }
    let vid = |k: u64| -> u32 {
        if k & TAG != 0 {
            n_lattice_used as u32 + cut_list.binary_search_by_key(&k, |c| c.0).expect("cut point present") as u32
        } else {
            remap[k as usize]
        }
    };
    lap("  mesh3d: numbering, cut points");

    // 5. cells (oriented), domain numbers, boundary facets and their numbers
    let mut cells: Vec<u32> = Vec::with_capacity((n_inner + parts.iter().map(|p| p.surf.len()).sum::<usize>()) * 4);
    for p in &parts {
        for t in &p.inner {
            cells.extend(t.iter().map(|&v| remap[v as usize]));
        }
        for t in &p.surf {
            cells.extend(t.iter().map(|&k| vid(k)));
        }
    }
    let pt = |v: u32| -> [f64; 3] { [points[3 * v as usize], points[3 * v as usize + 1], points[3 * v as usize + 2]] };
    let nc = cells.len() / 4;
    // faces on the outer surface: three vertices on it, used by one cell
    let mut keyed: Vec<[u32; 3]> = vec![];
    for c in 0..nc {
        let t = [cells[4 * c], cells[4 * c + 1], cells[4 * c + 2], cells[4 * c + 3]];
        if t.iter().filter(|&&v| vflag[v as usize]).count() >= 3 {
            for skip in 0..4 {
                let f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| t[i]).collect();
                if f.iter().all(|&v| vflag[v as usize]) {
                    let mut k = [f[0], f[1], f[2]];
                    k.sort_unstable();
                    keyed.push(k);
                }
            }
        }
    }
    keyed.sort_unstable();
    let mut facets: Vec<u32> = vec![];
    let mut i = 0;
    while i < keyed.len() {
        let mut j = i + 1;
        while j < keyed.len() && keyed[j] == keyed[i] {
            j += 1;
        }
        if j - i == 1 {
            facets.extend_from_slice(&keyed[i]);
        }
        i = j;
    }
    let field = &sf.field;
    let blocks = par_blocks(nc, 20_000, |r| {
        let mut fixed: Vec<[u32; 4]> = Vec::with_capacity(r.len());
        let mut labels: Vec<i32> = Vec::with_capacity(r.len());
        for c in r {
            let mut t = [cells[4 * c], cells[4 * c + 1], cells[4 * c + 2], cells[4 * c + 3]];
            let (p0, p1, p2, p3) = (pt(t[0]), pt(t[1]), pt(t[2]), pt(t[3]));
            let (a, b, d) = ([p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]], [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]], [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]]);
            let det = a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) + a[2] * (b[0] * d[1] - b[1] * d[0]);
            if det < 0.0 {
                t.swap(1, 2);
            }
            fixed.push(t);
            let cen = [(p0[0] + p1[0] + p2[0] + p3[0]) / 4.0, (p0[1] + p1[1] + p2[1] + p3[1]) / 4.0, (p0[2] + p1[2] + p2[2] + p3[2]) / 4.0];
            labels.push(field.domain_near(cen));
        }
        (fixed, labels)
    });
    let mut cells: Vec<u32> = Vec::with_capacity(nc * 4);
    let mut cell_labels: Vec<i32> = Vec::with_capacity(nc);
    for (f, l) in blocks {
        for t in f {
            cells.extend_from_slice(&t);
        }
        cell_labels.extend(l);
    }
    let nf = facets.len() / 3;
    let flabels: Vec<i32> = par_blocks(nf, 10_000, |r| r.map(|f| {
        let (p0, p1, p2) = (pt(facets[3 * f]), pt(facets[3 * f + 1]), pt(facets[3 * f + 2]));
        field.boundary_label([(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0])
    }).collect::<Vec<i32>>()).into_iter().flatten().collect();
    lap("  mesh3d: cells, labels, facets");
    let mut report = Report3d { lattice_vertices: nv, warped: n_warped, cut_points: n_cuts, cells: nc, refined: 0 };
    let mesh = SimplexMesh::new(3, points, cells, cell_labels, facets, flabels);
    lap("  mesh3d: SimplexMesh::new (edges)");
    // 4. local refinement inside bodies that ask for a finer element size
    let local: Vec<(usize, f64)> = spec.bodies.iter().enumerate().filter(|(_, b)| b.domain > 0).filter_map(|(i, b)| b.h.filter(|&hb| hb > 0.0 && hb < 0.999 * spec.h).map(|hb| (i, hb))).collect();
    if local.is_empty() {
        return Ok((mesh, report));
    }
    let bodies = &spec.bodies;
    let size = |p: [f64; 3]| -> f64 {
        let mut s = spec.h;
        for &(i, hb) in &local {
            if bodies[i].solid.sdf(p) <= 0.0 {
                s = s.min(hb);
            }
        }
        s
    };
    let field = Field { bodies };
    let project = |mid: [f64; 3], pa: [f64; 3], pb: [f64; 3], kind: crate::refine3d::EdgeKind| -> [f64; 3] {
        let len = len3([pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]]);
        let mut p = mid;
        if kind == crate::refine3d::EdgeKind::Boundary {
            p = newton(&|q| field.f(q), p, len);
        }
        // a curved face of a body that carries both ends of the edge
        for b in bodies {
            if b.solid.sdf(pa).abs() < 1e-7 * len && b.solid.sdf(pb).abs() < 1e-7 * len && b.solid.sdf(p).abs() > 1e-9 * len {
                p = newton(&|q| b.solid.sdf(q), p, len);
            }
        }
        if len3([p[0] - mid[0], p[1] - mid[1], p[2] - mid[2]]) > 0.5 * len {
            mid
        } else {
            p
        }
    };
    let (mesh, rr) = crate::refine3d::refine_to_size(&mesh, &size, &project, 1.3, 60_000_000 / 12)?;
    report.cells = rr.cells_after;
    report.refined = rr.bisections;
    lap("  mesh3d: local refinement");
    Ok((mesh, report))
}

/// Cuts a tetrahedron (given as four vertices) along surface `s` and keeps the part on `side`
/// (-1: inside the surface, +1: outside); vertices on the surface count for both sides.
fn clip_piece(sf: &Surfaces, s: usize, side: i8, t: &[V; 4], tol: f64, cuts: &mut Vec<(u64, [f64; 3], bool)>, out: &mut Vec<[V; 4]>) {
    let cl: Vec<i8> = t.iter().map(|w| w.cl[s]).collect();
    let (mut ins, mut outs, mut ons): (Vec<V>, Vec<V>, Vec<V>) = (vec![], vec![], vec![]);
    for k in 0..4 {
        if cl[k] == side {
            ins.push(t[k]);
        } else if cl[k] == -side {
            outs.push(t[k]);
        } else {
            ons.push(t[k]);
        }
    }
    if ins.is_empty() {
        return;
    }
    if outs.is_empty() {
        out.push(*t);
        return;
    }
    let mut cut = |a: V, b: V| -> V {
        let (lo, hi) = if a.id < b.id { (a, b) } else { (b, a) };
        let id = cut_id(lo.id, hi.id, s);
        let (_, p) = bisect_surface(sf, s, lo.p, hi.p);
        let on_outer = if s == 0 { true } else { a.surf && b.surf };
        cuts.push((id, p, on_outer));
        let mut cl = [0i8; MAX_SURFACES];
        for s2 in 0..sf.n() {
            cl[s2] = if s2 == s { 0 } else { sf.class(s2, p, tol) };
        }
        if s != 0 && on_outer {
            cl[0] = 0;
        }
        V { id, p, surf: on_outer, cl }
    };
    let mut tmp: Vec<[V; 4]> = vec![];
    match (ins.len(), outs.len()) {
        (1, 1) => {
            let c = cut(ins[0], outs[0]);
            tmp.push([ins[0], c, ons[0], ons[1]]);
        }
        (1, 2) => {
            let (c1, c2) = (cut(ins[0], outs[0]), cut(ins[0], outs[1]));
            tmp.push([ins[0], c1, c2, ons[0]]);
        }
        (1, 3) => {
            let (c1, c2, c3) = (cut(ins[0], outs[0]), cut(ins[0], outs[1]), cut(ins[0], outs[2]));
            tmp.push([ins[0], c1, c2, c3]);
        }
        (2, 1) => {
            let (c1, c2) = (cut(ins[0], outs[0]), cut(ins[1], outs[0]));
            let mut qt = vec![];
            quad_tris([ins[0], ins[1], c2, c1], &mut qt);
            for q in qt {
                tmp.push([ons[0], q[0], q[1], q[2]]);
            }
        }
        (2, 2) => {
            let (c11, c12, c21, c22) = (cut(ins[0], outs[0]), cut(ins[0], outs[1]), cut(ins[1], outs[0]), cut(ins[1], outs[1]));
            split_prism([ins[0], c11, c12, ins[1], c21, c22], &mut tmp);
        }
        (3, 1) => {
            let (c1, c2, c3) = (cut(ins[0], outs[0]), cut(ins[1], outs[0]), cut(ins[2], outs[0]));
            split_prism([ins[0], ins[1], ins[2], c1, c2, c3], &mut tmp);
        }
        _ => unreachable!(),
    }
    out.extend(tmp);
}
