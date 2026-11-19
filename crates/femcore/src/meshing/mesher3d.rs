//! Tetrahedral meshes of 3D models: a lattice fill that is fitted to the exact geometry of any model that can say which region a point is in.
//!
//! The geometry is reached through the small trait [`Domain3`] (`region_at`, region boxes, labels, a projection onto the boundary). The geom
//! adapter (`meshing::geomdrive3`) implements it with `Evaluated::classify` and `Complex::facet_project`, so the mesher never meets a
//! primitive: the same code meshes a sphere in a box, a torus, a layer stack or whatever geom can evaluate.
//!
//! The method is isosurface stuffing (Labelle and Shewchuk 2007, with the vertex warping and marching-tetrahedra cuts of Quartet): the
//! model is filled with the tetrahedra of a body-centred cubic lattice, which have good shape; a lattice vertex that lies within `alpha`
//! of an edge length of a region boundary is moved onto it, and the few tetrahedra that a boundary still crosses are cut along it. Every
//! step is local and independent, so the vertices are classified in parallel, tetrahedra are generated in parallel blocks, only the thin layer
//! next to the boundaries needs any bookkeeping, and cut points are identified by the pair of lattice vertices of their edge, so no search
//! structure is needed.
//!
//! The lattice is aligned with the bounding box of the meshed regions, so a face of the model that lies on the box is reproduced exactly;
//! curved surfaces are piecewise linear with a geometric error of order h^2 (ADR 18: no curved elements). Sharp creases that do not lie on a
//! lattice plane are *not* fitted in this version (the vertices near them are warped onto one face and the other face cuts the rest), so
//! the error there is of order `alpha h`.
//!
//! What the geometry has to provide is only the region at a point: the surfaces that are fitted are the boundaries of the regions
//! (membership, found by bisection on a lattice edge), and the solid, the union of the regions. Where two surfaces coincide (an interface
//! is the boundary of both of its regions) the points on it are classed as lying on both, from the two ends of the bisection bracket.
//!
//! Labels: a cell takes the domain number of its region (`Domain3::region_label`); the boundary facets are the faces of cells that have
//! the void on the other side and the faces between cells of different regions, and take the label of the geom face (`boundary_label`)
//! they lie on. Both are 1-based (ADR 32).

use crate::mesh::cells::{CellMesh, CellType};
use std::collections::HashMap;

/// The geometry the mesher fits to. Regions are numbered `0..n_regions()`; a point outside every region is in the void.
pub trait Domain3: Sync {
    /// Region at `p` (`None`: the void, which includes holes).
    fn region_at(&self, p: [f64; 3]) -> Option<u32>;
    fn n_regions(&self) -> usize;
    /// Axis-aligned box of a region (lower and upper corner).
    fn region_box(&self, r: u32) -> ([f64; 3], [f64; 3]);
    /// Element size the region asks for (finer than the global one means local refinement); `None`: the global size.
    fn region_size(&self, r: u32) -> Option<f64>;
    /// Domain number (1-based) of the cells of a region.
    fn region_label(&self, r: u32) -> i32;
    /// Boundary number (1-based) of the model face that a mesh facet with centroid `p` lies on, between the regions `a` and `b`
    /// (`None`: the void).
    fn boundary_label(&self, p: [f64; 3], a: Option<u32>, b: Option<u32>) -> i32;
    /// The position of a new vertex on the mesh edge `pa`-`pb`, both of which lie on a model face, moved onto that face (`None` when the
    /// two ends lie on no common face).
    fn project_edge(&self, mid: [f64; 3], pa: [f64; 3], pb: [f64; 3]) -> Option<[f64; 3]>;
    /// Whether the edge `pa`-`pb` lies along a crease: both ends lie on two different model faces that are not both flat. Such an edge is
    /// left straight by `curve3` (ADR 82); the default says no.
    fn edge_on_crease(&self, _pa: [f64; 3], _pb: [f64; 3]) -> bool {
        false
    }
    /// The element size a model face asks for around `p` (`None`: no face asks, the sizes of the regions apply). ADR 83.
    fn size_at(&self, _p: [f64; 3]) -> Option<f64> {
        None
    }
    /// Whether any face asks for a size (`size_at` can answer `Some`); the local refinement then runs even where no region is finer.
    fn asks_size_at(&self) -> bool {
        false
    }
}

/// Settings of the 3D mesher.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec3d {
    /// target element (lattice) size
    pub h: f64,
    /// a lattice vertex within this fraction of an edge length of a boundary is moved onto it (0.3 gives dihedral angles of at least about
    /// 10 degrees on smooth surfaces)
    pub alpha: f64,
    /// a region thinner than this many elements along an axis gets a finer lattice along that axis (at most 4 times finer than `h`), so
    /// that it holds at least this many layers; 0 turns it off
    pub min_layers: f64,
}

impl Spec3d {
    pub fn new(h: f64) -> Self {
        Spec3d { h, alpha: 0.3, min_layers: 2.0 }
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
    /// what the curving of the edges did (ADR 82); `None` when the mesh is straight-sided
    pub curve: Option<crate::meshing::curve3::CurveReport>,
}

const TAG: u64 = 1 << 63;

/// Largest number of cutting surfaces: the solid and the regions.
pub const MAX_SURFACES: usize = 64;

/// A vertex of a piece during cutting: its id (lattice number, or a tagged hash for a cut point), its position, its class with respect
/// to every cutting surface (-1 inside, 0 on, +1 outside) and the region it certainly lies in (`-2`: on a surface, not known).
#[derive(Clone, Copy, Debug)]
struct V {
    id: u64,
    p: [f64; 3],
    cl: [i8; MAX_SURFACES],
    reg: i32,
}

fn len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn det3(a: [f64; 3], b: [f64; 3], d: [f64; 3]) -> f64 {
    a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) + a[2] * (b[0] * d[1] - b[1] * d[0])
}

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

/// Splits a prism (triangles `v[0..3]` and `v[3..6]`, vertical edges v[i]-v[i+3]) into three tetrahedra by coning from its
/// lowest-numbered vertex; quadrilateral faces are cut along the diagonal from their lowest-numbered vertex so neighbours agree.
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

/// Runs `f(range)` on blocks of `0..n` in parallel and returns the results in order.
fn par_blocks<T: Send, F: Fn(std::ops::Range<usize>) -> T + Sync>(n: usize, min_block: usize, f: F) -> Vec<T> {
    let nt = crate::par::threads().min(n.div_ceil(min_block.max(1))).max(1);
    let chunk = n.div_ceil(nt).max(1);
    use rayon::prelude::*;
    (0..n.div_ceil(chunk)).into_par_iter().map(|t| f(t * chunk..((t + 1) * chunk).min(n))).collect()
}

/// The cutting surfaces: 0 is the boundary of the solid (the union of the regions); surface `1 + r` is the boundary of region `r` (only
/// when there are several regions). A point is inside a surface (-1) or outside it (+1), decided by the region at the point alone.
struct Surfaces<'a> {
    dom: &'a dyn Domain3,
    ns: usize,
}

/// A point found by bisection on a segment: the middle of the final bracket, and the two ends of the bracket (on either side).
struct Cut {
    t: f64,
    p: [f64; 3],
    inside: [f64; 3],
    outside: [f64; 3],
}

impl Surfaces<'_> {
    fn sign(&self, s: usize, region: i32) -> i8 {
        if s == 0 {
            if region >= 0 {
                -1
            } else {
                1
            }
        } else if region == (s - 1) as i32 {
            -1
        } else {
            1
        }
    }
    fn region(&self, p: [f64; 3]) -> i32 {
        self.dom.region_at(p).map_or(-1, |r| r as i32)
    }
    fn class(&self, s: usize, p: [f64; 3]) -> i8 {
        self.sign(s, self.region(p))
    }
    /// Point on the segment `p0`-`p1` where the side of surface `s` changes (bisection); `inside` is the end of the final bracket on the
    /// side of `p0`.
    fn bisect(&self, s: usize, p0: [f64; 3], p1: [f64; 3]) -> Cut {
        let c0 = self.class(s, p0);
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        let at = |t: f64| [p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1]), p0[2] + t * (p1[2] - p0[2])];
        for _ in 0..44 {
            let mid = 0.5 * (lo + hi);
            if self.class(s, at(mid)) == c0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t = 0.5 * (lo + hi);
        Cut { t, p: at(t), inside: at(lo), outside: at(hi) }
    }
    /// Classes of a point found on a bracket: zero for every surface whose side changes across it (the surface it was found on, and any
    /// other that coincides there, such as the second boundary of an interface), the side for the others.
    fn cut_classes(&self, c: &Cut) -> [i8; MAX_SURFACES] {
        let (ra, rb) = (self.region(c.inside), self.region(c.outside));
        let mut cl = [0i8; MAX_SURFACES];
        for (s, slot) in cl.iter_mut().enumerate().take(self.ns) {
            let (a, b) = (self.sign(s, ra), self.sign(s, rb));
            *slot = if a == b { a } else { 0 };
        }
        cl
    }
}

/// Meshes the model with tetrahedra. Cells carry `region_label`, boundary facets (outer faces and faces between regions) carry
/// `boundary_label`.
pub fn mesh_3d(dom: &dyn Domain3, spec: &Spec3d) -> Result<(CellMesh<3>, Report3d), String> {
    if !(spec.h > 0.0 && spec.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    if !(0.05..=0.45).contains(&spec.alpha) {
        return Err("alpha must be between 0.05 and 0.45".into());
    }
    let nreg = dom.n_regions();
    if nreg == 0 {
        return Err("the model has no region to mesh".into());
    }
    if nreg + 1 > MAX_SURFACES {
        return Err(format!("too many regions for fitted interfaces ({nreg}; the limit is {}): use fewer regions", MAX_SURFACES - 1));
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for r in 0..nreg as u32 {
        let (l, h) = dom.region_box(r);
        for k in 0..3 {
            lo[k] = lo[k].min(l[k]);
            hi[k] = hi[k].max(h[k]);
        }
    }
    // lattice spacing per axis: a thin region gets a finer lattice across its thickness
    let mut h_axis = [spec.h; 3];
    if spec.min_layers > 0.0 {
        for r in 0..nreg as u32 {
            let (l, h) = dom.region_box(r);
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
            return Err("the model has no extent in some direction".into());
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
    let ns = if nreg >= 2 { 1 + nreg } else { 1 };
    let sf = Surfaces { dom, ns };
    let hmin = hh[0].min(hh[1]).min(hh[2]);
    let alpha = spec.alpha;
    let mut tl = std::time::Instant::now();
    let mut lap = |what: &str| {
        crate::progress::lap(what, tl.elapsed().as_secs_f64());
        tl = std::time::Instant::now();
    };

    // 1. the region at every lattice vertex, and from it the class for every cutting surface: -1 inside, +1 outside
    let reg_lat: Vec<i32> = par_blocks(nv, 20_000, |r| r.map(|v| sf.region(lat.pos(v as u32))).collect::<Vec<i32>>()).into_iter().flatten().collect();
    let mut cls: Vec<Vec<i8>> = (0..ns).map(|s| reg_lat.iter().map(|&r| sf.sign(s, r)).collect()).collect();
    lap("  mesh3d: classify vertices");
    crate::progress::check()?;

    // 2. warp vertices that are within alpha of an edge length of a surface onto it
    let warps = par_blocks(nv, 20_000, |r| {
        let mut out: Vec<(u32, [f64; 3], [i8; MAX_SURFACES])> = vec![];
        let mut nb = Vec::with_capacity(14);
        for v in r {
            let pv = lat.pos(v as u32);
            lat.neighbours(v as u32, &mut nb);
            // nearest cut point on an edge to an opposite-side vertex, over all surfaces
            let mut best: Option<(f64, Cut)> = None;
            for &w in &nb {
                let cw_all = |s: usize| cls[s][w as usize];
                if (0..ns).all(|s| cw_all(s) == cls[s][v]) {
                    continue;
                }
                let pw = lat.pos(w);
                let elen = len3(sub3(pw, pv));
                for s in 0..ns {
                    let (cv, cw) = (cls[s][v], cw_all(s));
                    if cw != cv {
                        let pa = [pv[0] + alpha * (pw[0] - pv[0]), pv[1] + alpha * (pw[1] - pv[1]), pv[2] + alpha * (pw[2] - pv[2])];
                        if sf.class(s, pa) != cv {
                            let c = sf.bisect(s, pv, pw);
                            if c.t < alpha {
                                let d = c.t * elen;
                                if best.as_ref().is_none_or(|(bd, _)| d < *bd) {
                                    best = Some((d, c));
                                }
                            }
                        }
                    }
                }
            }
            if let Some((_, c)) = best {
                out.push((v as u32, c.p, sf.cut_classes(&c)));
            }
        }
        out
    });
    let mut warped: HashMap<u32, [f64; 3]> = HashMap::new();
    for w in warps {
        for (v, c, cl) in w {
            warped.insert(v, c);
            for (s, col) in cls.iter_mut().enumerate() {
                col[v as usize] = cl[s];
            }
        }
    }
    // a warp that puts all four vertices of a lattice tetrahedron on one surface would leave a flat tetrahedron: undo the largest of the
    // warps in it, and repeat
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
                        if !v.iter().any(|w| warped.contains_key(w)) {
                            continue;
                        }
                        let pp: Vec<[f64; 3]> = v.iter().map(|&w| warped.get(&w).copied().unwrap_or_else(|| lat.pos(w))).collect();
                        let vol = det3(sub3(pp[1], pp[0]), sub3(pp[2], pp[0]), sub3(pp[3], pp[0])).abs() / 6.0;
                        if vol < 2e-3 * hmin * hmin * hmin {
                            // the warped vertex that moved furthest
                            let mut best: Option<(f64, u32)> = None;
                            for &w in &v {
                                if let Some(p) = warped.get(&w) {
                                    let d = len3(sub3(*p, lat.pos(w)));
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
            for (s, col) in cls.iter_mut().enumerate() {
                col[w as usize] = sf.sign(s, reg_lat[w as usize]);
            }
        }
    }
    let n_warped = warped.len();
    lap("  mesh3d: warp");
    crate::progress::check()?;
    let pos_of = |v: u32| -> [f64; 3] { warped.get(&v).copied().unwrap_or_else(|| lat.pos(v)) };
    // a lattice vertex that was not moved lies strictly inside the region it was classified in
    let known_region = |v: u32| -> i32 {
        if warped.contains_key(&v) {
            -2
        } else {
            reg_lat[v as usize]
        }
    };

    // 3. generate tetrahedra: whole ones in the interior, cut pieces next to the surfaces
    struct Part {
        inner: Vec<([u32; 4], i32)>,
        surf: Vec<([u64; 4], i32)>,
        /// cut points made here: id and position
        cuts: Vec<(u64, [f64; 3])>,
    }
    let make_v = |v: u32| -> V {
        let mut cl = [0i8; MAX_SURFACES];
        for s in 0..ns {
            cl[s] = cls[s][v as usize];
        }
        V { id: v as u64, p: pos_of(v), cl, reg: known_region(v) }
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
                    // the solid: nothing inside, nothing to cut, or cut
                    let c0 = [cls[0][v[0] as usize], cls[0][v[1] as usize], cls[0][v[2] as usize], cls[0][v[3] as usize]];
                    if c0.iter().all(|&x| x >= 0) {
                        continue;
                    }
                    // does any surface cross this tetrahedron?
                    let mut straddles = c0.iter().any(|&x| x > 0);
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
                        let region = v.iter().map(|&w| known_region(w)).find(|&r| r >= 0).unwrap_or(-1);
                        part.inner.push((v, region));
                        continue;
                    }
                    pieces.clear();
                    pieces.push([make_v(v[0]), make_v(v[1]), make_v(v[2]), make_v(v[3])]);
                    for s in 0..ns {
                        next.clear();
                        for t in pieces.iter() {
                            let sides: &[i8] = if s == 0 { &[-1] } else { &[-1, 1] };
                            for &side in sides {
                                clip_piece(&sf, s, side, t, &mut part.cuts, &mut next);
                            }
                        }
                        std::mem::swap(&mut pieces, &mut next);
                        if pieces.is_empty() {
                            break;
                        }
                    }
                    for t in &pieces {
                        let region = t.iter().map(|w| w.reg).find(|&r| r >= 0).unwrap_or(-1);
                        if t.iter().all(|w| w.id & TAG == 0) {
                            part.inner.push(([t[0].id as u32, t[1].id as u32, t[2].id as u32, t[3].id as u32], region));
                        } else {
                            part.surf.push(([t[0].id, t[1].id, t[2].id, t[3].id], region));
                        }
                    }
                }
            }
        }
        part
    });
    lap("  mesh3d: tetrahedra");
    crate::progress::check()?;

    // 4. number the vertices: used lattice vertices, then the cut points
    let mut used = vec![false; nv];
    let mut cut_list: Vec<(u64, [f64; 3])> = vec![];
    for p in &parts {
        for (t, _) in &p.inner {
            for &v in t {
                used[v as usize] = true;
            }
        }
        for (t, _) in &p.surf {
            for &k in t {
                if k & TAG == 0 {
                    used[k as usize] = true;
                }
            }
        }
        cut_list.extend(p.cuts.iter().copied());
    }
    // the same cut point made by two blocks may differ in the last bits (bisected from either end): keep the smallest, so that the choice does
    // not depend on how the lattice was divided into blocks, that is on the number of threads (ADR 11)
    cut_list.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1[0].total_cmp(&b.1[0]).then_with(|| a.1[1].total_cmp(&b.1[1])).then_with(|| a.1[2].total_cmp(&b.1[2]))));
    cut_list.dedup_by_key(|c| c.0);
    let n_cuts = cut_list.len();
    let mut remap = vec![u32::MAX; nv];
    let mut points: Vec<[f64; 3]> = Vec::new();
    // a vertex may lie on a boundary or an interface: a lattice vertex that was moved or whose classes include a zero, and every cut point
    let mut vflag: Vec<bool> = vec![];
    let mut next_id = 0u32;
    for v in 0..nv {
        if used[v] {
            remap[v] = next_id;
            next_id += 1;
            points.push(pos_of(v as u32));
            vflag.push((0..ns).any(|s| cls[s][v] == 0));
        }
    }
    let n_lattice_used = next_id as usize;
    for (_, p) in &cut_list {
        points.push(*p);
        vflag.push(true);
    }
    let vid = |k: u64| -> u32 {
        if k & TAG != 0 {
            n_lattice_used as u32 + cut_list.binary_search_by_key(&k, |c| c.0).expect("cut point present") as u32
        } else {
            remap[k as usize]
        }
    };
    lap("  mesh3d: numbering, cut points");
    crate::progress::check()?;

    // 5. cells (oriented) with their regions
    let mut raw: Vec<([u32; 4], i32)> = Vec::with_capacity(parts.iter().map(|p| p.inner.len() + p.surf.len()).sum());
    // all whole tetrahedra, then all cut ones, each in lattice order: the order must not depend on how the lattice was divided into blocks
    // (that is, on the number of threads), or the numbering of cells would (ADR 11)
    for p in &parts {
        for (t, r) in &p.inner {
            raw.push(([remap[t[0] as usize], remap[t[1] as usize], remap[t[2] as usize], remap[t[3] as usize]], *r));
        }
    }
    for p in &parts {
        for (t, r) in &p.surf {
            raw.push(([vid(t[0]), vid(t[1]), vid(t[2]), vid(t[3])], *r));
        }
    }
    drop(parts);
    let pts = &points;
    let blocks = par_blocks(raw.len(), 20_000, |range| {
        let mut out: Vec<([u32; 4], i32)> = Vec::with_capacity(range.len());
        for c in range {
            let (mut t, mut region) = raw[c];
            let q = t.map(|v| pts[v as usize]);
            if det3(sub3(q[1], q[0]), sub3(q[2], q[0]), sub3(q[3], q[0])) < 0.0 {
                t.swap(1, 2);
            }
            if region < 0 {
                region = resolve_region(dom, &q);
            }
            if region >= 0 {
                out.push((t, region));
            }
        }
        out
    });
    let mut cells: Vec<([u32; 4], i32)> = blocks.into_iter().flatten().collect();
    drop(raw);

    // 5b. weld vertices that coincide to rounding. A surface of the model that lies on a lattice plane (a layer interface, a face of
    // the window of an absorbing layer) warps the lattice vertices onto it and also cuts the edges next to them, so a cut point can be
    // a few ulp away from a warped vertex; the cells between the two are slivers of zero volume. Welding the ends of every edge that is
    // shorter than `eps` removes them without touching a cell of ordinary size (ADR 52).
    let eps2 = (1e-8 * hmin).powi(2);
    let close: Vec<Vec<(u32, u32)>> = par_blocks(cells.len(), 20_000, |r| {
        let mut out = vec![];
        for c in r {
            let t = cells[c].0;
            for i in 0..4 {
                for j in i + 1..4 {
                    let (a, b) = (points[t[i] as usize], points[t[j] as usize]);
                    if (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2) < eps2 {
                        out.push((t[i].min(t[j]), t[i].max(t[j])));
                    }
                }
            }
        }
        out
    });
    if close.iter().any(|b| !b.is_empty()) {
        let mut rep: Vec<u32> = (0..points.len() as u32).collect();
        let find = |rep: &mut Vec<u32>, mut x: u32| {
            while rep[x as usize] != x {
                rep[x as usize] = rep[rep[x as usize] as usize];
                x = rep[x as usize];
            }
            x
        };
        for (a, b) in close.into_iter().flatten() {
            let (ra, rb) = (find(&mut rep, a), find(&mut rep, b));
            if ra != rb {
                rep[ra.max(rb) as usize] = ra.min(rb);
            }
        }
        for v in 0..points.len() as u32 {
            let r = find(&mut rep, v);
            if r != v && vflag[v as usize] {
                vflag[r as usize] = true;
            }
        }
        let mut kept = Vec::with_capacity(cells.len());
        for (t, r) in &cells {
            let t = t.map(|v| find(&mut rep, v));
            if (0..4).all(|i| (i + 1..4).all(|j| t[i] != t[j])) {
                kept.push((t, *r));
            }
        }
        cells = kept;
    }

    // 6. boundary facets: faces of cells with the void on the other side, and faces between cells of different regions. Such a face has
    // all three vertices on a surface.
    let mut keyed: Vec<([u32; 3], u32, u8)> = vec![];
    for (c, (t, _)) in cells.iter().enumerate() {
        if t.iter().filter(|&&v| vflag[v as usize]).count() < 3 {
            continue;
        }
        for skip in 0..4 {
            let f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| t[i]).collect();
            if f.iter().all(|&v| vflag[v as usize]) {
                let mut k = [f[0], f[1], f[2]];
                k.sort_unstable();
                keyed.push((k, c as u32, skip as u8));
            }
        }
    }
    keyed.sort_unstable();
    // (outward triple of the first cell, regions on its two sides)
    let mut facets: Vec<([u32; 3], Option<u32>, Option<u32>)> = vec![];
    let mut i = 0;
    while i < keyed.len() {
        let mut j = i + 1;
        while j < keyed.len() && keyed[j].0 == keyed[i].0 {
            j += 1;
        }
        let outward = |&(_, c, skip): &([u32; 3], u32, u8)| -> [u32; 3] {
            let t = cells[c as usize].0;
            let mut f = [0u32; 3];
            let mut m = 0;
            for (q, &v) in t.iter().enumerate() {
                if q != skip as usize {
                    f[m] = v;
                    m += 1;
                }
            }
            let (a, b, cc, d) = (points[f[0] as usize], points[f[1] as usize], points[f[2] as usize], points[t[skip as usize] as usize]);
            // the normal (b - a) x (c - a) must point away from the opposite vertex
            if det3(sub3(b, a), sub3(cc, a), sub3(d, a)) > 0.0 {
                f.swap(1, 2);
            }
            f
        };
        if j - i == 1 {
            let r = cells[keyed[i].1 as usize].1;
            facets.push((outward(&keyed[i]), Some(r as u32), None));
        } else if j - i == 2 {
            let (r1, r2) = (cells[keyed[i].1 as usize].1, cells[keyed[i + 1].1 as usize].1);
            if r1 != r2 {
                let first = if r1 < r2 { i } else { i + 1 };
                facets.push((outward(&keyed[first]), Some(r1.min(r2) as u32), Some(r1.max(r2) as u32)));
            }
        }
        i = j;
    }
    drop(keyed);
    let flabels: Vec<i32> = par_blocks(facets.len(), 5_000, |r| {
        r.map(|f| {
            let (tri, a, b) = &facets[f];
            let (p0, p1, p2) = (points[tri[0] as usize], points[tri[1] as usize], points[tri[2] as usize]);
            dom.boundary_label([(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0], *a, *b)
        })
        .collect::<Vec<i32>>()
    })
    .into_iter()
    .flatten()
    .collect();
    lap("  mesh3d: cells, labels, facets");
    crate::progress::check()?;

    // 7. the native mesh; vertices that no cell uses (cells dropped for having no region) are left out
    let mut new_id = vec![u32::MAX; points.len()];
    let mut mesh = CellMesh::<3>::new();
    for (t, _) in &cells {
        for &v in t {
            if new_id[v as usize] == u32::MAX {
                new_id[v as usize] = 0;
            }
        }
    }
    for (v, slot) in new_id.iter_mut().enumerate() {
        if *slot == 0 {
            *slot = mesh.points.len() as u32;
            mesh.points.push(points[v]);
        }
    }
    for (t, r) in &cells {
        mesh.push_cell(CellType::Tet, &t.map(|v| new_id[v as usize]), dom.region_label(*r as u32));
    }
    for ((tri, _, _), &l) in facets.iter().zip(&flabels) {
        mesh.push_facet(&tri.map(|v| new_id[v as usize]), l);
    }
    let mut report = Report3d { lattice_vertices: nv, warped: n_warped, cut_points: n_cuts, cells: mesh.n_cells(), refined: 0, curve: None };

    // 8. local refinement inside regions that ask for a finer element size
    let local: Vec<(u32, f64)> = (0..nreg as u32).filter_map(|r| dom.region_size(r).filter(|&hb| hb > 0.0 && hb < 0.999 * spec.h).map(|hb| (r, hb))).collect();
    let faces = dom.asks_size_at();
    if local.is_empty() && !faces {
        return Ok((mesh, report));
    }
    let size = |p: [f64; 3]| -> f64 {
        let s = match dom.region_at(p) {
            Some(r) => local.iter().find(|l| l.0 == r).map_or(spec.h, |l| l.1),
            None => spec.h,
        };
        // the sizes of faces (ADR 83): graded away from the face
        if faces { dom.size_at(p).map_or(s, |f| s.min(f)) } else { s }
    };
    let project = |mid: [f64; 3], pa: [f64; 3], pb: [f64; 3], _kind: super::refine3d::EdgeKind| -> [f64; 3] {
        let len = len3(sub3(pb, pa));
        match dom.project_edge(mid, pa, pb) {
            Some(q) if len3(sub3(q, mid)) <= 0.5 * len => q,
            _ => mid,
        }
    };
    let (mesh, rr) = super::refine3d::refine_to_size(&mesh, &size, &project, 1.3, 60_000_000 / 12)?;
    report.cells = rr.cells_after;
    report.refined = rr.bisections;
    lap("  mesh3d: local refinement");
    crate::progress::check()?;
    Ok((mesh, report))
}

/// The region of a piece none of whose vertices is known to lie in one: the region at the centroid, else at points between the centroid
/// and the corners. `-1` if all of them are in the void (a sliver outside the model).
fn resolve_region(dom: &dyn Domain3, q: &[[f64; 3]; 4]) -> i32 {
    let cen = [(q[0][0] + q[1][0] + q[2][0] + q[3][0]) / 4.0, (q[0][1] + q[1][1] + q[2][1] + q[3][1]) / 4.0, (q[0][2] + q[1][2] + q[2][2] + q[3][2]) / 4.0];
    if let Some(r) = dom.region_at(cen) {
        return r as i32;
    }
    for f in [0.25, 0.5, 0.75] {
        for p in q {
            let x = [cen[0] + f * (p[0] - cen[0]), cen[1] + f * (p[1] - cen[1]), cen[2] + f * (p[2] - cen[2])];
            if let Some(r) = dom.region_at(x) {
                return r as i32;
            }
        }
    }
    -1
}

/// Cuts a tetrahedron (given as four vertices) along surface `s` and keeps the part on `side` (-1: inside the surface, +1: outside);
/// vertices on the surface count for both sides.
fn clip_piece(sf: &Surfaces, s: usize, side: i8, t: &[V; 4], cuts: &mut Vec<(u64, [f64; 3])>, out: &mut Vec<[V; 4]>) {
    let (mut ins, mut outs, mut ons): (Vec<V>, Vec<V>, Vec<V>) = (vec![], vec![], vec![]);
    for k in 0..4 {
        let c = t[k].cl[s];
        if c == side {
            ins.push(t[k]);
        } else if c == -side {
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
        let c = sf.bisect(s, lo.p, hi.p);
        cuts.push((id, c.p));
        V { id, p: c.p, cl: sf.cut_classes(&c), reg: -2 }
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
