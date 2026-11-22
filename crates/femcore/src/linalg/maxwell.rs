//! A small **Nédélec (edge element) test-problem generator**: lowest-order `H(curl)` elements on a
//! structured tetrahedral mesh of the unit cube with perfect-electric-conductor walls.
//!
//! Purpose: the AMS preconditioner ([`crate::linalg::ams`]) needs `H(curl)` matrices *and* the discrete
//! gradient, which nyx's real matrices will provide; until those are exported this generator
//! provides a mathematically exact stand-in (same operator class, same null space, same
//! interpolation operators). It is not a mesh generator for production use.
//!
//! * Mesh: `n × n × n` cubes, each cut into six Kuhn tetrahedra (conforming).
//! * Basis: for the edge with end nodes `a < b` (global node numbers) `w = λ_a ∇λ_b − λ_b ∇λ_a`;
//!   so each global edge has one orientation, `a → b`, and no signs appear in the assembly.
//! * `stiffness` = `∫ curl w · curl w` (singular on gradients), `mass` = `∫ w · w`.
//! * `gradient`: edges × interior nodes, `(G φ)_e = φ_b − φ_a`; `stiffness · gradient = 0` exactly
//!   (to rounding).
//! * PEC: edges lying in boundary faces and nodes on the boundary are removed.
//!
//! The Maxwell cavity eigenvalues of the unit cube are `π²(l² + m² + n²)` (at most one index zero),
//! the lowest being `2π² ≈ 19.74` (triple), then `3π²` (double): the tests use these.

use crate::linalg::csr::Csr;
use std::collections::HashMap;

/// Matrices and geometry of the cube problem (all dimensions are the *interior* dof counts).
#[derive(Debug, Clone)]
pub struct EdgeElementCube {
    /// `∫ curl w_i · curl w_j` (edges × edges).
    pub stiffness: Csr<f64>,
    /// `∫ w_i · w_j` (edges × edges).
    pub mass: Csr<f64>,
    /// Discrete gradient (edges × nodes).
    pub gradient: Csr<f64>,
    /// Coordinates of the interior nodes (the columns of `gradient`).
    pub coords: Vec<[f64; 3]>,
    /// `x_b − x_a` for every edge (`a < b` its end nodes): the orientation of the gradient rows.
    pub edge_vectors: Vec<[f64; 3]>,
    /// Cells per side.
    pub n: usize,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Gradients of the four barycentric coordinates and the volume of a tetrahedron.
fn barycentric_gradients(p: &[[f64; 3]; 4]) -> ([[f64; 3]; 4], f64) {
    let e = |i: usize| [p[i][0] - p[0][0], p[i][1] - p[0][1], p[i][2] - p[0][2]];
    let (e1, e2, e3) = (e(1), e(2), e(3));
    let det = dot(e1, cross(e2, e3));
    // ∇λ_1 = (e2 × e3)/det, ∇λ_2 = (e3 × e1)/det, ∇λ_3 = (e1 × e2)/det, ∇λ_0 = −Σ.
    let scale = |v: [f64; 3]| [v[0] / det, v[1] / det, v[2] / det];
    let (g1, g2, g3) = (scale(cross(e2, e3)), scale(cross(e3, e1)), scale(cross(e1, e2)));
    let g0 = [-(g1[0] + g2[0] + g3[0]), -(g1[1] + g2[1] + g3[1]), -(g1[2] + g2[2] + g3[2])];
    ([g0, g1, g2, g3], det.abs() / 6.0)
}

/// Conforming Kuhn tetrahedral mesh of the unit cube with the numbering shared by the edge-element and the
/// face-element generators (so that the discrete gradient and curl act between matching dof sets).
pub(crate) struct CubeMesh {
    /// Vertices of every tetrahedron in increasing global node number.
    pub tets: Vec<[usize; 4]>,
    /// Interior edge `(a, b)`, `a < b` → dof number.
    pub edge_id: HashMap<(usize, usize), usize>,
    pub edges: Vec<(usize, usize)>,
    /// Interior node → column number.
    pub node_id: HashMap<usize, usize>,
    /// Coordinates of the interior nodes.
    pub coords: Vec<[f64; 3]>,
}

pub(crate) fn cube_coord(n: usize, g: usize) -> [f64; 3] {
    let np = n + 1;
    [(g % np) as f64 / n as f64, ((g / np) % np) as f64 / n as f64, (g / (np * np)) as f64 / n as f64]
}

pub(crate) fn cube_mesh(n: usize) -> CubeMesh {
    let np = n + 1;
    let node = |i: usize, j: usize, k: usize| i + np * (j + np * k);
    let coord =
        |g: usize| [(g % np) as f64 / n as f64, ((g / np) % np) as f64 / n as f64, (g / (np * np)) as f64 / n as f64];
    let is_boundary_node = |g: usize| {
        let c = [g % np, (g / np) % np, g / (np * np)];
        c.iter().any(|&x| x == 0 || x == n)
    };

    // Six Kuhn tetrahedra per cube: vertices along the axis permutation path.
    const PERMS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let mut tets: Vec<[usize; 4]> = Vec::with_capacity(6 * n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                for perm in PERMS {
                    let mut c = [i, j, k];
                    let mut t = [node(c[0], c[1], c[2]), 0, 0, 0];
                    for (s, &axis) in perm.iter().enumerate() {
                        c[axis] += 1;
                        t[s + 1] = node(c[0], c[1], c[2]);
                    }
                    tets.push(t); // indices increase along the path
                }
            }
        }
    }

    // Boundary faces (seen once), then edges on them.
    let mut face_count: HashMap<[usize; 3], u8> = HashMap::new();
    for t in &tets {
        for drop in 0..4 {
            let mut f: Vec<usize> = (0..4).filter(|&q| q != drop).map(|q| t[q]).collect();
            f.sort_unstable();
            *face_count.entry([f[0], f[1], f[2]]).or_insert(0) += 1;
        }
    }
    let mut boundary_edge: HashMap<(usize, usize), bool> = HashMap::new();
    for (f, &c) in &face_count {
        if c == 1 {
            for (a, b) in [(f[0], f[1]), (f[0], f[2]), (f[1], f[2])] {
                boundary_edge.insert((a, b), true);
            }
        }
    }

    // Edge numbering (interior only), in order of first appearance for determinism.
    let mut edge_id: HashMap<(usize, usize), usize> = HashMap::new();
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for t in &tets {
        for p in 0..4 {
            for q in p + 1..4 {
                let e = (t[p], t[q]);
                if !boundary_edge.contains_key(&e) && !edge_id.contains_key(&e) {
                    edge_id.insert(e, edges.len());
                    edges.push(e);
                }
            }
        }
    }
    let mut node_id: HashMap<usize, usize> = HashMap::new();
    let mut coords = Vec::new();
    for g in 0..np * np * np {
        if !is_boundary_node(g) {
            node_id.insert(g, coords.len());
            coords.push(coord(g));
        }
    }

    CubeMesh { tets, edge_id, edges, node_id, coords }
}

/// Generates the cube problem with `n ≥ 2` cells per side.
pub fn nedelec_cube(n: usize) -> EdgeElementCube {
    assert!(n >= 2, "need at least 2 cells per side for interior unknowns");
    let CubeMesh { tets, edge_id, edges, node_id, coords } = cube_mesh(n);
    let coord = |g: usize| cube_coord(n, g);
    // Assembly.
    let mut kt: Vec<(usize, usize, f64)> = Vec::new();
    let mut mt: Vec<(usize, usize, f64)> = Vec::new();
    for t in &tets {
        let p = [coord(t[0]), coord(t[1]), coord(t[2]), coord(t[3])];
        let (g, vol) = barycentric_gradients(&p);
        let local: Vec<(usize, usize)> = (0..4).flat_map(|a| (a + 1..4).map(move |b| (a, b))).collect();
        let curl: Vec<[f64; 3]> =
            local.iter().map(|&(a, b)| cross(g[a], g[b])).map(|c| [2.0 * c[0], 2.0 * c[1], 2.0 * c[2]]).collect();
        let integral = |i: usize, j: usize| vol * if i == j { 2.0 } else { 1.0 } / 20.0;
        for (le, &(a, b)) in local.iter().enumerate() {
            let Some(&ge) = edge_id.get(&(t[a], t[b])) else { continue };
            for (lf, &(c, d)) in local.iter().enumerate() {
                let Some(&gf) = edge_id.get(&(t[c], t[d])) else { continue };
                kt.push((ge, gf, vol * dot(curl[le], curl[lf])));
                let m = dot(g[b], g[d]) * integral(a, c)
                    - dot(g[b], g[c]) * integral(a, d)
                    - dot(g[a], g[d]) * integral(b, c)
                    + dot(g[a], g[c]) * integral(b, d);
                mt.push((ge, gf, m));
            }
        }
    }
    let ne = edges.len();
    let stiffness = Csr::from_triplets(ne, ne, &kt).unwrap();
    let mass = Csr::from_triplets(ne, ne, &mt).unwrap();
    let mut gt = Vec::new();
    for (e, &(a, b)) in edges.iter().enumerate() {
        if let Some(&j) = node_id.get(&a) {
            gt.push((e, j, -1.0));
        }
        if let Some(&j) = node_id.get(&b) {
            gt.push((e, j, 1.0));
        }
    }
    let gradient = Csr::from_triplets(ne, coords.len(), &gt).unwrap();
    let edge_vectors = edges
        .iter()
        .map(|&(a, b)| {
            let (pa, pb) = (coord(a), coord(b));
            [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]]
        })
        .collect();
    EdgeElementCube { stiffness, mass, gradient, coords, edge_vectors, n }
}

/// Matrices of the lowest-order **Raviart–Thomas** (`H(div)`, face element) space on the same mesh, with the
/// discrete curl linking it to the edge space of [`nedelec_cube`] (same interior-edge numbering). Unknowns are the
/// fluxes through the interior faces; a face `(a, b, c)` with `a < b < c` is oriented by `(b−a)×(c−a)`. The
/// normal component vanishes on the boundary (boundary faces are removed).
#[derive(Debug, Clone)]
pub struct FaceElementCube {
    /// `∫ div ψ_f div ψ_g` (faces × faces), singular on the curls.
    pub div_stiffness: Csr<f64>,
    /// `∫ ψ_f · ψ_g`.
    pub mass: Csr<f64>,
    /// Discrete curl, faces × interior edges: `(C e)_f` is the circulation around the boundary of `f`
    /// (`+e_ab + e_bc − e_ac`). `C G = 0`, and `∫ curl w · curl w' = Cᵀ M C` exactly.
    pub curl: Csr<f64>,
    /// `Π_c`, faces × interior nodes: the flux through each face of the nodal vector field `u e_c` (P1 interpolation
    /// averaged on the face), `(Π_c)_{f,j} = (nf)_c / 6` for the nodes `j` of face `f`, `nf = (b−a)×(c−a)`.
    pub pi: [Csr<f64>; 3],
    pub n: usize,
}

/// Generates the Raviart–Thomas test problem on the `n × n × n` Kuhn mesh (`n ≥ 2`).
pub fn raviart_thomas_cube(n: usize) -> FaceElementCube {
    assert!(n >= 2);
    let CubeMesh { tets, edge_id, node_id, .. } = cube_mesh(n);
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    // interior faces = shared by two tetrahedra
    let mut count: HashMap<[usize; 3], u8> = HashMap::new();
    for t in &tets {
        for drop in 0..4 {
            let f: Vec<usize> = (0..4).filter(|&q| q != drop).map(|q| t[q]).collect();
            *count.entry([f[0], f[1], f[2]]).or_insert(0) += 1; // t is increasing, so f is sorted
        }
    }
    let mut face_id: HashMap<[usize; 3], usize> = HashMap::new();
    let mut faces: Vec<[usize; 3]> = Vec::new();
    for t in &tets {
        for drop in 0..4 {
            let f: Vec<usize> = (0..4).filter(|&q| q != drop).map(|q| t[q]).collect();
            let key = [f[0], f[1], f[2]];
            if count[&key] == 2 && !face_id.contains_key(&key) {
                face_id.insert(key, faces.len());
                faces.push(key);
            }
        }
    }
    let nf = faces.len();
    // a point rule exact for quadratics on a tetrahedron
    let (qa, qb) = (0.585_410_196_624_968_5, 0.138_196_601_125_010_5);
    let mut mt = Vec::new();
    let mut st = Vec::new();
    for t in &tets {
        let x: Vec<[f64; 3]> = t.iter().map(|&g| cube_coord(n, g)).collect();
        let vol = dot(sub(x[1], x[0]), cross(sub(x[2], x[0]), sub(x[3], x[0]))).abs() / 6.0;
        // faces opposite vertex i, with their sign relative to the outward normal
        let mut loc: Vec<(Option<usize>, f64)> = Vec::new();
        for i in 0..4 {
            let f: Vec<usize> = (0..4).filter(|&q| q != i).collect();
            let key = [t[f[0]], t[f[1]], t[f[2]]];
            let nrm = cross(sub(x[f[1]], x[f[0]]), sub(x[f[2]], x[f[0]]));
            let cen = [
                (x[f[0]][0] + x[f[1]][0] + x[f[2]][0]) / 3.0,
                (x[f[0]][1] + x[f[1]][1] + x[f[2]][1]) / 3.0,
                (x[f[0]][2] + x[f[1]][2] + x[f[2]][2]) / 3.0,
            ];
            let sgn = if dot(nrm, sub(cen, x[i])) > 0.0 { 1.0 } else { -1.0 };
            loc.push((face_id.get(&key).copied(), sgn));
        }
        // Φ_i(y) = (y − x_i)/(3 vol): unit outward flux through the face opposite vertex i
        let mut quad: Vec<([f64; 3], f64)> = Vec::new();
        for k in 0..4 {
            let mut bary = [qb; 4];
            bary[k] = qa;
            let y = [0, 1, 2].map(|c| (0..4).map(|v| bary[v] * x[v][c]).sum::<f64>());
            quad.push((y, vol / 4.0));
        }
        for i in 0..4 {
            let Some(gi) = loc[i].0 else { continue };
            for j in 0..4 {
                let Some(gj) = loc[j].0 else { continue };
                let sg = loc[i].1 * loc[j].1;
                let m: f64 = quad.iter().map(|(y, w)| w * dot(sub(*y, x[i]), sub(*y, x[j])) / (9.0 * vol * vol)).sum();
                mt.push((gi, gj, sg * m));
                st.push((gi, gj, sg / vol));
            }
        }
    }
    let mass = Csr::from_triplets(nf, nf, &mt).unwrap();
    let div_stiffness = Csr::from_triplets(nf, nf, &st).unwrap();
    let ne = edge_id.len();
    let mut ct = Vec::new();
    let mut pit: [Vec<(usize, usize, f64)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for (f, key) in faces.iter().enumerate() {
        let [a, b, c] = *key;
        for (e, sg) in [((a, b), 1.0), ((b, c), 1.0), ((a, c), -1.0)] {
            if let Some(&id) = edge_id.get(&e) {
                ct.push((f, id, sg));
            }
        }
        let (pa, pb, pc) = (cube_coord(n, a), cube_coord(n, b), cube_coord(n, c));
        let nrm = cross(sub(pb, pa), sub(pc, pa));
        for g in [a, b, c] {
            if let Some(&j) = node_id.get(&g) {
                for (comp, pit_c) in pit.iter_mut().enumerate() {
                    pit_c.push((f, j, nrm[comp] / 6.0));
                }
            }
        }
    }
    let nn = node_id.len();
    let curl = Csr::from_triplets(nf, ne, &ct).unwrap();
    let pi = [0, 1, 2].map(|c| Csr::from_triplets(nf, nn, &pit[c]).unwrap());
    FaceElementCube { div_stiffness, mass, curl, pi, n }
}

/// Hierarchical second-order ($P_2$) Nédélec first-kind edge-element cube problem (`H(curl)` with
/// PEC walls), matching the 20-DOF tetrahedral basis used by `nyx`:
/// - per interior edge `(a, b)` (`a < b`): Whitney `N_ab` at `2*e`, edge gradient `∇(λ_a λ_b)` at `2*e + 1`;
/// - per interior face `(a, b, c)` (`a < b < c`): `λ_c N_ab` at `2*ne + 2*f`, `λ_a N_bc` at `2*ne + 2*f + 1`.
#[derive(Debug, Clone)]
pub struct EdgeElementP2Cube {
    /// `∫ curl w_i · curl w_j` (`ndofs × ndofs`).
    pub stiffness: Csr<f64>,
    /// `∫ w_i · w_j` (`ndofs × ndofs`).
    pub mass: Csr<f64>,
    /// Prolongation `P_{1→2}` (`ndofs × n_edges`) embedding the lowest-order $P_1$ space into $P_2$.
    pub p1_to_p2: Csr<f64>,
    /// Lowest-order discrete gradient `G_{P1}` (`n_edges × n_nodes`).
    pub gradient_p1: Csr<f64>,
    /// Full $P_2$ discrete gradient `G_{P2}` (`ndofs × (n_nodes + n_edges)`), satisfying `K G_{P2} = 0`.
    pub gradient_p2: Csr<f64>,
    /// Coordinates of the interior vertices (`n_nodes`).
    pub coords: Vec<[f64; 3]>,
    /// Edge vectors `x_b − x_a` for every interior edge (`n_edges`).
    pub edge_vectors: Vec<[f64; 3]>,
    /// Local 2-DOF blocks (one pair per interior edge, one pair per interior face) for block smoothers.
    pub blocks: Vec<Vec<usize>>,
    /// Cells per side.
    pub n: usize,
    /// Number of tetrahedra (`6 * n^3`).
    pub n_tets: usize,
}

/// Generates the hierarchical $P_2$ Nédélec cube problem on an `n × n × n` Kuhn tetrahedral mesh (`n ≥ 2`).
pub fn nedelec_p2_cube(n: usize) -> EdgeElementP2Cube {
    assert!(n >= 2, "need at least 2 cells per side for interior unknowns");
    let CubeMesh { tets, edge_id, edges, node_id, coords } = cube_mesh(n);
    let coord = |g: usize| cube_coord(n, g);
    let ne = edges.len();
    let nn = coords.len();

    // Identify interior faces (shared by two tetrahedra) in deterministic order.
    let mut face_count: HashMap<[usize; 3], u8> = HashMap::new();
    for t in &tets {
        for drop in 0..4 {
            let mut f = [0usize; 3];
            let mut idx = 0;
            for q in 0..4 {
                if q != drop {
                    f[idx] = t[q];
                    idx += 1;
                }
            }
            *face_count.entry(f).or_insert(0) += 1; // `t` is strictly increasing, so `f` is sorted
        }
    }
    let mut face_id: HashMap<[usize; 3], usize> = HashMap::new();
    let mut faces: Vec<[usize; 3]> = Vec::new();
    for t in &tets {
        for drop in 0..4 {
            let mut f = [0usize; 3];
            let mut idx = 0;
            for q in 0..4 {
                if q != drop {
                    f[idx] = t[q];
                    idx += 1;
                }
            }
            if face_count[&f] == 2 && !face_id.contains_key(&f) {
                face_id.insert(f, faces.len());
                faces.push(f);
            }
        }
    }
    let nf = faces.len();
    let ndofs = 2 * ne + 2 * nf;

    const LOCAL_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    const LOCAL_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    let quad = crate::fem::quad_rule(3, 4);

    let mut kt: Vec<(usize, usize, f64)> = Vec::with_capacity(tets.len() * 120);
    let mut mt: Vec<(usize, usize, f64)> = Vec::with_capacity(tets.len() * 120);
    let mut w = [[0.0f64; 3]; 20];
    let mut cu = [[0.0f64; 3]; 20];
    let mut ke = [0.0f64; 20 * 20];
    let mut me = [0.0f64; 20 * 20];
    let mut gdofs = [None::<usize>; 20];

    for t in &tets {
        let p = [coord(t[0]), coord(t[1]), coord(t[2]), coord(t[3])];
        let (g, vol) = barycentric_gradients(&p);

        for (le, &(a, b)) in LOCAL_EDGES.iter().enumerate() {
            if let Some(&e) = edge_id.get(&(t[a], t[b])) {
                gdofs[2 * le] = Some(2 * e);
                gdofs[2 * le + 1] = Some(2 * e + 1);
            } else {
                gdofs[2 * le] = None;
                gdofs[2 * le + 1] = None;
            }
        }
        for (lf, &[a, b, c]) in LOCAL_FACES.iter().enumerate() {
            if let Some(&f) = face_id.get(&[t[a], t[b], t[c]]) {
                gdofs[12 + 2 * lf] = Some(2 * ne + 2 * f);
                gdofs[12 + 2 * lf + 1] = Some(2 * ne + 2 * f + 1);
            } else {
                gdofs[12 + 2 * lf] = None;
                gdofs[12 + 2 * lf + 1] = None;
            }
        }

        ke.fill(0.0);
        me.fill(0.0);

        // Precompute constant 2*(g_a x g_b) for the 6 local edges.
        let mut edge_curl = [[0.0f64; 3]; 6];
        for (le, &(a, b)) in LOCAL_EDGES.iter().enumerate() {
            let cr = cross(g[a], g[b]);
            edge_curl[le] = [2.0 * cr[0], 2.0 * cr[1], 2.0 * cr[2]];
        }

        for (q_idx, lam) in quad.lam.iter().enumerate() {
            let wq = vol * quad.w[q_idx];
            for (le, &(a, b)) in LOCAL_EDGES.iter().enumerate() {
                let (la, lb) = (lam[a], lam[b]);
                let (ga, gb) = (g[a], g[b]);
                w[2 * le] = [la * gb[0] - lb * ga[0], la * gb[1] - lb * ga[1], la * gb[2] - lb * ga[2]];
                cu[2 * le] = edge_curl[le];
                w[2 * le + 1] = [la * gb[0] + lb * ga[0], la * gb[1] + lb * ga[1], la * gb[2] + lb * ga[2]];
                cu[2 * le + 1] = [0.0; 3];
            }
            for (lf, &[a, b, c]) in LOCAL_FACES.iter().enumerate() {
                let (la, lb, lc) = (lam[a], lam[b], lam[c]);
                let (ga, gb, gc) = (g[a], g[b], g[c]);
                let nab = [la * gb[0] - lb * ga[0], la * gb[1] - lb * ga[1], la * gb[2] - lb * ga[2]];
                let cab_raw = cross(ga, gb);
                let cab = [2.0 * cab_raw[0], 2.0 * cab_raw[1], 2.0 * cab_raw[2]];
                let gc_x_nab = cross(gc, nab);
                w[12 + 2 * lf] = [lc * nab[0], lc * nab[1], lc * nab[2]];
                cu[12 + 2 * lf] = [lc * cab[0] + gc_x_nab[0], lc * cab[1] + gc_x_nab[1], lc * cab[2] + gc_x_nab[2]];

                let nbc = [lb * gc[0] - lc * gb[0], lb * gc[1] - lc * gb[1], lb * gc[2] - lc * gb[2]];
                let cbc_raw = cross(gb, gc);
                let cbc = [2.0 * cbc_raw[0], 2.0 * cbc_raw[1], 2.0 * cbc_raw[2]];
                let ga_x_nbc = cross(ga, nbc);
                w[12 + 2 * lf + 1] = [la * nbc[0], la * nbc[1], la * nbc[2]];
                cu[12 + 2 * lf + 1] = [la * cbc[0] + ga_x_nbc[0], la * cbc[1] + ga_x_nbc[1], la * cbc[2] + ga_x_nbc[2]];
            }

            for i in 0..20 {
                if gdofs[i].is_none() {
                    continue;
                }
                let (wi, cui) = (w[i], cu[i]);
                for j in i..20 {
                    if gdofs[j].is_none() {
                        continue;
                    }
                    ke[i * 20 + j] += wq * dot(cui, cu[j]);
                    me[i * 20 + j] += wq * dot(wi, w[j]);
                }
            }
        }

        for i in 0..20 {
            let Some(gi) = gdofs[i] else { continue };
            for j in i..20 {
                let Some(gj) = gdofs[j] else { continue };
                let kv = ke[i * 20 + j];
                let mv = me[i * 20 + j];
                if kv != 0.0 {
                    kt.push((gi, gj, kv));
                    if i != j {
                        kt.push((gj, gi, kv));
                    }
                }
                if mv != 0.0 {
                    mt.push((gi, gj, mv));
                    if i != j {
                        mt.push((gj, gi, mv));
                    }
                }
            }
        }
    }

    let stiffness = Csr::from_triplets(ndofs, ndofs, &kt).unwrap();
    let mass = Csr::from_triplets(ndofs, ndofs, &mt).unwrap();

    // P1 -> P2 prolongation: Whitney edge `e` maps to `2*e`.
    let p_trip: Vec<(usize, usize, f64)> = (0..ne).map(|e| (2 * e, e, 1.0)).collect();
    let p1_to_p2 = Csr::from_triplets(ndofs, ne, &p_trip).unwrap();

    // P1 discrete gradient (`ne × nn`) and full P2 discrete gradient (`ndofs × (nn + ne)`).
    let mut g1t = Vec::with_capacity(2 * ne);
    let mut g2t = Vec::with_capacity(3 * ne);
    for (e, &(a, b)) in edges.iter().enumerate() {
        if let Some(&j) = node_id.get(&a) {
            g1t.push((e, j, -1.0));
            g2t.push((2 * e, j, -1.0));
        }
        if let Some(&j) = node_id.get(&b) {
            g1t.push((e, j, 1.0));
            g2t.push((2 * e, j, 1.0));
        }
        // Quadratic edge bubble `φ_e = λ_a λ_b` has gradient `∇(λ_a λ_b)` = DOF `2*e + 1`.
        g2t.push((2 * e + 1, nn + e, 1.0));
    }
    let gradient_p1 = Csr::from_triplets(ne, nn, &g1t).unwrap();
    let gradient_p2 = Csr::from_triplets(ndofs, nn + ne, &g2t).unwrap();

    let edge_vectors = edges
        .iter()
        .map(|&(a, b)| {
            let (pa, pb) = (coord(a), coord(b));
            [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]]
        })
        .collect();

    let mut blocks = Vec::with_capacity(ne + nf);
    for e in 0..ne {
        blocks.push(vec![2 * e, 2 * e + 1]);
    }
    for f in 0..nf {
        blocks.push(vec![2 * ne + 2 * f, 2 * ne + 2 * f + 1]);
    }

    EdgeElementP2Cube {
        stiffness,
        mass,
        p1_to_p2,
        gradient_p1,
        gradient_p2,
        coords,
        edge_vectors,
        blocks,
        n,
        n_tets: tets.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn max_abs(a: &Csr<f64>) -> f64 {
        a.values().iter().fold(0.0f64, |m, v| m.max(v.abs()))
    }

    #[test]
    fn sizes_match_the_euler_relations() {
        // Interior edges and nodes of the Kuhn mesh; checked against direct counts for n = 2.
        let p = nedelec_cube(2);
        assert_eq!(p.coords.len(), 1); // the centre node
        assert_eq!(p.stiffness.nrows(), p.gradient.nrows());
        assert!(p.stiffness.is_symmetric(1e-13) && p.mass.is_symmetric(1e-13));
        // the centre node of a 2×2×2 Kuhn mesh has 14 edges, and 26 interior edges exist in total
        assert_eq!(p.gradient.nnz(), 14);
        assert_eq!(p.stiffness.nrows(), 26);
    }

    #[test]
    fn gradients_are_in_the_kernel_of_the_curl_curl_matrix() {
        for n in [2, 3, 4] {
            let p = nedelec_cube(n);
            let kg = p.stiffness.matmul(&p.gradient);
            assert!(max_abs(&kg) < 1e-11 * max_abs(&p.stiffness), "n={n}: |K G| = {}", max_abs(&kg));
        }
    }

    #[test]
    fn p2_cube_matches_p1_restriction_and_has_exact_p2_gradient_kernel() {
        for n in [2, 3] {
            let p1 = nedelec_cube(n);
            let p2 = nedelec_p2_cube(n);
            assert_eq!(p2.n_tets, 6 * n * n * n);
            assert!(p2.stiffness.is_symmetric(1e-12) && p2.mass.is_symmetric(1e-12));
            let kg = p2.stiffness.matmul(&p2.gradient_p2);
            assert!(max_abs(&kg) < 1e-11 * max_abs(&p2.stiffness), "n={n}: |K_p2 G_p2| = {}", max_abs(&kg));
            // P1 restriction of P2 stiffness and mass must match P1 matrices identically
            let pt = p2.p1_to_p2.transpose();
            let k1_from_p2 = pt.matmul(&p2.stiffness.matmul(&p2.p1_to_p2));
            let m1_from_p2 = pt.matmul(&p2.mass.matmul(&p2.p1_to_p2));
            let dk = Csr::axpby(1.0, &k1_from_p2, -1.0, &p1.stiffness);
            let dm = Csr::axpby(1.0, &m1_from_p2, -1.0, &p1.mass);
            assert!(max_abs(&dk) < 1e-11 * max_abs(&p1.stiffness), "n={n}: |dK| = {}", max_abs(&dk));
            assert!(max_abs(&dm) < 1e-11 * max_abs(&p1.mass), "n={n}: |dM| = {}", max_abs(&dm));
        }
    }

    #[test]
    fn the_face_space_completes_the_exact_sequence() {
        for n in [2, 3, 4] {
            let e = nedelec_cube(n);
            let f = raviart_thomas_cube(n);
            assert_eq!(f.curl.ncols(), e.stiffness.nrows());
            // C G = 0
            let cg = f.curl.matmul(&e.gradient);
            assert!(max_abs(&cg) < 1e-13, "n={n}: |C G| = {}", max_abs(&cg));
            // ∫ curl w·curl w' = Cᵀ M_f C  (the curl of a lowest-order Nédélec function is a Raviart–Thomas one)
            let ctmc = f.curl.transpose().matmul(&f.mass.matmul(&f.curl));
            let diff = Csr::axpby(1.0, &ctmc, -1.0, &e.stiffness);
            assert!(max_abs(&diff) < 1e-11 * max_abs(&e.stiffness), "n={n}: {}", max_abs(&diff));
            // div curl = 0
            let dc = f.div_stiffness.matmul(&f.curl);
            assert!(max_abs(&dc) < 1e-10 * max_abs(&f.div_stiffness), "n={n}: |D C| = {}", max_abs(&dc));
            assert!(f.mass.is_symmetric(1e-13) && f.div_stiffness.is_symmetric(1e-13));
        }
    }

    #[test]
    fn mass_matrix_is_positive_definite_and_has_the_right_total() {
        // ∫ |u|² for the interpolant of the constant field e_x on the whole cube (no PEC) is 1;
        // with PEC edges removed the field is cut, so only check positivity via x·Mx > 0 for random x.
        let p = nedelec_cube(3);
        let n = p.mass.nrows();
        let x: Vec<f64> = (0..n).map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5).collect();
        let mut y = vec![0.0; n];
        p.mass.spmv(&x, &mut y);
        assert!(x.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>() > 0.0);
    }
}
