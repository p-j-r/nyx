// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Structured curvilinear PML mesh generators for Spherical topologies in 2D Axisymmetric and 3D:
//! - [`build_pml_semicircle_2d`]: extrudes semicircular boundary edges on `R = R_in` (`r >= 0`) radially
//!   into `N_R` concentric rings of [`CellType::Quad`], clamping polar endpoints (`theta = 0, pi`) strictly
//!   to IEEE `r = 0.0`, tagging newly extruded axial segments on `r = 0` with the `r = 0` symmetry boundary
//!   label (never outer PEC), and placing lockstep `CurvedQuad` mid-edge nodes on `R_k` (`k = 0..=N_R`).
//! - [`build_pml_sphere_3d`]: sweeps exterior surface triangles on `R = R_in` radially along unit sphere
//!   normals into `N_R` concentric shells of [`CellType::Prism`] (zero pyramids, zero corner sub-domains),
//!   with lockstep `CurvedPrism` mid-edge projection across shells `k = 0..=N_R`.

use crate::adapt::work::{facet_key, FacetKey, SimplexMesh};
use crate::mesh::cells::{canonical_quad_face, CellMesh, CellType};
use crate::meshing::pml_cyl::{graded_offsets, relax_against_phys_tets};
use std::collections::{BTreeMap, BTreeSet, HashMap};

fn orient_cell<const D: usize>(t: CellType, v: &mut [u32], pts: &[[f64; D]]) {
    let mut probe = CellMesh::<D>::new();
    probe.points = v.iter().map(|&i| pts[i as usize]).collect();
    let local: Vec<u32> = (0..v.len() as u32).collect();
    probe.push_cell(t, &local, 0);
    if crate::mesh::quality::cell_scaled_jacobian_min(&probe, 0) < 0.0 {
        let swaps: &[(usize, usize)] = match t {
            CellType::Tri => &[(1, 2)],
            CellType::Quad => &[(1, 3)],
            CellType::Tet => &[(2, 3)],
            CellType::Prism => &[(1, 2), (4, 5)],
            CellType::Hex => &[(1, 3), (5, 7)],
        };
        for &(i, j) in swaps {
            v.swap(i, j);
        }
    }
}

// ------------------------------------------------------------------ 2D Axisymmetric Semicircle

/// Specification for a 2D axisymmetric spherical semicircle PML (`R_in <= R <= R_in + thickness`, `r >= 0`).
#[derive(Clone, Debug)]
pub struct SpherePmlSpec2d {
    /// Axial center coordinate `z_0` on the symmetry axis `[0.0, z_0]`.
    pub center_z: f64,
    /// Inner spherical radius `R_in > 0`.
    pub r_in: f64,
    /// Radial shell thickness `d > 0`.
    pub thickness: f64,
    /// Number of concentric quad shells `N_R >= 1`.
    pub layers: usize,
    /// Geometric layer width ratio (`> 0`, `1.0` = uniform).
    pub grading: f64,
    /// Place lockstep `CurvedQuad` mid-edge nodes on exact radius `R_k` (`k = 0..=N_R`).
    pub curved: bool,
    /// Encoded domain label is `label_base + material` (`0 <= material < label_stride`).
    pub label_base: i32,
    pub label_stride: i32,
    /// Optional fallback `r = 0` symmetry boundary label near north pole `theta = 0` (`z > z_0`).
    pub axis_label_north: Option<i32>,
    /// Optional fallback `r = 0` symmetry boundary label near south pole `theta = pi` (`z < z_0`).
    pub axis_label_south: Option<i32>,
    /// Optional fallback outer semicircle boundary label at `R = R_in + thickness`.
    pub outer_label: Option<i32>,
    /// Relative tolerance (`tol * r_in`).
    pub tol: f64,
}

impl SpherePmlSpec2d {
    /// Creates a 2D axisymmetric semicircle PML spec centered at `[0.0, center_z]`.
    pub fn new(center_z: f64, r_in: f64, thickness: f64, layers: usize) -> Self {
        Self {
            center_z,
            r_in,
            thickness,
            layers,
            grading: 1.0,
            curved: false,
            label_base: 1000,
            label_stride: 1000,
            axis_label_north: None,
            axis_label_south: None,
            outer_label: None,
            tol: 1e-4,
        }
    }

    pub fn layer_offsets(&self) -> Vec<f64> {
        graded_offsets(self.thickness, self.layers, self.grading)
    }

    pub fn label(&self, material: i32) -> i32 {
        debug_assert!((0..self.label_stride).contains(&material));
        self.label_base + material
    }

    pub fn decode(&self, label: i32) -> Option<i32> {
        let rel = label.checked_sub(self.label_base)?;
        (0..self.label_stride).contains(&rel).then_some(rel)
    }
}

/// Output of [`build_pml_semicircle_2d`].
#[derive(Clone, Debug)]
pub struct SpherePmlMesh2d {
    /// Physical `Tri` cells (`0..n_physical_cells`) followed by `N_R` semicircular rings of `Quad` cells.
    pub mesh: CellMesh<2>,
    pub n_physical_cells: usize,
    pub n_physical_points: usize,
    pub spec: SpherePmlSpec2d,
}

/// Extrudes the outer semicircle `R = R_in` (`r >= 0`) of `phys` radially into `N_R` rings of [`CellType::Quad`],
/// clamping polar vertices strictly to `r = 0.0` and tagging extruded `r = 0` axial segments with the symmetry axis label.
pub fn build_pml_semicircle_2d(phys: &SimplexMesh<2>, spec: &SpherePmlSpec2d) -> Result<SpherePmlMesh2d, String> {
    if !(spec.r_in > 0.0 && spec.r_in.is_finite() && spec.center_z.is_finite()) {
        return Err("SpherePmlSpec2d: r_in must be positive and center_z finite".into());
    }
    if !(spec.thickness > 0.0 && spec.thickness.is_finite()) {
        return Err("SpherePmlSpec2d: thickness must be positive and finite".into());
    }
    if spec.layers == 0 {
        return Err("SpherePmlSpec2d: layers must be at least 1".into());
    }
    if !(spec.grading > 0.0 && spec.grading.is_finite()) {
        return Err("SpherePmlSpec2d: grading must be positive and finite".into());
    }
    if spec.label_stride <= 0 {
        return Err("SpherePmlSpec2d: label_stride must be positive".into());
    }

    let z0 = spec.center_z;
    let tol_abs = (spec.tol * spec.r_in).max(1e-12);
    let rad_of = |p: [f64; 2]| p[0].hypot(p[1] - z0);

    for (vi, p) in phys.points.iter().enumerate() {
        if p[0] < -tol_abs {
            return Err(format!(
                "physical vertex {vi} has negative radial coordinate r = {:.4e}; 2D axisymmetric spherical PML requires r >= 0",
                p[0]
            ));
        }
    }

    let mut phys_facet_label: HashMap<FacetKey, i32> = HashMap::new();
    for f in 0..phys.n_facets() {
        phys_facet_label.insert(facet_key(phys.facet(f)), phys.facet_label[f]);
    }

    let fc = phys.facet_cells();
    let mut ext_keys: Vec<(FacetKey, u32)> = fc
        .iter()
        .filter_map(|(&k, cs)| (cs.len() == 1).then_some((k, cs[0])))
        .collect();
    ext_keys.sort_unstable_by_key(|x| x.0);
    if ext_keys.is_empty() {
        return Err("physical mesh has no exterior boundary edges".into());
    }

    struct ArcEdge2 {
        u: u32,
        v: u32,
        mat: i32,
        bnd_label: Option<i32>,
    }

    let mut arc_edges = Vec::new();
    let mut arc_verts_set: BTreeSet<u32> = BTreeSet::new();
    let mut north_axis_best: Option<(f64, i32)> = None;
    let mut south_axis_best: Option<(f64, i32)> = None;

    for (fk, cell_idx) in ext_keys {
        let cf = phys.cell_facets(cell_idx as usize);
        let mut vs = [fk[0], fk[1]];
        for f in cf.iter().take(3) {
            if facet_key(&f[..2]) == fk {
                vs = [f[0], f[1]];
                break;
            }
        }
        let (pu, pv) = (phys.points[vs[0] as usize], phys.points[vs[1] as usize]);
        let on_axis = pu[0].abs() <= tol_abs && pv[0].abs() <= tol_abs;
        if on_axis {
            if let Some(&lab) = phys_facet_label.get(&fk) {
                let z_mid = 0.5 * (pu[1] + pv[1]);
                if north_axis_best.is_none_or(|(bz, _)| z_mid > bz) {
                    north_axis_best = Some((z_mid, lab));
                }
                if south_axis_best.is_none_or(|(bz, _)| z_mid < bz) {
                    south_axis_best = Some((z_mid, lab));
                }
            }
            continue;
        }
        let (ru, rv) = (rad_of(pu), rad_of(pv));
        if (ru - spec.r_in).abs() > tol_abs || (rv - spec.r_in).abs() > tol_abs {
            return Err(format!(
                "physical boundary edge {vs:?} has radii ({ru:.4e}, {rv:.4e}), expected semicircle R = {:.4e}",
                spec.r_in
            ));
        }
        let mat = phys.cell_label[cell_idx as usize];
        if !(0..spec.label_stride).contains(&mat) {
            return Err(format!("material label {mat} outside 0..label_stride"));
        }
        arc_verts_set.insert(vs[0]);
        arc_verts_set.insert(vs[1]);
        arc_edges.push(ArcEdge2 {
            u: vs[0],
            v: vs[1],
            mat,
            bnd_label: phys_facet_label.get(&fk).copied().or(spec.outer_label),
        });
    }

    if arc_edges.is_empty() {
        return Err("no exterior edges found on semicircle R = r_in".into());
    }

    let north_sym_label = north_axis_best.map(|x| x.1).or(spec.axis_label_north).or(south_axis_best.map(|x| x.1));
    let south_sym_label = south_axis_best.map(|x| x.1).or(spec.axis_label_south).or(north_axis_best.map(|x| x.1));

    let offsets = spec.layer_offsets();
    let mut out = CellMesh::<2>::from_simplex(phys);
    out.facet_offsets = vec![0];
    out.facet_conn.clear();
    out.facet_label.clear();

    // Clamp physical axis vertices strictly to exact IEEE `+0.0`.
    for p in &mut out.points {
        if p[0].abs() <= tol_abs {
            p[0] = 0.0;
        }
    }

    // Sort semicircle vertices from north pole (`theta = 0`) to south pole (`theta = pi`).
    let mut arc_verts: Vec<u32> = arc_verts_set.into_iter().collect();
    arc_verts.sort_by(|&a, &b| {
        let pa = out.points[a as usize];
        let pb = out.points[b as usize];
        let ta = pa[0].atan2(pa[1] - z0);
        let tb = pb[0].atan2(pb[1] - z0);
        ta.total_cmp(&tb).then_with(|| a.cmp(&b))
    });

    let mut shell_vid: HashMap<(u32, usize), u32> = HashMap::with_capacity(arc_verts.len() * (spec.layers + 1));
    let mut pole_verts: Vec<(u32, bool)> = Vec::new(); // (v, is_north)
    for &v in &arc_verts {
        shell_vid.insert((v, 0), v);
        let p0 = out.points[v as usize];
        if p0[0] == 0.0 {
            pole_verts.push((v, p0[1] >= z0));
        }
    }

    for k in 1..=spec.layers {
        let rk = spec.r_in + offsets[k];
        for &v in &arc_verts {
            let p0 = out.points[v as usize];
            let pt = if p0[0] == 0.0 {
                // Exact IEEE `0.0` clamping on `theta = 0` and `theta = pi`!
                let sgn = if p0[1] >= z0 { 1.0 } else { -1.0 };
                [0.0, z0 + sgn * rk]
            } else {
                let dz = p0[1] - z0;
                let r0 = p0[0].hypot(dz).max(1e-300);
                [rk * (p0[0] / r0), z0 + rk * (dz / r0)]
            };
            let id = out.points.len() as u32;
            out.points.push(pt);
            shell_vid.insert((v, k), id);
        }
    }

    let mut boundary_label_map: HashMap<FacetKey, i32> = phys_facet_label;
    for &(vp, is_north) in &pole_verts {
        let sym_opt = if is_north { north_sym_label } else { south_sym_label };
        if let Some(sym_lab) = sym_opt {
            for k in 0..spec.layers {
                let a = shell_vid[&(vp, k)];
                let b = shell_vid[&(vp, k + 1)];
                boundary_label_map.insert(facet_key(&[a, b]), sym_lab);
            }
        }
    }

    let mut edge_nodes: Vec<([u32; 2], [f64; 2])> = Vec::new();
    let mut quad_face_nodes: Vec<([u32; 4], [f64; 2])> = Vec::new();

    for ae in &arc_edges {
        let label = spec.label(ae.mat);
        let pu0 = out.points[ae.u as usize];
        let pv0 = out.points[ae.v as usize];
        let mr = 0.5 * (pu0[0] + pv0[0]);
        let mz = 0.5 * ((pu0[1] - z0) + (pv0[1] - z0));
        let mnrm = mr.hypot(mz);

        if spec.curved && mnrm > 1e-30 {
            let dir_r = (mr / mnrm).max(0.0);
            let dir_z = mz / mnrm;
            for k in 0..=spec.layers {
                let rk = spec.r_in + offsets[k];
                let uk = shell_vid[&(ae.u, k)];
                let vk = shell_vid[&(ae.v, k)];
                edge_nodes.push(([uk.min(vk), uk.max(vk)], [(rk * dir_r).max(0.0), z0 + rk * dir_z]));
            }
        }

        for k in 0..spec.layers {
            let u0 = shell_vid[&(ae.u, k)];
            let v0 = shell_vid[&(ae.v, k)];
            let u1 = shell_vid[&(ae.u, k + 1)];
            let v1 = shell_vid[&(ae.v, k + 1)];
            let mut nodes = [u0, v0, v1, u1];
            orient_cell(CellType::Quad, &mut nodes, &out.points);
            out.push_cell(CellType::Quad, &nodes, label);
            if spec.curved && mnrm > 1e-30 {
                let r_mid = spec.r_in + 0.5 * (offsets[k] + offsets[k + 1]);
                let dir_r = (mr / mnrm).max(0.0);
                let dir_z = mz / mnrm;
                quad_face_nodes.push((
                    canonical_quad_face(nodes),
                    [(r_mid * dir_r).max(0.0), z0 + r_mid * dir_z],
                ));
            }
            if k + 1 == spec.layers {
                if let Some(bl) = ae.bnd_label {
                    boundary_label_map.insert(facet_key(&[u1, v1]), bl);
                }
            }
        }
    }

    for f in out.boundary_faces_ordered() {
        let fk = facet_key(&f);
        if let Some(&bl) = boundary_label_map.get(&fk) {
            out.push_facet(&f, bl);
        }
    }

    if spec.curved {
        out.set_edge_nodes(edge_nodes);
        out.set_quad_face_nodes(quad_face_nodes);
    }

    Ok(SpherePmlMesh2d {
        mesh: out,
        n_physical_cells: phys.n_cells(),
        n_physical_points: phys.n_points(),
        spec: spec.clone(),
    })
}

// ------------------------------------------------------------------ 3D Concentric Spherical Shell

/// Specification for a 3D concentric spherical shell PML (`R_in <= R <= R_in + thickness`).
#[derive(Clone, Debug)]
pub struct SpherePmlSpec3d {
    /// Sphere center `[x_0, y_0, z_0]` in wafer coordinates.
    pub center: [f64; 3],
    /// Inner physical sphere radius `R_in > 0`.
    pub r_in: f64,
    /// Radial shell thickness `d > 0`.
    pub thickness: f64,
    /// Number of concentric `Prism` shells `N_R >= 1`.
    pub layers: usize,
    /// Geometric layer width ratio (`> 0`, `1.0` = uniform).
    pub grading: f64,
    /// Place lockstep `CurvedPrism` mid-edge nodes on exact shells `R_k` (`k = 0..=N_R`).
    pub curved: bool,
    /// Encoded cell label is `label_base + material` (`0 <= material < label_stride`).
    pub label_base: i32,
    pub label_stride: i32,
    /// Fallback outer surface facet label when `phys` has unlabelled exterior faces.
    pub outer_label: Option<i32>,
    /// Relative radial tolerance (`|R - r_in| <= tol * r_in`).
    pub tol: f64,
}

impl SpherePmlSpec3d {
    /// Creates a 3D spherical shell PML spec around `center`.
    pub fn new(center: [f64; 3], r_in: f64, thickness: f64, layers: usize) -> Self {
        Self {
            center,
            r_in,
            thickness,
            layers,
            grading: 1.0,
            curved: false,
            label_base: 1000,
            label_stride: 1000,
            outer_label: None,
            tol: 1e-3,
        }
    }

    pub fn layer_offsets(&self) -> Vec<f64> {
        graded_offsets(self.thickness, self.layers, self.grading)
    }

    pub fn label(&self, material: i32) -> i32 {
        debug_assert!((0..self.label_stride).contains(&material));
        self.label_base + material
    }

    pub fn decode(&self, label: i32) -> Option<i32> {
        let rel = label.checked_sub(self.label_base)?;
        (0..self.label_stride).contains(&rel).then_some(rel)
    }
}

/// Output of [`build_pml_sphere_3d`].
#[derive(Clone, Debug)]
pub struct SpherePmlMesh3d {
    /// Physical `Tet`s (`0..n_physical_cells`) followed by `N_R` concentric shells of `Prism`s (zero pyramids).
    pub mesh: CellMesh<3>,
    pub n_physical_cells: usize,
    pub n_physical_points: usize,
    pub spec: SpherePmlSpec3d,
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Sweeps the spherical boundary `R = R_in` of `phys` radially into `N_R` concentric shells of [`CellType::Prism`].
pub fn build_pml_sphere_3d(phys: &SimplexMesh<3>, spec: &SpherePmlSpec3d) -> Result<SpherePmlMesh3d, String> {
    build_pml_sphere_3d_with_edges(phys, spec, &[])
}

/// Sweeps the spherical boundary `R = R_in` of `phys` radially into `N_R` concentric shells of [`CellType::Prism`],
/// preserving interior `phys_edge_nodes` and placing lockstep `CurvedPrism` `edge_nodes` across `k = 0..=N_R`.
pub fn build_pml_sphere_3d_with_edges(
    phys: &SimplexMesh<3>,
    spec: &SpherePmlSpec3d,
    phys_edge_nodes: &[([u32; 2], [f64; 3])],
) -> Result<SpherePmlMesh3d, String> {
    if !(spec.r_in > 0.0 && spec.r_in.is_finite()) {
        return Err("SpherePmlSpec3d: r_in must be positive and finite".into());
    }
    if !(spec.thickness > 0.0 && spec.thickness.is_finite()) {
        return Err("SpherePmlSpec3d: thickness must be positive and finite".into());
    }
    if spec.layers == 0 {
        return Err("SpherePmlSpec3d: layers must be at least 1".into());
    }
    if !(spec.grading > 0.0 && spec.grading.is_finite()) {
        return Err("SpherePmlSpec3d: grading must be positive and finite".into());
    }
    if spec.label_stride <= 0 {
        return Err("SpherePmlSpec3d: label_stride must be positive".into());
    }

    let c0 = spec.center;
    let tol_abs = (spec.tol * spec.r_in).max(1e-11);

    let mut phys_facet_label: HashMap<FacetKey, i32> = HashMap::new();
    for f in 0..phys.n_facets() {
        phys_facet_label.insert(facet_key(phys.facet(f)), phys.facet_label[f]);
    }

    let fc = phys.facet_cells();
    let mut ext_keys: Vec<(FacetKey, u32)> = fc
        .iter()
        .filter_map(|(&k, cs)| (cs.len() == 1).then_some((k, cs[0])))
        .collect();
    ext_keys.sort_unstable_by_key(|x| x.0);
    if ext_keys.is_empty() {
        return Err("physical mesh has no exterior surface triangles".into());
    }

    struct STri3 {
        vs: [u32; 3],
        mat: i32,
        bnd_label: Option<i32>,
    }

    let mut stris = Vec::with_capacity(ext_keys.len());
    let mut surf_verts_set: BTreeSet<u32> = BTreeSet::new();
    let mut surf_edges_set: BTreeSet<(u32, u32)> = BTreeSet::new();

    for (fk, cell_idx) in ext_keys {
        let cf = phys.cell_facets(cell_idx as usize);
        let mut vs = [fk[0], fk[1], fk[2]];
        for f in cf.iter().take(4) {
            if facet_key(&f[..3]) == fk {
                vs = [f[0], f[1], f[2]];
                break;
            }
        }
        for &v in &vs {
            let rv = norm3(sub3(phys.points[v as usize], c0));
            if (rv - spec.r_in).abs() > tol_abs {
                return Err(format!(
                    "physical surface vertex {v} has radius {rv:.4e}, expected sphere R = {:.4e}",
                    spec.r_in
                ));
            }
            surf_verts_set.insert(v);
        }
        // Orient outward along `+(p0 - c0)`.
        let (p0, p1, p2) = (
            phys.points[vs[0] as usize],
            phys.points[vs[1] as usize],
            phys.points[vs[2] as usize],
        );
        if dot3(cross3(sub3(p1, p0), sub3(p2, p0)), sub3(p0, c0)) < 0.0 {
            vs.swap(1, 2);
        }
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            surf_edges_set.insert((vs[i].min(vs[j]), vs[i].max(vs[j])));
        }
        let mat = phys.cell_label[cell_idx as usize];
        if !(0..spec.label_stride).contains(&mat) {
            return Err(format!("material label {mat} outside 0..label_stride"));
        }
        stris.push(STri3 {
            vs,
            mat,
            bnd_label: phys_facet_label.get(&fk).copied().or(spec.outer_label),
        });
    }

    let offsets = spec.layer_offsets();
    let mut out = CellMesh::<3>::from_simplex(phys);
    out.facet_offsets = vec![0];
    out.facet_conn.clear();
    out.facet_label.clear();

    let surf_verts: Vec<u32> = surf_verts_set.into_iter().collect();
    let mut shell_vid: HashMap<(u32, usize), u32> = HashMap::with_capacity(surf_verts.len() * (spec.layers + 1));
    for &v in &surf_verts {
        shell_vid.insert((v, 0), v);
    }
    for k in 1..=spec.layers {
        let rk = spec.r_in + offsets[k];
        for &v in &surf_verts {
            let d = sub3(phys.points[v as usize], c0);
            let r0 = norm3(d).max(1e-300);
            let s = rk / r0;
            let id = out.points.len() as u32;
            out.points.push([c0[0] + s * d[0], c0[1] + s * d[1], c0[2] + s * d[2]]);
            shell_vid.insert((v, k), id);
        }
    }

    let mut outer_face_label: HashMap<FacetKey, i32> = HashMap::with_capacity(stris.len());
    for st in &stris {
        let label = spec.label(st.mat);
        for k in 0..spec.layers {
            let lo = st.vs.map(|v| shell_vid[&(v, k)]);
            let hi = st.vs.map(|v| shell_vid[&(v, k + 1)]);
            let mut nodes = [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]];
            orient_cell(CellType::Prism, &mut nodes, &out.points);
            out.push_cell(CellType::Prism, &nodes, label);
            if k + 1 == spec.layers {
                if let Some(bl) = st.bnd_label {
                    outer_face_label.insert(facet_key(&hi), bl);
                }
            }
        }
    }

    for f in out.boundary_faces_ordered() {
        let fk = facet_key(&f);
        if let Some(&bl) = outer_face_label.get(&fk).or_else(|| phys_facet_label.get(&fk)) {
            out.push_facet(&f, bl);
        }
    }

    // Lockstep spherical curving across shells `k = 0..=N_R`!
    let mut edge_node_map: BTreeMap<[u32; 2], [f64; 3]> = phys_edge_nodes
        .iter()
        .map(|&([a, b], p)| ([a.min(b), a.max(b)], p))
        .collect();

    if spec.curved {
        let mut surf_edge_rel0: BTreeMap<(u32, u32), [f64; 3]> = BTreeMap::new();
        let mut uncurved_cand: Vec<([u32; 2], [f64; 3])> = Vec::new();

        for &(u, v) in &surf_edges_set {
            let key = [u, v];
            if let Some(&q_exist) = edge_node_map.get(&key) {
                surf_edge_rel0.insert((u, v), sub3(q_exist, c0));
            } else {
                let (pu, pv) = (phys.points[u as usize], phys.points[v as usize]);
                let chord = norm3(sub3(pv, pu));
                let m_rel = [
                    0.5 * ((pu[0] - c0[0]) + (pv[0] - c0[0])),
                    0.5 * ((pu[1] - c0[1]) + (pv[1] - c0[1])),
                    0.5 * ((pu[2] - c0[2]) + (pv[2] - c0[2])),
                ];
                let mnrm = norm3(m_rel);
                if mnrm > 1e-30 && chord > 0.0 {
                    let s = spec.r_in / mnrm;
                    let q_rel = [s * m_rel[0], s * m_rel[1], s * m_rel[2]];
                    if norm3(sub3(q_rel, m_rel)) > 1e-9 * chord {
                        uncurved_cand.push((key, [c0[0] + q_rel[0], c0[1] + q_rel[1], c0[2] + q_rel[2]]));
                    }
                }
            }
        }

        if !uncurved_cand.is_empty() {
            let phys_tets: Vec<[u32; 4]> = (0..phys.n_cells())
                .map(|c| {
                    let cv = phys.cell(c);
                    [cv[0], cv[1], cv[2], cv[3]]
                })
                .collect();
            let _ = relax_against_phys_tets(&phys.points, &phys_tets, &edge_node_map, &mut uncurved_cand);
            for ([u, v], q0) in uncurved_cand {
                edge_node_map.insert([u, v], q0);
                surf_edge_rel0.insert((u, v), sub3(q0, c0));
            }
        }

        let mut quad_face_map: BTreeMap<[u32; 4], [f64; 3]> = BTreeMap::new();
        for (&(u, v), &rel0) in &surf_edge_rel0 {
            for k in 0..=spec.layers {
                let scale = (spec.r_in + offsets[k]) / spec.r_in;
                let uk = shell_vid[&(u, k)];
                let vk = shell_vid[&(v, k)];
                edge_node_map.insert(
                    [uk.min(vk), uk.max(vk)],
                    [c0[0] + scale * rel0[0], c0[1] + scale * rel0[1], c0[2] + scale * rel0[2]],
                );
                if k < spec.layers {
                    let u1 = shell_vid[&(u, k + 1)];
                    let v1 = shell_vid[&(v, k + 1)];
                    let s_mid = (spec.r_in + 0.5 * (offsets[k] + offsets[k + 1])) / spec.r_in;
                    quad_face_map.insert(
                        canonical_quad_face([uk, vk, v1, u1]),
                        [c0[0] + s_mid * rel0[0], c0[1] + s_mid * rel0[1], c0[2] + s_mid * rel0[2]],
                    );
                }
            }
        }
        out.set_quad_face_nodes(quad_face_map);
    }

    if !edge_node_map.is_empty() {
        out.set_edge_nodes(edge_node_map);
    }

    Ok(SpherePmlMesh3d {
        mesh: out,
        n_physical_cells: phys.n_cells(),
        n_physical_points: phys.n_points(),
        spec: spec.clone(),
    })
}
