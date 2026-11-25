// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Structured curvilinear PML mesh generators for Cylindrical topologies in 2D and 3D:
//! - [`build_pml_annulus_2d`]: sweeps circular boundary edges on `r = R_in` radially into `N_r` concentric
//!   rings of [`CellType::Quad`], preserving per-facet base materials and placing lockstep `CurvedQuad`
//!   mid-edge nodes on exact circles `r_k` (`k = 0..=N_r`).
//! - [`build_pml_cylinder_3d`]: sweeps barrel triangles on `r = R_in` radially into `N_r` layers of
//!   [`CellType::Prism`] (preserving `z`-stratified wafer materials), bottom/top endcap triangles along
//!   `-+axis` / `+axis` into `N_z` layers of [`CellType::Prism`], circular rim edge chains into `N_r x N_z`
//!   corner annuli of [`CellType::Hex`] (zero pyramids), interior `Op::Sheet` edges touching `r = R_in` into
//!   conductive `Quad` sheet facets, and lockstep cylindrical `edge_nodes` + barrel `quad_face_nodes`.

use crate::adapt::work::{facet_key, FacetKey, SimplexMesh};
use crate::fem::{quad_rule_collapsed, CurvedTet};
use crate::mesh::cells::{canonical_quad_face, CellMesh, CellType};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Geometric layer offsets `0.0 ..= thickness` (`layers + 1` entries) with width growth ratio `grading`.
pub fn graded_offsets(thickness: f64, layers: usize, grading: f64) -> Vec<f64> {
    let mut off = vec![0.0; layers + 1];
    if layers == 0 {
        return off;
    }
    let q = if layers > 1 {
        grading.powf(1.0 / (layers as f64 - 1.0))
    } else {
        1.0
    };
    let total: f64 = (0..layers).map(|k| q.powi(k as i32)).sum();
    let w0 = thickness / total;
    for k in 0..layers {
        off[k + 1] = off[k] + w0 * q.powi(k as i32);
    }
    off[layers] = thickness;
    off
}

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

// ------------------------------------------------------------------ 2D Cylindrical Annulus

/// Specification for a 2D planar concentric cylindrical annulus PML (`R_in <= r <= R_in + thickness`).
#[derive(Clone, Debug)]
pub struct CylPmlSpec2d {
    /// Center `[x_0, z_0]` in the 2D cross-section `(x, z)`.
    pub center: [f64; 2],
    /// Inner physical boundary radius `R_in > 0`.
    pub r_in: f64,
    /// Radial thickness `d > 0` of the PML ring.
    pub thickness: f64,
    /// Number of concentric quad rings `N_r >= 1`.
    pub layers: usize,
    /// Ratio of outermost to innermost radial layer width (`> 0`, `1.0` = uniform).
    pub grading: f64,
    /// Place lockstep `CurvedQuad` mid-edge nodes on exact circle `r_k` for every shell `k = 0..=N_r`.
    pub curved: bool,
    /// Encoded domain label is `label_base + material` (`0 <= material < label_stride`).
    pub label_base: i32,
    pub label_stride: i32,
    /// Fallback label for outer boundary facets when `phys` carries no boundary facet labels.
    pub outer_label: Option<i32>,
    /// Relative radial tolerance (`|r - r_in| <= tol * r_in`).
    pub tol: f64,
}

impl CylPmlSpec2d {
    /// Creates a 2D annulus PML spec with uniform grading (`1.0`) and `tol = 1e-4`.
    pub fn new(center: [f64; 2], r_in: f64, thickness: f64, layers: usize) -> Self {
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
            tol: 1e-4,
        }
    }

    /// Radial layer offsets `0.0 ..= thickness` (`layers + 1` values).
    pub fn layer_offsets(&self) -> Vec<f64> {
        graded_offsets(self.thickness, self.layers, self.grading)
    }

    /// Encodes PML domain label for physical `material`.
    pub fn label(&self, material: i32) -> i32 {
        debug_assert!((0..self.label_stride).contains(&material));
        self.label_base + material
    }

    /// Decodes PML domain label back to physical `material`, or `None` for physical cells.
    pub fn decode(&self, label: i32) -> Option<i32> {
        let rel = label.checked_sub(self.label_base)?;
        (0..self.label_stride).contains(&rel).then_some(rel)
    }
}

/// Output of [`build_pml_annulus_2d`].
#[derive(Clone, Debug)]
pub struct CylPmlMesh2d {
    /// Physical `Tri` cells (`0..n_physical_cells`) followed by `N_r` rings of `Quad` cells.
    pub mesh: CellMesh<2>,
    pub n_physical_cells: usize,
    pub n_physical_points: usize,
    pub spec: CylPmlSpec2d,
}

/// Extrudes the circular boundary `r = R_in` of `phys` radially into `N_r` concentric rings of [`CellType::Quad`].
pub fn build_pml_annulus_2d(phys: &SimplexMesh<2>, spec: &CylPmlSpec2d) -> Result<CylPmlMesh2d, String> {
    if !(spec.r_in > 0.0 && spec.r_in.is_finite()) {
        return Err("CylPmlSpec2d: r_in must be positive and finite".into());
    }
    if !(spec.thickness > 0.0 && spec.thickness.is_finite()) {
        return Err("CylPmlSpec2d: thickness must be positive and finite".into());
    }
    if spec.layers == 0 {
        return Err("CylPmlSpec2d: layers must be at least 1".into());
    }
    if !(spec.grading > 0.0 && spec.grading.is_finite()) {
        return Err("CylPmlSpec2d: grading must be positive and finite".into());
    }
    if spec.label_stride <= 0 {
        return Err("CylPmlSpec2d: label_stride must be positive".into());
    }

    let c0 = spec.center;
    let r_tol = (spec.tol * spec.r_in).max(1e-12);
    let rad_of = |p: [f64; 2]| (p[0] - c0[0]).hypot(p[1] - c0[1]);

    // Map physical boundary facet keys to their physical boundary labels.
    let mut phys_facet_label: HashMap<FacetKey, i32> = HashMap::new();
    for f in 0..phys.n_facets() {
        phys_facet_label.insert(facet_key(phys.facet(f)), phys.facet_label[f]);
    }

    // Collect exterior boundary edges of `phys` in deterministic order.
    let fc = phys.facet_cells();
    let mut ext_keys: Vec<(FacetKey, u32)> = fc
        .iter()
        .filter_map(|(&k, cs)| (cs.len() == 1).then_some((k, cs[0])))
        .collect();
    ext_keys.sort_unstable_by_key(|x| x.0);
    if ext_keys.is_empty() {
        return Err("physical mesh has no exterior boundary edges".into());
    }

    struct BEdge2 {
        u: u32,
        v: u32,
        mat: i32,
        bnd_label: Option<i32>,
    }

    let mut bedges = Vec::with_capacity(ext_keys.len());
    let mut rim_verts_set: BTreeSet<u32> = BTreeSet::new();
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
        let (ru, rv) = (rad_of(pu), rad_of(pv));
        if (ru - spec.r_in).abs() > r_tol || (rv - spec.r_in).abs() > r_tol {
            return Err(format!(
                "physical boundary edge {:?} has radii ({ru:.4e}, {rv:.4e}), expected circle r = {:.4e}; cylindrical annulus PML requires a circular physical boundary",
                vs, spec.r_in
            ));
        }
        let mat = phys.cell_label[cell_idx as usize];
        if !(0..spec.label_stride).contains(&mat) {
            return Err(format!("material label {mat} outside 0..label_stride"));
        }
        rim_verts_set.insert(vs[0]);
        rim_verts_set.insert(vs[1]);
        bedges.push(BEdge2 {
            u: vs[0],
            v: vs[1],
            mat,
            bnd_label: phys_facet_label.get(&fk).copied().or(spec.outer_label),
        });
    }

    let offsets = spec.layer_offsets();
    let mut out = CellMesh::<2>::from_simplex(phys);
    out.facet_offsets = vec![0];
    out.facet_conn.clear();
    out.facet_label.clear();

    // Create shell vertices deterministically in `(layer_k, angle_order)` order.
    let mut rim_verts: Vec<u32> = rim_verts_set.into_iter().collect();
    rim_verts.sort_by(|&a, &b| {
        let pa = phys.points[a as usize];
        let pb = phys.points[b as usize];
        let ta = (pa[1] - c0[1]).atan2(pa[0] - c0[0]);
        let tb = (pb[1] - c0[1]).atan2(pb[0] - c0[0]);
        ta.total_cmp(&tb).then_with(|| a.cmp(&b))
    });

    let mut shell_vid: HashMap<(u32, usize), u32> = HashMap::with_capacity(rim_verts.len() * (spec.layers + 1));
    for &v in &rim_verts {
        shell_vid.insert((v, 0), v);
    }
    for k in 1..=spec.layers {
        let rk = spec.r_in + offsets[k];
        for &v in &rim_verts {
            let p0 = phys.points[v as usize];
            let dx = p0[0] - c0[0];
            let dz = p0[1] - c0[1];
            let r0 = dx.hypot(dz).max(1e-300);
            let id = out.points.len() as u32;
            out.points.push([c0[0] + rk * (dx / r0), c0[1] + rk * (dz / r0)]);
            shell_vid.insert((v, k), id);
        }
    }

    let mut outer_edge_label: HashMap<FacetKey, i32> = HashMap::with_capacity(bedges.len());
    let mut edge_nodes: Vec<([u32; 2], [f64; 2])> = Vec::new();
    let mut quad_face_nodes: Vec<([u32; 4], [f64; 2])> = Vec::new();

    for be in &bedges {
        let label = spec.label(be.mat);
        let pu0 = phys.points[be.u as usize];
        let pv0 = phys.points[be.v as usize];
        let mx = 0.5 * ((pu0[0] - c0[0]) + (pv0[0] - c0[0]));
        let mz = 0.5 * ((pu0[1] - c0[1]) + (pv0[1] - c0[1]));
        let mnrm = mx.hypot(mz);

        if spec.curved && mnrm > 1e-30 {
            let dir = [mx / mnrm, mz / mnrm];
            for k in 0..=spec.layers {
                let rk = spec.r_in + offsets[k];
                let uk = shell_vid[&(be.u, k)];
                let vk = shell_vid[&(be.v, k)];
                edge_nodes.push(([uk.min(vk), uk.max(vk)], [c0[0] + rk * dir[0], c0[1] + rk * dir[1]]));
            }
        }

        for k in 0..spec.layers {
            let u0 = shell_vid[&(be.u, k)];
            let v0 = shell_vid[&(be.v, k)];
            let u1 = shell_vid[&(be.u, k + 1)];
            let v1 = shell_vid[&(be.v, k + 1)];
            let mut nodes = [u0, v0, v1, u1];
            orient_cell(CellType::Quad, &mut nodes, &out.points);
            out.push_cell(CellType::Quad, &nodes, label);
            if spec.curved && mnrm > 1e-30 {
                let r_mid = spec.r_in + 0.5 * (offsets[k] + offsets[k + 1]);
                let dir = [mx / mnrm, mz / mnrm];
                quad_face_nodes.push((
                    canonical_quad_face(nodes),
                    [c0[0] + r_mid * dir[0], c0[1] + r_mid * dir[1]],
                ));
            }
            if k + 1 == spec.layers {
                if let Some(bl) = be.bnd_label {
                    outer_edge_label.insert(facet_key(&[u1, v1]), bl);
                }
            }
        }
    }

    for f in out.boundary_faces_ordered() {
        let fk = facet_key(&f);
        if let Some(&bl) = outer_edge_label.get(&fk).or_else(|| phys_facet_label.get(&fk)) {
            out.push_facet(&f, bl);
        }
    }

    if spec.curved {
        out.set_edge_nodes(edge_nodes);
        out.set_quad_face_nodes(quad_face_nodes);
    }

    Ok(CylPmlMesh2d {
        mesh: out,
        n_physical_cells: phys.n_cells(),
        n_physical_points: phys.n_points(),
        spec: spec.clone(),
    })
}

// ------------------------------------------------------------------ 3D Cylindrical Shell + Endcaps + Rim Annuli

/// Canonical sub-domain roles of a 3D Cylindrical PML.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CylRegionRole {
    /// Radial barrel shell (`r_in <= r <= r_in + r_thickness`, `z_min <= zeta <= z_max`).
    RadialShell = 0,
    /// Lower axial endcap (`r <= r_in`, `z_min - z_min_thickness <= zeta <= z_min`).
    EndCapMin = 1,
    /// Upper axial endcap (`r <= r_in`, `z_max <= zeta <= z_max + z_max_thickness`).
    EndCapMax = 2,
    /// Lower corner rim annulus (`r_in <= r <= r_in + r_thickness`, `z_min - z_min_thickness <= zeta <= z_min`).
    CornerAnnulusMin = 3,
    /// Upper corner rim annulus (`r_in <= r <= r_in + r_thickness`, `z_max <= zeta <= z_max + z_max_thickness`).
    CornerAnnulusMax = 4,
}

impl CylRegionRole {
    pub fn from_index(idx: i32) -> Option<Self> {
        match idx {
            0 => Some(CylRegionRole::RadialShell),
            1 => Some(CylRegionRole::EndCapMin),
            2 => Some(CylRegionRole::EndCapMax),
            3 => Some(CylRegionRole::CornerAnnulusMin),
            4 => Some(CylRegionRole::CornerAnnulusMax),
            _ => None,
        }
    }
}

/// Specification for a 3D cylindrical PML (barrel `Prism`s + endcap `Prism`s + corner rim `Hex`es).
#[derive(Clone, Debug)]
pub struct CylPmlSpec3d {
    /// Point on the cylinder axis in wafer coordinates `(x, y, z)`.
    pub center: [f64; 3],
    /// Unit cylinder axis vector `a_hat` in wafer coordinates `(x, y, z)` (typically `[0, 0, 1]`).
    pub axis: [f64; 3],
    /// Barrel inner radius `R_in > 0`.
    pub r_in: f64,
    /// Physical axial coordinate interval `[z_min, z_max]` (`x . a_hat` along `axis`).
    pub z_range: [f64; 2],
    /// Radial PML thickness (`0.0` if radial barrel PML is inactive).
    pub r_thickness: f64,
    /// Number of radial PML layers `N_r` (`0` if radial barrel PML is inactive).
    pub r_layers: usize,
    /// Bottom endcap PML thickness (`0.0` if `-axis` endcap is inactive).
    pub z_min_thickness: f64,
    /// Number of bottom endcap PML layers `N_{z-}` (`0` if inactive).
    pub z_min_layers: usize,
    /// Top endcap PML thickness (`0.0` if `+axis` endcap is inactive).
    pub z_max_thickness: f64,
    /// Number of top endcap PML layers `N_{z+}` (`0` if inactive).
    pub z_max_layers: usize,
    /// Geometric layer growth ratio (`> 0`, `1.0` = uniform).
    pub grading: f64,
    /// Place lockstep cylindrical `edge_nodes` and barrel `quad_face_nodes`.
    pub curved: bool,
    /// Encoded cell label is `label_base + (role as i32) * label_stride + material`.
    pub label_base: i32,
    pub label_stride: i32,
    /// Fallback outer boundary facet label when `phys` has unlabelled exterior faces.
    pub outer_label: Option<i32>,
    /// Relative tolerance (`tol * max(r_in, z_max - z_min)`).
    pub tol: f64,
}

impl CylPmlSpec3d {
    /// Creates a full 3D cylindrical PML with identical `thickness` and `layers` on the barrel and both endcaps.
    pub fn uniform(center: [f64; 3], axis: [f64; 3], r_in: f64, z_range: [f64; 2], thickness: f64, layers: usize) -> Self {
        Self {
            center,
            axis,
            r_in,
            z_range,
            r_thickness: thickness,
            r_layers: layers,
            z_min_thickness: thickness,
            z_min_layers: layers,
            z_max_thickness: thickness,
            z_max_layers: layers,
            grading: 1.0,
            curved: false,
            label_base: 1000,
            label_stride: 1000,
            outer_label: None,
            tol: 1e-4,
        }
    }

    #[inline]
    pub fn radial_active(&self) -> bool {
        self.r_layers > 0 && self.r_thickness > 0.0
    }

    #[inline]
    pub fn z_min_active(&self) -> bool {
        self.z_min_layers > 0 && self.z_min_thickness > 0.0
    }

    #[inline]
    pub fn z_max_active(&self) -> bool {
        self.z_max_layers > 0 && self.z_max_thickness > 0.0
    }

    pub fn r_offsets(&self) -> Vec<f64> {
        graded_offsets(self.r_thickness, self.r_layers, self.grading)
    }

    pub fn z_min_offsets(&self) -> Vec<f64> {
        graded_offsets(self.z_min_thickness, self.z_min_layers, self.grading)
    }

    pub fn z_max_offsets(&self) -> Vec<f64> {
        graded_offsets(self.z_max_thickness, self.z_max_layers, self.grading)
    }

    pub fn label(&self, role: CylRegionRole, material: i32) -> i32 {
        debug_assert!((0..self.label_stride).contains(&material));
        self.label_base + (role as i32) * self.label_stride + material
    }

    pub fn decode(&self, label: i32) -> Option<(CylRegionRole, i32)> {
        let rel = label.checked_sub(self.label_base)?;
        if rel < 0 {
            return None;
        }
        let role = CylRegionRole::from_index(rel / self.label_stride)?;
        let mat = rel % self.label_stride;
        Some((role, mat))
    }
}

/// Output of [`build_pml_cylinder_3d`].
#[derive(Clone, Debug)]
pub struct CylPmlMesh3d {
    /// Physical `Tet`s (`0..n_physical_cells`) followed by barrel/endcap `Prism`s and corner rim `Hex`es.
    pub mesh: CellMesh<3>,
    pub n_physical_cells: usize,
    pub n_physical_points: usize,
    /// Interior conductive `Quad` sheet facets (`(quad_verts, sheet_label)`) extruded radially from `r = R_in`.
    pub extruded_sheets: Vec<([u32; 4], i32)>,
    pub spec: CylPmlSpec3d,
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
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

/// Builds a 3D structured cylindrical PML around `phys` (no pre-existing curved edges or interior sheets).
pub fn build_pml_cylinder_3d(phys: &SimplexMesh<3>, spec: &CylPmlSpec3d) -> Result<CylPmlMesh3d, String> {
    build_pml_cylinder_3d_with_edges_and_sheets(phys, spec, &[], &[])
}

/// Builds a 3D structured cylindrical PML around `phys`, preserving `phys_edge_nodes` on interior `Tet`s,
/// extruding interior `Op::Sheet` edges touching `r = R_in` into conductive `Quad` facets, and placing
/// lockstep cylindrical `edge_nodes` + curved barrel `quad_face_nodes`.
pub fn build_pml_cylinder_3d_with_edges_and_sheets(
    phys: &SimplexMesh<3>,
    spec: &CylPmlSpec3d,
    phys_edge_nodes: &[([u32; 2], [f64; 3])],
    interior_sheets: &[([u32; 3], i32)],
) -> Result<CylPmlMesh3d, String> {
    if !(spec.r_in > 0.0 && spec.r_in.is_finite()) {
        return Err("CylPmlSpec3d: r_in must be positive and finite".into());
    }
    if !(spec.grading > 0.0 && spec.grading.is_finite()) {
        return Err("CylPmlSpec3d: grading must be positive and finite".into());
    }
    if spec.label_stride <= 0 {
        return Err("CylPmlSpec3d: label_stride must be positive".into());
    }
    let ax_len = norm3(spec.axis);
    if !(ax_len > 1e-12 && ax_len.is_finite()) {
        return Err("CylPmlSpec3d: axis must be non-zero".into());
    }
    let a_hat = [spec.axis[0] / ax_len, spec.axis[1] / ax_len, spec.axis[2] / ax_len];
    let (z_min, z_max) = (spec.z_range[0], spec.z_range[1]);
    if !(z_max > z_min && z_min.is_finite() && z_max.is_finite()) {
        return Err("CylPmlSpec3d: z_range must have z_max > z_min".into());
    }

    let c_dot_a = dot3(spec.center, a_hat);
    // Axis point at world axial coordinate `zeta = x . a_hat`:
    let axis_pt = |zeta: f64| -> [f64; 3] {
        let s = zeta - c_dot_a;
        [
            spec.center[0] + s * a_hat[0],
            spec.center[1] + s * a_hat[1],
            spec.center[2] + s * a_hat[2],
        ]
    };
    let decompose = |p: [f64; 3]| -> (f64, [f64; 3], f64) {
        let zeta = dot3(p, a_hat);
        let c_ax = axis_pt(zeta);
        let r_vec = sub3(p, c_ax);
        let r = norm3(r_vec);
        (zeta, r_vec, r)
    };

    let char_len = spec.r_in.max(z_max - z_min);
    let tol_abs = (spec.tol * char_len).max(1e-11);

    let mut phys_facet_label: HashMap<FacetKey, i32> = HashMap::new();
    for f in 0..phys.n_facets() {
        phys_facet_label.insert(facet_key(phys.facet(f)), phys.facet_label[f]);
    }

    // Classify all exterior triangles of `phys`.
    let fc = phys.facet_cells();
    let mut ext_keys: Vec<(FacetKey, u32)> = fc
        .iter()
        .filter_map(|(&k, cs)| (cs.len() == 1).then_some((k, cs[0])))
        .collect();
    ext_keys.sort_unstable_by_key(|x| x.0);
    if ext_keys.is_empty() {
        return Err("physical mesh has no exterior boundary facets".into());
    }

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum PatchKind {
        Barrel,
        CapMin,
        CapMax,
        InactiveOther,
    }

    struct BTri3 {
        vs: [u32; 3],
        kind: PatchKind,
        mat: i32,
        bnd_label: Option<i32>,
    }

    let mut btris = Vec::with_capacity(ext_keys.len());
    for (fk, cell_idx) in ext_keys {
        let cf = phys.cell_facets(cell_idx as usize);
        let mut vs = [fk[0], fk[1], fk[2]];
        for f in cf.iter().take(4) {
            if facet_key(&f[..3]) == fk {
                vs = [f[0], f[1], f[2]];
                break;
            }
        }
        // Ensure outward orientation from owning tet `cell_idx`.
        let cell_vs = phys.cell(cell_idx as usize);
        let opp = cell_vs.iter().copied().find(|v| !vs.contains(v)).unwrap_or(vs[0]);
        let (p0, p1, p2, po) = (
            phys.points[vs[0] as usize],
            phys.points[vs[1] as usize],
            phys.points[vs[2] as usize],
            phys.points[opp as usize],
        );
        let n_tri = cross3(sub3(p1, p0), sub3(p2, p0));
        if dot3(n_tri, sub3(p0, po)) < 0.0 {
            vs.swap(1, 2);
        }
        let n_out = cross3(sub3(phys.points[vs[1] as usize], p0), sub3(phys.points[vs[2] as usize], p0));
        let n_len = norm3(n_out).max(1e-300);
        let n_unit = [n_out[0] / n_len, n_out[1] / n_len, n_out[2] / n_len];
        let dec = vs.map(|v| decompose(phys.points[v as usize]));

        let on_zmin = dec.iter().all(|d| (d.0 - z_min).abs() <= tol_abs) && dot3(n_unit, a_hat) < -0.5;
        let on_zmax = dec.iter().all(|d| (d.0 - z_max).abs() <= tol_abs) && dot3(n_unit, a_hat) > 0.5;
        let on_barrel = dec.iter().all(|d| {
            (d.2 - spec.r_in).abs() <= tol_abs && d.0 >= z_min - tol_abs && d.0 <= z_max + tol_abs
        }) && dot3(n_unit, a_hat).abs() < 0.5;

        let kind = if on_zmin {
            PatchKind::CapMin
        } else if on_zmax {
            PatchKind::CapMax
        } else if on_barrel {
            PatchKind::Barrel
        } else if spec.radial_active() && spec.z_min_active() && spec.z_max_active() {
            return Err(format!(
                "physical boundary triangle {vs:?} does not lie on cylinder r = {:.4e} or endcaps [{z_min:.4e}, {z_max:.4e}]",
                spec.r_in
            ));
        } else {
            PatchKind::InactiveOther
        };

        let mat = phys.cell_label[cell_idx as usize];
        if !(0..spec.label_stride).contains(&mat) {
            return Err(format!("material label {mat} outside 0..label_stride"));
        }
        btris.push(BTri3 {
            vs,
            kind,
            mat,
            bnd_label: phys_facet_label.get(&fk).copied().or(spec.outer_label),
        });
    }

    let r_off = spec.r_offsets();
    let z_min_off = spec.z_min_offsets();
    let z_max_off = spec.z_max_offsets();

    let mut out = CellMesh::<3>::from_simplex(phys);
    out.facet_offsets = vec![0];
    out.facet_conn.clear();
    out.facet_label.clear();

    // Shifted vertex cache keyed by `(base_v, kr, kz_signed)` (`kz_signed < 0` for CapMin, `> 0` for CapMax).
    let mut shifted_cache: HashMap<(u32, usize, isize), u32> = HashMap::new();
    let shifted_pos = |base_v: u32, kr: usize, kz: isize| -> [f64; 3] {
        let p0 = phys.points[base_v as usize];
        if kr == 0 && kz == 0 {
            return p0;
        }
        let (zeta0, r_vec0, r0) = decompose(p0);
        let zeta = if kz == 0 {
            zeta0
        } else if kz < 0 {
            z_min - z_min_off[(-kz) as usize]
        } else {
            z_max + z_max_off[kz as usize]
        };
        let r_vec = if kr == 0 {
            r_vec0
        } else {
            let scale = (spec.r_in + r_off[kr]) / r0.max(1e-300);
            [r_vec0[0] * scale, r_vec0[1] * scale, r_vec0[2] * scale]
        };
        let c_ax = axis_pt(zeta);
        [c_ax[0] + r_vec[0], c_ax[1] + r_vec[1], c_ax[2] + r_vec[2]]
    };

    // Pre-create vertices deterministically in `(layer_depth, base_v)` order.
    let mut barrel_verts: BTreeSet<u32> = BTreeSet::new();
    let mut cap_min_verts: BTreeSet<u32> = BTreeSet::new();
    let mut cap_max_verts: BTreeSet<u32> = BTreeSet::new();
    for bt in &btris {
        match bt.kind {
            PatchKind::Barrel => barrel_verts.extend(bt.vs),
            PatchKind::CapMin => cap_min_verts.extend(bt.vs),
            PatchKind::CapMax => cap_max_verts.extend(bt.vs),
            PatchKind::InactiveOther => {}
        }
    }
    let rim_min_verts: Vec<u32> = barrel_verts.intersection(&cap_min_verts).copied().collect();
    let rim_max_verts: Vec<u32> = barrel_verts.intersection(&cap_max_verts).copied().collect();

    let mut ensure_vertex = |out: &mut CellMesh<3>, base_v: u32, kr: usize, kz: isize| -> u32 {
        if kr == 0 && kz == 0 {
            shifted_cache.insert((base_v, 0, 0), base_v);
            return base_v;
        }
        if let Some(&id) = shifted_cache.get(&(base_v, kr, kz)) {
            return id;
        }
        let id = out.points.len() as u32;
        out.points.push(shifted_pos(base_v, kr, kz));
        shifted_cache.insert((base_v, kr, kz), id);
        id
    };

    if spec.radial_active() {
        for kr in 1..=spec.r_layers {
            for &v in &barrel_verts {
                ensure_vertex(&mut out, v, kr, 0);
            }
        }
    }
    if spec.z_min_active() {
        for j in 1..=spec.z_min_layers {
            let kz = -(j as isize);
            for &v in &cap_min_verts {
                ensure_vertex(&mut out, v, 0, kz);
            }
        }
    }
    if spec.z_max_active() {
        for j in 1..=spec.z_max_layers {
            let kz = j as isize;
            for &v in &cap_max_verts {
                ensure_vertex(&mut out, v, 0, kz);
            }
        }
    }
    if spec.radial_active() && spec.z_min_active() {
        for j in 1..=spec.z_min_layers {
            let kz = -(j as isize);
            for kr in 1..=spec.r_layers {
                for &v in &rim_min_verts {
                    ensure_vertex(&mut out, v, kr, kz);
                }
            }
        }
    }
    if spec.radial_active() && spec.z_max_active() {
        for j in 1..=spec.z_max_layers {
            let kz = j as isize;
            for kr in 1..=spec.r_layers {
                for &v in &rim_max_verts {
                    ensure_vertex(&mut out, v, kr, kz);
                }
            }
        }
    }

    let mut outer_face_label: HashMap<FacetKey, i32> = HashMap::new();

    // 1. Barrel triangles -> `N_r` layers of `Prism` (preserving `z`-stratified materials).
    // Also record barrel edges `(u.min(v), u.max(v)) -> (u, v, mat, bnd_label)` to drive rim `Hex` rings!
    type BarrelEdgeItem = (u32, u32, i32, Option<i32>);
    type RimEdgeItem = (u32, u32, i32, Option<i32>, Option<i32>);
    let mut barrel_edge_info: BTreeMap<(u32, u32), Vec<BarrelEdgeItem>> = BTreeMap::new();
    for bt in &btris {
        if bt.kind != PatchKind::Barrel {
            continue;
        }
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            let (u, v) = (bt.vs[i], bt.vs[j]);
            barrel_edge_info
                .entry((u.min(v), u.max(v)))
                .or_default()
                .push((u, v, bt.mat, bt.bnd_label));
        }
        if !spec.radial_active() {
            continue;
        }
        let label = spec.label(CylRegionRole::RadialShell, bt.mat);
        for kr in 0..spec.r_layers {
            let lo = bt.vs.map(|v| ensure_vertex(&mut out, v, kr, 0));
            let hi = bt.vs.map(|v| ensure_vertex(&mut out, v, kr + 1, 0));
            let mut nodes = [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]];
            orient_cell(CellType::Prism, &mut nodes, &out.points);
            out.push_cell(CellType::Prism, &nodes, label);
            if kr + 1 == spec.r_layers {
                if let Some(bl) = bt.bnd_label {
                    outer_face_label.insert(facet_key(&hi), bl);
                }
            }
            // If an endcap is INACTIVE while radial IS active, the exposed annulus quads on `z_min` / `z_max` are exterior!
        }
    }

    // 2. Bottom & Top EndCap triangles -> `N_z` layers of `Prism`.
    // Also collect cap rim edges touching `Barrel` to build CornerAnnulus `Hex`es!
    let mut cap_min_edge_label: HashMap<(u32, u32), Option<i32>> = HashMap::new();
    let mut cap_max_edge_label: HashMap<(u32, u32), Option<i32>> = HashMap::new();
    for bt in &btris {
        match bt.kind {
            PatchKind::CapMin => {
                for (i, j) in [(0, 1), (1, 2), (2, 0)] {
                    let (u, v) = (bt.vs[i], bt.vs[j]);
                    cap_min_edge_label.insert((u.min(v), u.max(v)), bt.bnd_label);
                }
                if spec.z_min_active() {
                    let label = spec.label(CylRegionRole::EndCapMin, bt.mat);
                    for j in 0..spec.z_min_layers {
                        let kz0 = -(j as isize);
                        let kz1 = -((j + 1) as isize);
                        let lo = bt.vs.map(|v| ensure_vertex(&mut out, v, 0, kz0));
                        let hi = bt.vs.map(|v| ensure_vertex(&mut out, v, 0, kz1));
                        let mut nodes = [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]];
                        orient_cell(CellType::Prism, &mut nodes, &out.points);
                        out.push_cell(CellType::Prism, &nodes, label);
                        if j + 1 == spec.z_min_layers {
                            if let Some(bl) = bt.bnd_label {
                                outer_face_label.insert(facet_key(&hi), bl);
                            }
                        }
                    }
                }
            }
            PatchKind::CapMax => {
                for (i, j) in [(0, 1), (1, 2), (2, 0)] {
                    let (u, v) = (bt.vs[i], bt.vs[j]);
                    cap_max_edge_label.insert((u.min(v), u.max(v)), bt.bnd_label);
                }
                if spec.z_max_active() {
                    let label = spec.label(CylRegionRole::EndCapMax, bt.mat);
                    for j in 0..spec.z_max_layers {
                        let kz0 = j as isize;
                        let kz1 = (j + 1) as isize;
                        let lo = bt.vs.map(|v| ensure_vertex(&mut out, v, 0, kz0));
                        let hi = bt.vs.map(|v| ensure_vertex(&mut out, v, 0, kz1));
                        let mut nodes = [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]];
                        orient_cell(CellType::Prism, &mut nodes, &out.points);
                        out.push_cell(CellType::Prism, &nodes, label);
                        if j + 1 == spec.z_max_layers {
                            if let Some(bl) = bt.bnd_label {
                                outer_face_label.insert(facet_key(&hi), bl);
                            }
                        }
                    }
                }
            }
            _ => {}
        }
    }

    // Handle outer labels on exposed side quads when ONLY radial OR ONLY axial is active:
    let mut rim_min_edges: Vec<RimEdgeItem> = Vec::new();
    let mut rim_max_edges: Vec<RimEdgeItem> = Vec::new();
    for (&ekey, b_list) in &barrel_edge_info {
        if let Some(&cap_bl) = cap_min_edge_label.get(&ekey) {
            let &(u, v, mat, bar_bl) = &b_list[0];
            rim_min_edges.push((u, v, mat, bar_bl, cap_bl));
        }
        if let Some(&cap_bl) = cap_max_edge_label.get(&ekey) {
            let &(u, v, mat, bar_bl) = &b_list[0];
            rim_max_edges.push((u, v, mat, bar_bl, cap_bl));
        }
    }

    if spec.radial_active() && !spec.z_min_active() {
        for &(u, v, _, _, cap_bl) in &rim_min_edges {
            if let Some(bl) = cap_bl.or(spec.outer_label) {
                for kr in 0..spec.r_layers {
                    let q = [
                        ensure_vertex(&mut out, u, kr, 0),
                        ensure_vertex(&mut out, v, kr, 0),
                        ensure_vertex(&mut out, v, kr + 1, 0),
                        ensure_vertex(&mut out, u, kr + 1, 0),
                    ];
                    outer_face_label.insert(facet_key(&q), bl);
                }
            }
        }
    }
    if spec.radial_active() && !spec.z_max_active() {
        for &(u, v, _, _, cap_bl) in &rim_max_edges {
            if let Some(bl) = cap_bl.or(spec.outer_label) {
                for kr in 0..spec.r_layers {
                    let q = [
                        ensure_vertex(&mut out, u, kr, 0),
                        ensure_vertex(&mut out, v, kr, 0),
                        ensure_vertex(&mut out, v, kr + 1, 0),
                        ensure_vertex(&mut out, u, kr + 1, 0),
                    ];
                    outer_face_label.insert(facet_key(&q), bl);
                }
            }
        }
    }
    if !spec.radial_active() && spec.z_min_active() {
        for &(u, v, _, bar_bl, _) in &rim_min_edges {
            if let Some(bl) = bar_bl.or(spec.outer_label) {
                for j in 0..spec.z_min_layers {
                    let kz0 = -(j as isize);
                    let kz1 = -((j + 1) as isize);
                    let q = [
                        ensure_vertex(&mut out, u, 0, kz0),
                        ensure_vertex(&mut out, v, 0, kz0),
                        ensure_vertex(&mut out, v, 0, kz1),
                        ensure_vertex(&mut out, u, 0, kz1),
                    ];
                    outer_face_label.insert(facet_key(&q), bl);
                }
            }
        }
    }
    if !spec.radial_active() && spec.z_max_active() {
        for &(u, v, _, bar_bl, _) in &rim_max_edges {
            if let Some(bl) = bar_bl.or(spec.outer_label) {
                for j in 0..spec.z_max_layers {
                    let kz0 = j as isize;
                    let kz1 = (j + 1) as isize;
                    let q = [
                        ensure_vertex(&mut out, u, 0, kz0),
                        ensure_vertex(&mut out, v, 0, kz0),
                        ensure_vertex(&mut out, v, 0, kz1),
                        ensure_vertex(&mut out, u, 0, kz1),
                    ];
                    outer_face_label.insert(facet_key(&q), bl);
                }
            }
        }
    }

    // 3. Rim circular edge chains -> `N_r x N_z` CornerAnnulus `Hex`es!
    if spec.radial_active() && spec.z_min_active() {
        for &(u, v, mat, bar_bl, cap_bl) in &rim_min_edges {
            let label = spec.label(CylRegionRole::CornerAnnulusMin, mat);
            for j in 0..spec.z_min_layers {
                let kz0 = -(j as isize);
                let kz1 = -((j + 1) as isize);
                for kr in 0..spec.r_layers {
                    let b0 = ensure_vertex(&mut out, u, kr, kz0);
                    let b1 = ensure_vertex(&mut out, u, kr + 1, kz0);
                    let b2 = ensure_vertex(&mut out, v, kr + 1, kz0);
                    let b3 = ensure_vertex(&mut out, v, kr, kz0);
                    let t0 = ensure_vertex(&mut out, u, kr, kz1);
                    let t1 = ensure_vertex(&mut out, u, kr + 1, kz1);
                    let t2 = ensure_vertex(&mut out, v, kr + 1, kz1);
                    let t3 = ensure_vertex(&mut out, v, kr, kz1);
                    let mut nodes = [b0, b1, b2, b3, t0, t1, t2, t3];
                    orient_cell(CellType::Hex, &mut nodes, &out.points);
                    out.push_cell(CellType::Hex, &nodes, label);
                    if kr + 1 == spec.r_layers {
                        if let Some(bl) = bar_bl.or(spec.outer_label) {
                            outer_face_label.insert(facet_key(&[b1, b2, t2, t1]), bl);
                        }
                    }
                    if j + 1 == spec.z_min_layers {
                        if let Some(bl) = cap_bl.or(spec.outer_label) {
                            outer_face_label.insert(facet_key(&[t0, t1, t2, t3]), bl);
                        }
                    }
                }
            }
        }
    }

    if spec.radial_active() && spec.z_max_active() {
        for &(u, v, mat, bar_bl, cap_bl) in &rim_max_edges {
            let label = spec.label(CylRegionRole::CornerAnnulusMax, mat);
            for j in 0..spec.z_max_layers {
                let kz0 = j as isize;
                let kz1 = (j + 1) as isize;
                for kr in 0..spec.r_layers {
                    let b0 = ensure_vertex(&mut out, u, kr, kz0);
                    let b1 = ensure_vertex(&mut out, u, kr + 1, kz0);
                    let b2 = ensure_vertex(&mut out, v, kr + 1, kz0);
                    let b3 = ensure_vertex(&mut out, v, kr, kz0);
                    let t0 = ensure_vertex(&mut out, u, kr, kz1);
                    let t1 = ensure_vertex(&mut out, u, kr + 1, kz1);
                    let t2 = ensure_vertex(&mut out, v, kr + 1, kz1);
                    let t3 = ensure_vertex(&mut out, v, kr, kz1);
                    let mut nodes = [b0, b1, b2, b3, t0, t1, t2, t3];
                    orient_cell(CellType::Hex, &mut nodes, &out.points);
                    out.push_cell(CellType::Hex, &nodes, label);
                    if kr + 1 == spec.r_layers {
                        if let Some(bl) = bar_bl.or(spec.outer_label) {
                            outer_face_label.insert(facet_key(&[b1, b2, t2, t1]), bl);
                        }
                    }
                    if j + 1 == spec.z_max_layers {
                        if let Some(bl) = cap_bl.or(spec.outer_label) {
                            outer_face_label.insert(facet_key(&[t0, t1, t2, t3]), bl);
                        }
                    }
                }
            }
        }
    }

    // 4. Radial extrusion of interior `Op::Sheet` edges touching `r = R_in` into conductive `Quad` facets!
    let mut extruded_sheets: Vec<([u32; 4], i32)> = Vec::new();
    if spec.radial_active() && !interior_sheets.is_empty() {
        // Collect all interior Quad faces of `out` (`count == 2`) for topology verification.
        let mut quad_counts: HashMap<[u32; 4], u32> = HashMap::new();
        for c in phys.n_cells()..out.n_cells() {
            let vs = out.cell(c);
            for f in out.types[c].faces() {
                if f.len() == 4 {
                    let q = canonical_quad_face([
                        vs[f[0] as usize],
                        vs[f[1] as usize],
                        vs[f[2] as usize],
                        vs[f[3] as usize],
                    ]);
                    *quad_counts.entry(q).or_insert(0) += 1;
                }
            }
        }

        let mut seen_sheet_edges: BTreeSet<(u32, u32, i32)> = BTreeSet::new();
        for &(tri, sheet_lab) in interior_sheets {
            let (p0, p1, p2) = (
                phys.points[tri[0] as usize],
                phys.points[tri[1] as usize],
                phys.points[tri[2] as usize],
            );
            let n_sheet = cross3(sub3(p1, p0), sub3(p2, p0));
            for (i, j) in [(0, 1), (1, 2), (2, 0)] {
                let (u, v) = (tri[i], tri[j]);
                let ekey = (u.min(v), u.max(v));
                if !barrel_edge_info.contains_key(&ekey) || !seen_sheet_edges.insert((ekey.0, ekey.1, sheet_lab)) {
                    continue;
                }
                for kr in 0..spec.r_layers {
                    let u0 = ensure_vertex(&mut out, u, kr, 0);
                    let v0 = ensure_vertex(&mut out, v, kr, 0);
                    let v1 = ensure_vertex(&mut out, v, kr + 1, 0);
                    let u1 = ensure_vertex(&mut out, u, kr + 1, 0);
                    let mut q = [u0, v0, v1, u1];
                    if quad_counts.get(&canonical_quad_face(q)).copied() == Some(2) {
                        let n_q = cross3(
                            sub3(out.points[q[1] as usize], out.points[q[0] as usize]),
                            sub3(out.points[q[3] as usize], out.points[q[0] as usize]),
                        );
                        if dot3(n_q, n_sheet) < 0.0 {
                            q = [v0, u0, u1, v1];
                        }
                        extruded_sheets.push((q, sheet_lab));
                    }
                }
            }
        }
    }

    // 5. Exterior boundary facets of `out`.
    for f in out.boundary_faces_ordered() {
        let fk = facet_key(&f);
        let all_phys = f.iter().all(|&v| (v as usize) < phys.n_points());
        let bl_opt = if all_phys {
            phys_facet_label.get(&fk).copied().or_else(|| outer_face_label.get(&fk).copied())
        } else {
            outer_face_label.get(&fk).copied().or(spec.outer_label)
        };
        if let Some(bl) = bl_opt {
            out.push_facet(&f, bl);
        }
    }

    // 6. Lockstep cylindrical curving (`edge_nodes` + barrel `quad_face_nodes`).
    let mut edge_node_map: BTreeMap<[u32; 2], [f64; 3]> = phys_edge_nodes
        .iter()
        .map(|&([a, b], p)| ([a.min(b), a.max(b)], p))
        .collect();

    if spec.curved {
        // Precompute interior physical tets touching barrel edges so `k = 0` never inverts an interior cut-tet.
        let mut barrel_edge_r0_mid: BTreeMap<(u32, u32), [f64; 3]> = BTreeMap::new();
        let mut uncurved_barrel_cand: Vec<([u32; 2], [f64; 3])> = Vec::new();

        for &(u, v) in barrel_edge_info.keys() {
            let key = [u, v];
            let (pu, pv) = (phys.points[u as usize], phys.points[v as usize]);
            let (_, ru_vec, _) = decompose(pu);
            let (_, rv_vec, _) = decompose(pv);
            let chord_perp = norm3(sub3(rv_vec, ru_vec));
            if chord_perp <= tol_abs {
                continue;
            }
            let sum_r = [0.5 * (ru_vec[0] + rv_vec[0]), 0.5 * (ru_vec[1] + rv_vec[1]), 0.5 * (ru_vec[2] + rv_vec[2])];
            let sum_nrm = norm3(sum_r);
            if sum_nrm <= 1e-30 {
                continue;
            }
            if let Some(&q_existing) = edge_node_map.get(&key) {
                let (_, r_vec_exist, _) = decompose(q_existing);
                barrel_edge_r0_mid.insert((u, v), r_vec_exist);
            } else {
                let z_mid = 0.5 * (dot3(pu, a_hat) + dot3(pv, a_hat));
                let c_ax = axis_pt(z_mid);
                let r_arc = [
                    spec.r_in * sum_r[0] / sum_nrm,
                    spec.r_in * sum_r[1] / sum_nrm,
                    spec.r_in * sum_r[2] / sum_nrm,
                ];
                let off = norm3(sub3(r_arc, sum_r));
                if off > 1e-9 * chord_perp {
                    uncurved_barrel_cand.push((key, [c_ax[0] + r_arc[0], c_ax[1] + r_arc[1], c_ax[2] + r_arc[2]]));
                }
            }
        }

        if !uncurved_barrel_cand.is_empty() {
            // Relax candidate `k = 0` barrel midpoints against physical interior `Tet`s ONLY (if any fail `MIN_SCALED_JACOBIAN`).
            let phys_tets: Vec<[u32; 4]> = (0..phys.n_cells())
                .map(|c| {
                    let cv = phys.cell(c);
                    [cv[0], cv[1], cv[2], cv[3]]
                })
                .collect();
            let _ = relax_against_phys_tets(&phys.points, &phys_tets, &edge_node_map, &mut uncurved_barrel_cand);
            for ([u, v], q0) in uncurved_barrel_cand {
                edge_node_map.insert([u, v], q0);
                let (_, r_vec0, _) = decompose(q0);
                barrel_edge_r0_mid.insert((u, v), r_vec0);
            }
        }

        let mut quad_face_map: BTreeMap<[u32; 4], [f64; 3]> = BTreeMap::new();
        let rim_min_set: BTreeSet<(u32, u32)> = rim_min_edges.iter().map(|&(u, v, ..)| (u.min(v), u.max(v))).collect();
        let rim_max_set: BTreeSet<(u32, u32)> = rim_max_edges.iter().map(|&(u, v, ..)| (u.min(v), u.max(v))).collect();

        // Lockstep radial propagation across `kr = 0..=N_r` (and axial `kz` on rim chains)!
        for (&(u, v), &r_vec0) in &barrel_edge_r0_mid {
            let max_kr = if spec.radial_active() { spec.r_layers } else { 0 };
            for kr in 0..=max_kr {
                let scale = (spec.r_in + r_off[kr]) / spec.r_in;
                let r_vec_k = [r_vec0[0] * scale, r_vec0[1] * scale, r_vec0[2] * scale];

                // Active axial slab `kz = 0`:
                if let (Some(&uk), Some(&vk)) = (shifted_cache.get(&(u, kr, 0)), shifted_cache.get(&(v, kr, 0))) {
                    let z_mid = 0.5 * (dot3(out.points[uk as usize], a_hat) + dot3(out.points[vk as usize], a_hat));
                    let c_ax = axis_pt(z_mid);
                    edge_node_map.insert(
                        [uk.min(vk), uk.max(vk)],
                        [c_ax[0] + r_vec_k[0], c_ax[1] + r_vec_k[1], c_ax[2] + r_vec_k[2]],
                    );
                }

                // Bottom rim `kz = -1 ..= -N_{z-}`:
                if spec.z_min_active() && rim_min_set.contains(&(u, v)) {
                    for j in 1..=spec.z_min_layers {
                        let kz = -(j as isize);
                        if let (Some(&uk), Some(&vk)) = (shifted_cache.get(&(u, kr, kz)), shifted_cache.get(&(v, kr, kz)))
                        {
                            let z_val = z_min - z_min_off[j];
                            let c_ax = axis_pt(z_val);
                            edge_node_map.insert(
                                [uk.min(vk), uk.max(vk)],
                                [c_ax[0] + r_vec_k[0], c_ax[1] + r_vec_k[1], c_ax[2] + r_vec_k[2]],
                            );
                        }
                        // Curved barrel `Quad` face on `[kz_prev, kz]` at radius `r_k`!
                        let kz_prev = -((j - 1) as isize);
                        if let (Some(&u_a), Some(&v_a), Some(&u_b), Some(&v_b)) = (
                            shifted_cache.get(&(u, kr, kz_prev)),
                            shifted_cache.get(&(v, kr, kz_prev)),
                            shifted_cache.get(&(u, kr, kz)),
                            shifted_cache.get(&(v, kr, kz)),
                        ) {
                            let z_mid = z_min - 0.5 * (z_min_off[j - 1] + z_min_off[j]);
                            let c_ax = axis_pt(z_mid);
                            quad_face_map.insert(
                                canonical_quad_face([u_a, v_a, v_b, u_b]),
                                [c_ax[0] + r_vec_k[0], c_ax[1] + r_vec_k[1], c_ax[2] + r_vec_k[2]],
                            );
                        }
                    }
                }

                // Top rim `kz = 1 ..= N_{z+}`:
                if spec.z_max_active() && rim_max_set.contains(&(u, v)) {
                    for j in 1..=spec.z_max_layers {
                        let kz = j as isize;
                        if let (Some(&uk), Some(&vk)) = (shifted_cache.get(&(u, kr, kz)), shifted_cache.get(&(v, kr, kz)))
                        {
                            let z_val = z_max + z_max_off[j];
                            let c_ax = axis_pt(z_val);
                            edge_node_map.insert(
                                [uk.min(vk), uk.max(vk)],
                                [c_ax[0] + r_vec_k[0], c_ax[1] + r_vec_k[1], c_ax[2] + r_vec_k[2]],
                            );
                        }
                        // Curved barrel `Quad` face on `[kz_prev, kz]` at radius `r_k`!
                        let kz_prev = (j - 1) as isize;
                        if let (Some(&u_a), Some(&v_a), Some(&u_b), Some(&v_b)) = (
                            shifted_cache.get(&(u, kr, kz_prev)),
                            shifted_cache.get(&(v, kr, kz_prev)),
                            shifted_cache.get(&(u, kr, kz)),
                            shifted_cache.get(&(v, kr, kz)),
                        ) {
                            let z_mid = z_max + 0.5 * (z_max_off[j - 1] + z_max_off[j]);
                            let c_ax = axis_pt(z_mid);
                            quad_face_map.insert(
                                canonical_quad_face([u_a, v_a, v_b, u_b]),
                                [c_ax[0] + r_vec_k[0], c_ax[1] + r_vec_k[1], c_ax[2] + r_vec_k[2]],
                            );
                        }
                    }
                }
            }
        }

        out.set_quad_face_nodes(quad_face_map);
    }

    if !edge_node_map.is_empty() {
        out.set_edge_nodes(edge_node_map);
    }

    Ok(CylPmlMesh3d {
        mesh: out,
        n_physical_cells: phys.n_cells(),
        n_physical_points: phys.n_points(),
        extruded_sheets,
        spec: spec.clone(),
    })
}

/// Relaxes newly proposed `k = 0` interface edge midpoints against incident physical `Tet`s only,
/// keeping any pre-existing interior `fixed_nodes` fixed.
pub(crate) fn relax_against_phys_tets(
    pts: &[[f64; 3]],
    phys_tets: &[[u32; 4]],
    fixed_nodes: &BTreeMap<[u32; 2], [f64; 3]>,
    cand: &mut Vec<([u32; 2], [f64; 3])>,
) -> usize {
    const TET_E: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    const MAX_HALVINGS: u8 = 4;
    if cand.is_empty() || phys_tets.is_empty() {
        return 0;
    }
    cand.sort_by_key(|x| x.0);
    cand.dedup_by_key(|x| x.0);
    let idx_map: HashMap<[u32; 2], usize> = cand.iter().enumerate().map(|(i, &(e, _))| (e, i)).collect();
    let mut level = vec![0u8; cand.len()];

    let touched_tets: Vec<(usize, [Option<usize>; 6])> = phys_tets
        .iter()
        .enumerate()
        .filter_map(|(ti, t)| {
            let mut ids = [None; 6];
            let mut any = false;
            for (k, &(i, j)) in TET_E.iter().enumerate() {
                let pr = [t[i].min(t[j]), t[i].max(t[j])];
                if let Some(&ci) = idx_map.get(&pr) {
                    ids[k] = Some(ci);
                    any = true;
                }
            }
            any.then_some((ti, ids))
        })
        .collect();
    if touched_tets.is_empty() {
        return 0;
    }

    let rule = quad_rule_collapsed(6);
    let mut samples: Vec<[f64; 4]> = vec![
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    for &(i, j) in &TET_E {
        let mut l = [0.0; 4];
        l[i] = 0.5;
        l[j] = 0.5;
        samples.push(l);
    }
    samples.extend(rule.lam.iter().map(|l| [l[0], l[1], l[2], l[3]]));

    let scale_of = |lv: u8| if lv > MAX_HALVINGS { 0.0 } else { 0.5_f64.powi(lv as i32) };
    let tet_ratio = |ti: usize, ids: &[Option<usize>; 6], level: &[u8]| -> f64 {
        let t = phys_tets[ti];
        let v = [pts[t[0] as usize], pts[t[1] as usize], pts[t[2] as usize], pts[t[3] as usize]];
        let mut x = CurvedTet::straight(v).x;
        for (k, &(i, j)) in TET_E.iter().enumerate() {
            let pr = [t[i].min(t[j]), t[i].max(t[j])];
            if let Some(ci) = ids[k] {
                let s = scale_of(level[ci]);
                for d in 0..3 {
                    let mid = 0.5 * (v[i][d] + v[j][d]);
                    x[4 + k][d] = mid + s * (cand[ci].1[d] - mid);
                }
            } else if let Some(&q_fix) = fixed_nodes.get(&pr) {
                x[4 + k] = q_fix;
            }
        }
        let ct = CurvedTet::new(x);
        let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
        if det_s == 0.0 || !det_s.is_finite() {
            return 1.0;
        }
        samples
            .iter()
            .map(|l| ct.jacobian(l).1 / det_s)
            .fold(f64::INFINITY, f64::min)
    };

    let mut relaxed_count = 0usize;
    for _ in 0..=(MAX_HALVINGS as usize + 1) {
        let mut bump = BTreeSet::new();
        for &(ti, ref ids) in &touched_tets {
            if tet_ratio(ti, ids, &level) <= crate::meshing::curve3::MIN_SCALED_JACOBIAN {
                for ci in ids.iter().flatten() {
                    if level[*ci] <= MAX_HALVINGS {
                        bump.insert(*ci);
                    }
                }
            }
        }
        if bump.is_empty() {
            break;
        }
        for ci in bump {
            level[ci] += 1;
        }
    }

    let mut kept = Vec::with_capacity(cand.len());
    for (ci, &(pr, q)) in cand.iter().enumerate() {
        let s = scale_of(level[ci]);
        if level[ci] > 0 {
            relaxed_count += 1;
        }
        if s > 0.0 {
            let (pa, pb) = (pts[pr[0] as usize], pts[pr[1] as usize]);
            let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            kept.push((
                pr,
                [
                    mid[0] + s * (q[0] - mid[0]),
                    mid[1] + s * (q[1] - mid[1]),
                    mid[2] + s * (q[2] - mid[2]),
                ],
            ));
        }
    }
    *cand = kept;
    relaxed_count
}
