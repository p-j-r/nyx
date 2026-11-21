// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-Rust 3D P2 slicer (`Slice3d`), convex-box cutaway (`Cutaway3d`), marching-tetrahedra isosurfaces (`Isosurface3d`), region mesh view
//! (`Mesh3d`), and zero-copy Bloch tiler for `nyx-analysis::figure`.
//!
//! Every 10-node P2 tetrahedron (`v0..v3, e01, e02, e03, e12, e13, e23`) is subdivided into eight conformal positive-volume linear sub-tetrahedra
//! ([`P2_SUBTETS`]) and evaluated directly on per-cell un-averaged `CellField` entries (`cf.re[c * nloc + a]`, `cf.im[c * nloc + a]`) without
//! cross-cell vertex welding, preserving sharp normal-field discontinuities across material interfaces.

use crate::canvas::{Anchor, Canvas, Item, Rgba, INK, MUTED};
use crate::figure::{
    bytes, clip, draw_ticks, fit_aspect, length_unit, mode_lines, si, Cmap, Figure, FigureOptions, Frame, Kind,
    ModeNote, OUTLINE, PANEL_BG,
};
use crate::raster::{raster_field, raster_field_3d, text_width, ColouredTri, Tri3d};
use femcore::fieldexpr::{Fun, Op, Tree};
use femcore::formula::{self, Named, Names};
use femcore::linalg::csr::C;
use femcore::mesh::SimplexMesh;
use nyx_physics::{Error, FieldFrame, Scalar, Solution, SolveResult};
use std::collections::HashMap;

/// Exact 8-subtet P2 split in nyx local node order (`v0..v3, e01, e02, e03, e12, e13, e23`).
/// Every sub-tetrahedron has `+1/8` the signed volume of the parent straight tetrahedron.
pub const P2_SUBTETS: [[usize; 4]; 8] =
    [[0, 4, 5, 6], [4, 1, 7, 8], [5, 7, 2, 9], [6, 8, 9, 3], [4, 5, 6, 9], [4, 6, 8, 9], [4, 8, 7, 9], [4, 7, 5, 9]];

const P1_SUBTETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

/// Local edges of a 4-node linear sub-tetrahedron.
const SUBTET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// Faces of a parent tetrahedron `(v_a, v_b, v_c, opp_v)` and their three P2 edge-midpoint indices `(e_ab, e_ac, e_bc)`.
const TET_FACES: [([usize; 3], usize, [usize; 3]); 4] =
    [([1, 2, 3], 0, [7, 8, 9]), ([0, 2, 3], 1, [5, 6, 9]), ([0, 1, 3], 2, [4, 6, 8]), ([0, 1, 2], 3, [4, 5, 7])];

/// Slice plane specification for `Kind::Slice3d`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlaneSpec {
    Axis { axis: char, at: f64 },
    General { origin: [f64; 3], normal: [f64; 3] },
}

/// Optional axis-aligned cutaway bounds in metres; unset bounds stay open.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ClipBox {
    pub xmin: Option<f64>,
    pub xmax: Option<f64>,
    pub ymin: Option<f64>,
    pub ymax: Option<f64>,
    pub zmin: Option<f64>,
    pub zmax: Option<f64>,
}

/// Edge overlay mode for 3D figures.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum EdgeMode {
    None,
    #[default]
    Regions,
    Mesh,
}

impl EdgeMode {
    pub fn parse(s: &str) -> Option<EdgeMode> {
        match s.to_ascii_lowercase().as_str() {
            "none" => Some(EdgeMode::None),
            "regions" => Some(EdgeMode::Regions),
            "mesh" => Some(EdgeMode::Mesh),
            _ => None,
        }
    }
}

/// Orthographic 3D turntable camera in the wafer frame (`z` up).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera3d {
    pub azimuth_deg: f64,
    pub elevation_deg: f64,
}

impl Default for Camera3d {
    fn default() -> Self {
        Camera3d { azimuth_deg: -52.0, elevation_deg: 24.0 }
    }
}

impl Camera3d {
    /// Orthonormal `(right, up, view_dir)` basis where `view_dir` points from the camera into the scene (larger `dot(p, view_dir)` is farther).
    pub fn basis(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let az = self.azimuth_deg.to_radians();
        let el = self.elevation_deg.to_radians();
        let (sa, ca) = az.sin_cos();
        let (se, ce) = el.sin_cos();
        let right = [-sa, ca, 0.0];
        let up = [-ca * se, -sa * se, ce];
        let view_dir = [-ca * ce, -sa * ce, -se];
        (right, up, view_dir)
    }
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n > 1e-30 {
        [a[0] / n, a[1] / n, a[2] / n]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Signed volume `det([b - a, c - a, d - a]) / 6` of a linear tetrahedron `[a, b, c, d]`.
pub fn tet_signed_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    dot3(sub3(b, a), cross3(sub3(c, a), sub3(d, a))) / 6.0
}

/// Decodes a PNG byte slice into `(width, height, rgba_bytes)` for verification and pixel inspection.
pub fn decode_png(bytes: &[u8]) -> Result<(u32, u32, Vec<u8>), String> {
    let decoder = png::Decoder::new(std::io::Cursor::new(bytes));
    let mut reader = decoder.read_info().map_err(|e| format!("PNG header: {e}"))?;
    let size = reader.output_buffer_size().ok_or_else(|| "PNG buffer size overflow".to_string())?;
    let mut buf = vec![0u8; size];
    let info = reader.next_frame(&mut buf).map_err(|e| format!("PNG frame: {e}"))?;
    buf.truncate(info.buffer_size());
    Ok((info.width, info.height, buf))
}

// ---------------------------------------------------------------------------------------------------------- 3D mesh topology

#[derive(Clone, Copy)]
struct FeatureSeg {
    cell: u32,
    la: u8,
    lb: u8,
    lm: u8,
}

struct Mesh3dTopo {
    lo: [f64; 3],
    hi: [f64; 3],
    geom10: Vec<[[f64; 3]; 10]>,
    /// For `c * 4 + f`, `u32::MAX` if exterior face, else neighbour cell index `c2`.
    face_neighbor: Vec<u32>,
    /// Exterior feature edges (sharp dihedral, domain border, or facet-label border) on the unclipped surface.
    ext_feature_edges: Vec<FeatureSeg>,
}

impl Mesh3dTopo {
    fn build(mesh: &SimplexMesh) -> Mesh3dTopo {
        let nc = mesh.n_cells();
        let (lo, hi) = mesh.bbox();
        let curved = mesh.has_curved();
        let cell_ids: Vec<usize> = (0..nc).collect();
        let geom_chunks = femcore::par::map_chunks(&cell_ids, 512, |chunk| {
            let mut out = Vec::with_capacity(chunk.len());
            for &c in chunk {
                if curved {
                    out.push(mesh.cell_geometry_nodes(c));
                } else {
                    let v = mesh.cell(c);
                    let p = |i: usize| -> [f64; 3] {
                        let b = v[i] as usize * 3;
                        [mesh.points[b], mesh.points[b + 1], mesh.points[b + 2]]
                    };
                    let (p0, p1, p2, p3) = (p(0), p(1), p(2), p(3));
                    let mid =
                        |a: [f64; 3], b: [f64; 3]| [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])];
                    out.push([
                        p0,
                        p1,
                        p2,
                        p3,
                        mid(p0, p1),
                        mid(p0, p2),
                        mid(p0, p3),
                        mid(p1, p2),
                        mid(p1, p3),
                        mid(p2, p3),
                    ]);
                }
            }
            out
        });
        let mut geom10 = Vec::with_capacity(nc);
        for g in geom_chunks {
            geom10.extend_from_slice(&g);
        }

        let sort3 = |mut k: [u32; 3]| -> [u32; 3] {
            if k[0] > k[1] {
                k.swap(0, 1);
            }
            if k[1] > k[2] {
                k.swap(1, 2);
            }
            if k[0] > k[1] {
                k.swap(0, 1);
            }
            k
        };

        let mut half_faces: Vec<([u32; 3], u32)> = Vec::with_capacity(nc * 4);
        for c in 0..nc {
            let v = mesh.cell(c);
            for (f, &(fv, _, _)) in TET_FACES.iter().enumerate() {
                half_faces.push((sort3([v[fv[0]], v[fv[1]], v[fv[2]]]), (c * 4 + f) as u32));
            }
        }
        half_faces.sort_unstable_by_key(|x| x.0);

        let mut face_neighbor = vec![u32::MAX; nc * 4];
        let mut i = 0;
        while i < half_faces.len() {
            if i + 1 < half_faces.len() && half_faces[i].0 == half_faces[i + 1].0 {
                let (cf0, cf1) = (half_faces[i].1 as usize, half_faces[i + 1].1 as usize);
                face_neighbor[cf0] = (cf1 / 4) as u32;
                face_neighbor[cf1] = (cf0 / 4) as u32;
                i += 2;
            } else {
                i += 1;
            }
        }

        // Collect exterior half-edges to find sharp exterior creases / region borders
        struct ExtHalfEdge {
            key: (u32, u32),
            cell: u32,
            la: u8,
            lb: u8,
            lm: u8,
            dom: i32,
            normal: [f64; 3],
        }
        let mut ext_edges: Vec<ExtHalfEdge> = Vec::new();
        for c in 0..nc {
            let v = mesh.cell(c);
            let dom = mesh.cell_labels[c];
            let g = &geom10[c];
            for (f, &(fv, opp, mids)) in TET_FACES.iter().enumerate() {
                if face_neighbor[c * 4 + f] != u32::MAX {
                    continue;
                }
                let mut n = cross3(sub3(g[fv[1]], g[fv[0]]), sub3(g[fv[2]], g[fv[0]]));
                if dot3(n, sub3(g[opp], g[fv[0]])) > 0.0 {
                    n = scale3(n, -1.0);
                }
                let n = normalize3(n);
                for &(ea, eb, em) in &[(fv[0], fv[1], mids[0]), (fv[0], fv[2], mids[1]), (fv[1], fv[2], mids[2])] {
                    let (va, vb) = (v[ea], v[eb]);
                    ext_edges.push(ExtHalfEdge {
                        key: (va.min(vb), va.max(vb)),
                        cell: c as u32,
                        la: ea as u8,
                        lb: eb as u8,
                        lm: em as u8,
                        dom,
                        normal: n,
                    });
                }
            }
        }
        ext_edges.sort_unstable_by_key(|e| e.key);
        let mut ext_feature_edges = Vec::new();
        let mut j = 0;
        while j < ext_edges.len() {
            let mut k = j + 1;
            while k < ext_edges.len() && ext_edges[k].key == ext_edges[j].key {
                k += 1;
            }
            let e0 = &ext_edges[j];
            let sharp = if k - j != 2 {
                true
            } else {
                let e1 = &ext_edges[j + 1];
                e0.dom != e1.dom || dot3(e0.normal, e1.normal) < 0.906
            };
            if sharp {
                ext_feature_edges.push(FeatureSeg { cell: e0.cell, la: e0.la, lb: e0.lb, lm: e0.lm });
            }
            j = k;
        }

        Mesh3dTopo { lo, hi, geom10, face_neighbor, ext_feature_edges }
    }
}

// ----------------------------------------------------------------------------------------------- mode phase alignment & formula

#[derive(Clone)]
struct Fam3d {
    prefix: String,
    axes: Vec<char>,
    fields: Vec<usize>,
}

fn families_3d(sol: &Solution) -> Vec<Fam3d> {
    let mut out: Vec<Fam3d> = Vec::new();
    for (i, f) in sol.cell_fields.iter().enumerate() {
        let mut ch = f.name.chars();
        let (Some(axis), prefix) = (ch.next_back(), ch.as_str()) else { continue };
        if !matches!(axis, 'x' | 'y' | 'z') || prefix.is_empty() {
            continue;
        }
        match out.iter_mut().find(|v| v.prefix == prefix) {
            Some(v) => {
                v.axes.push(axis);
                v.fields.push(i);
            }
            None => out.push(Fam3d { prefix: prefix.to_string(), axes: vec![axis], fields: vec![i] }),
        }
    }
    out
}

/// Prepared field evaluator for one solution: keeps zero-copy references to `Solution::cell_fields` (or nodal fallback) together with the
/// mode's phase rotation `rot` and per-field normalization factor `scale`.
struct EvalMode<'a> {
    note: ModeNote,
    nloc: usize,
    order: usize,
    tree: Tree,
    expr_label: String,
    /// Per leaf slot >= 3: `(cell_field_index, phase_and_scale_multiplier)`; slots 0..3 are `(x, y, z)`.
    field_slots: Vec<(&'a [f64], &'a [f64], C)>,
    fallback_re: Vec<Vec<f64>>,
    periodic_phase: f64,
    vmin: f64,
    vmax: f64,
    signed: bool,
    abs_max: f64,
}

struct FormulaBinder<'a> {
    sol: &'a Solution,
    fams: &'a [Fam3d],
    field_map: HashMap<String, usize>,
    slots: Vec<usize>,
}

impl<'a> FormulaBinder<'a> {
    fn new(sol: &'a Solution, fams: &'a [Fam3d]) -> Self {
        let mut field_map = HashMap::new();
        for (i, cf) in sol.cell_fields.iter().enumerate() {
            field_map.insert(cf.name.clone(), i);
        }
        Self { sol, fams, field_map, slots: Vec::new() }
    }

    fn leaf_for_cf(&mut self, cf_idx: usize) -> Tree {
        if let Some(pos) = self.slots.iter().position(|&x| x == cf_idx) {
            Tree::Leaf(3 + pos)
        } else {
            let pos = self.slots.len();
            self.slots.push(cf_idx);
            Tree::Leaf(3 + pos)
        }
    }
}

impl Names for FormulaBinder<'_> {
    fn name(&mut self, name: &str) -> Result<Option<Named>, String> {
        match name {
            "x" => return Ok(Some(Named::Scalar(Tree::Leaf(0)))),
            "y" => return Ok(Some(Named::Scalar(Tree::Leaf(1)))),
            "z" => return Ok(Some(Named::Scalar(Tree::Leaf(2)))),
            _ => {}
        }
        if let Some(&cf_idx) = self.field_map.get(name) {
            return Ok(Some(Named::Scalar(self.leaf_for_cf(cf_idx))));
        }
        // Vector family name (`E`, `H`, `u`, ...): bind present components and zero-fill missing Cartesian axes when a family is partial
        if let Some(fam) = self.fams.iter().find(|f| f.prefix == name) {
            let mut comps = Vec::with_capacity(3);
            for (axis_i, axis_ch) in ['x', 'y', 'z'].iter().enumerate() {
                if let Some(p) = fam.axes.iter().position(|a| a == axis_ch) {
                    let cf_idx = fam.fields[p];
                    comps.push((axis_i, self.leaf_for_cf(cf_idx)));
                } else {
                    comps.push((axis_i, Tree::num(0.0)));
                }
            }
            return Ok(Some(Named::Vector(comps)));
        }
        // Component of a present family whose orthogonal axis was omitted from a synthetic solution (e.g. `Ey` when `Ex` exists)
        if name.len() >= 2 {
            let mut ch = name.chars();
            if let Some(axis) = ch.next_back() {
                let prefix = ch.as_str();
                if matches!(axis, 'x' | 'y' | 'z') && self.fams.iter().any(|f| f.prefix == prefix) {
                    return Ok(Some(Named::Scalar(Tree::num(0.0))));
                }
            }
        }
        Ok(None)
    }

    fn material(&mut self, name: &str) -> Result<Tree, String> {
        Err(format!("material('{name}') is not available in figure expressions"))
    }

    fn region(&mut self, name: &str) -> Result<Tree, String> {
        Err(format!("region('{name}') is not available in figure expressions"))
    }

    fn known(&self) -> String {
        let mut k = vec!["x".to_string(), "y".to_string(), "z".to_string()];
        for cf in &self.sol.cell_fields {
            k.push(cf.name.clone());
        }
        for f in self.fams {
            k.push(f.prefix.clone());
        }
        k.join(", ")
    }
}

/// Rewrites `abs(<vec>)` -> `norm(<vec>)` and `abs2(<vec>)` -> `norm2(<vec>)` when `<vec>` evaluates to a vector in `formula::parse`.
fn rewrite_vector_abs(text: &str, sol: &Solution, fams: &[Fam3d]) -> String {
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len() + 8);
    let mut i = 0;
    while i < bytes.len() {
        let prev_ident = i > 0 && (bytes[i - 1].is_ascii_alphanumeric() || bytes[i - 1] == b'_');
        if !prev_ident && (bytes[i..].starts_with(b"abs2") || bytes[i..].starts_with(b"abs")) {
            let fn_len = if bytes[i..].starts_with(b"abs2")
                && !bytes.get(i + 4).is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
            {
                4
            } else if bytes[i..].starts_with(b"abs")
                && !bytes.get(i + 3).is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
            {
                3
            } else {
                0
            };
            if fn_len > 0 {
                let mut p = i + fn_len;
                while p < bytes.len() && (bytes[p] as char).is_whitespace() {
                    p += 1;
                }
                if p < bytes.len() && bytes[p] == b'(' {
                    // Find matching ')'
                    let mut depth = 1usize;
                    let mut q = p + 1;
                    let mut in_quote: Option<u8> = None;
                    while q < bytes.len() && depth > 0 {
                        let ch = bytes[q];
                        if let Some(qt) = in_quote {
                            if ch == qt {
                                in_quote = None;
                            }
                        } else if ch == b'"' || ch == b'\'' || ch == b'`' {
                            in_quote = Some(ch);
                        } else if ch == b'(' {
                            depth += 1;
                        } else if ch == b')' {
                            depth -= 1;
                        }
                        q += 1;
                    }
                    if depth == 0 {
                        let inner_raw = &text[p + 1..q - 1];
                        let inner = rewrite_vector_abs(inner_raw, sol, fams);
                        let mut probe = FormulaBinder::new(sol, fams);
                        let is_vec = match formula::parse(&inner, &mut probe) {
                            Err(msg) => msg.contains("is a vector"),
                            Ok(_) => false,
                        };
                        if is_vec {
                            out.push_str(if fn_len == 4 { "norm2(" } else { "norm(" });
                        } else {
                            out.push_str(&text[i..i + fn_len]);
                            out.push('(');
                        }
                        out.push_str(&inner);
                        out.push(')');
                        i = q;
                        continue;
                    }
                }
            }
        }
        let ch = text[i..].chars().next().unwrap();
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn prepare_mode<'a>(
    mesh: &SimplexMesh,
    topo: &Mesh3dTopo,
    result: &'a SolveResult,
    k: usize,
    o: &FigureOptions,
    shift: [f64; 3],
    repeat: usize,
    mode_periodic_phase: f64,
) -> Result<EvalMode<'a>, Error> {
    let sol = &result.solutions[k];
    let order = if result.order == 2 { 2 } else { 1 };
    let nloc = mesh.nloc(order);
    let nc = mesh.n_cells();

    // Optional fallback when a scalar 3D solve (e.g. dfs V/T) stores nodal `fields` without `cell_fields`
    let mut fallback_re: Vec<Vec<f64>> = Vec::new();
    let mut fallback_names: Vec<String> = Vec::new();
    if sol.cell_fields.is_empty() {
        let flat = sol.flat_fields();
        let mut dofs = Vec::new();
        for (fname, _, vals) in flat {
            if vals.len() == mesh.n_dofs(order) {
                let mut cell_vals = Vec::with_capacity(nc * nloc);
                for c in 0..nc {
                    mesh.cell_dofs(c, order, &mut dofs);
                    for &d in &dofs {
                        cell_vals.push(vals[d]);
                    }
                }
                fallback_names.push(fname);
                fallback_re.push(cell_vals);
            }
        }
    }

    let fams = families_3d(sol);
    let frame = FieldFrame::of(result);
    let prop_axis = frame.as_ref().map(|f| ['x', 'y', 'z'][f.propagation]).unwrap_or('y');

    // Core cell mask: `o.core` if explicitly requested, else non-PML physical cells (falling back to all cells if all are marked outside)
    let mut core_cell = vec![false; nc];
    if let Some(core_doms) = &o.core {
        for (c, slot) in core_cell.iter_mut().enumerate() {
            *slot = core_doms.contains(&mesh.cell_labels[c]);
        }
        if !core_cell.iter().any(|&b| b) {
            return Err(Error::invalid("no cell of the 3D mesh carries the domains named as the core of the mode"));
        }
    } else {
        let outside = frame.as_ref().map(|f| f.outside.as_slice()).unwrap_or(&[]);
        for (c, slot) in core_cell.iter_mut().enumerate() {
            *slot = !outside.contains(&mesh.cell_labels[c]);
        }
        if !core_cell.iter().any(|&b| b) {
            core_cell.fill(true);
        }
    }

    // Validate cell_field lengths
    for cf in &sol.cell_fields {
        if cf.re.len() < nc * nloc || (!cf.im.is_empty() && cf.im.len() < nc * nloc) {
            return Err(Error::invalid(format!(
                "cell field '{}' has {} entries, fewer than {} cells x {} nodes",
                cf.name,
                cf.re.len(),
                nc,
                nloc
            )));
        }
    }

    let (field_name, theta, flipped, peak, quadrature, complex, fam_peaks): (
        String,
        f64,
        bool,
        f64,
        Option<f64>,
        bool,
        HashMap<usize, f64>,
    ) = if !fams.is_empty() {
        let primary = match &o.family {
            Some(p) => fams
                .iter()
                .find(|f| &f.prefix == p)
                .ok_or_else(|| Error::invalid(format!("the solution ('{}') has no vector family '{p}'", sol.label)))?,
            None => &fams[0],
        };
        let complex = primary.fields.iter().any(|&fi| !sol.cell_fields[fi].im.is_empty());
        let along = if primary.axes.len() >= 2 { primary.axes.iter().position(|&a| a == prop_axis) } else { None };
        let in_plane: Vec<usize> = (0..primary.fields.len()).filter(|&i| Some(i) != along).collect();

        // 1. Phase alignment over un-averaged core cell nodes
        let mut theta = 0.0f64;
        if complex {
            let (mut sr, mut si_acc) = (0.0f64, 0.0f64);
            for c in (0..nc).filter(|&c| core_cell[c]) {
                let base = c * nloc;
                for a in 0..nloc {
                    let idx = base + a;
                    for &ci in &in_plane {
                        let cf = &sol.cell_fields[primary.fields[ci]];
                        let re = cf.re[idx];
                        let im = if cf.im.is_empty() { 0.0 } else { cf.im[idx] };
                        sr += re * re - im * im;
                        si_acc += 2.0 * re * im;
                    }
                }
            }
            if sr.hypot(si_acc) > 1e-300 {
                theta = -0.5 * si_acc.atan2(sr);
            }
        }
        let (mut ct, mut st) = (theta.cos(), theta.sin());

        // 2. Sign alignment: largest in-plane component at node of largest in-plane amplitude is positive
        let mut best_amp = -1.0f64;
        let mut best_idx = 0usize;
        for c in (0..nc).filter(|&c| core_cell[c]) {
            let base = c * nloc;
            for a in 0..nloc {
                let idx = base + a;
                let mut amp = 0.0f64;
                for &ci in &in_plane {
                    let cf = &sol.cell_fields[primary.fields[ci]];
                    let re = cf.re[idx];
                    let im = if cf.im.is_empty() { 0.0 } else { cf.im[idx] };
                    amp += re * re + im * im;
                }
                if amp > best_amp {
                    best_amp = amp;
                    best_idx = idx;
                }
            }
        }
        let mut lead_re = 0.0f64;
        let mut lead_abs = -1.0f64;
        for &ci in &in_plane {
            let cf = &sol.cell_fields[primary.fields[ci]];
            let re = cf.re[best_idx];
            let im = if cf.im.is_empty() { 0.0 } else { cf.im[best_idx] };
            let r_rot = re * ct - im * st;
            if r_rot.abs() > lead_abs {
                lead_abs = r_rot.abs();
                lead_re = r_rot;
            }
        }
        let mut flipped = false;
        if lead_re < 0.0 {
            flipped = true;
            theta += if theta > 0.0 { -std::f64::consts::PI } else { std::f64::consts::PI };
            ct = -ct;
            st = -st;
        }

        // 3. Peak per vector family over core cells (un-averaged) & quadrature of primary
        let mut fam_peaks = HashMap::new();
        let mut primary_peak = 0.0f64;
        for fam in &fams {
            let mut pk2 = 0.0f64;
            for c in (0..nc).filter(|&c| core_cell[c]) {
                let base = c * nloc;
                for a in 0..nloc {
                    let idx = base + a;
                    let mut m2 = 0.0f64;
                    for &fi in &fam.fields {
                        let cf = &sol.cell_fields[fi];
                        let re = cf.re[idx];
                        let im = if cf.im.is_empty() { 0.0 } else { cf.im[idx] };
                        m2 += re * re + im * im;
                    }
                    if m2 > pk2 {
                        pk2 = m2;
                    }
                }
            }
            let pk = pk2.sqrt();
            if fam.prefix == primary.prefix {
                primary_peak = pk;
            }
            for &fi in &fam.fields {
                fam_peaks.insert(fi, pk);
            }
        }

        let quadrature = along.filter(|_| complex).and_then(|ci| {
            let cf = &sol.cell_fields[primary.fields[ci]];
            let (mut er, mut ei) = (0.0f64, 0.0f64);
            for c in (0..nc).filter(|&c| core_cell[c]) {
                let base = c * nloc;
                for a in 0..nloc {
                    let idx = base + a;
                    let re = cf.re[idx];
                    let im = if cf.im.is_empty() { 0.0 } else { cf.im[idx] };
                    let r_rot = re * ct - im * st;
                    let i_rot = re * st + im * ct;
                    er += r_rot * r_rot;
                    ei += i_rot * i_rot;
                }
            }
            (er + ei > 0.0).then_some(er / (er + ei))
        });

        (primary.prefix.clone(), theta, flipped, primary_peak, quadrature, complex, fam_peaks)
    } else if !sol.cell_fields.is_empty() {
        let cf = &sol.cell_fields[0];
        let complex = !cf.im.is_empty();
        let mut pk2 = 0.0f64;
        for c in (0..nc).filter(|&c| core_cell[c]) {
            for a in 0..nloc {
                let idx = c * nloc + a;
                let re = cf.re[idx];
                let im = if complex { cf.im[idx] } else { 0.0 };
                pk2 = pk2.max(re * re + im * im);
            }
        }
        let pk = pk2.sqrt();
        let mut fp = HashMap::new();
        for i in 0..sol.cell_fields.len() {
            fp.insert(i, pk);
        }
        (cf.name.clone(), 0.0, false, pk, None, complex, fp)
    } else if !fallback_re.is_empty() {
        let pk = fallback_re[0].iter().fold(0.0f64, |m, x| m.max(x.abs()));
        (fallback_names[0].clone(), 0.0, false, pk, None, false, HashMap::new())
    } else {
        return Err(Error::invalid(format!("the solution ('{}') has no fields to draw", sol.label)));
    };

    // Build expression Tree
    let empty_slice: &[f64] = &[];
    let (tree, expr_label, field_slots) = if !sol.cell_fields.is_empty() {
        let mut binder = FormulaBinder::new(sol, &fams);
        let (tree, label) = match &o.expr {
            Some(raw) => {
                let norm_expr = rewrite_vector_abs(raw, sol, &fams);
                let t = formula::parse(&norm_expr, &mut binder)
                    .map_err(|e| Error::invalid(format!("the formula '{raw}' does not read: {e}")))?;
                (t, raw.clone())
            }
            None => {
                if let Some(primary) = fams.iter().find(|f| f.prefix == field_name) {
                    let mut terms =
                        primary.fields.iter().map(|&fi| Tree::Call(Fun::Abs2, Box::new(binder.leaf_for_cf(fi))));
                    let first = terms.next().expect("at least one component");
                    let sum = terms.fold(first, |acc, t| Tree::Bin(Op::Add, Box::new(acc), Box::new(t)));
                    (Tree::Call(Fun::RealSqrt, Box::new(sum)), format!("|{field_name}|"))
                } else {
                    let leaf = binder.leaf_for_cf(0);
                    (leaf, field_name.clone())
                }
            }
        };
        let bloch_complex = repeat > 1 && mode_periodic_phase != 0.0;
        let mut leaf_complex = vec![false, false, false];
        let mut slots = Vec::with_capacity(binder.slots.len());
        let (ct, st) = if complex {
            (theta.cos(), theta.sin())
        } else if flipped {
            (-1.0, 0.0)
        } else {
            (1.0, 0.0)
        };
        for &fi in &binder.slots {
            let cf = &sol.cell_fields[fi];
            let pk = fam_peaks.get(&fi).copied().unwrap_or(1.0);
            let inv_pk = if pk > 0.0 { 1.0 / pk } else { 1.0 };
            let mult = C::new(ct * inv_pk, st * inv_pk);
            let is_c = !cf.im.is_empty() || bloch_complex;
            leaf_complex.push(is_c);
            let im_ref: &[f64] = if cf.im.is_empty() { empty_slice } else { &cf.im };
            slots.push((cf.re.as_slice(), im_ref, mult));
        }
        tree.check(&leaf_complex).map_err(|m| Error::invalid(format!("the formula '{label}': {m}")))?;
        if !tree.is_real_valued(&leaf_complex) {
            return Err(Error::invalid(format!(
                "the formula '{label}' is complex-valued; wrap it in re(...), im(...), abs(...) or abs2(...)"
            )));
        }
        (tree, label, slots)
    } else {
        // Fallback scalar nodal field
        let label = o.expr.clone().unwrap_or_else(|| fallback_names[0].clone());
        let idx = fallback_names.iter().position(|n| n == &label).unwrap_or(0);
        if idx != 0 {
            fallback_re.swap(0, idx);
        }
        let inv_pk = if peak > 0.0 { 1.0 / peak } else { 1.0 };
        let slots = vec![(empty_slice, empty_slice, C::new(inv_pk, 0.0))];
        (Tree::Leaf(3), label, slots)
    };

    let mut em = EvalMode {
        note: ModeNote {
            mode: k + 1,
            field: field_name,
            phase_deg: theta.to_degrees(),
            flipped,
            peak,
            quadrature,
            complex,
        },
        nloc,
        order,
        tree,
        expr_label,
        field_slots,
        fallback_re,
        periodic_phase: mode_periodic_phase,
        vmin: 0.0,
        vmax: 1.0,
        signed: false,
        abs_max: 1.0,
    };

    // Parallel range scan across core cells & periodic copies (zero intermediate field allocation)
    let bloch_rots: Vec<C> = (0..repeat)
        .map(|m| {
            let ph = -(m as f64) * mode_periodic_phase;
            C::new(ph.cos(), ph.sin())
        })
        .collect();
    let cell_ids: Vec<usize> = (0..nc).collect();
    let ranges = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
        let mut lo_v = f64::INFINITY;
        let mut hi_v = f64::NEG_INFINITY;
        let mut vals = [0.0f64; 10];
        for &c in chunk {
            if !core_cell[c] {
                continue;
            }
            for (m, &brot) in bloch_rots.iter().enumerate() {
                let sh = scale3(shift, m as f64);
                em.eval_cell(c, &topo.geom10[c], sh, brot, &mut vals);
                for &v in &vals[..nloc] {
                    if v.is_finite() {
                        lo_v = lo_v.min(v);
                        hi_v = hi_v.max(v);
                    }
                }
            }
        }
        (lo_v, hi_v)
    });
    let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for (a, b) in ranges {
        vmin = vmin.min(a);
        vmax = vmax.max(b);
    }
    if !(vmin <= vmax && vmin.is_finite() && vmax.is_finite()) {
        vmin = 0.0;
        vmax = 1.0;
    }
    let abs_max = vmin.abs().max(vmax.abs()).max(1e-30);
    // Signed Coolwarm symmetric around 0 when expr evaluates both negative and positive values (ignore sub-1e-12 rounding dust on abs/abs2)
    let signed = vmin < -1e-12 * abs_max && vmax > 1e-12 * abs_max;
    em.vmin = vmin;
    em.vmax = vmax;
    em.signed = signed;
    em.abs_max = abs_max;
    Ok(em)
}

impl EvalMode<'_> {
    #[inline]
    fn eval_cell(&self, cell: usize, geom: &[[f64; 3]; 10], shift: [f64; 3], bloch_rot: C, out: &mut [f64; 10]) {
        let nloc = self.nloc;
        let base = cell * nloc;
        let mut stack_leaf = [C::new(0.0, 0.0); 16];
        let total_leaves = 3 + self.field_slots.len();
        if total_leaves <= 16 {
            for a in 0..nloc {
                let p = add3(geom[a], shift);
                stack_leaf[0] = C::new(p[0], 0.0);
                stack_leaf[1] = C::new(p[1], 0.0);
                stack_leaf[2] = C::new(p[2], 0.0);
                let idx = base + a;
                if self.fallback_re.is_empty() {
                    for (s, &(re_s, im_s, mult)) in self.field_slots.iter().enumerate() {
                        let r0 = re_s[idx];
                        let i0 = if im_s.is_empty() { 0.0 } else { im_s[idx] };
                        let z0 = C::new(r0 * mult.re - i0 * mult.im, r0 * mult.im + i0 * mult.re);
                        stack_leaf[3 + s] = C::new(
                            z0.re * bloch_rot.re - z0.im * bloch_rot.im,
                            z0.re * bloch_rot.im + z0.im * bloch_rot.re,
                        );
                    }
                } else {
                    let r0 = self.fallback_re[0][idx] * self.field_slots[0].2.re;
                    stack_leaf[3] = C::new(r0 * bloch_rot.re, r0 * bloch_rot.im);
                }
                out[a] = self.tree.eval(&stack_leaf[..total_leaves]).re;
            }
        } else {
            let mut heap_leaf = vec![C::new(0.0, 0.0); total_leaves];
            for a in 0..nloc {
                let p = add3(geom[a], shift);
                heap_leaf[0] = C::new(p[0], 0.0);
                heap_leaf[1] = C::new(p[1], 0.0);
                heap_leaf[2] = C::new(p[2], 0.0);
                let idx = base + a;
                for (s, &(re_s, im_s, mult)) in self.field_slots.iter().enumerate() {
                    let r0 = re_s[idx];
                    let i0 = if im_s.is_empty() { 0.0 } else { im_s[idx] };
                    let z0 = C::new(r0 * mult.re - i0 * mult.im, r0 * mult.im + i0 * mult.re);
                    heap_leaf[3 + s] = C::new(
                        z0.re * bloch_rot.re - z0.im * bloch_rot.im,
                        z0.re * bloch_rot.im + z0.im * bloch_rot.re,
                    );
                }
                out[a] = self.tree.eval(&heap_leaf).re;
            }
        }
    }

    #[inline]
    fn colour(&self, v: f64, base_cmap: Cmap) -> [f32; 3] {
        if self.signed {
            Cmap::Coolwarm.at(0.5 + 0.5 * (v / self.abs_max))
        } else if self.vmax > 0.0 {
            base_cmap.at(v / self.vmax)
        } else if self.vmax > self.vmin {
            base_cmap.at((v - self.vmin) / (self.vmax - self.vmin))
        } else {
            base_cmap.at(0.0)
        }
    }
}

// ----------------------------------------------------------------------------------------- sub-tet slicing & clipping primitives

#[derive(Clone, Copy)]
struct PV {
    p: [f64; 3],
    v: f64,
}

#[inline]
fn interp_pv(a: PV, b: PV, da: f64, db: f64) -> PV {
    let denom = da - db;
    let t = if denom.abs() > 1e-300 { (da / denom).clamp(0.0, 1.0) } else { 0.5 };
    PV { p: lerp3(a.p, b.p, t), v: a.v + (b.v - a.v) * t }
}

/// Slices a single linear tetrahedron `tet` at `d == 0` (with Simulation-of-Simplicity `eps > 0` so `d == 0` never stalls),
/// appending 0, 1, or 2 triangles oriented with `plane_normal` into `out`.
fn slice_subtet(tet: [PV; 4], mut d: [f64; 4], eps: f64, plane_normal: [f64; 3], out: &mut Vec<[PV; 3]>) {
    for di in &mut d {
        if di.abs() <= eps {
            *di = eps;
        }
    }
    let mask =
        ((d[0] > 0.0) as u8) | (((d[1] > 0.0) as u8) << 1) | (((d[2] > 0.0) as u8) << 2) | (((d[3] > 0.0) as u8) << 3);
    if mask == 0 || mask == 15 {
        return;
    }
    let mut pts = [PV { p: [0.0; 3], v: 0.0 }; 4];
    let mut n = 0usize;
    for &(i, j) in &SUBTET_EDGES {
        if (d[i] > 0.0) != (d[j] > 0.0) {
            pts[n] = interp_pv(tet[i], tet[j], d[i], d[j]);
            n += 1;
        }
    }
    let orient = |mut tri: [PV; 3]| -> [PV; 3] {
        let nrm = cross3(sub3(tri[1].p, tri[0].p), sub3(tri[2].p, tri[0].p));
        if dot3(nrm, plane_normal) < 0.0 {
            tri.swap(1, 2);
        }
        tri
    };
    if n == 3 {
        out.push(orient([pts[0], pts[1], pts[2]]));
    } else if n == 4 {
        // Order the 4 quad vertices cyclically around their centroid in the plane orthogonal to `plane_normal`
        let cen = scale3(add3(add3(pts[0].p, pts[1].p), add3(pts[2].p, pts[3].p)), 0.25);
        let u_ax = normalize3(sub3(pts[0].p, cen));
        let v_ax = cross3(plane_normal, u_ax);
        let mut angs: [(f64, PV); 4] = std::array::from_fn(|i| {
            let r = sub3(pts[i].p, cen);
            (dot3(r, v_ax).atan2(dot3(r, u_ax)), pts[i])
        });
        angs.sort_by(|a, b| a.0.total_cmp(&b.0));
        out.push(orient([angs[0].1, angs[1].1, angs[2].1]));
        out.push(orient([angs[0].1, angs[2].1, angs[3].1]));
    }
}

/// Active half-space `{ p : sign * p[axis] >= sign * bound }` (`sign = +1.0` for `min`, `-1.0` for `max`).
#[derive(Clone, Copy)]
struct HalfSpace {
    axis: usize,
    sign: f64,
    bound: f64,
    /// Outward normal of the exposed cut face (`-sign` along `axis`).
    out_normal: [f64; 3],
}

impl HalfSpace {
    #[inline]
    fn signed_inside(&self, p: [f64; 3]) -> f64 {
        self.sign * (p[self.axis] - self.bound)
    }
}

fn active_half_spaces(cb: &ClipBox, tot_lo: [f64; 3], tot_hi: [f64; 3]) -> (Vec<HalfSpace>, bool) {
    let span = (0..3).map(|a| tot_hi[a] - tot_lo[a]).fold(0.0f64, f64::max).max(1e-30);
    let tol = 1e-9 * span;
    let nudge = 1e-11 * span;
    let mut planes = Vec::new();
    let pairs = [(0, cb.xmin, cb.xmax), (1, cb.ymin, cb.ymax), (2, cb.zmin, cb.zmax)];
    for (axis, lo_opt, hi_opt) in pairs {
        if let (Some(a), Some(b)) = (lo_opt, hi_opt) {
            if a >= b {
                return (Vec::new(), true);
            }
        }
        if let Some(lo_v) = lo_opt {
            if lo_v >= tot_hi[axis] - tol {
                return (Vec::new(), true);
            }
            if lo_v > tot_lo[axis] + tol {
                let mut n = [0.0; 3];
                n[axis] = -1.0;
                planes.push(HalfSpace { axis, sign: 1.0, bound: lo_v + nudge, out_normal: n });
            }
        }
        if let Some(hi_v) = hi_opt {
            if hi_v <= tot_lo[axis] + tol {
                return (Vec::new(), true);
            }
            if hi_v < tot_hi[axis] - tol {
                let mut n = [0.0; 3];
                n[axis] = 1.0;
                planes.push(HalfSpace { axis, sign: -1.0, bound: hi_v - nudge, out_normal: n });
            }
        }
    }
    (planes, false)
}

/// Sutherland-Hodgman clip of a 3D triangle against `planes` (skipping `skip_idx` when clipping a cut face on plane `skip_idx`).
fn clip_poly_against_planes(tri: [PV; 3], planes: &[HalfSpace], skip_idx: Option<usize>, out: &mut Vec<[PV; 3]>) {
    let mut poly: Vec<PV> = vec![tri[0], tri[1], tri[2]];
    let mut next: Vec<PV> = Vec::with_capacity(8);
    for (pi, pl) in planes.iter().enumerate() {
        if Some(pi) == skip_idx {
            continue;
        }
        next.clear();
        let m = poly.len();
        if m < 3 {
            return;
        }
        for i in 0..m {
            let cur = poly[i];
            let prv = poly[(i + m - 1) % m];
            let dc = pl.signed_inside(cur.p);
            let dp = pl.signed_inside(prv.p);
            let cur_in = dc >= 0.0;
            let prv_in = dp >= 0.0;
            if cur_in != prv_in {
                next.push(interp_pv(prv, cur, dp, dc));
            }
            if cur_in {
                next.push(cur);
            }
        }
        std::mem::swap(&mut poly, &mut next);
    }
    if poly.len() >= 3 {
        for i in 1..poly.len() - 1 {
            out.push([poly[0], poly[i], poly[i + 1]]);
        }
    }
}

/// Clips a 3D line segment `[a, b]` against `planes`.
fn clip_seg3d(mut a: [f64; 3], mut b: [f64; 3], planes: &[HalfSpace]) -> Option<([f64; 3], [f64; 3])> {
    for pl in planes {
        let da = pl.signed_inside(a);
        let db = pl.signed_inside(b);
        if da < 0.0 && db < 0.0 {
            return None;
        }
        if (da < 0.0) != (db < 0.0) {
            let t = (da / (da - db)).clamp(0.0, 1.0);
            let p = lerp3(a, b, t);
            if da < 0.0 {
                a = p;
            } else {
                b = p;
            }
        }
    }
    Some((a, b))
}

// ----------------------------------------------------------------------------------------------- plane frame for Slice3d

struct SliceFrame {
    origin: [f64; 3],
    normal: [f64; 3],
    e_u: [f64; 3],
    e_v: [f64; 3],
    at: Option<f64>,
    subtitle: String,
}

impl SliceFrame {
    fn from_options(o: &FigureOptions, tot_lo: [f64; 3], tot_hi: [f64; 3]) -> Result<SliceFrame, Error> {
        let spec =
            o.plane.unwrap_or_else(|| PlaneSpec::Axis { axis: 'z', at: o.at.unwrap_or(0.5 * (tot_lo[2] + tot_hi[2])) });
        match spec {
            PlaneSpec::Axis { axis, at } => {
                if !at.is_finite() {
                    return Err(Error::invalid("plane `at` coordinate must be finite"));
                }
                match axis.to_ascii_lowercase() {
                    'x' => Ok(SliceFrame {
                        origin: [at, 0.0, 0.0],
                        normal: [1.0, 0.0, 0.0],
                        e_u: [0.0, 1.0, 0.0],
                        e_v: [0.0, 0.0, 1.0],
                        at: Some(at),
                        subtitle: format!("slice x = {} (axes y, z)", si(at, "m")),
                    }),
                    'y' => Ok(SliceFrame {
                        origin: [0.0, at, 0.0],
                        normal: [0.0, 1.0, 0.0],
                        e_u: [1.0, 0.0, 0.0],
                        e_v: [0.0, 0.0, 1.0],
                        at: Some(at),
                        subtitle: format!("slice y = {} (axes x, z)", si(at, "m")),
                    }),
                    'z' => Ok(SliceFrame {
                        origin: [0.0, 0.0, at],
                        normal: [0.0, 0.0, 1.0],
                        e_u: [1.0, 0.0, 0.0],
                        e_v: [0.0, 1.0, 0.0],
                        at: Some(at),
                        subtitle: format!("slice z = {} (axes x, y)", si(at, "m")),
                    }),
                    other => Err(Error::invalid(format!("plane axis must be 'x', 'y' or 'z', not '{other}'"))),
                }
            }
            PlaneSpec::General { origin, normal } => {
                if origin.iter().chain(normal.iter()).any(|x| !x.is_finite()) {
                    return Err(Error::invalid("plane `origin` and `normal` must be finite"));
                }
                let nlen = norm3(normal);
                if nlen <= 1e-30 {
                    return Err(Error::invalid("plane `normal` must be non-zero"));
                }
                let n = scale3(normal, 1.0 / nlen);
                let up_ref = if n[2].abs() < 0.99 { [0.0, 0.0, 1.0] } else { [0.0, 1.0, 0.0] };
                let e_v = normalize3(sub3(up_ref, scale3(n, dot3(up_ref, n))));
                let e_u = normalize3(cross3(e_v, n));
                Ok(SliceFrame {
                    origin,
                    normal: n,
                    e_u,
                    e_v,
                    at: None,
                    subtitle: format!("slice n = ({:.2}, {:.2}, {:.2})", n[0], n[1], n[2]),
                })
            }
        }
    }

    #[inline]
    fn dist(&self, p: [f64; 3]) -> f64 {
        dot3(sub3(p, self.origin), self.normal)
    }

    #[inline]
    fn project_2d(&self, p: [f64; 3]) -> [f64; 2] {
        let r = sub3(p, self.origin);
        [dot3(r, self.e_u), dot3(r, self.e_v)]
    }
}

// ----------------------------------------------------------------------------------------- main 3D figure entry point

fn mode_bloch_phase(
    result: &SolveResult,
    k: usize,
    lambda: f64,
    fallback_phase: f64,
    base_phase_1: Option<f64>,
) -> f64 {
    let exact_scalar = |name: &str| -> Option<f64> {
        result
            .scalars
            .iter()
            .find(|(n, _)| n == name)
            .and_then(|(_, v)| if let Scalar::Num(x) = v { Some(*x) } else { None })
    };
    let mode_phase_of = |m_1based: usize| -> Option<f64> {
        exact_scalar(&format!("phase {m_1based} (rad)"))
            .or_else(|| exact_scalar(&format!("k_y re {m_1based} (rad/m)")).map(|ky| ky * lambda))
    };
    // If fallback_phase was explicitly set by caller to something different from mode 1's scalar phase, honour caller's override
    if fallback_phase != 0.0 && base_phase_1.is_none_or(|p1| (fallback_phase - p1).abs() > 1e-12) {
        return fallback_phase;
    }
    mode_phase_of(k + 1).unwrap_or(fallback_phase)
}

/// Renders `Slice3d`, `Cutaway3d`, `Isosurface3d`, or `Mesh3d` on a 3D `SimplexMesh`.
pub fn figure_3d_on(mesh: &SimplexMesh, result: &SolveResult, o: &FigureOptions) -> Result<Figure, Error> {
    if mesh.dim != 3 {
        return Err(Error::invalid(format!("`{}` requires a 3D mesh (this mesh is {}D)", o.kind.as_str(), mesh.dim)));
    }
    let repeat = o.repeat_periodic.max(1);
    let topo = Mesh3dTopo::build(mesh);
    let shift = if repeat > 1 { o.periodic_shift.unwrap_or([0.0, topo.hi[1] - topo.lo[1], 0.0]) } else { [0.0; 3] };
    let lambda = norm3(shift);

    // Tiled bounding box
    let mut tot_lo = topo.lo;
    let mut tot_hi = topo.hi;
    let max_shift = scale3(shift, (repeat - 1) as f64);
    for a in 0..3 {
        tot_lo[a] = tot_lo[a].min(topo.lo[a] + max_shift[a]);
        tot_hi[a] = tot_hi[a].max(topo.hi[a] + max_shift[a]);
    }

    if o.kind == Kind::Mesh3d {
        if let Some(m) = &o.modes {
            if let Some(&bad) = m.iter().find(|&&k| k >= result.solutions.len()) {
                return Err(Error::invalid(format!(
                    "mode {} is not in the result, which has {}",
                    bad + 1,
                    result.solutions.len()
                )));
            }
        }
        return render_mesh_3d(mesh, &topo, result, o, shift, repeat, tot_lo, tot_hi);
    }

    if result.solutions.is_empty() {
        return Err(Error::invalid("the result has no solutions to draw"));
    }
    let picks: Vec<usize> = match &o.modes {
        Some(m) => m.clone(),
        None => (0..result.solutions.len().min(o.max_modes.max(1))).collect(),
    };
    if picks.is_empty() {
        return Err(Error::invalid("`modes` selects no modes to draw"));
    }
    if let Some(&bad) = picks.iter().find(|&&k| k >= result.solutions.len()) {
        return Err(Error::invalid(format!(
            "mode {} is not in the result, which has {}",
            bad + 1,
            result.solutions.len()
        )));
    }

    let base_phase_1 = {
        let p1 = result
            .scalars
            .iter()
            .find(|(n, _)| n == "phase 1 (rad)")
            .and_then(|(_, v)| if let Scalar::Num(x) = v { Some(*x) } else { None })
            .or_else(|| {
                result.scalars.iter().find(|(n, _)| n == "k_y re 1 (rad/m)").and_then(|(_, v)| {
                    if let Scalar::Num(x) = v {
                        Some(*x * lambda)
                    } else {
                        None
                    }
                })
            });
        p1
    };

    let mut modes = Vec::with_capacity(picks.len());
    for &k in &picks {
        let m_phase = mode_bloch_phase(result, k, lambda, o.periodic_phase, base_phase_1);
        modes.push((k, prepare_mode(mesh, &topo, result, k, o, shift, repeat, m_phase)?));
    }

    match o.kind {
        Kind::Slice3d => render_slice_3d(mesh, &topo, result, o, &modes, shift, repeat, tot_lo, tot_hi),
        Kind::Cutaway3d => render_cutaway_3d(mesh, &topo, result, o, &modes, shift, repeat, tot_lo, tot_hi),
        Kind::Isosurface3d => render_isosurface_3d(mesh, &topo, result, o, &modes, shift, repeat, tot_lo, tot_hi),
        _ => unreachable!(),
    }
}

// ----------------------------------------------------------------------------------------- 1. Slice3d implementation

#[allow(clippy::too_many_arguments)]
fn render_slice_3d(
    mesh: &SimplexMesh,
    topo: &Mesh3dTopo,
    result: &SolveResult,
    o: &FigureOptions,
    modes: &[(usize, EvalMode<'_>)],
    shift: [f64; 3],
    repeat: usize,
    tot_lo: [f64; 3],
    tot_hi: [f64; 3],
) -> Result<Figure, Error> {
    let mut sf = SliceFrame::from_options(o, tot_lo, tot_hi)?;

    // Check distance range of tiled mesh corners against slice plane; nudge inward by 1e-11 * span if sitting right on an outer mesh face
    let mut d_min = f64::INFINITY;
    let mut d_max = f64::NEG_INFINITY;
    for ix in 0..2 {
        for iy in 0..2 {
            for iz in 0..2 {
                let p = [
                    if ix == 0 { tot_lo[0] } else { tot_hi[0] },
                    if iy == 0 { tot_lo[1] } else { tot_hi[1] },
                    if iz == 0 { tot_lo[2] } else { tot_hi[2] },
                ];
                let d = sf.dist(p);
                d_min = d_min.min(d);
                d_max = d_max.max(d);
            }
        }
    }
    let span = (0..3).map(|a| tot_hi[a] - tot_lo[a]).fold(0.0f64, f64::max).max(1e-30);
    let eps = 1e-13 * span;
    if d_min > 1e-7 * span || d_max < -1e-7 * span {
        return Err(Error::invalid("the slice plane does not intersect the 3D mesh"));
    }
    if d_min >= -1e-10 * span && d_max > 1e-9 * span {
        sf.origin = add3(sf.origin, scale3(sf.normal, 1e-11 * span));
    } else if d_max <= 1e-10 * span && d_min < -1e-9 * span {
        sf.origin = add3(sf.origin, scale3(sf.normal, -1e-11 * span));
    }

    // Extract crisp 2D region/exterior vector outline segments on the slice plane once
    let nc = mesh.n_cells();
    let cell_ids: Vec<usize> = (0..nc).collect();
    let outline_segs: Vec<[[f64; 2]; 2]> = if o.show_edges == EdgeMode::Regions {
        let seg_chunks = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
            let mut segs = Vec::new();
            for &c in chunk {
                let dom_c = mesh.cell_labels[c];
                let g0 = &topo.geom10[c];
                for m in 0..repeat {
                    let sh = scale3(shift, m as f64);
                    for (f, &(fv, _, mids)) in TET_FACES.iter().enumerate() {
                        let nb = topo.face_neighbor[c * 4 + f];
                        let draw_face = if nb == u32::MAX {
                            true
                        } else {
                            let c2 = nb as usize;
                            // Emit interior inter-region interface once (from lower cell index)
                            c < c2 && dom_c != mesh.cell_labels[c2]
                        };
                        if !draw_face {
                            continue;
                        }
                        let sub_tris: &[[usize; 3]] = if result.order == 2 {
                            &[
                                [fv[0], mids[0], mids[1]],
                                [mids[0], fv[1], mids[2]],
                                [mids[1], mids[2], fv[2]],
                                [mids[0], mids[2], mids[1]],
                            ]
                        } else {
                            &[[fv[0], fv[1], fv[2]]]
                        };
                        for st in sub_tris {
                            let p = [add3(g0[st[0]], sh), add3(g0[st[1]], sh), add3(g0[st[2]], sh)];
                            let mut d = [sf.dist(p[0]), sf.dist(p[1]), sf.dist(p[2])];
                            for di in &mut d {
                                if di.abs() <= eps {
                                    *di = eps;
                                }
                            }
                            let mut cross_pts = [[0.0; 2]; 3];
                            let mut np = 0usize;
                            for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                                if (d[a] > 0.0) != (d[b] > 0.0) && np < 3 {
                                    let t = (d[a] / (d[a] - d[b])).clamp(0.0, 1.0);
                                    cross_pts[np] = sf.project_2d(lerp3(p[a], p[b], t));
                                    np += 1;
                                }
                            }
                            if np == 2 {
                                segs.push([cross_pts[0], cross_pts[1]]);
                            }
                        }
                    }
                }
            }
            segs
        });
        seg_chunks.into_iter().flatten().collect()
    } else {
        Vec::new()
    };

    // Extract sliced 2D field triangles per mode in parallel
    type ModeSlice = (Vec<ColouredTri>, [f64; 4]);
    let mut mode_slices: Vec<ModeSlice> = Vec::with_capacity(modes.len());
    for (_, em) in modes {
        let bloch_rots: Vec<C> = (0..repeat)
            .map(|m| {
                let ph = -(m as f64) * em.periodic_phase;
                C::new(ph.cos(), ph.sin())
            })
            .collect();
        let subtets: &[[usize; 4]] = if em.order == 2 { &P2_SUBTETS } else { &P1_SUBTETS };
        let tri_chunks = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
            let mut tris: Vec<ColouredTri> = Vec::new();
            let mut bb = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
            let mut vals = [0.0f64; 10];
            let mut p_sh = [[0.0; 3]; 10];
            let mut d10 = [0.0f64; 10];
            let mut cut_buf = Vec::with_capacity(4);

            for &c in chunk {
                let g0 = &topo.geom10[c];
                for (m, &brot) in bloch_rots.iter().enumerate() {
                    let sh = scale3(shift, m as f64);
                    // Fast reject on the 4 vertices (or 10 nodes if curved)
                    let n_check = if mesh.has_curved() { em.nloc } else { 4 };
                    let (mut c_dmin, mut c_dmax) = (f64::INFINITY, f64::NEG_INFINITY);
                    for a in 0..n_check {
                        let pa = add3(g0[a], sh);
                        p_sh[a] = pa;
                        let da = sf.dist(pa);
                        d10[a] = da;
                        c_dmin = c_dmin.min(da);
                        c_dmax = c_dmax.max(da);
                    }
                    if c_dmin > eps * 4.0 || c_dmax < -eps * 4.0 {
                        continue;
                    }
                    for a in n_check..em.nloc {
                        let pa = add3(g0[a], sh);
                        p_sh[a] = pa;
                        d10[a] = sf.dist(pa);
                    }
                    em.eval_cell(c, g0, sh, brot, &mut vals);
                    for st in subtets {
                        let tet = [
                            PV { p: p_sh[st[0]], v: vals[st[0]] },
                            PV { p: p_sh[st[1]], v: vals[st[1]] },
                            PV { p: p_sh[st[2]], v: vals[st[2]] },
                            PV { p: p_sh[st[3]], v: vals[st[3]] },
                        ];
                        let dt = [d10[st[0]], d10[st[1]], d10[st[2]], d10[st[3]]];
                        cut_buf.clear();
                        slice_subtet(tet, dt, eps, sf.normal, &mut cut_buf);
                        for tri in &cut_buf {
                            let q = [sf.project_2d(tri[0].p), sf.project_2d(tri[1].p), sf.project_2d(tri[2].p)];
                            for qi in &q {
                                bb[0] = bb[0].min(qi[0]);
                                bb[1] = bb[1].max(qi[0]);
                                bb[2] = bb[2].min(qi[1]);
                                bb[3] = bb[3].max(qi[1]);
                            }
                            let cols =
                                [em.colour(tri[0].v, o.cmap), em.colour(tri[1].v, o.cmap), em.colour(tri[2].v, o.cmap)];
                            tris.push((q, cols));
                        }
                    }
                }
            }
            (tris, bb)
        });
        let mut all_tris = Vec::new();
        let mut bb = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
        for (t, b) in tri_chunks {
            all_tris.extend(t);
            bb[0] = bb[0].min(b[0]);
            bb[1] = bb[1].max(b[1]);
            bb[2] = bb[2].min(b[2]);
            bb[3] = bb[3].max(b[3]);
        }
        mode_slices.push((all_tris, bb));
    }

    if mode_slices.iter().all(|(t, _)| t.is_empty()) {
        return Err(Error::invalid("the slice plane intersected no cells of the 3D mesh"));
    }

    // Window & layout
    let win = match o.window {
        Some(w) => w,
        None => {
            let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
            for (_, mb) in &mode_slices {
                b[0] = b[0].min(mb[0]);
                b[1] = b[1].max(mb[1]);
                b[2] = b[2].min(mb[2]);
                b[3] = b[3].max(mb[3]);
            }
            let pad = 0.04 * (b[1] - b[0]).max(b[3] - b[2]).max(1e-12);
            [b[0] - pad, b[1] + pad, b[2] - pad, b[3] + pad]
        }
    };
    if !(win[1] > win[0] && win[3] > win[2]) {
        return Err(Error::invalid("the window of the figure is empty"));
    }
    let aspect = ((win[3] - win[2]) / (win[1] - win[0])).clamp(0.25, 1.6);
    let (pw, ph) = (o.panel_width, o.panel_width * aspect);
    let win = fit_aspect(win, aspect);

    let n_panels = modes.len();
    let ncols = o.columns.clamp(1, 8).min(n_panels);
    let nrows = n_panels.div_ceil(ncols);
    let (left, lab, gap, head_h, top_title, foot) = (14.0, 42.0, 14.0, 46.0, 54.0, 68.0);
    let stride = lab + pw + gap;
    let (_, uname) = length_unit(win[1] - win[0]);
    let key_note = format!(
        "{} (un-averaged P2 slice, axes in {uname}{})",
        modes[0].1.expr_label,
        if repeat > 1 { format!(", tiled {repeat}x") } else { String::new() }
    );
    let width =
        (left + ncols as f64 * stride + 6.0).max(left + lab + 180.0 + 12.0 + text_width(&key_note, 11.0) + 16.0);
    let row_h = head_h + ph + 30.0;
    let height = top_title + nrows as f64 * row_h + foot;
    let mut cv = Canvas::new(width, height);
    cv.bold(left, 30.0, &format!("{}: {} ({})", result.physics, result.study, sf.subtitle), 18.0, INK, Anchor::Start);
    let (unit, _) = length_unit(win[1] - win[0]);

    for (idx, ((k, em), (tris, _))) in modes.iter().zip(mode_slices.iter()).enumerate() {
        let (r, c0) = (idx / ncols, idx % ncols);
        let x = left + lab + c0 as f64 * stride;
        let y = top_title + r as f64 * row_h + head_h;
        let fr = Frame { x, y, w: pw, h: ph, win };
        cv.push(Item::Rect { x, y, w: pw, h: ph, fill: Some(PANEL_BG), stroke: None });

        let (iw, ih) = ((pw * 2.0).ceil() as u32, (ph * 2.0).ceil() as u32);
        let to_img = |q: [f64; 2]| -> [f64; 2] {
            let f = fr.map(q);
            [(f[0] - x) * 2.0, (f[1] - y) * 2.0]
        };
        let img_tris: Vec<ColouredTri> =
            tris.iter().map(|(q, col)| ([to_img(q[0]), to_img(q[1]), to_img(q[2])], *col)).collect();
        cv.push(Item::Image { x, y, w: pw, h: ph, px: iw, py: ih, rgba: raster_field(iw, ih, &img_tris) });

        if o.show_edges == EdgeMode::Mesh {
            for (q, _) in tris {
                for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                    if let Some(pts) = clip(fr.map(q[a]), fr.map(q[b]), [x, y, x + pw, y + ph]) {
                        cv.push(Item::Line { pts: pts.to_vec(), color: [148, 163, 184, 110], width: 0.5 });
                    }
                }
            }
        } else if o.show_edges == EdgeMode::Regions {
            for seg in &outline_segs {
                if let Some(pts) = clip(fr.map(seg[0]), fr.map(seg[1]), [x, y, x + pw, y + ph]) {
                    cv.push(Item::Line { pts: pts.to_vec(), color: OUTLINE, width: 1.0 });
                }
            }
        }

        cv.push(Item::Rect { x, y, w: pw, h: ph, fill: None, stroke: Some(([71, 85, 105, 255], 0.8)) });
        draw_ticks(&mut cv, &fr, unit);
        if em.signed {
            cv.text(x + pw - 4.0, y + 12.0, &format!("+-{:.2}", em.abs_max), 10.5, WHITE_INK, Anchor::End);
        }

        let (head, sub) = mode_lines(result, *k);
        let hy = top_title + r as f64 * row_h;
        cv.bold(x, hy + 16.0, &head, 14.0, INK, Anchor::Start);
        if !sub.is_empty() {
            cv.text(x, hy + 31.0, &sub, 11.0, MUTED, Anchor::Start);
        }
    }

    draw_colorbar(
        &mut cv,
        left + lab,
        top_title + nrows as f64 * row_h + 14.0,
        if modes[0].1.signed { Cmap::Coolwarm } else { o.cmap },
        modes[0].1.signed,
        modes[0].1.abs_max,
        modes[0].1.vmax,
        &key_note,
    );

    let notes = modes.iter().map(|(_, m)| m.note.clone()).collect();
    Ok(Figure { canvas: cv, modes: notes, window: win, at: sf.at })
}

const WHITE_INK: Rgba = [241, 245, 249, 240];

// ------------------------------------------------------------------------------- 3D camera projection & shared helpers

struct Proj3d {
    center: [f64; 3],
    right: [f64; 3],
    up: [f64; 3],
    view_dir: [f64; 3],
    radius: f64,
    win: [f64; 4],
}

impl Proj3d {
    fn new(
        cam: Camera3d,
        lo: [f64; 3],
        hi: [f64; 3],
        user_win: Option<[f64; 4]>,
        aspect_clamp: (f64, f64),
    ) -> (Proj3d, f64) {
        let (right, up, view_dir) = cam.basis();
        let center = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
        let diag = norm3(sub3(hi, lo)).max(1e-30);
        let radius = 0.5 * diag;
        let mut bb = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
        for ix in 0..2 {
            for iy in 0..2 {
                for iz in 0..2 {
                    let p = [
                        if ix == 0 { lo[0] } else { hi[0] },
                        if iy == 0 { lo[1] } else { hi[1] },
                        if iz == 0 { lo[2] } else { hi[2] },
                    ];
                    let r = sub3(p, center);
                    let u = dot3(r, right);
                    let v = dot3(r, up);
                    bb[0] = bb[0].min(u);
                    bb[1] = bb[1].max(u);
                    bb[2] = bb[2].min(v);
                    bb[3] = bb[3].max(v);
                }
            }
        }
        let raw_win = user_win.unwrap_or_else(|| {
            let pad = 0.08 * (bb[1] - bb[0]).max(bb[3] - bb[2]).max(1e-12);
            [bb[0] - pad, bb[1] + pad, bb[2] - pad, bb[3] + pad]
        });
        let aspect =
            ((raw_win[3] - raw_win[2]) / (raw_win[1] - raw_win[0]).max(1e-30)).clamp(aspect_clamp.0, aspect_clamp.1);
        let win = fit_aspect(raw_win, aspect);
        (Proj3d { center, right, up, view_dir, radius, win }, aspect)
    }

    /// Maps world 3D point `p` into image pixel coordinates `(px, py)` of an `iw x ih` image plus normalized depth `z` in `[0.05, 0.95]`.
    #[inline]
    fn map_px(&self, p: [f64; 3], iw: f64, ih: f64) -> [f64; 3] {
        let r = sub3(p, self.center);
        let u = dot3(r, self.right);
        let v = dot3(r, self.up);
        let d = dot3(r, self.view_dir);
        let px = (u - self.win[0]) / (self.win[1] - self.win[0]) * iw;
        let py = (self.win[3] - v) / (self.win[3] - self.win[2]) * ih;
        let z = 0.5 + 0.45 * (d / self.radius.max(1e-30));
        [px, py, z]
    }
}

fn draw_triad(cv: &mut Canvas, cam: Camera3d, ox: f64, oy: f64, len: f64) {
    let (right, up, _) = cam.basis();
    let axes = [
        ([1.0, 0.0, 0.0], "x", [248, 113, 113, 240]),
        ([0.0, 1.0, 0.0], "y", [74, 222, 128, 240]),
        ([0.0, 0.0, 1.0], "z", [96, 165, 250, 240]),
    ];
    for (ax, name, col) in axes {
        let dx = dot3(ax, right) * len;
        let dy = -dot3(ax, up) * len;
        let end = [ox + dx, oy + dy];
        cv.push(Item::Line { pts: vec![[ox, oy], end], color: col, width: 1.5 });
        cv.text(ox + dx * 1.32, oy + dy * 1.32 + 3.5, name, 9.5, col, Anchor::Middle);
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_colorbar(cv: &mut Canvas, x0: f64, ky: f64, cmap: Cmap, signed: bool, abs_max: f64, vmax: f64, label: &str) {
    let kw = 180.0;
    let steps = 40;
    for i in 0..steps {
        let t = i as f64 / (steps - 1) as f64;
        cv.push(Item::Rect {
            x: x0 + i as f64 * kw / steps as f64,
            y: ky,
            w: kw / steps as f64 + 0.6,
            h: 10.0,
            fill: Some(bytes(cmap.at(t))),
            stroke: None,
        });
    }
    let (l_txt, r_txt) = if signed {
        (format!("-{abs_max:.2}"), format!("+{abs_max:.2}"))
    } else {
        ("0".to_string(), format!("{:.2}", if vmax > 0.0 { vmax } else { 1.0 }))
    };
    cv.text(x0, ky + 24.0, &l_txt, 10.5, MUTED, Anchor::Start);
    cv.text(x0 + kw, ky + 24.0, &r_txt, 10.5, MUTED, Anchor::End);
    cv.text(x0 + kw + 12.0, ky + 9.0, label, 11.0, MUTED, Anchor::Start);
}

// ------------------------------------------------------------------------------- marching-tets smooth isosurface extraction

type QuantEdge = ([i64; 3], [i64; 3], [f64; 3], [f64; 3]);

struct IsoTriRaw {
    p: [[f64; 3]; 3],
    shell: u16,
    dom: i32,
}

#[allow(clippy::too_many_arguments)]
fn extract_isosurface_tris(
    mesh: &SimplexMesh,
    topo: &Mesh3dTopo,
    em: &EvalMode<'_>,
    levels: &[f64],
    planes: &[HalfSpace],
    shift: [f64; 3],
    repeat: usize,
    base_cmap: Cmap,
    proj: &Proj3d,
    iw: f64,
    ih: f64,
    alpha: f32,
) -> Vec<Tri3d> {
    if levels.is_empty() {
        return Vec::new();
    }
    let nc = mesh.n_cells();
    let span = norm3(sub3(topo.hi, topo.lo)).max(1e-30);
    let val_eps = 1e-12 * em.abs_max.max(1e-30);
    let bloch_rots: Vec<C> = (0..repeat)
        .map(|m| {
            let ph = -(m as f64) * em.periodic_phase;
            C::new(ph.cos(), ph.sin())
        })
        .collect();
    let subtets: &[[usize; 4]] = if em.order == 2 { &P2_SUBTETS } else { &P1_SUBTETS };
    let cell_ids: Vec<usize> = (0..nc).collect();

    let raw_chunks = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
        let mut out: Vec<IsoTriRaw> = Vec::new();
        let mut vals = [0.0f64; 10];
        let mut p_sh = [[0.0; 3]; 10];
        let mut slice_buf = Vec::with_capacity(4);
        let mut clip_buf = Vec::with_capacity(6);

        for &c in chunk {
            let dom = mesh.cell_labels[c];
            let g0 = &topo.geom10[c];
            for (m, &brot) in bloch_rots.iter().enumerate() {
                let sh = scale3(shift, m as f64);
                // Optional ClipBox reject
                if !planes.is_empty() {
                    let mut outside = false;
                    for pl in planes {
                        if (0..4).all(|vi| pl.signed_inside(add3(g0[vi], sh)) < 0.0) && !mesh.has_curved() {
                            outside = true;
                            break;
                        }
                    }
                    if outside {
                        continue;
                    }
                }
                em.eval_cell(c, g0, sh, brot, &mut vals);
                let (mut vmin_c, mut vmax_c) = (f64::INFINITY, f64::NEG_INFINITY);
                for &v in &vals[..em.nloc] {
                    vmin_c = vmin_c.min(v);
                    vmax_c = vmax_c.max(v);
                }
                if !levels.iter().any(|&lv| lv >= vmin_c - val_eps && lv <= vmax_c + val_eps) {
                    continue;
                }
                for a in 0..em.nloc {
                    p_sh[a] = add3(g0[a], sh);
                }
                for (si, &lv) in levels.iter().enumerate() {
                    if lv < vmin_c - val_eps || lv > vmax_c + val_eps {
                        continue;
                    }
                    for st in subtets {
                        let d = [vals[st[0]] - lv, vals[st[1]] - lv, vals[st[2]] - lv, vals[st[3]] - lv];
                        // Estimate scalar gradient direction in sub-tet so normals point consistently along -grad(v) (outward from high-|v| lobe)
                        let tet = [
                            PV { p: p_sh[st[0]], v: lv },
                            PV { p: p_sh[st[1]], v: lv },
                            PV { p: p_sh[st[2]], v: lv },
                            PV { p: p_sh[st[3]], v: lv },
                        ];
                        let sign_lv = if lv < 0.0 { 1.0 } else { -1.0 };
                        let grad_approx = add3(
                            add3(
                                scale3(sub3(tet[1].p, tet[0].p), d[1] - d[0]),
                                scale3(sub3(tet[2].p, tet[0].p), d[2] - d[0]),
                            ),
                            scale3(sub3(tet[3].p, tet[0].p), d[3] - d[0]),
                        );
                        let out_dir = normalize3(scale3(grad_approx, sign_lv));
                        slice_buf.clear();
                        slice_subtet(tet, d, val_eps, out_dir, &mut slice_buf);
                        if planes.is_empty() {
                            for tri in &slice_buf {
                                out.push(IsoTriRaw { p: [tri[0].p, tri[1].p, tri[2].p], shell: si as u16, dom });
                            }
                        } else {
                            for tri in &slice_buf {
                                clip_buf.clear();
                                clip_poly_against_planes(*tri, planes, None, &mut clip_buf);
                                for ct in &clip_buf {
                                    out.push(IsoTriRaw { p: [ct[0].p, ct[1].p, ct[2].p], shell: si as u16, dom });
                                }
                            }
                        }
                    }
                }
            }
        }
        out
    });

    let raw: Vec<IsoTriRaw> = raw_chunks.into_iter().flatten().collect();
    if raw.is_empty() {
        return Vec::new();
    }

    // Smooth vertex normals per (shell, domain_label) without welding field values across domains
    let inv_q = 1e6 / span;
    let qkey = |shell: u16, dom: i32, p: [f64; 3]| -> (u16, i32, i64, i64, i64) {
        (shell, dom, (p[0] * inv_q).round() as i64, (p[1] * inv_q).round() as i64, (p[2] * inv_q).round() as i64)
    };
    let mut n_acc: HashMap<(u16, i32, i64, i64, i64), [f64; 3]> = HashMap::with_capacity(raw.len() * 2);
    for t in &raw {
        let fnrm = cross3(sub3(t.p[1], t.p[0]), sub3(t.p[2], t.p[0]));
        for k in 0..3 {
            let e = n_acc.entry(qkey(t.shell, t.dom, t.p[k])).or_insert([0.0; 3]);
            e[0] += fnrm[0];
            e[1] += fnrm[1];
            e[2] += fnrm[2];
        }
    }

    // Key light from upper-left of camera + headlight fill
    let key_light =
        normalize3(add3(add3(scale3(proj.right, -0.35), scale3(proj.up, 0.55)), scale3(proj.view_dir, -0.75)));
    let head_light = scale3(proj.view_dir, -1.0);

    let shell_cols: Vec<[f32; 3]> = levels.iter().map(|&lv| em.colour(lv, base_cmap)).collect();
    raw.iter()
        .map(|t| {
            let base = shell_cols[t.shell as usize];
            let mut c = [[0.0f32; 3]; 3];
            let mut p2d = [[0.0f64; 3]; 3];
            let fallback_n = normalize3(cross3(sub3(t.p[1], t.p[0]), sub3(t.p[2], t.p[0])));
            for k in 0..3 {
                p2d[k] = proj.map_px(t.p[k], iw, ih);
                let vn = n_acc.get(&qkey(t.shell, t.dom, t.p[k])).copied().map(normalize3).unwrap_or(fallback_n);
                let ndl = dot3(vn, key_light).abs();
                let ndh = dot3(vn, head_light).abs();
                let spec = ndl.powi(12) * 0.18;
                let shade = (0.34 + 0.52 * ndl + 0.14 * ndh) as f32;
                c[k] = [
                    (base[0] * shade + spec as f32).clamp(0.0, 1.0),
                    (base[1] * shade + spec as f32).clamp(0.0, 1.0),
                    (base[2] * shade + spec as f32).clamp(0.0, 1.0),
                ];
            }
            Tri3d { p: p2d, c, alpha }
        })
        .collect()
}

// ------------------------------------------------------------------------------- 2. Cutaway3d implementation

#[allow(clippy::too_many_arguments)]
fn render_cutaway_3d(
    mesh: &SimplexMesh,
    topo: &Mesh3dTopo,
    result: &SolveResult,
    o: &FigureOptions,
    modes: &[(usize, EvalMode<'_>)],
    shift: [f64; 3],
    repeat: usize,
    tot_lo: [f64; 3],
    tot_hi: [f64; 3],
) -> Result<Figure, Error> {
    let default_cb = ClipBox { zmax: Some(0.5 * (tot_lo[2] + tot_hi[2])), ..Default::default() };
    let cb = o.clip.as_ref().unwrap_or(&default_cb);
    let (planes, empty_clip) = active_half_spaces(cb, tot_lo, tot_hi);
    if empty_clip {
        return Err(Error::invalid("the `clip` box excludes the entire 3D mesh"));
    }

    let (proj, aspect) = Proj3d::new(o.camera, tot_lo, tot_hi, o.window, (0.35, 1.4));
    let (pw, ph) = (o.panel_width, o.panel_width * aspect);
    let (iw, ih) = ((pw * 2.0).ceil() as u32, (ph * 2.0).ceil() as u32);
    let (iw_f, ih_f) = (f64::from(iw), f64::from(ih));

    let nc = mesh.n_cells();
    let span = norm3(sub3(tot_hi, tot_lo)).max(1e-30);
    let cut_eps = 1e-13 * span;
    let inv_q = 1e6 / span;
    let qpt = |p: [f64; 3]| -> [i64; 3] {
        [(p[0] * inv_q).round() as i64, (p[1] * inv_q).round() as i64, (p[2] * inv_q).round() as i64]
    };

    // Precompute 3D feature/region edges (clipped to `planes`) when `show_edges == EdgeMode::Regions`
    let region_edges_3d: Vec<([f64; 3], [f64; 3], Rgba, f64)> = if o.show_edges == EdgeMode::Regions {
        let mut segs: Vec<([f64; 3], [f64; 3], Rgba, f64)> = Vec::new();
        let stroke: Rgba = [203, 213, 225, 215];
        // 1. Clipped exterior sharp/region feature edges
        for fe in &topo.ext_feature_edges {
            let g0 = &topo.geom10[fe.cell as usize];
            let parts: &[(usize, usize)] = if result.order == 2 && mesh.has_curved() {
                &[(fe.la as usize, fe.lm as usize), (fe.lm as usize, fe.lb as usize)]
            } else {
                &[(fe.la as usize, fe.lb as usize)]
            };
            for m in 0..repeat {
                let sh = scale3(shift, m as f64);
                for &(ia, ib) in parts {
                    if let Some((pa, pb)) = clip_seg3d(add3(g0[ia], sh), add3(g0[ib], sh), &planes) {
                        segs.push((proj.map_px(pa, iw_f, ih_f), proj.map_px(pb, iw_f, ih_f), stroke, 1.4));
                    }
                }
            }
        }
        // 2. Cut-plane interface & rim edges (exterior + inter-region faces intersected with each active clip plane)
        if !planes.is_empty() {
            let cell_ids: Vec<usize> = (0..nc).collect();
            let cut_edge_chunks = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
                let mut local = Vec::new();
                for &c in chunk {
                    let dom_c = mesh.cell_labels[c];
                    let g0 = &topo.geom10[c];
                    for m in 0..repeat {
                        let sh = scale3(shift, m as f64);
                        for (f, &(fv, _, mids)) in TET_FACES.iter().enumerate() {
                            let nb = topo.face_neighbor[c * 4 + f];
                            let keep_face = if nb == u32::MAX {
                                true
                            } else {
                                let c2 = nb as usize;
                                c < c2 && dom_c != mesh.cell_labels[c2]
                            };
                            if !keep_face {
                                continue;
                            }
                            let sub_tris: &[[usize; 3]] = if result.order == 2 {
                                &[
                                    [fv[0], mids[0], mids[1]],
                                    [mids[0], fv[1], mids[2]],
                                    [mids[1], mids[2], fv[2]],
                                    [mids[0], mids[2], mids[1]],
                                ]
                            } else {
                                &[[fv[0], fv[1], fv[2]]]
                            };
                            for (pi, pl) in planes.iter().enumerate() {
                                for st in sub_tris {
                                    let p = [add3(g0[st[0]], sh), add3(g0[st[1]], sh), add3(g0[st[2]], sh)];
                                    let mut d =
                                        [pl.signed_inside(p[0]), pl.signed_inside(p[1]), pl.signed_inside(p[2])];
                                    for di in &mut d {
                                        if di.abs() <= cut_eps {
                                            *di = cut_eps;
                                        }
                                    }
                                    let mut cp = [[0.0; 3]; 3];
                                    let mut np = 0usize;
                                    for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                                        if (d[a] > 0.0) != (d[b] > 0.0) && np < 3 {
                                            let t = (d[a] / (d[a] - d[b])).clamp(0.0, 1.0);
                                            cp[np] = lerp3(p[a], p[b], t);
                                            np += 1;
                                        }
                                    }
                                    if np == 2 {
                                        let other_planes: Vec<HalfSpace> = planes
                                            .iter()
                                            .enumerate()
                                            .filter(|&(j, _)| j != pi)
                                            .map(|(_, &x)| x)
                                            .collect();
                                        if let Some((pa, pb)) = clip_seg3d(cp[0], cp[1], &other_planes) {
                                            local.push((
                                                proj.map_px(pa, iw_f, ih_f),
                                                proj.map_px(pb, iw_f, ih_f),
                                                stroke,
                                                1.5,
                                            ));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                local
            });
            for ch in cut_edge_chunks {
                segs.extend(ch);
            }
        }
        segs
    } else {
        Vec::new()
    };

    let n_panels = modes.len();
    let ncols = o.columns.clamp(1, 8).min(n_panels);
    let nrows = n_panels.div_ceil(ncols);
    let (left, gap, head_h, top_title, foot) = (18.0, 16.0, 46.0, 54.0, 68.0);
    let stride = pw + gap;
    let key_note = format!(
        "{} (3D cutaway{}{})",
        modes[0].1.expr_label,
        if o.isosurfaces.is_some() { " + isosurfaces" } else { "" },
        if repeat > 1 { format!(", tiled {repeat}x") } else { String::new() }
    );
    let width = (left + ncols as f64 * stride + 6.0).max(left + 180.0 + 12.0 + text_width(&key_note, 11.0) + 16.0);
    let row_h = head_h + ph + 20.0;
    let height = top_title + nrows as f64 * row_h + foot;
    let mut cv = Canvas::new(width, height);
    cv.bold(left, 30.0, &format!("{}: {} (3D cutaway)", result.physics, result.study), 18.0, INK, Anchor::Start);

    let cell_ids: Vec<usize> = (0..nc).collect();
    for (idx, (k, em)) in modes.iter().enumerate() {
        let bloch_rots: Vec<C> = (0..repeat)
            .map(|m| {
                let ph = -(m as f64) * em.periodic_phase;
                C::new(ph.cos(), ph.sin())
            })
            .collect();
        let subtets: &[[usize; 4]] = if em.order == 2 { &P2_SUBTETS } else { &P1_SUBTETS };
        let collect_mesh_edges = o.show_edges == EdgeMode::Mesh;

        let chunks = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
            let mut opaque: Vec<Tri3d> = Vec::new();
            let mut raw_edges: Vec<QuantEdge> = Vec::new();
            let mut vals = [0.0f64; 10];
            let mut p_sh = [[0.0; 3]; 10];
            let mut clip_buf = Vec::with_capacity(8);
            let mut cut_buf = Vec::with_capacity(4);

            for &c in chunk {
                let g0 = &topo.geom10[c];
                let has_ext = (0..4).any(|f| topo.face_neighbor[c * 4 + f] == u32::MAX);

                for (m, &brot) in bloch_rots.iter().enumerate() {
                    let sh = scale3(shift, m as f64);
                    // Classify cell against `planes`
                    let n_chk = if mesh.has_curved() { em.nloc } else { 4 };
                    for a in 0..n_chk {
                        p_sh[a] = add3(g0[a], sh);
                    }
                    let mut outside = false;
                    let mut straddles = false;
                    for pl in &planes {
                        let (mut any_in, mut any_out) = (false, false);
                        for &pt in p_sh.iter().take(n_chk) {
                            if pl.signed_inside(pt) >= 0.0 {
                                any_in = true;
                            } else {
                                any_out = true;
                            }
                        }
                        if !any_in {
                            outside = true;
                            break;
                        }
                        if any_out {
                            straddles = true;
                        }
                    }
                    if outside || (!has_ext && !straddles) {
                        continue;
                    }
                    for a in n_chk..em.nloc {
                        p_sh[a] = add3(g0[a], sh);
                    }
                    em.eval_cell(c, g0, sh, brot, &mut vals);

                    let mut emit_tri = |tri: [PV; 3]| {
                        let p2 = [
                            proj.map_px(tri[0].p, iw_f, ih_f),
                            proj.map_px(tri[1].p, iw_f, ih_f),
                            proj.map_px(tri[2].p, iw_f, ih_f),
                        ];
                        let cols =
                            [em.colour(tri[0].v, o.cmap), em.colour(tri[1].v, o.cmap), em.colour(tri[2].v, o.cmap)];
                        opaque.push(Tri3d { p: p2, c: cols, alpha: 1.0 });
                        if collect_mesh_edges {
                            for (ia, ib) in [(0, 1), (1, 2), (2, 0)] {
                                let (qa, qb) = (qpt(tri[ia].p), qpt(tri[ib].p));
                                if qa <= qb {
                                    raw_edges.push((qa, qb, p2[ia], p2[ib]));
                                } else {
                                    raw_edges.push((qb, qa, p2[ib], p2[ia]));
                                }
                            }
                        }
                    };

                    // 1. Exterior mesh faces (unclipped if fully inside, else Sutherland-Hodgman clipped)
                    if has_ext {
                        for (f, &(fv, _, mids)) in TET_FACES.iter().enumerate() {
                            if topo.face_neighbor[c * 4 + f] != u32::MAX {
                                continue;
                            }
                            let sub_tris: &[[usize; 3]] = if em.order == 2 {
                                &[
                                    [fv[0], mids[0], mids[1]],
                                    [mids[0], fv[1], mids[2]],
                                    [mids[1], mids[2], fv[2]],
                                    [mids[0], mids[2], mids[1]],
                                ]
                            } else {
                                &[[fv[0], fv[1], fv[2]]]
                            };
                            for st in sub_tris {
                                let base_tri = [
                                    PV { p: p_sh[st[0]], v: vals[st[0]] },
                                    PV { p: p_sh[st[1]], v: vals[st[1]] },
                                    PV { p: p_sh[st[2]], v: vals[st[2]] },
                                ];
                                if !straddles {
                                    emit_tri(base_tri);
                                } else {
                                    clip_buf.clear();
                                    clip_poly_against_planes(base_tri, &planes, None, &mut clip_buf);
                                    for ct in &clip_buf {
                                        emit_tri(*ct);
                                    }
                                }
                            }
                        }
                    }

                    // 2. Exposed interior cut faces on each straddled ClipBox plane
                    if straddles {
                        for (pi, pl) in planes.iter().enumerate() {
                            for st in subtets {
                                let tet = [
                                    PV { p: p_sh[st[0]], v: vals[st[0]] },
                                    PV { p: p_sh[st[1]], v: vals[st[1]] },
                                    PV { p: p_sh[st[2]], v: vals[st[2]] },
                                    PV { p: p_sh[st[3]], v: vals[st[3]] },
                                ];
                                let dt = [
                                    pl.signed_inside(tet[0].p),
                                    pl.signed_inside(tet[1].p),
                                    pl.signed_inside(tet[2].p),
                                    pl.signed_inside(tet[3].p),
                                ];
                                cut_buf.clear();
                                slice_subtet(tet, dt, cut_eps, pl.out_normal, &mut cut_buf);
                                for ct_tri in &cut_buf {
                                    if planes.len() == 1 {
                                        emit_tri(*ct_tri);
                                    } else {
                                        clip_buf.clear();
                                        clip_poly_against_planes(*ct_tri, &planes, Some(pi), &mut clip_buf);
                                        for clipped in &clip_buf {
                                            emit_tri(*clipped);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            (opaque, raw_edges)
        });

        let mut opaque = Vec::new();
        let mut mesh_edge_raw = Vec::new();
        for (op, ed) in chunks {
            opaque.extend(op);
            mesh_edge_raw.extend(ed);
        }

        let edges_3d: Vec<([f64; 3], [f64; 3], Rgba, f64)> = match o.show_edges {
            EdgeMode::None => Vec::new(),
            EdgeMode::Regions => region_edges_3d.clone(),
            EdgeMode::Mesh => {
                mesh_edge_raw.sort_unstable_by_key(|e| (e.0, e.1));
                mesh_edge_raw.dedup_by_key(|e| (e.0, e.1));
                let stroke: Rgba = [148, 163, 184, 150];
                mesh_edge_raw.into_iter().map(|(_, _, a, b)| (a, b, stroke, 1.0)).collect()
            }
        };

        let translucent = if let Some(levels) = &o.isosurfaces {
            extract_isosurface_tris(mesh, topo, em, levels, &planes, shift, repeat, o.cmap, &proj, iw_f, ih_f, 0.52)
        } else {
            Vec::new()
        };

        let rgba = raster_field_3d(iw, ih, &opaque, &translucent, &edges_3d);
        let (r, c0) = (idx / ncols, idx % ncols);
        let x = left + c0 as f64 * stride;
        let y = top_title + r as f64 * row_h + head_h;
        cv.push(Item::Rect { x, y, w: pw, h: ph, fill: Some(PANEL_BG), stroke: None });
        cv.push(Item::Image { x, y, w: pw, h: ph, px: iw, py: ih, rgba });
        draw_triad(&mut cv, o.camera, x + 26.0, y + ph - 24.0, 15.0);
        cv.push(Item::Rect { x, y, w: pw, h: ph, fill: None, stroke: Some(([71, 85, 105, 255], 0.8)) });
        if em.signed {
            cv.text(x + pw - 6.0, y + 14.0, &format!("+-{:.2}", em.abs_max), 10.5, WHITE_INK, Anchor::End);
        }

        let (head, sub) = mode_lines(result, *k);
        let hy = top_title + r as f64 * row_h;
        cv.bold(x, hy + 16.0, &head, 14.0, INK, Anchor::Start);
        if !sub.is_empty() {
            cv.text(x, hy + 31.0, &sub, 11.0, MUTED, Anchor::Start);
        }
    }

    draw_colorbar(
        &mut cv,
        left,
        top_title + nrows as f64 * row_h + 10.0,
        if modes[0].1.signed { Cmap::Coolwarm } else { o.cmap },
        modes[0].1.signed,
        modes[0].1.abs_max,
        modes[0].1.vmax,
        &key_note,
    );

    let notes = modes.iter().map(|(_, m)| m.note.clone()).collect();
    Ok(Figure { canvas: cv, modes: notes, window: proj.win, at: None })
}

// ------------------------------------------------------------------------------- 3. Isosurface3d implementation

#[allow(clippy::too_many_arguments)]
fn render_isosurface_3d(
    mesh: &SimplexMesh,
    topo: &Mesh3dTopo,
    result: &SolveResult,
    o: &FigureOptions,
    modes: &[(usize, EvalMode<'_>)],
    shift: [f64; 3],
    repeat: usize,
    tot_lo: [f64; 3],
    tot_hi: [f64; 3],
) -> Result<Figure, Error> {
    let (planes, empty_clip) = match &o.clip {
        Some(cb) => active_half_spaces(cb, tot_lo, tot_hi),
        None => (Vec::new(), false),
    };
    if empty_clip {
        return Err(Error::invalid("the `clip` box excludes the entire 3D mesh"));
    }

    let (proj, aspect) = Proj3d::new(o.camera, tot_lo, tot_hi, o.window, (0.35, 1.4));
    let (pw, ph) = (o.panel_width, o.panel_width * aspect);
    let (iw, ih) = ((pw * 2.0).ceil() as u32, (ph * 2.0).ceil() as u32);
    let (iw_f, ih_f) = (f64::from(iw), f64::from(ih));

    let wire_edges: Vec<([f64; 3], [f64; 3], Rgba, f64)> = match o.show_edges {
        EdgeMode::None => Vec::new(),
        EdgeMode::Regions | EdgeMode::Mesh => {
            let mut segs = Vec::new();
            let stroke: Rgba = [148, 163, 184, 130];
            for fe in &topo.ext_feature_edges {
                let g0 = &topo.geom10[fe.cell as usize];
                for m in 0..repeat {
                    let sh = scale3(shift, m as f64);
                    if let Some((pa, pb)) =
                        clip_seg3d(add3(g0[fe.la as usize], sh), add3(g0[fe.lb as usize], sh), &planes)
                    {
                        segs.push((proj.map_px(pa, iw_f, ih_f), proj.map_px(pb, iw_f, ih_f), stroke, 1.1));
                    }
                }
            }
            segs
        }
    };

    let n_panels = modes.len();
    let ncols = o.columns.clamp(1, 8).min(n_panels);
    let nrows = n_panels.div_ceil(ncols);
    let (left, gap, head_h, top_title, foot) = (18.0, 16.0, 46.0, 54.0, 68.0);
    let stride = pw + gap;

    let first_levels = o.isosurfaces.clone().unwrap_or_else(|| {
        if modes[0].1.signed {
            vec![-0.5 * modes[0].1.abs_max, 0.5 * modes[0].1.abs_max]
        } else {
            let sc = if modes[0].1.vmax > 0.0 { modes[0].1.vmax } else { modes[0].1.abs_max };
            vec![0.35 * sc, 0.65 * sc]
        }
    });
    let lvl_txt: Vec<String> = first_levels.iter().map(|v| format!("{v:.2}")).collect();
    let key_note = format!(
        "{} isosurfaces at [{}]{}",
        modes[0].1.expr_label,
        lvl_txt.join(", "),
        if repeat > 1 { format!(" (tiled {repeat}x)") } else { String::new() }
    );
    let width = (left + ncols as f64 * stride + 6.0).max(left + 180.0 + 12.0 + text_width(&key_note, 11.0) + 16.0);
    let row_h = head_h + ph + 20.0;
    let height = top_title + nrows as f64 * row_h + foot;
    let mut cv = Canvas::new(width, height);
    cv.bold(left, 30.0, &format!("{}: {} (3D P2 isosurfaces)", result.physics, result.study), 18.0, INK, Anchor::Start);

    for (idx, (k, em)) in modes.iter().enumerate() {
        let levels = o.isosurfaces.clone().unwrap_or_else(|| {
            if em.signed {
                vec![-0.5 * em.abs_max, 0.5 * em.abs_max]
            } else {
                let sc = if em.vmax > 0.0 { em.vmax } else { em.abs_max };
                vec![0.35 * sc, 0.65 * sc]
            }
        });
        let translucent =
            extract_isosurface_tris(mesh, topo, em, &levels, &planes, shift, repeat, o.cmap, &proj, iw_f, ih_f, 0.58);
        let rgba = raster_field_3d(iw, ih, &[], &translucent, &wire_edges);

        let (r, c0) = (idx / ncols, idx % ncols);
        let x = left + c0 as f64 * stride;
        let y = top_title + r as f64 * row_h + head_h;
        cv.push(Item::Rect { x, y, w: pw, h: ph, fill: Some(PANEL_BG), stroke: None });
        cv.push(Item::Image { x, y, w: pw, h: ph, px: iw, py: ih, rgba });
        draw_triad(&mut cv, o.camera, x + 26.0, y + ph - 24.0, 15.0);
        cv.push(Item::Rect { x, y, w: pw, h: ph, fill: None, stroke: Some(([71, 85, 105, 255], 0.8)) });

        let (head, sub) = mode_lines(result, *k);
        let hy = top_title + r as f64 * row_h;
        cv.bold(x, hy + 16.0, &head, 14.0, INK, Anchor::Start);
        if !sub.is_empty() {
            cv.text(x, hy + 31.0, &sub, 11.0, MUTED, Anchor::Start);
        }
    }

    draw_colorbar(
        &mut cv,
        left,
        top_title + nrows as f64 * row_h + 10.0,
        if modes[0].1.signed { Cmap::Coolwarm } else { o.cmap },
        modes[0].1.signed,
        modes[0].1.abs_max,
        modes[0].1.vmax,
        &key_note,
    );

    let notes = modes.iter().map(|(_, m)| m.note.clone()).collect();
    Ok(Figure { canvas: cv, modes: notes, window: proj.win, at: None })
}

// ------------------------------------------------------------------------------- 4. Mesh3d implementation

const REGION_PALETTE: [[f32; 3]; 8] = [
    [0.23, 0.51, 0.96],
    [0.08, 0.72, 0.65],
    [0.96, 0.62, 0.04],
    [0.94, 0.33, 0.31],
    [0.66, 0.33, 0.97],
    [0.13, 0.77, 0.37],
    [0.93, 0.28, 0.60],
    [0.39, 0.45, 0.55],
];

#[allow(clippy::too_many_arguments)]
fn render_mesh_3d(
    mesh: &SimplexMesh,
    topo: &Mesh3dTopo,
    result: &SolveResult,
    o: &FigureOptions,
    shift: [f64; 3],
    repeat: usize,
    tot_lo: [f64; 3],
    tot_hi: [f64; 3],
) -> Result<Figure, Error> {
    let (planes, empty_clip) = match &o.clip {
        Some(cb) => active_half_spaces(cb, tot_lo, tot_hi),
        None => (Vec::new(), false),
    };
    if empty_clip {
        return Err(Error::invalid("the `clip` box excludes the entire 3D mesh"));
    }

    let (proj, aspect) = Proj3d::new(o.camera, tot_lo, tot_hi, o.window, (0.35, 1.4));
    let (pw, ph) = (o.panel_width, o.panel_width * aspect);
    let (iw, ih) = ((pw * 2.0).ceil() as u32, (ph * 2.0).ceil() as u32);
    let (iw_f, ih_f) = (f64::from(iw), f64::from(ih));

    let doms = mesh.domain_labels();
    let dom_col = |lab: i32| -> [f32; 3] {
        let pos = doms.binary_search(&lab).unwrap_or((lab.unsigned_abs() as usize) % REGION_PALETTE.len());
        REGION_PALETTE[pos % REGION_PALETTE.len()]
    };

    let key_light =
        normalize3(add3(add3(scale3(proj.right, -0.35), scale3(proj.up, 0.55)), scale3(proj.view_dir, -0.75)));
    let span = norm3(sub3(tot_hi, tot_lo)).max(1e-30);
    let cut_eps = 1e-13 * span;
    let inv_q = 1e6 / span;
    let qpt = |p: [f64; 3]| -> [i64; 3] {
        [(p[0] * inv_q).round() as i64, (p[1] * inv_q).round() as i64, (p[2] * inv_q).round() as i64]
    };

    let nc = mesh.n_cells();
    let cell_ids: Vec<usize> = (0..nc).collect();
    let collect_mesh_edges = o.show_edges == EdgeMode::Mesh;
    let curved = mesh.has_curved();

    let chunks = femcore::par::map_chunks(&cell_ids, 256, |chunk| {
        let mut opaque: Vec<Tri3d> = Vec::new();
        let mut raw_edges: Vec<QuantEdge> = Vec::new();
        let mut p_sh = [[0.0; 3]; 10];
        let mut clip_buf = Vec::with_capacity(8);
        let mut cut_buf = Vec::with_capacity(4);

        for &c in chunk {
            let dom = mesh.cell_labels[c];
            let base_c = dom_col(dom);
            let g0 = &topo.geom10[c];
            let has_ext = (0..4).any(|f| topo.face_neighbor[c * 4 + f] == u32::MAX);

            for m in 0..repeat {
                let sh = scale3(shift, m as f64);
                let n_chk = if curved { 10 } else { 4 };
                for a in 0..n_chk {
                    p_sh[a] = add3(g0[a], sh);
                }
                let mut outside = false;
                let mut straddles = false;
                for pl in &planes {
                    let (mut any_in, mut any_out) = (false, false);
                    for &pt in p_sh.iter().take(n_chk) {
                        if pl.signed_inside(pt) >= 0.0 {
                            any_in = true;
                        } else {
                            any_out = true;
                        }
                    }
                    if !any_in {
                        outside = true;
                        break;
                    }
                    if any_out {
                        straddles = true;
                    }
                }
                if outside || (!has_ext && !straddles) {
                    continue;
                }
                for a in n_chk..10 {
                    p_sh[a] = add3(g0[a], sh);
                }

                let mut emit_face = |tri: [PV; 3]| {
                    let fnrm = normalize3(cross3(sub3(tri[1].p, tri[0].p), sub3(tri[2].p, tri[0].p)));
                    let shade = (0.52 + 0.48 * dot3(fnrm, key_light).abs()) as f32;
                    let col = [
                        (base_c[0] * shade).clamp(0.0, 1.0),
                        (base_c[1] * shade).clamp(0.0, 1.0),
                        (base_c[2] * shade).clamp(0.0, 1.0),
                    ];
                    let p2 = [
                        proj.map_px(tri[0].p, iw_f, ih_f),
                        proj.map_px(tri[1].p, iw_f, ih_f),
                        proj.map_px(tri[2].p, iw_f, ih_f),
                    ];
                    opaque.push(Tri3d { p: p2, c: [col; 3], alpha: 1.0 });
                    if collect_mesh_edges {
                        for (ia, ib) in [(0, 1), (1, 2), (2, 0)] {
                            let (qa, qb) = (qpt(tri[ia].p), qpt(tri[ib].p));
                            if qa <= qb {
                                raw_edges.push((qa, qb, p2[ia], p2[ib]));
                            } else {
                                raw_edges.push((qb, qa, p2[ib], p2[ia]));
                            }
                        }
                    }
                };

                if has_ext {
                    for (f, &(fv, _, mids)) in TET_FACES.iter().enumerate() {
                        if topo.face_neighbor[c * 4 + f] != u32::MAX {
                            continue;
                        }
                        // Use 4 sub-tris only when curved, else the 1 parent tet face so mesh edges show true tet edges!
                        let sub_tris: &[[usize; 3]] = if curved && result.order == 2 {
                            &[
                                [fv[0], mids[0], mids[1]],
                                [mids[0], fv[1], mids[2]],
                                [mids[1], mids[2], fv[2]],
                                [mids[0], mids[2], mids[1]],
                            ]
                        } else {
                            &[[fv[0], fv[1], fv[2]]]
                        };
                        for st in sub_tris {
                            let base_tri = [
                                PV { p: p_sh[st[0]], v: 0.0 },
                                PV { p: p_sh[st[1]], v: 0.0 },
                                PV { p: p_sh[st[2]], v: 0.0 },
                            ];
                            if !straddles {
                                emit_face(base_tri);
                            } else {
                                clip_buf.clear();
                                clip_poly_against_planes(base_tri, &planes, None, &mut clip_buf);
                                for ct in &clip_buf {
                                    emit_face(*ct);
                                }
                            }
                        }
                    }
                }

                if straddles {
                    for (pi, pl) in planes.iter().enumerate() {
                        let tet = [
                            PV { p: p_sh[0], v: 0.0 },
                            PV { p: p_sh[1], v: 0.0 },
                            PV { p: p_sh[2], v: 0.0 },
                            PV { p: p_sh[3], v: 0.0 },
                        ];
                        let dt = [
                            pl.signed_inside(tet[0].p),
                            pl.signed_inside(tet[1].p),
                            pl.signed_inside(tet[2].p),
                            pl.signed_inside(tet[3].p),
                        ];
                        cut_buf.clear();
                        slice_subtet(tet, dt, cut_eps, pl.out_normal, &mut cut_buf);
                        for ct_tri in &cut_buf {
                            if planes.len() == 1 {
                                emit_face(*ct_tri);
                            } else {
                                clip_buf.clear();
                                clip_poly_against_planes(*ct_tri, &planes, Some(pi), &mut clip_buf);
                                for clipped in &clip_buf {
                                    emit_face(*clipped);
                                }
                            }
                        }
                    }
                }
            }
        }
        (opaque, raw_edges)
    });

    let mut opaque = Vec::new();
    let mut mesh_edge_raw = Vec::new();
    for (op, ed) in chunks {
        opaque.extend(op);
        mesh_edge_raw.extend(ed);
    }

    let edges_3d: Vec<([f64; 3], [f64; 3], Rgba, f64)> = match o.show_edges {
        EdgeMode::None => Vec::new(),
        EdgeMode::Mesh => {
            mesh_edge_raw.sort_unstable_by_key(|e| (e.0, e.1));
            mesh_edge_raw.dedup_by_key(|e| (e.0, e.1));
            let stroke: Rgba = [15, 23, 42, 200];
            mesh_edge_raw.into_iter().map(|(_, _, a, b)| (a, b, stroke, 1.1)).collect()
        }
        EdgeMode::Regions => {
            let mut segs = Vec::new();
            let stroke: Rgba = [15, 23, 42, 220];
            for fe in &topo.ext_feature_edges {
                let g0 = &topo.geom10[fe.cell as usize];
                for m in 0..repeat {
                    let sh = scale3(shift, m as f64);
                    if let Some((pa, pb)) =
                        clip_seg3d(add3(g0[fe.la as usize], sh), add3(g0[fe.lb as usize], sh), &planes)
                    {
                        segs.push((proj.map_px(pa, iw_f, ih_f), proj.map_px(pb, iw_f, ih_f), stroke, 1.4));
                    }
                }
            }
            segs
        }
    };

    let rgba = raster_field_3d(iw, ih, &opaque, &[], &edges_3d);
    let (left, head_h, top_title, foot) = (18.0, 28.0, 48.0, 56.0);
    let summary = format!(
        "{} tets, {} domains{}",
        nc * repeat,
        doms.len(),
        if repeat > 1 { format!(" ({repeat}x periodic)") } else { String::new() }
    );
    let width = (left + pw + 18.0).max(left + text_width(&summary, 12.0) + 180.0);
    let height = top_title + head_h + ph + foot;
    let mut cv = Canvas::new(width, height);
    cv.bold(left, 30.0, &format!("{}: 3D mesh ({summary})", result.physics), 18.0, INK, Anchor::Start);

    let y = top_title + head_h;
    cv.push(Item::Rect { x: left, y, w: pw, h: ph, fill: Some(PANEL_BG), stroke: None });
    cv.push(Item::Image { x: left, y, w: pw, h: ph, px: iw, py: ih, rgba });
    draw_triad(&mut cv, o.camera, left + 26.0, y + ph - 24.0, 15.0);
    cv.push(Item::Rect { x: left, y, w: pw, h: ph, fill: None, stroke: Some(([71, 85, 105, 255], 0.8)) });

    // Region legend swatches
    let mut lx = left;
    let ly = y + ph + 16.0;
    for &d in doms.iter().take(8) {
        let c = bytes(dom_col(d));
        cv.push(Item::Rect { x: lx, y: ly, w: 12.0, h: 12.0, fill: Some(c), stroke: Some(([71, 85, 105, 255], 0.6)) });
        let txt = format!("domain {d}");
        cv.text(lx + 16.0, ly + 10.0, &txt, 10.5, MUTED, Anchor::Start);
        lx += 24.0 + text_width(&txt, 10.5);
    }

    Ok(Figure { canvas: cv, modes: Vec::new(), window: proj.win, at: None })
}
