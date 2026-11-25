// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Absorbing layers across the two halves of a project: the geometry cuts or extrudes a region into the physical part and one piece per
//! side or shell (`Pml` feature), and the physics stretches the coordinates in each piece (`pml` or `infinite_element` condition).
//! `add_pml_conditions` writes one condition per piece from the geometry, with joint multi-material wavelength auto-calibration so
//! adjacent regions sharing a PML slab/shell use identical coordinate stretching.

use crate::engine::Engine;
use femcore::material::{isotropic_stiffness, rotate_m6_global, M6};
use femcore::mesh::cells::CellType;
use femcore::meshing::physics_size::backward_wave_mpml_ratio;
use femcore::pml::auto_scale_params;
use geom::model::{PmlCoordKind, PmlGrowMode, PmlRole};
use geom::sim::PmlRegion;
use nyx_physics::{BcSpec, Error, MaterialAssign, Mesh};
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};

fn dir_name(s: i8) -> &'static str {
    match s {
        0 => "off",
        s if s < 0 => "minus",
        _ => "plus",
    }
}

fn coord_name(c: PmlCoordKind) -> &'static str {
    match c {
        PmlCoordKind::Auto | PmlCoordKind::Cartesian => "cartesian",
        PmlCoordKind::Cylindrical => "cylindrical",
        PmlCoordKind::Spherical => "spherical",
    }
}

fn mode_name(m: PmlGrowMode) -> &'static str {
    match m {
        PmlGrowMode::Inward => "inward",
        PmlGrowMode::Outward => "outward",
    }
}

fn cartesian_piece_bounds(r: &PmlRegion, d_actual: f64, dim2: bool) -> ([f64; 3], [f64; 3]) {
    let mut lo = r.window.lo;
    let mut hi = r.window.hi;
    for a in 0..3 {
        if dim2 && a == 1 {
            lo[1] = 0.0;
            hi[1] = 0.0;
            continue;
        }
        match (r.mode, r.sides[a]) {
            (PmlGrowMode::Outward, s) if s > 0 => {
                lo[a] = r.window.hi[a];
                hi[a] = r.window.hi[a] + d_actual;
            }
            (PmlGrowMode::Outward, s) if s < 0 => {
                lo[a] = r.window.lo[a] - d_actual;
                hi[a] = r.window.lo[a];
            }
            (PmlGrowMode::Inward, s) if s > 0 => {
                lo[a] = r.window.hi[a];
                hi[a] = r.window.hi[a] + d_actual;
            }
            (PmlGrowMode::Inward, s) if s < 0 => {
                lo[a] = r.window.lo[a] - d_actual;
                hi[a] = r.window.lo[a];
            }
            _ => {}
        }
    }
    (lo, hi)
}

fn piece_outward_normals(r: &PmlRegion, dim2: bool) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    match r.coord {
        PmlCoordKind::Auto | PmlCoordKind::Cartesian => {
            if r.sides[0] != 0 {
                out.push([r.sides[0] as f64, 0.0, 0.0]);
            }
            if !dim2 && r.sides[1] != 0 {
                out.push([0.0, r.sides[1] as f64, 0.0]);
            }
            if r.sides[2] != 0 {
                out.push([0.0, 0.0, r.sides[2] as f64]);
            }
        }
        PmlCoordKind::Cylindrical => {
            let ax = r.axis;
            let an = (ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]).sqrt();
            let a_hat = if an > 1e-30 {
                [ax[0] / an, ax[1] / an, ax[2] / an]
            } else if dim2 {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, 1.0]
            };
            if matches!(r.role, PmlRole::EndCap | PmlRole::CornerAnnulus) || r.sides[2] != 0 {
                let s = if r.sides[2] < 0 { -1.0 } else { 1.0 };
                out.push([s * a_hat[0], s * a_hat[1], s * a_hat[2]]);
            }
            if matches!(r.role, PmlRole::RadialShell | PmlRole::CornerAnnulus) || r.sides[0] != 0 {
                for k in 0..8 {
                    let th = 2.0 * std::f64::consts::PI * (k as f64) / 8.0;
                    let (s, c) = th.sin_cos();
                    if dim2 {
                        out.push([c, 0.0, s]);
                    } else {
                        out.push([c, s, 0.0]);
                    }
                }
            }
        }
        PmlCoordKind::Spherical => {
            out.extend_from_slice(&[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ]);
        }
    }
    if out.is_empty() {
        out.push([1.0, 0.0, 0.0]);
    }
    out
}

impl Engine {
    /// Adds a `pml` (or `infinite_element` for static studies) condition for every absorbing piece of the geometry to `physics`
    /// (replacing existing absorbing conditions on those pieces). When `strength <= 0.0`, parameters are jointly auto-calibrated
    /// across adjacent materials sharing each PML slab/shell from the study's wavelengths.
    pub fn add_pml_conditions(&mut self, physics: &str, strength: f64) -> Result<Value, Error> {
        let s_opt = (strength > 0.0).then_some(strength);
        self.add_pml_conditions_with(physics, s_opt, None, None, None)
    }

    /// Full absorbing-condition generator supporting explicit or auto-calibrated `strength`, `profile` (`"auto"`, `"kappa_polynomial"`,
    /// `"rational"`, `"exponential"`), `auto_scale`, and `target_reflection`.
    pub fn add_pml_conditions_with(
        &mut self,
        physics: &str,
        strength: Option<f64>,
        profile: Option<&str>,
        auto_scale: Option<bool>,
        target_reflection: Option<f64>,
    ) -> Result<Value, Error> {
        let p = crate::physics_for(&self.registry, physics)?;
        let spec_before = self.project.spec(physics).cloned();
        let study_kind = spec_before.as_ref().map(|s| s.study.kind.as_str()).unwrap_or("");
        let use_infinite_element = physics == "dfs"
            || (physics == "sms" && study_kind == "static")
            || (!p.bc_kinds().iter().any(|k| k.kind == "pml")
                && p.bc_kinds().iter().any(|k| k.kind == "infinite_element"));
        let bc_kind = if use_infinite_element { "infinite_element" } else { "pml" };
        let Some(bc_desc) = p.bc_kinds().into_iter().find(|k| k.kind == bc_kind) else {
            return Err(Error::invalid(format!("{} has no absorbing layers (no `{bc_kind}` condition)", p.name())));
        };
        let allowed_keys: BTreeSet<String> = bc_desc.params.into_iter().map(|pd| pd.name).collect();

        let ev = femcore::meshing::geomdrive::evaluate(&self.project.geometry).map_err(Error::Invalid)?;
        let labels = ev.label_map().map_err(|e| Error::invalid(e.to_string()))?;
        if ev.pml_regions().is_empty() {
            return Err(Error::invalid("the geometry has no absorbing layers; add them with add_pml first"));
        }

        let dim2 = self.project.geometry.dim == 2;
        let want_auto = auto_scale.unwrap_or_else(|| strength.is_none_or(|s| s <= 0.0));
        let r0 = target_reflection.unwrap_or(1e-6).clamp(1e-15, 0.5);
        let is_eigen =
            !matches!(study_kind, "driven" | "static" | "electrostatic" | "heat_stationary" | "magnetostatic");
        let resolved_profile = match profile.map(str::trim).filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("auto"))
        {
            Some(pr) => pr.to_ascii_lowercase(),
            None if use_infinite_element => "rational".to_string(),
            None if !is_eigen => "rational".to_string(),
            None => "kappa_polynomial".to_string(),
        };
        let is_rational = resolved_profile == "rational";

        // Collect per-domain wavelengths when a problem spec exists for this physics.
        let lib = self.library().ok();
        let domain_names: Vec<&str> = labels.domains.iter().map(|d| d.name.as_str()).collect();
        let expanded_wavelengths: BTreeMap<String, f64> = spec_before
            .as_ref()
            .map(|s| p.sizing(s, lib.as_ref()).expand_wildcards(&domain_names).wavelengths)
            .unwrap_or_default();

        let region_wavelength = |r: &PmlRegion| -> Option<f64> {
            let dom_name = &labels.domains[r.region as usize].name;
            if let Some(&w) = expanded_wavelengths.get(dom_name) {
                if w.is_finite() && w > 0.0 {
                    return Some(w);
                }
            }
            if let Some(parent_idx) = r.parent {
                if let Some(parent_dom) = labels.domains.get(parent_idx as usize) {
                    if let Some(&w) = expanded_wavelengths.get(&parent_dom.name) {
                        if w.is_finite() && w > 0.0 {
                            return Some(w);
                        }
                    }
                }
            }
            None
        };

        // Group adjacent materials sharing the same PML slab/shell `(coord, sides, role)` so their
        // coordinate stretching parameters `(strength, order, profile, kappa_max/curvature, bounds)` are identical.
        type GroupKey = (u8, [i8; 3], u8);
        let group_key = |r: &PmlRegion| -> GroupKey {
            let c_code = match r.coord {
                PmlCoordKind::Auto | PmlCoordKind::Cartesian => 0,
                PmlCoordKind::Cylindrical => 1,
                PmlCoordKind::Spherical => 2,
            };
            let r_code = r.role as u8;
            (c_code, r.sides, r_code)
        };

        let mut group_lambdas: BTreeMap<GroupKey, (f64, f64)> = BTreeMap::new();
        for r in ev.pml_regions() {
            let key = group_key(r);
            let entry = group_lambdas.entry(key).or_insert((f64::INFINITY, 0.0));
            if let Some(w) = region_wavelength(r) {
                entry.0 = entry.0.min(w);
                entry.1 = entry.1.max(w);
            }
        }
        let all_min =
            expanded_wavelengths.values().copied().filter(|&w| w.is_finite() && w > 0.0).fold(f64::INFINITY, f64::min);
        let all_max =
            expanded_wavelengths.values().copied().filter(|&w| w.is_finite() && w > 0.0).fold(0.0_f64, f64::max);

        let mut pieces: Vec<(String, Value)> = Vec::with_capacity(ev.pml_regions().len());
        for r in ev.pml_regions() {
            let name = labels.domains[r.region as usize].name.clone();
            let y_dir = if dim2 { "off" } else { dir_name(r.sides[1]) };
            let key = group_key(r);
            let (mut lmin, mut lmax) = group_lambdas.get(&key).copied().unwrap_or((f64::INFINITY, 0.0));
            if !lmin.is_finite() || lmax <= 0.0 {
                if all_min.is_finite() && all_max > 0.0 {
                    lmin = all_min;
                    lmax = all_max;
                } else {
                    let fallback = 2.0 * r.thickness.abs().max(1e-6);
                    lmin = fallback;
                    lmax = fallback;
                }
            }
            let d_opt = (r.thickness > 0.0 && r.thickness.is_finite()).then_some(r.thickness);
            let cal = auto_scale_params(lmin, lmax, d_opt, r.grading, r0, is_eigen, is_rational);
            let d_actual = d_opt.unwrap_or(cal.d_auto);
            let strength_val = match strength {
                Some(s) if s > 0.0 && !want_auto => s,
                _ => cal.psi_max,
            };

            let is_legacy_cartesian_inward = !use_infinite_element
                && !want_auto
                && profile.is_none()
                && target_reflection.is_none()
                && matches!(r.coord, PmlCoordKind::Auto | PmlCoordKind::Cartesian)
                && r.mode == PmlGrowMode::Inward;

            if is_legacy_cartesian_inward {
                pieces.push((
                    name,
                    json!({
                        "x": dir_name(r.sides[0]),
                        "y": y_dir,
                        "z": dir_name(r.sides[2]),
                        "strength": strength_val,
                        "order": r.grading,
                        "re": 1.0
                    }),
                ));
                continue;
            }

            let mut m = Map::new();
            let mut insert_if_allowed = |k: &str, v: Value| {
                if allowed_keys.contains(k) {
                    m.insert(k.to_string(), v);
                }
            };

            insert_if_allowed("coord", json!(coord_name(r.coord)));
            insert_if_allowed("x", json!(dir_name(r.sides[0])));
            insert_if_allowed("y", json!(y_dir));
            insert_if_allowed("z", json!(dir_name(r.sides[2])));

            match r.coord {
                PmlCoordKind::Auto | PmlCoordKind::Cartesian => {
                    if r.mode == PmlGrowMode::Outward || want_auto {
                        let (blo, bhi) = cartesian_piece_bounds(r, d_actual, dim2);
                        insert_if_allowed("bounds_lo", json!(blo));
                        insert_if_allowed("bounds_hi", json!(bhi));
                    }
                }
                PmlCoordKind::Cylindrical => {
                    let radial_on = matches!(r.role, PmlRole::RadialShell | PmlRole::CornerAnnulus)
                        || r.sides[0] != 0
                        || (!dim2 && r.sides[1] != 0);
                    let axial_side = if matches!(r.role, PmlRole::EndCap | PmlRole::CornerAnnulus) || r.sides[2] != 0 {
                        r.sides[2]
                    } else {
                        0
                    };
                    let z_in = if axial_side >= 0 { r.z_range[1] } else { r.z_range[0] };
                    insert_if_allowed("center", json!(r.center));
                    insert_if_allowed("axis", json!(r.axis));
                    insert_if_allowed("r_in", json!(r.r_in));
                    insert_if_allowed("r_thickness", json!(d_actual));
                    insert_if_allowed("thickness", json!(d_actual));
                    insert_if_allowed("radial_on", json!(radial_on));
                    insert_if_allowed("axial", json!(dir_name(axial_side)));
                    insert_if_allowed("z_in", json!(z_in));
                    insert_if_allowed("z_thickness", json!(d_actual));
                }
                PmlCoordKind::Spherical => {
                    insert_if_allowed("center", json!(r.center));
                    insert_if_allowed("r_in", json!(r.r_in));
                    insert_if_allowed("thickness", json!(d_actual));
                }
            }

            if use_infinite_element {
                insert_if_allowed("profile", json!(resolved_profile));
            } else {
                insert_if_allowed("strength", json!(strength_val));
                insert_if_allowed("order", json!(r.grading));
                insert_if_allowed("re", json!(1.0));
                insert_if_allowed("profile", json!(resolved_profile));
                if is_rational {
                    insert_if_allowed("curvature", json!(cal.rational_curvature()));
                } else {
                    insert_if_allowed("kappa_max", json!(cal.kappa_max));
                }
                if want_auto {
                    insert_if_allowed("auto_scale", json!(true));
                    insert_if_allowed("target_reflection", json!(r0));
                }
            }

            pieces.push((name, Value::Object(m)));
        }

        let spec = self.project.need_spec(physics)?;
        spec.bcs.retain(|b| {
            !(matches!(b.kind.as_str(), "pml" | "infinite_element") && pieces.iter().any(|(n, _)| *n == b.on))
        });
        for (name, params) in &pieces {
            spec.bcs.push(BcSpec { on: name.clone(), kind: bc_kind.into(), params: params.clone() });
        }
        Ok(json!({
            "kind": bc_kind,
            "profile": resolved_profile,
            "auto_scale": want_auto,
            "added": pieces.iter().map(|(n, _)| n.clone()).collect::<Vec<_>>()
        }))
    }

    /// Pre-solve and post-mesh absorbing-layer diagnostics: uncut element counts, minimum PML scaled Jacobian,
    /// layer count warnings, normalized penetration depth `sDist` bounds `[0.0, 1.0]`, and Christoffel
    /// backward-wave check (`has_backward_wave`, `min_projection`, `recommended_mpml_ratio`).
    pub fn pml_diagnostics(&mut self) -> Result<Value, Error> {
        let ev = femcore::meshing::geomdrive::evaluate(&self.project.geometry).map_err(Error::Invalid)?;
        let labels = ev.label_map().map_err(|e| Error::invalid(e.to_string()))?;
        let dim2 = self.project.geometry.dim == 2;

        if self.model().is_none() && !self.project.geometry.features.is_empty() {
            let _ = self.mesh_now();
        }

        let mut element_counts: BTreeMap<String, usize> = BTreeMap::new();
        let mut min_pml_scaled_jacobian: Option<f64> = None;
        let mut pml_warnings: Vec<String> = Vec::new();

        if let Some(m) = self.model() {
            if let Some(rep) = m.report() {
                element_counts = rep.element_counts.clone();
                min_pml_scaled_jacobian = rep.min_pml_scaled_jacobian;
                pml_warnings = rep.pml_warnings.clone();
            }
            if element_counts.is_empty() {
                let counts = match m.mesh().as_ref() {
                    Mesh::D2(cm) => cm.counts(),
                    Mesh::D3(cm) => cm.counts(),
                };
                for (ct, cnt) in counts {
                    if cnt > 0 {
                        let key = match ct {
                            CellType::Tri => "tri",
                            CellType::Quad => "quad",
                            CellType::Tet => "tet",
                            CellType::Prism => "prism",
                            CellType::Hex => "hex",
                        };
                        element_counts.insert(key.into(), cnt);
                    }
                }
            }
        }
        if pml_warnings.is_empty() {
            for r in ev.pml_regions() {
                let n = (r.layers.round() as usize).max(1);
                if r.layers > 0.0 && n < 5 {
                    pml_warnings.push(format!(
                        "absorbing layer uses {n} (< 5) layer(s); at least 5-6 layers are recommended to suppress discretization reflection"
                    ));
                    break;
                }
            }
        }

        // Resolve elastic materials (if `sms` is loaded) for Christoffel backward-wave stability diagnostics.
        let lib = self.library().ok();
        let sms_spec = self.project.spec("sms").cloned();

        let find_elastic = |r: &PmlRegion| -> Option<(M6, f64)> {
            let spec = sms_spec.as_ref()?;
            let try_dom = |reg_idx: u32| -> Option<(M6, f64)> {
                let d_name = &labels.domains.get(reg_idx as usize)?.name;
                let (_, assign) = spec.materials.iter().find(|(region, _)| {
                    *region == d_name
                        || region
                            .strip_suffix('*')
                            .is_some_and(|base| d_name == base || d_name.starts_with(&format!("{base}.")))
                })?;
                match assign {
                    MaterialAssign::Named(mat_name) => {
                        let ms = lib.as_ref()?.get(mat_name)?;
                        let rho = ms.density?;
                        let el = ms.elastic.as_ref()?;
                        Some((rotate_m6_global(&el.c, [0.0, 0.0, 0.0]), rho))
                    }
                    MaterialAssign::Oriented { material, euler_deg } => {
                        let ms = lib.as_ref()?.get(material)?;
                        let rho = ms.density?;
                        let el = ms.elastic.as_ref()?;
                        Some((rotate_m6_global(&el.c, *euler_deg), rho))
                    }
                    MaterialAssign::Props(map) => {
                        let rho = *map.get("rho")?;
                        let e = *map.get("E")?;
                        let nu = *map.get("nu")?;
                        Some((isotropic_stiffness(e, nu), rho))
                    }
                    MaterialAssign::PropsTensor(_) => None,
                }
            };
            try_dom(r.region).or_else(|| r.parent.and_then(try_dom))
        };

        let mut global_has_backward = false;
        let mut global_min_proj = 1.0_f64;
        let mut global_mpml_ratio = 0.0_f64;
        let mut domain_reports = Vec::with_capacity(ev.pml_regions().len());

        for r in ev.pml_regions() {
            let name = labels.domains[r.region as usize].name.clone();
            let mut dom_backward = false;
            let mut dom_min_proj = 1.0_f64;
            let mut dom_mpml = 0.0_f64;
            if let Some((c, rho)) = find_elastic(r) {
                for n in piece_outward_normals(r, dim2) {
                    let chk = backward_wave_mpml_ratio(&c, rho, n);
                    dom_backward |= chk.has_backward_wave;
                    dom_min_proj = dom_min_proj.min(chk.min_projection);
                    dom_mpml = dom_mpml.max(chk.recommended_mpml_ratio);
                }
            }
            global_has_backward |= dom_backward;
            global_min_proj = global_min_proj.min(dom_min_proj);
            global_mpml_ratio = global_mpml_ratio.max(dom_mpml);

            domain_reports.push(json!({
                "name": name,
                "coord": coord_name(r.coord),
                "mode": mode_name(r.mode),
                "thickness": r.thickness,
                "layers": r.layers,
                "layer_grading": r.layer_grading,
                "s_dist_bounds": [0.0, 1.0],
                "has_backward_wave": dom_backward,
                "min_projection": dom_min_proj,
                "recommended_mpml_ratio": dom_mpml,
            }));
        }

        let s_dist_bounds = if ev.pml_regions().is_empty() { json!(null) } else { json!([0.0, 1.0]) };

        Ok(json!({
            "element_counts": element_counts,
            "min_pml_scaled_jacobian": min_pml_scaled_jacobian,
            "pml_warnings": pml_warnings,
            "s_dist_bounds": s_dist_bounds,
            "sDist_bounds": s_dist_bounds,
            "has_backward_wave": global_has_backward,
            "min_projection": global_min_proj,
            "recommended_mpml_ratio": global_mpml_ratio,
            "domains": domain_reports,
        }))
    }

    /// The physics solved when none is named: the active one and, before it, those of the project that export what it imports (a heat
    /// result feeding thermal expansion).
    pub fn default_targets(&self) -> Vec<String> {
        let Some(active) =
            self.project.active.clone().or_else(|| self.project.physics.first().map(|s| s.physics.clone()))
        else {
            return vec![];
        };
        let mut out = vec![active.clone()];
        if let Some(p) = self.registry.get(&active) {
            for imp in p.imports() {
                for s in &self.project.physics {
                    if s.physics != active
                        && !out.contains(&s.physics)
                        && self.registry.get(&s.physics).is_some_and(|q| q.exports().iter().any(|e| e.role == imp.role))
                    {
                        out.push(s.physics.clone());
                    }
                }
            }
        }
        out
    }
}
