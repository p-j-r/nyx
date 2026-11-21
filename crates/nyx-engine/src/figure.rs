// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The `figure` command (ADR 99): a figure of the modes or 3D fields of the latest result of a physics, written as SVG or PNG. The window, a script
//! and the MCP server all call it; the drawing is `nyx_analysis::figure`, which names no physics.

use crate::engine::Engine;
use crate::meshing::MeshedModel;
use nyx_analysis::figure::{mode_figure_on, Camera3d, ClipBox, Cmap, EdgeMode, FigureOptions, Kind, PlaneSpec};
use nyx_physics::{Error, LabelResolver, Scalar, SolveResult};
use serde_json::{json, Value};
use std::path::Path;

fn numbers(a: &Value, k: &str) -> Result<Option<Vec<f64>>, Error> {
    match a.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(v)) => v
            .iter()
            .map(|x| {
                x.as_f64()
                    .filter(|n| n.is_finite())
                    .ok_or_else(|| Error::invalid(format!("`{k}` must be a list of finite numbers")))
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some),
        Some(_) => Err(Error::invalid(format!("`{k}` must be a list of numbers"))),
    }
}

fn whole(a: &Value, k: &str) -> Result<Option<usize>, Error> {
    match a.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => {
            v.as_u64().map(|x| Some(x as usize)).ok_or_else(|| Error::invalid(format!("`{k}` must be a whole number")))
        }
    }
}

fn real(a: &Value, k: &str) -> Result<Option<f64>, Error> {
    match a.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => v
            .as_f64()
            .filter(|n| n.is_finite())
            .map(Some)
            .ok_or_else(|| Error::invalid(format!("`{k}` must be a finite number"))),
    }
}

fn opt_real_val(v: Option<&Value>, what: &str) -> Result<Option<f64>, Error> {
    match v {
        None | Some(Value::Null) => Ok(None),
        Some(x) => x
            .as_f64()
            .filter(|n| n.is_finite())
            .map(Some)
            .ok_or_else(|| Error::invalid(format!("`{what}` must be a finite number or null"))),
    }
}

fn vec3(v: &Value, what: &str) -> Result<[f64; 3], Error> {
    let arr = v.as_array().ok_or_else(|| Error::invalid(format!("`{what}` must be [x, y, z]")))?;
    if arr.len() != 3 {
        return Err(Error::invalid(format!("`{what}` must have 3 numbers [x, y, z]")));
    }
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[i] = arr[i]
            .as_f64()
            .filter(|n| n.is_finite())
            .ok_or_else(|| Error::invalid(format!("`{what}` entries must be finite numbers")))?;
    }
    Ok(out)
}

fn parse_axis_char(s: &str) -> Result<(char, usize), Error> {
    match s.trim().to_ascii_lowercase().as_str() {
        "x" => Ok(('x', 0)),
        "y" => Ok(('y', 1)),
        "z" => Ok(('z', 2)),
        other => Err(Error::invalid(format!("plane `axis` must be 'x', 'y' or 'z', not '{other}'"))),
    }
}

fn parse_plane(args: &Value, at_fallback: Option<f64>, lo: [f64; 3], hi: [f64; 3]) -> Result<Option<PlaneSpec>, Error> {
    match args.get("plane") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) => {
            let (axis, idx) = parse_axis_char(s)?;
            let at = at_fallback.unwrap_or(0.5 * (lo[idx] + hi[idx]));
            Ok(Some(PlaneSpec::Axis { axis, at }))
        }
        Some(Value::Object(m)) => {
            if let Some(ax_val) = m.get("axis") {
                let s = ax_val.as_str().ok_or_else(|| Error::invalid("`plane.axis` must be 'x', 'y' or 'z'"))?;
                let (axis, idx) = parse_axis_char(s)?;
                let at = opt_real_val(m.get("at"), "plane.at")?
                    .unwrap_or_else(|| at_fallback.unwrap_or(0.5 * (lo[idx] + hi[idx])));
                return Ok(Some(PlaneSpec::Axis { axis, at }));
            }
            let orig_val = m
                .get("origin")
                .or_else(|| m.get("point"))
                .ok_or_else(|| Error::invalid("`plane` object requires either `axis` or `origin` and `normal`"))?;
            let norm_val = m
                .get("normal")
                .ok_or_else(|| Error::invalid("general `plane` object requires `normal: [nx, ny, nz]`"))?;
            let origin = vec3(orig_val, "plane.origin")?;
            let normal = vec3(norm_val, "plane.normal")?;
            if normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2) <= 1e-30 {
                return Err(Error::invalid("`plane.normal` must be non-zero"));
            }
            Ok(Some(PlaneSpec::General { origin, normal }))
        }
        Some(_) => Err(Error::invalid("`plane` must be `{axis, at}` or `{origin, normal}`")),
    }
}

fn parse_clip(args: &Value) -> Result<Option<ClipBox>, Error> {
    match args.get("clip") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(arr)) => {
            if arr.len() != 6 {
                return Err(Error::invalid("`clip` array must be [xmin, xmax, ymin, ymax, zmin, zmax]"));
            }
            let cb = ClipBox {
                xmin: opt_real_val(arr.first(), "clip[0]")?,
                xmax: opt_real_val(arr.get(1), "clip[1]")?,
                ymin: opt_real_val(arr.get(2), "clip[2]")?,
                ymax: opt_real_val(arr.get(3), "clip[3]")?,
                zmin: opt_real_val(arr.get(4), "clip[4]")?,
                zmax: opt_real_val(arr.get(5), "clip[5]")?,
            };
            check_clip_bounds(&cb)?;
            Ok(Some(cb))
        }
        Some(Value::Object(m)) => {
            let mut cb = ClipBox {
                xmin: opt_real_val(m.get("xmin").or_else(|| m.get("x_min")), "clip.xmin")?,
                xmax: opt_real_val(m.get("xmax").or_else(|| m.get("x_max")), "clip.xmax")?,
                ymin: opt_real_val(m.get("ymin").or_else(|| m.get("y_min")), "clip.ymin")?,
                ymax: opt_real_val(m.get("ymax").or_else(|| m.get("y_max")), "clip.ymax")?,
                zmin: opt_real_val(m.get("zmin").or_else(|| m.get("z_min")), "clip.zmin")?,
                zmax: opt_real_val(m.get("zmax").or_else(|| m.get("z_max")), "clip.zmax")?,
            };
            for (k, lo_ref, hi_ref) in [
                ("x", &mut cb.xmin, &mut cb.xmax),
                ("y", &mut cb.ymin, &mut cb.ymax),
                ("z", &mut cb.zmin, &mut cb.zmax),
            ] {
                if let Some(v) = m.get(k) {
                    let pair = v.as_array().ok_or_else(|| Error::invalid(format!("`clip.{k}` must be [min, max]")))?;
                    if pair.len() != 2 {
                        return Err(Error::invalid(format!("`clip.{k}` must have 2 elements [min, max]")));
                    }
                    if let Some(a) = opt_real_val(pair.first(), k)? {
                        *lo_ref = Some(a);
                    }
                    if let Some(b) = opt_real_val(pair.get(1), k)? {
                        *hi_ref = Some(b);
                    }
                }
            }
            check_clip_bounds(&cb)?;
            Ok(Some(cb))
        }
        Some(_) => Err(Error::invalid("`clip` must be an object `{xmin, xmax, ymin, ymax, zmin, zmax}`")),
    }
}

fn check_clip_bounds(cb: &ClipBox) -> Result<(), Error> {
    for (name, lo, hi) in [("x", cb.xmin, cb.xmax), ("y", cb.ymin, cb.ymax), ("z", cb.zmin, cb.zmax)] {
        if let (Some(a), Some(b)) = (lo, hi) {
            if a >= b {
                return Err(Error::invalid(format!("`clip`: {name}min ({a}) must be smaller than {name}max ({b})")));
            }
        }
    }
    Ok(())
}

fn parse_camera(args: &Value) -> Result<Camera3d, Error> {
    let mut cam = Camera3d::default();
    match args.get("camera") {
        None | Some(Value::Null) => Ok(cam),
        Some(Value::Array(arr)) => {
            if arr.len() != 2 {
                return Err(Error::invalid("`camera` array must be [azimuth_deg, elevation_deg]"));
            }
            cam.azimuth_deg = arr[0]
                .as_f64()
                .filter(|n| n.is_finite())
                .ok_or_else(|| Error::invalid("`camera[0]` must be a finite number"))?;
            cam.elevation_deg = arr[1]
                .as_f64()
                .filter(|n| n.is_finite())
                .ok_or_else(|| Error::invalid("`camera[1]` must be a finite number"))?;
            Ok(cam)
        }
        Some(Value::Object(m)) => {
            if let Some(az) = opt_real_val(m.get("azimuth_deg").or_else(|| m.get("azimuth")), "camera.azimuth_deg")? {
                cam.azimuth_deg = az;
            }
            if let Some(el) =
                opt_real_val(m.get("elevation_deg").or_else(|| m.get("elevation")), "camera.elevation_deg")?
            {
                cam.elevation_deg = el;
            }
            Ok(cam)
        }
        Some(_) => Err(Error::invalid("`camera` must be `{azimuth_deg, elevation_deg}`")),
    }
}

#[allow(clippy::too_many_arguments)]
fn resolve_periodic_shift_and_phase(
    engine: &Engine,
    physics: &str,
    model: &MeshedModel,
    solver_mesh: &femcore::mesh::SimplexMesh,
    result: &SolveResult,
    first_mode_1based: usize,
    explicit_shift: Option<[f64; 3]>,
    explicit_phase: Option<f64>,
) -> ([f64; 3], f64) {
    let shift = explicit_shift.unwrap_or_else(|| {
        let pairs = model.evaluated().periodic_pairs();
        let chosen = pairs
            .iter()
            .find(|p| {
                let n = (p.translation[0].powi(2) + p.translation[1].powi(2) + p.translation[2].powi(2)).sqrt();
                n > 0.0 && p.translation[1].abs() > 0.9 * n
            })
            .or_else(|| {
                pairs
                    .iter()
                    .find(|p| p.translation[0].powi(2) + p.translation[1].powi(2) + p.translation[2].powi(2) > 0.0)
            });
        if let Some(pair) = chosen {
            let t = pair.translation;
            // Orient along + of dominant axis so y-periodic tiling advances in +y
            let dom = if t[1].abs() >= t[0].abs() && t[1].abs() >= t[2].abs() {
                t[1]
            } else if t[0].abs() >= t[2].abs() {
                t[0]
            } else {
                t[2]
            };
            if dom < 0.0 {
                [-t[0], -t[1], -t[2]]
            } else {
                t
            }
        } else {
            let (lo, hi) = solver_mesh.bbox();
            [0.0, (hi[1] - lo[1]).max(0.0), 0.0]
        }
    });

    if let Some(ph) = explicit_phase {
        return (shift, ph);
    }

    let lambda = (shift[0].powi(2) + shift[1].powi(2) + shift[2].powi(2)).sqrt();
    let exact_scalar = |name: &str| -> Option<f64> {
        result
            .scalars
            .iter()
            .find(|(n, _)| n == name)
            .and_then(|(_, v)| if let Scalar::Num(x) = v { Some(*x) } else { None })
    };

    // 1. SolveResult scalar `phase N (rad)`
    if let Some(ph) = exact_scalar(&format!("phase {first_mode_1based} (rad)")) {
        return (shift, ph);
    }
    // 2. SolveResult scalar `k_y re N (rad/m) * Lambda`
    if let Some(ky) = exact_scalar(&format!("k_y re {first_mode_1based} (rad/m)")) {
        return (shift, ky * lambda);
    }
    // 3. Periodic BC `phase` on the physics problem
    if let Some(spec) = engine.project.spec(physics) {
        let facet_centroid = |labels: &[i32]| -> Option<[f64; 3]> {
            let mut sum = [0.0; 3];
            let mut cnt = 0usize;
            for f in 0..solver_mesh.n_facets() {
                if labels.contains(&solver_mesh.facet_labels[f]) {
                    for &v in solver_mesh.facet(f) {
                        let p = solver_mesh.point(v as usize);
                        for k in 0..solver_mesh.dim.min(3) {
                            sum[k] += p[k];
                        }
                        cnt += 1;
                    }
                }
            }
            (cnt > 0).then(|| [sum[0] / cnt as f64, sum[1] / cnt as f64, sum[2] / cnt as f64])
        };

        let mut best_bc_phase: Option<(f64, f64)> = None;
        for bc in spec.bcs.iter().filter(|b| b.kind == "periodic") {
            let raw_phase = bc.params.get("phase").and_then(Value::as_f64).unwrap_or(0.0);
            let dest_sel = bc.params.get("dest").and_then(Value::as_str);
            if let (Ok(src_l), Some(Ok(dst_l))) = (model.boundaries(&bc.on), dest_sel.map(|d| model.boundaries(d))) {
                if let (Some(cs), Some(cd)) = (facet_centroid(&src_l), facet_centroid(&dst_l)) {
                    let d = [cd[0] - cs[0], cd[1] - cs[1], cd[2] - cs[2]];
                    let dot = d[0] * shift[0] + d[1] * shift[1] + d[2] * shift[2];
                    let align = dot.abs();
                    let signed_ph = if dot < 0.0 { -raw_phase } else { raw_phase };
                    if best_bc_phase.is_none_or(|(best_align, _)| align > best_align) {
                        best_bc_phase = Some((align, signed_ph));
                    }
                    continue;
                }
            }
            if best_bc_phase.is_none() {
                best_bc_phase = Some((0.0, raw_phase));
            }
        }
        if let Some((_, ph)) = best_bc_phase {
            return (shift, ph);
        }
    }

    (shift, 0.0)
}

impl Engine {
    /// Draws the modes or 3D fields of `physics` and writes the figure to `path` (`.svg` or `.png`).
    pub fn figure(&self, args: &Value) -> Result<Value, Error> {
        let physics = args
            .get("physics")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::invalid("`physics` is required and must be text"))?;
        let path = args
            .get("path")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::invalid("`path` is required and must be text"))?;
        let ext = Path::new(path).extension().and_then(|e| e.to_str()).map(str::to_ascii_lowercase).unwrap_or_default();
        if ext != "svg" && ext != "png" {
            return Err(Error::invalid(format!("the file of a figure ends in .svg or .png, not '{path}'")));
        }
        let kind_name = args.get("kind").and_then(Value::as_str).unwrap_or("gallery");
        let kind = Kind::parse(kind_name).ok_or_else(|| {
            Error::invalid(format!(
                "there is no figure kind '{kind_name}' (gallery, components, slice3d, cutaway3d, isosurface3d, mesh3d)"
            ))
        })?;
        let mut o = FigureOptions { kind, ..Default::default() };
        if let Some(m) = numbers(args, "modes")? {
            // the modes are numbered from 1, as the study numbers them
            if m.iter().any(|x| *x < 1.0 || x.fract() != 0.0) {
                return Err(Error::invalid("`modes` are whole numbers from 1"));
            }
            o.modes = Some(m.iter().map(|x| *x as usize - 1).collect());
        }
        if let Some(n) = whole(args, "max_modes")? {
            o.max_modes = n;
        }
        if let Some(n) = whole(args, "columns")? {
            o.columns = n;
        }
        if let Some(c) = args.get("cmap").and_then(Value::as_str) {
            o.cmap = Cmap::parse(c).ok_or_else(|| {
                Error::invalid(format!("there is no colour map '{c}' (viridis, cividis, inferno, coolwarm)"))
            })?;
        }
        o.warp = real(args, "warp")?;
        if let Some(b) = args.get("arrows").and_then(Value::as_bool) {
            o.arrows = b;
        }
        if let Some(w) = numbers(args, "window")? {
            if w.len() != 4 {
                return Err(Error::invalid("`window` is [x0, x1, z0, z1] in metres"));
            }
            o.window = Some([w[0], w[1], w[2], w[3]]);
        }
        o.at = real(args, "at")?;
        if let Some(c) = numbers(args, "core")? {
            o.core = Some(c.iter().map(|x| *x as i32).collect());
        }
        o.family = args.get("family").and_then(Value::as_str).map(str::to_string);
        if let Some(w) = real(args, "panel_width")? {
            o.panel_width = w.clamp(60.0, 2000.0);
        }
        match args.get("expr") {
            None | Some(Value::Null) => {}
            Some(Value::String(s)) => o.expr = Some(s.clone()),
            Some(_) => return Err(Error::invalid("`expr` must be a formula string, e.g. \"re(Ex)\"")),
        }
        if let Some(iso) = numbers(args, "isosurfaces")? {
            if iso.is_empty() {
                return Err(Error::invalid("`isosurfaces` must not be empty"));
            }
            o.isosurfaces = Some(iso);
        }
        if let Some(rep) = whole(args, "repeat_periodic")? {
            if rep == 0 {
                return Err(Error::invalid("`repeat_periodic` must be at least 1"));
            }
            o.repeat_periodic = rep;
        }
        o.show_edges = match args.get("show_edges") {
            None | Some(Value::Null) => {
                if kind == Kind::Mesh3d {
                    EdgeMode::Mesh
                } else {
                    EdgeMode::Regions
                }
            }
            Some(Value::String(s)) => EdgeMode::parse(s)
                .ok_or_else(|| Error::invalid(format!("there is no `show_edges` mode '{s}' (none, regions, mesh)")))?,
            Some(Value::Bool(b)) => {
                if *b {
                    EdgeMode::Mesh
                } else {
                    EdgeMode::None
                }
            }
            Some(_) => return Err(Error::invalid("`show_edges` must be 'none', 'regions' or 'mesh'")),
        };
        o.camera = parse_camera(args)?;
        o.clip = parse_clip(args)?;
        let explicit_shift = match args.get("periodic_shift") {
            None | Some(Value::Null) => None,
            Some(v) => Some(vec3(v, "periodic_shift")?),
        };
        let explicit_phase = real(args, "periodic_phase")?;
        let scale = real(args, "scale")?.unwrap_or(2.0);

        let (model, result) = self.analysis_inputs(physics)?;
        let solver_mesh = model.solver_mesh();
        let (lo, hi) = solver_mesh.bbox();
        o.plane = parse_plane(args, o.at, lo, hi)?;

        if o.repeat_periodic > 1 {
            let first_mode_1based = o.modes.as_ref().and_then(|m| m.first().copied()).unwrap_or(0) + 1;
            let (shift, phase) = resolve_periodic_shift_and_phase(
                self,
                physics,
                &model,
                &solver_mesh,
                &result,
                first_mode_1based,
                explicit_shift,
                explicit_phase,
            );
            o.periodic_shift = Some(shift);
            o.periodic_phase = phase;
        } else {
            o.periodic_shift = explicit_shift;
            if let Some(ph) = explicit_phase {
                o.periodic_phase = ph;
            }
        }

        let fig = mode_figure_on(&solver_mesh, &result, &o)?;
        let (bytes, w, h): (Vec<u8>, f64, f64) = if ext == "svg" {
            (fig.canvas.to_svg().map_err(Error::invalid)?.into_bytes(), fig.canvas.width, fig.canvas.height)
        } else {
            (
                fig.canvas.to_png(scale).map_err(Error::invalid)?,
                (fig.canvas.width * scale).ceil(),
                (fig.canvas.height * scale).ceil(),
            )
        };
        // written under another name and renamed, so a figure that is half written is never under its own name
        let target = Path::new(path);
        let mut tmp = target.as_os_str().to_owned();
        tmp.push(".partial");
        let tmp = std::path::PathBuf::from(tmp);
        std::fs::write(&tmp, &bytes).map_err(|e| Error::invalid(format!("{path}: {e}")))?;
        std::fs::rename(&tmp, target).map_err(|e| {
            let _ = std::fs::remove_file(&tmp);
            Error::invalid(format!("{path}: {e}"))
        })?;
        let modes: Vec<Value> = fig
            .modes
            .iter()
            .map(|m| {
                json!({
                    "mode": m.mode,
                    "field": m.field,
                    "phase_deg": m.phase_deg,
                    "flipped": m.flipped,
                    "peak": m.peak,
                    "quadrature": m.quadrature,
                    "complex": m.complex
                })
            })
            .collect();
        Ok(json!({
            "path": path,
            "format": ext,
            "bytes": bytes.len(),
            "width": w,
            "height": h,
            "physics": physics,
            "kind": kind_name,
            "modes": modes,
            "window": fig.window,
            "section_y": fig.at,
            "stale": self.is_stale(physics),
            "age_s": self.result_age(physics),
        }))
    }
}
