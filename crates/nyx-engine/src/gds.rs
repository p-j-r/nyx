// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! GDSII layout import on `Engine`: reads a binary `.gds` file via `geom::gds`, adds merged/split
//! `Polygons` and conductive `Sheet` features at plane `z`, and returns sheet names, per-layer
//! mappings, chip bounding box, and matched `TEXT` pins.

use crate::engine::Engine;
use crate::geometry;
use geom::gds::{import_gds_bytes, GdsImportOptions};
use nyx_physics::Error;
use serde_json::{json, Map, Value};
use std::path::PathBuf;

const ALLOWED_FIELDS: &[&str] =
    &["path", "layers", "z", "split_islands", "pin_layer", "simplify_tol", "cell", "layer_bbox", "lateral_pad"];

const CHIP_STACK_ALLOWED_FIELDS: &[&str] = &[
    "bbox",
    "lateral_pad",
    "z",
    "t_sub",
    "substrate_thickness",
    "h_vac_top",
    "t_vac_top",
    "vac_top_thickness",
    "h_vac_bot",
    "t_vac_bot",
    "vac_bot_thickness",
    "vac_bot_box",
    "substrate_box",
    "vac_top_box",
    "vac_bot_region",
    "substrate_region",
    "vac_top_region",
    "overlay",
];

fn parse_layer_pair(v: &Value, what: &str) -> Result<(i16, i16), Error> {
    let arr = v
        .as_array()
        .filter(|a| a.len() == 2)
        .ok_or_else(|| Error::invalid(format!("import_gds: `{what}` must be a 2-integer pair `[layer, datatype]`")))?;
    let parse_i16 = |x: &Value, label: &str| -> Result<i16, Error> {
        let n = x
            .as_i64()
            .ok_or_else(|| Error::invalid(format!("import_gds: `{what}` {label} must be an integer, got {x}")))?;
        i16::try_from(n).map_err(|_| {
            Error::invalid(format!(
                "import_gds: `{what}` {label} {n} is outside valid GDSII i16 range (-32768..=32767)"
            ))
        })
    };
    Ok((parse_i16(&arr[0], "layer")?, parse_i16(&arr[1], "datatype/texttype")?))
}

fn free_feature_name(list: &geom::model::FeatureList, base: &str) -> String {
    if list.feature_by_name(base).is_none() {
        return base.to_string();
    }
    for k in 1usize.. {
        let cand = format!("{base}_{k}");
        if list.feature_by_name(&cand).is_none() {
            return cand;
        }
    }
    unreachable!()
}

fn point_in_ring_even_odd(ring: &[[f64; 2]], p: [f64; 2]) -> bool {
    let n = ring.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0] {
            inside = !inside;
        }
    }
    inside
}

impl Engine {
    /// Imports GDSII layers as 2D conductive `Sheet` features at height `z` in a 3D project.
    pub fn import_gds(&mut self, args: &Value) -> Result<Value, Error> {
        let obj = args
            .as_object()
            .ok_or_else(|| Error::invalid(format!("import_gds: arguments must be a JSON object, got {args}")))?;

        for k in obj.keys() {
            if !ALLOWED_FIELDS.contains(&k.as_str()) {
                return Err(Error::invalid(format!(
                    "import_gds: unknown field '{k}' (expected: {})",
                    ALLOWED_FIELDS.join(", ")
                )));
            }
        }

        if self.project.geometry.dim != 3 {
            return Err(Error::invalid(format!(
                "import_gds requires a 3D project (`project_new` with `dim: 3`), but current project has dim = {}",
                self.project.geometry.dim
            )));
        }

        let raw_path = obj
            .get("path")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| Error::invalid("import_gds: `path` is required and must be a non-empty file path"))?;

        let resolved_path: PathBuf = match &self.path_root {
            Some(root) => crate::path_guard::confine(root, raw_path)
                .map_err(|why| Error::invalid(format!("import_gds: `path`: {why}")))?,
            None => PathBuf::from(raw_path),
        };

        let layers_obj = obj.get("layers").and_then(Value::as_object).filter(|m| !m.is_empty()).ok_or_else(|| {
            Error::invalid(
                "import_gds: `layers` is required and must be a non-empty map `{ \"<name>\": [layer, datatype], ... }`",
            )
        })?;

        let mut layers = Vec::with_capacity(layers_obj.len());
        for (name, spec) in layers_obj {
            geom::model::check_name("layer", name).map_err(|e| Error::invalid(format!("import_gds: {e}")))?;
            let pair = parse_layer_pair(spec, &format!("layers.{name}"))?;
            layers.push((name.clone(), pair));
        }

        let z = match obj.get("z") {
            None | Some(Value::Null) => 0.0,
            Some(v) => v
                .as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| Error::invalid("import_gds: `z` must be a finite number in metres"))?,
        };

        let split_islands = match obj.get("split_islands") {
            None | Some(Value::Null) => true,
            Some(v) => v.as_bool().ok_or_else(|| Error::invalid("import_gds: `split_islands` must be a boolean"))?,
        };

        let pin_layer = match obj.get("pin_layer") {
            None | Some(Value::Null) => None,
            Some(v) => Some(parse_layer_pair(v, "pin_layer")?),
        };

        let simplify_tol = match obj.get("simplify_tol") {
            None | Some(Value::Null) => 0.0,
            Some(v) => v
                .as_f64()
                .filter(|&x| x.is_finite() && x >= 0.0)
                .ok_or_else(|| Error::invalid("import_gds: `simplify_tol` must be a non-negative finite number"))?,
        };

        let cell = match obj.get("cell") {
            None | Some(Value::Null) => None,
            Some(v) => {
                let s = v.as_str().ok_or_else(|| Error::invalid("import_gds: `cell` must be a string"))?;
                if s.trim().is_empty() {
                    return Err(Error::invalid("import_gds: `cell` must not be empty"));
                }
                Some(s.to_string())
            }
        };

        let layer_bbox = match obj.get("layer_bbox") {
            None | Some(Value::Null) => true,
            Some(v) => v.as_bool().ok_or_else(|| Error::invalid("import_gds: `layer_bbox` must be a boolean"))?,
        };

        let lateral_pad = match obj.get("lateral_pad") {
            None | Some(Value::Null) => 0.0,
            Some(v) => v
                .as_f64()
                .filter(|&x| x.is_finite() && x >= 0.0)
                .ok_or_else(|| Error::invalid("import_gds: `lateral_pad` must be a non-negative finite number"))?,
        };

        let bytes = std::fs::read(&resolved_path)
            .map_err(|e| Error::invalid(format!("import_gds: failed to read '{}': {e}", resolved_path.display())))?;

        let gds_opts = GdsImportOptions { layers, cell, split_islands, pin_layer, simplify_tol };
        let imported = import_gds_bytes(&bytes, &gds_opts).map_err(|e| Error::invalid(format!("import_gds: {e}")))?;

        let mut copy = self.project.geometry.clone();
        for sh in &imported.sheets {
            if copy.feature_by_name(&sh.name).is_some() {
                return Err(Error::invalid(format!(
                    "import_gds: a feature named '{}' already exists in the project",
                    sh.name
                )));
            }

            // Detect whether a merged (`split_islands == false`) sheet contains an island nested inside
            // another island's hole on the same layer; if so, union per-island Polygons features cleanly.
            let nested_in_hole = sh.islands.len() > 1
                && (0..sh.islands.len()).any(|i| {
                    let pt = sh.islands[i].0[0];
                    (0..sh.islands.len())
                        .any(|j| i != j && sh.islands[j].1.iter().any(|h| point_in_ring_even_odd(h, pt)))
                });

            let poly_input = if nested_in_hole {
                let mut part_names = Vec::with_capacity(sh.islands.len());
                for (idx, (outer, holes)) in sh.islands.iter().enumerate() {
                    let pname = free_feature_name(&copy, &format!("poly_{}_{idx}", sh.name));
                    geometry::add_feature(
                        &mut copy,
                        "Polygons",
                        Some(&pname),
                        &json!({
                            "polygons": [outer],
                            "holes": holes,
                            "simplify_tol": simplify_tol,
                            "min_hole_area": 0.0,
                        }),
                    )?;
                    part_names.push(pname);
                }
                let uname = free_feature_name(&copy, &format!("poly_{}", sh.name));
                geometry::add_feature(&mut copy, "Union", Some(&uname), &json!({ "ids": part_names }))?;
                uname
            } else {
                let pname = free_feature_name(&copy, &format!("poly_{}", sh.name));
                geometry::add_feature(
                    &mut copy,
                    "Polygons",
                    Some(&pname),
                    &json!({
                        "polygons": sh.outers,
                        "holes": sh.holes,
                        "simplify_tol": simplify_tol,
                        "min_hole_area": 0.0,
                    }),
                )?;
                pname
            };

            geometry::add_feature(&mut copy, "Sheet", Some(&sh.name), &json!({ "input": poly_input, "z": z }))?;
        }

        copy.validate().map_err(|e| Error::invalid(format!("import_gds: {e}")))?;
        self.project.geometry = copy;

        let sheet_names: Vec<String> = imported.sheets.iter().map(|s| s.name.clone()).collect();
        let mut layers_json = Map::new();
        for (k, v) in imported.layers {
            layers_json.insert(k, json!(v));
        }
        let mut pins_json = Map::new();
        for (pin, sh_name) in imported.pins {
            pins_json.insert(pin, json!(sh_name));
        }
        let unpadded_bbox = if layer_bbox { imported.layer_bbox } else { imported.cell_bbox };
        let bbox = if lateral_pad > 0.0 {
            [
                unpadded_bbox[0] - lateral_pad,
                unpadded_bbox[1] - lateral_pad,
                unpadded_bbox[2] + lateral_pad,
                unpadded_bbox[3] + lateral_pad,
            ]
        } else {
            unpadded_bbox
        };

        Ok(json!({
            "sheets": sheet_names,
            "layers": Value::Object(layers_json),
            "bbox": bbox,
            "unpadded_bbox": unpadded_bbox,
            "padded_bbox": bbox,
            "lateral_pad": lateral_pad,
            "pins": Value::Object(pins_json),
        }))
    }

    /// Pure `geom` feature-list helper that constructs a 3D chip dielectric stack `(VacBot, Si, VacTop)`
    /// over a 2D `[xmin, ymin, xmax, ymax]` bounding box (explicitly supplied or inferred from existing
    /// 2D/sheet features in `self.project.geometry`), expanded laterally by `lateral_pad >= 0.0`.
    pub fn add_chip_stack(&mut self, args: &Value) -> Result<Value, Error> {
        let obj = args
            .as_object()
            .ok_or_else(|| Error::invalid(format!("add_chip_stack: arguments must be a JSON object, got {args}")))?;

        for k in obj.keys() {
            if !CHIP_STACK_ALLOWED_FIELDS.contains(&k.as_str()) {
                return Err(Error::invalid(format!(
                    "add_chip_stack: unknown field '{k}' (expected: {})",
                    CHIP_STACK_ALLOWED_FIELDS.join(", ")
                )));
            }
        }

        if self.project.geometry.dim != 3 {
            return Err(Error::invalid(format!(
                "add_chip_stack requires a 3D project (`project_new` with `dim: 3`), but current project has dim = {}",
                self.project.geometry.dim
            )));
        }

        let lateral_pad = match obj.get("lateral_pad") {
            None | Some(Value::Null) => 0.0,
            Some(v) => v
                .as_f64()
                .filter(|&x| x.is_finite() && x >= 0.0)
                .ok_or_else(|| Error::invalid("add_chip_stack: `lateral_pad` must be a non-negative finite number"))?,
        };

        let z = match obj.get("z") {
            None | Some(Value::Null) => 0.0,
            Some(v) => v
                .as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| Error::invalid("add_chip_stack: `z` must be a finite number in metres"))?,
        };

        let t_sub_val = obj.get("t_sub").or_else(|| obj.get("substrate_thickness"));
        let t_sub = match t_sub_val {
            None | Some(Value::Null) => 350.0e-6,
            Some(v) => v
                .as_f64()
                .filter(|&x| x.is_finite() && x > 0.0)
                .ok_or_else(|| Error::invalid("add_chip_stack: `t_sub` must be a positive finite number in metres"))?,
        };

        let h_vac_top_val =
            obj.get("h_vac_top").or_else(|| obj.get("t_vac_top")).or_else(|| obj.get("vac_top_thickness"));
        let h_vac_top = match h_vac_top_val {
            None | Some(Value::Null) => 3.0 * t_sub,
            Some(v) => v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid("add_chip_stack: `h_vac_top` must be a positive finite number in metres")
            })?,
        };

        let h_vac_bot_val =
            obj.get("h_vac_bot").or_else(|| obj.get("t_vac_bot")).or_else(|| obj.get("vac_bot_thickness"));
        let h_vac_bot = match h_vac_bot_val {
            None | Some(Value::Null) => h_vac_top,
            Some(v) => v.as_f64().filter(|&x| x.is_finite() && x >= 0.0).ok_or_else(|| {
                Error::invalid("add_chip_stack: `h_vac_bot` must be a non-negative finite number in metres")
            })?,
        };

        let unpadded_bbox: [f64; 4] = match obj.get("bbox") {
            Some(v) if !v.is_null() => {
                let arr = v
                    .as_array()
                    .filter(|a| a.len() == 4)
                    .ok_or_else(|| Error::invalid("add_chip_stack: `bbox` must be `[xmin, ymin, xmax, ymax]`"))?;
                let mut b = [0.0; 4];
                for i in 0..4 {
                    b[i] = arr[i]
                        .as_f64()
                        .filter(|x| x.is_finite())
                        .ok_or_else(|| Error::invalid("add_chip_stack: `bbox` entries must be finite numbers"))?;
                }
                if b[0] >= b[2] || b[1] >= b[3] {
                    return Err(Error::invalid(format!(
                        "add_chip_stack: `bbox` must satisfy xmin < xmax and ymin < ymax, got {b:?}"
                    )));
                }
                b
            }
            _ => infer_xy_bbox(&self.project.geometry)?,
        };

        let padded_bbox = if lateral_pad > 0.0 {
            [
                unpadded_bbox[0] - lateral_pad,
                unpadded_bbox[1] - lateral_pad,
                unpadded_bbox[2] + lateral_pad,
                unpadded_bbox[3] + lateral_pad,
            ]
        } else {
            unpadded_bbox
        };

        let str_opt = |key: &str, default: &str| -> Result<String, Error> {
            match obj.get(key) {
                None | Some(Value::Null) => Ok(default.to_string()),
                Some(v) => {
                    let s = v
                        .as_str()
                        .filter(|s| !s.trim().is_empty())
                        .ok_or_else(|| Error::invalid(format!("add_chip_stack: `{key}` must be a non-empty string")))?;
                    geom::model::check_name(key, s).map_err(|e| Error::invalid(format!("add_chip_stack: {e}")))?;
                    Ok(s.to_string())
                }
            }
        };

        let vac_bot_box = str_opt("vac_bot_box", "vac_bot")?;
        let substrate_box = str_opt("substrate_box", "si_sub")?;
        let vac_top_box = str_opt("vac_top_box", "vac_top")?;
        let vac_bot_region = str_opt("vac_bot_region", "VacBot")?;
        let substrate_region = str_opt("substrate_region", "Si")?;
        let vac_top_region = str_opt("vac_top_region", "VacTop")?;

        let add_overlay = match obj.get("overlay") {
            None | Some(Value::Null) => true,
            Some(v) => v.as_bool().ok_or_else(|| Error::invalid("add_chip_stack: `overlay` must be a boolean"))?,
        };

        let [xmin, ymin, xmax, ymax] = padded_bbox;
        let z_sub_top = z;
        let z_sub_bot = z - t_sub;
        let z_vac_top = z + h_vac_top;
        let z_vac_bot = z_sub_bot - h_vac_bot;

        let mut copy = self.project.geometry.clone();
        let mut boxes = Vec::new();
        let mut regions = Vec::new();

        if h_vac_bot > 0.0 {
            geometry::add_feature(
                &mut copy,
                "Box",
                Some(&vac_bot_box),
                &json!({
                    "min": [xmin, ymin, z_vac_bot],
                    "max": [xmax, ymax, z_sub_bot],
                }),
            )?;
            boxes.push(vac_bot_box.clone());
            regions.push(json!({ "name": vac_bot_region, "input": vac_bot_box }));
        }

        geometry::add_feature(
            &mut copy,
            "Box",
            Some(&substrate_box),
            &json!({
                "min": [xmin, ymin, z_sub_bot],
                "max": [xmax, ymax, z_sub_top],
            }),
        )?;
        boxes.push(substrate_box.clone());
        regions.push(json!({ "name": substrate_region, "input": substrate_box }));

        geometry::add_feature(
            &mut copy,
            "Box",
            Some(&vac_top_box),
            &json!({
                "min": [xmin, ymin, z_sub_top],
                "max": [xmax, ymax, z_vac_top],
            }),
        )?;
        boxes.push(vac_top_box.clone());
        regions.push(json!({ "name": vac_top_region, "input": vac_top_box }));

        if add_overlay {
            if let Some(existing_ov) =
                copy.features.iter().find(|f| matches!(f.op, geom::model::Op::Overlay { .. })).map(|f| f.name.clone())
            {
                geometry::update_feature(&mut copy, &existing_ov, &json!({ "regions": regions }))?;
            } else {
                geometry::add_feature(&mut copy, "Overlay", None, &json!({ "regions": regions }))?;
            }
        }

        copy.validate().map_err(|e| Error::invalid(format!("add_chip_stack: {e}")))?;
        self.project.geometry = copy;

        let region_names: Vec<String> = if h_vac_bot > 0.0 {
            vec![vac_bot_region, substrate_region, vac_top_region]
        } else {
            vec![substrate_region, vac_top_region]
        };

        Ok(json!({
            "bbox": padded_bbox,
            "unpadded_bbox": unpadded_bbox,
            "padded_bbox": padded_bbox,
            "lateral_pad": lateral_pad,
            "z_range": [z_vac_bot, z_vac_top],
            "boxes": boxes,
            "regions": region_names,
        }))
    }
}

fn infer_xy_bbox(list: &geom::model::FeatureList) -> Result<[f64; 4], Error> {
    use geom::model::Op;
    let mut env = geom::expr::Env::default();
    for p in &list.params {
        let v = env.eval(&p.value, &p.name).map_err(|e| Error::invalid(format!("add_chip_stack: {e}")))?;
        env.vars.push((p.name.clone(), v));
    }
    let val = |e: &geom::expr::Expr| -> Result<f64, Error> {
        env.eval(e, "coordinate").map_err(|err| Error::invalid(format!("add_chip_stack: {err}")))
    };
    let (mut xmin, mut ymin) = (f64::INFINITY, f64::INFINITY);
    let (mut xmax, mut ymax) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut include = |x: f64, y: f64| {
        if x.is_finite() && y.is_finite() {
            xmin = xmin.min(x);
            xmax = xmax.max(x);
            ymin = ymin.min(y);
            ymax = ymax.max(y);
        }
    };
    for f in &list.features {
        if f.suppressed {
            continue;
        }
        match &f.op {
            Op::Rect { min, max } | Op::RoundedRect { min, max, .. } => {
                include(val(&min[0])?, val(&min[1])?);
                include(val(&max[0])?, val(&max[1])?);
            }
            Op::Circle { center, r } => {
                let (cx, cy, rv) = (val(&center[0])?, val(&center[1])?, val(r)?.abs());
                include(cx - rv, cy - rv);
                include(cx + rv, cy + rv);
            }
            Op::Ellipse { center, radii, .. } => {
                let (cx, cy) = (val(&center[0])?, val(&center[1])?);
                let rm = val(&radii[0])?.abs().max(val(&radii[1])?.abs());
                include(cx - rm, cy - rm);
                include(cx + rm, cy + rm);
            }
            Op::Polygon { pts, .. } => {
                for p in pts {
                    include(val(&p[0])?, val(&p[1])?);
                }
            }
            Op::Polygons { polygons, .. } => {
                for ring in polygons {
                    for p in ring {
                        include(val(&p[0])?, val(&p[1])?);
                    }
                }
            }
            Op::Box { min, max } => {
                include(val(&min[0])?, val(&min[1])?);
                include(val(&max[0])?, val(&max[1])?);
            }
            _ => {}
        }
    }
    if !xmin.is_finite() || !ymin.is_finite() || xmin >= xmax || ymin >= ymax {
        return Err(Error::invalid(
            "add_chip_stack: `bbox` was not provided and no 2D/sheet feature bounding box could be inferred from the project",
        ));
    }
    Ok([xmin, ymin, xmax, ymax])
}
