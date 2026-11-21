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
    &["path", "layers", "z", "split_islands", "pin_layer", "simplify_tol", "cell", "layer_bbox"];

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
        let bbox = if layer_bbox { imported.layer_bbox } else { imported.cell_bbox };

        Ok(json!({
            "sheets": sheet_names,
            "layers": Value::Object(layers_json),
            "bbox": bbox,
            "pins": Value::Object(pins_json),
        }))
    }
}
