// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The `figure` command (ADR 99): a figure of the modes of the latest result of a physics, written as SVG or PNG. The window, a script and the MCP
//! server all call it; the drawing is `nyx_analysis::figure`, which names no physics.

use crate::engine::Engine;
use nyx_analysis::figure::{mode_figure, Cmap, FigureOptions, Kind};
use nyx_physics::Error;
use serde_json::{json, Value};
use std::path::Path;

fn numbers(a: &Value, k: &str) -> Result<Option<Vec<f64>>, Error> {
    match a.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(v)) => v
            .iter()
            .map(|x| x.as_f64().ok_or_else(|| Error::invalid(format!("`{k}` must be a list of numbers"))))
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
        Some(v) => v.as_f64().map(Some).ok_or_else(|| Error::invalid(format!("`{k}` must be a number"))),
    }
}

impl Engine {
    /// Draws the modes of `physics` and writes the figure to `path` (`.svg` or `.png`).
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
        let kind = match args.get("kind").and_then(Value::as_str).unwrap_or("gallery") {
            "gallery" => Kind::Gallery,
            "components" => Kind::Components,
            other => return Err(Error::invalid(format!("there is no figure kind '{other}' (gallery, components)"))),
        };
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
        let scale = real(args, "scale")?.unwrap_or(2.0);

        let (model, result) = self.analysis_inputs(physics)?;
        let fig = mode_figure(model.mesh(), &result, &o)?;
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
            .map(|m| json!({ "mode": m.mode, "field": m.field, "phase_deg": m.phase_deg, "flipped": m.flipped, "peak": m.peak, "quadrature": m.quadrature, "complex": m.complex }))
            .collect();
        Ok(json!({
            "path": path, "format": ext, "bytes": bytes.len(), "width": w, "height": h, "physics": physics,
            "kind": args.get("kind").and_then(Value::as_str).unwrap_or("gallery"),
            "modes": modes, "window": fig.window, "section_y": fig.at,
            "stale": self.is_stale(physics), "age_s": self.result_age(physics),
        }))
    }
}
