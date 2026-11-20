// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Absorbing layers across the two halves of a project: the geometry cuts a region into the physical part and one piece per side of the
//! window (`Pml` feature), the physics stretches the coordinates in each piece (`pml` condition). `add_pml_conditions` writes one condition
//! per piece from the geometry, so the two cannot disagree.

use crate::engine::Engine;
use nyx_physics::{BcSpec, Error};
use serde_json::{json, Value};

impl Engine {
    /// Adds a `pml` condition for every absorbing piece of the geometry to `physics` (replacing the conditions it had on those pieces).
    /// `strength` is the imaginary stretch at the outer edge; the order of the profile is the grading of the `Pml` feature.
    pub fn add_pml_conditions(&mut self, physics: &str, strength: f64) -> Result<Value, Error> {
        let p = crate::physics_for(&self.registry, physics)?;
        if !p.bc_kinds().iter().any(|k| k.kind == "pml") {
            return Err(Error::invalid(format!("{} has no absorbing layers (no `pml` condition)", p.name())));
        }
        let ev = femcore::meshing::geomdrive::evaluate(&self.project.geometry).map_err(Error::Invalid)?;
        let labels = ev.label_map().map_err(|e| Error::invalid(e.to_string()))?;
        let dim2 = self.project.geometry.dim == 2;
        let dir = |s: i8| match s {
            0 => "off",
            s if s < 0 => "minus",
            _ => "plus",
        };
        let pieces: Vec<(String, Value)> = ev
            .pml_regions()
            .iter()
            .map(|r| {
                let name = labels.domains[r.region as usize].name.clone();
                // a 2D section is the x-z plane: the layers are written along x and z (y is out of the plane)
                let y = if dim2 { "off" } else { dir(r.sides[1]) };
                (
                    name,
                    json!({ "x": dir(r.sides[0]), "y": y, "z": dir(r.sides[2]), "strength": strength, "order": r.grading, "re": 1.0 }),
                )
            })
            .collect();
        if pieces.is_empty() {
            return Err(Error::invalid("the geometry has no absorbing layers; add them with add_pml first"));
        }
        let spec = self.project.need_spec(physics)?;
        spec.bcs.retain(|b| !(b.kind == "pml" && pieces.iter().any(|(n, _)| *n == b.on)));
        for (name, params) in &pieces {
            spec.bcs.push(BcSpec { on: name.clone(), kind: "pml".into(), params: params.clone() });
        }
        Ok(json!({ "added": pieces.iter().map(|(n, _)| n.clone()).collect::<Vec<_>>() }))
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
