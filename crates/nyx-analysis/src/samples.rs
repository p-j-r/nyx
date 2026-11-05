//! The modes of a `SolveResult` as weighted samples of the field, which is what the polarization analysis integrates.
//!
//! A mode study of `ems` states its frame (`nyx_physics::FieldFrame`: the lateral, normal and propagation axes, the handedness, the absorbing
//! domains to leave out) and gives, per mode, the complex electric field at the nodes of every cell as the cell fields `Ex`, `Ey`, `Ez`
//! (wafer axes). This module reads them back:
//!
//! * a **cross-section** (2D) is sampled at the 24 points of `femcore::optical::fine_rule` of every cell outside the absorbing domains, with
//!   the field interpolated from the nodes (P1 or P2 shape functions); the weight is the rule weight times the cell area. The mode fields are
//!   piecewise quadratic, so the interpolation is exact and the Stokes sums are those of the solver's own field;
//! * a **cavity** (3D) is sampled at the four vertices of every cell, each standing for a quarter of the cell volume.

use crate::modes::{AnalysisOptions, ModeInput};
use crate::polarization::{Cx, Sample, TimeConvention};
use femcore::fem::{geometry, shape};
use femcore::fields::local_bary;
use femcore::optical::fine_rule;
use nyx_physics::{discretise, CellField, Error, FieldFrame, Mesh, Solution, SolveResult};

/// The modes of a result, ready for `analyze_modes`, and the options that go with them.
pub struct PolarizationInput {
    pub modes: Vec<ModeInput>,
    pub options: AnalysisOptions,
}

/// Whether a result has what the polarization analysis needs: a stated frame and, in every solution, the three complex field components as
/// cell fields.
pub fn applies(result: &SolveResult) -> bool {
    let Some(frame) = FieldFrame::of(result) else {
        return false;
    };
    !result.solutions.is_empty()
        && result
            .solutions
            .iter()
            .all(|s| frame.names().iter().all(|a| s.cell_fields.iter().any(|f| f.name == format!("E{a}"))))
}

/// Reads the modes of `result` (solved on `mesh`) as weighted samples.
pub fn polarization_input(mesh: &Mesh, result: &SolveResult) -> Result<PolarizationInput, Error> {
    let frame = FieldFrame::of(result).ok_or_else(|| {
        Error::invalid("the result states no field frame (scalar \"field frame\"): it has no modes to analyse")
    })?;
    let disc = discretise(mesh)?;
    let dim = disc.dim;
    if dim != result.dim as usize {
        return Err(Error::invalid(format!("the result is {}D but the mesh is {dim}D", result.dim)));
    }
    let geo = geometry(&disc);
    let skipped: Vec<bool> = (0..disc.n_cells()).map(|c| frame.outside.contains(&disc.cell_labels[c])).collect();
    let cavity = dim == 3;
    let mut modes = vec![];
    for (j, sol) in result.solutions.iter().enumerate() {
        let comp = |axis: usize| -> Result<&CellField, Error> {
            let name = format!("E{}", ["x", "y", "z"][axis]);
            sol.cell_fields.iter().find(|f| f.name == name).ok_or_else(|| {
                let have: Vec<&str> = sol.cell_fields.iter().map(|f| f.name.as_str()).collect();
                Error::invalid(format!(
                    "mode {} ('{}') has no cell field '{name}'; it has: {}",
                    j + 1,
                    sol.label,
                    have.join(", ")
                ))
            })
        };
        let (lat, nor, prop) = (comp(frame.lateral)?, comp(frame.normal)?, comp(frame.propagation)?);
        let nloc = lat.nloc;
        if nor.nloc != nloc || prop.nloc != nloc || nloc == 0 {
            return Err(Error::invalid(format!(
                "the field components of mode {} have different numbers of nodes per cell",
                j + 1
            )));
        }
        let ncell = disc.n_cells();
        for f in [lat, nor, prop] {
            if f.re.len() != ncell * nloc || (!f.im.is_empty() && f.im.len() != ncell * nloc) {
                return Err(Error::invalid(format!(
                    "cell field '{}' of mode {} has {} values for {ncell} cells of {nloc} nodes (is the result from this mesh?)",
                    f.name,
                    j + 1,
                    f.re.len()
                )));
            }
        }
        let order = match (dim, nloc) {
            (2, 3) | (3, 4) => 1,
            (2, 6) | (3, 10) => 2,
            _ => return Err(Error::invalid(format!("{nloc} nodes per cell is not a P1 or P2 cell of a {dim}D mesh"))),
        };
        let at = |f: &CellField, base: usize, n: &[f64]| -> Cx {
            let mut v = [0.0, 0.0];
            for (a, na) in n.iter().enumerate() {
                v[0] += na * f.re[base + a];
                if !f.im.is_empty() {
                    v[1] += na * f.im[base + a];
                }
            }
            v
        };
        let mut samples = vec![];
        if cavity {
            for c in (0..ncell).filter(|&c| !skipped[c]) {
                for a in 0..4 {
                    let i = c * nloc + a;
                    let one = |f: &CellField| -> Cx { [f.re[i], if f.im.is_empty() { 0.0 } else { f.im[i] }] };
                    samples.push(Sample { lat: one(lat), nor: one(nor), prop: one(prop), w: geo.vol[c] / 4.0 });
                }
            }
        } else {
            let (pts, wts) = fine_rule();
            let shapes: Vec<Vec<f64>> = pts.iter().map(|p| shape(p, order).0).collect();
            debug_assert_eq!(shapes[0].len(), local_bary(2, order).len());
            for c in (0..ncell).filter(|&c| !skipped[c]) {
                let base = c * nloc;
                for (n, w) in shapes.iter().zip(&wts) {
                    samples.push(Sample {
                        lat: at(lat, base, n),
                        nor: at(nor, base, n),
                        prop: at(prop, base, n),
                        w: w * geo.vol[c],
                    });
                }
            }
        }
        modes.push(ModeInput { label: format!("mode {}", j + 1), neff: eigenvalue_of(sol, cavity)?, samples });
    }
    let mut options = AnalysisOptions::default();
    if frame.left_handed {
        options.class.time = TimeConvention::MinusI;
    }
    if cavity {
        let fmax = modes.iter().map(|m| m.neff.0.abs()).fold(0.0, f64::max);
        options.value_label = "f (THz)".into();
        options.degeneracy_tol = 1e-6 * fmax;
    }
    Ok(PolarizationInput { modes, options })
}

/// The number shown beside a mode: the effective index of a waveguide mode, the frequency in THz of a cavity mode.
fn eigenvalue_of(sol: &Solution, cavity: bool) -> Result<(f64, f64), Error> {
    if cavity {
        let [re, im] = sol
            .eigenvalue
            .or_else(|| sol.value.map(|v| [v, 0.0]))
            .ok_or_else(|| Error::invalid("a cavity mode has no frequency"))?;
        Ok((re / 1e12, im / 1e12))
    } else {
        let [re, im] = sol.eigenvalue.ok_or_else(|| Error::invalid("a mode has no effective index"))?;
        Ok((re, im))
    }
}
