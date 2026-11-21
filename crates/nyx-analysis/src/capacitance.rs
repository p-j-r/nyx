// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-terminal capacitance (`capacitance_matrix`) and inductance (`inductance_matrix`) analyses of `dfs` results.
//!
//! `CapacitanceMatrixAnalysis` reads `CapacitanceMatrixData` from `SolveResult::scalars` (produced by an electrostatic
//! study of `dfs` with one or more `terminal` boundary conditions) and returns the Maxwell capacitance matrix (F), the
//! mutual/SPICE capacitance matrix (F), and the charging energy matrix `E_C / h = (e^2 / 2h) C_M^(-1)` (Hz).
//!
//! `InductanceMatrixAnalysis` reads `InductanceMatrixData` from `SolveResult::scalars` (produced by a magnetostatic
//! study of `dfs`) and returns the multi-terminal inductance matrix `L` (H, or H/m for 2D `cross_section` models) and
//! the active magnetostatic formulation (`"cross_section"`, `"planar_stream"`, or `"mixed_3d"`).

use crate::Analysis;
use nyx_physics::{CapacitanceMatrixData, Error, InductanceMatrixData, Mesh, SolveResult};
use serde_json::{json, Value};

/// Extracts the typed [`CapacitanceMatrixData`] payload from `result`.
pub fn capacitance_matrix(result: &SolveResult) -> Result<CapacitanceMatrixData, Error> {
    CapacitanceMatrixData::of(result).ok_or_else(|| {
        Error::invalid(
            "the result carries no capacitance matrix (run an electrostatic study of dfs with terminal conditions)",
        )
    })
}

/// True when `result` carries a multi-terminal capacitance matrix payload.
pub fn capacitance_matrix_applies(result: &SolveResult) -> bool {
    CapacitanceMatrixData::of(result).is_some()
}

/// Multi-terminal Maxwell, mutual, and charging-energy capacitance matrix analysis behind the [`Analysis`] trait.
pub struct CapacitanceMatrixAnalysis;

impl Analysis for CapacitanceMatrixAnalysis {
    fn id(&self) -> &'static str {
        "capacitance_matrix"
    }

    fn name(&self) -> &'static str {
        "Multi-terminal capacitance and charging energy matrix"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        capacitance_matrix_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        let cap = capacitance_matrix(result)?;
        Ok(json!({
            "terminals": cap.terminals,
            "maxwell_matrix_F": cap.maxwell_matrix_f,
            "mutual_matrix_F": cap.mutual_matrix_f,
            "ec_matrix_Hz": cap.ec_matrix_hz,
        }))
    }
}

/// Extracts the typed [`InductanceMatrixData`] payload from `result`.
pub fn inductance_matrix(result: &SolveResult) -> Result<InductanceMatrixData, Error> {
    InductanceMatrixData::of(result)
        .ok_or_else(|| Error::invalid("the result carries no inductance matrix (run a magnetostatic study of dfs)"))
}

/// True when `result` carries a magnetostatic inductance matrix payload.
pub fn inductance_matrix_applies(result: &SolveResult) -> bool {
    InductanceMatrixData::of(result).is_some()
}

/// Multi-terminal magnetostatic inductance matrix analysis behind the [`Analysis`] trait.
pub struct InductanceMatrixAnalysis;

impl Analysis for InductanceMatrixAnalysis {
    fn id(&self) -> &'static str {
        "inductance_matrix"
    }

    fn name(&self) -> &'static str {
        "Multi-terminal inductance matrix"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        inductance_matrix_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        let ind = inductance_matrix(result)?;
        Ok(json!({
            "terminals": ind.terminals,
            "inductance_matrix_H": ind.inductance_matrix_h,
            "formulation": ind.formulation,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capacitance_and_inductance_analyses_round_trip_to_json() {
        let sm = femcore::testkit::rect(1.0, 0.5, 2, 2);
        let mesh = Mesh::from_solver_mesh(&sm).unwrap();

        let cap = CapacitanceMatrixData {
            terminals: vec!["pad1".into(), "pad2".into()],
            maxwell_matrix_f: vec![vec![80e-15, -30e-15], vec![-30e-15, 85e-15]],
            mutual_matrix_f: vec![vec![50e-15, 30e-15], vec![30e-15, 55e-15]],
            ec_matrix_hz: vec![vec![2.8e8, 1.0e8], vec![1.0e8, 2.6e8]],
        };
        let ind = InductanceMatrixData {
            terminals: vec!["loop1".into()],
            inductance_matrix_h: vec![vec![150e-12]],
            formulation: "mixed_3d".into(),
        };
        let res = SolveResult {
            physics: "dfs".into(),
            study: "electrostatic".into(),
            scalars: vec![cap.scalar(), ind.scalar()],
            ..Default::default()
        };

        assert!(CapacitanceMatrixAnalysis.applies_to(&res));
        assert!(InductanceMatrixAnalysis.applies_to(&res));

        let cap_json = CapacitanceMatrixAnalysis.run(&mesh, &res).unwrap();
        assert_eq!(cap_json["terminals"], json!(["pad1", "pad2"]));
        assert_eq!(cap_json["maxwell_matrix_F"][0][0].as_f64().unwrap(), 80e-15);
        assert_eq!(cap_json["mutual_matrix_F"][0][1].as_f64().unwrap(), 30e-15);
        assert_eq!(cap_json["ec_matrix_Hz"][1][1].as_f64().unwrap(), 2.6e8);

        let ind_json = InductanceMatrixAnalysis.run(&mesh, &res).unwrap();
        assert_eq!(ind_json["terminals"], json!(["loop1"]));
        assert_eq!(ind_json["inductance_matrix_H"][0][0].as_f64().unwrap(), 150e-12);
        assert_eq!(ind_json["formulation"], json!("mixed_3d"));

        let empty = SolveResult::default();
        assert!(!CapacitanceMatrixAnalysis.applies_to(&empty));
        assert!(!InductanceMatrixAnalysis.applies_to(&empty));
        assert!(CapacitanceMatrixAnalysis.run(&mesh, &empty).is_err());
        assert!(InductanceMatrixAnalysis.run(&mesh, &empty).is_err());
    }
}
