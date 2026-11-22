// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cross-crate superconducting-qubit (`cQED`) data contracts stored inside `SolveResult::scalars` (following `FieldFrame`, ADR 56).
//!
//! Physics crates (`nyx-dfs`, `nyx-ems`) never depend on `nyx-analysis` and `nyx-analysis` never depends on them (`layering.rs`).
//! Serialising these plain structs into `SolveResult::scalars` via `serde_json` (`float_roundtrip`) preserves that layering while
//! carrying exact multi-port capacitances, inductances, TLS interface metadata, Craig–Bampton reduced-order models, and driven
//! S/Z/Y parameters.

use crate::result::{Scalar, SolveResult};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Scalar key carrying `CapacitanceMatrixData` on a multi-terminal `dfs` electrostatic result.
pub const CAPACITANCE_KEY: &str = "capacitance matrix";
/// Scalar key carrying `InductanceMatrixData` on a `dfs` magnetostatic result.
pub const INDUCTANCE_KEY: &str = "inductance matrix";
/// Scalar key carrying `TlsMeta` on `dfs` or `ems` results for dielectric TLS surface participation.
pub const TLS_META_KEY: &str = "tls metadata";
/// Scalar key carrying `CircuitRom` on `ems` `cavity` results with `lumped_jj` or `surface_impedance` sheets.
pub const CIRCUIT_ROM_KEY: &str = "circuit rom";
/// Scalar key carrying `DrivenPortData` on `ems` `driven` results.
pub const DRIVEN_PORT_KEY: &str = "driven port matrix";

fn encode_json<T: Serialize>(key: &str, value: &T) -> (String, Scalar) {
    let text = serde_json::to_string(value).expect("cQED payload is plain serialisable data");
    (key.to_string(), Scalar::Text(text))
}

fn decode_json<T: for<'de> Deserialize<'de>>(key: &str, result: &SolveResult) -> Option<T> {
    result.scalars.iter().find(|(n, _)| n == key).and_then(|(_, v)| match v {
        Scalar::Text(t) => serde_json::from_str(t).ok(),
        Scalar::Num(_) => None,
    })
}

/// Multi-terminal Maxwell, SPICE/mutual capacitance matrices (F) and charging energy matrix `E_C / h` (Hz).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CapacitanceMatrixData {
    pub terminals: Vec<String>,
    pub maxwell_matrix_f: Vec<Vec<f64>>,
    pub mutual_matrix_f: Vec<Vec<f64>>,
    pub ec_matrix_hz: Vec<Vec<f64>>,
}

impl CapacitanceMatrixData {
    pub fn of(result: &SolveResult) -> Option<Self> {
        decode_json(CAPACITANCE_KEY, result)
    }

    pub fn scalar(&self) -> (String, Scalar) {
        encode_json(CAPACITANCE_KEY, self)
    }
}

/// Multi-terminal inductance matrix `L` (H, or H/m for 2D cross-sections) from `dfs` magnetostatics.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct InductanceMatrixData {
    pub terminals: Vec<String>,
    pub inductance_matrix_h: Vec<Vec<f64>>,
    pub formulation: String,
}

impl InductanceMatrixData {
    pub fn of(result: &SolveResult) -> Option<Self> {
        decode_json(INDUCTANCE_KEY, result)
    }

    pub fn scalar(&self) -> (String, Scalar) {
        encode_json(INDUCTANCE_KEY, self)
    }
}

/// Native surface-oxide TLS loss layer parameters (`thickness_m`, `eps_r`, `tan_delta`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TlsInterfaceSpec {
    pub thickness_m: f64,
    pub eps_r: f64,
    pub tan_delta: f64,
}

/// Material-resolved `[interfaces.{ma,ms,sa}]` TLS surface-loss parameters.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TlsInterfaces {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ma: Option<TlsInterfaceSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ms: Option<TlsInterfaceSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sa: Option<TlsInterfaceSpec>,
}

impl TlsInterfaces {
    pub fn is_empty(&self) -> bool {
        self.ma.is_none() && self.ms.is_none() && self.sa.is_none()
    }
}

/// Domain permittivities and metallized boundary/sheet annotations used by `tls_participation`.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TlsMeta {
    /// 1-based domain label -> isotropic relative permittivity `eps_r`.
    pub domain_eps_r: BTreeMap<i32, f64>,
    /// Optional 1-based domain label -> assigned region name from `ProblemSpec::materials`.
    #[serde(default)]
    pub domain_names: BTreeMap<i32, String>,
    /// 1-based boundary labels that act as metallic surfaces (`ground`, `potential`, `terminal`, `pec`, `surface_impedance`, etc.).
    #[serde(default)]
    pub conductor_labels: Vec<i32>,
    /// 1-based boundary labels that correspond to interior zero-thickness 2D sheets.
    #[serde(default)]
    pub sheet_labels: Vec<i32>,
    /// 1-based boundary label -> physical metal film thickness (m), if supplied (e.g. from `surface_impedance`).
    #[serde(default)]
    pub sheet_thickness_m: BTreeMap<i32, f64>,
    /// 1-based boundary label -> exterior wall flag (if known).
    #[serde(default)]
    pub pec_labels: Vec<i32>,
    /// Optional `[interfaces.{ma,ms,sa}]` parameters resolved from the problem's material library.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interfaces: Option<TlsInterfaces>,
}

impl TlsMeta {
    pub fn of(result: &SolveResult) -> Option<Self> {
        decode_json(TLS_META_KEY, result)
    }

    pub fn scalar(&self) -> (String, Scalar) {
        encode_json(TLS_META_KEY, self)
    }
}

/// Metadata describing one lumped Josephson / nonlinear port inside `CircuitRom`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RomPort {
    pub name: String,
    pub kind: String,
    pub l_j_h: Option<f64>,
    pub e_j_hz: Option<f64>,
    pub c_j_f: f64,
    pub e_l_hz: Option<f64>,
    pub phi_ext_rad: f64,
    pub n_g: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub squid_asymmetry: Option<f64>,
    pub z0_ohm: Option<f64>,
    pub direction: [f64; 3],
    pub formulation: String,
}

/// Exact energy partition and terminal phasors of one converged `ems` eigenmode.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModeEnergyRow {
    pub freq_hz: f64,
    pub q_factor: Option<f64>,
    pub u_elec_j: f64,
    pub u_mag_vol_j: f64,
    pub u_kin_ibc_j: f64,
    pub u_ind_jj_j: Vec<f64>,
    /// Complex port voltage phasor `(re, im)` in V at `U_elec = 1 J` normalization (or scaled consistently with `u_elec_j`).
    pub v_mj: Vec<[f64; 2]>,
    /// Complex port inductive current phasor `(re, im)` in A.
    pub i_mj: Vec<[f64; 2]>,
}

/// Multi-port Craig–Bampton passive reduced-order model `(K_r, C_r, M_r, B_r)` + exact modal energies.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CircuitRom {
    pub ports: Vec<RomPort>,
    pub r_static: usize,
    pub r_dyn: usize,
    /// Reduced passive stiffness matrix `K_r = Q^T (K + K_IBC) Q` (`r x r`).
    pub k_r: Vec<Vec<f64>>,
    /// Reduced loss matrix `C_r = Q^T C_IBC Q` (`r x r`, empty if lossless).
    #[serde(default)]
    pub c_r: Vec<Vec<f64>>,
    /// Reduced passive mass matrix `M_r = Q^T M Q` (`r x r`).
    pub m_r: Vec<Vec<f64>>,
    /// Reduced terminal incidence matrix `B_r = Q^T B` (`r x J`).
    pub b_r: Vec<Vec<f64>>,
    /// Passive 3D terminal capacitance matrix `C_port` (F, `J x J`, excluding parallel `C_j`).
    pub c_port_f: Vec<Vec<f64>>,
    pub modes: Vec<ModeEnergyRow>,
}

impl CircuitRom {
    pub fn of(result: &SolveResult) -> Option<Self> {
        decode_json(CIRCUIT_ROM_KEY, result)
    }

    pub fn scalar(&self) -> (String, Scalar) {
        encode_json(CIRCUIT_ROM_KEY, self)
    }
}

/// Multi-frequency scattering, impedance, and admittance matrices from `ems` `driven` study.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DrivenPortData {
    pub freqs_hz: Vec<f64>,
    pub ports: Vec<String>,
    /// `[re, im]` reference port impedance at each frequency or per port (`n_freqs x P`, or `1 x P`).
    pub z0_ohm: Vec<Vec<[f64; 2]>>,
    /// `n_freqs x P x P` complex `[re, im]`.
    pub s_matrix: Vec<Vec<Vec<[f64; 2]>>>,
    /// `n_freqs x P x P` complex `[re, im]` in Ohms.
    pub z_matrix: Vec<Vec<Vec<[f64; 2]>>>,
    /// `n_freqs x P x P` complex `[re, im]` in Siemens.
    pub y_matrix: Vec<Vec<Vec<[f64; 2]>>>,
    pub hdm_solves: usize,
    pub rom_dim: usize,
    pub anchor_freqs_hz: Vec<f64>,
    pub max_residual: f64,
}

impl DrivenPortData {
    pub fn of(result: &SolveResult) -> Option<Self> {
        decode_json(DRIVEN_PORT_KEY, result)
    }

    pub fn scalar(&self) -> (String, Scalar) {
        encode_json(DRIVEN_PORT_KEY, self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cqed_payloads_round_trip_through_solve_result_scalars() {
        let cap = CapacitanceMatrixData {
            terminals: vec!["Q1_pad1".into(), "Q1_pad2".into()],
            maxwell_matrix_f: vec![vec![60e-15, -20e-15], vec![-20e-15, 65e-15]],
            mutual_matrix_f: vec![vec![40e-15, 20e-15], vec![20e-15, 45e-15]],
            ec_matrix_hz: vec![vec![3.5e8, 1.1e8], vec![1.1e8, 3.2e8]],
        };
        let tls = TlsMeta {
            domain_eps_r: [(1, 1.0), (2, 11.45)].into(),
            domain_names: [(1, "vacuum".into()), (2, "silicon".into())].into(),
            conductor_labels: vec![3, 4],
            sheet_labels: vec![4],
            sheet_thickness_m: [(4, 100e-9)].into(),
            pec_labels: vec![3],
            interfaces: Some(TlsInterfaces {
                ma: Some(TlsInterfaceSpec { thickness_m: 3.0e-9, eps_r: 9.8, tan_delta: 1.5e-3 }),
                ms: Some(TlsInterfaceSpec { thickness_m: 2.0e-9, eps_r: 9.8, tan_delta: 1.5e-3 }),
                sa: Some(TlsInterfaceSpec { thickness_m: 2.0e-9, eps_r: 3.9, tan_delta: 1.5e-3 }),
            }),
        };
        let rom = CircuitRom {
            ports: vec![RomPort {
                name: "Q1_JJ".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(8.15e-9),
                e_j_hz: Some(20.06e9),
                c_j_f: 2.0e-15,
                e_l_hz: Some(1.0e9),
                phi_ext_rad: std::f64::consts::PI,
                n_g: 0.0,
                squid_asymmetry: Some(0.1),
                z0_ohm: None,
                direction: [1.0, 0.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 1,
            r_dyn: 1,
            k_r: vec![vec![1.0, 0.0], vec![0.0, 2.0]],
            c_r: vec![],
            m_r: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            b_r: vec![vec![0.5], vec![0.1]],
            c_port_f: vec![vec![19.37e-15]],
            modes: vec![ModeEnergyRow {
                freq_hz: 6.5e9,
                q_factor: None,
                u_elec_j: 1.0,
                u_mag_vol_j: 0.8,
                u_kin_ibc_j: 0.05,
                u_ind_jj_j: vec![0.15],
                v_mj: vec![[1.2, 0.0]],
                i_mj: vec![[0.0, -3.5e-3]],
            }],
        };
        let res = SolveResult { scalars: vec![cap.scalar(), tls.scalar(), rom.scalar()], ..Default::default() };
        assert_eq!(CapacitanceMatrixData::of(&res), Some(cap));
        assert_eq!(TlsMeta::of(&res), Some(tls));
        assert_eq!(CircuitRom::of(&res), Some(rom));
    }
}
