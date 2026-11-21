# cQED result payloads: code

C4 level 4 (code) of the component **cQED result payloads** (`physics.cqed`) in the container **nyx-physics** (crate `nyx_physics`). [Component view of the container](../../3-components/physics.md) · [All code pages](../README.md)

> Typed superconducting quantum-circuit payloads (capacitance and inductance matrices, TLS geometry metadata, passive CircuitRom, driven port parameters).

**Decisions:** [ADR 103](../../../adr/0103-multi-rhs-capacitance-and-gauge-stabilized-magnetostatics.md), [ADR 104](../../../adr/0104-bcs-surface-impedance-and-lumped-josephson-ports.md), [ADR 105](../../../adr/0105-kicq-epr-craig-bampton-quantization-and-tls-participation.md), [ADR 106](../../../adr/0106-driven-frequency-domain-waveports-and-hermite-prom.md)  
**Size:** 1 file, 254 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-physics/src/cqed.rs` | 254 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Result and frame](./result.md) | result | Scalar, SolveResult |

## Used by

| Component | Modules | Items |
|---|---|---|
| [analysis: Capacitance and inductance matrices](../analysis/capacitance.md) | cqed | CapacitanceMatrixData, InductanceMatrixData |
| [analysis: cQED quantization and TLS participation](../analysis/cqed.md) | cqed | CircuitRom, ModeEnergyRow, RomPort, TlsMeta |
| [dfs: Electrostatics kernel](../dfs/electro.md) | cqed | CapacitanceMatrixData, TlsMeta |
| [dfs: Magnetostatics kernel](../dfs/magneto.md) | cqed | InductanceMatrixData |
| [ems: Driven frequency-domain and PROM](../ems/driven.md) | cqed | DrivenPortData |
| [ems: Thin-film IBC and Josephson ROM](../ems/ibc.md) | cqed | CircuitRom, ModeEnergyRow, RomPort, TlsMeta |
| [Physics trait and registry](./trait.md) | cqed | CAPACITANCE_KEY, CIRCUIT_ROM_KEY, CapacitanceMatrixData, CircuitRom, DRIVEN_PORT_KEY, DrivenPortData, INDUCTANCE_KEY, InductanceMatrixData +4 |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class CapacitanceMatrixData {
      <<struct>>
      +terminals Vec~String~
      +maxwell_matrix_f Vec
      +mutual_matrix_f Vec
      +ec_matrix_hz Vec
      +of(result) Option~Self~
      +scalar() tuple
  }
  class InductanceMatrixData {
      <<struct>>
      +terminals Vec~String~
      +inductance_matrix_h Vec
      +formulation String
      +of(result) Option~Self~
      +scalar() tuple
  }
  class TlsMeta {
      <<struct>>
      +domain_eps_r BTreeMap
      +domain_names BTreeMap
      +conductor_labels Vec~i32~
      +sheet_labels Vec~i32~
      +sheet_thickness_m BTreeMap
      +pec_labels Vec~i32~
      +of(result) Option~Self~
      +scalar() tuple
  }
  class RomPort {
      <<struct>>
      +name String
      +kind String
      +l_j_h Option~f64~
      +e_j_hz Option~f64~
      +c_j_f f64
      +e_l_hz Option~f64~
      +phi_ext_rad f64
      +n_g f64
      +more_3_fields
  }
  class ModeEnergyRow {
      <<struct>>
      +freq_hz f64
      +q_factor Option~f64~
      +u_elec_j f64
      +u_mag_vol_j f64
      +u_kin_ibc_j f64
      +u_ind_jj_j Vec~f64~
      +v_mj Vec
      +i_mj Vec
  }
  class CircuitRom {
      <<struct>>
      +ports Vec~RomPort~
      +r_static usize
      +r_dyn usize
      +k_r Vec
      +c_r Vec
      +m_r Vec
      +b_r Vec
      +c_port_f Vec
      +more_1_fields
      +of(result) Option~Self~
      +scalar() tuple
  }
  class DrivenPortData {
      <<struct>>
      +freqs_hz Vec~f64~
      +ports Vec~String~
      +z0_ohm Vec
      +s_matrix Vec
      +z_matrix Vec
      +y_matrix Vec
      +hdm_solves usize
      +rom_dim usize
      +more_2_fields
      +of(result) Option~Self~
      +scalar() tuple
  }
  class fns_cqed_rs {
      <<module cqed.rs>>
      +CAPACITANCE_KEY const
      +INDUCTANCE_KEY const
      +TLS_META_KEY const
      +CIRCUIT_ROM_KEY const
      +more_1_items
  }
  CircuitRom --> RomPort
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| const | `CAPACITANCE_KEY` | `cqed.rs:15` | `&str` —  Scalar key carrying `CapacitanceMatrixData` on a multi-terminal `dfs` electrostatic result. |
| const | `INDUCTANCE_KEY` | `cqed.rs:17` | `&str` —  Scalar key carrying `InductanceMatrixData` on a `dfs` magnetostatic result. |
| const | `TLS_META_KEY` | `cqed.rs:19` | `&str` —  Scalar key carrying `TlsMeta` on `dfs` or `ems` results for dielectric TLS surface participation. |
| const | `CIRCUIT_ROM_KEY` | `cqed.rs:21` | `&str` —  Scalar key carrying `CircuitRom` on `ems` `cavity` results with `lumped_jj` or `surface_impedance` sheets. |
| const | `DRIVEN_PORT_KEY` | `cqed.rs:23` | `&str` —  Scalar key carrying `DrivenPortData` on `ems` `driven` results. |
| struct | `CapacitanceMatrixData` | `cqed.rs:39` | fields: `terminals`, `maxwell_matrix_f`, `mutual_matrix_f`, `ec_matrix_hz` —  Multi-terminal Maxwell, SPICE/mutual capacitance matrices (F) and charging energy matrix `E_C / h` (Hz). |
| struct | `InductanceMatrixData` | `cqed.rs:58` | fields: `terminals`, `inductance_matrix_h`, `formulation` —  Multi-terminal inductance matrix `L` (H, or H/m for 2D cross-sections) from `dfs` magnetostatics. |
| struct | `TlsMeta` | `cqed.rs:76` | fields: `domain_eps_r`, `domain_names`, `conductor_labels`, `sheet_labels`, `sheet_thickness_m`, `pec_labels` —  Domain permittivities and metallized boundary/sheet annotations used by `tls_participation`. |
| struct | `RomPort` | `cqed.rs:108` | fields: `name`, `kind`, `l_j_h`, `e_j_hz`, `c_j_f`, `e_l_hz`, `phi_ext_rad`, `n_g`, `z0_ohm`, `direction`, `formulation` —  Metadata describing one lumped Josephson / nonlinear port inside `CircuitRom`. |
| struct | `ModeEnergyRow` | `cqed.rs:124` | fields: `freq_hz`, `q_factor`, `u_elec_j`, `u_mag_vol_j`, `u_kin_ibc_j`, `u_ind_jj_j`, `v_mj`, `i_mj` —  Exact energy partition and terminal phasors of one converged `ems` eigenmode. |
| struct | `CircuitRom` | `cqed.rs:139` | fields: `ports`, `r_static`, `r_dyn`, `k_r`, `c_r`, `m_r`, `b_r`, `c_port_f`, `modes` —  Multi-port Craig–Bampton passive reduced-order model `(K_r, C_r, M_r, B_r)` + exact modal energies. |
| struct | `DrivenPortData` | `cqed.rs:169` | fields: `freqs_hz`, `ports`, `z0_ohm`, `s_matrix`, `z_matrix`, `y_matrix`, `hdm_solves`, `rom_dim`, `anchor_freqs_hz`, `max_residual` —  Multi-frequency scattering, impedance, and admittance matrices from `ems` `driven` study. |
| fn | `CapacitanceMatrixData::of` | `cqed.rs:47` | `pub fn of(result: &SolveResult) -> Option<Self>` |
| fn | `CapacitanceMatrixData::scalar` | `cqed.rs:51` | `pub fn scalar(&self) -> (String, Scalar)` |
| fn | `InductanceMatrixData::of` | `cqed.rs:65` | `pub fn of(result: &SolveResult) -> Option<Self>` |
| fn | `InductanceMatrixData::scalar` | `cqed.rs:69` | `pub fn scalar(&self) -> (String, Scalar)` |
| fn | `TlsMeta::of` | `cqed.rs:97` | `pub fn of(result: &SolveResult) -> Option<Self>` |
| fn | `TlsMeta::scalar` | `cqed.rs:101` | `pub fn scalar(&self) -> (String, Scalar)` |
| fn | `CircuitRom::of` | `cqed.rs:158` | `pub fn of(result: &SolveResult) -> Option<Self>` |
| fn | `CircuitRom::scalar` | `cqed.rs:162` | `pub fn scalar(&self) -> (String, Scalar)` |
| fn | `DrivenPortData::of` | `cqed.rs:187` | `pub fn of(result: &SolveResult) -> Option<Self>` |
| fn | `DrivenPortData::scalar` | `cqed.rs:191` | `pub fn scalar(&self) -> (String, Scalar)` |
