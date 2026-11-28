# 145. SIMD structured subpixel FDTD seeded by finite-element eigenmodes

Date: 2026-11-27

## Status

Accepted (FR-42 Track B). Extends ADR 1, ADR 3, ADR 9, ADR 11, ADR 12, ADR 13, ADR 19, ADR 25, ADR 45, ADR 82, ADR 106, ADR 131, ADR 133, and ADR 138.

## Context

Broadband photonic, plasmonic, and microwave scattering problems — octave-spanning transmission/reflection spectra, free-space metasurface diffraction, ultrashort pulse propagation, and multi-mode cavity ringdowns — complement frequency-domain finite-element sweeps (`Study::Driven`, ADR 106) and eigenmode expansion (`Study::Eme`, ADR 140) when solved directly in the time domain on a structured Cartesian Yee lattice:

1. **Structured SoA throughput and GPU-readiness:** Matrix-free leapfrog Yee updates over contiguous 64-byte-aligned Structure-of-Arrays (SoA) slabs `(Ex, Ey, Ez, Hx, Hy, Hz)` achieve high cache-line utilization and auto-vectorization on CPU SIMD lanes while providing a flat 1D buffer layout ready for future portable GPU compute backends without changing physics interfaces.
2. **Staircasing vs. curved CAD interfaces:** Naive cell-wise permittivity assignment on a rectilinear grid degrades convergence from $O(\Delta x^2)$ to $O(\Delta x)$ at slanted or curved dielectric interfaces and curved metallic walls. Sampling exact `geom` CSG/conic geometry (and curved `SimplexMesh` / `CurvedTet` fallback domains) across Yee edge dual rods and dual faces enables Farjadpour–Oskooi–Joannopoulos anisotropic inverse-permittivity tensor smoothing (`PolarizedAveraging`) and stabilized Dey–Mittra fractional edge/area PEC updates (`PECConformal`).
3. **Hybrid FEM-seeded modal ports:** Finite-difference cross-sectional mode solvers suffer from staircasing, split-step staggered eigenvector mismatch, and poor accuracy on curved or anisotropic waveguides. Because `nyx-ems` already provides high-order curved/anisotropic finite-element waveguide mode solvers (`modes` and `modes_complex`), solving the 2D port cross-section with `nyx-ems::modes` and projecting the continuous modal $(\mathbf{E}, \mathbf{H})$ profile onto staggered Yee tangential edges and faces eliminates modal injection artifacts and yields broadband S-parameters directly in `DrivenPortData` format.
4. **Broadband open boundaries, dispersion, NFFFT compatibility, and high-$Q$ ringdown:** Practical time-domain electromagnetics requires convolutional stretched-coordinate `CPML` (and adiabatic `Absorber` layers), multi-pole Drude/Lorentz auxiliary differential equations (ADEs) compatible with `femcore::material::Dispersion::Oscillators`, complex Bloch-periodic unit-cell wraps, dispersion-matched 1D auxiliary Total-Field/Scattered-Field (`TFSF`) propagators, on-the-fly discrete Fourier transform (`DFT`) monitors producing `Ex..Hz` `CellField`s consumable by `nyx-analysis::farfield` (`FarFieldAnalysis`), and harmonic-inversion (`ResonanceFinder`) to extract $(f_0, Q)$ from short ringdown traces without waiting for full $1/e$ decay of high-$Q$ cavities.

## Decision

* **Keep `"ems"` as the single electromagnetic physics crate (`ADR 12`, `ADR 25`):**
  * Add `Study::Fdtd(Box<FdtdStudyConfig>)` (`"fdtd"`, `Dims::Both`) to `nyx-ems`, implemented in `crates/nyx-ems/src/fdtd/` as a three-stage ADR 9 kernel (`assemble_fdtd`, `solve_fdtd`, `fdtd_fields`, plus convenience driver `run_fdtd`). Re-export `geom` from `femcore` (`femcore::geom`) so `nyx-ems` evaluates exact parametric CAD geometry while preserving strict workspace dependency layering (`crates/nyx-engine/tests/layering.rs`).
* **64-byte cache-aligned SoA slabs and wavelength-graded rectilinear `YeeGrid` (`fdtd::grid`):**
  * Store electric and magnetic components in 64-byte-aligned `AlignedSlab` buffers indexed in `(x, y, z)` contiguous `z`-pencils (`FdtdDim::TwoD` supporting `Pol2d::{TmEy, TeHy, Full}` in wafer `(x, z)` coordinates and `FdtdDim::ThreeD` in full `(x, y, z)` wafer coordinates).
  * Build smoothly graded rectilinear primal/dual coordinates (`AutoGridSpec` from `wavelength_min_m`, `min_steps_per_wave`, and `max_grading_ratio`, refined by `LayerRefinementSpec`, `GridOverrideBox`, and exact source/monitor `snap_planes`) directly from `geom` (`YeeGrid::from_geom`) or `Mesh` (`YeeGrid::from_mesh`).
* **Subpixel Farjadpour polarized dielectric averaging and Dey–Mittra conformal PEC (`fdtd::subpixel`):**
  * Query `FdtdGeometryOracle` (exact `geom::eval::Evaluated` CSG/curve containment and `closest_point` clearance, or BVH-accelerated straight/curved `SimplexMesh` point location) with $O(1)$ interior cell-clearance fast paths so subpixel quadrature runs only on true interface cells.
  * Assemble $\mathbf{W} = \mathbf{P}\langle\boldsymbol{\varepsilon}^{-1}\rangle\mathbf{P} + (\mathbf{I} - \mathbf{P})\langle\boldsymbol{\varepsilon}\rangle^{-1}(\mathbf{I} - \mathbf{P})$ along interface normals (`SubpixelMode::Polarized`, alongside `Volumetric` and `Staircased`), storing diagonal updates in aligned slabs and sparse off-diagonal interface/anisotropic couplings in `OffDiagEntry` lists.
  * Compute unblocked primal edge fractions $l_{e}$ and dual face open-area fractions $A_{f}$ (`conformal_pec` with sliver clamp `dm_min_edge_frac` and automatic CFL scaling) for curved PEC boundaries.
* **Drude/Lorentz ADEs, shared `[nonlinear]` Yee Ampere updates, Roden–Gedney `CPML`, `Absorber`, and Bloch periodicity (`fdtd::ade_pml`, `fdtd::stepper`):**
  * Advance second-order auxiliary polarization currents for arbitrary multi-pole Drude ($\omega_0 = 0$) and Lorentz ($\omega_0 > 0$) susceptibilities (`DispersiveDomainSpec`, automatically populated from library `[optical]` `Dispersion::Oscillators` when present or specified inline), and evaluate format-2 `[nonlinear]` (`ADR 143`) instantaneous $\chi^{(2)}$/QPM poling, $\chi^{(3)}$ Kerr, and semi-analytic TPA + free-carrier (`drude_carrier_step`) updates synchronously across non-PML Yee edges.
  * Provide six-face `FdtdBoundaries` supporting recursive-convolution stretched-coordinate `Cpml` ($\kappa, \sigma, \alpha$ polynomial grading), adiabatic conductivity `Absorber`, `Pec`, `Pmc`, midplane `Symmetry`, and real/complex Bloch `Periodic { bloch_k }`.
  * Partition outer `x`-slabs into fixed-size chunks via `femcore::par`, guaranteeing bitwise-identical (`0 ULP`) fields and monitor spectra across `NYX_THREADS=1` vs `8` (`ADR 11`).
* **FEM-seeded `ModeSource` / `ModeMonitor`, dispersion-matched `TFSF`, DFT monitors, and Prony `ResonanceFinder` (`fdtd::modes_bridge`, `fdtd::sources_monitors`):**
  * Extract a conforming 2D cross-sectional `SimplexMesh` on each `ModeSource` / `ModeMonitor` plane, solve guided eigenmodes via `nyx-ems::modes` / `modes_complex`, normalize discrete Yee cross-power $\frac{1}{2}\operatorname{Re}\sum (\mathbf{E}_t \times \mathbf{H}_t^*)\cdot\hat{\mathbf{n}}\,\Delta A = 1\,\mathrm{W}$, and inject/project unidirectional modal waves (`ModeMonitor` emitting `DrivenPortData` S/Z/Y matrices).
  * Drive `PlaneWave` and closed 6-face `Tfsf` boxes from a step-matched 1D auxiliary Yee propagator (`Aux1dYeeLine`) whose primal/dual steps along the propagation axis match the 3D grid identically, achieving $< 10^{-12}$ empty-box scattered-field leakage, alongside `GaussianBeam` and trilinearly weighted `PointDipole` sources.
  * Accumulate running $\exp(-j\omega t)\,\mathrm{d}t$ phasors (matching `nyx-ems` $\exp(+j\omega t)$ time-harmonic conventions) across `FluxMonitor`, Rayleigh-harmonic `DiffractionMonitor`, and volumetric `FieldMonitor` (projecting complex `Ex..Hz` `CellField`s onto `SimplexMesh` for direct downstream `nyx-analysis::farfield` Stratton–Chu NFFFT evaluation), plus `FieldTimeMonitor` coupled to a pure-Rust SVD Matrix-Pencil / Prony `find_resonances_prony` harmonic-inversion solver reporting $(f_0, Q)$.

## Consequences

* `physics_ems` advances additively (`StudyDesc` `"fdtd"`), `physics_ems_result` adds the `"fdtd"` result shape, and `kernel_stages` exports `assemble_fdtd`, `solve_fdtd`, `fdtd_fields`, and `run_fdtd`.
* Zero external C/Fortran or Python dependencies are introduced (`ADR 13`); all non-finite inputs and outputs are rejected (`ADR 131`), and VVUQ spatial/temporal discretization budgets (`ADR 133`) are populated on every `Study::Fdtd` result.

## Code

`crates/femcore/src/lib.rs`, `crates/nyx-ems/src/fdtd/{mod,grid,subpixel,ade_pml,modes_bridge,sources_monitors,stepper}.rs`, `crates/nyx-ems/src/{problem,sensitivity,lib}.rs`, `crates/nyx-ems/tests/{fr42_fdtd,contracts,finite_guard}.rs`, `crates/nyx-engine/tests/phase3_integration.rs`.
