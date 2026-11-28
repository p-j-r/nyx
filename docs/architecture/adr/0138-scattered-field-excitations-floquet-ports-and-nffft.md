# ADR 139: Scattered-field excitations, Floquet diffraction ports, and Stratton-Chu NFFFT / beam overlap

- Status: Accepted (FR-36)
- Date: 2026-11-27
- Scope: `nyx-ems::{excitation, floquet_port, driven, problem, lib}`, `nyx-analysis::{farfield, beam, lib}`, `nyx-engine::analyse`

## Context

The driven EM study (`ems` `driven`, ADR 106/127) only knew circuit-style excitations: lumped sheet ports and 2D boundary-mode
waveports. Free-space photonics (scatterers, metasurfaces, gratings, focusing optics, fiber coupling) needs three more things:

1. incident fields that are not guided modes of a port (plane waves, focused beams, point dipoles), solved without having to
   resolve the source in the mesh;
2. exact radiation conditions and per-order bookkeeping on the faces of periodic unit cells (diffraction orders);
3. post-solve quantities defined at infinity or on remote planes (radiation patterns, cross-sections, focal spots, coupling
   efficiency into a beam or a fiber mode).

## Decision

### 1. Scattered-field formulation (`nyx-ems::excitation`, `driven`)

The total field is split as `E = E_inc + E_scat`, where `E_inc` solves Maxwell's equations in a homogeneous background
`(eps_bg, mu_bg)`. The driven operator `A(omega)` is unchanged; the right-hand side is

    <w, b_scat> = ∫_{Ω_scat} [ k0² w·(eps_r − eps_bg) E_inc − (curl w)·(mu_r⁻¹ − mu_bg⁻¹)(curl E_inc) ] dV

over the contrast cells (domains named by the condition's `on`, non-PML, with non-zero contrast). It vanishes identically in the
background and in every PML cell, so PMLs and diffraction ports only ever see outgoing scattered waves. The integrand is evaluated
with the same quadrature as the system (degree-4 Keast on straight cells, the ADR 82 curved rule on `CurvedTet` cells) and is
accumulated in fixed 1024-cell chunks merged in chunk order (bitwise identical for any `NYX_THREADS`).

Excitations are analytical (`Excitation`, `exp(+iωt)`, SI, wafer frame), each with closed-form `E` and `curl E`:

| kind | model | exactness |
|---|---|---|
| `plane_wave` | `A e^{iφ} p e^{−i k k̂·(r−r0)}` | exact |
| `gaussian_beam` | regular complex-source/sink beam `curl curl (p j0(kR))/k²`, `R² = (r − r_f + i b k̂)²`, `b = k w0²/2`, evaluated with `e^{−kb}`-scaled spherical Bessel ratios (no overflow, no branch cut) | exact Maxwell solution; `paraxial` option is approximate |
| `astigmatic_gaussian_beam` | elliptical paraxial Gaussian with longitudinal correction, `curl` by 2nd-order central differences | approximate (paraxial) |
| `thin_lens_beam` | vectorial angular spectrum of an aplanatic lens of numerical aperture NA (Gauss-Legendre × trapezoid plane-wave superposition, optional Gaussian pupil filling) | exact for the discrete spectrum |
| `dipole` | Hertzian electric / magnetic point dipole radiating in the background; must lie outside the contrast cells (checked) | exact |

Beams are normalized so that `p̂*·E(focus) = A e^{iφ}`.

Several excitations are solved as extra right-hand-side columns on the same factorization as the port columns (multi-RHS on
`DirectSolver` / `BlrLu` / AMS through the shared `femcore::linalg::prom::SparseDrivenSolver`). The scattered-field load is not an
affine (`γ(ω) b`) function of frequency, so it cannot enter the Hermite/Galerkin PROM image space of ADR 129: when `prom_tol` is
set, ports keep the PROM sweep and the excitation columns are solved at full order on one factorization per frequency (counted in
`hdm_solves`). Stage 3 reconstructs total fields `E = E_inc + E_scat`, `H = i curl E/(ω μ0 μ_r)` (cell fields `Ex..Hz`) and, with
`keep_scattered`, the scattered `Esx..Esz`, `Hsx..Hsz` (`H_s = H − H_inc`). Powers per excitation and frequency:

* `P_abs = −(ω ε0/2) ∫ Im(E*·ε_r E) dV + (1/(2 ω μ0)) ∫ Im(curl E*·μ_r⁻¹ curl E) dV` (non-PML cells),
* `P_ext = −(ω ε0/2) Im ∫_{Ω_scat} E_inc*·Δε E dV − (ω μ0/2) Im ∫ Δμ H_inc*·H dV` (work of the incident field on the induced
  currents), `P_sca = P_ext − P_abs` (optical theorem),
* `P_sca_pml = Im(x_s^H (K_pml − k0² M_pml) x_s)/(2 ω μ0)`, an independent check through the PML,

and, for plane waves, `σ = P/I_inc`. They are multiplied by the midplane-symmetry factor (`symmetry`, now also accepted by the
driven study: 2 per mirror plane). Results: JSON scalar `scattering` (all frequencies) plus numeric scalars
`<name>: sigma_ext (m^2)` etc. for single-frequency runs.

### 2. Floquet diffraction ports (`nyx-ems::floquet_port`, `driven`)

`diffraction_port` (boundary) on a planar, axis-normal, full unit-cell face of a transversely Floquet-periodic cell expands the
tangential field in orthonormal Rayleigh harmonics `e_mn^σ = û_mn^σ e^{−i k_t,mn·(r − r0)}/√A`,
`k_t,mn = k_t0 + 2π(m/Λ1, n/Λ2)`, `û^TM = k_t/|k_t|`, `û^TE = â_n × û^TM` (shared global axis so both ports of a cell use the same
basis). `k_t0` follows from the periodic pairs (`k_t0·s = φ`, matching `x_dest = e^{−iφ} x_src`). For an outgoing harmonic the
exact DtN map is `n × μ_r⁻¹ curl E = γ E_t`, `γ^TE = iβ/μ_r`, `γ^TM = i k²/(β μ_r)`, `β = √(k² − |k_t|²)`, `Im β ≤ 0`. Each
harmonic adds the rank-one affine term `γ_h(ω) ĉ_h t̂_hᵀ` (`ĉ = D Pᴴ c`, `t̂ = D Pᵀ c̄ = conj(ĉ)`), with analytic `γ'`, `γ''`
for the Hermite PROM. Because the term is not complex-symmetric, `PromOperator::is_symmetric` is false and the receive functional
`port_b_recv_hat` differs from the drive (both already supported by `femcore::linalg::prom`).

Every harmonic that propagates at `f_max` becomes an S-matrix channel `name:(m,n)TE|TM` of the unified port list (unit-power
waves: excitation `2γ√(2/Re Y)` with `Y = γ/(iωμ0)`, response `S_pq = t̂_pᵀ ŷ_q √(Re Y_p/2) − δ_pq`), so reflection and transmission
efficiencies are `R_mn = |S_{(m,n),in}|²`, `T_mn = |S_{(m,n)',in}|²` inside the existing `DrivenPortData` (no schema change);
rows/columns of channels below cutoff at a given frequency are zero. Evanescent harmonics (orders up to `m_max`, `n_max`,
automatic = all propagating at `f_max` plus one ring) only enter the DtN. Channel metadata are in the JSON scalar
`diffraction_channels`. The DtN block is dense on the face DOFs (`n_face²` entries per harmonic): adequate for unit cells; an
augmented low-rank formulation is the documented upgrade path.

### 3. Far field and overlap (`nyx-analysis::farfield`)

A pure post-solve analysis on `SolveResult` + `Mesh` (no dependency on `nyx-ems`, ADR 25 layering): Stratton-Chu equivalent
currents `J_s = n × H`, `M_s = −n × E` on a closed facet set, box or planar aperture; far field on angular grids (radiation
intensity, directivity, gain, bistatic RCS), k-space grids, and near-zone-valid focal planes (PSF, Strehl against the Airy
reference); vectorial four-field overlap `overlap_beam` against Gaussian, thin-lens and LP fiber modes. Both are registered as
engine analyses (`farfield`, `overlap_beam`).

## Consequences

* The driven study accepts problems without lumped/wave ports (excitation-only or diffraction-port-only runs); the representative
  port solution is emitted only when ports exist, followed by one solution per excitation (all frequencies when
  `n_freq × n_exc ≤ 8`, else midband).
* New BC kinds (`diffraction_port`, `plane_wave`, `gaussian_beam`, `astigmatic_gaussian_beam`, `thin_lens_beam`, `dipole`) and the
  driven `symmetry` parameter are additive descriptor changes (contract exemplars and `VERSIONS.md` updated with them).
* Mixed prism/hex meshes reject diffraction ports and excitations with a clear error (pure tetrahedral meshes, straight or
  curved, are supported).
* Port faces must be straight (curved adjacent cells are refused, as for waveports) and the port medium homogeneous, isotropic and
  lossless.
* The optional stratified (layered-background) incident field is not part of this decision; the contrast is always taken against
  one homogeneous background.

## Verification

`crates/nyx-ems/tests/fr36_scattering.rs` and the unit tests of `excitation.rs`: analytic `curl` vs finite differences for every
excitation, complex-source beam normalization / Rayleigh-range law / no overflow; periodic slab Fresnel `R`, `T` (TE, TM, normal
and 30° incidence) through diffraction ports against the Airy formulas, power conservation and reciprocity; lamellar grating with
±1 orders (power conservation over all Rayleigh channels, mirrored-order reciprocity); scattered-field plane wave on the periodic
slab (scattered amplitudes vs Airy, optical theorem); Lorenz-Mie extinction/scattering of dielectric and lossy spheres on curved
tetrahedra with a spherical PML; bitwise determinism `NYX_THREADS=1` vs `8`. The far-field analyses carry their own tests in
`nyx-analysis` (dipole directivity 1.5, aperture PSF vs Airy, Gaussian self-overlap and waist-mismatch closed form).
