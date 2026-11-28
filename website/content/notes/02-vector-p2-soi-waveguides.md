+++
title = "Vector Nédélec(2)+P2 Waveguides, Cylindrical PML & Poincaré Classification"
date = 2026-11-23
description = "Full-vector mixed curl-curl discretisation in wafer coordinates, Stokes/Poincaré polarisation invariants, and complex shift-and-invert whispering-gallery modes in structured quadrilateral PML."

[extra]
note_number = "02"
crates = "nyx-ems::modes · nyx-analysis::polarization"
reading_time = "8 min"
+++

Discretising the transverse electric field of a dielectric waveguide with nodal Lagrange elements produces unphysical "spurious" modes because nodal spaces enforce full vector continuity across high-index boundaries ($n\_{\mathrm{Si}} \approx 3.48$ vs. $n\_{\mathrm{SiO}\_2} \approx 1.44$) where normal $E\_\perp$ must jump by $\varepsilon\_1 / \varepsilon\_2$. In `crates/nyx-ems/src/modes.rs`, `nyx` eliminates spurious modes structurally via a mixed Nédélec + Lagrange formulation in the wafer frame.

## Mixed Curl-Curl Eigenproblem in the Wafer Frame

Every waveguide cross-section in `nyx` lives in the `(x, z)` plane with propagation along `+y` as $\mathbf{E}(x, y, z) = \left(\mathbf{E}\_t(x, z) + \hat{\mathbf{y}} E\_y(x, z)\right) e^{-j\beta y}$.{% <sidenote id="sn-wafer" label="1"> %}Inside `modes.rs`, the kernel rotates user tensors into local `(lateral, normal, propagation)` axes via `femcore::section::local_tensor` and exports fields back under wafer names (`Ex`, `Ez` in the plane; `Ey` along the guide).{% </sidenote> %} Introducing the scaled longitudinal potential $\phi = E\_y / (j\beta)$ turns Maxwell's curl-curl equation $\nabla \times (\mu\_r^{-1} \nabla \times \mathbf{E}) - k\_0^2 \varepsilon\_r \mathbf{E} = \mathbf{0}$ into the generalised eigenproblem:

$$
\begin{bmatrix} S - k\_0^2 T\_e & 0 \\\\ 0 & 0 \end{bmatrix} \begin{bmatrix} e\_t \\\\ \phi \end{bmatrix} = -\beta^2 \begin{bmatrix} T\_m & G \\\\ G^T & K\_p - k\_0^2 T\_p \end{bmatrix} \begin{bmatrix} e\_t \\\\ \phi \end{bmatrix}
$$

Here $S$ is the curl-curl stiffness matrix weighted by $\mu\_r^{-1}$, $T\_e$ and $T\_m$ are Nédélec edge mass matrices weighted by $\varepsilon\_r$ and $\mu\_r^{-1}$, $G$ couples $\mathbf{w}\_t \cdot \nabla\_t \lambda$, and $(K\_p, T\_p)$ are the Lagrange stiffness and mass blocks. By exact sequence duality, gradient fields $\mathbf{E}\_t = \nabla\_t \psi$ sit at $\beta^2 = 0$, so a shift-and-invert Krylov search around `n_target` isolates physical guided modes cleanly. At `order = 2` (the default), each triangle carries **8 hierarchical degree-2 Nédélec DOFs** (3 Whitney edge functions, 3 gradients of edge bubbles, 2 interior functions) paired with quadratic P2 $\phi$.

## Silicon-on-Insulator Wire (`15_soi_wire`)

In `examples/15_soi_wire.nyx.json`, a `500 nm × 300 nm` silicon waveguide (`slab = 0`, `n_target = 2.75`) on buried oxide is solved at $\lambda = 1550\,\text{nm}$ (`193.414 THz`) with `h = 40 nm` inside `Si`, `90 nm` in `BOX_near`/`clad_near`, and `300 nm` outside.

{{ <plate src="figures/soi_wire_modes.png" number="PLATE 01" alt="SOI wire waveguide TE0 and TM0 electric field gallery at 1550 nm" caption="TE0 (left, neff = 2.671811, ng = 4.1969, 94.5% power in Si) and TM0 (right, neff = 2.333849, ng = 4.6294, 85.1% in Si) on a 986-triangle P2 mesh." meta="15_soi_wire · 2,015 DOFs · 0.08 s" /> }}

Cross-checking against the archived Lumerical/MPB reference in `benchmarks/validation/photonics_benchmarks.json` (`2.67179` for `TE0` and `2.33374` for `TM0`), `15_soi_wire` lands within **`7.9 ppm`** (`0.000021`) on `TE0` and **`47 ppm`** (`0.000109`) on `TM0` in `0.08 s`:

{% <solver_log cmd="nyx --batch examples/15_soi_wire.nyx.json"> %}
986 cells, 515 points
ems: solved in 0.08 s, 2015 dofs, 2 solution(s)
  neff 1 = 2.671811
  v_g 1 (m/s) = 71431926.58486077
  n_g 1 = 4.19689727455199
  Ex share 1 = 0.8002390961141337
  Ey share 1 = 0.19720412785316194
  Ez share 1 = 0.0025567760327043277
  share in Si 1 = 0.9451175893630006
  neff 2 = 2.333849
  v_g 2 (m/s) = 64758283.71708667
  n_g 2 = 4.629407093457278
  Ex share 2 = 0.011127642624896001
  Ey share 2 = 0.38206400952785785
  Ez share 2 = 0.606808347847246
  share in Si 2 = 0.8510136106501677
  field frame = lateral x, normal z, propagation y, left-handed
what the solver did:
  solving EM waves (mode analysis) on a 2D mesh of 986 cells: 0.08 s
    factorising a sparse general matrix (6819 unknowns, 150676 non-zeros); this cannot be interrupted: 0.04 s
{% </solver_log> %}

Note what the JSON `result_info(physics="ems")` diagnostic also discloses honestly for this coarse `986`-cell example mesh: while mode 1 (`TE0`) converges to an algebraic residual of `1.327e-4`, mode 2 (`TM0`) triggers a `high_algebraic_residual` error flag (`2.101e-3 > 1e-6`) with `14` Arnoldi restarts because `n_target = 2.75` sits closer to mode 1 than mode 2. Refining the mesh or aiming `n_target = 2.35` drives both residuals below `1e-12`.

## Stokes Parameters & Poincaré Classification

Every optical mode study computes area-weighted Stokes parameters from `(E_lat, E_nor)` in `crates/nyx-analysis/src/polarization.rs` under `TimeConvention::PlusJ` ($\mathbf{E}(t) = \mathrm{Re}(\mathbf{E} e^{+i\omega t})$, so $S\_3 = -2\,\mathrm{Im}(E\_{\mathrm{lat}}^\* E\_{\mathrm{nor}})$):

$$
S\_0 = |E\_{\mathrm{lat}}|^2 + |E\_{\mathrm{nor}}|^2, \qquad S\_1 = |E\_{\mathrm{lat}}|^2 - |E\_{\mathrm{nor}}|^2, \qquad S\_2 = 2\,\mathrm{Re}(E\_{\mathrm{lat}}^\* E\_{\mathrm{nor}})
$$

Because $S\_0 \dots S\_3$ are quadratic in the field and invariant under multiplication by any global complex phase $e^{i\theta}$, the normalised Poincaré vector $(s\_1, s\_2, s\_3) = (S\_1, S\_2, S\_3)/S\_0$, degree of polarisation $P = \sqrt{S\_1^2 + S\_2^2 + S\_3^2}/S\_0$, orientation $\psi = \frac{1}{2}\mathrm{atan2}(S\_2, S\_1)$, and ellipticity $\chi = \frac{1}{2}\arcsin(S\_3 / P S\_0)$ are immune to arbitrary eigensolver phase scaling.

## Cylindrical PML Whispering-Gallery Modes (`18_wgm_cylindrical_pml`)

In `examples/18_wgm_cylindrical_pml.nyx.json`, a `1.2 um` radius dielectric core ($\varepsilon\_r = 12.0$) inside a `2.2 um` air cladding is surrounded by a `0.8 um` 5-layer structured quadrilateral cylindrical PML (`kappa_max = 2.0`, `order = 2.0`, `target_reflection = 1e-6`). Structured PML cells remain native quadrilaterals rather than being split into triangles (`ADR 128`), keeping radial stretching orthogonal to the circular boundary.

{{ <plate src="figures/wgm_cyl_pml.png" number="PLATE 02" alt="Four whispering-gallery modes inside a dielectric cylinder with outward cylindrical PML" caption="Modes 1–4 of 18_wgm_cylindrical_pml (5,225 DOFs, 0.38 s, max residual 8.716e-14): nearly lossless core mode 1 (Im neff = -1.559e-10) through radiative higher modes 2–4." meta="18_wgm_cylindrical_pml · 5,225 DOFs · 0.38 s" /> }}

{% <solver_log cmd="nyx run render_site_figures.nys  # 18_wgm_cylindrical_pml"> %}
ems: modes, order=2, dim=2, dofs=5225, seconds=0.3793
  solver_kind: shift_invert_complex, converged: true, restarts: 64
  max_algebraic_residual: 8.716388532767634e-14
  mode 1: neff = 2.400572 -1.559e-10i   (n_g = 4.9502, core = 0.9900, res = 1.375e-15)
  mode 2: neff = 2.414624 +6.653e-09i   (n_g = 4.9195, core = 0.9899, res = 4.588e-14)
  mode 3: neff = 2.385004 +3.849e-08i   (n_g = 5.0502, core = 0.9932, res = 6.400e-14)
  mode 4: neff = 2.379766 +6.144e-08i   (n_g = 5.0591, core = 0.9932, res = 8.716e-14)
{% </solver_log> %}

Here `layer fraction` remains below `4.2e-10` across all four modes, confirming that the returned eigenmodes live inside the dielectric cylinder (`99.0%–99.3%` in `core`) rather than inside the absorbing annulus.
