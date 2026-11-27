# Photonic use cases: 300 nm SOI wires, Bragg gratings, and the polarization analysis

Everything here is in `femcore::usecases::{soi, bragg}`, in `nyx-analysis` (`polarization`, `modes`; they were `femcore::{analysis, polarization}` until phase 5) and in the application
(Examples 16 to 19, the **Analysis** tab). All numbers below were produced by this code on the machine it was written on; the tests named in each
section re-check the ones that can be checked against a formula.

## The polarization analysis (Stokes parameters, Poincare sphere)

The analysis is a function of a result (phase 5, ADR 56): `nyx_analysis::polarization(mesh, result)`. The mode studies of `ems` put the complex field of every mode in their cell fields and state their frame in the scalar `field frame`, and the analysis reads that: the 2D mode solver (real and complex, with PML or a tilted crystal), the 3D
eigenmode solver (cavities and Floquet cells), the Brillouin study of a straight guide (its optical pump and Stokes modes) and the Brillouin study of a ring
(whispering-gallery modes). It is computed from the **complex fields the solver found** at quadrature points, not from the plotted real fields (those lose
the phase between Ex and Ey, so they cannot give s2 and s3 for hybrid or leaky modes).

* Per mode: integrated Stokes vector, point `(s1, s2, s3)` on the unit sphere, degree of polarization `P`, ellipticity and orientation angles, the share of the
  energy in Ex (`TE`) and in the longitudinal component (`Ez`), and a class (quasi-TE, quasi-TM, linear +-45, elliptical, circular, hybrid).
* Axes: x lateral (quasi-TE), y the surface normal (quasi-TM), z the propagation axis. Wafer-frame projects keep this convention for the analysis (it
  is the local frame of the mode solvers, `section::local_tensor`). For 3D problems the propagation axis is the axis a Floquet pair is periodic along (else y in the wafer frame, z
  otherwise); for rings (E_r, E_z, E_theta) are (lateral, normal, propagation). The sense of rotation is `s3 > 0` counter-clockwise for a wave
  `exp(-i k z)`; left-handed axis triples (rings, wafer-frame 3D) are handled by flipping the sign convention.
* **Degenerate modes** (effective indices within 5e-4, or frequencies within 1e-6 relative) form a group: the eigensolver returns an arbitrary mixture of
  them, so their own Stokes points mean nothing. The analysis reports the group Stokes vector, which is independent of the mixture, and the *resolved* pair
  (the two modes of definite x/y polarization). Example: the square 300 nm wire in oxide has an exactly degenerate TE/TM pair; the solver returned two
  mixed modes with `s = (+0.58, -0.66, 0)` and the resolved pair is quasi-TE and quasi-TM.
* In the application: the **Analysis** tab appears next to Results whenever the result carries an analysis. It draws the sphere with one numbered marker per mode
  (the line from the centre to the marker ends on the surface for a fully polarized mode, inside it for a hybrid one; dashed for degenerate modes), lists the
  modes (click one to show its field, for electromagnetic mode studies), the degenerate groups, the table, and exports CSV and SVG. In a batch run the table is
  part of the printed log, and `nyx --script "example 15; compute; tab analysis; shot FILE; quit"` takes a picture of the tab (scripts count the examples from 0, the menu from 1).

Tests: `polarization` and `analysis` unit tests (from the session that wrote them), `tests/usecases.rs`, the PML case in `tests/photonics.rs`, the cavity case in
`tests/emw3d.rs`, and in `tests/brillouin.rs` a ring of 20 mm radius against the straight guide (same TE share and longitudinal share to 0.01).

## Use case 1: modes of a wire on 300 nm SOI (example 16, `femcore::usecases::soi`, `cargo run --release -p femcore --example soi_modes`)

Silicon (n = 3.4757) of 300 nm thickness and 500 nm width on 2 um of oxide (n = 1.4440) under an oxide cladding, 1550 nm, perfectly conducting walls 1 um
from the core, shape mesh with 40 nm elements in the silicon (1312 cells, 9077 unknowns), about 0.1 s for two modes. This is an assumption about what "300 nm SOI" means (a 300 nm
device layer, strip wire); width, thickness, rib slab, side-wall angles (left and right separately), cladding and mesh sizes are fields of `SoiWire`.

| mode | neff | TE share | P | s1 | class |
|---|---|---|---|---|---|
| 1 | 2.67179 | 0.989 | 0.978 | +0.978 | quasi-TE (x) |
| 2 | 2.33374 | 0.033 | 0.933 | -0.933 | quasi-TM (y) |

Checks: both lie below the effective indices of the infinite 300 nm slab (TE 3.0489, TM 2.6158, from the transcendental equation, `slab_neff`) and approach the TE value
as the wire widens (error 0.096 at 1 um, falls by more than half at 2 um); the TE index changes by less than 1.5e-3 between 50 nm and 30 nm elements. Width sweep
(`soi_modes -- sweep`): TE0 2.0467 / 2.4635 / 2.6718 / 2.7862 / 2.8553 / 2.9003 / 2.9534 and TM0 2.0466 / 2.2248 / 2.3337 / 2.4037 / 2.4508 / 2.4840 / 2.5263 for 300 / 400 / 500 / 600 / 700 / 800 / 1000 nm.

* **Which modes you get:** the solver returns the `n_modes` modes nearest `n_target` (2.75 by default). Two modes take 0.1 s; a third that sits near the oxide index
  (here a hybrid at 1.63) takes about 15 s because the box modes of the closed domain cluster next to it. For wires much wider than 1 um put `n_target` just above the
  index of the wanted mode, or the solver returns higher lateral orders.
* **Hybrid polarization needs a broken mirror symmetry** (example 17, `soi_modes -- hybrid`): a 300 nm wide wire (TE and TM nearly degenerate) with the left wall sloped to 75 degrees gives
  modes with `s = (-0.586, +0.641, 0)` and `(+0.583, -0.647, 0)`, P = 0.87, TE shares 0.21 and 0.79: two linear states tilted by about +-45 degrees, nearly antipodal on the sphere.
  A side-wall angle applied to both walls keeps s2 = 0.

## Use case 2: Bragg gratings, symmetric and asymmetric (examples 18 and 19, `femcore::usecases::bragg`, `cargo run --release -p femcore --example bragg_gratings`)

One period of a sidewall grating as a 3D cell: the 500 nm x 300 nm wire has rectangular teeth 80 nm deep with a 50 % duty cycle. The 2D cross-section of the widest part is extruded
through two layers (teeth present, absent), the domain numbers of the layers choose the material of each tooth, the two ends carry identical meshes and are a Floquet pair with phase pi (the
band edge), and the other walls are perfect conductors 0.45 to 0.6 um from the wire. Kinds: `Symmetric` (teeth on both sides), `Asymmetric` (one side only; an
interpretation of "asymmetric", the other common one is `HalfPeriodShifted`: right teeth shifted by half a period, which couples TE0 to the odd mode), `Uniform` (no teeth, the reference).
The period is `lambda / (2 n_eff)` of the wire of mean width (280 nm; `matched_period`).

Results with `BraggGrating::coarse()` (second-order elements, 6696 cells, 44,838 unknowns, about 40 s each):

| grating | band-edge modes (THz) | gap | stop band | kappa from the gap | polarization of the edge modes |
|---|---|---|---|---|---|
| uniform | 198.03, 198.03 (degenerate) | 0 | | | quasi-TE, group G1 |
| symmetric | 194.680, 195.415 | 0.377 % | 5.8 nm at 1537 nm | 30 /mm | quasi-TE, TE 0.994 / 0.974, s2 = +0.001 |
| asymmetric (one side) | 196.236, 196.644 | 0.208 % | 3.2 nm at 1526 nm | 17 /mm | quasi-TE, TE 0.991 / 0.980, s2 = +-0.001 |

The quasi-TM band edges follow at 210.7 / 211.1 THz (symmetric). The coupling strength is `kappa = pi n_g (f_high - f_low) / c` with the group index of the unperturbed wire of mean width
(3.95, central differences of the 2D mode solver). The finer mesh (period 300 nm, 128,880 unknowns, 283 s) gives 185.625 / 186.457 THz (gap 0.45 %).

What to read from them, and what not to:

* The **one-sided grating has about 0.55 times the gap** of the two-sided one (the test asserts 0.3 to 0.8), as the halved perturbation suggests.
* The **band-edge modes of the one-sided grating are not rotated**: |s2| stays at 0.001 for the 500 nm wire (and below 0.02 on the 300 nm wire at its own matched period). TE0 and TM0
  differ by 0.34 in effective index, so the asymmetry hardly mixes them at one frequency. A polarization tilt comes from an asymmetric *cross-section* (use case 1, example 17), not from this corrugation.
* The **coupled-mode estimate** `kappa = 2 |n_wide - n_narrow| sin(pi duty) / lambda` from the effective indices of the two sections (206 /mm, 124 /mm) is
  **6 to 7 times larger than the FEM value**. It assumes sections much longer than the mode; here they are 140 nm, far shorter than the field, so the estimate is an upper bound
  (the test asserts `kappa < estimate`) and not a prediction.
* **Accuracy:** the 3D Floquet solver was checked on a structure with an exact answer (`tests/bragg_stack.rs`: a stack of two dielectrics seen as plane waves,
  band edges within 1 %, gap within 2 %, edges converging with refinement). The coarse mesh puts absolute frequencies about 1 % low (the uniform wire: 198.03 THz against 200.4 THz
  from the 2D index) and the gap, a small difference of two such numbers, is good to about 20 % (0.38 % coarse against 0.45 % fine). Use the fine settings (`BraggGrating::default()`, minutes) for numbers you want to quote.
* The factorisation of the complex Floquet matrix is most of the time (234 s of 283 s on the fine mesh). The block low-rank LU of `solvers::blr` is the port author's candidate drop-in for it; it is not wired in.

## What was not done

* The new solvers (`femcore::solvers`) and the adaptive mesher are not used by these use cases; the use cases run on the existing solvers.
* No substrate, no finite oxide thickness with leakage, no absorbing layers in the SOI example (perfectly conducting walls at a distance; the guided modes do not touch them).
* The grating has rectangular teeth on a structured layer mesh; no rounded corners, no apodization, no chirp, no transmission/reflection spectrum of a finite grating (the stop band and kappa are what the Floquet cell gives).
