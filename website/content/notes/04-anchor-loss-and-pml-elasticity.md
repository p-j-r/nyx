+++
title = "Elastic Perfectly Matched Layers, Axisymmetric Regularity & MEMS Anchor Loss"
date = 2026-11-25
description = "Curvilinear complex coordinate stretching across Cartesian, cylindrical, and spherical topologies, Christoffel backward-wave checks, and Bindel–Govindjee MEMS anchor loss."

[extra]
note_number = "04"
crates = "femcore::pml · nyx-sms::pml · nyx-sms::axisym"
reading_time = "8 min"
+++

A mechanical resonator anchored to a semi-infinite wafer loses energy by radiating elastic waves through its support post into the bulk substrate. Truncating the substrate with a rigid clamp reflects those outgoing waves back into the resonator, destroying the quality factor $Q$; truncating with a simple viscous dashpot (`lrb`) absorbs waves cleanly only near normal incidence and turns the eigenproblem quadratic in $\omega$.

In `crates/femcore/src/pml.rs` and `crates/nyx-sms/src/pml.rs`, `nyx` surrounds the finite substrate with a **curvilinear complex-coordinate-stretched Perfectly Matched Layer (PML)** across Cartesian, cylindrical, and spherical topologies, yielding a complex linear generalised eigenproblem $K(\tilde{\mathbf{x}})\mathbf{u} = \omega^2 M(\tilde{\mathbf{x}})\mathbf{u}$ whose complex frequencies $f = f_{\mathrm{re}} + i f_{\mathrm{im}}$ directly measure radiation damping via $Q = \mathrm{Re}(f) / (2|\mathrm{Im}(f)|)$.

## Complex Coordinate Stretching & Transformation Pullback

Under `nyx`'s $e^{+i\omega t}$ engineering time convention, outgoing waves decay inside the absorbing layer when real spatial coordinates are continued into the fourth quadrant with a complex stretch function along each active normal direction:

$$
s(t) = \kappa(t) - i \psi(t), \qquad \kappa(t) \ge 1, \quad \psi(t) \ge 0, \quad t \in [0, 1]
$$

where $t = 0$ marks the physical interface and $t = 1$ the clamped outer PML wall (`OUTER_T_TOL = 1e-5`). At every quadrature point, `femcore::pml::pml_jacobian` forms the $3 \times 3$ metric Jacobian $J_s = \partial \tilde{\mathbf{x}} / \partial \mathbf{x}$, its analytical inverse $J_s^{-1}$, and its complex determinant $\det J_s$. Because variational energy integration uses the **non-conjugated complex-symmetric bilinear form**, strain gradients pull back via $J_s^{-T}$ and volume elements scale by $\det J_s$ without complex conjugation.{% sidenote(id="sn-pml-order", label="1") %}Polynomial grading (`order = 2.0`, `kappa_max = 2.0`) avoids reflection at the physical-to-PML interface where $s(0) = 1$, while real evanescent stretching $\kappa(t) > 1$ accelerates decay of non-propagating near-field tails.{% end %}

## Cantilever Anchor Loss (`11_anchor_loss`)

In `examples/11_anchor_loss.nyx.json`, a `20 mm × 300 mm` elastic beam ($E = 200\,\text{GPa}$, $\nu = 0.3$, $\rho = 7850\,\text{kg/m}^3$, plane strain, P2) rests on a `600 mm × 350 mm` substrate flanked on `-x`, `+x`, and `-z` by 5-layer Cartesian PML wings (`thickness = 0.2 m`, `strength = 3.0`). Solving `4,832` DOFs takes `0.086 s` (`9` restarts, max algebraic residual `6.255e-12`).

{{ plate(src="figures/anchor_loss_modes.png", number="PLATE 03", alt="Six elastic-PML eigenmodes of a beam anchored to a PML substrate", caption="Modes 1–6 of 11_anchor_loss (4,832 DOFs, 0.086 s): flexural beam modes 1, 2, and 6 (Q = 14,444.6, 1,239.9, and 119.3) contrasted with strongly damped substrate-PML modes 3–5 (Q = 1.46–1.61).", meta="11_anchor_loss · sms::eigenfrequency · 4,832 DOFs") }}

{% solver_log(cmd="nyx run render_site_figures.nys  # 11_anchor_loss") %}
sms: eigenfrequency, order=2, dim=2, dofs=4832, seconds=0.0858
  solver_kind: shift_invert_complex, converged: true, restarts: 9
  max_algebraic_residual: 6.254677404565312e-12
  mode 1: 183.703 Hz + 6.359e-3i Hz   (Q = 14444.57, layer_frac = 2.714e-07, P_pml = 553.37 W)
  mode 2: 1132.207 Hz + 0.4566i Hz    (Q =  1239.95, layer_frac = 9.086e-05, P_pml = 251.08 kW)
  mode 3: 2286.157 Hz + 710.22i Hz    (Q =     1.61, layer_frac = 2.715e-01, P_pml = 36.66 GW)
  mode 4: 2410.962 Hz + 827.02i Hz    (Q =     1.46, layer_frac = 3.537e-01, P_pml = 105.59 GW)
  mode 5: 2762.168 Hz + 925.32i Hz    (Q =     1.49, layer_frac = 5.202e-01, P_pml = 87.59 GW)
  mode 6: 3096.439 Hz + 12.979i Hz    (Q =   119.28, layer_frac = 5.753e-03, P_pml = 20.48 MW)
  unstable_modes = 0
{% end %}

Notice how `layer fraction` — the share of kinetic energy inside the PML domains — cleanly distinguishes physical resonator modes (`2.71e-7` on mode 1, `9.09e-5` on mode 2, `5.75e-3` on mode 6) from continuous-spectrum Berenger substrate modes (`0.27–0.52` on modes 3–5, `Q < 2`).

## Christoffel Backward-Wave Stabilisation & Axisymmetric Disks

In anisotropic crystals, slowness surfaces can exhibit concave regions where phase and group velocities point in opposite directions along a PML normal ($\mathbf{v}_p \cdot \mathbf{v}_g < 0$). Left unchecked, standard stretching turns a backward wave into an exponentially *growing* mode (`Im f < 0`, `unstable_modes > 0`). Before assembly, `nyx-sms::pml` calls `backward_wave_mpml_ratio` (`femcore::meshing::physics_size`) on the wafer-rotated stiffness tensor `domain_wafer_stiffness` to detect backward-wave cones and apply real stretching stabilization (`PmlStabilization`).

For bodies of revolution (`sms::axisymmetric_eigenfrequency`, `crates/nyx-sms/src/axisym.rs`), displacement varies as $e^{-i m \theta}$ in $(r, z, \theta)$. Whenever a mesh touches the symmetry axis $r = 0$, regularity forces exact Dirichlet constraints (`axis r = 0: N dofs, K held at zero`):

- $m = 0$ (torsional / radial-axial): $u_r = u_\theta = 0$ at $r = 0$.
- $|m| = 1$ (flexural): $u_z = 0$ at $r = 0$.
- $|m| \ge 2$ (whispering-gallery / higher azimuthal): $u_r = u_z = u_\theta = 0$ at $r = 0$.

Cylindrical PML wings stretch $r \to \tilde{r}$ on $+x$, $+z$, and $-z$ (never towards $r = 0$), scaling $1/r$ hoop-strain terms and volume weights by $\tilde{r}/r$. In `crates/nyx-engine/tests/validation_mems_piezo.rs`, this axisymmetric PML kernel is verified against Bindel & Govindjee (2005) on a polysilicon radial-contour disk-on-post ($R = 6\,\mu\text{m}$, $f_0 = 450.155\,\text{MHz}$) for both closed-form frequency and PML thickness/strength invariance.
