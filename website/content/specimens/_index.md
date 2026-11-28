+++
title = "Specimens & Field Atlas"
description = "Verified finite-element specimens across 2D vector photonics, cylindrical PML resonators, elastic MEMS anchor loss, and 3D open cavities."
template = "specimens.html"
+++

Every specimen in this atlas is drawn directly by `nyx.figure(...)` from the unmodified example projects in `examples/`, using SI metres in the **wafer frame**: `x` and `y` span the wafer plane, `z` is the surface normal, and `+y` is the propagation axis of a waveguide (`exp(-j beta y)`). A 2D cross-section therefore lies in the `(x, z)` plane, and all result fields (`Ex`, `Ey`, `Ez`, `ux`, `uz`) are exported in wafer coordinates.

## 2D Vector Photonics & Cylindrical PML

- **`15_soi_wire` (`ems::modes`, P2):** A `500 nm × 300 nm` silicon core (`w = 5e-7 m`, `t = 3e-7 m`, `slab = 0.0`) embedded in silica cladding and buried oxide, solved at `lambda = 1550 nm` (`freq = 193.414 THz`) with local near-field mesh grading (`h = 40 nm` in `Si`, `90 nm` in `BOX_near`/`clad_near`, `300 nm` outer). Hierarchical degree-2 Nédélec transverse edge elements paired with P2 longitudinal `phi` yield `neff 1 = 2.671811` (`n_g = 4.1969`, `94.5%` power in `Si`) and `neff 2 = 2.333849` (`n_g = 4.6294`, `85.1%` in `Si`) in `0.08 s` across `986` triangles (`2,015` DOFs).
- **`18_wgm_cylindrical_pml` (`ems::modes`, P2):** A `1.2 um` radius dielectric core (`eps_r = 12.0`) inside a `2.2 um` air cladding ringed by an `800 nm` 5-layer structured quadrilateral cylindrical PML (`kappa_max = 2.0`, `order = 2.0`, `target_reflection = 1e-6`). The complex shift-and-invert solver resolves four whispering-gallery modes (`5,225` DOFs) in `0.38 s` down to a maximum algebraic residual of `8.716e-14`.

## Elastic Anchor Loss & 3D Open Cavities

- **`11_anchor_loss` (`sms::eigenfrequency`, P2 plane strain):** A `20 mm × 300 mm` elastic beam (`E = 200 GPa`, `nu = 0.3`, `rho = 7850 kg/m^3`) anchored onto a `600 mm × 350 mm` substrate flanked on `-x`, `+x`, and `-z` by 5-layer Cartesian PML wings (`thickness = 0.2 m`, `kappa_max = 2.0`, `strength = 3.0`). Across `4,832` DOFs (`0.086 s`, max residual `6.255e-12`), `nyx-sms` separates high-`Q` cantilever flexure (`183.703 Hz`, `Q = 14,444.6`, `layer fraction = 2.71e-7`) from substrate-radiating modes (`2.286–2.762 kHz`, `Q = 1.46–1.61`, `layer fraction = 0.27–0.52`).
- **`20_cavity_3d_spherical_pml` (`ems::cavity`):** An open `25 mm` radius air sphere wrapped in a `10 mm` 5-layer spherical PML shell (`r_in = 0.025 m`, `strength = 16.480`). `projected_shift_invert_complex` projects out zero-frequency gradient modes and solves `3,622` DOFs in `0.68 s` (max residual `1.054e-8`), recovering the triply-degenerate fundamental leaky resonance at `1.642024e9 + 1.700781e9i Hz` (`Q = 0.483`) and the higher mode at `2.824384e9 + 2.557435e9i Hz` (`Q = 0.552`). Both `cutaway3d` and `slice3d` plates evaluate expressions on un-averaged 8-subtet P2 `CellFields`, preserving exact normal-field jumps across material boundaries.
