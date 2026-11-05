Optical modes of a 2D cross-section (x-z plane; the mode travels along y) with physics ems, study "modes".

- Study parameters: freq (Hz; c / wavelength), n_target (the solver returns the modes whose effective index is nearest it; put it just below the highest index of the cores, for silicon on oxide about 2.8 to 3.0), n_modes (keep it small; a mode near the cladding index is slow), evanescent.
- Materials come from the library at the study wavelength: set_material {physics:"ems", region:"Si", material:"silicon"}. Dispersion models are evaluated at c / freq, so the wavelength you set is the one used. A lossy material is refused by the real-mode solver; use an absorbing layer (add_pml) for loss.
- The outer boundary is a perfect conductor unless you add a pml condition; keep it at least a micron from the core.
- Mesh: h about 0.1 um with adapt_scale 0.4 resolves a 500 x 300 nm silicon wire; set_mesh {physics_sizing:true} sizes elements from the wavelength in each material. Repeat with a finer mesh to check the number (the TE index of the reference wire changes by about 1e-3 between 50 nm and 30 nm elements).
- result_info: each solution has eigenvalue [re, im] = effective index. Real materials give real indices; loss in dB/cm = 8.686 * k0 * Im(neff) / 100 with k0 = 2 pi / wavelength in 1/m.
- analyse {analysis:"polarization", physics:"ems"}: class per mode from the Stokes parameters of the complex field: quasi-TE (field mostly along x), quasi-TM (mostly along z) or hybrid; te_fraction, degree of polarization, degenerate groups. "TE0" means the quasi-TE mode with the highest index.
- Exactly degenerate modes (a square wire) come back as an arbitrary mixture and their classes are not reliable. Make the width differ from the thickness, or say so.
- Reference: a 500 nm x 300 nm silicon wire in silica at 1550 nm has TE0 about 2.67 and TM0 about 2.33 (silicon 3.4757, silica 1.444); wider wires approach the slab index.
