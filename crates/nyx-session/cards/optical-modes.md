Optical modes of a 2D cross-section (x-z plane; the mode travels along y).

- The solver returns the n_modes modes whose effective index is nearest the search index (n_target). With n_target omitted nyx places it just above the analytic slab index of the silicon layer, an upper bound for a wire of that thickness, so the fundamental mode is found first.
- Real materials give real neff. A complex neff would mean loss; loss in dB/cm = 8.686 * k0 * Im(neff) / 100 with k0 = 2 pi / wavelength in 1/m.
- Class comes from the Stokes parameters of the complex field: quasi-TE = field mostly along x, quasi-TM = mostly along z, hybrid = polarization changes across the mode. "TE0" means the quasi-TE mode with the highest neff.
- Exactly degenerate modes (a square wire) come back as an arbitrary mixture; their classes are then not reliable. Make the width differ from the thickness, or ask the user.
- Defaults you may see in the ledger: 1550 nm, silicon 3.4757, SiO2 1.4440 (as in the repository's SOI use case), perfectly conducting window 1 um from the core, no silicon substrate under the oxide.
- Mesh: 40 nm elements in the silicon, 90 nm near it, 300 nm far (the repository's converged setting: the TE index changes by less than 1.5e-3 between 50 nm and 30 nm). check_run mesh_converge repeats the solve on a 40 percent finer mesh.
- Reference (documented in the repository for an oxide-clad 500 nm x 300 nm wire on 2 um oxide): TE0 2.6718, TM0 2.3337. Wider wires approach the slab index.
- A third mode close to the cladding index can take about 15 s; keep n_modes small.
- Wide wires (about 1.5 um and more at 300 nm thick) have many TE modes crowded near the TM slab index; a "TM0" can then be a TE/TM mixture above the TM slab bound. recipe_run reports status "unverified" when slab_bound or mesh_converge fails; say so instead of quoting the number.
