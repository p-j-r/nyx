# Golden changes

Golden files are the numerical oracle (ADR 10). A change is listed here with the reason; none is made to hide a difference.

## `solid_kz_pml_z.json` (not changed)

The same effect (the library converges one more mode: 5 solutions against 4) appears in `solid_kz_pml_z`. `sms_kz` (out-of-plane wave number) is being removed from the repository by another thread, so the golden is left as it is and the test `golden_solid_kz_sapphire_with_a_pml_along_z` is expected to fail until it goes with the model.

**Check made after adopting it**: with `NYX_EIG_TOL=1e-13` the study returns 3 modes again, like the old file. The 4th pair is therefore marginal: it meets the eigensolver tolerance at the kernel's 1e-10 and not at 1e-13, so whether it is listed depends on the tolerance (the same dependence as in the old solver, at a different point). This does not show that it is spurious; it shows that the mode count of this PML study is tolerance-sensitive (gate G2 option B: loose tolerances, so the count may change with the tolerance chosen).
