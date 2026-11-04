# Golden changes

Golden files are the numerical oracle (ADR 10). A change is listed here with the reason; none is made to hide a difference.

## 2026-11-04, `emw_pml_x.json` (phase 3, linear algebra; decided by paul: adopt the extra mode for now)

* **What changed**: the file had 3 modes, it has 4. The new 4th mode is `neff 4 = 1.730216 -5.076e-6i` (Q = Re/(2|Im|) = 1.704e5, quasi-TM (z), TE 0.000, Ey 0.188, P 1.000). Modes 1 to 3 are the same modes: `neff` 1.4393617 -2.791e-5i, 1.5762346 +7.334e-6i, 1.6630504 -5.812e-5i agree to 1e-12 relative in the real part (their tiny imaginary parts, 1e-5 to 1e-6, differ in the 7th digit, so Q 2 moves by 4e-7 relative); the fields of modes 1 to 3 differ by rounding only.
* **Why**: the study asks for 4 modes, and the old file holds 3. `ems_pml::solve_mode_pairs_complex` drops a pair whose own error estimate is above its threshold, and its doc says so ("pairs the Krylov process could not converge are dropped"). The old solver therefore delivered only three acceptable pairs; I did not trace which of its limits (the 400-vector basis cap, the tolerance, the start block) left the fourth unconverged. The library's block Krylov-Schur (`solvers::driver::eigs_complex`, ADR 33) converges four pairs at the same tolerance (1e-10), and the fourth is the next eigenvalue of the list (0.067 in `neff` above mode 3), a genuine eigenvalue of the pencil. It is kept for now (paul's decision); whether the old file was an artefact of the old solver or the new solver adds a spurious pair has not been established independently (a check with a much tighter tolerance, or the mode field inside the PML, would settle it).
* **How to undo**: `NYX_OLD_LINALG` no longer exists (the old solver is deleted); restore the file from the commit before `Phase 3: adopt the extra converged mode` and lower the study's `n_modes` to 3.

## `solid_kz_pml_z.json` (not changed)

The same effect (the library converges one more mode: 5 solutions against 4) appears in `solid_kz_pml_z`. `sms_kz` (out-of-plane wave number) is being removed from the repository by another thread, so the golden is left as it is and the test `golden_solid_kz_sapphire_with_a_pml_along_z` is expected to fail until it goes with the model.

**Check made after adopting it**: with `NYX_EIG_TOL=1e-13` the study returns 3 modes again, like the old file. The 4th pair is therefore marginal: it meets the eigensolver tolerance at the kernel's 1e-10 and not at 1e-13, so whether it is listed depends on the tolerance (the same dependence as in the old solver, at a different point). This does not show that it is spurious; it shows that the mode count of this PML study is tolerance-sensitive (gate G2 option B: loose tolerances, so the count may change with the tolerance chosen).
