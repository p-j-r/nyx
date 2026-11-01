nyx: finite-element photonics. Use it for numbers (effective indices, modes, polarization); do not estimate them from memory.

Conventions
- SI units everywhere (metres, hertz). Wafer frame: x lateral, y propagation, z surface normal. A 2D model is the x-z cross-section.
- Regions and boundaries have names (Si, BOX, clad, "BOX|Si"); you never see or pass integer labels.
- Every result carries nyx.{session, revision, project_hash}. If project_hash is not what you last saw, something changed: call session_status.

How to work
1. One-line question about a waveguide mode: call recipe_run (recipe wg_modes_2d). It runs geometry, mesh, problem, solve, mode selection and checks and returns an answer block.
2. If a result says status "needs_input", ask the user that one question (it names the slot, why it matters and suggested values), then call again with the value. Use use_defaults:true only if the user says to assume.
3. To deviate (other mesh sizes, other wavelength, more modes) use the step tools in order: geom_template -> mesh_make -> problem_define -> solve_run -> result_modes -> check_run. Each takes the handle of the previous step (g_, m_, p_, r_).
4. Cheap tools: get_guide, session_status, result_modes. A solve of a 2D section takes about a second; check_run with mesh_converge solves once more.

Answering
- Quote the answer block: the number, the assumptions (assumed[]), and the checks. Say which values the user did not give.
- Do not report a mode whose class does not match what was asked; result_modes warns if a higher-index mode of another class exists.

More: get_guide topic "optical-modes" (reading mode tables), "troubleshooting" (errors).
