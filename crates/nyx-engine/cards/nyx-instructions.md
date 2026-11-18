nyx: finite-element photonics and phononics. Use it for numbers (effective indices, modes, frequencies, fields); do not estimate them from memory.

Conventions
- SI units everywhere (metres, hertz, kelvin). Wafer frame: x lateral, y propagation, z surface normal. A 2D model is the x-z cross-section.
- Regions and boundaries have names (Si, BOX, "BOX|Si", clad.top); you never see or pass integer labels. Conditions use selectors over the geometry or those names.
- Every result carries nyx.{session, revision, project_hash}. If project_hash is not the one you last saw, the project changed: call project_info.

How to work (one project is open; every tool is a command)
1. Geometry: add_rect / add_circle / add_polygon ... (list_operations shows all with fields), then add_overlay to name the regions in painter order (a later region wins where they overlap). geometry_names lists regions and boundaries. get_guide topic "geometry".
2. Physics: list_physics, describe_physics, then add_physics (ems = EM waves, sms = solid mechanics, dfs = diffusion) with its study. Materials: list_materials, then set_material per region (a library name, a name placed at an angle, or numbers). add_material writes your own as TOML.
3. Conditions: add_condition (kind and parameters from describe_physics). Mesh: set_mesh (h, adapt_scale, physics_sizing), mesh_build.
4. solve_run (a 2D section takes about a second). convergence_check repeats it on a finer mesh and compares: run it before you quote a number, with a `tolerance`: the result then answers `validated: checked` (a result nobody has checked is `exploratory`, one whose check would not fit the machine is `check_not_affordable`, a stale one is never checked; every answer that carries `stale` carries `validated`). result_info lists solutions, fields and scalars; result_field reads values (sampled, node-averaged: for a number, do not add them up yourself); analyse gives the polarization of optical modes. For numbers from a field use the reductions, which run in the engine on the cells the solver made: integrate (an integral of a field or of abs2(Ex) over a region), extremum, probe and line_cut (values at points), overlap and confinement (the share of |E|^2 in a region). A field is named by a cell field of result_info (`cell_fields`, e.g. Ex) and combined with field_op; `revolution: true` weights a body of revolution by 2 pi r.
5. store_save binds the project to a workspace folder; every later solve is a row in its run table, so the person sees your work in the dashboard.

Answering
- Quote the number with its assumptions: the geometry you built, the materials, the mesh size, the wavelength. Say which values the person did not give.
- Do not report a mode whose class does not match what was asked: check analyse. A result is stale when what it was made from changed after it: its own problem, the materials, the mesh or a physics it takes fields from (result_info, result_field, analyse and list_analyses say `stale`; do not quote a stale number).

More: get_guide topics "geometry", "optical-modes" (reading mode tables), "troubleshooting" (errors).
