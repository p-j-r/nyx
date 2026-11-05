nyx: finite-element photonics and phononics. Use it for numbers (effective indices, modes, frequencies, fields); do not estimate them from memory.

Conventions
- SI units everywhere (metres, hertz, kelvin). Wafer frame: x lateral, y propagation, z surface normal. A 2D model is the x-z cross-section.
- Regions and boundaries have names (Si, BOX, "BOX|Si", clad.top); you never see or pass integer labels. Conditions use selectors over the geometry or those names.
- Every result carries nyx.{session, revision, project_hash}. If project_hash is not the one you last saw, the project changed: call project_info.

How to work (one project is open; every tool is a command)
1. Geometry: add_rect / add_circle / add_polygon ... (list_operations shows all with fields), then add_overlay to name the regions in painter order (a later region wins where they overlap). geometry_names lists regions and boundaries. get_guide topic "geometry".
2. Physics: list_physics, describe_physics, then add_physics (ems = EM waves, sms = solid mechanics, dfs = diffusion) with its study. Materials: list_materials, then set_material per region (a library name, a name placed at an angle, or numbers). add_material writes your own as TOML.
3. Conditions: add_condition (kind and parameters from describe_physics). Mesh: set_mesh (h, adapt_scale, physics_sizing), mesh_build.
4. solve_run (a 2D section takes about a second). result_info lists solutions, fields and scalars; result_field reads values; analyse gives the polarization of optical modes.
5. store_save binds the project to a workspace folder; every later solve is a row in its run table, so the person sees your work in the dashboard.

Answering
- Quote the number with its assumptions: the geometry you built, the materials, the mesh size, the wavelength. Say which values the person did not give.
- Do not report a mode whose class does not match what was asked: check analyse. A result is stale when the project changed after it (result_info says so).

More: get_guide topics "geometry", "optical-modes" (reading mode tables), "troubleshooting" (errors).
