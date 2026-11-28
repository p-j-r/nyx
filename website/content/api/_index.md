+++
title = "Engine API & 118 Commands"
description = "Single versioned JSON-Schema command catalog driving the desktop window, Monty scripts, headless batch CLI, and MCP server."
template = "api.html"
+++

Every operation in `nyx` — from creating a parametric wafer-frame rectangle or importing a GDSII layout to building an adaptive mesh, submitting a sweep, running a Craig–Bampton circuit quantisation, or exporting an HDF5 results file — is an **engine command** registered in `Engine::catalog()` (`crates/nyx-engine`).

The current catalog holds **118 commands** (`b3:f4c059caae7ee3100bbf472db2c560bd`), each backed by a JSON Schema and a generated Markdown page in `docs/api/`. Because the egui desktop window, `.nys` Monty scripts (`nyx-script`), JSON-lines command lists (`nyx run`), and the MCP server (`nyx-mcp`) all dispatch through the same `Engine::call` boundary, an analysis run interactively in the window is bit-for-bit reproducible from a headless script.

## Inspecting the Catalog from the Terminal

```sh
nyx api list                 # print all 118 commands and their access tier
nyx api show figure          # full schema, defaults, and example for `figure`
nyx api show solve_run --json
```

Each command declares an explicit access tier (`read`, `edit`, or `expensive`) so dry-run checkers (`nyx check FILE`) and AI hosts can distinguish side-effect-free queries from mesh builds, solves, and disk exports.

## Command Groups

- **Parametric Geometry & GDSII (`38` commands):** `add_rect`, `add_circle`, `add_ellipse`, `add_polygon`, `add_polygons`, `add_rounded_rect`, `add_box`, `add_cylinder`, `add_cone`, `add_sphere`, `add_torus`, `add_extrude`, `add_revolve`, `add_union`, `add_difference`, `add_intersect`, `add_cut`, `add_transform`, `add_pattern`, `add_overlay`, `add_chip_stack`, `add_sheet`, `add_pml`, `add_absorbing_layer`, `add_periodic`, `add_size`, `add_tag`, `import_gds`, `geometry_names`, `set_geometry`, `set_param`, `remove_param`, `rename_feature`, `remove_feature`, `suppress_feature`, `tag_feature`, `list_features`, `list_operations`.
- **Physics, Materials & Mesh (`24` commands):** `add_physics`, `remove_physics`, `set_active`, `set_problem`, `set_study`, `set_physics_params`, `set_coupled_loop`, `add_condition`, `remove_condition`, `add_pml_conditions`, `add_material`, `remove_material`, `set_material`, `describe_material`, `list_materials`, `evaluate_material`, `bcs_surface_impedance`, `list_remote_materials`, `search_remote_materials`, `fetch_material`, `import_material`, `set_mesh`, `mesh_build`, `check_project`.
- **Solves, Jobs, Sweeps & Optimization (`24` commands):** `solve_run`, `convergence_check`, `eme_sweep`, `sensitivity_run`, `optimize_run`, `sweep_run`, `sweep_table`, `sweep_cancel`, `job_list`, `job_status`, `job_wait`, `job_cancel`, `poll_jobs`, `queue_submit`, `queue_list`, `queue_get`, `queue_wait`, `queue_hold`, `queue_release`, `queue_pause`, `queue_resume`, `queue_set_priority`, `queue_cancel`, `pml_diagnostics`.
- **Fields, Analyses, Figures & HDF5 Storage (`30` commands):** `result_info`, `result_field`, `field`, `field_op`, `probe`, `line_cut`, `integrate`, `extremum`, `confinement`, `overlap`, `analyse`, `list_analyses`, `figure`, `export_results`, `open_results`, `results_file_info`, `project_new`, `project_open`, `project_load`, `project_save`, `project_info`, `project_json`, `project_script`, `store_open`, `store_save`, `store_runs`, `describe_physics`, `list_physics`, `list_commands`, `get_guide`.
