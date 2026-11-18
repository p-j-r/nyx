# nyx engine commands

Generated from `Engine::catalog()` by `nyx_engine::apiref` and checked by `crates/nyx-engine/tests/api_docs.rs`: **do not edit**. After a deliberate change of a command, rewrite with `NYX_WRITE_API=1 cargo test --release -p nyx-engine --test api_docs` (and the interface version, `engine_commands` in `docs/architecture/contracts/VERSIONS.md`, changes with it).

Catalog hash: `b3:f3346c8659d805385492b2c6ba82432b` · 99 commands. `nyx api list` and `nyx api show <command>` print the same text.

| Command | Access | What it does |
|---|---|---|
| [`add_box`](add_box.md) | edit | Box |
| [`add_circle`](add_circle.md) | edit | Circle |
| [`add_condition`](add_condition.md) | edit | Add a boundary or domain condition |
| [`add_cone`](add_cone.md) | edit | Cone frustum |
| [`add_cut`](add_cut.md) | edit | Cut by a plane |
| [`add_cylinder`](add_cylinder.md) | edit | Cylinder |
| [`add_difference`](add_difference.md) | edit | Difference |
| [`add_ellipse`](add_ellipse.md) | edit | Ellipse |
| [`add_extrude`](add_extrude.md) | edit | Extrude |
| [`add_intersect`](add_intersect.md) | edit | Intersection |
| [`add_material`](add_material.md) | edit | Add a material to the project |
| [`add_overlay`](add_overlay.md) | edit | Regions (overlay) |
| [`add_pattern`](add_pattern.md) | edit | Pattern |
| [`add_periodic`](add_periodic.md) | edit | Periodic pair |
| [`add_physics`](add_physics.md) | edit | Add a physics to the project |
| [`add_pml`](add_pml.md) | edit | Absorbing layers |
| [`add_pml_conditions`](add_pml_conditions.md) | edit | Add the absorbing-layer conditions |
| [`add_polygon`](add_polygon.md) | edit | Polygon |
| [`add_rect`](add_rect.md) | edit | Rectangle |
| [`add_revolve`](add_revolve.md) | edit | Revolve |
| [`add_rounded_rect`](add_rounded_rect.md) | edit | Rounded rectangle |
| [`add_size`](add_size.md) | edit | Mesh size |
| [`add_sphere`](add_sphere.md) | edit | Sphere |
| [`add_tag`](add_tag.md) | edit | Tag |
| [`add_torus`](add_torus.md) | edit | Torus |
| [`add_transform`](add_transform.md) | edit | Transform |
| [`add_union`](add_union.md) | edit | Union |
| [`analyse`](analyse.md) | read | Analyse a result |
| [`check_project`](check_project.md) | read | Check the project |
| [`confinement`](confinement.md) | read | Fraction of a field in a region |
| [`convergence_check`](convergence_check.md) | expensive | Check mesh convergence |
| [`describe_material`](describe_material.md) | read | Describe a material |
| [`describe_physics`](describe_physics.md) | read | Describe a physics |
| [`export_results`](export_results.md) | expensive | Write the results to an HDF5 file |
| [`extremum`](extremum.md) | read | Largest or smallest value of an expression |
| [`field`](field.md) | read | A handle to a field |
| [`field_op`](field_op.md) | read | Arithmetic on handles |
| [`geometry_names`](geometry_names.md) | read | Evaluate the geometry and list its names |
| [`get_guide`](get_guide.md) | read | Read a nyx guide |
| [`integrate`](integrate.md) | read | Integrate an expression over a region |
| [`job_cancel`](job_cancel.md) | expensive | Cancel a job |
| [`job_list`](job_list.md) | read | List the jobs |
| [`job_status`](job_status.md) | read | Status of a job |
| [`job_wait`](job_wait.md) | expensive | Wait for a job |
| [`line_cut`](line_cut.md) | read | Values of an expression along a line |
| [`list_analyses`](list_analyses.md) | read | List the analyses |
| [`list_commands`](list_commands.md) | read | List the commands |
| [`list_features`](list_features.md) | read | List the features |
| [`list_materials`](list_materials.md) | read | List the materials |
| [`list_operations`](list_operations.md) | read | List the geometry operations |
| [`list_physics`](list_physics.md) | read | List the physics of this build |
| [`mesh_build`](mesh_build.md) | expensive | Mesh the geometry |
| [`open_results`](open_results.md) | edit | Open an HDF5 results file |
| [`overlap`](overlap.md) | read | Overlap of two fields |
| [`poll_jobs`](poll_jobs.md) | read | Apply finished jobs |
| [`probe`](probe.md) | read | Value of an expression at points |
| [`project_info`](project_info.md) | read | Project summary |
| [`project_json`](project_json.md) | read | Project as JSON |
| [`project_load`](project_load.md) | edit | Load a project from JSON |
| [`project_new`](project_new.md) | edit | New project |
| [`project_open`](project_open.md) | edit | Open a project file |
| [`project_save`](project_save.md) | expensive | Save the project to a file |
| [`queue_cancel`](queue_cancel.md) | expensive | Cancel a queued or running job |
| [`queue_get`](queue_get.md) | read | One job of the queue |
| [`queue_hold`](queue_hold.md) | expensive | Hold a queued job |
| [`queue_list`](queue_list.md) | read | List the queue |
| [`queue_pause`](queue_pause.md) | expensive | Pause the queue |
| [`queue_release`](queue_release.md) | expensive | Release a held job |
| [`queue_resume`](queue_resume.md) | expensive | Resume the queue |
| [`queue_set_priority`](queue_set_priority.md) | expensive | Change the priority of a queued job |
| [`queue_submit`](queue_submit.md) | expensive | Submit work to the queue |
| [`queue_wait`](queue_wait.md) | expensive | Wait for a job, at most a while |
| [`remove_condition`](remove_condition.md) | edit | Remove a condition |
| [`remove_feature`](remove_feature.md) | edit | Remove a feature |
| [`remove_material`](remove_material.md) | edit | Remove a project material |
| [`remove_param`](remove_param.md) | edit | Remove a model parameter |
| [`remove_physics`](remove_physics.md) | edit | Remove a physics from the project |
| [`rename_feature`](rename_feature.md) | edit | Rename a feature |
| [`result_field`](result_field.md) | read | Read a field |
| [`result_info`](result_info.md) | read | Describe the results |
| [`results_file_info`](results_file_info.md) | read | Describe an HDF5 results file |
| [`set_active`](set_active.md) | edit | Set the active physics |
| [`set_geometry`](set_geometry.md) | edit | Replace the geometry |
| [`set_material`](set_material.md) | edit | Assign a material to a region |
| [`set_mesh`](set_mesh.md) | edit | Set the mesh settings |
| [`set_param`](set_param.md) | edit | Set a model parameter |
| [`set_physics_params`](set_physics_params.md) | edit | Set the options of a physics |
| [`set_problem`](set_problem.md) | edit | Replace the problem of a physics |
| [`set_study`](set_study.md) | edit | Set the study of a physics |
| [`solve_run`](solve_run.md) | expensive | Solve |
| [`store_open`](store_open.md) | edit | Open a workflow |
| [`store_runs`](store_runs.md) | read | Run table |
| [`store_save`](store_save.md) | expensive | Save into a workspace |
| [`suppress_feature`](suppress_feature.md) | edit | Suppress or restore a feature |
| [`sweep_cancel`](sweep_cancel.md) | expensive | Cancel the points of a sweep that have not finished |
| [`sweep_run`](sweep_run.md) | expensive | Run a sweep: one recorded run per point |
| [`sweep_table`](sweep_table.md) | read | The table of a sweep |
| [`tag_feature`](tag_feature.md) | edit | Tag a feature |
| [`update_feature`](update_feature.md) | edit | Change a feature |
