An error is JSON with an "error" code and a message that says what to change.

- invalid: the project or a command argument is wrong; the message names the region, condition or parameter and lists the valid choices. Fix that call; do not guess.
- unsupported: the physics does not offer it (a study, a condition, a property), or a material is lossy for a real-mode solver. The message says what it does offer.
- solver: the solver failed (singular system, no convergence). Check the boundary conditions (a static problem needs a fixed boundary), the mesh and the materials; the stage log of the run is in result_info (scalars) and in the run table.
- non_finite: NaN or infinity in the geometry or a material. Look for a zero-size feature or a value written without its unit.
- cancelled: a job was cancelled with job_cancel.
- A region name that does not exist: call geometry_names; names are those of add_overlay, and a region with absorbing layers also has pieces (use region*).
- A selector that selects nothing: use a boundary name from geometry_names, or check the coordinates in the selector.
- "no mesh" or "stale": mesh_build or solve_run again after an edit; result_info marks a result stale when the project changed.
- A material not found: list_materials. The error lists the names; add_material writes one of your own.
- Before quoting a number: convergence_check {physics, factor, tolerance} solves again on a finer mesh and compares every eigenvalue and scalar; do not quote an index that is not converged.
- Long solves: solve_run {wait:false} returns a job id; job_status shows the stage and fraction done; job_wait collects it.
- The server can be restarted by the host; open the stored workflow again with store_open {folder} (store_save returned the folder).
