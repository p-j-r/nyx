An error is JSON with an "error" code and a message that says what to change.

- invalid: the project or a command argument is wrong; the message names the region, condition or parameter and lists the valid choices. Fix that call; do not guess.
- unsupported: the physics does not offer it (a study, a condition, a property), or a material is lossy for a real-mode solver. The message says what it does offer.
- solver: the solver failed (singular system, no convergence). Check the boundary conditions (a static problem needs a fixed boundary), the mesh and the materials; the stage log of the run is in result_info (scalars) and in the run table.
- non_finite: NaN or infinity in the geometry or a material. Look for a zero-size feature or a value written without its unit.
- cancelled: a job was cancelled with job_cancel.
- A region name that does not exist: call geometry_names; names are those of add_overlay, and a region with absorbing layers also has pieces (use region*).
- A selector that selects nothing: use a boundary name from geometry_names, or check the coordinates in the selector.
- "no mesh" or "stale": mesh_build or solve_run again after an edit; result_info marks a result stale when what it was made from changed (an edit of another physics does not stale it).
- A material not found: list_materials. The error lists the names; add_material writes one of your own.
- Before quoting a number: convergence_check {physics, factor, tolerance} solves again on a finer mesh and compares every eigenvalue and scalar; do not quote an index that is not converged.
- Long solves: solve_run {wait:false} returns a job id; job_status shows the stage, the seconds in it (`stage_seconds`) and `cancel_requested`; job_wait collects it. The same holds for mesh_build, convergence_check, export_results and open_results with {wait:false}. A cancel is seen at the next check: a factorisation has to finish first, so `cancel_requested` can stay true for as long as that stage runs.
- The server can be restarted by the host; open the stored workflow again with store_open {folder} (store_save returned the folder).
- A long wait: one script, not many calls (each costs a model turn). wait_s is at most 90; the answer "running" carries the job (stage, seconds, estimate); script_wait continues, script_wait cancel stops a solve that is too large.
- A mesh that is too large is found by solving it. mesh_build answers `cost`: aim at 2000 to 4000 cells for a 2D mode with order 2 (points_per_wavelength is 6 by default; h 200 to 250 nm and adapt_scale 0.85 to 1 are a start; both a high points_per_wavelength and a small adapt_scale multiply the cells). A crystal with off-diagonal permittivity (a tilted sapphire) or an absorbing layer runs the complex solver, about 4 times the work.
- convergence_check {dry_run:true} builds both meshes and tells the cells and the expected seconds (null until this machine has solved a few); {max_cells, max_seconds} refuse the check before any solve. Factor 0.6 gives an adaptive 2D mesh 6 to 8 times the cells.
- A lossless guide on a tilted crystal has an imaginary part of 1e-10 and Q of 1e10 that are rounding: convergence_check marks those rows `noise` and leaves them out of `converged`; do not quote that Q.
- A script is refused when a loop variable is given different literal lists in turn: use a new name for each loop.
