# `convergence_check`

**Check mesh convergence** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the physics on 1 comparison mesh (2-grid default) or 2 comparison meshes (3-grid Roache GCI when `grids: 3`) with element sizes scaled by `factor` (`< 1` refines, `> 1` coarsens, or inverted coarsening ladder `[1.0, 1/factor, (1/factor)^2]` when `coarsen: true`; the active project is never mutated) and compares every eigenvalue, value and scalar across grids. By default `scale_sizes` is true so non-suppressed `add_size` rules (`h`, `max_deviation`, `min_layers`) and `gap_cells` scale proportionally on comparison copies. With `tolerance` each row says pass or fail, `converged` summarises, and the verdict is kept (`validated` becomes `checked` with tier `gci_verified` / `two_grid_checked` or stays `exploratory` / `non_convergent`). Sub-solves are isolated through `nyx-worker` and admitted against machine RAM when available. When mode ordering changes between meshes, `match: "nearest"` pairs mutual nearest eigenvalues and `track` filters candidate modes by target (`f`, `n_b` or `n_eff`) and/or physical mode signatures. With `extrapolate` (`"richardson_log_spd"` or `"richardson_linear"`), computes Richardson extrapolation (`extrapolate_order` defaults to 1.0 on 3D/sheet capacitance matrices, else `2 * order`) and promotes extrapolated capacitance matrices when `apply_extrapolated` is true (tagging `SpatialErrorEstimate.extrapolated_scalars_only = true`).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `apply_extrapolated` | true or false | no |  | when true, promote the extrapolated capacitance matrix onto the stored latest result of physics and mark extrapolated_scalars_only in error_budget, default false |
| `coarsen` | true or false | no |  | when true and factor < 1, invert factor into a coarsening ladder [1.0, 1/factor, (1/factor)^2] relative to the present mesh when memory headroom is tight, default false |
| `coarsen_factors` | list of number | no |  | explicit coarsening element-size multipliers (> 1.0) relative to the present mesh (1 entry for grids=2, 2 entries for grids=3) |
| `dry_run` | true or false | no |  | only build the comparison mesh(es) and answer their cells and what a solve of each is expected to cost; nothing is solved |
| `extrapolate` | text: `richardson_log_spd`, `richardson_linear` | no |  | optional two-grid Richardson extrapolation mode: richardson_log_spd (SPD-preserving log-space Maxwell capacitance extrapolation) or richardson_linear |
| `extrapolate_order` | number | no |  | asymptotic convergence order p > 0 for Richardson extrapolation and GCI theoretical order; default 1.0 for 3D/sheet capacitance matrices, else 2 * order |
| `factor` | number | no |  | comparison/present element-size ratio, positive and not 1 (< 1 refines, > 1 coarsens), default 0.6 |
| `force` | true or false | no |  | start convergence sub-solves even when machine memory estimate flags low headroom, default false |
| `grids` | integer: `2`, `3` | no |  | number of mesh levels to solve: 2 (default fast check) or 3 (Roache GCI with observed order and safety factor) |
| `isolate` | true or false | no |  | override worker process isolation for convergence sub-solves (default uses nyx-worker when configured) |
| `match` | text: `index`, `nearest` | no |  | how the solutions of the meshes are paired: index (default) compares solution i with solution i; nearest pairs each mode with the nearest eigenvalue of the other mesh in the complex plane and lists the modes without a partner |
| `max_cells` | number | no |  | refuse, before any solve, a check whose comparison mesh has more cells than this (an adaptive 2D mesh has 6 to 8 times as many cells at factor 0.6) |
| `max_seconds` | number | no |  | refuse, before any solve, a check estimated (from the solves this machine has run; no refusal without an estimate) to take longer than this for the comparison mesh |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `quantities` | list of text | no |  | non-empty list of quantity names to check; only matching rows count toward converged, worst_relative, and worst_quantity. each item: quantity name or scalar name to include in the verdict (e.g. "eigenvalue[0].re", "scalar C_0_0 (F)", or "C_0_0 (F)") |
| `scale_sizes` | true or false | no |  | scale non-suppressed add_size features (h, max_deviation, min_layers) proportionally with factor on comparison copies, default true |
| `tolerance` | number | no |  | largest accepted relative difference (or 3-grid GCI), e.g. 1e-3; omit to only report |
| `track` | object | no |  | follow mode(s) matching a target eigenvalue (f, n_b or n_eff) and/or physical mode-signature filters across meshes |
| `track.component` | text | no |  | field component name (e.g. Ex, Ey, Ez) whose share must be >= min_component |
| `track.f` | number | no |  | target eigenvalue real part (e.g. frequency in Hz) |
| `track.max_layer_fraction` | number | no |  | maximum absorbing-layer energy fraction (`layer fraction N`), in [0, 1] |
| `track.min_component` | number | no |  | minimum component share in [0, 1] |
| `track.min_share` | number | no |  | minimum energy share in `region`, in [0, 1] |
| `track.n_b` | number | no |  | target Bloch index real part |
| `track.n_eff` | number | no |  | target effective index real part |
| `track.region` | text | no |  | region name whose energy share (`share in <region> N`) must be >= min_share |
| `wait` | true or false | no |  | wait for the answer, default true; false starts a job that job_list, job_status and job_cancel see and job_wait answers |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"convergence_check","args":{"physics":"dfs"}}
```
