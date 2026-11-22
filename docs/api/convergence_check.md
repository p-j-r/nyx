# `convergence_check`

**Check mesh convergence** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the physics on a second mesh with element sizes multiplied by `factor` (`< 1` refines, `> 1` coarsens; the active project is never mutated) and compares every eigenvalue, value and scalar between the coarse and fine grids. By default `scale_sizes` is true so non-suppressed `add_size` rules (`h`, `max_deviation`, `min_layers`) and `gap_cells` scale proportionally with `factor` on the comparison copy. With `tolerance` each row says pass or fail, `converged` summarises, and the verdict is kept (`validated` becomes `checked` or stays `exploratory`). When mode ordering changes between meshes, `match: "nearest"` pairs mutual nearest eigenvalues and `track` filters candidate modes on both meshes by target (`f`, `n_b` or `n_eff`) and/or physical mode signatures (`region`/`min_share`, `component`/`min_component`, `max_layer_fraction`). With `extrapolate` (`"richardson_log_spd"` or `"richardson_linear"`), computes two-grid Richardson extrapolation (`extrapolate_order` defaults to 1.0 on 3D/sheet capacitance matrices, else `2 * order`) and optionally promotes the extrapolated capacitance matrix onto the stored result when `apply_extrapolated` is true.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `apply_extrapolated` | true or false | no |  | when true, promote the extrapolated capacitance matrix onto the stored latest result of physics, default false |
| `dry_run` | true or false | no |  | only build both meshes and answer their cells and what a solve of each is expected to cost; nothing is solved |
| `extrapolate` | text: `richardson_log_spd`, `richardson_linear` | no |  | optional two-grid Richardson extrapolation mode: richardson_log_spd (SPD-preserving log-space Maxwell capacitance extrapolation) or richardson_linear |
| `extrapolate_order` | number | no |  | asymptotic convergence order p > 0 for Richardson extrapolation; default 1.0 for 3D/sheet capacitance matrices, else 2 * order |
| `factor` | number | no |  | comparison/present element-size ratio, positive and not 1 (< 1 refines, > 1 coarsens), default 0.6 |
| `match` | text: `index`, `nearest` | no |  | how the solutions of the two meshes are paired: index (default) compares solution i with solution i; nearest pairs each mode with the nearest eigenvalue of the other mesh in the complex plane and lists the modes without a partner |
| `max_cells` | number | no |  | refuse, before any solve, a check whose comparison mesh has more cells than this (an adaptive 2D mesh has 6 to 8 times as many cells at factor 0.6) |
| `max_seconds` | number | no |  | refuse, before any solve, a check estimated (from the solves this machine has run; no refusal without an estimate) to take longer than this for the comparison mesh |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `quantities` | list of text | no |  | non-empty list of quantity names to check; only matching rows count toward converged, worst_relative, and worst_quantity. each item: quantity name or scalar name to include in the verdict (e.g. "eigenvalue[0].re", "scalar C_0_0 (F)", or "C_0_0 (F)") |
| `scale_sizes` | true or false | no |  | scale non-suppressed add_size features (h, max_deviation, min_layers) proportionally with factor on the comparison copy, default true |
| `tolerance` | number | no |  | largest accepted relative difference, e.g. 1e-3; omit to only report |
| `track` | object | no |  | follow mode(s) matching a target eigenvalue (f, n_b or n_eff) and/or physical mode-signature filters across both meshes |
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
