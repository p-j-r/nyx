# `optimize_run`

**Run gradient-based shape, material, or topology optimization** · access `expensive`: starts a mesh build or a solve, or writes files

Optimizes CAD geometry parameters, scalar material properties, EME propagation segment lengths (`eme_segment_length`), or unstructured Helmholtz-filtered cellwise topology densities using `lbfgs`, `ncg` (Polak-Ribiere+), or constrained `mma` (method of moving asymptotes) driven by exact discrete adjoint / analytical Redheffer gradients, optional Heaviside `beta` continuation stages, and modal eigenvector tracking (ADR 140, ADR 142). Stores the final solved result (attaching `topology_density` when present) and optionally writes a JSON checkpoint to `checkpoint_path`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `algorithm` | text: `lbfgs`, `ncg`, `mma` | no |  | gradient optimizer: lbfgs (default), ncg (Polak-Ribiere+ nonlinear conjugate gradients), or mma (method of moving asymptotes) |
| `apply` | true or false | no |  | apply optimized parameters and store the final solved result on the active project, default true |
| `checkpoint_path` | text | no |  | optional JSON file path for writing optimization history and final design checkpoint |
| `constraints` | list of object | no |  | inequality constraints c_i(x) <= 0 enforced by mma. each item: one optimization constraint |
| `constraints[].kind` | text: `volume_fraction` | yes |  | constraint kind |
| `constraints[].max_fraction` | number | yes |  | upper bound on physical volume fraction in (0, 1] |
| `constraints[].min_fraction` | number | no |  | optional lower bound on physical volume fraction in [0, max_fraction] |
| `constraints[].param_index` | integer | no |  | 0-based index of the topology_density entry in params (default 0) |
| `continuation_schedule` | list of object | no |  | optional multi-stage Heaviside beta continuation schedule [{ beta, max_iters }]. each item: one Heaviside projection sharpness stage |
| `continuation_schedule[].beta` | number | yes |  | Heaviside sharpness beta (>= 0) for this stage |
| `continuation_schedule[].max_iters` | integer | yes |  | maximum optimizer iterations in this stage (>= 1) |
| `ftol` | number | no |  | relative objective change threshold (> 0) for lbfgs and ncg, default 1e-12 |
| `gtol` | number | no |  | relative gradient norm convergence threshold (> 0) for lbfgs and ncg, default 1e-6 |
| `max_iters` | integer | no |  | maximum optimizer iterations per stage (>= 1), default 50 |
| `move_limit` | number | no |  | fractional per-iteration box move limit (> 0) for mma, default 0.2 |
| `objective` | object | yes |  | scalar objective differentiated by discrete adjoint or modal sensitivity (`eigenfrequency`, `effective_index`, `quality_factor`, `transmission`, `splitter_balance`, `capacitance`, `compliance`, `target_temperature`, or `field_power`) |
| `objective.domain` | text | no |  | region name or pattern for target_temperature or field_power |
| `objective.freq_index` | integer | no |  | optional 0-based driven frequency index (omit to average across all sweep frequencies) |
| `objective.in_port` | integer | no |  | 0-based excited input port index for transmission, splitter_balance, or field_power (default 0) |
| `objective.inverse` | true or false | no |  | for quality_factor: differentiate loss tangent 1/Q instead of Q, default false |
| `objective.kind` | text: `eigenfrequency`, `effective_index`, `quality_factor`, `transmission`, `splitter_balance`, `capacitance`, `compliance`, `target_temperature`, `field_power` | yes |  | objective family |
| `objective.mode` | integer | no |  | 0-based eigenmode index (default 0) |
| `objective.out_port` | integer | no |  | 0-based output port index for transmission (default 1) |
| `objective.out_ports` | list of integer | no |  | 0-based output port indices for splitter_balance |
| `objective.p_norm` | number | no |  | optional smooth L_p peak norm exponent (>= 1) for target_temperature |
| `objective.scale_hz` | number | no |  | optional normalization scale in Hz (> 0) for target_hz mismatch |
| `objective.target` | number | no |  | optional target transmission or compliance value (switches to squared mismatch) |
| `objective.target_f` | number | no |  | optional target capacitance in farads (switches to squared mismatch) |
| `objective.target_hz` | number | no |  | optional target resonant frequency in Hz (switches objective to normalized squared mismatch) |
| `objective.target_k` | number | no |  | optional target temperature in kelvin for squared tracking error |
| `objective.target_neff` | number | no |  | optional target effective or Bloch index (switches objective to 0.5 * (neff - target_neff)^2) |
| `objective.targets` | list of number | no |  | target power fractions per output port in splitter_balance |
| `objective.term_i` | integer | no |  | 0-based row terminal index for capacitance C_{ij} (default 0) |
| `objective.term_j` | integer | no |  | 0-based column terminal index for capacitance C_{ij} (default 0) |
| `optimize` | object | no |  | optional OptimizeSpec object overriding default optimizer tolerances, constraints, or continuation stages |
| `params` | list of object | yes |  | design parameters differentiated or optimized (`param_var`, `material_scalar`, `eme_segment_length`, or `topology_density`). each item: one scalar CAD, scalar material, EME segment-length, or cellwise topology-density design parameter |
| `params[].beta` | number | no |  | Heaviside tanh projection sharpness beta (>= 0, default 1.0) |
| `params[].densities` | list of number | no |  | optional explicit per-cell raw densities in [0, 1] |
| `params[].domain` | text | no |  | design region name or pattern (for topology_density) |
| `params[].eta` | number | no |  | Heaviside projection threshold eta in [0, 1] (default 0.5) |
| `params[].filter_radius` | number | no |  | unstructured Helmholtz PDE filter radius l_f in metres (>= 0, default 0) |
| `params[].initial_density` | number | no |  | uniform initial raw density in [0, 1] when densities is omitted |
| `params[].kind` | text: `param_var`, `material_scalar`, `eme_segment_length`, `topology_density` | yes |  | design parameter category |
| `params[].max` | number | no |  | optional upper box bound |
| `params[].min` | number | no |  | optional lower box bound |
| `params[].name` | text | no |  | geometry parameter name in FeatureList::params (for param_var) |
| `params[].plasmon_damping` | number | no |  | intermediate extinction floor for plasmon-safe metallic optical interpolation (default 0.0) |
| `params[].property` | text | no |  | scalar material property name such as eps_r, k, E, nu, rho, Q (for material_scalar or topology_density) |
| `params[].ramp_q` | number | no |  | optional RAMP concavity parameter q (>= 0); overrides simp_penalty when given |
| `params[].region` | text | no |  | material region name (for material_scalar) |
| `params[].segment` | integer | no |  | 0-based EME segment index for eme_segment_length (default 0) |
| `params[].simp_penalty` | number | no |  | SIMP power-law exponent p (>= 1, default 3.0) |
| `params[].val_max` | number | no |  | solid endpoint scalar/real property value at rho_bar = 1 (default 1.0) |
| `params[].val_max_im` | number | no |  | solid imaginary permittivity endpoint at rho_bar = 1 (default 0.0) |
| `params[].val_min` | number | no |  | void endpoint scalar/real property value at rho_bar = 0 (default 1e-3) |
| `params[].val_min_im` | number | no |  | void imaginary permittivity endpoint at rho_bar = 0 (default 0.0) |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `remesh` | true or false | no |  | rebuild the CAD mesh at every shape evaluation instead of smooth harmonic mesh morphing, default false |
| `sense` | text: `minimize`, `maximize` | no |  | optimization direction: minimize (default) or maximize |
| `shape_step_rel` | number | no |  | optional relative CAD step size for harmonic boundary velocity extraction, default 1e-6 |
| `xtol` | number | no |  | maximum design variable step threshold (> 0) for mma, default 1e-6 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"optimize_run","args":{"objective":{"kind":"eigenfrequency"},"params":[{"kind":"param_var"}]}}
```
