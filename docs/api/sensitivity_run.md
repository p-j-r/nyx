# `sensitivity_run`

**Compute discrete adjoint sensitivities** · access `expensive`: starts a mesh build or a solve, or writes files

Evaluates the scalar objective `J(p)` and its exact discrete adjoint, modal, or analytical Redheffer EME sensitivity gradient `dJ/dp` across named CAD parameters (`param_var`, via harmonic mesh-velocity extension), scalar material properties (`material_scalar`), EME segment lengths (`eme_segment_length`, via `dP_c/dL_c = -j beta_c P_c`), and Helmholtz-filtered SIMP/RAMP/optical cellwise densities (`topology_density`) without finite-difference solve loops (ADR 140, ADR 142).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `include_densities` | true or false | no |  | include projected physical densities in the reply when topology parameters are present, default true |
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
| `shape_step_rel` | number | no |  | optional relative CAD step size for harmonic boundary velocity extraction, default 1e-6 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"sensitivity_run","args":{"objective":{"kind":"eigenfrequency"},"params":[{"kind":"param_var"}]}}
```
