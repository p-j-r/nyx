# `set_coupled_loop`

**Configure two-way Anderson coupled loop** · access `edit`: changes the project

Configures or disables an explicit opt-in two-way Anderson-accelerated self-consistent coupled loop (`femcore::linalg::snes::anderson`, default history depth `m = 5`) for photothermal (`ems` <-> `dfs` `heat_stationary`) or electrothermal (`dfs` `conduction` `auxiliary` <-> `dfs` `heat_stationary`) convergence (`||T^{k+1} - T^k||_inf / max(1, ||T^k||_inf) <= rtol`). Pass `enabled: false` to clear.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `anderson_m` | integer | no |  | Anderson type-II history depth m (0..=64; 0 is Picard), default 5 |
| `auxiliary` | object or null | no |  | optional auxiliary source ProblemSpec (e.g. dfs conduction problem for intra-physics electrothermal V <-> T coupling); null clears auxiliary |
| `enabled` | true or false | no |  | when false, clears coupled_loop so solve_run enforces acyclic one-way chains; default true |
| `max_iters` | integer | no |  | maximum outer Anderson iterations (>= 1), default 25 |
| `rtol` | number | no |  | relative L_inf temperature convergence tolerance (> 0), default 1e-6 |
| `source_scale` | number | no |  | non-negative multiplier applied to exported heat_source before thermal assembly (e.g. input optical power in W), default 1.0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_coupled_loop","args":{}}
```
