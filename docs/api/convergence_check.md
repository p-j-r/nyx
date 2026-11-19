# `convergence_check`

**Check mesh convergence** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the physics again on a finer mesh (element size times `factor`; the project is not changed) and compares every eigenvalue, value and scalar with the present result. With `tolerance` each row says pass or fail and `converged` summarises. With a `tolerance` the verdict is kept for this problem: the answer has `validated` (the result is then `checked` or stays `exploratory`) and the run record that solved it gets the verdict. Without a tolerance nothing is recorded. Use it before quoting a number: a mode index that moves by more than your tolerance is not converged. When the order of the modes changes between the meshes (complex modes of an open resonator), `match` "nearest" pairs each mode with the nearest eigenvalue of the finer mesh and `track` follows one mode; the answer then lists the pairs and the modes without a partner (the present mesh is then always solved again, never taken from an earlier request). Only the global element size and its rules are made finer: a size set with add_size stays as it is, so refine those by hand and solve again to test them.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `factor` | number | no |  | finer/present element size, between 0 and 1, default 0.6 |
| `match` | text: `index`, `nearest` | no |  | how the solutions of the two meshes are paired: index (default) compares solution i with solution i; nearest pairs each mode with the nearest eigenvalue of the finer mesh in the complex plane and lists the modes without a partner |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `tolerance` | number | no |  | largest accepted relative difference, e.g. 1e-3; omit to only report |
| `track` | object | no |  | follow one mode: { f } is a value (Hz for a vibration study); the mode whose eigenvalue is nearest it is compared with the nearest eigenvalue of the finer mesh |
| `track.f` | number | yes |  |  |
| `wait` | true or false | no |  | wait for the answer, default true; false starts a job that job_list, job_status and job_cancel see and job_wait answers |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"convergence_check","args":{"physics":"dfs"}}
```
