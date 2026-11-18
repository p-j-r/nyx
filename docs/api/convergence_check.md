# `convergence_check`

**Check mesh convergence** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the physics again on a finer mesh (element size times `factor`; the project is not changed) and compares every eigenvalue, value and scalar with the present result. With `tolerance` each row says pass or fail and `converged` summarises. With a `tolerance` the verdict is kept for this problem: the answer has `validated` (the result is then `checked` or stays `exploratory`) and the run record that solved it gets the verdict. Without a tolerance nothing is recorded. Use it before quoting a number: a mode index that moves by more than your tolerance is not converged.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `factor` | number | no |  | finer/present element size, between 0 and 1, default 0.6 |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `tolerance` | number | no |  | largest accepted relative difference, e.g. 1e-3; omit to only report |
| `wait` | true or false | no |  | wait for the answer, default true; false starts a job that job_list, job_status and job_cancel see and job_wait answers |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"convergence_check","args":{"physics":"dfs"}}
```
