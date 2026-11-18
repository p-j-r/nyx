# `field`

**A handle to a field** · access `read`: reads the project or the results

A handle to a field of a solution, to compute with (`integrate`, `extremum`, `probe`, `line_cut`, `overlap`, `confinement`, `field_op`). A handle is a small value, not the field: it names the physics, the solution and the field, and is read from the stored result each time it is used. A name is a per-cell field of the solution (`Ex`, `uz`, ...: complex when the solve is) or else a flat field of `result_info` (`T`, `|E|`: real). A vector name (`E` when `Ex` and `Ez` exist) answers one handle per component.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | field name, e.g. Ex, or a vector name such as E |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `solution` | integer | no |  | solution number from 0, default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"field","args":{"name":"<name>"}}
```
