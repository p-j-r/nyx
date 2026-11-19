# `set_problem`

**Replace the problem of a physics** · access `edit`: changes the project

Replaces the whole problem (study, order, materials, conditions, options) of a physics by one given as JSON, after checking it against the descriptors. Adds the physics if the project has no problem for it yet.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `problem` | object | yes |  | a problem: { physics, order, study, materials, bcs, params } |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_problem","args":{"physics":"dfs","problem":{}}}
```
