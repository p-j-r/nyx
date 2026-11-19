# `add_pattern`

**Pattern** · access `edit`: changes the project

n copies of a feature, in a line or around an axis. [build; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `input` | text or integer | yes |  | the feature to repeat. the name of an earlier feature (or its id) |
| `kind` | object | yes |  | Linear { step } or Polar { axis, center, total_deg }. {"Linear": {"step": [dx, dy, dz]}} \| {"Polar": {axis, center, total_deg}} |
| `n` | number or text | yes |  | number of copies. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_pattern","args":{"input":"<input>","kind":{},"n":0.0}}
```
