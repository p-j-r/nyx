# `set_param`

**Set a model parameter** · access `edit`: changes the project

Sets (or creates) a named parameter that expressions in the features can use. A number, or an expression over other parameters.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | parameter name |
| `value` | number or text | yes |  | a number in metres (or the unit the features give it), or an expression over other parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_param","args":{"name":"<name>","value":0.0}}
```
