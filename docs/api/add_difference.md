# `add_difference`

**Difference** · access `edit`: changes the project

The feature `a` with the features `b` cut away. [build; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `a` | text or integer | yes |  | the feature kept. the name of an earlier feature (or its id) |
| `b` | list of text or integer | yes |  | the features removed |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_difference","args":{"a":"<a>","b":[]}}
```
