# `add_intersect`

**Intersection** · access `edit`: changes the project

The common part of features. [build; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `ids` | list of text or integer | yes |  | the features to intersect |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_intersect","args":{"ids":[]}}
```
