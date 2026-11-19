# `add_transform`

**Transform** · access `edit`: changes the project

Moves, rotates, mirrors or scales a feature; the steps apply in order. [build; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `input` | text or integer | yes |  | the feature to transform. the name of an earlier feature (or its id) |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `steps` | list of object | yes |  | list of steps in the serde form of TStep. each item: {"Translate": [dx, dy, dz]} \| {"Rotate": {axis, deg, about}} \| {"Mirror": {axis, at}} \| {"Scale": {factor, about}} |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_transform","args":{"input":"<input>","steps":[]}}
```
