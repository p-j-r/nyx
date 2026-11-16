# `set_geometry`

**Replace the geometry** · access `edit`: changes the project

Replaces the whole feature list by one given as JSON (what the window's geometry editor writes); it is checked first.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `geometry` | object | yes |  | a feature list: { dim, params, features } |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_geometry","args":{"geometry":{}}}
```
