# `add_tag`

**Tag** · access `edit`: changes the project

Gives a selection a name (a boundary group) to refer to later. [naming; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | a name for the feature; default the kind and a number (`rect1`) |
| `sel` | text | yes |  | what to tag. selector text, e.g. `outer & plane(x = 0)` |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_tag","args":{"name":"<name>","sel":"<sel>"}}
```
