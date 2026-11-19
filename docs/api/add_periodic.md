# `add_periodic`

**Periodic pair** · access `edit`: changes the project

Declares two boundary sets as periodic images of each other. [boundary; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `sel_a` | text | yes |  | first boundary set. selector text, e.g. `outer & plane(x = 0)` |
| `sel_b` | text | yes |  | second boundary set. selector text, e.g. `outer & plane(x = 0)` |
| `translation` | list of 3 (number or text) | yes |  | the translation from A to B. each item: a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_periodic","args":{"sel_a":"<sel_a>","sel_b":"<sel_b>","translation":[0.0,0.0,0.0]}}
```
