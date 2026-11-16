# `suppress_feature`

**Suppress or restore a feature** · access `edit`: changes the project

A suppressed feature stays in the list and is skipped.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | feature name |
| `suppressed` | true or false | no |  | true (default) to suppress, false to restore |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"suppress_feature","args":{"name":"<name>"}}
```
