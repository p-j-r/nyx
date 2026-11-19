# `tag_feature`

**Tag a feature** · access `edit`: changes the project

Adds (or removes) a tag of a feature.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | feature name |
| `remove` | true or false | no |  | true removes the tag |
| `tag` | text | yes |  | tag |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"tag_feature","args":{"name":"<name>","tag":"<tag>"}}
```
