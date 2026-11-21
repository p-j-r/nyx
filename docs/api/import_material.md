# `import_material`

**Import an external material record** · access `edit`: changes the project

Converts raw payload text from an open database (`rii` YAML, `aflow` JSON, `mp` JSON or `numbat` JSONC) into validated TOML, optionally extending a library material and adding it to the project or writing `<name>.toml` into `folder`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `add_to_project` | true or false | no |  | add the converted material to the active project (default true) |
| `content` | text | yes |  | raw YAML, JSON or JSONC payload text from the external database |
| `description` | text | no |  | optional one-line comparative summary written into the material description field |
| `extends` | text | no |  | optional base material name in the library to inherit unprovided properties from |
| `folder` | text | no |  | optional destination folder to write <name>.toml into |
| `format` | text: `rii`, `aflow`, `mp`, `numbat` | yes |  | source database format: rii, aflow, mp or numbat |
| `name` | text | yes |  | material name without slashes or commas |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"import_material","args":{"content":"<content>","format":"rii","name":"<name>"}}
```
