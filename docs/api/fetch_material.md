# `fetch_material`

**Download a material from an open database** · access `expensive`: starts a mesh build or a solve, or writes files

Downloads a preset or database entry over HTTPS from `refractiveindex.info`, `AFLOW`, `Materials Project` or `NumBAT`, converts it into validated TOML, and optionally adds it to the project or writes `<name>.toml` into `folder`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `add_to_project` | true or false | no |  | add the downloaded material to the active project (default true) |
| `api_key` | text | no |  | optional Materials Project API key (falls back to MP_API_KEY environment variable; passed via curl stdin, never stored) |
| `description` | text | no |  | optional one-line comparative summary written into the material description field |
| `extends` | text | no |  | optional base material name in the library to inherit unprovided properties from |
| `folder` | text | no |  | optional destination folder to write <name>.toml into |
| `name` | text | no |  | material name in TOML; defaults to the preset id or derived entry stem |
| `source` | text | yes |  | curated preset id (from `list_remote_materials`), `catalog`, `rii`, `aflow`, `mp` or `numbat` |
| `target` | text | no |  | preset id, database entry key, Materials Project id (e.g. mp-149), AFLOW query, NumBAT stem or https URL; omit when `source` is a preset id |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"fetch_material","args":{"source":"<source>"}}
```
