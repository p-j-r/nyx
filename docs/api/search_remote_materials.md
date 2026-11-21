# `search_remote_materials`

**Search curated and live open material databases** · access `read`: reads the project or the results

Searches the curated remote material presets by keyword, formula, domain or source database, and optionally queries live AFLOW and Materials Project HTTPS endpoints.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `api_key` | text | no |  | optional Materials Project API key for live mp queries (falls back to MP_API_KEY environment variable) |
| `domain` | text: `all`, `photonics`, `phononics`, `superconducting` | no |  | physics domain filter, default all |
| `live` | true or false | no |  | also query live AFLOW and Materials Project HTTPS endpoints, default false |
| `query` | text | yes |  | search keyword, chemical formula or material identifier (e.g. Si, AlN, LiNbO3, mp-149) |
| `source` | text: `all`, `rii`, `aflow`, `mp`, `numbat` | no |  | database filter, default all |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"search_remote_materials","args":{"query":"<query>"}}
```
