# `list_remote_materials`

**List downloadable remote materials** · access `read`: reads the project or the results

Curated presets from open material databases (`refractiveindex.info`, `AFLOW`, `NumBAT`) that `fetch_material` can download and convert on demand.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `domain` | text: `all`, `photonics`, `phononics`, `superconducting` | no |  | physics domain filter, default all |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"list_remote_materials","args":{}}
```
