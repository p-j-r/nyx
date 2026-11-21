# `get_guide`

**Read a nyx guide** · access `read`: reads the project or the results

Read a short guide. Call with topic nyx-instructions before using other nyx tools. Other topics: geometry, optical-modes (reading mode tables), cqed (superconducting qubit capacitance, inductance, EPR, Hamiltonian and TLS participation), troubleshooting (errors).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `topic` | text: `nyx-instructions`, `geometry`, `optical-modes`, `cqed`, `troubleshooting` | no |  | which card; omit for nyx-instructions |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"get_guide","args":{}}
```
