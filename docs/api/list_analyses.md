# `list_analyses`

**List the analyses** · access `read`: reads the project or the results

The analyses of results (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator), which solved physics each applies to and, per physics, whether its result is `stale`.

## Arguments

None.

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"list_analyses","args":{}}
```
