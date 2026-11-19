# `poll_jobs`

**Apply finished jobs** · access `read`: reads the project or the results

Applies every job that has finished since the last call and returns what each came to.

## Arguments

None.

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"poll_jobs","args":{}}
```
