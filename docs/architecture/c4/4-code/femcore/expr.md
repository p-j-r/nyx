# Expressions: code

C4 level 4 (code) of the component **Expressions** (`femcore.expr`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> A small safe expression evaluator for user input; a copy of geom's kept equal by a corpus test.

**Decisions:** [ADR 65](../../../adr/0065-closing-the-restructure-what-the-structure-replaced.md), [ADR 2](../../../adr/0002-geometry-as-named-feature-list.md)  
**Size:** 1 file, 286 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/expr.rs` | 286 | A small, safe expression evaluator for user input such as `293.15 + 50*exp(-((x-0.5)^2)/0.01)` (initial temperature fields, for example). Numbers, `+ - * / ^`, parentheses, the functions `sin cos tan exp log ln log10 sqrt abs sinh cosh t... |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Frame and sections](./section.md) | section | wafer_point |

## Used by

| Component | Modules | Items |
|---|---|---|
| [dfs: Heat kernel](../dfs/heat.md) | expr | eval_at |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Expr {
      <<struct>>
      +source String
      +parse(text) Result
      +eval(x, y, z) f64
  }
  class fns_expr_rs {
      <<module expr.rs>>
      +eval_at(text, coords, dim) Result
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `Expr` | `expr.rs:34` | fields: `source` |
| fn | `eval_at` | `expr.rs:256` | `pub fn eval_at(text: &str, coords: &[f64], dim: usize) -> Result<Vec<f64>, String>` —  Evaluates `text` at dof coordinates (flat, `dim` per point) in the wafer axes: `x`, `y`, `z`. On a 2D cross-section the mesh is the x-z |
| fn | `Expr::parse` | `expr.rs:199` | `pub fn parse(text: &str) -> Result<Expr, String>` |
| fn | `Expr::eval` | `expr.rs:212` | `pub fn eval(&self, x: f64, y: f64, z: f64) -> f64` |
