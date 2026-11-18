# Expressions: code

C4 level 4 (code) of the component **Expressions** (`femcore.expr`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> A small safe expression evaluator for user input (a copy of geom's kept equal by a corpus test), and the typed tree over complex field values that the reductions of results evaluate, and the text formulas parsed into it (ADR 90).

**Decisions:** [ADR 65](../../../adr/0065-closing-the-restructure-what-the-structure-replaced.md), [ADR 2](../../../adr/0002-geometry-as-named-feature-list.md), [ADR 90](../../../adr/0090-whole-array-operations-are-engine-commands-over-a-typed-tree.md)  
**Size:** 3 files, 1562 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/expr.rs` | 286 | A small, safe expression evaluator for user input such as `293.15 + 50*exp(-((x-0.5)^2)/0.01)` (initial temperature fields, for example). Numbers, `+ - * / ^`, parentheses, the functions `sin cos tan exp log ln log10 sqrt abs sinh cosh t... |
| `crates/femcore/src/fieldexpr.rs` | 436 | A typed expression tree over complex field values, evaluated at one point at a time (productivity plan N1 and N2, ADR 90). |
| `crates/femcore/src/formula.rs` | 840 | Text formulas over fields, coordinates and material values (productivity plan N2, ADR 90): `abs(Ex)^2 * eps_r`, `where(x > 0, T, 0)`, `dot(E, conj(E))`. The text is parsed once into a `fieldexpr::Tree`, the same tree that handle arithmet... |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Matrices and vectors](./lacore.md) | linalg::csr::assembly | C |
| [Frame and sections](./section.md) | section | wafer_point |

## Used by

| Component | Modules | Items |
|---|---|---|
| [dfs: Heat kernel](../dfs/heat.md) | expr | eval_at |
| [engine: Commands](../engine/commands.md) | fieldexpr, formula | Fun, MAX_POWER, Named, Names, Op, Tree |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace expr_rs {
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
  }
  namespace fieldexpr_rs {
    class Op {
        <<enumeration>>
        Add
        Sub
        Mul
        Div
        And
        Or
    }
    class Cmp {
        <<enumeration>>
        Lt
        Le
        Gt
        Ge
        Eq
        Ne
    }
    class Fun {
        <<enumeration>>
        Abs
        Abs2
        Re
        Im
        Conj
        Sqrt
        Exp
        Sin
        Cos
        Arg
        RealSqrt
        Log
        more_8_variants
        +name() str
    }
    class Fun2 {
        <<enumeration>>
        Min
        Max
        Pow
        Atan2
        +name() str
    }
    class Tree {
        <<enumeration>>
        Num
        Leaf
        Neg
        Bin
        Call
        Call2
        Pow
        Cmp
        Where
        +num(x) Tree
        +leaves() usize
        +eval(leaf) C
        +is_real_valued(leaf_is_complex) bool
        +check(leaf_is_complex) Result
    }
    class fns_fieldexpr_rs {
        <<module fieldexpr.rs>>
        +MAX_POWER const
    }
  }
  namespace formula_rs {
    class Named {
        <<enumeration>>
        Scalar
        Vector
    }
    class Names {
        <<interface>>
        +name(name) Result
        +material(name) Result
        +region(name) Result
        +known() String
    }
    class fns_formula_rs {
        <<module formula.rs>>
        +parse(text, names) Result
    }
  }
  Tree --> Op
  Tree --> Fun
  Tree --> Fun2
  Tree --> Cmp
  Named --> Tree
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `Expr` | `expr.rs:34` | fields: `source` |
| fn | `eval_at` | `expr.rs:256` | `pub fn eval_at(text: &str, coords: &[f64], dim: usize) -> Result<Vec<f64>, String>` —  Evaluates `text` at dof coordinates (flat, `dim` per point) in the wafer axes: `x`, `y`, `z`. On a 2D cross-section the mesh is the x-z |
| fn | `Expr::parse` | `expr.rs:199` | `pub fn parse(text: &str) -> Result<Expr, String>` |
| fn | `Expr::eval` | `expr.rs:212` | `pub fn eval(&self, x: f64, y: f64, z: f64) -> f64` |
| enum | `Op` | `fieldexpr.rs:18` | variants: `Add`, `Sub`, `Mul`, `Div`, `And`, `Or` —  Binary operations. `And` and `Or` read the truth of the real parts (nonzero is true) and answer 1 or 0. |
| enum | `Cmp` | `fieldexpr.rs:29` | variants: `Lt`, `Le`, `Gt`, `Ge`, `Eq`, `Ne` —  Comparisons of two real values; the answer is 1 or 0 (false when either side is NaN). |
| enum | `Fun` | `fieldexpr.rs:40` | variants: `Abs`, `Abs2`, `Re`, `Im`, `Conj`, `Sqrt`, `Exp`, `Sin`, `Cos`, `Arg`, `RealSqrt`, `Log`, ... —  Functions of one value. |
| enum | `Fun2` | `fieldexpr.rs:112` | variants: `Min`, `Max`, `Pow`, `Atan2` —  Functions of two real values. |
| enum | `Tree` | `fieldexpr.rs:133` | variants: `Num`, `Leaf`, `Neg`, `Bin`, `Call`, `Call2`, `Pow`, `Cmp`, `Where` |
| const | `MAX_POWER` | `fieldexpr.rs:149` | `u32` —  Largest integer power a tree may carry. |
| fn | `Fun::name` | `fieldexpr.rs:76` | `pub fn name(self) -> &'static str` —  The name a formula uses. |
| fn | `Fun2::name` | `fieldexpr.rs:122` | `pub fn name(self) -> &'static str` |
| fn | `Tree::num` | `fieldexpr.rs:180` | `pub fn num(x: f64) -> Tree` |
| fn | `Tree::leaves` | `fieldexpr.rs:185` | `pub fn leaves(&self) -> usize` —  The number of leaf slots the tree reads (one more than the largest slot). |
| fn | `Tree::eval` | `fieldexpr.rs:196` | `pub fn eval(&self, leaf: &[C]) -> C` —  The value at one point; `leaf[i]` holds slot `i`. |
| fn | `Tree::is_real_valued` | `fieldexpr.rs:283` | `pub fn is_real_valued(&self, leaf_is_complex: &[bool]) -> bool` —  Whether the value is real for these leaves: `true` when the tree ends in a function that answers a real value, or reads no imaginary |
| fn | `Tree::check` | `fieldexpr.rs:300` | `pub fn check(&self, leaf_is_complex: &[bool]) -> Result<(), String>` —  Refuses what has no meaning for the leaves' values: an order, a comparison, a logical operation, `log`, `tan`, a non-integer power |
| enum | `Named` | `formula.rs:26` | variants: `Scalar`, `Vector` —  What a name stands for. |
| trait | `Names` | `formula.rs:33` | methods: `name`, `material`, `region`, `known` —  The caller's side of a formula: what its names mean. |
| fn | `parse` | `formula.rs:78` | `pub fn parse(text: &str, names: &mut dyn Names) -> Result<Tree, String>` —  Parses `text` into a scalar tree; `names` says what its names are. |
