# 90. Whole-array operations on results are engine commands over a typed expression tree evaluated in Rust

Date: 2026-11-17

## Status

Accepted as step S3 (plan item N1) of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3b). **The number 90 is the one the plan reserves for "per-element numerics are formulas evaluated in Rust" and is provisional**: take the next free ADR number at merge and renumber this file, the index row and the code comments. This record covers level N1 (whole-array operations) and the tree that levels N2 to N4 build on; N2 (text formulas, material values, `where`), N3 (tracing script functions) and N4 (a code generator) amend it when they land.

## Context

A researcher or an agent wants numbers computed from a field: the integral of `|E|^2` over the core, the confinement factor, the largest stress, the field along a line, the overlap of two modes. Before this step the only way was `result_field`, which returns up to 2000 sampled node-averaged values, and a loop of arithmetic written by the caller. In an interpreted script such a loop is 10 to 100 times slower than Rust and, worse, it reads values that were averaged over the cells and sampled, so the number it gets is not the integral of the field the solver made.

## Decision

* **Eight commands, one implementation** (`nyx_engine::reduce`): `field`, `field_op`, `integrate`, `extremum`, `probe`, `line_cut`, `overlap` and `confinement`. They are ordinary engine commands (read access), so the MCP server lists them as tools, a script (S4) calls them, and the window's Reduce form in the Results tab calls the very same commands (`App::reduce_now`; the test compares the line it shows with the command's number).
* **A handle is a value, not an id.** `field` returns `{"field": "Ex", "physics": "ems", "solution": 0}`; arithmetic builds `{"op": "mul", "args": [a, b]}` (`add sub mul div neg abs abs2 re im conj sqrt pow`); numbers are `{"num": x}` or `{"num": [re, im]}`; a bare string is a field of the request's `physics` and `solution`. The engine keeps no handle table, so a handle has no lifetime and cannot go stale in the engine: it is read against the stored result each time and the answer says `stale` and `age_s` (as `result_field` does). The tree is `femcore::fieldexpr::Tree` (complex values, leaves are numbered slots, every operation written out so the order of the floating-point operations is fixed by that file).
* **A leaf reads the field the solver made.** A name is first looked up among the per-cell fields of the solution (`cell_fields`: `Ex`, `uz`, ..., complex when the solve is, not averaged over cells), then among the node-averaged flat fields of `result_info` (`T`, `|E|`: real, one value per dof). Both are Lagrange P1 or P2 interpolants on the cells of the mesh the result was solved on. **A stored result keeps that mesh** (`Stored.model`), so a result is read on its own cells whatever the project's mesh is now (a changed element size makes the result stale, not unreadable).
* **The rule (what N2, N3 and N4 must reproduce bit for bit).** Integrals use the degree-4 symmetric rule of `femcore::fem::quad_rule` (6 points on a triangle, 14 on a tetrahedron), exact for the product of two P2 fields; an integrand of higher degree carries the quadrature error. At each point every leaf is `sum_a N_a(q) v_a` (local node order, starting from 0.0) and the tree is evaluated on those values. A cell adds `vol * sum_q (w_q * f_q)`; cells are summed in blocks of 2048 consecutive cell numbers in cell order and the blocks in block order, so the sum has the same bits for any thread count (`tests/reduce.rs` runs it on 1, 3 and the pool's threads). The reference loops in `tests/reduce.rs` are written on the raw arrays of the result and agree **to the bit** with the commands; the closed forms (polynomials of degree 2 integrated by Gauss-Legendre over boxes) agree to 1e-9.
* **Measure.** The plane measure (area in 2D, volume in 3D) is the default. `revolution: true` multiplies the integrand at each quadrature point by `2 pi r` (r is mesh axis 0, the mesh must be 2D) for bodies of revolution. It is an explicit argument, echoed in every answer, because the engine does not know which physics is axisymmetric; when a result's study name says `axisym` and the argument is absent, the answer carries a note.
* **Maxima are nodal.** `extremum` takes the largest or smallest real value at the nodes of the cells of a region (vertices, and edge midpoints of a P2 mesh), not over the interpolant; it says so, and returns the cell, the local node and the point. A complex expression has no order and is refused with the advice to wrap it in `abs`, `abs2`, `re` or `im`.
* **`probe` and `line_cut`** find the cell of a point with a uniform grid of the cells' boxes and interpolate in that cell. Points are in mesh axes (a 2D mesh: x then z of the wafer; 3D: x, y, z). A point outside the mesh answers `inside: false` and a null value, not an error; a point on the face between two cells is read in the cell with the lowest number (the fields are not averaged, so the two cells can differ by the discretisation error).
* **`overlap` and `confinement` have no material weight.** `overlap` returns the complex integral of `sum_k conj(a_k) b_k` over a region, the two norms and the coefficient `|overlap|^2 / (norm_a norm_b)` (0 to 1, Cauchy-Schwarz); `a` and `b` are a field or a list of components. `confinement` returns the integral of `sum_k |c_k|^2` over a region divided by the same over the whole mesh. The energy or power fraction needs the permittivity at each point, which is N2 (material values in formulas); until then the commands say what they computed in `meaning`.
* **Where the tree lives.** The plan says N1's handle arithmetic builds the formula tree of N2, extending `femcore::expr`. That evaluator is real-valued and has no fields; making it complex and giving it field leaves is the substance of N2. So S3 adds the typed tree `femcore::fieldexpr` (complex, field leaves) and N2 adds a text parser that produces the same tree and moves `expr::Expr` (initial values from `x y z`) onto it. Until then there are two small evaluators in `femcore`; the gate for N2 is that a formula and a handle give the same bits.

## Not in this decision

* Material values (`eps`, `rho`) at each point, text formulas, `where`, vectors and dot products (N2); tracing a script function into the tree (N3); a code generator (N4).
* Convergence of an integral under refinement (a reduction reads one result on one mesh).
* Reductions in a worker process: they read the result the parent holds and take about as long as reading it.
* Maxima of the interpolant between nodes; surface and line integrals along boundaries (`integrate` is over cells).

## Consequences

* `engine_commands` gains eight commands and `result_info` gains `cell_fields` (additive); `mcp_tools` gains eight tools (additive). Provisional versions 1.12 and 2.8.
* `Stored` holds the mesh model of its result, so a result keeps its mesh alive until it is replaced.
* `nyx check` treats the eight commands as readers of a result.

## Code

* `crates/femcore/src/fieldexpr.rs` (the tree and its evaluator), `crates/nyx-engine/src/reduce.rs` (the commands), `crates/nyx-engine/src/engine.rs` (`Stored.model`), `crates/app/src/tabs.rs` (the Reduce form).
* Tests: `crates/nyx-engine/tests/reduce.rs` (closed forms, loops to the bit, thread counts, a worker result, a real mode and a real thermal field), `crates/app/src/app.rs` (the form shows the engine's numbers).
