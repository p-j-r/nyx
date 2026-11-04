# solid3d: stuffed tetrahedral meshes of solids (kept on purpose)

`mesher3d.rs` (isosurface stuffing of sphere, box, cylinder and extruded-outline bodies, with local refinement) and
`refine3d.rs` (Rivara longest-edge bisection, used only by `mesher3d`). They were `femcore::meshing::{mesher3d, refine3d}`
until phase 2 of the restructure (R05), where they moved here; paul asked to keep the code for now because it may have
been meant for the meshing path and never properly used: no library, app or MCP code calls them today.

* What they do: `mesher3d` is an isosurface-stuffing tetrahedral mesher for primitive solids, `refine3d` is Rivara refinement.
* Used by: `tests/emw3d.rs` (PEC sphere against the closed form) and `examples/cavity_bench.rs` (baseline benchmark),
  both through `#[path = "common/solid3d/mod.rs"] mod solid3d;`.
* Their own 15 tests (`tests/mesher3d.rs`) and `examples/mesh3d_bench.rs` were deleted in commit `1de99a33`; `jj file show -r 1de99a33~ crates/femcore/tests/mesher3d.rs` brings them back.
* To make them library code again: move the two files back to `crates/femcore/src/meshing/`, add `pub mod mesher3d; pub mod refine3d;` to `meshing/mod.rs`, replace `femcore::` by `crate::` in them.
