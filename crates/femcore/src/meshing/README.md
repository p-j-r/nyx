# Moved: mesher3d and refine3d

The 3D solid mesher (`mesher3d`) and its refinement (`refine3d`) are no longer modules of `femcore::meshing`. They are kept, unchanged apart from paths, as a test fixture in `crates/femcore/tests/common/solid3d/` (see the README there). Production 3D meshes come from geom models through `mesher::extrude_y`.
