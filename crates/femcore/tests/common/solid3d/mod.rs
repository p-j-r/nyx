//! Test fixture: tetrahedral meshes of solids built from primitives (sphere, box, cylinder), used by the sphere-cavity
//! test and benchmark. Not part of the library (R05, phase 2): production 3D meshes come from `geom` models and
//! `meshing::mesher::extrude_y`. Include with `#[path = "common/solid3d/mod.rs"] mod solid3d;`.
#![allow(dead_code)]
pub mod mesher3d;
pub mod refine3d;
