//! Validation harness: a small, independent scalar Helmholtz finite-element solver with PML.
//!
//! It exists to test the meshes this crate produces without depending on nyx: P1 on triangles and
//! Q1 on quadrilaterals, complex coordinate stretching in the PML regions described by
//! [`crate::meshers::adapt::pml::PmlMesh`], and a direct banded solver. It is deliberately simple and not
//! tuned; it is not a replacement for the nyx solvers.

pub mod band;
pub mod cplx;
pub mod helmholtz;
