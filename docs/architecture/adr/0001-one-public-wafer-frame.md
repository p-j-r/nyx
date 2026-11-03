# 1. One public wafer frame, propagation along +y

Date: 2026-11-02

## Status

Accepted (recorded retrospectively on 2026-11-03 from the history and from `docs/FRAME_UNIFICATION_STATUS.md` and `docs/KERNEL_PORT_STATUS.md`).

## Context

Before this decision nyx had several coordinate conventions at once: 2D legacy models in "model axes", a `Frame` field on `Materials`, an `Axes` toggle in the window, wrong 2D labels (y where z was meant), and a positive out-of-plane wave number `k` that meant a wave towards the viewer (-y). A frame adapter translated between the public frame and the kernels' frame, and every new field name had to be placed by a 22-string table.

## Decision

There is one public frame, the **wafer frame**: x lateral, y the direction of propagation, z the surface normal, right-handed. A 2D cross-section is the x-z plane with y into the page. A wave along +y is `exp(-j beta y)` and has positive `k` (`PROP_SIGN = +1`). The time convention is `e^{j omega t}`.

* `femcore::section` is the one place that knows axis numbers: `plane_tensor`, `local_tensor`, `LOCAL_AXES`, `wafer_point`, `PLANE_VOIGT`, field names. Nobody else permutes axes.
* Solid kernels use the Voigt order (x, z, y); mode solvers use the local triple (lateral, normal, propagation), which is left-handed and only reverses the reading of the polarization sense.
* The ring (axisymmetric) kernels use their local wafer frame (r, theta, z) = (x, y, z); `m > 0` is a wave along +theta = +y.
* 2D PML is written along x and z. A y layer on a 2D mesh is refused, with no alias, because y means y in 3D files.
* `Frame`, `Materials::frame`, the adapter, `rename_wafer_fields` and the legacy front-ends were deleted rather than kept as shims: nothing was released.

## Consequences

* All kernels name result fields in wafer axes straight from the kernel (`Ex`, `Ez`, `Ey`, `ux`, `uz`, `uy`, `sxz`, `syz`); a 2D field the table cannot place fails a debug assertion in `solve`.
* Positive `k` changed meaning (cosmetic for elastic frequencies, which are even in `k`); neff did not change, only the sign of `Ey`.
* geom is wafer-native too (`crates/geom/tests/wafer.rs`, `sweep_axes.rs`), so a geom model needs no translation on the way to a mesh.
* The seven failing golden tests (rough spot R22) concern kernels ported under this decision; they fail by machine noise, not by frame errors.

## Code

`crates/femcore/src/section.rs`, `solve.rs`, `ems*.rs`, `sms*.rs`, `sbs*.rs`, `dfs.rs`, `ems_pml.rs::cell_stretches`; `crates/geom/src/stack.rs`, `topo3.rs`. History: `4eab426`, `1d51d96`, `dde3caa`, `f2cb8b7`, `35cceba`, `dcb8ff2` to `a0bb492`.

## Contract

`contracts/frames` (`crates/femcore/tests/contracts.rs`, tests `frames_*`) pins the constants and name tables of `section`.
