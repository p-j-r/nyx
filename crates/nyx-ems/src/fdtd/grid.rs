// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cache-aligned Structure-of-Arrays (SoA) Yee field storage and non-uniform graded rectilinear grid generation (`AutoGrid`,
//! `LayerRefinementSpec`, and `GridOverrideBox`) from `geom` models or `Mesh` domains (ADR 145).

use femcore::geom::eval::Evaluated;
use femcore::geom::labels::LabelMap;
use femcore::geom::model::RegionKind;
use nyx_physics::Mesh;
use serde::{Deserialize, Serialize};

/// Speed of light in vacuum (m/s).
pub const C0: f64 = 299_792_458.0;
/// Vacuum permeability (H/m).
pub const MU0: f64 = 4.0 * std::f64::consts::PI * 1e-7;
/// Vacuum permittivity (F/m).
pub const EPS0: f64 = 1.0 / (MU0 * C0 * C0);

/// 64-byte cache-line / AVX-512 aligned block of eight `f64` lanes, structured for SIMD vectorization
/// and direct Phase 4 GPU (`CubeCL`) contiguous buffer mapping.
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Align64(pub [f64; 8]);

impl Default for Align64 {
    fn default() -> Self {
        Self([0.0; 8])
    }
}

/// Contiguous 64-byte-aligned 1D buffer backing a 3D `(nx + 1, ny + 1, nz_pad)` Yee field component slab.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignedSlab {
    blocks: Vec<Align64>,
    len: usize,
}

impl AlignedSlab {
    /// Allocates a zero-initialized 64-byte-aligned buffer of at least `len` `f64` elements (rounded up to a multiple of 8).
    pub fn zeros(len: usize) -> Self {
        let n_blocks = len.div_ceil(8).max(1);
        Self { blocks: vec![Align64::default(); n_blocks], len: n_blocks * 8 }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Pointer alignment check (always a multiple of 64 bytes).
    #[inline]
    pub fn is_cache_aligned(&self) -> bool {
        (self.blocks.as_ptr() as usize).is_multiple_of(64)
    }

    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        // SAFETY: `Align64` is `#[repr(C, align(64))]` containing `[f64; 8]` with size `64 == 8 * size_of::<f64>()`
        // and no padding, so `blocks.len() * 8` contiguous `f64` values are valid for shared reads.
        unsafe { std::slice::from_raw_parts(self.blocks.as_ptr().cast::<f64>(), self.blocks.len() * 8) }
    }

    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        // SAFETY: `Align64` is `#[repr(C, align(64))]` containing `[f64; 8]` with size `64 == 8 * size_of::<f64>()`
        // and no padding, so `blocks.len() * 8` contiguous `f64` values are valid for exclusive writes.
        unsafe { std::slice::from_raw_parts_mut(self.blocks.as_mut_ptr().cast::<f64>(), self.blocks.len() * 8) }
    }

    pub fn fill_zero(&mut self) {
        self.as_mut_slice().fill(0.0);
    }
}

/// Polarization decoupling mode on a 2D wafer `(x, z)` cross-section (`∂/∂y = 0`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Pol2d {
    /// Out-of-plane electric field `(Ey, Hx, Hz)` (`TM_y` / s-polarized in the `x-z` plane).
    #[default]
    TmEy,
    /// In-plane electric field `(Ex, Ez, Hy)` (`TE_y` / p-polarized in the `x-z` plane).
    TeHy,
    /// Full 6-component 2.5D update `(Ex, Ey, Ez, Hx, Hy, Hz)` on `(nx, 1, nz)` for anisotropic/coupled media.
    Full,
}

impl Pol2d {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "tm" | "tm_ey" | "ey" | "" => Ok(Self::TmEy),
            "te" | "te_hy" | "hy" => Ok(Self::TeHy),
            "full" | "both" => Ok(Self::Full),
            other => Err(format!("unknown 2D FDTD polarization '{other}' (expected 'tm_ey', 'te_hy', or 'full')")),
        }
    }
}

/// Spatial dimensionality of an FDTD simulation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "dim", rename_all = "snake_case")]
pub enum FdtdDim {
    TwoD { polarization: Pol2d },
    ThreeD,
}

/// Automatic wavelength-aware non-uniform grid sizing specification.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AutoGridSpec {
    /// Shortest vacuum wavelength `lambda_min` (m) resolved by the simulation.
    pub wavelength_min_m: f64,
    /// Minimum grid steps per local medium wavelength `lambda_min / n_max` (default `20.0`).
    #[serde(default = "default_steps_per_wave")]
    pub min_steps_per_wave: f64,
    /// Maximum adjacent cell size ratio `Delta_{i+1} / Delta_i` across graded transitions (default `1.35`).
    #[serde(default = "default_max_grading")]
    pub max_grading_ratio: f64,
    /// Optional global floor on cell step size (m) to guard against degenerate micro-features.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dl_min_m: Option<f64>,
}

fn default_steps_per_wave() -> f64 {
    20.0
}

fn default_max_grading() -> f64 {
    1.35
}

impl Default for AutoGridSpec {
    fn default() -> Self {
        Self {
            wavelength_min_m: 1.5e-6,
            min_steps_per_wave: default_steps_per_wave(),
            max_grading_ratio: default_max_grading(),
            dl_min_m: None,
        }
    }
}

/// Planar layer / interface refinement specification along a wafer axis (`0 = x`, `1 = y`, `2 = z`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LayerRefinementSpec {
    /// Axis index (`0 = x`, `1 = y`, `2 = z`).
    pub axis: usize,
    /// `[lo, hi]` coordinate interval (m) along `axis` to refine and snap interfaces to.
    pub bounds_m: [f64; 2],
    /// Maximum step size (m) inside `[lo, hi]`.
    pub max_step_m: f64,
    /// Optional finer step size (m) applied within `edge_band_m` of the interval endpoints (corner/interface refinement).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edge_step_m: Option<f64>,
    /// Half-width (m) around `bounds_m[0]` and `bounds_m[1]` where `edge_step_m` applies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edge_band_m: Option<f64>,
}

/// Axis-aligned 3D override box imposing local maximum Yee cell steps `(dx, dy, dz)`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GridOverrideBox {
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    /// Target maximum step `[dx, dy, dz]` (m) on each axis (`None` leaves that axis unaffected).
    pub dl: [Option<f64>; 3],
}

/// Full grid generation specification (`AutoGrid`, explicit uniform steps, or explicit vertex coordinates).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct GridSpec {
    /// Automatic wavelength-graded grid parameters (used unless `uniform_dl` or `explicit_vertices` overrides an axis).
    #[serde(default)]
    pub auto: AutoGridSpec,
    /// Optional uniform cell size `[dx, dy, dz]` (m) overriding `auto` when positive.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uniform_dl: Option<[f64; 3]>,
    /// Optional explicit monotonically increasing primal vertex coordinates `[vx, vy, vz]` (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub explicit_vertices: Option<[Vec<f64>; 3]>,
    /// Optional simulation bounding box `([x0, y0, z0], [x1, y1, z1])` overriding the `geom` / `Mesh` extent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounds_override: Option<([f64; 3], [f64; 3])>,
    /// Layer refinement specs along individual wafer axes.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub layer_refinements: Vec<LayerRefinementSpec>,
    /// Local 3D refinement override boxes.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub override_boxes: Vec<GridOverrideBox>,
    /// Extra interior planes `(axis, position_m)` that primal grid vertices should snap to (e.g. port/monitor/TFSF planes).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub snap_planes: Vec<(usize, f64)>,
}

impl GridSpec {
    pub fn check_finite(&self) -> Result<(), String> {
        if !(self.auto.wavelength_min_m.is_finite() && self.auto.wavelength_min_m > 0.0) {
            return Err(format!(
                "FDTD wavelength_min_m must be positive and finite (got {})",
                self.auto.wavelength_min_m
            ));
        }
        if !(self.auto.min_steps_per_wave.is_finite() && self.auto.min_steps_per_wave >= 2.0) {
            return Err(format!(
                "FDTD min_steps_per_wave must be finite and >= 2.0 (got {})",
                self.auto.min_steps_per_wave
            ));
        }
        if !(self.auto.max_grading_ratio.is_finite() && self.auto.max_grading_ratio > 1.0) {
            return Err(format!(
                "FDTD max_grading_ratio must be finite and > 1.0 (got {})",
                self.auto.max_grading_ratio
            ));
        }
        if let Some(dmin) = self.auto.dl_min_m {
            if !(dmin.is_finite() && dmin > 0.0) {
                return Err(format!("FDTD dl_min_m must be positive and finite (got {dmin})"));
            }
        }
        if let Some(udl) = self.uniform_dl {
            if udl.iter().any(|&d| !d.is_finite() || d < 0.0) {
                return Err("FDTD uniform_dl entries must be finite and non-negative".into());
            }
        }
        if let Some((lo, hi)) = self.bounds_override {
            if lo.iter().chain(&hi).any(|v| !v.is_finite()) {
                return Err("FDTD bounds_override contains NaN or infinite coordinates".into());
            }
        }
        for lr in &self.layer_refinements {
            if lr.axis > 2
                || !lr.bounds_m[0].is_finite()
                || !lr.bounds_m[1].is_finite()
                || !(lr.max_step_m.is_finite() && lr.max_step_m > 0.0)
                || lr.edge_step_m.is_some_and(|s| !s.is_finite() || s <= 0.0)
                || lr.edge_band_m.is_some_and(|b| !b.is_finite() || b < 0.0)
            {
                return Err("FDTD layer_refinement contains invalid or non-finite parameters".into());
            }
        }
        for ob in &self.override_boxes {
            if ob.lo.iter().chain(&ob.hi).any(|v| !v.is_finite())
                || ob.dl.iter().flatten().any(|&d| !d.is_finite() || d <= 0.0)
            {
                return Err("FDTD override_box contains invalid or non-finite parameters".into());
            }
        }
        for &(ax, pos) in &self.snap_planes {
            if ax > 2 || !pos.is_finite() {
                return Err("FDTD snap_plane contains invalid axis or non-finite coordinate".into());
            }
        }
        Ok(())
    }
}

/// Non-uniform staggered rectilinear Yee grid in wafer coordinates `(x, y, z)`.
///
/// For 2D simulations in the `(x, z)` wafer cross-section, `ny == 1` with `vy = [-0.5, 0.5]` and `dim = FdtdDim::TwoD`.
/// Every 3D SoA field slab has shape `(nx + 1, ny + 1, nz_pad)` where `nz_pad = ((nz + 1) + 7) & !7`, guaranteeing
/// 64-byte cache-line alignment at the start of every `(ix, iy)` z-pencil.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct YeeGrid {
    pub dim: FdtdDim,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub nz_pad: usize,
    pub stride_x: usize,
    pub stride_y: usize,
    /// Primal grid vertices (`nx + 1`, `ny + 1`, `nz + 1`) in metres.
    pub vx: Vec<f64>,
    pub vy: Vec<f64>,
    pub vz: Vec<f64>,
    /// Primal cell edge lengths `dp_a[i] = v_a[i + 1] - v_a[i]` (`nx`, `ny`, `nz`).
    pub dp_x: Vec<f64>,
    pub dp_y: Vec<f64>,
    pub dp_z: Vec<f64>,
    /// Reciprocal primal edge lengths `1 / dp_a[i]`.
    pub inv_dp_x: Vec<f64>,
    pub inv_dp_y: Vec<f64>,
    pub inv_dp_z: Vec<f64>,
    /// Cell centers (dual vertices) `c_a[i] = 0.5 * (v_a[i] + v_a[i + 1])` (`nx`, `ny`, `nz`).
    pub cx: Vec<f64>,
    pub cy: Vec<f64>,
    pub cz: Vec<f64>,
    /// Dual edge lengths `dd_a[i]` (`nx + 1`, `ny + 1`, `nz + 1`) associated with primal vertices `v_a[i]`.
    pub dd_x: Vec<f64>,
    pub dd_y: Vec<f64>,
    pub dd_z: Vec<f64>,
    /// Reciprocal dual edge lengths `1 / dd_a[i]`.
    pub inv_dd_x: Vec<f64>,
    pub inv_dd_y: Vec<f64>,
    pub inv_dd_z: Vec<f64>,
}

/// One region refractive-index interval `(lo, hi, n_max)` along a coordinate axis used by `AutoGrid`.
#[derive(Clone, Copy, Debug)]
pub struct AxisIndexBox {
    pub lo: f64,
    pub hi: f64,
    pub n_max: f64,
}

impl YeeGrid {
    /// Builds a `YeeGrid` from explicit primal vertex arrays `(vx, vy, vz)` and `periodic` flags `[px, py, pz]`.
    pub fn from_vertices(
        dim: FdtdDim,
        vx: Vec<f64>,
        vy: Vec<f64>,
        vz: Vec<f64>,
        periodic: [bool; 3],
    ) -> Result<Self, String> {
        let check_axis = |name: &str, v: &[f64]| -> Result<(), String> {
            if v.len() < 2 {
                return Err(format!("YeeGrid axis '{name}' requires at least 2 primal vertices (got {})", v.len()));
            }
            for (i, &x) in v.iter().enumerate() {
                if !x.is_finite() {
                    return Err(format!("YeeGrid axis '{name}' vertex {i} is non-finite ({x})"));
                }
            }
            for w in v.windows(2) {
                if w[1] <= w[0] {
                    return Err(format!(
                        "YeeGrid axis '{name}' vertices must be strictly ascending ({:.6e} >= {:.6e})",
                        w[0], w[1]
                    ));
                }
            }
            Ok(())
        };
        check_axis("x", &vx)?;
        check_axis("y", &vy)?;
        check_axis("z", &vz)?;

        let (nx, ny, nz) = (vx.len() - 1, vy.len() - 1, vz.len() - 1);
        if matches!(dim, FdtdDim::TwoD { .. }) && ny != 1 {
            return Err(format!("2D YeeGrid in wafer (x, z) plane requires ny == 1 (got ny = {ny})"));
        }

        let build_1d = |v: &[f64], is_per: bool| {
            let n = v.len() - 1;
            let dp: Vec<f64> = (0..n).map(|i| v[i + 1] - v[i]).collect();
            let inv_dp: Vec<f64> = dp.iter().map(|&d| 1.0 / d).collect();
            let c: Vec<f64> = (0..n).map(|i| 0.5 * (v[i] + v[i + 1])).collect();
            let mut dd = vec![0.0; n + 1];
            for i in 1..n {
                dd[i] = c[i] - c[i - 1];
            }
            if is_per {
                let wrap = 0.5 * (dp[0] + dp[n - 1]);
                dd[0] = wrap;
                dd[n] = wrap;
            } else {
                dd[0] = dp[0];
                dd[n] = dp[n - 1];
            }
            let inv_dd: Vec<f64> = dd.iter().map(|&d| 1.0 / d).collect();
            (dp, inv_dp, c, dd, inv_dd)
        };

        let (dp_x, inv_dp_x, cx, dd_x, inv_dd_x) = build_1d(&vx, periodic[0]);
        let (dp_y, inv_dp_y, cy, dd_y, inv_dd_y) = build_1d(&vy, periodic[1] || matches!(dim, FdtdDim::TwoD { .. }));
        let (dp_z, inv_dp_z, cz, dd_z, inv_dd_z) = build_1d(&vz, periodic[2]);

        let nz_pad = ((nz + 1) + 7) & !7;
        let stride_y = nz_pad;
        let stride_x = (ny + 1) * nz_pad;

        Ok(Self {
            dim,
            nx,
            ny,
            nz,
            nz_pad,
            stride_x,
            stride_y,
            vx,
            vy,
            vz,
            dp_x,
            dp_y,
            dp_z,
            inv_dp_x,
            inv_dp_y,
            inv_dp_z,
            cx,
            cy,
            cz,
            dd_x,
            dd_y,
            dd_z,
            inv_dd_x,
            inv_dd_y,
            inv_dd_z,
        })
    }

    /// Generates a wavelength-aware graded `YeeGrid` directly from an evaluated `geom` model (`Evaluated` + `LabelMap`)
    /// and a per-domain maximum refractive index map (`domain_n_max`).
    pub fn from_geom(
        ev: &Evaluated,
        labels: &LabelMap,
        pol_2d: Pol2d,
        spec: &GridSpec,
        domain_n_max: &std::collections::BTreeMap<i32, f64>,
        periodic: [bool; 3],
    ) -> Result<Self, String> {
        spec.check_finite()?;
        let is_2d = ev.list.dim == 2;
        let dim = if is_2d { FdtdDim::TwoD { polarization: pol_2d } } else { FdtdDim::ThreeD };
        let b = ev.bounds();
        let (lo, hi) = if let Some((blo, bhi)) = spec.bounds_override {
            (blo, bhi)
        } else if is_2d {
            ([b.lo[0], -0.5, b.lo[2]], [b.hi[0], 0.5, b.hi[2]])
        } else {
            (b.lo, b.hi)
        };

        let mut axis_boxes: [Vec<AxisIndexBox>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        let mut auto_snaps: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];

        if ev.complex.is_some() {
            for (ri, reg) in ev.regions().iter().enumerate() {
                if matches!(reg.kind, RegionKind::Hole | RegionKind::Void) || reg.bbox.is_empty() {
                    continue;
                }
                let dom_id = labels.domains.get(ri).map(|d| d.id).unwrap_or((ri + 1) as i32);
                let n_reg = domain_n_max.get(&dom_id).copied().unwrap_or(1.0).max(1.0);
                for ax in 0..3 {
                    if is_2d && ax == 1 {
                        continue;
                    }
                    let (rlo, rhi) = (reg.bbox.lo[ax], reg.bbox.hi[ax]);
                    if rhi > rlo {
                        axis_boxes[ax].push(AxisIndexBox { lo: rlo, hi: rhi, n_max: n_reg });
                        auto_snaps[ax].push(rlo);
                        auto_snaps[ax].push(rhi);
                    }
                }
            }
            if let Some(st) = ev.layer_stack() {
                for l in &st.layers {
                    auto_snaps[2].push(l.z0);
                    auto_snaps[2].push(l.z1);
                }
            }
        }

        Self::build_graded_axes(dim, lo, hi, spec, &axis_boxes, &auto_snaps, periodic)
    }

    /// Generates a wavelength-aware graded `YeeGrid` from a `nyx_physics::Mesh` and per-domain maximum refractive indices.
    pub fn from_mesh(
        mesh: &Mesh,
        pol_2d: Pol2d,
        spec: &GridSpec,
        domain_n_max: &std::collections::BTreeMap<i32, f64>,
        periodic: [bool; 3],
    ) -> Result<Self, String> {
        spec.check_finite()?;
        let simplex = nyx_physics::discretise(mesh).map_err(|e| e.to_string())?;
        let is_2d = simplex.dim == 2;
        let dim = if is_2d { FdtdDim::TwoD { polarization: pol_2d } } else { FdtdDim::ThreeD };
        let (mlo, mhi) = simplex.bbox();
        let (lo, hi) = if let Some((blo, bhi)) = spec.bounds_override {
            (blo, bhi)
        } else if is_2d {
            ([mlo[0], -0.5, mlo[1]], [mhi[0], 0.5, mhi[1]])
        } else {
            ([mlo[0], mlo[1], mlo[2]], [mhi[0], mhi[1], mhi[2]])
        };

        let mut dom_bounds: std::collections::BTreeMap<i32, ([f64; 3], [f64; 3])> = std::collections::BTreeMap::new();
        for ci in 0..simplex.n_cells() {
            let lab = simplex.cell_labels[ci];
            let entry = dom_bounds.entry(lab).or_insert(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]));
            for &v in simplex.cell(ci) {
                let p = simplex.point(v as usize);
                let wp = if is_2d { [p[0], 0.0, p[1]] } else { [p[0], p[1], p[2]] };
                for k in 0..3 {
                    entry.0[k] = entry.0[k].min(wp[k]);
                    entry.1[k] = entry.1[k].max(wp[k]);
                }
            }
        }

        let mut axis_boxes: [Vec<AxisIndexBox>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        let mut auto_snaps: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for (lab, (dlo, dhi)) in dom_bounds {
            let n_reg = domain_n_max.get(&lab).copied().unwrap_or(1.0).max(1.0);
            for ax in 0..3 {
                if is_2d && ax == 1 {
                    continue;
                }
                if dhi[ax] > dlo[ax] {
                    axis_boxes[ax].push(AxisIndexBox { lo: dlo[ax], hi: dhi[ax], n_max: n_reg });
                    auto_snaps[ax].push(dlo[ax]);
                    auto_snaps[ax].push(dhi[ax]);
                }
            }
        }

        Self::build_graded_axes(dim, lo, hi, spec, &axis_boxes, &auto_snaps, periodic)
    }

    fn build_graded_axes(
        dim: FdtdDim,
        lo: [f64; 3],
        hi: [f64; 3],
        spec: &GridSpec,
        axis_boxes: &[Vec<AxisIndexBox>; 3],
        auto_snaps: &[Vec<f64>; 3],
        periodic: [bool; 3],
    ) -> Result<Self, String> {
        let is_2d = matches!(dim, FdtdDim::TwoD { .. });
        let mut verts: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];

        for ax in 0..3 {
            if is_2d && ax == 1 {
                verts[1] = vec![-0.5, 0.5];
                continue;
            }
            if !hi[ax].is_finite() || !lo[ax].is_finite() || hi[ax] <= lo[ax] {
                return Err(format!(
                    "FDTD domain bounds along axis {ax} must satisfy lo < hi (got [{:.6e}, {:.6e}])",
                    lo[ax], hi[ax]
                ));
            }
            if let Some(ref exp) = spec.explicit_vertices {
                if exp[ax].len() >= 2 {
                    verts[ax] = exp[ax].clone();
                    continue;
                }
            }
            if let Some(udl) = spec.uniform_dl {
                if udl[ax] > 0.0 {
                    let span = hi[ax] - lo[ax];
                    let n = (span / udl[ax]).ceil().max(2.0) as usize;
                    verts[ax] = (0..=n).map(|i| lo[ax] + span * (i as f64) / (n as f64)).collect();
                    continue;
                }
            }
            verts[ax] = generate_graded_1d_axis(ax, lo[ax], hi[ax], spec, &axis_boxes[ax], &auto_snaps[ax])?;
        }

        Self::from_vertices(dim, verts[0].clone(), verts[1].clone(), verts[2].clone(), periodic)
    }

    /// Flat index into any `(nx + 1, ny + 1, nz_pad)` SoA slab.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.stride_x + iy * self.stride_y + iz
    }

    /// Total number of `f64` slots per SoA component slab.
    #[inline]
    pub fn slab_len(&self) -> usize {
        (self.nx + 1) * (self.ny + 1) * self.nz_pad
    }

    /// Total number of Yee primal cells `nx * ny * nz`.
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Smallest primal cell step across all active axes (m).
    pub fn min_step_m(&self) -> f64 {
        let mx = self.dp_x.iter().copied().fold(f64::INFINITY, f64::min);
        let mz = self.dp_z.iter().copied().fold(f64::INFINITY, f64::min);
        if matches!(self.dim, FdtdDim::TwoD { .. }) {
            mx.min(mz)
        } else {
            let my = self.dp_y.iter().copied().fold(f64::INFINITY, f64::min);
            mx.min(my).min(mz)
        }
    }

    /// Largest primal cell step across all active axes (m).
    pub fn max_step_m(&self) -> f64 {
        let mx = self.dp_x.iter().copied().fold(0.0_f64, f64::max);
        let mz = self.dp_z.iter().copied().fold(0.0_f64, f64::max);
        if matches!(self.dim, FdtdDim::TwoD { .. }) {
            mx.max(mz)
        } else {
            let my = self.dp_y.iter().copied().fold(0.0_f64, f64::max);
            mx.max(my).max(mz)
        }
    }

    /// Simulation domain `[lo, hi]` bounding box in wafer coordinates.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        ([self.vx[0], self.vy[0], self.vz[0]], [self.vx[self.nx], self.vy[self.ny], self.vz[self.nz]])
    }

    /// Computes the Courant–Friedrichs–Lewy (CFL) stable time step `dt` (s) given the minimum local background/high-frequency
    /// refractive index `n_min_cell(ix, iy, iz)` (or global `n_min >= 1.0`) and `courant` safety factor in `(0, 1)`.
    pub fn cfl_dt(&self, courant: f64, min_refractive_index: f64, conformal_pec_scale: f64) -> f64 {
        let s = courant.clamp(0.01, 0.995) * conformal_pec_scale.clamp(0.1, 1.0);
        let n_min = min_refractive_index.max(1.0);
        let is_2d = matches!(self.dim, FdtdDim::TwoD { .. });
        let mut max_inv_ds2 = 0.0_f64;
        for ix in 0..self.nx {
            let idx2 = self.inv_dp_x[ix] * self.inv_dp_x[ix];
            if is_2d {
                for iz in 0..self.nz {
                    let idz2 = self.inv_dp_z[iz] * self.inv_dp_z[iz];
                    max_inv_ds2 = max_inv_ds2.max(idx2 + idz2);
                }
            } else {
                for iy in 0..self.ny {
                    let idxy2 = idx2 + self.inv_dp_y[iy] * self.inv_dp_y[iy];
                    for iz in 0..self.nz {
                        let tot = idxy2 + self.inv_dp_z[iz] * self.inv_dp_z[iz];
                        max_inv_ds2 = max_inv_ds2.max(tot);
                    }
                }
            }
        }
        s * n_min / (C0 * max_inv_ds2.max(1e-30).sqrt())
    }

    /// Finds the primal vertex index `i` in `0..=n` on `axis` (`0=x, 1=y, 2=z`) nearest `pos`.
    pub fn nearest_vertex(&self, axis: usize, pos: f64) -> usize {
        let v = self.axis_vertices(axis);
        let mut best = 0usize;
        let mut best_d = f64::INFINITY;
        for (i, &x) in v.iter().enumerate() {
            let d = (x - pos).abs();
            if d < best_d {
                best_d = d;
                best = i;
            }
        }
        best
    }

    /// Finds the cell index `i` in `0..n` on `axis` containing `pos` (clamped to `[0, n - 1]`).
    pub fn locate_cell(&self, axis: usize, pos: f64) -> usize {
        let v = self.axis_vertices(axis);
        let n = v.len() - 1;
        if pos <= v[0] {
            return 0;
        }
        if pos >= v[n] {
            return n - 1;
        }
        match v.binary_search_by(|x| x.total_cmp(&pos)) {
            Ok(i) => i.min(n - 1),
            Err(i) => i.saturating_sub(1).min(n - 1),
        }
    }

    #[inline]
    pub fn axis_vertices(&self, axis: usize) -> &[f64] {
        match axis {
            0 => &self.vx,
            1 => &self.vy,
            _ => &self.vz,
        }
    }

    #[inline]
    pub fn axis_centers(&self, axis: usize) -> &[f64] {
        match axis {
            0 => &self.cx,
            1 => &self.cy,
            _ => &self.cz,
        }
    }

    #[inline]
    pub fn axis_dp(&self, axis: usize) -> &[f64] {
        match axis {
            0 => &self.dp_x,
            1 => &self.dp_y,
            _ => &self.dp_z,
        }
    }

    #[inline]
    pub fn axis_dd(&self, axis: usize) -> &[f64] {
        match axis {
            0 => &self.dd_x,
            1 => &self.dd_y,
            _ => &self.dd_z,
        }
    }
}

/// Generates 1D smoothly graded primal vertices on `[lo, hi]` honoring `AutoGridSpec`, `LayerRefinementSpec`,
/// `GridOverrideBox`, and snapped structural planes.
fn generate_graded_1d_axis(
    axis: usize,
    lo: f64,
    hi: f64,
    spec: &GridSpec,
    boxes: &[AxisIndexBox],
    auto_snaps: &[f64],
) -> Result<Vec<f64>, String> {
    let span = hi - lo;
    let base_dl = (spec.auto.wavelength_min_m / spec.auto.min_steps_per_wave).max(span * 1e-5);
    let dl_floor = spec.auto.dl_min_m.unwrap_or(base_dl * 1e-2).max(span * 1e-6);
    let r_max = spec.auto.max_grading_ratio.clamp(1.02, 2.0);

    // Collect snapped breakpoints inside [lo, hi]
    let tol_snap = 1e-7 * span.max(1e-12);
    let mut snaps = vec![lo, hi];
    let push_clamped = |snaps: &mut Vec<f64>, x: f64| {
        if x > lo + tol_snap && x < hi - tol_snap {
            snaps.push(x);
        }
    };
    for &s in auto_snaps {
        push_clamped(&mut snaps, s);
    }
    for &(ax, pos) in &spec.snap_planes {
        if ax == axis {
            push_clamped(&mut snaps, pos);
        }
    }
    for lr in &spec.layer_refinements {
        if lr.axis == axis {
            let (a, b) = (lr.bounds_m[0].min(lr.bounds_m[1]), lr.bounds_m[0].max(lr.bounds_m[1]));
            push_clamped(&mut snaps, a);
            push_clamped(&mut snaps, b);
            if let Some(band) = lr.edge_band_m.filter(|&eb| eb > 0.0) {
                push_clamped(&mut snaps, a - band);
                push_clamped(&mut snaps, a + band);
                push_clamped(&mut snaps, b - band);
                push_clamped(&mut snaps, b + band);
            }
        }
    }
    for ob in &spec.override_boxes {
        if ob.dl[axis].is_some() {
            push_clamped(&mut snaps, ob.lo[axis].min(ob.hi[axis]));
            push_clamped(&mut snaps, ob.lo[axis].max(ob.hi[axis]));
        }
    }
    snaps.sort_by(|a, b| a.total_cmp(b));
    let mut knots = Vec::with_capacity(snaps.len());
    for x in snaps {
        if knots.last().is_none_or(|&prev: &f64| (x - prev).abs() > tol_snap.max(0.25 * dl_floor)) {
            knots.push(x);
        }
    }
    if let Some(last) = knots.last_mut() {
        *last = hi;
    }
    if knots.len() < 2 {
        knots = vec![lo, hi];
    }

    // Target step size on each sub-interval `[knots[m], knots[m + 1]]`
    let n_seg = knots.len() - 1;
    let mut seg_dl = vec![base_dl; n_seg];
    for m in 0..n_seg {
        let mid = 0.5 * (knots[m] + knots[m + 1]);
        let mut want = base_dl;
        for bx in boxes {
            if mid >= bx.lo - tol_snap && mid <= bx.hi + tol_snap {
                want = want.min(spec.auto.wavelength_min_m / (bx.n_max.max(1.0) * spec.auto.min_steps_per_wave));
            }
        }
        for lr in &spec.layer_refinements {
            if lr.axis == axis {
                let (a, b) = (lr.bounds_m[0].min(lr.bounds_m[1]), lr.bounds_m[0].max(lr.bounds_m[1]));
                if mid >= a - tol_snap && mid <= b + tol_snap {
                    want = want.min(lr.max_step_m);
                }
                if let (Some(es), Some(eb)) = (lr.edge_step_m, lr.edge_band_m) {
                    if (mid - a).abs() <= eb + tol_snap || (mid - b).abs() <= eb + tol_snap {
                        want = want.min(es);
                    }
                }
            }
        }
        for ob in &spec.override_boxes {
            if let Some(odl) = ob.dl[axis] {
                let (a, b) = (ob.lo[axis].min(ob.hi[axis]), ob.lo[axis].max(ob.hi[axis]));
                if mid >= a - tol_snap && mid <= b + tol_snap {
                    want = want.min(odl);
                }
            }
        }
        seg_dl[m] = want.max(dl_floor);
    }

    // Forward/backward Lipschitz smoothing across interval targets so boundary steps match smoothly
    let mut knot_dl = vec![base_dl; knots.len()];
    for k in 0..knots.len() {
        let mut d = f64::INFINITY;
        if k > 0 {
            d = d.min(seg_dl[k - 1]);
        }
        if k < n_seg {
            d = d.min(seg_dl[k]);
        }
        knot_dl[k] = d.max(dl_floor);
    }
    let lip = r_max - 1.0;
    for k in 1..knots.len() {
        let dist = knots[k] - knots[k - 1];
        knot_dl[k] = knot_dl[k].min(knot_dl[k - 1] + lip * dist);
    }
    for k in (0..knots.len() - 1).rev() {
        let dist = knots[k + 1] - knots[k];
        knot_dl[k] = knot_dl[k].min(knot_dl[k + 1] + lip * dist);
    }

    // Discretize each `[knots[m], knots[m + 1]]` with symmetric geometric grading from `knot_dl[m]`, `seg_dl[m]`, `knot_dl[m + 1]`
    let mut out = vec![lo];
    for m in 0..n_seg {
        let (a, b) = (knots[m], knots[m + 1]);
        let w = b - a;
        let d_left = knot_dl[m].min(seg_dl[m]).min(w);
        let d_right = knot_dl[m + 1].min(seg_dl[m]).min(w);
        let d_mid = seg_dl[m].min(w);
        let sub = fill_graded_interval(a, b, d_left, d_mid, d_right, r_max);
        for &x in sub.iter().skip(1) {
            out.push(x);
        }
    }
    if out.len() < 3 {
        let mid = 0.5 * (lo + hi);
        out = vec![lo, mid, hi];
    }
    Ok(out)
}

/// Fills `[a, b]` with strictly ascending vertices starting near step `d_left`, plateauing at `d_mid`,
/// and ending near `d_right`, with adjacent step ratio bounded by `r_max`.
fn fill_graded_interval(a: f64, b: f64, d_left: f64, d_mid: f64, d_right: f64, r_max: f64) -> Vec<f64> {
    let span = b - a;
    if span <= 1.05 * d_mid.max(d_left.min(d_right)) {
        let n = (span / d_mid.min(d_left).min(d_right).max(span * 1e-4)).ceil().max(1.0) as usize;
        return (0..=n).map(|i| if i == n { b } else { a + span * (i as f64) / (n as f64) }).collect();
    }

    let mut left_steps = Vec::new();
    let mut right_steps = Vec::new();
    let mut cur_l = d_left;
    let mut cur_r = d_right;
    let mut used = 0.0;

    // Grow symmetrically from both ends toward `d_mid` while total < span
    while (cur_l < d_mid * (1.0 - 1e-6) || cur_r < d_mid * (1.0 - 1e-6)) && used + cur_l + cur_r <= span {
        if cur_l <= cur_r && cur_l < d_mid {
            left_steps.push(cur_l);
            used += cur_l;
            cur_l = (cur_l * r_max).min(d_mid);
        } else if cur_r < d_mid {
            right_steps.push(cur_r);
            used += cur_r;
            cur_r = (cur_r * r_max).min(d_mid);
        } else {
            break;
        }
    }

    let rem = (span - used).max(0.0);
    let n_mid = if rem > 0.0 { (rem / d_mid).ceil().max(1.0) as usize } else { 0 };
    let mut all_steps = Vec::with_capacity(left_steps.len() + n_mid + right_steps.len());
    all_steps.extend_from_slice(&left_steps);
    if n_mid > 0 {
        let dm = rem / (n_mid as f64);
        for _ in 0..n_mid {
            all_steps.push(dm);
        }
    }
    for &s in right_steps.iter().rev() {
        all_steps.push(s);
    }
    if all_steps.is_empty() {
        return vec![a, b];
    }
    // Rescale steps to sum identically to `span`
    let sum: f64 = all_steps.iter().sum();
    let scale = span / sum;
    let mut verts = Vec::with_capacity(all_steps.len() + 1);
    verts.push(a);
    let mut acc = a;
    for (i, &s) in all_steps.iter().enumerate() {
        if i + 1 == all_steps.len() {
            verts.push(b);
        } else {
            acc += s * scale;
            verts.push(acc);
        }
    }
    verts
}

/// Cache-aligned Structure-of-Arrays (SoA) electromagnetic field state `(Ex, Ey, Ez, Hx, Hy, Hz)`.
#[derive(Clone, Debug)]
pub struct YeeFields {
    pub ex: AlignedSlab,
    pub ey: AlignedSlab,
    pub ez: AlignedSlab,
    pub hx: AlignedSlab,
    pub hy: AlignedSlab,
    pub hz: AlignedSlab,
}

impl YeeFields {
    pub fn new(grid: &YeeGrid) -> Self {
        let n = grid.slab_len();
        Self {
            ex: AlignedSlab::zeros(n),
            ey: AlignedSlab::zeros(n),
            ez: AlignedSlab::zeros(n),
            hx: AlignedSlab::zeros(n),
            hy: AlignedSlab::zeros(n),
            hz: AlignedSlab::zeros(n),
        }
    }

    pub fn fill_zero(&mut self) {
        self.ex.fill_zero();
        self.ey.fill_zero();
        self.ez.fill_zero();
        self.hx.fill_zero();
        self.hy.fill_zero();
        self.hz.fill_zero();
    }

    /// Interpolates `(Ex, Ey, Ez)` and `(Hx, Hy, Hz)` at an arbitrary wafer point `p = [x, y, z]`.
    pub fn sample_eh(&self, grid: &YeeGrid, p: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let ex = sample_staggered(grid, self.ex.as_slice(), [true, false, false], p);
        let ey = sample_staggered(grid, self.ey.as_slice(), [false, true, false], p);
        let ez = sample_staggered(grid, self.ez.as_slice(), [false, false, true], p);
        let hx = sample_staggered(grid, self.hx.as_slice(), [false, true, true], p);
        let hy = sample_staggered(grid, self.hy.as_slice(), [true, false, true], p);
        let hz = sample_staggered(grid, self.hz.as_slice(), [true, true, false], p);
        ([ex, ey, ez], [hx, hy, hz])
    }
}

/// Trilinear (or 2D bilinear when `grid.dim` is `TwoD`) interpolation of one staggered Yee component
/// (`is_center[a] == true` -> sampled on `grid.axis_centers(a)`, `false` -> sampled on `grid.axis_vertices(a)`).
pub fn sample_staggered(grid: &YeeGrid, slab: &[f64], is_center: [bool; 3], p: [f64; 3]) -> f64 {
    let bracket = |coords: &[f64], x: f64| -> (usize, usize, f64) {
        let m = coords.len();
        if m <= 1 || x <= coords[0] {
            return (0, 0, 0.0);
        }
        if x >= coords[m - 1] {
            return (m - 1, m - 1, 0.0);
        }
        let hi = match coords.binary_search_by(|v| v.total_cmp(&x)) {
            Ok(i) => return (i, i, 0.0),
            Err(i) => i.clamp(1, m - 1),
        };
        let lo = hi - 1;
        let denom = coords[hi] - coords[lo];
        let t = if denom.abs() > 1e-30 { ((x - coords[lo]) / denom).clamp(0.0, 1.0) } else { 0.0 };
        (lo, hi, t)
    };

    let (ix0, ix1, tx) = bracket(if is_center[0] { &grid.cx } else { &grid.vx }, p[0]);
    let (iz0, iz1, tz) = bracket(if is_center[2] { &grid.cz } else { &grid.vz }, p[2]);

    if matches!(grid.dim, FdtdDim::TwoD { .. }) {
        let v00 = slab[grid.idx(ix0, 0, iz0)];
        let v10 = slab[grid.idx(ix1, 0, iz0)];
        let v01 = slab[grid.idx(ix0, 0, iz1)];
        let v11 = slab[grid.idx(ix1, 0, iz1)];
        return (1.0 - tx) * ((1.0 - tz) * v00 + tz * v01) + tx * ((1.0 - tz) * v10 + tz * v11);
    }

    let (iy0, iy1, ty) = bracket(if is_center[1] { &grid.cy } else { &grid.vy }, p[1]);
    let c00 = (1.0 - tx) * slab[grid.idx(ix0, iy0, iz0)] + tx * slab[grid.idx(ix1, iy0, iz0)];
    let c01 = (1.0 - tx) * slab[grid.idx(ix0, iy0, iz1)] + tx * slab[grid.idx(ix1, iy0, iz1)];
    let c10 = (1.0 - tx) * slab[grid.idx(ix0, iy1, iz0)] + tx * slab[grid.idx(ix1, iy1, iz0)];
    let c11 = (1.0 - tx) * slab[grid.idx(ix0, iy1, iz1)] + tx * slab[grid.idx(ix1, iy1, iz1)];
    let c0 = (1.0 - ty) * c00 + ty * c10;
    let c1 = (1.0 - ty) * c01 + ty * c11;
    (1.0 - tz) * c0 + tz * c1
}
