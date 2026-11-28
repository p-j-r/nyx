// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Conformal subpixel dielectric smoothing (`PolarizedAveraging`, Farjadpour et al. 2006) and stabilized Dey–Mittra
//! conformal PEC (`PECConformal`) edge/area fractions sampled from `geom` models or curved/simplex `Mesh` domains (ADR 145).

use super::grid::{AlignedSlab, FdtdDim, YeeGrid};
use femcore::fem::CurvedTet;
use femcore::geom::eval::Evaluated;
use femcore::geom::labels::LabelMap;
use femcore::geom::model::RegionKind;
use femcore::mesh::SimplexMesh;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Dielectric interface smoothing strategy on Yee electric edges.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SubpixelMode {
    /// Point-sampled staircased permittivity at each Yee edge midpoint (`O(Delta x)` interface error).
    Staircased,
    /// Isotropic volume-weighted arithmetic permittivity average `<eps>` over each dual control box.
    Volumetric,
    /// Tensorial Farjadpour–Oskooi–Joannopoulos polarized averaging:
    /// `W = P <eps^-1> P + (I - P) <eps>^-1 (I - P)` (inverse-permittivity normal + arithmetic tangential, `O(Delta x^2)`).
    #[default]
    Polarized,
}

impl SubpixelMode {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "polarized" | "conformal" | "" => Ok(Self::Polarized),
            "volumetric" | "volume" => Ok(Self::Volumetric),
            "staircased" | "none" | "off" => Ok(Self::Staircased),
            other => Err(format!(
                "unknown FDTD subpixel mode '{other}' (expected 'polarized', 'volumetric', or 'staircased')"
            )),
        }
    }
}

/// Constitutive parameters of a single domain in the FDTD background/high-frequency table.
#[derive(Clone, Copy, Debug)]
pub struct DomainMedium {
    /// Symmetric positive-definite high-frequency relative permittivity tensor `eps_inf` in wafer axes `(x, y, z)`.
    pub eps_inf: [[f64; 3]; 3],
    /// Inverse of `eps_inf`.
    pub inv_eps_inf: [[f64; 3]; 3],
    /// Effective isotropic electric conductivity `sigma_e` (S/m) per principal wafer axis (`[sx, sy, sz]`).
    pub sigma_e: [f64; 3],
    /// Relative magnetic permeability `mu_r` (`[mux, muy, muz]`).
    pub mu_r: [f64; 3],
    /// Whether this domain is treated as a volumetric perfect electric conductor (`PEC`).
    pub is_pec: bool,
}

impl Default for DomainMedium {
    fn default() -> Self {
        Self {
            eps_inf: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            inv_eps_inf: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            sigma_e: [0.0; 3],
            mu_r: [1.0; 3],
            is_pec: false,
        }
    }
}

impl DomainMedium {
    pub fn from_eps_tensor(eps_re: [[f64; 3]; 3], sigma_e: [f64; 3], mu_r: [f64; 3], is_pec: bool) -> Self {
        // Symmetrize real permittivity and clamp diagonal away from zero for unconditional SPD invertibility
        let mut s = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] = 0.5 * (eps_re[i][j] + eps_re[j][i]);
            }
            if s[i][i] < 1e-6 {
                s[i][i] = 1.0;
            }
        }
        let inv = invert_sym3(&s).unwrap_or([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        Self {
            eps_inf: s,
            inv_eps_inf: inv,
            sigma_e: [sigma_e[0].max(0.0), sigma_e[1].max(0.0), sigma_e[2].max(0.0)],
            mu_r: [mu_r[0].max(1e-6), mu_r[1].max(1e-6), mu_r[2].max(1e-6)],
            is_pec,
        }
    }
}

/// Off-diagonal coupling entry for anisotropic or tilted-interface `PolarizedAveraging` (`E = W D`).
#[derive(Clone, Copy, Debug)]
pub struct OffDiagEntry {
    /// Flat SoA index `(ix, iy, iz)` of the target electric field component (`Ex`, `Ey`, or `Ez`).
    pub idx: usize,
    pub ix: u32,
    pub iy: u32,
    pub iz: u32,
    /// Off-diagonal inverse permittivity entries `(W_{a,b}, W_{a,c})` coupling orthogonal `D_b`, `D_c` to `E_a`.
    pub w_ab: f64,
    pub w_ac: f64,
}

/// Pre-assembled cell-edge and cell-face material/conformal coefficients on the Yee grid.
#[derive(Clone, Debug)]
pub struct YeeMaterialCoeffs {
    /// Diagonal inverse relative permittivity `W_{xx}, W_{yy}, W_{zz}` at `Ex, Ey, Ez` Yee edges.
    pub inv_eps_x: AlignedSlab,
    pub inv_eps_y: AlignedSlab,
    pub inv_eps_z: AlignedSlab,
    /// Electric conductivity `sigma_e` (S/m) at `Ex, Ey, Ez` Yee edges.
    pub sigma_ex: AlignedSlab,
    pub sigma_ey: AlignedSlab,
    pub sigma_ez: AlignedSlab,
    /// Dominant 1-based domain label at each `Ex, Ey, Ez` Yee edge (`0` for PEC/void), used to bind ADE poles.
    pub dom_ex: Vec<i32>,
    pub dom_ey: Vec<i32>,
    pub dom_ez: Vec<i32>,
    /// Subpixel volume fraction in `(0, 1]` of the dispersive domain bound to `dom_ex/ey/ez` on interface edges.
    pub disp_frac_x: Vec<f64>,
    pub disp_frac_y: Vec<f64>,
    pub disp_frac_z: Vec<f64>,
    /// Inverse relative permeability `1 / mu_r` (multiplied by Dey–Mittra face area ratios when active) at `Hx, Hy, Hz` faces.
    pub inv_mu_x: AlignedSlab,
    pub inv_mu_y: AlignedSlab,
    pub inv_mu_z: AlignedSlab,
    /// Dey–Mittra unblocked primal edge length fractions `l_eff / dl` in `[0, 1]` (`0.0` = clamped PEC edge, `1.0` = full dielectric).
    pub pec_len_x: AlignedSlab,
    pub pec_len_y: AlignedSlab,
    pub pec_len_z: AlignedSlab,
    /// Whether any edge has non-trivial Dey–Mittra PEC conformal fraction (`0 < l < 1` or interior `l == 0`).
    pub has_conformal_pec: bool,
    /// CFL reduction factor `<= 1.0` induced by stabilized Dey–Mittra cut faces.
    pub conformal_cfl_scale: f64,
    /// Sparse off-diagonal inverse-permittivity couplings at `Ex`, `Ey`, `Ez` edges (populated only where `|W_{ab}| > 1e-11`).
    pub offdiag_x: Vec<OffDiagEntry>,
    pub offdiag_y: Vec<OffDiagEntry>,
    pub offdiag_z: Vec<OffDiagEntry>,
    /// Minimum effective principal refractive index observed across all active dielectric cells (for CFL calculation).
    pub min_refractive_index: f64,
}

impl YeeMaterialCoeffs {
    #[inline]
    pub fn has_offdiag(&self) -> bool {
        !self.offdiag_x.is_empty() || !self.offdiag_y.is_empty() || !self.offdiag_z.is_empty()
    }
}

/// Geometry oracle queried during Stage 1 subpixel dielectric & Dey–Mittra PEC assembly.
pub enum FdtdGeometryOracle {
    Geom {
        ev: Arc<Evaluated>,
        region_to_domain: Vec<Option<i32>>,
        pec_domains: BTreeSet<i32>,
        outside_is_pec: bool,
        periodic: [bool; 3],
        bounds: ([f64; 3], [f64; 3]),
    },
    Simplex2d {
        bvh: MeshGridLocator2d,
        pec_domains: BTreeSet<i32>,
        outside_is_pec: bool,
        periodic: [bool; 3],
        bounds: ([f64; 3], [f64; 3]),
    },
    Simplex3d {
        bvh: MeshGridLocator3d,
        pec_domains: BTreeSet<i32>,
        outside_is_pec: bool,
        periodic: [bool; 3],
        bounds: ([f64; 3], [f64; 3]),
    },
}

/// Result of classifying a point `p = [x, y, z]` in wafer space.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointClass {
    Domain(i32),
    Pec,
    Outside,
}

impl FdtdGeometryOracle {
    pub fn from_geom(
        ev: Arc<Evaluated>,
        labels: &LabelMap,
        pec_domains: BTreeSet<i32>,
        outside_is_pec: bool,
        periodic: [bool; 3],
        bounds: ([f64; 3], [f64; 3]),
    ) -> Self {
        let nr = ev.scene().regions.len();
        let mut region_to_domain = vec![None; nr];
        if ev.complex.is_some() {
            for (ri, d) in labels.domains.iter().enumerate() {
                if !matches!(d.kind, RegionKind::Hole | RegionKind::Void) {
                    if let Some(reg) = ev.regions().get(ri) {
                        if reg.scene_index < nr {
                            region_to_domain[reg.scene_index] = Some(d.id);
                        }
                    }
                }
            }
        } else {
            for (i, slot) in region_to_domain.iter_mut().enumerate() {
                *slot = Some((i + 1) as i32);
            }
        }
        Self::Geom { ev, region_to_domain, pec_domains, outside_is_pec, periodic, bounds }
    }

    pub fn from_simplex(
        mesh: &SimplexMesh,
        pec_domains: BTreeSet<i32>,
        pec_facet_labels: &BTreeSet<i32>,
        outside_is_pec: bool,
        periodic: [bool; 3],
        bounds: ([f64; 3], [f64; 3]),
    ) -> Self {
        if mesh.dim == 2 {
            Self::Simplex2d {
                bvh: MeshGridLocator2d::new(mesh, pec_facet_labels),
                pec_domains,
                outside_is_pec,
                periodic,
                bounds,
            }
        } else {
            Self::Simplex3d {
                bvh: MeshGridLocator3d::new(mesh, pec_facet_labels),
                pec_domains,
                outside_is_pec,
                periodic,
                bounds,
            }
        }
    }

    #[inline]
    fn map_query_point(
        mut p: [f64; 3],
        periodic: [bool; 3],
        bounds: ([f64; 3], [f64; 3]),
        clamp_nonperiodic: bool,
    ) -> [f64; 3] {
        let (lo, hi) = bounds;
        for a in 0..3 {
            let span = hi[a] - lo[a];
            if span <= 1e-30 {
                continue;
            }
            if periodic[a] {
                p[a] = lo[a] + (p[a] - lo[a]).rem_euclid(span);
                let eps = 1e-11 * span;
                p[a] = p[a].clamp(lo[a] + eps, hi[a] - eps);
            } else if clamp_nonperiodic {
                let eps = 1e-11 * span;
                p[a] = p[a].clamp(lo[a] + eps, hi[a] - eps);
            }
        }
        p
    }

    /// Lower bound on distance from `p` to any material or PEC interface (up to `margin`), enabling O(1) bulk skipping.
    #[inline]
    pub fn homogeneous_clearance(&self, p: [f64; 3], margin: f64) -> f64 {
        match self {
            Self::Geom { ev, periodic, bounds, .. } => {
                let q = Self::map_query_point(p, *periodic, *bounds, false);
                ev.sdf_bound(q, margin).abs()
            }
            Self::Simplex2d { bvh, periodic, bounds, outside_is_pec, .. } => {
                let q = Self::map_query_point(p, *periodic, *bounds, !*outside_is_pec);
                bvh.interface_clearance([q[0], q[2]], margin)
            }
            Self::Simplex3d { bvh, periodic, bounds, outside_is_pec, .. } => {
                let q = Self::map_query_point(p, *periodic, *bounds, !*outside_is_pec);
                bvh.interface_clearance(q, margin)
            }
        }
    }

    /// Classifies `p = [x, y, z]` into a 1-based domain label, volumetric/exterior PEC, or outside background.
    #[inline]
    pub fn classify(&self, p: [f64; 3]) -> PointClass {
        match self {
            Self::Geom { ev, region_to_domain, pec_domains, outside_is_pec, periodic, bounds } => {
                let q = Self::map_query_point(p, *periodic, *bounds, false);
                match ev.classify(q) {
                    Some(rid) => {
                        let dom = if ev.complex.is_some() {
                            ev.regions()
                                .get(rid as usize)
                                .and_then(|r| region_to_domain.get(r.scene_index).copied().flatten())
                                .unwrap_or(rid as i32 + 1)
                        } else {
                            region_to_domain.get(rid as usize).copied().flatten().unwrap_or(rid as i32 + 1)
                        };
                        if pec_domains.contains(&dom) {
                            PointClass::Pec
                        } else {
                            PointClass::Domain(dom)
                        }
                    }
                    None => {
                        if *outside_is_pec {
                            PointClass::Pec
                        } else {
                            PointClass::Outside
                        }
                    }
                }
            }
            Self::Simplex2d { bvh, pec_domains, outside_is_pec, periodic, bounds } => {
                let q = Self::map_query_point(p, *periodic, *bounds, !*outside_is_pec);
                match bvh.locate_domain([q[0], q[2]]) {
                    Some(dom) => {
                        if pec_domains.contains(&dom) {
                            PointClass::Pec
                        } else {
                            PointClass::Domain(dom)
                        }
                    }
                    None => {
                        if *outside_is_pec {
                            PointClass::Pec
                        } else {
                            PointClass::Outside
                        }
                    }
                }
            }
            Self::Simplex3d { bvh, pec_domains, outside_is_pec, periodic, bounds } => {
                let q = Self::map_query_point(p, *periodic, *bounds, !*outside_is_pec);
                match bvh.locate_domain(q) {
                    Some(dom) => {
                        if pec_domains.contains(&dom) {
                            PointClass::Pec
                        } else {
                            PointClass::Domain(dom)
                        }
                    }
                    None => {
                        if *outside_is_pec {
                            PointClass::Pec
                        } else {
                            PointClass::Outside
                        }
                    }
                }
            }
        }
    }

    /// Checks whether the straight primal segment `a -> b` crosses an explicit zero-thickness PEC facet sheet.
    #[inline]
    pub fn segment_hits_pec_sheet(&self, a: [f64; 3], b: [f64; 3]) -> bool {
        match self {
            Self::Geom { .. } => false,
            Self::Simplex2d { bvh, .. } => bvh.segment_hits_pec([a[0], a[2]], [b[0], b[2]]),
            Self::Simplex3d { bvh, .. } => bvh.segment_hits_pec(a, b),
        }
    }
}

/// Assembles all Yee dielectric (`PolarizedAveraging` / `Volumetric` / `Staircased`) and Dey–Mittra `PECConformal`
/// material arrays deterministically via fixed-chunk parallel loops (`femcore::par`).
#[allow(clippy::too_many_arguments)]
pub fn assemble_yee_materials(
    grid: &YeeGrid,
    oracle: &FdtdGeometryOracle,
    media: &BTreeMap<i32, DomainMedium>,
    default_medium: DomainMedium,
    disp_domains: &BTreeSet<i32>,
    subpixel: SubpixelMode,
    enable_conformal_pec: bool,
    dm_min_edge_frac: f64,
) -> YeeMaterialCoeffs {
    let n_slab = grid.slab_len();
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let dm_cut = dm_min_edge_frac.clamp(0.01, 0.45);

    let medium_of = |cls: PointClass| -> (i32, DomainMedium) {
        match cls {
            PointClass::Domain(d) => {
                let m = media.get(&d).copied().unwrap_or(default_medium);
                (if m.is_pec { 0 } else { d }, m)
            }
            PointClass::Pec => (0, DomainMedium { is_pec: true, ..default_medium }),
            PointClass::Outside => (0, default_medium),
        }
    };

    // Subgrid quadrature nodes in [-0.5, 0.5] (multiples of 12 so 2:1, 3:1, and 4:1 graded-grid interfaces are exact!)
    let m_sub: usize = if subpixel == SubpixelMode::Staircased {
        1
    } else if is_2d {
        24
    } else {
        12
    };
    let sub_offsets: Vec<f64> = (0..m_sub).map(|k| -0.5 + (k as f64 + 0.5) / (m_sub as f64)).collect();

    // Assemble per-x-slab in deterministic parallel chunks
    let run_par = femcore::par::threads() > 1 && grid.nx >= 8;
    let nx_plus = grid.nx + 1;
    let ny_plus = grid.ny + 1;

    struct SlabOut {
        inv_ex: Vec<f64>,
        inv_ey: Vec<f64>,
        inv_ez: Vec<f64>,
        sig_ex: Vec<f64>,
        sig_ey: Vec<f64>,
        sig_ez: Vec<f64>,
        dom_ex: Vec<i32>,
        dom_ey: Vec<i32>,
        dom_ez: Vec<i32>,
        df_x: Vec<f64>,
        df_y: Vec<f64>,
        df_z: Vec<f64>,
        pec_lx: Vec<f64>,
        pec_ly: Vec<f64>,
        pec_lz: Vec<f64>,
        inv_mx: Vec<f64>,
        inv_my: Vec<f64>,
        inv_mz: Vec<f64>,
        off_x: Vec<OffDiagEntry>,
        off_y: Vec<OffDiagEntry>,
        off_z: Vec<OffDiagEntry>,
        has_conf_pec: bool,
        min_n: f64,
    }

    let slabs: Vec<SlabOut> = femcore::par::chunk_map(nx_plus, run_par, |ix| {
        let slab_sz = ny_plus * grid.nz_pad;
        let mut out = SlabOut {
            inv_ex: vec![1.0; slab_sz],
            inv_ey: vec![1.0; slab_sz],
            inv_ez: vec![1.0; slab_sz],
            sig_ex: vec![0.0; slab_sz],
            sig_ey: vec![0.0; slab_sz],
            sig_ez: vec![0.0; slab_sz],
            dom_ex: vec![0; slab_sz],
            dom_ey: vec![0; slab_sz],
            dom_ez: vec![0; slab_sz],
            df_x: vec![1.0; slab_sz],
            df_y: vec![1.0; slab_sz],
            df_z: vec![1.0; slab_sz],
            pec_lx: vec![1.0; slab_sz],
            pec_ly: vec![1.0; slab_sz],
            pec_lz: vec![1.0; slab_sz],
            inv_mx: vec![1.0; slab_sz],
            inv_my: vec![1.0; slab_sz],
            inv_mz: vec![1.0; slab_sz],
            off_x: Vec::new(),
            off_y: Vec::new(),
            off_z: Vec::new(),
            has_conf_pec: false,
            min_n: f64::INFINITY,
        };

        let eval_e_edge = |axis: usize,
                           p_mid: [f64; 3],
                           box_size: [f64; 3],
                           seg_a: [f64; 3],
                           seg_b: [f64; 3]|
         -> ([[f64; 3]; 3], f64, i32, f64, f64, f64) {
            // Check explicit zero-thickness PEC facet sheet first
            if oracle.segment_hits_pec_sheet(seg_a, seg_b) {
                return ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 0.0, 0, 0.0, 0.0, f64::INFINITY);
            }
            let half_diag = 0.5
                * (box_size[0] * box_size[0]
                    + if is_2d { 0.0 } else { box_size[1] * box_size[1] }
                    + box_size[2] * box_size[2])
                    .sqrt();

            let (dom_c, med_c) = medium_of(oracle.classify(p_mid));
            let is_homog = oracle.homogeneous_clearance(p_mid, half_diag * 1.05) >= half_diag;
            if is_homog || (subpixel == SubpixelMode::Staircased && !enable_conformal_pec) {
                if med_c.is_pec {
                    return (med_c.inv_eps_inf, 0.0, 0, 0.0, 0.0, f64::INFINITY);
                }
                let n_loc = min_principal_n(&med_c);
                return (med_c.inv_eps_inf, med_c.sigma_e[axis], dom_c, 1.0, 1.0, n_loc);
            }

            // Compute Dey–Mittra 1D edge open fraction along `seg_a -> seg_b` when conformal PEC is active
            let mut open_len_frac = 1.0;
            let mut first_open_med: Option<(i32, DomainMedium)> = None;
            if enable_conformal_pec {
                let n_line = 16usize;
                let mut open_cnt = 0usize;
                for s in 0..n_line {
                    let t = (s as f64 + 0.5) / (n_line as f64);
                    let q = [
                        seg_a[0] + t * (seg_b[0] - seg_a[0]),
                        seg_a[1] + t * (seg_b[1] - seg_a[1]),
                        seg_a[2] + t * (seg_b[2] - seg_a[2]),
                    ];
                    let (dq, mq) = medium_of(oracle.classify(q));
                    if !mq.is_pec {
                        open_cnt += 1;
                        if first_open_med.is_none() {
                            first_open_med = Some((dq, mq));
                        }
                    }
                }
                open_len_frac = (open_cnt as f64) / (n_line as f64);
                if open_len_frac < dm_cut {
                    return ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 0.0, 0, 0.0, 0.0, f64::INFINITY);
                }
                if subpixel == SubpixelMode::Staircased {
                    let (d_open, m_open) = if !med_c.is_pec {
                        (dom_c, med_c)
                    } else if let Some(pair) = first_open_med {
                        pair
                    } else {
                        return ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 0.0, 0, 0.0, 0.0, f64::INFINITY);
                    };
                    let n_loc = min_principal_n(&m_open);
                    return (m_open.inv_eps_inf, m_open.sigma_e[axis], d_open, 1.0, open_len_frac, n_loc);
                }
            } else if med_c.is_pec {
                return (med_c.inv_eps_inf, 0.0, 0, 0.0, 0.0, f64::INFINITY);
            }

            // Subpixel quadrature over dual control box around `p_mid`
            let mut eps_avg = [[0.0; 3]; 3];
            let mut inv_eps_avg = [[0.0; 3]; 3];
            let mut sig_sum = 0.0;
            let mut grad_trace = [0.0; 3];
            let mut dom_counts: BTreeMap<i32, usize> = BTreeMap::new();
            let mut n_diel = 0usize;
            let mut min_n_box = f64::INFINITY;

            let ny_sub = if is_2d { 1 } else { m_sub };
            for &ox in &sub_offsets {
                let x = p_mid[0] + ox * box_size[0];
                for iy_s in 0..ny_sub {
                    let oy = if is_2d { 0.0 } else { sub_offsets[iy_s] };
                    let y = p_mid[1] + oy * box_size[1];
                    for &oz in &sub_offsets {
                        let z = p_mid[2] + oz * box_size[2];
                        let (d_q, m_q) = medium_of(oracle.classify([x, y, z]));
                        if m_q.is_pec {
                            continue;
                        }
                        n_diel += 1;
                        *dom_counts.entry(d_q).or_insert(0) += 1;
                        min_n_box = min_n_box.min(min_principal_n(&m_q));
                        sig_sum += m_q.sigma_e[axis];
                        let tr = m_q.eps_inf[0][0] + m_q.eps_inf[1][1] + m_q.eps_inf[2][2];
                        grad_trace[0] += ox * tr;
                        grad_trace[1] += oy * tr;
                        grad_trace[2] += oz * tr;
                        for a in 0..3 {
                            for b in 0..3 {
                                eps_avg[a][b] += m_q.eps_inf[a][b];
                                inv_eps_avg[a][b] += m_q.inv_eps_inf[a][b];
                            }
                        }
                    }
                }
            }

            if n_diel == 0 {
                return ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 0.0, 0, 0.0, 0.0, f64::INFINITY);
            }
            let inv_n = 1.0 / (n_diel as f64);
            for a in 0..3 {
                for b in 0..3 {
                    eps_avg[a][b] *= inv_n;
                    inv_eps_avg[a][b] *= inv_n;
                }
            }
            let (dom_best, disp_frac) = if let Some((&d_disp, &c_disp)) =
                dom_counts.iter().filter(|&(d, _)| disp_domains.contains(d)).max_by_key(|&(_, &c)| c)
            {
                (d_disp, (c_disp as f64) * inv_n)
            } else {
                let db = dom_counts.into_iter().max_by_key(|&(_, c)| c).map_or(dom_c, |(d, _)| d);
                (db, 1.0)
            };
            let inv_of_avg = invert_sym3(&eps_avg).unwrap_or(inv_eps_avg);

            let w_eff = if subpixel == SubpixelMode::Volumetric {
                inv_of_avg
            } else {
                // Estimate unit interface normal n_hat from physical gradient of tr(eps): grad_x = grad_trace[k] / box_size[k]
                let gx = if box_size[0] > 0.0 { grad_trace[0] / box_size[0] } else { 0.0 };
                let gy = if !is_2d && box_size[1] > 0.0 { grad_trace[1] / box_size[1] } else { 0.0 };
                let gz = if box_size[2] > 0.0 { grad_trace[2] / box_size[2] } else { 0.0 };
                let gnorm = (gx * gx + gy * gy + gz * gz).sqrt();
                if gnorm <= 1e-11 {
                    inv_of_avg
                } else {
                    let n_hat = [gx / gnorm, gy / gnorm, gz / gnorm];
                    polarized_inverse_eps(&inv_eps_avg, &inv_of_avg, n_hat)
                }
            };

            (w_eff, sig_sum * inv_n, dom_best, disp_frac, open_len_frac, min_n_box)
        };

        let dual_mid_x = if ix == 0 || ix == grid.nx { grid.vx[ix] } else { 0.5 * (grid.cx[ix - 1] + grid.cx[ix]) };
        for iy in 0..ny_plus {
            let y_v = grid.vy[iy];
            let y_c = if iy < grid.ny { grid.cy[iy] } else { grid.vy[iy] };
            let dy_p = if iy < grid.ny { grid.dp_y[iy] } else { grid.dp_y[grid.ny - 1] };
            let dy_d = grid.dd_y[iy];
            let dual_mid_y = if is_2d {
                0.0
            } else if iy == 0 || iy == grid.ny {
                y_v
            } else {
                0.5 * (grid.cy[iy - 1] + grid.cy[iy])
            };

            for iz in 0..=grid.nz {
                let local = iy * grid.stride_y + iz;
                let g_idx = grid.idx(ix, iy, iz);
                let z_v = grid.vz[iz];
                let z_c = if iz < grid.nz { grid.cz[iz] } else { grid.vz[iz] };
                let dz_p = if iz < grid.nz { grid.dp_z[iz] } else { grid.dp_z[grid.nz - 1] };
                let dz_d = grid.dd_z[iz];
                let dual_mid_z = if iz == 0 || iz == grid.nz { z_v } else { 0.5 * (grid.cz[iz - 1] + grid.cz[iz]) };

                // Ex at (cx[ix], dual_mid_y, dual_mid_z) for ix < nx
                if ix < grid.nx && (!is_2d || iy == 0) {
                    let p = [grid.cx[ix], dual_mid_y, dual_mid_z];
                    let bsz = [grid.dp_x[ix], dy_d, dz_d];
                    let sa = [grid.vx[ix], if is_2d { 0.0 } else { y_v }, z_v];
                    let sb = [grid.vx[ix + 1], if is_2d { 0.0 } else { y_v }, z_v];
                    let (w, sig, dom, dfrac, lf, mn) = eval_e_edge(0, p, bsz, sa, sb);
                    out.inv_ex[local] = w[0][0];
                    out.sig_ex[local] = sig;
                    out.dom_ex[local] = dom;
                    out.df_x[local] = dfrac;
                    out.pec_lx[local] = lf;
                    out.min_n = out.min_n.min(mn);
                    if lf > 0.0 && lf < 1.0 {
                        out.has_conf_pec = true;
                    }
                    if lf > 0.0 && (w[0][1].abs() > 1e-11 || w[0][2].abs() > 1e-11) {
                        out.off_x.push(OffDiagEntry {
                            idx: g_idx,
                            ix: ix as u32,
                            iy: iy as u32,
                            iz: iz as u32,
                            w_ab: w[0][1],
                            w_ac: w[0][2],
                        });
                    }
                }

                // Ey at (dual_mid_x, cy[iy], dual_mid_z) for iy < ny
                if iy < grid.ny {
                    let p = [dual_mid_x, if is_2d { 0.0 } else { y_c }, dual_mid_z];
                    let bsz = [grid.dd_x[ix], dy_p, dz_d];
                    let sa = [grid.vx[ix], if is_2d { 0.0 } else { grid.vy[iy] }, z_v];
                    let sb = [grid.vx[ix], if is_2d { 0.0 } else { grid.vy[iy + 1] }, z_v];
                    let (w, sig, dom, dfrac, lf, mn) = eval_e_edge(1, p, bsz, sa, sb);
                    out.inv_ey[local] = w[1][1];
                    out.sig_ey[local] = sig;
                    out.dom_ey[local] = dom;
                    out.df_y[local] = dfrac;
                    out.pec_ly[local] = lf;
                    out.min_n = out.min_n.min(mn);
                    if lf > 0.0 && lf < 1.0 {
                        out.has_conf_pec = true;
                    }
                    if lf > 0.0 && (w[1][0].abs() > 1e-11 || w[1][2].abs() > 1e-11) {
                        out.off_y.push(OffDiagEntry {
                            idx: g_idx,
                            ix: ix as u32,
                            iy: iy as u32,
                            iz: iz as u32,
                            w_ab: w[1][0],
                            w_ac: w[1][2],
                        });
                    }
                }

                // Ez at (dual_mid_x, dual_mid_y, cz[iz]) for iz < nz
                if iz < grid.nz && (!is_2d || iy == 0) {
                    let p = [dual_mid_x, dual_mid_y, z_c];
                    let bsz = [grid.dd_x[ix], dy_d, dz_p];
                    let sa = [grid.vx[ix], p[1], grid.vz[iz]];
                    let sb = [grid.vx[ix], p[1], grid.vz[iz + 1]];
                    let (w, sig, dom, dfrac, lf, mn) = eval_e_edge(2, p, bsz, sa, sb);
                    out.inv_ez[local] = w[2][2];
                    out.sig_ez[local] = sig;
                    out.dom_ez[local] = dom;
                    out.df_z[local] = dfrac;
                    out.pec_lz[local] = lf;
                    out.min_n = out.min_n.min(mn);
                    if lf > 0.0 && lf < 1.0 {
                        out.has_conf_pec = true;
                    }
                    if lf > 0.0 && (w[2][0].abs() > 1e-11 || w[2][1].abs() > 1e-11) {
                        out.off_z.push(OffDiagEntry {
                            idx: g_idx,
                            ix: ix as u32,
                            iy: iy as u32,
                            iz: iz as u32,
                            w_ab: w[2][0],
                            w_ac: w[2][1],
                        });
                    }
                }

                // Magnetic face inverse permeabilities (and Dey–Mittra open face area ratios)
                // Hx at (vx[ix], cy[iy], cz[iz]) for iy < ny, iz < nz
                if iy < grid.ny && iz < grid.nz {
                    let p = [grid.vx[ix], if is_2d { 0.0 } else { y_c }, z_c];
                    let (_, m) = medium_of(oracle.classify(p));
                    let mut im = 1.0 / m.mu_r[0];
                    if !is_2d
                        && enable_conformal_pec
                        && oracle.homogeneous_clearance(p, 0.6 * dy_p.hypot(dz_p)) < 0.55 * dy_p.hypot(dz_p)
                    {
                        let af = sample_face_open_fraction(
                            oracle,
                            [grid.vx[ix], grid.vy[iy], grid.vz[iz]],
                            [0.0, dy_p, 0.0],
                            [0.0, 0.0, dz_p],
                        );
                        if af < dm_cut * dm_cut {
                            im = 0.0;
                        } else if af < 1.0 {
                            out.has_conf_pec = true;
                            im /= af.max(0.25);
                        }
                    } else if m.is_pec {
                        im = 0.0;
                    }
                    out.inv_mx[local] = im;
                }
                // Hy at (cx[ix], vy[iy], cz[iz]) for ix < nx, iz < nz
                if ix < grid.nx && iz < grid.nz && (!is_2d || iy == 0) {
                    let p = [grid.cx[ix], if is_2d { 0.0 } else { y_v }, z_c];
                    let (_, m) = medium_of(oracle.classify(p));
                    let mut im = 1.0 / m.mu_r[1];
                    if enable_conformal_pec
                        && oracle.homogeneous_clearance(p, 0.6 * (grid.dp_x[ix].hypot(dz_p)))
                            < 0.55 * (grid.dp_x[ix].hypot(dz_p))
                    {
                        let af = sample_face_open_fraction(
                            oracle,
                            [grid.vx[ix], p[1], grid.vz[iz]],
                            [grid.dp_x[ix], 0.0, 0.0],
                            [0.0, 0.0, dz_p],
                        );
                        if af < dm_cut * dm_cut {
                            im = 0.0;
                        } else if af < 1.0 {
                            out.has_conf_pec = true;
                            im /= af.max(0.25);
                        }
                    } else if m.is_pec {
                        im = 0.0;
                    }
                    out.inv_my[local] = im;
                }
                // Hz at (cx[ix], cy[iy], vz[iz]) for ix < nx, iy < ny
                if ix < grid.nx && iy < grid.ny {
                    let p = [grid.cx[ix], if is_2d { 0.0 } else { y_c }, z_v];
                    let (_, m) = medium_of(oracle.classify(p));
                    let mut im = 1.0 / m.mu_r[2];
                    if !is_2d
                        && enable_conformal_pec
                        && oracle.homogeneous_clearance(p, 0.6 * grid.dp_x[ix].hypot(dy_p))
                            < 0.55 * grid.dp_x[ix].hypot(dy_p)
                    {
                        let af = sample_face_open_fraction(
                            oracle,
                            [grid.vx[ix], grid.vy[iy], z_v],
                            [grid.dp_x[ix], 0.0, 0.0],
                            [0.0, dy_p, 0.0],
                        );
                        if af < dm_cut * dm_cut {
                            im = 0.0;
                        } else if af < 1.0 {
                            out.has_conf_pec = true;
                            im /= af.max(0.25);
                        }
                    } else if m.is_pec {
                        im = 0.0;
                    }
                    out.inv_mz[local] = im;
                }
            }
        }
        out
    });

    let mut coeffs = YeeMaterialCoeffs {
        inv_eps_x: AlignedSlab::zeros(n_slab),
        inv_eps_y: AlignedSlab::zeros(n_slab),
        inv_eps_z: AlignedSlab::zeros(n_slab),
        sigma_ex: AlignedSlab::zeros(n_slab),
        sigma_ey: AlignedSlab::zeros(n_slab),
        sigma_ez: AlignedSlab::zeros(n_slab),
        dom_ex: vec![0; n_slab],
        dom_ey: vec![0; n_slab],
        dom_ez: vec![0; n_slab],
        disp_frac_x: vec![1.0; n_slab],
        disp_frac_y: vec![1.0; n_slab],
        disp_frac_z: vec![1.0; n_slab],
        inv_mu_x: AlignedSlab::zeros(n_slab),
        inv_mu_y: AlignedSlab::zeros(n_slab),
        inv_mu_z: AlignedSlab::zeros(n_slab),
        pec_len_x: AlignedSlab::zeros(n_slab),
        pec_len_y: AlignedSlab::zeros(n_slab),
        pec_len_z: AlignedSlab::zeros(n_slab),
        has_conformal_pec: false,
        conformal_cfl_scale: 1.0,
        offdiag_x: Vec::new(),
        offdiag_y: Vec::new(),
        offdiag_z: Vec::new(),
        min_refractive_index: f64::INFINITY,
    };

    let slab_sz = ny_plus * grid.nz_pad;
    for (ix, s) in slabs.into_iter().enumerate() {
        let base = ix * grid.stride_x;
        coeffs.inv_eps_x.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.inv_ex);
        coeffs.inv_eps_y.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.inv_ey);
        coeffs.inv_eps_z.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.inv_ez);
        coeffs.sigma_ex.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.sig_ex);
        coeffs.sigma_ey.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.sig_ey);
        coeffs.sigma_ez.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.sig_ez);
        coeffs.dom_ex[base..base + slab_sz].copy_from_slice(&s.dom_ex);
        coeffs.dom_ey[base..base + slab_sz].copy_from_slice(&s.dom_ey);
        coeffs.dom_ez[base..base + slab_sz].copy_from_slice(&s.dom_ez);
        coeffs.disp_frac_x[base..base + slab_sz].copy_from_slice(&s.df_x);
        coeffs.disp_frac_y[base..base + slab_sz].copy_from_slice(&s.df_y);
        coeffs.disp_frac_z[base..base + slab_sz].copy_from_slice(&s.df_z);
        coeffs.pec_len_x.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.pec_lx);
        coeffs.pec_len_y.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.pec_ly);
        coeffs.pec_len_z.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.pec_lz);
        coeffs.inv_mu_x.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.inv_mx);
        coeffs.inv_mu_y.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.inv_my);
        coeffs.inv_mu_z.as_mut_slice()[base..base + slab_sz].copy_from_slice(&s.inv_mz);
        coeffs.offdiag_x.extend(s.off_x);
        coeffs.offdiag_y.extend(s.off_y);
        coeffs.offdiag_z.extend(s.off_z);
        coeffs.has_conformal_pec |= s.has_conf_pec;
        coeffs.min_refractive_index = coeffs.min_refractive_index.min(s.min_n);
    }

    if !coeffs.min_refractive_index.is_finite() {
        coeffs.min_refractive_index = 1.0;
    }
    if coeffs.has_conformal_pec {
        // Dey–Mittra topological consistency: any primal E edge bordering a dead PEC magnetic face (`inv_mu == 0`,
        // where `H == 0` permanently) must also be clamped to PEC (`pec_len = 0, inv_eps = 0`), otherwise a zeroed
        // magnetic cell would act as a spurious PMC (`H_t = 0`) boundary instead of PEC!
        let imx = coeffs.inv_mu_x.as_slice();
        let imy = coeffs.inv_mu_y.as_slice();
        let imz = coeffs.inv_mu_z.as_slice();
        for ix in 0..=grid.nx {
            for iy in 0..ny_plus {
                for iz in 0..=grid.nz {
                    let g = grid.idx(ix, iy, iz);
                    if ix < grid.nx && iz > 0 && iz < grid.nz && (is_2d || (iy > 0 && iy < grid.ny)) {
                        let dead_y = imy[g] <= 0.0 || imy[g - 1] <= 0.0;
                        let dead_z = !is_2d && (imz[g] <= 0.0 || imz[g - grid.stride_y] <= 0.0);
                        if dead_y || dead_z {
                            coeffs.pec_len_x.as_mut_slice()[g] = 0.0;
                            coeffs.inv_eps_x.as_mut_slice()[g] = 0.0;
                        }
                    }
                    if !is_2d && iy < grid.ny && ix > 0 && ix < grid.nx && iz > 0 && iz < grid.nz {
                        let dead_x = imx[g] <= 0.0 || imx[g - 1] <= 0.0;
                        let dead_z = imz[g] <= 0.0 || imz[g - grid.stride_x] <= 0.0;
                        if dead_x || dead_z {
                            coeffs.pec_len_y.as_mut_slice()[g] = 0.0;
                            coeffs.inv_eps_y.as_mut_slice()[g] = 0.0;
                        }
                    }
                    if iz < grid.nz && ix > 0 && ix < grid.nx && (is_2d || (iy > 0 && iy < grid.ny)) {
                        let dead_y = imy[g] <= 0.0 || imy[g - grid.stride_x] <= 0.0;
                        let dead_x = !is_2d && (imx[g] <= 0.0 || imx[g - grid.stride_y] <= 0.0);
                        if dead_y || dead_x {
                            coeffs.pec_len_z.as_mut_slice()[g] = 0.0;
                            coeffs.inv_eps_z.as_mut_slice()[g] = 0.0;
                        }
                    }
                }
            }
        }
        coeffs.conformal_cfl_scale = 0.65;
    }
    coeffs
}

fn sample_face_open_fraction(oracle: &FdtdGeometryOracle, origin: [f64; 3], u_vec: [f64; 3], v_vec: [f64; 3]) -> f64 {
    let n = 6usize;
    let mut open = 0usize;
    for i in 0..n {
        let su = (i as f64 + 0.5) / (n as f64);
        for j in 0..n {
            let sv = (j as f64 + 0.5) / (n as f64);
            let p = [
                origin[0] + su * u_vec[0] + sv * v_vec[0],
                origin[1] + su * u_vec[1] + sv * v_vec[1],
                origin[2] + su * u_vec[2] + sv * v_vec[2],
            ];
            if oracle.classify(p) != PointClass::Pec {
                open += 1;
            }
        }
    }
    (open as f64) / ((n * n) as f64)
}

/// Computes `W = P <eps^-1> P + (I - P) <eps>^-1 (I - P)` where `P = n n^T`.
#[inline]
pub fn polarized_inverse_eps(
    avg_inv_eps: &[[f64; 3]; 3],
    inv_avg_eps: &[[f64; 3]; 3],
    n_hat: [f64; 3],
) -> [[f64; 3]; 3] {
    let mut p = [[0.0; 3]; 3];
    let mut q = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] = n_hat[i] * n_hat[j];
            q[i][j] = (if i == j { 1.0 } else { 0.0 }) - p[i][j];
        }
    }
    let pap = mat3_triple(&p, avg_inv_eps, &p);
    let qbq = mat3_triple(&q, inv_avg_eps, &q);
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = pap[i][j] + qbq[i][j];
        }
    }
    out
}

#[inline]
fn mat3_triple(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3], c_mat: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut ab = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            ab[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = ab[i][0] * c_mat[0][j] + ab[i][1] * c_mat[1][j] + ab[i][2] * c_mat[2][j];
        }
    }
    out
}

#[inline]
pub fn invert_sym3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    if !(det.is_finite() && det.abs() > 1e-30) {
        return None;
    }
    let inv_d = 1.0 / det;
    let c11 = m[0][0] * m[2][2] - m[0][2] * m[2][0];
    let c12 = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]);
    let c22 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    Some([
        [c00 * inv_d, c01 * inv_d, c02 * inv_d],
        [c01 * inv_d, c11 * inv_d, c12 * inv_d],
        [c02 * inv_d, c12 * inv_d, c22 * inv_d],
    ])
}

#[inline]
fn min_principal_n(m: &DomainMedium) -> f64 {
    let d_min = m.eps_inf[0][0].min(m.eps_inf[1][1]).min(m.eps_inf[2][2]).max(1.0);
    let mu_min = m.mu_r[0].min(m.mu_r[1]).min(m.mu_r[2]).max(1.0);
    (d_min * mu_min).sqrt()
}

// ------------------------------------------------------------------ Fast Binned Mesh Point/Interface Locators (2D & 3D)

pub struct MeshGridLocator2d {
    lo: [f64; 2],
    size: [f64; 2],
    dims: [usize; 2],
    bins: Vec<Vec<u32>>,
    tris: Vec<([f64; 2], [f64; 2], [f64; 2], i32)>,
    interface_segs: Vec<([f64; 2], [f64; 2], bool)>,
}

impl MeshGridLocator2d {
    pub fn new(mesh: &SimplexMesh, pec_facet_labels: &BTreeSet<i32>) -> Self {
        let (mlo, mhi) = mesh.bbox();
        let lo = [mlo[0], mlo[1]];
        let ext = [(mhi[0] - mlo[0]).max(1e-30), (mhi[1] - mlo[1]).max(1e-30)];
        let s = (ext[0] * ext[1] / (mesh.n_cells().max(1) as f64)).sqrt().max(1e-30);
        let dims = [((ext[0] / s).ceil() as usize).clamp(1, 256), ((ext[1] / s).ceil() as usize).clamp(1, 256)];
        let size = [ext[0] / (dims[0] as f64), ext[1] / (dims[1] as f64)];
        let mut bins = vec![Vec::new(); dims[0] * dims[1]];
        let mut tris = Vec::with_capacity(mesh.n_cells());
        let mut edge_doms: BTreeMap<(u32, u32), Vec<i32>> = BTreeMap::new();

        for ci in 0..mesh.n_cells() {
            let cv = mesh.cell(ci);
            let p0 = [mesh.point(cv[0] as usize)[0], mesh.point(cv[0] as usize)[1]];
            let p1 = [mesh.point(cv[1] as usize)[0], mesh.point(cv[1] as usize)[1]];
            let p2 = [mesh.point(cv[2] as usize)[0], mesh.point(cv[2] as usize)[1]];
            let lab = mesh.cell_labels[ci];
            tris.push((p0, p1, p2, lab));
            for (u, w) in [(cv[0], cv[1]), (cv[1], cv[2]), (cv[2], cv[0])] {
                edge_doms.entry((u.min(w), u.max(w))).or_default().push(lab);
            }
            let bx0 = (((p0[0].min(p1[0]).min(p2[0]) - lo[0]) / size[0]).floor() as isize)
                .clamp(0, dims[0] as isize - 1) as usize;
            let bx1 = (((p0[0].max(p1[0]).max(p2[0]) - lo[0]) / size[0]).floor() as isize)
                .clamp(0, dims[0] as isize - 1) as usize;
            let bz0 = (((p0[1].min(p1[1]).min(p2[1]) - lo[1]) / size[1]).floor() as isize)
                .clamp(0, dims[1] as isize - 1) as usize;
            let bz1 = (((p0[1].max(p1[1]).max(p2[1]) - lo[1]) / size[1]).floor() as isize)
                .clamp(0, dims[1] as isize - 1) as usize;
            for ix in bx0..=bx1 {
                for iz in bz0..=bz1 {
                    bins[ix * dims[1] + iz].push(ci as u32);
                }
            }
        }

        let mut pec_edges = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if pec_facet_labels.contains(&mesh.facet_labels[f]) {
                let fv = mesh.facet(f);
                pec_edges.insert((fv[0].min(fv[1]), fv[0].max(fv[1])));
            }
        }

        let mut interface_segs = Vec::new();
        for ((u, w), doms) in edge_doms {
            let is_pec = pec_edges.contains(&(u, w));
            let is_iface = (doms.len() == 2 && doms[0] != doms[1]) || is_pec;
            if is_iface {
                let pa = [mesh.point(u as usize)[0], mesh.point(u as usize)[1]];
                let pb = [mesh.point(w as usize)[0], mesh.point(w as usize)[1]];
                interface_segs.push((pa, pb, is_pec));
            }
        }

        Self { lo, size, dims, bins, tris, interface_segs }
    }

    pub fn locate_domain(&self, q: [f64; 2]) -> Option<i32> {
        let ix = ((q[0] - self.lo[0]) / self.size[0]).floor() as isize;
        let iz = ((q[1] - self.lo[1]) / self.size[1]).floor() as isize;
        if ix < 0 || iz < 0 || ix >= self.dims[0] as isize || iz >= self.dims[1] as isize {
            return None;
        }
        for &ci in &self.bins[(ix as usize) * self.dims[1] + (iz as usize)] {
            let (p0, p1, p2, lab) = self.tris[ci as usize];
            let det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
            if det.abs() <= 1e-30 {
                continue;
            }
            let l1 = ((q[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (q[1] - p0[1])) / det;
            let l2 = ((p1[0] - p0[0]) * (q[1] - p0[1]) - (q[0] - p0[0]) * (p1[1] - p0[1])) / det;
            let l0 = 1.0 - l1 - l2;
            if l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9 {
                return Some(lab);
            }
        }
        None
    }

    pub fn interface_clearance(&self, q: [f64; 2], margin: f64) -> f64 {
        let mut best2 = margin * margin;
        for &(a, b, _) in &self.interface_segs {
            let lo_x = a[0].min(b[0]) - margin;
            let hi_x = a[0].max(b[0]) + margin;
            let lo_z = a[1].min(b[1]) - margin;
            let hi_z = a[1].max(b[1]) + margin;
            if q[0] < lo_x || q[0] > hi_x || q[1] < lo_z || q[1] > hi_z {
                continue;
            }
            let d2 = point_seg_dist2_2d(q, a, b);
            if d2 < best2 {
                best2 = d2;
            }
        }
        best2.sqrt()
    }

    pub fn segment_hits_pec(&self, a: [f64; 2], b: [f64; 2]) -> bool {
        for &(pa, pb, is_pec) in &self.interface_segs {
            if is_pec && segments_intersect_2d(a, b, pa, pb) {
                return true;
            }
        }
        false
    }
}

fn point_seg_dist2_2d(q: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (dx, dz) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dz * dz;
    if l2 <= 1e-30 {
        return (q[0] - a[0]).powi(2) + (q[1] - a[1]).powi(2);
    }
    let t = (((q[0] - a[0]) * dx + (q[1] - a[1]) * dz) / l2).clamp(0.0, 1.0);
    (q[0] - (a[0] + t * dx)).powi(2) + (q[1] - (a[1] + t * dz)).powi(2)
}

fn segments_intersect_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let cross = |o: [f64; 2], p: [f64; 2], q: [f64; 2]| (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0]);
    let d1 = cross(c, d, a);
    let d2 = cross(c, d, b);
    let d3 = cross(a, b, c);
    let d4 = cross(a, b, d);
    ((d1 <= 0.0 && d2 >= 0.0) || (d1 >= 0.0 && d2 <= 0.0)) && ((d3 <= 0.0 && d4 >= 0.0) || (d3 >= 0.0 && d4 <= 0.0))
}

pub struct MeshGridLocator3d {
    lo: [f64; 3],
    size: [f64; 3],
    dims: [usize; 3],
    bins: Vec<Vec<u32>>,
    tets: Vec<([f64; 3], [f64; 3], [f64; 3], [f64; 3], Option<Box<CurvedTet>>, i32)>,
    interface_tris: Vec<([f64; 3], [f64; 3], [f64; 3], [f64; 3], [f64; 3], bool)>,
}

impl MeshGridLocator3d {
    pub fn new(mesh: &SimplexMesh, pec_facet_labels: &BTreeSet<i32>) -> Self {
        let (mlo, mhi) = mesh.bbox();
        let lo = [mlo[0], mlo[1], mlo[2]];
        let ext = [(mhi[0] - mlo[0]).max(1e-30), (mhi[1] - mlo[1]).max(1e-30), (mhi[2] - mlo[2]).max(1e-30)];
        let s = (ext[0] * ext[1] * ext[2] / (mesh.n_cells().max(1) as f64)).cbrt().max(1e-30);
        let dims = [0, 1, 2].map(|k| ((ext[k] / s).ceil() as usize).clamp(1, 96));
        let size = [0, 1, 2].map(|k| ext[k] / (dims[k] as f64));
        let mut bins = vec![Vec::new(); dims[0] * dims[1] * dims[2]];
        let mut tets = Vec::with_capacity(mesh.n_cells());
        let mut face_doms: BTreeMap<[u32; 3], Vec<i32>> = BTreeMap::new();
        const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

        for ci in 0..mesh.n_cells() {
            let cv = mesh.cell(ci);
            let pt = |u: u32| {
                let p = mesh.point(u as usize);
                [p[0], p[1], p[2]]
            };
            let (p0, p1, p2, p3) = (pt(cv[0]), pt(cv[1]), pt(cv[2]), pt(cv[3]));
            let curved = mesh.cell_is_curved(ci).then(|| Box::new(CurvedTet::new(mesh.cell_geometry_nodes(ci))));
            let lab = mesh.cell_labels[ci];
            for f in FACES {
                let mut fk = [cv[f[0]], cv[f[1]], cv[f[2]]];
                fk.sort_unstable();
                face_doms.entry(fk).or_default().push(lab);
            }
            let mut clo = p0;
            let mut chi = p0;
            for q in [p1, p2, p3] {
                for k in 0..3 {
                    clo[k] = clo[k].min(q[k]);
                    chi[k] = chi[k].max(q[k]);
                }
            }
            if curved.is_some() {
                let pad = 0.15 * ((chi[0] - clo[0]).max(chi[1] - clo[1]).max(chi[2] - clo[2]));
                for k in 0..3 {
                    clo[k] -= pad;
                    chi[k] += pad;
                }
            }
            let b0 = [0, 1, 2]
                .map(|k| (((clo[k] - lo[k]) / size[k]).floor() as isize).clamp(0, dims[k] as isize - 1) as usize);
            let b1 = [0, 1, 2]
                .map(|k| (((chi[k] - lo[k]) / size[k]).floor() as isize).clamp(0, dims[k] as isize - 1) as usize);
            for ix in b0[0]..=b1[0] {
                for iy in b0[1]..=b1[1] {
                    for iz in b0[2]..=b1[2] {
                        bins[(ix * dims[1] + iy) * dims[2] + iz].push(ci as u32);
                    }
                }
            }
            tets.push((p0, p1, p2, p3, curved, lab));
        }

        let mut pec_faces = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if pec_facet_labels.contains(&mesh.facet_labels[f]) {
                let fv = mesh.facet(f);
                let mut fk = [fv[0], fv[1], fv[2]];
                fk.sort_unstable();
                pec_faces.insert(fk);
            }
        }

        let mut interface_tris = Vec::new();
        for (fk, doms) in face_doms {
            let is_pec = pec_faces.contains(&fk);
            let is_iface = (doms.len() == 2 && doms[0] != doms[1]) || is_pec;
            if is_iface {
                let pt = |u: u32| {
                    let p = mesh.point(u as usize);
                    [p[0], p[1], p[2]]
                };
                let (a, b, c) = (pt(fk[0]), pt(fk[1]), pt(fk[2]));
                let tlo = [a[0].min(b[0]).min(c[0]), a[1].min(b[1]).min(c[1]), a[2].min(b[2]).min(c[2])];
                let thi = [a[0].max(b[0]).max(c[0]), a[1].max(b[1]).max(c[1]), a[2].max(b[2]).max(c[2])];
                interface_tris.push((a, b, c, tlo, thi, is_pec));
            }
        }

        Self { lo, size, dims, bins, tets, interface_tris }
    }

    pub fn locate_domain(&self, q: [f64; 3]) -> Option<i32> {
        let ix = ((q[0] - self.lo[0]) / self.size[0]).floor() as isize;
        let iy = ((q[1] - self.lo[1]) / self.size[1]).floor() as isize;
        let iz = ((q[2] - self.lo[2]) / self.size[2]).floor() as isize;
        if ix < 0
            || iy < 0
            || iz < 0
            || ix >= self.dims[0] as isize
            || iy >= self.dims[1] as isize
            || iz >= self.dims[2] as isize
        {
            return None;
        }
        let bin = &self.bins[((ix as usize) * self.dims[1] + (iy as usize)) * self.dims[2] + (iz as usize)];
        for &ci in bin {
            let (p0, p1, p2, p3, ref curved, lab) = self.tets[ci as usize];
            if let Some(ct) = curved {
                if point_in_curved_tet(ct, p0, p1, p2, p3, q) {
                    return Some(lab);
                }
            } else if let Some(lam) = tet_bary(p0, p1, p2, p3, q) {
                if lam.iter().all(|&l| l >= -1e-9) {
                    return Some(lab);
                }
            }
        }
        None
    }

    pub fn interface_clearance(&self, q: [f64; 3], margin: f64) -> f64 {
        let mut best2 = margin * margin;
        for &(a, b, c, tlo, thi, _) in &self.interface_tris {
            if q[0] < tlo[0] - margin
                || q[0] > thi[0] + margin
                || q[1] < tlo[1] - margin
                || q[1] > thi[1] + margin
                || q[2] < tlo[2] - margin
                || q[2] > thi[2] + margin
            {
                continue;
            }
            let d2 = point_tri_dist2_3d(q, a, b, c);
            if d2 < best2 {
                best2 = d2;
            }
        }
        best2.sqrt()
    }

    pub fn segment_hits_pec(&self, a: [f64; 3], b: [f64; 3]) -> bool {
        let slo = [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])];
        let shi = [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])];
        for &(p0, p1, p2, tlo, thi, is_pec) in &self.interface_tris {
            if !is_pec
                || shi[0] < tlo[0] - 1e-12
                || slo[0] > thi[0] + 1e-12
                || shi[1] < tlo[1] - 1e-12
                || slo[1] > thi[1] + 1e-12
                || shi[2] < tlo[2] - 1e-12
                || slo[2] > thi[2] + 1e-12
            {
                continue;
            }
            if seg_hits_tri_3d(a, b, p0, p1, p2) {
                return true;
            }
        }
        false
    }
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn tet_bary(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], q: [f64; 3]) -> Option<[f64; 4]> {
    let det6 = |a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]| dot3(sub3(b, a), cross3(sub3(c, a), sub3(d, a)));
    let v = det6(p0, p1, p2, p3);
    if v.abs() <= 1e-300 {
        return None;
    }
    Some([det6(q, p1, p2, p3) / v, det6(p0, q, p2, p3) / v, det6(p0, p1, q, p3) / v, det6(p0, p1, p2, q) / v])
}

fn point_in_curved_tet(ct: &CurvedTet, p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], q: [f64; 3]) -> bool {
    let Some(mut lam) = tet_bary(p0, p1, p2, p3, q) else { return false };
    if lam.iter().any(|&l| l < -0.35 || l > 1.35) {
        return false;
    }
    for _ in 0..8 {
        let pos = ct.position(&lam);
        let res = sub3(q, pos);
        if dot3(res, res) < 1e-24 {
            break;
        }
        let d = ct.dx_dlam(&lam);
        let j = [sub3(d[0], d[3]), sub3(d[1], d[3]), sub3(d[2], d[3])];
        let c0 = cross3(j[1], j[2]);
        let det = dot3(j[0], c0);
        if det.abs() < 1e-30 {
            break;
        }
        let c1 = cross3(j[2], j[0]);
        let c2 = cross3(j[0], j[1]);
        let dl0 = dot3(c0, res) / det;
        let dl1 = dot3(c1, res) / det;
        let dl2 = dot3(c2, res) / det;
        lam[0] += dl0;
        lam[1] += dl1;
        lam[2] += dl2;
        lam[3] = 1.0 - lam[0] - lam[1] - lam[2];
    }
    lam.iter().all(|&l| l >= -1e-7)
}

fn point_tri_dist2_3d(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ap = sub3(p, a);
    let d1 = dot3(ab, ap);
    let d2 = dot3(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return dot3(ap, ap);
    }
    let bp = sub3(p, b);
    let d3 = dot3(ab, bp);
    let d4 = dot3(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return dot3(bp, bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let q = [a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]];
        let d = sub3(p, q);
        return dot3(d, d);
    }
    let cp = sub3(p, c);
    let d5 = dot3(ab, cp);
    let d6 = dot3(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return dot3(cp, cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let q = [a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]];
        let d = sub3(p, q);
        return dot3(d, d);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let bc = sub3(c, b);
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let q = [b[0] + w * bc[0], b[1] + w * bc[1], b[2] + w * bc[2]];
        let d = sub3(p, q);
        return dot3(d, d);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let q = [a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w];
    let d = sub3(p, q);
    dot3(d, d)
}

fn seg_hits_tri_3d(a: [f64; 3], b: [f64; 3], p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> bool {
    let dir = sub3(b, a);
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    let h = cross3(dir, e2);
    let det = dot3(e1, h);
    if det.abs() < 1e-24 {
        return false;
    }
    let inv_det = 1.0 / det;
    let s = sub3(a, p0);
    let u = inv_det * dot3(s, h);
    if !(-1e-7..=1.0 + 1e-7).contains(&u) {
        return false;
    }
    let q = cross3(s, e1);
    let v = inv_det * dot3(dir, q);
    if v < -1e-7 || u + v > 1.0 + 1e-7 {
        return false;
    }
    let t = inv_det * dot3(e2, q);
    (-1e-7..=1.0 + 1e-7).contains(&t)
}
