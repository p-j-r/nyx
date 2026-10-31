//! Use case: Bragg gratings in a 300 nm SOI wire, symmetric and asymmetric.
//!
//! A silicon wire (`soi::SoiWire` cross-section) whose side walls carry rectangular teeth of depth `tooth` with period
//! `period` and duty cycle `duty` along the propagation axis z. One period is meshed as a 3D cell: the 2D cross-section of the
//! widest part is extruded through two layers (teeth present, teeth absent), and the domain numbers of the layers pick the material
//! of each tooth, so the two ends of the cell carry identical meshes and a Floquet-periodic pair (`Periodic`) closes it. All other
//! walls are perfect conductors, far enough from the wire for the guided field to have died out.
//!
//! * [`GratingKind::Symmetric`]: teeth on both sides, in phase (mirror symmetric about the plane x = 0). The fundamental TE mode couples
//!   to its own backward wave. Its standing waves are even and odd about the tooth plane and stay quasi-TE.
//! * [`GratingKind::Asymmetric`]: teeth on one side only. The mirror symmetry is gone: the band-edge modes are hybrid, with a
//!   longitudinal part and a TM share, and their polarization on the Poincare sphere moves away from the x pole.
//! * [`GratingKind::HalfPeriodShifted`]: teeth on both sides, the right ones shifted by half a period (antisymmetric corrugation). It
//!   does not couple TE0 to itself at first order; it couples it to the odd-symmetry mode (a contra-directional coupler).
//!
//! The Floquet phase of `pi` over one period is the band edge, where the forward and backward waves of the unperturbed wire have the same
//! frequency and the grating splits them: the two lowest eigenfrequencies are the edges of the stop band. Their gap gives the coupling
//! strength `kappa = pi n_g (f_high - f_low) / c`, with `n_g` the group index of the wire; a coupled-mode estimate
//! `kappa ~ 2 |n_wide - n_narrow| sin(pi duty) / lambda` from the two effective indices of the unperturbed sections is the check.
//!
//! The mesh is in the model frame with the propagation axis along z, x lateral and y the surface normal: the polarization analysis in the
//! result uses exactly these axes.

use super::soi::{Cladding, SoiWire};
use crate::mesh::SimplexMesh;
use crate::mesher::{extrude, mesh_2d, Layer, Part, Shape, Spec2d};
use crate::problem::*;

const C0: f64 = 299_792_458.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GratingKind {
    Symmetric,
    Asymmetric,
    HalfPeriodShifted,
    /// no teeth at all: the unperturbed wire, the reference for the band edge (`f = c / (2 period n_eff)`)
    Uniform,
}

impl GratingKind {
    pub fn name(self) -> &'static str {
        match self {
            GratingKind::Symmetric => "symmetric (teeth on both sides)",
            GratingKind::Asymmetric => "asymmetric (teeth on one side)",
            GratingKind::HalfPeriodShifted => "antisymmetric (right teeth shifted by half a period)",
            GratingKind::Uniform => "uniform wire (no teeth)",
        }
    }
}

#[derive(Clone, Debug)]
pub struct BraggGrating {
    pub kind: GratingKind,
    pub wavelength: f64,
    /// width of the wire between the teeth
    pub width: f64,
    pub thickness: f64,
    /// depth of a tooth
    pub tooth: f64,
    pub period: f64,
    /// fraction of the period that carries teeth
    pub duty: f64,
    pub n_si: f64,
    pub n_sio2: f64,
    pub cladding: Cladding,
    /// oxide around the widest part of the wire, laterally and above and below
    pub margin: f64,
    /// Floquet phase over one period (pi: the band edge)
    pub phase: f64,
    pub n_modes: usize,
    pub order: usize,
    /// element size in the silicon and in the teeth, near the wire, far from it
    pub h_core: f64,
    pub h_near: f64,
    pub h_far: f64,
    /// element layers along the propagation axis in the part with teeth and in the part without
    pub layers: (usize, usize),
}

impl Default for BraggGrating {
    fn default() -> Self {
        BraggGrating {
            kind: GratingKind::Symmetric,
            wavelength: 1550e-9,
            width: 500e-9,
            thickness: 300e-9,
            tooth: 80e-9,
            period: 300e-9,
            duty: 0.5,
            n_si: 3.4757,
            n_sio2: 1.4440,
            cladding: Cladding::Oxide,
            margin: 0.6e-6,
            phase: std::f64::consts::PI,
            n_modes: 4,
            order: 2,
            h_core: 45e-9,
            h_near: 110e-9,
            h_far: 0.35e-6,
            layers: (4, 4),
        }
    }
}

/// Domain numbers: the layer with teeth uses these, the layer without them the same plus `NO_TEETH`.
pub const OXIDE_BOX: i32 = 1;
pub const CORE: i32 = 2;
pub const CLADDING: i32 = 3;
pub const LEFT_TOOTH: i32 = 4;
pub const RIGHT_TOOTH: i32 = 5;
pub const NO_TEETH: i32 = 10;

/// Stop band of a grating from the two band-edge modes.
#[derive(Clone, Debug)]
pub struct BandEdges {
    /// frequencies (Hz) of the two lowest modes at the Floquet phase `pi`
    pub f_low: f64,
    pub f_high: f64,
    /// `(f_high - f_low) / mean`
    pub gap_fraction: f64,
    /// coupling coefficient (1/m) from the gap and the group index of the unperturbed wire
    pub kappa: f64,
    pub group_index: f64,
    /// width of the stop band in wavelength (m)
    pub stopband: f64,
    /// Bragg wavelength (m) at the middle of the gap
    pub bragg_wavelength: f64,
    /// coupled-mode estimate of `kappa` (1/m) from the effective indices of the two sections
    pub kappa_estimate: f64,
}

impl BraggGrating {
    /// A coarse mesh and a period matched to 1550 nm: a few seconds to a minute instead of minutes, frequencies good to a few tenths of a percent.
    pub fn coarse() -> Self {
        BraggGrating { period: 280e-9, margin: 0.45e-6, h_core: 65e-9, h_near: 170e-9, h_far: 0.45e-6, layers: (3, 3), ..Default::default() }
    }

    pub fn symmetric() -> Self {
        BraggGrating::default()
    }
    pub fn asymmetric() -> Self {
        BraggGrating { kind: GratingKind::Asymmetric, ..Default::default() }
    }

    /// The unperturbed wire of width `width + extra` as a 2D mode problem.
    fn wire(&self, extra: f64) -> SoiWire {
        SoiWire { wavelength: self.wavelength, width: self.width + extra, thickness: self.thickness, n_si: self.n_si, n_sio2: self.n_sio2, cladding: self.cladding, n_modes: 2, ..Default::default() }
    }

    /// Width of the wire averaged over the period.
    pub fn mean_width(&self) -> f64 {
        let teeth = match self.kind {
            GratingKind::Symmetric | GratingKind::HalfPeriodShifted => 2.0,
            GratingKind::Asymmetric => 1.0,
            GratingKind::Uniform => 0.0,
        };
        self.width + teeth * self.tooth * self.duty
    }

    /// The cross-section of the widest part: domains 1 box, 2 core, 3 cladding, 4 and 5 left and right tooth. x lateral, y normal.
    pub fn cross_section(&self) -> Spec2d {
        let (w, t, d) = (self.width, self.thickness, self.tooth);
        let half = 0.5 * w + d + self.margin;
        let (y0, y1) = (-self.margin - 0.1e-6, t + self.margin);
        let rect = |xa: f64, ya: f64, xb: f64, yb: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: xa, y0: ya, x1: xb, y1: yb }, domain: dom, h };
        let zone = 0.5 * w + d + 0.3e-6;
        Spec2d {
            parts: vec![
                rect(-half, y0, half, y1, CLADDING, None),
                rect(-half, y0, half, 0.0, OXIDE_BOX, None),
                rect(-zone, -0.3e-6, zone, 0.0, OXIDE_BOX, Some(self.h_near)),
                rect(-zone, t, zone, t + 0.3e-6, CLADDING, Some(self.h_near)),
                rect(-0.5 * w, 0.0, 0.5 * w, t, CORE, Some(self.h_core)),
                rect(-0.5 * w - d, 0.0, -0.5 * w, t, LEFT_TOOTH, Some(self.h_core)),
                rect(0.5 * w, 0.0, 0.5 * w + d, t, RIGHT_TOOTH, Some(self.h_core)),
            ],
            h: self.h_far,
            min_angle_deg: 25.0,
        }
    }

    /// The extrusion layers: teeth present (domains as they are), then absent (domain numbers plus `NO_TEETH`).
    pub fn layers(&self) -> Vec<Layer> {
        vec![
            Layer { thickness: self.duty * self.period, n: self.layers.0, domain_offset: 0 },
            Layer { thickness: (1.0 - self.duty) * self.period, n: self.layers.1, domain_offset: NO_TEETH },
        ]
    }

    /// The mesh of one period and the boundary numbers of its two ends (`source` z = 0, `dest` z = period).
    pub fn mesh(&self) -> Result<(SimplexMesh, i32, i32), String> {
        if !(self.duty > 0.0 && self.duty < 1.0) || self.period <= 0.0 {
            return Err("the duty cycle must lie between 0 and 1 and the period must be positive".into());
        }
        let m2 = mesh_2d(&self.cross_section())?;
        let nb = m2.facet_labels.iter().copied().max().unwrap_or(0);
        let m3 = extrude(&m2, &self.layers(), 0.0)?;
        Ok((m3, nb + 1, nb + 2))
    }

    /// Relative permittivity of each domain of the cell.
    fn permittivities(&self) -> Vec<(i32, f64)> {
        let (si, ox) = (self.n_si * self.n_si, self.n_sio2 * self.n_sio2);
        let top = match self.cladding {
            Cladding::Oxide => ox,
            Cladding::Air => 1.0,
        };
        // (left tooth with teeth, right tooth with teeth, left without, right without)
        let (l_with, r_with, l_without, r_without) = match self.kind {
            GratingKind::Symmetric => (si, si, ox, ox),
            GratingKind::Asymmetric => (si, ox, ox, ox),
            GratingKind::HalfPeriodShifted => (si, ox, ox, si),
            GratingKind::Uniform => (ox, ox, ox, ox),
        };
        let mut v = vec![];
        for off in [0, NO_TEETH] {
            v.push((OXIDE_BOX + off, ox));
            v.push((CORE + off, si));
            v.push((CLADDING + off, top));
        }
        v.extend([(LEFT_TOOTH, l_with), (RIGHT_TOOTH, r_with), (LEFT_TOOTH + NO_TEETH, l_without), (RIGHT_TOOTH + NO_TEETH, r_without)]);
        v
    }

    /// The eigenmode problem of the cell: Floquet pair on the two ends, modes near `f_search` (default: just below the
    /// frequency at which the unperturbed wire of mean width satisfies the Bragg condition for this period).
    pub fn spec(&self, source: i32, dest: i32, f_search: f64) -> ProblemSpec {
        let mut s = ProblemSpec::new(Physics::Emw);
        s.order = self.order;
        for (d, eps) in self.permittivities() {
            s.materials.set(d, "eps_r", eps);
            s.materials.set(d, "mu_r", 1.0);
        }
        s.periodic = vec![Periodic { source: vec![source], dest: vec![dest], phase: self.phase }];
        s.study = Study::Eigenfrequency { n_modes: self.n_modes, f_search };
        s
    }

    /// The Bragg frequency of the unperturbed wire of mean width: `f = c / (2 period n_eff)` for a phase of `pi`.
    pub fn bragg_frequency_estimate(&self) -> Result<f64, String> {
        let n = fundamental_neff(&self.wire(self.mean_width() - self.width))?;
        Ok(C0 * (self.phase / self.period) / (2.0 * std::f64::consts::PI * n))
    }

    /// The period at which the unperturbed wire of mean width is Bragg-matched at `wavelength`: `lambda / (2 n_eff)`.
    pub fn matched_period(&self) -> Result<f64, String> {
        Ok(self.wavelength / (2.0 * fundamental_neff(&self.wire(self.mean_width() - self.width))?))
    }

    /// Meshes and solves the cell at the configured Floquet phase.
    pub fn solve(&self) -> Result<SolveResult, String> {
        let (mesh, src, dst) = self.mesh()?;
        let f0 = 0.95 * self.bragg_frequency_estimate()?;
        crate::solve::solve(&mesh, &self.spec(src, dst, f0), &mut SolveState::default())
    }

    /// The stop band: solves the cell at the band edge (phase `pi`) and measures the splitting of the two lowest modes.
    pub fn band_edges(&self) -> Result<BandEdges, String> {
        let at_edge = BraggGrating { phase: std::f64::consts::PI, n_modes: self.n_modes.max(2), ..self.clone() };
        let r = at_edge.solve()?;
        at_edge.band_edges_of(&r)
    }

    /// The stop band from a solution of the cell at the band edge (`solve` with the phase `pi`).
    pub fn band_edges_of(&self, r: &SolveResult) -> Result<BandEdges, String> {
        let mut f: Vec<f64> = r.solutions.iter().filter_map(|s| s.value).collect();
        f.sort_by(|a, b| a.total_cmp(b));
        if f.len() < 2 {
            return Err("fewer than two band-edge modes were found".into());
        }
        let (f_low, f_high) = (f[0], f[1]);
        let mean = 0.5 * (f_low + f_high);
        let wire = self.wire(self.mean_width() - self.width);
        let ng = group_index(&wire)?;
        let kappa = std::f64::consts::PI * ng * (f_high - f_low) / C0;
        // coupled-mode estimate from the effective indices of the two sections of the teeth
        let extra = match self.kind {
            GratingKind::Symmetric | GratingKind::HalfPeriodShifted => 2.0 * self.tooth,
            GratingKind::Asymmetric => self.tooth,
            GratingKind::Uniform => 0.0,
        };
        let delta = (fundamental_neff(&self.wire(extra))? - fundamental_neff(&self.wire(0.0))?).abs();
        let kappa_estimate = 2.0 * delta * (std::f64::consts::PI * self.duty).sin() / self.wavelength;
        Ok(BandEdges { f_low, f_high, gap_fraction: (f_high - f_low) / mean, kappa, group_index: ng, stopband: C0 / mean * (f_high - f_low) / mean, bragg_wavelength: C0 / mean, kappa_estimate })
    }
}

/// Effective index of the fundamental quasi-TE mode of a wire (the mode with the largest x share among the modes found).
pub fn fundamental_neff(w: &SoiWire) -> Result<f64, String> {
    let r = w.solve()?;
    let a = r.analysis.as_ref().ok_or("the mode study carries no polarization analysis")?;
    let best = a
        .rows
        .iter()
        .enumerate()
        .filter_map(|(i, row)| row.report.as_ref().map(|m| (i, m.te_fraction)))
        .filter(|&(_, te)| te > 0.5)
        .max_by(|x, y| r.solutions[x.0].neff.unwrap().0.total_cmp(&r.solutions[y.0].neff.unwrap().0))
        .ok_or("no quasi-TE mode was found")?;
    Ok(r.solutions[best.0].neff.unwrap().0)
}

/// Group index `n_g = n_eff - lambda dn_eff/dlambda` of the fundamental quasi-TE mode, by central differences over +-1 % in wavelength.
pub fn group_index(w: &SoiWire) -> Result<f64, String> {
    let rel = 0.01;
    let (a, b) = (SoiWire { wavelength: w.wavelength * (1.0 - rel), ..w.clone() }, SoiWire { wavelength: w.wavelength * (1.0 + rel), ..w.clone() });
    let (na, nb, n0) = (fundamental_neff(&a)?, fundamental_neff(&b)?, fundamental_neff(w)?);
    let dn_dl = (nb - na) / (w.wavelength * 2.0 * rel);
    Ok(n0 - w.wavelength * dn_dl)
}
