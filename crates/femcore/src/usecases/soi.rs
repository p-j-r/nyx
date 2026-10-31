//! Use case: the optical modes of a silicon wire on 300 nm silicon-on-insulator at 1550 nm.
//!
//! A silicon core of thickness `thickness` (300 nm by default) and `width` on a buried oxide (BOX), under an oxide or air
//! cladding. The cross-section is meshed with the shape mesher, the modes come from the full-vector solver of `emw` (Nedelec
//! degree 2 + P2, perfectly conducting walls far from the core) and every solve carries the Stokes/Poincare analysis of
//! `analysis`: the quasi-TE and quasi-TM modes sit at the two ends of the equator of the sphere, hybrid modes (sloped
//! side walls, rib slab, near the TE-TM crossing) move off the poles of the x/y axis and lose degree of polarization.
//!
//! Domains: 1 buried oxide, 2 silicon (core and rib slab), 3 upper cladding. Wafer axes: x lateral, the 2D mesh y is the surface normal.
//!
//! Refractive indices at 1550 nm: Si 3.4757 (Li, 1980), SiO2 1.4440 (Malitson, 1965).
//!
//! The thin-slab limit gives an analytic reference: for a very wide wire the quasi-TE mode tends to the TE0 mode of the
//! symmetric slab and the quasi-TM mode to its TM0 mode (`slab_neff`).

use crate::mesh::SimplexMesh;
use crate::mesher::{mesh_2d, Part, Shape, Spec2d};
use crate::problem::*;

const C0: f64 = 299_792_458.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cladding {
    /// silicon dioxide above the core
    Oxide,
    /// air above the core
    Air,
}

/// Polarization of a slab mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlabPol {
    /// field parallel to the film (quasi-TE of a wide wire)
    Te,
    /// field with a component normal to the film (quasi-TM)
    Tm,
}

/// Effective index of mode `m` of a symmetric slab of index `n_core` and thickness `t` in a cladding `n_clad`:
/// `k_x t = m pi + 2 atan(p gamma / k_x)` with `p = 1` (TE) or `(n_core / n_clad)^2` (TM). `None` if the slab does not guide the mode.
pub fn slab_neff(n_core: f64, n_clad: f64, t: f64, wavelength: f64, pol: SlabPol, m: usize) -> Option<f64> {
    let k0 = 2.0 * std::f64::consts::PI / wavelength;
    let p = match pol {
        SlabPol::Te => 1.0,
        SlabPol::Tm => (n_core / n_clad).powi(2),
    };
    let f = |n: f64| {
        let kx = k0 * (n_core * n_core - n * n).max(0.0).sqrt();
        let gamma = k0 * (n * n - n_clad * n_clad).max(0.0).sqrt();
        kx * t - m as f64 * std::f64::consts::PI - 2.0 * (p * gamma / kx.max(1e-300)).atan()
    };
    let (mut lo, mut hi) = (n_clad * (1.0 + 1e-12), n_core * (1.0 - 1e-12));
    if f(lo) <= 0.0 || f(hi) >= 0.0 {
        return None;
    }
    // f falls with n: positive at the cladding index (mode present), negative near the core index
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// The geometry, materials and study of the use case.
#[derive(Clone, Debug)]
pub struct SoiWire {
    pub wavelength: f64,
    /// width of the core at its foot
    pub width: f64,
    /// silicon thickness (the SOI layer)
    pub thickness: f64,
    /// thickness of a silicon slab under the core, spanning the whole width of the domain (0: a strip wire)
    pub slab: f64,
    /// side wall angle from the horizontal in degrees (90: vertical walls); the core is narrower at the top. This is the angle of the right wall,
    /// and of the left one too unless `sidewall_left_deg` is set.
    pub sidewall_deg: f64,
    /// a different angle for the left wall: the wire loses its mirror symmetry and the modes become hybrid (their Stokes vectors leave the x/y axis)
    pub sidewall_left_deg: Option<f64>,
    pub n_si: f64,
    pub n_sio2: f64,
    pub cladding: Cladding,
    /// oxide between the core and the walls of the domain, laterally and above
    pub margin: f64,
    /// thickness of the buried oxide
    pub box_thickness: f64,
    pub n_modes: usize,
    /// the `n_modes` modes nearest this effective index are returned (the shift of the eigensolver). The default suits wires of 400 to 700 nm; for
    /// a much wider wire put it just above the index of the mode you want (the modes of higher lateral order crowd in between), and keep `n_modes`
    /// small: modes close to the cladding index take the iteration a long time.
    pub n_target: f64,
    /// element size in the silicon, in the oxide near the core and far away
    pub h_core: f64,
    pub h_near: f64,
    pub h_far: f64,
}

impl Default for SoiWire {
    fn default() -> Self {
        SoiWire {
            wavelength: 1550e-9,
            width: 500e-9,
            thickness: 300e-9,
            slab: 0.0,
            sidewall_deg: 90.0,
            sidewall_left_deg: None,
            n_si: 3.4757,
            n_sio2: 1.4440,
            cladding: Cladding::Oxide,
            margin: 1.0e-6,
            box_thickness: 2.0e-6,
            n_modes: 2,
            n_target: 2.75,
            h_core: 40e-9,
            h_near: 90e-9,
            h_far: 0.3e-6,
        }
    }
}

impl SoiWire {
    /// The shape specification of the cross-section: x lateral, y the surface normal, the silicon starts at y = 0.
    pub fn shapes(&self) -> Spec2d {
        let (w, t) = (self.width, self.thickness);
        let half = 0.5 * w + self.margin;
        let (y0, y1) = (-self.box_thickness, t + self.margin);
        let rect = |xa: f64, ya: f64, xb: f64, yb: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: xa, y0: ya, x1: xb, y1: yb }, domain: d, h };
        let (right, left) = (self.sidewall_deg, self.sidewall_left_deg.unwrap_or(self.sidewall_deg));
        // the top corners move inwards by (height of the core above the slab) / tan(angle)
        let inset = |deg: f64| if deg < 90.0 { (t - self.slab) / deg.to_radians().tan() } else { 0.0 };
        let core = if inset(right) == 0.0 && inset(left) == 0.0 {
            Shape::Rect { x0: -0.5 * w, y0: self.slab, x1: 0.5 * w, y1: t }
        } else {
            Shape::Polygon(vec![[-0.5 * w, self.slab], [0.5 * w, self.slab], [0.5 * w - inset(right), t], [-0.5 * w + inset(left), t]])
        };
        let zone = 0.5 * w + 0.4e-6;
        let mut parts = vec![
            rect(-half, y0, half, y1, 3, None),
            rect(-half, y0, half, 0.0, 1, None),
            // finer cells where the field lives, on both sides of the oxide interface
            rect(-zone, -0.4e-6, zone, 0.0, 1, Some(self.h_near)),
            rect(-zone, t, zone, t + 0.4e-6, 3, Some(self.h_near)),
        ];
        if self.slab > 0.0 {
            parts.push(rect(-half, 0.0, half, self.slab, 2, Some(self.h_near)));
        }
        parts.push(Part { shape: core, domain: 2, h: Some(self.h_core) });
        Spec2d { parts, h: self.h_far, min_angle_deg: 25.0 }
    }

    pub fn mesh(&self) -> Result<SimplexMesh, String> {
        mesh_2d(&self.shapes())
    }

    /// Optical problem: wafer axes, permittivity per domain, `n_modes` modes near the index of the silicon.
    pub fn spec(&self) -> ProblemSpec {
        let mut s = ProblemSpec::new(Physics::Emw);
        s.order = 2;
        s.materials.frame = Frame::Wafer;
        let n_top = match self.cladding {
            Cladding::Oxide => self.n_sio2,
            Cladding::Air => 1.0,
        };
        for (d, n) in [(1, self.n_sio2), (2, self.n_si), (3, n_top)] {
            s.materials.set(d, "eps_r", n * n);
            s.materials.set(d, "mu_r", 1.0);
        }
        s.study = Study::Modes { freq: C0 / self.wavelength, n_target: self.n_target, evanescent: false, n_modes: self.n_modes };
        s
    }

    /// Meshes and solves; the result carries the polarization analysis (`SolveResult::analysis`).
    pub fn solve(&self) -> Result<SolveResult, String> {
        let mesh = self.mesh()?;
        crate::solve::solve(&mesh, &self.spec(), &mut SolveState::default())
    }

    /// The same silicon thickness as an infinite slab: reference for very wide wires.
    pub fn slab_reference(&self, pol: SlabPol) -> Option<f64> {
        slab_neff(self.n_si, self.n_sio2, self.thickness, self.wavelength, pol, 0)
    }
}
