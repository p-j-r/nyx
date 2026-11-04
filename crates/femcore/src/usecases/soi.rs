//! Use case: the optical modes of a silicon wire on 300 nm silicon-on-insulator at 1550 nm.
//!
//! A silicon core of thickness `thickness` (300 nm by default) and `width` on a buried oxide (BOX), under an oxide or air
//! cladding. The cross-section is a geom model (exact lines, named regions) meshed by `geomdrive`, the modes come from the full-vector solver of `emw` (Nedelec
//! degree 2 + P2, perfectly conducting walls far from the core) and every solve carries the Stokes/Poincare analysis of
//! `analysis`: the quasi-TE and quasi-TM modes sit at the two ends of the equator of the sphere, hybrid modes (sloped
//! side walls, rib slab, near the TE-TM crossing) move off the poles of the x/y axis and lose degree of polarization.
//!
//! Regions (by name; the integers are the solver's business): `BOX` buried oxide, `Si` core, `Si_slab` rib slab, `clad` upper cladding, and
//! the `_near` zones of finer elements. Wafer axes: x lateral, the cross-section is the x-z plane, z the surface normal, y the propagation.
//!
//! Refractive indices at 1550 nm: Si 3.4757 (Li, 1980), SiO2 1.4440 (Malitson, 1965).
//!
//! The thin-slab limit gives an analytic reference: for a very wide wire the quasi-TE mode tends to the TE0 mode of the
//! symmetric slab and the quasi-TM mode to its TM0 mode (`slab_neff`).

use crate::meshing::geomdrive::MeshedGeom;
use crate::problem::*;
use geom::model::Model;
use geom::select::Selector;

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
    /// The cross-section as a geom model in the wafer frame: x lateral, z the surface normal, the silicon starts at z = 0. Regions are
    /// `BOX`, `Si` (the core), `Si_slab`, `clad`, plus `BOX_near` and `clad_near`, the same oxides with finer elements on both sides of
    /// the oxide interface.
    pub fn model(&self) -> Model {
        let mut m = Model::new_2d();
        m.param("w", self.width)
            .param("t", self.thickness)
            .param("slab", self.slab)
            .param("margin", self.margin)
            .param("boxt", self.box_thickness)
            .param("zone", 0.5 * self.width + 0.4e-6);
        let (right, left) = (self.sidewall_deg, self.sidewall_left_deg.unwrap_or(self.sidewall_deg));
        m.param("ar", right).param("al", left);
        let clad = m.rect("clad", ("-w/2-margin", "-boxt"), ("w/2+margin", "t+margin"));
        let bx = m.rect("box", ("-w/2-margin", "-boxt"), ("w/2+margin", 0.0));
        let box_near = m.rect("box_near", ("-zone", -0.4e-6), ("zone", 0.0));
        let clad_near = m.rect("clad_near", ("-zone", "t"), ("zone", "t+0.4e-6"));
        // the top corners move inwards by (height of the core above the slab) / tan(angle)
        let vertical = right >= 90.0 && left >= 90.0;
        let core = if vertical {
            m.rect("core", ("-w/2", "slab"), ("w/2", "t"))
        } else {
            let inset = |a: &str| if (if a == "ar" { right } else { left }) < 90.0 { format!("(t-slab)/tan({a}*pi/180)") } else { "0".to_string() };
            let pts: [(String, String); 4] = [
                ("-w/2".into(), "slab".into()),
                ("w/2".into(), "slab".into()),
                (format!("w/2-{}", inset("ar")), "t".into()),
                (format!("-w/2+{}", inset("al")), "t".into()),
            ];
            m.polygon_named("core", &pts, &["bottom", "right", "top", "left"])
        };
        let mut regions = vec![("clad", clad), ("BOX", bx), ("BOX_near", box_near), ("clad_near", clad_near)];
        let mut near = "region(BOX_near) | region(clad_near)".to_string();
        if self.slab > 0.0 {
            let slab = m.rect("slab_rect", ("-w/2-margin", 0.0), ("w/2+margin", "slab"));
            regions.push(("Si_slab", slab));
            near.push_str(" | region(Si_slab)");
        }
        regions.push(("Si", core));
        m.overlay(&regions);
        m.size("size_near", Selector::parse(&near).expect("a valid selector"), Some(self.h_near.into()), None, None);
        m.size("size_core", Selector::parse("region(Si)").expect("a valid selector"), Some(self.h_core.into()), None, None);
        m
    }

    /// The cross-section meshed from the exact model: the mesh, the names of its domains and boundaries, and the evaluation.
    pub fn mesh(&self) -> Result<MeshedGeom, String> {
        let ev = crate::meshing::geomdrive::evaluate(&self.model().list)?;
        crate::meshing::mesh_evaluated(ev, &crate::meshing::MeshPlan::simple(self.h_far))
    }

    /// Optical problem: wafer axes, permittivity per region (by name), `n_modes` modes near the index of the silicon.
    pub fn spec(&self, g: &MeshedGeom) -> ProblemSpec {
        let mut s = ProblemSpec::new(Physics::Emw);
        s.order = 2;
        let n_top = match self.cladding {
            Cladding::Oxide => self.n_sio2,
            Cladding::Air => 1.0,
        };
        for (regions, n) in [("BOX, BOX_near", self.n_sio2), ("Si, Si_slab", self.n_si), ("clad, clad_near", n_top)] {
            g.set_material(&mut s, regions, "eps_r", n * n);
            g.set_material(&mut s, regions, "mu_r", 1.0);
        }
        s.study = Study::Modes { freq: C0 / self.wavelength, n_target: self.n_target, evanescent: false, n_modes: self.n_modes };
        s
    }

    /// Meshes and solves; the result carries the polarization analysis (`SolveResult::analysis`).
    pub fn solve(&self) -> Result<SolveResult, String> {
        let g = self.mesh()?;
        crate::solve::solve(&g.mesh, &self.spec(&g), &mut SolveState::default())
    }

    /// The same silicon thickness as an infinite slab: reference for very wide wires.
    pub fn slab_reference(&self, pol: SlabPol) -> Option<f64> {
        slab_neff(self.n_si, self.n_sio2, self.thickness, self.wavelength, pol, 0)
    }
}
