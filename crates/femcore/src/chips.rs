//! The silicon-on-sapphire rib waveguides of the four ring chips, as Brillouin problems in the
//! wafer frame: x lateral, y the direction of propagation, z the surface normal; the cross-section is the x-z plane.
//!
//! Regions (by name): `sapphire`, `slab`, `rib`, `air`, and the pieces of the absorbing layer around the physical region
//! (`sapphire.pml_xmin`, `sapphire.pml_zmin`, `sapphire.pml_xmin_zmin`, `slab.pml_xmax`, `air.pml_zmax`, ...). The physical region is 2 um (3 um for the 500TE chip) wide on each
//! side of the rib, 2 um of sapphire below, 1.2 um of air above, with 1 um absorbing layers around it.
//!
//! The crystals are oriented as follows: the silicon axes are `[p, q, z]` and the sapphire axes
//! `[a, c x a, c]` with `p = (sin psi, -cos psi, 0)` (the projection of the sapphire c axis on the wafer; at `psi = 0` it
//! lies along the waveguide), `q = (cos psi, sin psi, 0)` the in-plane a axis and `c = sin(57.6 deg) p + cos(57.6 deg) z`
//! (r-plane sapphire). The wafer frame is (lateral x, propagation y, normal z), right-handed; a wave along `-y` is as good as one along `+y`: the
//! elastic and optical modes have the same frequencies and gains.
//! The sapphire permittivity is the uniaxial tensor `diag(n_o^2, n_o^2, n_e^2)` of the same crystal.

use crate::geomdrive::{mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use crate::material::*;
use crate::mesh::SimplexMesh;
use crate::problem::*;
use geom::model::Model;
use geom::select::Selector;

#[derive(Clone, Copy, Debug)]
pub struct ChipSpec {
    pub name: &'static str,
    /// rib height (total thickness of the silicon) in um
    pub htot: f64,
    /// etch depth in um (the slab is `htot - etch` thick)
    pub etch: f64,
    /// rib width in um
    pub width: f64,
    /// half width of the physical region in um
    pub xin: f64,
    /// effective index of the design mode (the search target for the optical modes)
    pub n_eff: f64,
    pub polarisation: &'static str,
}

pub const CHIPS: [ChipSpec; 4] = [
    ChipSpec { name: "300TE", htot: 0.3, etch: 0.17, width: 0.5, xin: 2.0, n_eff: 2.784, polarisation: "TE" },
    ChipSpec { name: "500TE", htot: 0.5, etch: 0.15, width: 1.0, xin: 3.0, n_eff: 3.2223, polarisation: "TE" },
    ChipSpec { name: "300TM", htot: 0.3, etch: 0.235, width: 0.5, xin: 2.0, n_eff: 2.3677, polarisation: "TM" },
    ChipSpec { name: "500TM", htot: 0.5, etch: 0.325, width: 0.5, xin: 2.0, n_eff: 2.9401, polarisation: "TM" },
];

pub fn chip_by_name(name: &str) -> Option<ChipSpec> {
    CHIPS.iter().copied().find(|c| c.name.eq_ignore_ascii_case(name))
}

/// Optical constants at 1550 nm.
pub const N_SILICON: f64 = 3.476;
pub const N_SAPPHIRE_ORDINARY: f64 = 1.7462;
pub const N_SAPPHIRE_EXTRAORDINARY: f64 = 1.7383;

/// The cross-section as a geom model in the wafer frame (x lateral, z the surface normal, lengths in m): sapphire, the silicon slab and
/// rib, and air, with a 1 um absorbing layer around the physical region (`air`, `sapphire`, `slab`, `rib` and their pieces
/// `<region>.pml_xmin`, `.pml_zmin`, `.pml_xmin_zmin`, ...). `scale` multiplies the element sizes.
pub fn chip_geometry(c: &ChipSpec, scale: f64) -> Model {
    let um = 1e-6;
    let (w, htot) = (c.width * um, c.htot * um);
    let hslab = (c.htot - c.etch) * um;
    let (xin, hin, lp, ha) = (c.xin * um, 2.0 * um, 1.0 * um, 1.2 * um);
    let x = xin + lp;
    let (h_si, h_sap, h_pml) = (0.028 * scale * um, 0.05 * scale * um, 0.07 * scale * um);
    let ztop = htot + ha;
    let mut m = Model::new_2d();
    let air = m.rect("air_box", (-x, -hin - lp), (x, ztop + lp));
    let sapphire = m.rect("sapphire_box", (-x, -hin - lp), (x, 0.0));
    let rib = m.rect("rib_box", (-w / 2.0, hslab), (w / 2.0, htot));
    let window = m.rect("window", (-xin, -hin), (xin, ztop));
    // a fully etched strip (`etch == htot`) has no slab
    let has_slab = hslab > 1e-12;
    if has_slab {
        let slab = m.rect("slab_box", (-x, 0.0), (x, hslab));
        m.overlay(&[("air", air), ("sapphire", sapphire), ("slab", slab), ("rib", rib)]);
    } else {
        m.overlay(&[("air", air), ("sapphire", sapphire), ("rib", rib)]);
    }
    // the absorbing layer: pieces of every region outside the window, as many elements across as the old fixed size gave
    m.pml("pml", window, &[], lp, lp / h_pml, 2.0);
    let si = if has_slab { "region(slab) | region(rib)" } else { "region(rib)" };
    m.size("size_si", Selector::parse(si).expect("a valid selector"), Some(h_si.into()), None, None);
    m.size("size_sapphire", Selector::parse("region(sapphire)").expect("a valid selector"), Some(h_sap.into()), None, None);
    m
}

/// The mesh of a chip: global element size 0.12 um.
pub fn chip_mesh(c: &ChipSpec, scale: f64) -> Result<MeshedGeom, String> {
    let ev = crate::geomdrive::evaluate(&chip_geometry(c, scale).list)?;
    mesh_geom_2d(ev, &GeomMeshOpts::new(0.12e-6))
}

fn cross(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
}

fn columns(c1: [f64; 3], c2: [f64; 3], c3: [f64; 3]) -> M3 {
    [[c1[0], c2[0], c3[0]], [c1[1], c2[1], c3[1]], [c1[2], c2[2], c3[2]]]
}

/// Orientation matrices (columns = crystal axes in the wafer frame) of silicon and sapphire at the azimuth `psi`.
pub fn crystal_rotations(psi_deg: f64) -> (M3, M3) {
    let (ps, b) = (psi_deg.to_radians(), 57.6f64.to_radians());
    let p = [ps.sin(), -ps.cos(), 0.0];
    let q = [ps.cos(), ps.sin(), 0.0];
    let z = [0.0, 0.0, 1.0];
    let c = [b.sin() * p[0], b.sin() * p[1], b.cos()];
    (columns(p, q, z), columns(q, cross(c, q), c))
}

/// The Brillouin problem of a chip on its mesh (backward scattering of the design mode with itself): `psi_deg` is the crystal azimuth,
/// `n_modes` and `f_search` (Hz) set the elastic search. Materials, absorbing layers, inactive air and the clamped outer walls are
/// written by region and boundary names.
pub fn chip_problem(c: &ChipSpec, psi_deg: f64, n_modes: usize, f_search: f64, g: &MeshedGeom) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Brillouin);
    s.order = 2;
    let (sapphire, silicon, air) = (g.domains("sapphire*"), g.domains("slab*, rib*"), g.domains("air*"));
    s.inactive = air.clone();
    let (r_si, r_sap) = crystal_rotations(psi_deg);
    let (eul_si, eul_sap) = (euler_from_rotation(&r_si), euler_from_rotation(&r_sap));
    for d in g.labels.domains.iter().map(|d| d.id) {
        s.materials.set(d, "E", 100e9);
        s.materials.set(d, "nu", 0.25);
        s.materials.set(d, "rho", 1.0);
    }
    for &d in &sapphire {
        s.materials.aniso.insert(
            d,
            Anisotropy {
                stiffness: Some(Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 }),
                permittivity: Some(Tensor::Diagonal { x: N_SAPPHIRE_ORDINARY.powi(2), y: N_SAPPHIRE_ORDINARY.powi(2), z: N_SAPPHIRE_EXTRAORDINARY.powi(2) }),
                euler_deg: eul_sap,
                ..Default::default()
            },
        );
        s.materials.set(d, "rho", 3980.0);
    }
    for &d in &silicon {
        s.materials.aniso.insert(d, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }), euler_deg: eul_si, ..Default::default() });
        s.materials.set(d, "rho", 2329.0);
    }
    // the absorbing layers of the model, except in the inactive air
    s.materials.pml.extend(g.pml_materials(3.0).into_iter().filter(|(d, _)| !air.contains(d)));
    // the clamped walls: the outer boundary of the sapphire and the slab (not the air above)
    let hslab = (c.htot - c.etch) * 1e-6;
    let outer = g.select(&format!("outer & in_box([-1, 0, -1], [1, 0, {}], within)", hslab + 1e-12)).expect("a valid selector");
    s.bcs = vec![Bc { on: None, labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s.study = Study::Eigenfrequency { n_modes, f_search };
    let b = s.brillouin.as_mut().unwrap();
    b.wavelength = 1.55e-6;
    b.n_target = c.n_eff;
    b.n_optical_modes = 1;
    b.backward = true;
    b.pump_modes = vec![1];
    b.stokes_modes = vec![1];
    for &d in &silicon {
        b.refractive_index.insert(d, N_SILICON);
    }
    for &d in &sapphire {
        b.refractive_index.insert(d, N_SAPPHIRE_ORDINARY);
    }
    for &d in &air {
        b.refractive_index.insert(d, 1.0);
    }
    for d in g.domains("slab, rib") {
        b.photoelastic.insert(d, PhotoSpec { kind: PhotoKind::Cubic { p11: -0.094, p12: 0.017, p44: -0.051 }, euler_deg: [0.0; 3], follow_stiffness: true });
    }
    for d in g.domains("sapphire") {
        b.photoelastic.insert(d, PhotoSpec { kind: PhotoKind::Trigonal { p11: -0.23, p12: 0.03, p13: 0.02, p14: 0.0, p33: -0.2, p44: -0.1 }, euler_deg: [0.0; 3], follow_stiffness: true });
    }
    b.moving = g.domains("slab, rib");
    s
}

/// Mesh and problem of a chip (use `chip_mesh` and `chip_problem` to keep the names of regions and boundaries).
pub fn chip_model(c: &ChipSpec, psi_deg: f64, scale: f64, n_modes: usize, f_search: f64) -> Result<(SimplexMesh, ProblemSpec), String> {
    let g = chip_mesh(c, scale)?;
    let spec = chip_problem(c, psi_deg, n_modes, f_search, &g);
    Ok((g.mesh, spec))
}
