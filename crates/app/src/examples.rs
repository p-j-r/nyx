//! Ready-made projects (geometry, materials, boundary conditions, studies) to start from.
//!
//! Every geometry is a `geom` model. Regions and boundaries are looked up by name (`MeshedGeom::domains`, `select`) on the mesh the
//! model gives, so no example types a region or boundary number. A 3D example is a 2D x-z section swept along y, the propagation axis.

use crate::colormap::Colormap;
use crate::state::*;
use femcore::geomdrive::MeshedGeom;
use femcore::mesher::Layer;
use femcore::problem::*;
use geom::model::Model;
use geom::select::Selector;

fn base(geometry: GeometrySpec, active: Physics) -> Project {
    let mut specs: Vec<ProblemSpec> = Physics::ALL.iter().map(|p| ProblemSpec::new(*p)).collect();
    for s in &mut specs {
        for (k, v) in [("eps_r", 1.0), ("mu_r", 1.0), ("k", 1.0), ("rho", 1000.0), ("Cp", 500.0), ("E", 200e9), ("nu", 0.3), ("alpha", 12e-6)] {
            s.materials.default.insert(k.into(), v);
        }
        s.materials.default.insert("rho".into(), if s.physics == Physics::Solid { 7850.0 } else { 1000.0 });
    }
    // everything is written in the wafer frame: x lateral, y propagation, z surface normal; a 2D mesh is the x-z cross-section
    for s in &mut specs {
        s.materials.frame = Frame::Wafer;
    }
    Project { version: 1, geometry, source: Source::Generated, active, specs, colormap: Colormap::Viridis }
}

/// The model meshed once, for the names of its regions and boundaries.
fn named(m: Model, h: f64) -> (GeometrySpec, MeshedGeom) {
    let g = GeometrySpec::model(m.list, h);
    let mg = g.build_named().expect("the example geometry meshes");
    (g, mg)
}

/// A 2D model swept along y through one layer of `len` and `n` element layers: the geometry, the names, and the boundary numbers of the
/// end at y = 0 and the far end (the side walls keep their 2D numbers).
fn swept(m: Model, h: f64, len: f64, n: usize) -> (GeometrySpec, MeshedGeom, i32, i32) {
    let (g, mg) = named(m, h);
    let nb = mg.mesh.facet_labels.iter().copied().max().unwrap_or(0);
    let g = GeometrySpec::model_extruded(g.list, h, vec![Layer { thickness: len, n, domain_offset: 0 }]);
    (g, mg, nb + 1, nb + 2)
}

fn sel(g: &MeshedGeom, rule: &str) -> Vec<i32> {
    g.select(rule).unwrap_or_else(|e| panic!("selector '{rule}': {e}"))
}

fn bc(labels: &[i32], kind: BcKind) -> Bc {
    Bc { labels: labels.to_vec(), kind }
}

fn set(p: &mut Project, phys: Physics, props: &[(&str, f64)]) {
    for (k, v) in props {
        p.spec_mut(phys).materials.default.insert((*k).into(), *v);
    }
}

fn fixed_all() -> BcKind {
    BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] }
}

/// A model of one rectangle `plate` from the origin to `(lx, lz)`.
fn plate(lx: f64, lz: f64) -> Model {
    let mut m = Model::new_2d();
    let r = m.rect("plate", (0.0, 0.0), (lx, lz));
    m.overlay(&[("plate", r)]);
    m
}

fn selector(s: &str) -> Selector {
    Selector::parse(s).expect("a valid selector")
}

pub const NAMES: [&str; 19] = [
    "Parallel-plate capacitor with two dielectrics (electrostatics, 2D)",
    "Cantilever beam (solid mechanics, static + eigenfrequencies, 2D)",
    "Heated plate with convective cooling (heat transfer, 2D, transient)",
    "Thermal expansion of a bar (heat + solid, 2D)",
    "Block cantilever (solid mechanics, 3D, a square section swept along y)",
    "WR-90 rectangular waveguide modes (EM waves, 2D cross-section)",
    "Dielectric-loaded waveguide (EM waves, 2D cross-section)",
    "Sapphire cantilever, anisotropic crystal (solid mechanics, 3D eigenfrequencies)",
    "Dielectric slab guide with PML ends (EM waves, complex modes)",
    "Plate with a hole (solid mechanics, stress concentration)",
    "Beam on a PML substrate: anchor loss (solid eigenfrequencies, complex frequency and Q)",
    "Plate with a drilled hole in 3D, the section swept along y (solid mechanics, stress concentration)",
    "Silicon strip on sapphire at 1550 nm: optical modes (EM waves) and acoustic modes with out-of-plane wave number (solid mechanics)",
    "Rectangular cavity: 3D electromagnetic eigenmodes (40 x 30 x 20 mm, perfect conductor walls)",
    "Brillouin gain of a silicon rib on sapphire (300TE ring chip, backward, 1550 nm): optical mode, phonons, photoelastic and moving-boundary gain",
    "Silicon wire on 300 nm SOI at 1550 nm: TE0 and TM0 modes with Stokes parameters and Poincare sphere (EM waves, 2D cross-section)",
    "300 nm SOI wire with a sloped left wall: hybrid modes tilted on the Poincare sphere (EM waves, 2D cross-section)",
    "Symmetric Bragg grating in a 300 nm SOI wire: band-edge modes of a 3D Floquet cell (teeth on both sides)",
    "Asymmetric Bragg grating in a 300 nm SOI wire: band-edge modes of a 3D Floquet cell (teeth on one side)",
];

/// An empty project to build a simulation from: the chosen physics active, default materials, no boundary conditions and a unit
/// plate (one region `plate`, x-z) or, in 3D, the unit cube made by sweeping it along y. Wafer axes.
pub fn blank(physics: Physics, three_d: bool) -> Project {
    let three_d = three_d && physics != Physics::Brillouin;
    let g = if three_d { GeometrySpec::plate(1.0, 1.0, 0.125, Some(1.0)) } else { GeometrySpec::plate(1.0, 1.0, 0.05, None) };
    let mut p = base(g, physics);
    if three_d && physics == Physics::Emw {
        p.spec_mut(physics).study = Study::Eigenfrequency { n_modes: 6, f_search: 10e9 };
    }
    p
}

pub fn by_index(i: usize) -> Project {
    match i {
        0 => capacitor(),
        1 => cantilever(),
        2 => heated_plate(),
        3 => thermal_bar(),
        4 => block_cantilever(),
        5 => waveguide(),
        6 => loaded_waveguide(),
        7 => sapphire_cantilever(),
        8 => pml_slab(),
        9 => plate_with_hole(),
        10 => anchor_loss(),
        11 => drilled_plate(),
        12 => sos_strip(),
        13 => box_cavity(),
        14 => brillouin_rib(),
        15 => soi_wire(),
        16 => soi_hybrid_wire(),
        17 => bragg_symmetric(),
        _ => bragg_asymmetric(),
    }
}

pub fn capacitor() -> Project {
    let mut m = Model::new_2d();
    let a = m.rect("left_box", (0.0, 0.0), (0.5, 0.5));
    let b = m.rect("right_box", (0.5, 0.0), (1.0, 0.5));
    m.overlay(&[("left", a), ("right", b)]);
    let (g, mg) = named(m, 0.025);
    let mut p = base(g, Physics::Electrostatics);
    let es = p.spec_mut(Physics::Electrostatics);
    for (name, eps) in [("left", 2.0), ("right", 8.0)] {
        mg.set_material(es, name, "eps_r", eps);
    }
    es.bcs = vec![bc(&sel(&mg, "outer & plane(x = 0)"), BcKind::Potential(1.0)), bc(&sel(&mg, "outer & plane(x = 1)"), BcKind::Ground)];
    p
}

pub fn cantilever() -> Project {
    let (g, mg) = named(plate(1.0, 0.1), 0.02);
    let mut p = base(g, Physics::Solid);
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&sel(&mg, "outer & plane(x = 0)"), fixed_all()), bc(&sel(&mg, "outer & plane(x = 1)"), BcKind::Load([0.0, 0.0, -1e5]))];
    s.plane = Plane::Stress;
    p
}

pub fn heated_plate() -> Project {
    let (g, mg) = named(plate(1.0, 0.2), 0.025);
    let mut p = base(g, Physics::Heat);
    set(&mut p, Physics::Heat, &[("k", 50.0), ("rho", 7800.0), ("Cp", 470.0)]);
    let h = p.spec_mut(Physics::Heat);
    h.bcs = vec![bc(&sel(&mg, "outer & plane(x = 0)"), BcKind::Temperature(500.0)), bc(&sel(&mg, "outer - plane(x = 0)"), BcKind::Convection { h: 25.0, text: 293.15 })];
    h.study = Study::Transient { t_end: 20000.0, n_out: 20, t0: "293.15".into() };
    p
}

pub fn thermal_bar() -> Project {
    let (g, mg) = named(plate(1.0, 0.1), 0.02);
    let mut p = base(g, Physics::Solid);
    set(&mut p, Physics::Heat, &[("k", 160.0)]);
    set(&mut p, Physics::Solid, &[("E", 70e9), ("nu", 0.33), ("alpha", 23e-6), ("rho", 2700.0)]);
    let (left, right) = (sel(&mg, "outer & plane(x = 0)"), sel(&mg, "outer & plane(x = 1)"));
    p.spec_mut(Physics::Heat).bcs = vec![bc(&left, BcKind::Temperature(293.15)), bc(&right, BcKind::Temperature(393.15))];
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&left, BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] })];
    s.thermal = Some(ThermalExpansion { source: ThermalSource::FromHeat, t_ref: 293.15 });
    s.plane = Plane::Stress;
    p
}

/// A 0.1 m square bar 1 m long: its section swept along y, clamped at y = 0 and loaded across at the far end.
pub fn block_cantilever() -> Project {
    let (g, _mg, start, end) = swept(plate(0.1, 0.1), 0.035, 1.0, 16);
    let mut p = base(g, Physics::Solid);
    p.spec_mut(Physics::Solid).bcs = vec![bc(&[start], fixed_all()), bc(&[end], BcKind::Load([0.0, 0.0, -1e5]))];
    p
}

pub fn waveguide() -> Project {
    let (g, _) = named(plate(22.86e-3, 10.16e-3), 1.0e-3);
    let mut p = base(g, Physics::Emw);
    p.spec_mut(Physics::Emw).study = Study::Modes { freq: 20e9, n_target: 0.9, evanescent: false, n_modes: 4 };
    p
}

pub fn loaded_waveguide() -> Project {
    let (lx, lz) = (22.86e-3, 10.16e-3);
    let mut m = Model::new_2d();
    let air = m.rect("air_box", (0.0, 0.0), (lx, lz));
    let diel = m.rect("dielectric_box", (0.0, 0.0), (11.43e-3, lz));
    m.overlay(&[("air", air), ("dielectric", diel)]);
    let (g, mg) = named(m, 0.7e-3);
    let mut p = base(g, Physics::Emw);
    let e = p.spec_mut(Physics::Emw);
    mg.set_material(e, "air", "eps_r", 1.0);
    mg.set_material(e, "dielectric", "eps_r", 4.0);
    e.study = Study::Modes { freq: 8e9, n_target: 1.7, evanescent: false, n_modes: 3 };
    p
}

/// A sapphire beam 1 mm long, 0.2 mm wide and 0.1 mm thick, clamped at y = 0: its eigenfrequencies with the full anisotropic stiffness.
pub fn sapphire_cantilever() -> Project {
    // 1 mm long, 0.2 mm wide (x), 0.1 mm thick (z), swept along y
    let (g, mg, start, _end) = swept(plate(0.2e-3, 0.1e-3), 0.05e-3, 1.0e-3, 20);
    let mut p = base(g, Physics::Solid);
    let (st, eps, k, alpha, rho) = femcore::material::sapphire();
    let an = femcore::material::Anisotropy { stiffness: Some(st), conductivity: Some(k), permittivity: Some(eps), expansion: Some(alpha), euler_deg: [0.0, 90.0, 0.0] };
    let plate_domains = mg.domains("plate");
    for s in &mut p.specs {
        for &d in &plate_domains {
            s.materials.aniso.insert(d, an.clone());
        }
        s.materials.default.insert("rho".into(), rho);
    }
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&[start], fixed_all())];
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    p
}

/// A 40 mm x 10 mm window of air with a dielectric slab in the middle and an absorbing layer at each end.
pub fn pml_slab() -> Project {
    let (lx, lz) = (40e-3, 10e-3);
    let mut m = Model::new_2d();
    let air = m.rect("air_box", (0.0, 0.0), (lx, lz));
    let slab = m.rect("slab_box", (15e-3, 0.0), (25e-3, lz));
    let left = m.rect("pml_left_box", (0.0, 0.0), (8e-3, lz));
    let right = m.rect("pml_right_box", (32e-3, 0.0), (lx, lz));
    m.overlay(&[("air", air), ("slab", slab), ("pml_left", left), ("pml_right", right)]);
    let (g, mg) = named(m, 0.8e-3);
    let mut p = base(g, Physics::Emw);
    for s in &mut p.specs {
        for (name, eps) in [("air", 1.0), ("slab", 4.0), ("pml_left", 1.0), ("pml_right", 1.0)] {
            mg.set_material(s, name, "eps_r", eps);
            mg.set_material(s, name, "mu_r", 1.0);
        }
        for (name, dir) in [("pml_left", PmlDir::Minus), ("pml_right", PmlDir::Plus)] {
            for d in mg.domains(name) {
                s.materials.pml.insert(d, Pml { x: dir, y: PmlDir::Off, z: PmlDir::Off, strength: 4.0, order: 2.0, re: 1.0 });
            }
        }
    }
    p.spec_mut(Physics::Emw).study = Study::Modes { freq: 12e9, n_target: 1.9, evanescent: false, n_modes: 6 };
    p
}

/// A 0.4 x 0.2 plate with a hole of radius 0.04 at the centre, pulled along x.
pub fn plate_with_hole() -> Project {
    let mut m = Model::new_2d();
    let r = m.rect("plate_box", (0.0, 0.0), (0.4, 0.2));
    let hole = m.circle("hole", (0.2, 0.1), 0.04);
    let p = m.difference("plate", r, &[hole]);
    m.overlay(&[("plate", p)]);
    m.size("size_hole", selector("kind(circle)"), Some(0.004.into()), None, None);
    let (g, mg) = named(m, 0.012);
    let mut p = base(g, Physics::Solid);
    let s = p.spec_mut(Physics::Solid);
    s.plane = Plane::Stress;
    s.bcs = vec![
        bc(&sel(&mg, "outer & plane(x = 0)"), BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] }),
        bc(&sel(&mg, "outer & plane(z = 0)"), BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] }),
        bc(&sel(&mg, "outer & plane(x = 0.4)"), BcKind::Load([1e6, 0.0, 0.0])),
    ];
    p
}

/// The substrate (a 1 x 0.5 block with absorbing layers on the left, right and bottom) carrying a thin beam on its top face.
pub fn anchor_loss() -> Project {
    let mut m = Model::new_2d();
    let sub = m.rect("substrate_box", (0.0, 0.0), (1.0, 0.5));
    let pl = m.rect("pml_left_box", (0.0, 0.15), (0.2, 0.5));
    let pr = m.rect("pml_right_box", (0.8, 0.15), (1.0, 0.5));
    let pb = m.rect("pml_bottom_box", (0.2, 0.0), (0.8, 0.15));
    let pbl = m.rect("pml_bottom_left_box", (0.0, 0.0), (0.2, 0.15));
    let pbr = m.rect("pml_bottom_right_box", (0.8, 0.0), (1.0, 0.15));
    let beam = m.rect("beam_box", (0.49, 0.5), (0.51, 0.8));
    m.overlay(&[("substrate", sub), ("pml_left", pl), ("pml_right", pr), ("pml_bottom", pb), ("pml_bottom_left", pbl), ("pml_bottom_right", pbr), ("beam", beam)]);
    let (g, mg) = named(m, 0.03);
    let mut p = base(g, Physics::Solid);
    let pieces = [
        ("pml_left", PmlDir::Minus, PmlDir::Off),
        ("pml_right", PmlDir::Plus, PmlDir::Off),
        ("pml_bottom", PmlDir::Off, PmlDir::Minus),
        ("pml_bottom_left", PmlDir::Minus, PmlDir::Minus),
        ("pml_bottom_right", PmlDir::Plus, PmlDir::Minus),
    ];
    for s in &mut p.specs {
        for (name, x, y) in pieces {
            for d in mg.domains(name) {
                s.materials.pml.insert(d, Pml { x, y, z: PmlDir::Off, strength: 3.0, order: 2.0, re: 1.0 });
            }
        }
    }
    // clamped on the outer edge of the layers: every boundary below the substrate surface
    let outer = sel(&mg, "(outer & in_box([-1, 0, -1], [2, 0, 0.5000000001], within)) - plane(z = 0.5)");
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&outer, fixed_all())];
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    s.plane = Plane::Strain;
    p
}

/// A plate 0.4 x 0.2 x 0.04 with a drilled hole (axis along y, through the thickness): the section with the hole swept along y.
pub fn drilled_plate() -> Project {
    let mut m = Model::new_2d();
    let r = m.rect("plate_box", (0.0, 0.0), (0.4, 0.2));
    let hole = m.circle("hole", (0.2, 0.1), 0.04);
    let p = m.difference("plate", r, &[hole]);
    m.overlay(&[("plate", p)]);
    m.size("size_hole", selector("kind(circle)"), Some(0.01.into()), None, None);
    let (g, mg, _start, _end) = swept(m, 0.02, 0.04, 4);
    let mut p = base(g, Physics::Solid);
    p.spec_mut(Physics::Solid).bcs = vec![bc(&sel(&mg, "outer & plane(x = 0)"), fixed_all()), bc(&sel(&mg, "outer & plane(x = 0.4)"), BcKind::Load([1e6, 0.0, 0.0]))];
    p
}

/// The cross-section of a silicon strip (450 nm x 230 nm) on a sapphire substrate under air, at 1550 nm:
/// optical modes with PML on all sides (regions `air`, `sapphire`, `silicon` and the absorbing pieces `<region>.pml_xmin`, ...) and
/// acoustic modes of the solid (silicon and sapphire, rotated crystal, PML in the sapphire) at a small
/// out-of-plane wave number (forward Brillouin scattering).
pub fn sos_strip() -> Project {
    use femcore::material::*;
    let um = 1e-6;
    let (xin, hin, ztop, lp) = (6.0 * um, 10.0 * um, 3.0 * um, 1.0 * um);
    let (w_si, h_si, zone) = (450e-9, 230e-9, 1.5 * um);
    let x = xin + lp;
    let mut m = Model::new_2d();
    let air = m.rect("air_box", (-x, -hin - lp), (x, ztop + lp));
    let sapphire = m.rect("sapphire_box", (-x, -hin - lp), (x, 0.0));
    let sap_near = m.rect("sapphire_near_box", (-zone, -zone), (zone, 0.0));
    let air_near = m.rect("air_near_box", (-zone, 0.0), (zone, zone));
    let si = m.rect("silicon_box", (-w_si / 2.0, 0.0), (w_si / 2.0, h_si));
    let window = m.rect("window", (-xin, -hin), (xin, ztop));
    m.overlay(&[("air", air), ("sapphire", sapphire), ("sapphire_near", sap_near), ("air_near", air_near), ("silicon", si)]);
    m.pml("pml", window, &[], lp, 3.0, 2.0);
    m.size("size_near", selector("region(sapphire_near) | region(air_near)"), Some((0.08 * um).into()), None, None);
    m.size("size_si", selector("region(silicon)"), Some(40e-9.into()), None, None);
    let (g, mg) = named(m, 0.4 * um);
    let mut p = base(g, Physics::Emw);
    let (sapphire, air, silicon) = (mg.domains("sapphire*, sapphire_near*"), mg.domains("air*, air_near*"), mg.domains("silicon"));
    let pml = mg.pml_materials(3.0);
    // optical: sapphire n = 1.746, silicon 3.45, air 1
    {
        let e = p.spec_mut(Physics::Emw);
        for (ds, n) in [(&sapphire, 1.746), (&silicon, 3.45), (&air, 1.0)] {
            for &d in ds {
                e.materials.set(d, "eps_r", n * n);
                e.materials.set(d, "mu_r", 1.0);
            }
        }
        e.materials.pml.extend(pml.clone());
        e.study = Study::Modes { freq: 299792458.0 / 1550e-9, n_target: 2.4, evanescent: false, n_modes: 4 };
    }
    // acoustic: silicon strip and sapphire block (rotated), the layers above the surface are not solid
    {
        let so = p.spec_mut(Physics::Solid);
        so.inactive = air.clone();
        for d in sapphire.iter().chain(&silicon).chain(&air) {
            so.materials.set(*d, "E", 100e9);
            so.materials.set(*d, "nu", 0.25);
            so.materials.set(*d, "rho", 3980.0);
        }
        let saph = Anisotropy { stiffness: Some(Stiffness::Trigonal { c11: 497.73e9, c12: 162.8e9, c13: 116.0e9, c14: -21.9e9, c33: 500.9e9, c44: 146.8e9 }), euler_deg: femcore::material::euler_from_rotation(&femcore::chips::crystal_rotations(0.0).1), ..Default::default() };
        for &d in &sapphire {
            so.materials.aniso.insert(d, saph.clone());
        }
        for &d in &silicon {
            so.materials.aniso.insert(d, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 166e9, c12: 64e9, c44: 79e9 }), euler_deg: [0.0; 3], ..Default::default() });
            so.materials.set(d, "rho", 2329.0);
        }
        so.materials.pml.extend(pml.into_iter().filter(|(d, _)| sapphire.contains(d)));
        // clamped: the outer sides and bottom of the sapphire block
        let outer = sel(&mg, "outer & in_box([-1, 0, -1], [1, 0, 1e-12], within)");
        so.bcs = vec![Bc { labels: outer, kind: fixed_all() }];
        so.study = Study::Eigenfrequency { n_modes: 8, f_search: 9.47e9 };
        so.out_of_plane = Some(2.3172 * 2.0 * std::f64::consts::PI * 9.47e9 / 299792458.0);
    }
    p
}

/// A 40 x 30 x 20 mm air-filled cavity with perfectly conducting walls (the 40 x 20 mm section swept 30 mm along y); the modes are TE101
/// at 6.2457 GHz, then 8.3795, 9.0076 (twice) GHz and so on (f = c/2 sqrt((m/a)^2 + (n/b)^2 + (p/d)^2)).
pub fn box_cavity() -> Project {
    let (g, _mg, _s, _e) = swept(plate(0.04, 0.02), 0.005, 0.03, 6);
    let mut p = base(g, Physics::Emw);
    let e = p.spec_mut(Physics::Emw);
    e.materials.default.insert("eps_r".into(), 1.0);
    e.materials.default.insert("mu_r".into(), 1.0);
    e.study = Study::Eigenfrequency { n_modes: 6, f_search: 6.0e9 };
    p
}

/// Backward Brillouin scattering in the rib waveguide of the 300TE ring chip (0.5 um wide, 0.3 um high rib on a
/// 0.13 um slab, silicon on r-plane sapphire, 1550 nm), in the wafer frame (x lateral, y propagation, z normal; the
/// cross-section is the x-z plane): the optical TE0 mode, the phonons at q = 2 beta and the photoelastic and
/// moving-boundary gain, with the full anisotropic stiffness, photoelastic and permittivity tensors of the crystals.
/// The geometry is a `geom` model (`femcore::chips::chip_geometry`: regions air, sapphire, slab, rib and their absorbing layers).
pub fn brillouin_rib() -> Project {
    let chip = femcore::chips::CHIPS[0];
    let g = femcore::chips::chip_mesh(&chip, 1.0).expect("the chip mesh");
    let spec = femcore::chips::chip_problem(&chip, 0.0, 24, 18.35e9, &g);
    let geometry = GeometrySpec::model(femcore::chips::chip_geometry(&chip, 1.0).list, 0.12e-6);
    let mut p = base(geometry, Physics::Brillouin);
    *p.spec_mut(Physics::Brillouin) = spec;
    p
}

/// A project for the optical modes of a `SoiWire`: its model and its mode study, in the wafer frame.
fn soi_project(wire: &femcore::usecases::soi::SoiWire) -> Project {
    let g = wire.mesh().expect("the wire mesh");
    let geometry = GeometrySpec::model(wire.model().list, wire.h_far);
    let mut p = base(geometry, Physics::Emw);
    *p.spec_mut(Physics::Emw) = wire.spec(&g);
    p
}

/// Modes of a 500 nm x 300 nm silicon wire in oxide on 300 nm silicon-on-insulator at 1550 nm (`femcore::usecases::soi`): the quasi-TE and
/// quasi-TM modes with their Stokes parameters and Poincare points in the Analysis tab. Regions `BOX`, `Si`, `clad` (and `_near` zones).
pub fn soi_wire() -> Project {
    soi_project(&femcore::usecases::soi::SoiWire::default())
}

/// The same 300 nm SOI film with a 300 nm wide wire whose left wall is sloped to 75 degrees: the quasi-TE/quasi-TM pair of the square core is
/// mixed into two hybrid modes, tilted by about 45 degrees on the Poincare sphere.
pub fn soi_hybrid_wire() -> Project {
    soi_project(&femcore::usecases::soi::SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), ..Default::default() })
}

/// One period of a Bragg grating in the SOI wire as a 3D Floquet cell at the band edge (phase pi), with the coarse mesh of
/// `BraggGrating::coarse` (about 40 s): the two lowest modes bound the stop band. The x-z cross-section is a `geom` model swept along y,
/// the propagation axis.
fn bragg_project(grating: &femcore::usecases::bragg::BraggGrating) -> Project {
    let cell = grating.mesh().expect("the grating mesh");
    let f_search = 0.95 * grating.bragg_frequency_estimate().expect("the Bragg frequency of the wire");
    let geometry = GeometrySpec::model_extruded(grating.cross_section().list, grating.h_far, grating.layers());
    let mut p = base(geometry, Physics::Emw);
    *p.spec_mut(Physics::Emw) = grating.spec(&cell, f_search);
    p.active = Physics::Emw;
    p
}

/// Symmetric Bragg grating: 80 nm deep teeth on both sides of a 500 nm x 300 nm wire, 280 nm period, 50 % duty (`femcore::usecases::bragg`).
pub fn bragg_symmetric() -> Project {
    bragg_project(&femcore::usecases::bragg::BraggGrating { kind: femcore::usecases::bragg::GratingKind::Symmetric, ..femcore::usecases::bragg::BraggGrating::coarse() })
}

/// Asymmetric Bragg grating: the same teeth on one side only, so the mirror symmetry of the wire is gone and the stop band is narrower.
pub fn bragg_asymmetric() -> Project {
    bragg_project(&femcore::usecases::bragg::BraggGrating { kind: femcore::usecases::bragg::GratingKind::Asymmetric, ..femcore::usecases::bragg::BraggGrating::coarse() })
}
