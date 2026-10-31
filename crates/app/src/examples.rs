//! Ready-made projects (geometry, materials, boundary conditions, studies) to start from.

use crate::colormap::Colormap;
use crate::state::*;
use femcore::problem::*;

fn base(geometry: GeometrySpec, active: Physics) -> Project {
    let mut specs: Vec<ProblemSpec> = Physics::ALL.iter().map(|p| ProblemSpec::new(*p)).collect();
    for s in &mut specs {
        for (k, v) in [("eps_r", 1.0), ("mu_r", 1.0), ("k", 1.0), ("rho", 1000.0), ("Cp", 500.0), ("E", 200e9), ("nu", 0.3), ("alpha", 12e-6)] {
            s.materials.default.insert(k.into(), v);
        }
        s.materials.default.insert("rho".into(), if s.physics == Physics::Solid { 7850.0 } else { 1000.0 });
    }
    Project { version: 1, geometry, source: Source::Generated, active, specs, colormap: Colormap::Viridis }
}

/// Turns a project written for the 2D model frame (mesh x-y, z out of plane) into the wafer frame (mesh x-z, y out of
/// plane): constraints, loads and crystal orientations keep their physical meaning.
fn to_wafer(mut p: Project) -> Project {
    let is_2d = p.geometry.geometry.dim() == 2;
    for s in &mut p.specs {
        s.materials.frame = Frame::Wafer;
        if !is_2d {
            continue;
        }
        // model components (1, 2, 3) = (x, z, -y): wafer = (c1, -c3, c2)
        for bc in &mut s.bcs {
            match &mut bc.kind {
                BcKind::Fixed { fixed, values } => {
                    *fixed = [fixed[0], fixed[2], fixed[1]];
                    *values = [values[0], -values[2], values[1]];
                }
                BcKind::Load(f) => *f = [f[0], -f[2], f[1]],
                _ => {}
            }
        }
        for a in s.materials.aniso.values_mut() {
            // wafer orientation = S^-1 R_model
            let r = femcore::material::rotation(a.euler_deg);
            let st = femcore::problem::WAFER_TO_MODEL;
            let sinv = [[st[0][0], st[1][0], st[2][0]], [st[0][1], st[1][1], st[2][1]], [st[0][2], st[1][2], st[2][2]]];
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        m[i][j] += sinv[i][k] * r[k][j];
                    }
                }
            }
            a.euler_deg = femcore::material::euler_from_rotation(&m);
        }
    }
    p
}

fn bc(labels: &[i32], kind: BcKind) -> Bc {
    Bc { labels: labels.to_vec(), kind }
}

fn set(p: &mut Project, phys: Physics, props: &[(&str, f64)]) {
    for (k, v) in props {
        p.spec_mut(phys).materials.default.insert((*k).into(), *v);
    }
}

pub const NAMES: [&str; 19] = [
    "Parallel-plate capacitor with two dielectrics (electrostatics, 2D)",
    "Cantilever beam (solid mechanics, static + eigenfrequencies, 2D)",
    "Heated plate with convective cooling (heat transfer, 2D, transient)",
    "Thermal expansion of a bar (heat + solid, 2D)",
    "Block cantilever (solid mechanics, 3D)",
    "WR-90 rectangular waveguide modes (EM waves, 2D cross-section)",
    "Dielectric-loaded waveguide (EM waves, 2D cross-section)",
    "Sapphire cantilever, anisotropic crystal (solid mechanics, 3D eigenfrequencies)",
    "Dielectric slab guide with PML ends (EM waves, complex modes)",
    "Plate with a hole, meshed from shapes (solid mechanics, stress concentration)",
    "Beam on a PML substrate: anchor loss (solid eigenfrequencies, complex frequency and Q)",
    "Plate with a drilled hole in 3D, from the 3D mesher (solid mechanics, stress concentration)",
    "Silicon strip on sapphire at 1550 nm: optical modes (EM waves) and acoustic modes with out-of-plane wave number (solid mechanics)",
    "Rectangular cavity: 3D electromagnetic eigenmodes (40 x 30 x 20 mm, perfect conductor walls)",
    "Brillouin gain of a silicon rib on sapphire (300TE ring chip, backward, 1550 nm): optical mode, phonons, photoelastic and moving-boundary gain",
    "Silicon wire on 300 nm SOI at 1550 nm: TE0 and TM0 modes with Stokes parameters and Poincare sphere (EM waves, 2D cross-section)",
    "300 nm SOI wire with a sloped left wall: hybrid modes tilted on the Poincare sphere (EM waves, 2D cross-section)",
    "Symmetric Bragg grating in a 300 nm SOI wire: band-edge modes of a 3D Floquet cell (teeth on both sides)",
    "Asymmetric Bragg grating in a 300 nm SOI wire: band-edge modes of a 3D Floquet cell (teeth on one side)",
];

/// An empty project to build a simulation from: the chosen physics active, default materials, no boundary conditions, a
/// unit rectangle or block (`structured`) or a one-shape/one-body geometry for the shape and solid meshers, wafer axes.
pub fn blank(physics: Physics, three_d: bool, structured: bool) -> Project {
    let three_d = three_d && physics != Physics::Brillouin;
    let g = match (three_d, structured) {
        (false, true) => GeometrySpec::rectangle(1.0, 1.0, 20, 20),
        (false, false) => GeometrySpec::shapes(),
        (true, true) => GeometrySpec::block(1.0, 1.0, 1.0, 8, 8, 8),
        (true, false) => GeometrySpec::solids(),
    };
    let mut p = base(g, physics);
    if three_d && physics == Physics::Emw {
        p.spec_mut(physics).study = Study::Eigenfrequency { n_modes: 6, f_search: 10e9 };
    }
    to_wafer(p)
}

pub fn by_index(i: usize) -> Project {
    let p = by_index_model(i);
    // the 3D examples and the Brillouin example are described in the wafer frame already
    if i == 14 {
        let mut p = p;
        for s in &mut p.specs {
            s.materials.frame = Frame::Wafer;
        }
        return p;
    }
    to_wafer(p)
}

fn by_index_model(i: usize) -> Project {
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
    let mut g = GeometrySpec::rectangle(1.0, 0.5, 40, 20);
    g.domains = vec![Region { bbox: vec![0.5, 0.0, 1.0, 0.5], label: 2 }];
    let mut p = base(g, Physics::Electrostatics);
    let es = p.spec_mut(Physics::Electrostatics);
    es.materials.set(1, "eps_r", 2.0);
    es.materials.set(2, "eps_r", 8.0);
    es.bcs = vec![bc(&[1], BcKind::Potential(1.0)), bc(&[4], BcKind::Ground)];
    p
}

pub fn cantilever() -> Project {
    let mut p = base(GeometrySpec::rectangle(1.0, 0.1, 60, 6), Physics::Solid);
    p.spec_mut(Physics::Solid).bcs = vec![bc(&[1], BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] }), bc(&[4], BcKind::Load([0.0, -1e5, 0.0]))];
    p.spec_mut(Physics::Solid).plane = Plane::Stress;
    p
}

pub fn heated_plate() -> Project {
    let mut p = base(GeometrySpec::rectangle(1.0, 0.2, 40, 8), Physics::Heat);
    set(&mut p, Physics::Heat, &[("k", 50.0), ("rho", 7800.0), ("Cp", 470.0)]);
    let h = p.spec_mut(Physics::Heat);
    h.bcs = vec![bc(&[1], BcKind::Temperature(500.0)), bc(&[2, 3, 4], BcKind::Convection { h: 25.0, text: 293.15 })];
    h.study = Study::Transient { t_end: 20000.0, n_out: 20, t0: "293.15".into() };
    p
}

pub fn thermal_bar() -> Project {
    let mut p = base(GeometrySpec::rectangle(1.0, 0.1, 50, 5), Physics::Solid);
    set(&mut p, Physics::Heat, &[("k", 160.0)]);
    set(&mut p, Physics::Solid, &[("E", 70e9), ("nu", 0.33), ("alpha", 23e-6), ("rho", 2700.0)]);
    p.spec_mut(Physics::Heat).bcs = vec![bc(&[1], BcKind::Temperature(293.15)), bc(&[4], BcKind::Temperature(393.15))];
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&[1], BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] })];
    s.thermal = Some(ThermalExpansion { source: ThermalSource::FromHeat, t_ref: 293.15 });
    s.plane = Plane::Stress;
    p
}

pub fn block_cantilever() -> Project {
    let mut p = base(GeometrySpec::block(1.0, 0.1, 0.1, 16, 3, 3), Physics::Solid);
    p.spec_mut(Physics::Solid).bcs = vec![bc(&[1], BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] }), bc(&[6], BcKind::Load([0.0, 0.0, -1e5]))];
    p
}

pub fn waveguide() -> Project {
    let mut p = base(GeometrySpec::rectangle(22.86e-3, 10.16e-3, 24, 12), Physics::Emw);
    p.spec_mut(Physics::Emw).study = Study::Modes { freq: 20e9, n_target: 0.9, evanescent: false, n_modes: 4 };
    p
}

pub fn loaded_waveguide() -> Project {
    let mut g = GeometrySpec::rectangle(22.86e-3, 10.16e-3, 36, 12);
    g.domains = vec![Region { bbox: vec![0.0, 0.0, 11.43e-3, 10.16e-3], label: 2 }];
    let mut p = base(g, Physics::Emw);
    let e = p.spec_mut(Physics::Emw);
    e.materials.set(1, "eps_r", 1.0);
    e.materials.set(2, "eps_r", 4.0);
    e.study = Study::Modes { freq: 8e9, n_target: 1.7, evanescent: false, n_modes: 3 };
    p
}

pub fn sapphire_cantilever() -> Project {
    let mut p = base(GeometrySpec::block(1.0e-3, 0.2e-3, 0.1e-3, 20, 4, 2), Physics::Solid);
    let (st, eps, k, alpha, rho) = femcore::material::sapphire();
    let an = femcore::material::Anisotropy { stiffness: Some(st), conductivity: Some(k), permittivity: Some(eps), expansion: Some(alpha), euler_deg: [0.0, 90.0, 0.0] };
    for s in &mut p.specs {
        s.materials.aniso.insert(1, an.clone());
        s.materials.default.insert("rho".into(), rho);
    }
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&[1], BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] })];
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    p
}

pub fn pml_slab() -> Project {
    let (lx, ly) = (40e-3, 10e-3);
    let mut g = GeometrySpec::rectangle(lx, ly, 80, 8);
    g.domains = vec![
        Region { bbox: vec![15e-3, 0.0, 25e-3, ly], label: 2 },
        Region { bbox: vec![0.0, 0.0, 8e-3, ly], label: 3 },
        Region { bbox: vec![32e-3, 0.0, lx, ly], label: 4 },
    ];
    let mut p = base(g, Physics::Emw);
    for s in &mut p.specs {
        s.materials.set(1, "eps_r", 1.0);
        s.materials.set(2, "eps_r", 4.0);
        s.materials.set(3, "eps_r", 1.0);
        s.materials.set(4, "eps_r", 1.0);
        s.materials.pml.insert(3, Pml { x: PmlDir::Minus, y: PmlDir::Off, z: PmlDir::Off, strength: 4.0, order: 2.0, re: 1.0 });
        s.materials.pml.insert(4, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 4.0, order: 2.0, re: 1.0 });
    }
    p.spec_mut(Physics::Emw).study = Study::Modes { freq: 12e9, n_target: 1.9, evanescent: false, n_modes: 6 };
    p
}

pub fn plate_with_hole() -> Project {
    use femcore::mesher::*;
    let spec = Spec2d {
        parts: vec![
            Part { shape: Shape::Rect { x0: 0.0, y0: 0.0, x1: 0.4, y1: 0.2 }, domain: 1, h: None },
            Part { shape: Shape::Circle { cx: 0.2, cy: 0.1, r: 0.04 }, domain: 0, h: Some(0.004) },
        ],
        h: 0.012,
        min_angle_deg: 25.0,
    };
    let g = GeometrySpec { geometry: Geometry::Shapes { spec, extrude: None }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Solid);
    let s = p.spec_mut(Physics::Solid);
    s.plane = Plane::Stress;
    s.bcs = vec![bc(&[1], BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] }), bc(&[2], BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] }), bc(&[4], BcKind::Load([1e6, 0.0, 0.0]))];
    p
}

pub fn anchor_loss() -> Project {
    use femcore::mesher::*;
    let rect = |x0, y0, x1, y1, d| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: d, h: None };
    let spec = Spec2d {
        parts: vec![
            rect(0.0, 0.0, 1.0, 0.5, 2),
            rect(0.0, 0.15, 0.2, 0.5, 3),
            rect(0.8, 0.15, 1.0, 0.5, 4),
            rect(0.2, 0.0, 0.8, 0.15, 5),
            rect(0.0, 0.0, 0.2, 0.15, 6),
            rect(0.8, 0.0, 1.0, 0.15, 7),
            rect(0.49, 0.5, 0.51, 0.8, 1),
        ],
        h: 0.03,
        min_angle_deg: 25.0,
    };
    // the outer edge of the layers: every boundary facet below the substrate surface
    let mut outer: Vec<i32> = vec![];
    if let Ok(m) = mesh_2d(&spec) {
        for f in 0..m.n_facets() {
            let (a, b) = (m.point(m.facet(f)[0] as usize), m.point(m.facet(f)[1] as usize));
            if (a[1] + b[1]) / 2.0 < 0.5 - 1e-9 {
                outer.push(m.facet_labels[f]);
            }
        }
    }
    outer.sort();
    outer.dedup();
    let g = GeometrySpec { geometry: Geometry::Shapes { spec, extrude: None }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Solid);
    for s in &mut p.specs {
        let st = |x, y| Pml { x, y, z: PmlDir::Off, strength: 3.0, order: 2.0, re: 1.0 };
        s.materials.pml.insert(3, st(PmlDir::Minus, PmlDir::Off));
        s.materials.pml.insert(4, st(PmlDir::Plus, PmlDir::Off));
        s.materials.pml.insert(5, st(PmlDir::Off, PmlDir::Minus));
        s.materials.pml.insert(6, st(PmlDir::Minus, PmlDir::Minus));
        s.materials.pml.insert(7, st(PmlDir::Plus, PmlDir::Minus));
    }
    let s = p.spec_mut(Physics::Solid);
    s.bcs = vec![bc(&outer, BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] })];
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    s.plane = Plane::Strain;
    p
}

pub fn drilled_plate() -> Project {
    use femcore::mesher3d::*;
    let spec = Spec3d {
        bodies: vec![
            Body3d { solid: Solid::Box { min: [0.0, 0.0, 0.0], max: [0.4, 0.2, 0.04] }, domain: 1, h: None },
            Body3d { solid: Solid::Cylinder { base: [0.2, 0.1, -0.01], axis: 2, r: 0.04, len: 0.06 }, domain: 0, h: None },
        ],
        h: 0.01,
        alpha: 0.3,
        snap_features: true, min_layers: 2.0
    };
    let g = GeometrySpec { geometry: Geometry::Solids { spec }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Solid);
    let s = p.spec_mut(Physics::Solid);
    // faces: 1 x=0, 2 x=0.4, 7 the wall of the hole
    s.bcs = vec![bc(&[1], BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] }), bc(&[2], BcKind::Load([1e6, 0.0, 0.0]))];
    p
}

/// The cross-section of a silicon strip (450 nm x 230 nm) on a sapphire substrate under air, at 1550 nm:
/// optical modes with PML on all sides (domain 3 air, 1 sapphire, 2 silicon, 4..11 the absorbing layers) and
/// acoustic modes of the solid (silicon and sapphire, rotated crystal, PML in the sapphire) at a small
/// out-of-plane wave number (forward Brillouin scattering).
pub fn sos_strip() -> Project {
    use femcore::material::*;
    use femcore::mesher::*;
    let um = 1e-6;
    let (w_sapph, h_sapph, h_air, w_pml) = (12.0 * um, 10.0 * um, 3.0 * um, 1.0 * um);
    let (w_si, h_si) = (450e-9, 230e-9);
    let (x0, x1, y0, y1) = (-w_sapph / 2.0, w_sapph / 2.0, -h_sapph / 2.0, h_sapph / 2.0 + h_air);
    let ys = h_sapph / 2.0;
    let rect = |xa: f64, ya: f64, xb: f64, yb: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: xa, y0: ya, x1: xb, y1: yb }, domain: d, h };
    let zone = 1.5 * um;
    let spec2 = Spec2d {
        parts: vec![
            rect(x0 - w_pml, y0 - w_pml, x1 + w_pml, y1 + w_pml, 3, None),
            rect(x0, y0, x1, ys, 1, None),
            rect(-zone, ys - zone, zone, ys, 1, Some(0.08 * um)),
            rect(-zone, ys, zone, ys + zone, 3, Some(0.08 * um)),
            rect(-w_si / 2.0, ys, w_si / 2.0, ys + h_si, 2, Some(40e-9)),
            // absorbing layers: sapphire side below the substrate surface, air above it
            rect(x0 - w_pml, y0, x0, ys, 4, None),
            rect(x1, y0, x1 + w_pml, ys, 5, None),
            rect(x0, y0 - w_pml, x1, y0, 6, None),
            rect(x0, y1, x1, y1 + w_pml, 7, None),
            rect(x0 - w_pml, y0 - w_pml, x0, y0, 8, None),
            rect(x1, y0 - w_pml, x1 + w_pml, y0, 9, None),
            rect(x0 - w_pml, ys, x0, y1, 12, None),
            rect(x1, ys, x1 + w_pml, y1, 13, None),
            rect(x0 - w_pml, y1, x0, y1 + w_pml, 10, None),
            rect(x1, y1, x1 + w_pml, y1 + w_pml, 11, None),
        ],
        h: 0.4 * um,
        min_angle_deg: 25.0,
    };
    // boundaries of the solid: the outer sides and bottom of the sapphire block
    let mut outer: Vec<i32> = vec![];
    if let Ok(m) = mesh_2d(&spec2) {
        for f in 0..m.n_facets() {
            let (pa, pb) = (m.point(m.facet(f)[0] as usize), m.point(m.facet(f)[1] as usize));
            let (mx, my) = (0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]));
            let side = (mx - (x0 - w_pml)).abs() < 1e-9 || (mx - (x1 + w_pml)).abs() < 1e-9 || (my - (y0 - w_pml)).abs() < 1e-9;
            if side && my < ys + 1e-9 {
                outer.push(m.facet_labels[f]);
            }
        }
    }
    outer.sort();
    outer.dedup();
    let g = GeometrySpec { geometry: Geometry::Shapes { spec: spec2, extrude: None }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Emw);
    let pml = |x, y, strength| Pml { x, y, z: PmlDir::Off, strength, order: 2.0, re: 1.0 };
    let absorbing = |strength: f64| -> Vec<(i32, Pml)> {
        vec![
            (4, pml(PmlDir::Minus, PmlDir::Off, strength)),
            (5, pml(PmlDir::Plus, PmlDir::Off, strength)),
            (6, pml(PmlDir::Off, PmlDir::Minus, strength)),
            (7, pml(PmlDir::Off, PmlDir::Plus, strength)),
            (8, pml(PmlDir::Minus, PmlDir::Minus, strength)),
            (9, pml(PmlDir::Plus, PmlDir::Minus, strength)),
            (10, pml(PmlDir::Minus, PmlDir::Plus, strength)),
            (11, pml(PmlDir::Plus, PmlDir::Plus, strength)),
            (12, pml(PmlDir::Minus, PmlDir::Off, strength)),
            (13, pml(PmlDir::Plus, PmlDir::Off, strength)),
        ]
    };
    // optical: sapphire n = 1.746, silicon 3.45, air 1
    {
        let e = p.spec_mut(Physics::Emw);
        for (d, n) in [(1, 1.746), (2, 3.45), (3, 1.0), (4, 1.746), (5, 1.746), (6, 1.746), (7, 1.0), (8, 1.746), (9, 1.746), (10, 1.0), (11, 1.0), (12, 1.0), (13, 1.0)] {
            e.materials.set(d, "eps_r", n * n);
            e.materials.set(d, "mu_r", 1.0);
        }
        for (d, pm) in absorbing(3.0) {
            e.materials.pml.insert(d, pm);
        }
        e.study = Study::Modes { freq: 299792458.0 / 1550e-9, n_target: 2.4, evanescent: false, n_modes: 4 };
    }
    // acoustic: silicon strip and sapphire block (rotated), the layers above the surface are not solid
    {
        let so = p.spec_mut(Physics::Solid);
        so.inactive = vec![3, 7, 10, 11, 12, 13];
        for d in 1..=13 {
            so.materials.set(d, "E", 100e9);
            so.materials.set(d, "nu", 0.25);
            so.materials.set(d, "rho", 3980.0);
        }
        let saph = Anisotropy { stiffness: Some(Stiffness::Trigonal { c11: 497.73e9, c12: 162.8e9, c13: 116.0e9, c14: -21.9e9, c33: 500.9e9, c44: 146.8e9 }), euler_deg: femcore::problem::euler_to_model_frame(femcore::material::euler_from_rotation(&femcore::chips::crystal_rotations(0.0).1)), ..Default::default() };
        for d in [1, 4, 5, 6, 8, 9] {
            so.materials.aniso.insert(d, saph.clone());
        }
        so.materials.aniso.insert(2, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 166e9, c12: 64e9, c44: 79e9 }), euler_deg: [0.0; 3], ..Default::default() });
        so.materials.set(2, "rho", 2329.0);
        for (d, pm) in absorbing(3.0) {
            if [4, 5, 6, 8, 9].contains(&d) {
                so.materials.pml.insert(d, pm);
            }
        }
        so.bcs = vec![Bc { labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        so.study = Study::Eigenfrequency { n_modes: 8, f_search: 9.47e9 };
        so.out_of_plane = Some(2.3172 * 2.0 * std::f64::consts::PI * 9.47e9 / 299792458.0);
    }
    p
}

/// A 40 x 30 x 20 mm air-filled cavity with perfectly conducting walls; the modes are TE101 at 6.2457 GHz,
/// then 8.3795, 9.0076 (twice) GHz and so on (f = c/2 sqrt((m/a)^2 + (n/b)^2 + (p/d)^2)).
pub fn box_cavity() -> Project {
    let g = GeometrySpec::block(0.04, 0.03, 0.02, 8, 6, 4);
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
/// See `femcore::chips` for the layers (domains 1 sapphire, 2 slab, 3 rib, 4 air, 11-15, 21-22, 41-45 absorbing layers).
pub fn brillouin_rib() -> Project {
    let chip = femcore::chips::CHIPS[0];
    let (_, spec) = femcore::chips::chip_model(&chip, 0.0, 1.0, 24, 18.35e9).expect("the chip mesh");
    let g = GeometrySpec { geometry: Geometry::Shapes { spec: femcore::chips::chip_shapes(&chip, 1.0), extrude: None }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Brillouin);
    *p.spec_mut(Physics::Brillouin) = spec;
    p
}

/// A project for the optical modes of a `SoiWire`: its shape mesh and its mode study, in the wafer frame.
fn soi_project(wire: &femcore::usecases::soi::SoiWire) -> Project {
    let g = GeometrySpec { geometry: Geometry::Shapes { spec: wire.shapes(), extrude: None }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Emw);
    *p.spec_mut(Physics::Emw) = wire.spec();
    p
}

/// Modes of a 500 nm x 300 nm silicon wire in oxide on 300 nm silicon-on-insulator at 1550 nm (`femcore::usecases::soi`): the quasi-TE and
/// quasi-TM modes with their Stokes parameters and Poincare points in the Analysis tab. Domains 1 buried oxide, 2 silicon, 3 cladding.
pub fn soi_wire() -> Project {
    soi_project(&femcore::usecases::soi::SoiWire::default())
}

/// The same 300 nm SOI film with a 300 nm wide wire whose left wall is sloped to 75 degrees: the quasi-TE/quasi-TM pair of the square core is
/// mixed into two hybrid modes, tilted by about 45 degrees on the Poincare sphere.
pub fn soi_hybrid_wire() -> Project {
    soi_project(&femcore::usecases::soi::SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), ..Default::default() })
}

/// One period of a Bragg grating in the SOI wire as a 3D Floquet cell at the band edge (phase pi), with the coarse mesh of
/// `BraggGrating::coarse` (about 40 s): the two lowest modes bound the stop band. z is the propagation axis, x lateral, y normal.
fn bragg_project(grating: &femcore::usecases::bragg::BraggGrating) -> Project {
    let (_, src, dst) = grating.mesh().expect("the grating mesh");
    let f_search = 0.95 * grating.bragg_frequency_estimate().expect("the Bragg frequency of the wire");
    let g = GeometrySpec { geometry: Geometry::Shapes { spec: grating.cross_section(), extrude: Some(grating.layers()) }, domains: vec![], boundaries: vec![] };
    let mut p = base(g, Physics::Emw);
    *p.spec_mut(Physics::Emw) = grating.spec(src, dst, f_search);
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
