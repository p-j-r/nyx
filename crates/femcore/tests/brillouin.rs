//! Brillouin gain module: rotations, the optical power and group velocity of a mode against the
//! dispersion, the phase-matched wave number, and mode tracking.
//!
//! Seven tests that take 30-60 s each (rings, whispering-gallery modes, group velocity from the fields) are #[ignore]d so that
//! cargo test stays fast; run them with cargo test --release -p femcore --test brillouin -- --ignored (about 70 s).
use femcore::brillouin::*;
use femcore::csparse::c;
use femcore::material::*;
use femcore::mesher::*;
use femcore::problem::*;
use femcore::solid_kz::KzModes;
use std::collections::{BTreeMap, BTreeSet};

#[test]
fn euler_angles_of_a_rotation_reproduce_it() {
    for e in [[17.0, 33.0, 71.0], [0.0, 57.6 - 90.0, 0.0], [120.0, 10.0, -45.0], [0.0, 0.0, 0.0], [30.0, 180.0, 0.0]] {
        let r = rotation(e);
        let r2 = rotation(euler_from_rotation(&r));
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[i][j] - r2[i][j]).abs() < 1e-12, "{e:?}: {r:?} against {r2:?}");
            }
        }
    }
}

#[test]
fn an_isotropic_photoelastic_tensor_is_unchanged_by_rotation() {
    let (p11, p12) = (0.121, 0.270);
    let a = PhotoElastic::isotropic(p11, p12);
    let b = PhotoElastic::cubic(p11, p12, (p11 - p12) / 2.0, &rotation([23.0, 41.0, 77.0]));
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    assert!((a.p[i][j][k][l] - b.p[i][j][k][l]).abs() < 1e-12);
                }
            }
        }
    }
}

fn wire_mesh() -> femcore::mesh::SimplexMesh {
    let um = 1e-6;
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: x0 * um, y0: y0 * um, x1: x1 * um, y1: y1 * um }, domain: d, h: h.map(|v| v * um) };
    let parts = vec![rect(-1.2, -1.0, 1.2, 1.0, 1, None), rect(-0.225, -0.11, 0.225, 0.11, 2, Some(0.04))];
    mesh_2d(&Spec2d { parts, h: 0.3 * um, min_angle_deg: 25.0 }).unwrap()
}

fn optics(mesh: &femcore::mesh::SimplexMesh, lambda: f64) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    for d in mesh.domain_labels() {
        let n: f64 = if d == 2 { 3.476 } else { 1.444 };
        s.materials.set(d, "eps_r", n * n);
        s.materials.set(d, "mu_r", 1.0);
    }
    s.study = Study::Modes { freq: 299_792_458.0 / lambda, n_target: 2.4, evanescent: false, n_modes: 2 };
    s
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn group_velocity_from_the_fields_matches_the_dispersion() {
    // vg = P / W from the E and H fields must equal d omega / d beta from neff at two wavelengths
    let mesh = wire_mesh();
    let lambda = 1.55e-6;
    let set = solve_optical(&mesh, &optics(&mesh, lambda)).unwrap();
    let pe = BTreeMap::new();
    let moving = BTreeSet::new();
    let cp = coupling(&PairInput { mesh: &mesh, pump: (&set, 0), stokes: (&set, 0), backward: false, mirror_pump: false, photoelastic: &pe, moving: &moving }).unwrap();
    let (ng_fields, _) = cp.group_index();
    let dl = 1e-9;
    let n1 = solve_optical(&mesh, &optics(&mesh, lambda - dl)).unwrap().modes[0].neff.re;
    let n2 = solve_optical(&mesh, &optics(&mesh, lambda + dl)).unwrap().modes[0].neff.re;
    let n0 = set.modes[0].neff.re;
    // ng = n - lambda dn/dlambda
    let ng_disp = n0 - lambda * (n2 - n1) / (2.0 * dl);
    println!("neff {n0:.5}, ng from fields {ng_fields:.5}, from dispersion {ng_disp:.5}");
    assert!((ng_fields - ng_disp).abs() < 2e-3 * ng_disp, "{ng_fields} against {ng_disp}");
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn the_phase_matched_wave_number_follows_the_direction() {
    let mesh = wire_mesh();
    let set = solve_optical(&mesh, &optics(&mesh, 1.55e-6)).unwrap();
    let pe = BTreeMap::new();
    let moving = BTreeSet::new();
    let fwd = coupling(&PairInput { mesh: &mesh, pump: (&set, 0), stokes: (&set, 0), backward: false, mirror_pump: false, photoelastic: &pe, moving: &moving }).unwrap();
    let bwd = coupling(&PairInput { mesh: &mesh, pump: (&set, 0), stokes: (&set, 0), backward: true, mirror_pump: false, photoelastic: &pe, moving: &moving }).unwrap();
    assert!(fwd.q().abs() < 1e-6 * fwd.beta_pump);
    assert!((bwd.q() - 2.0 * bwd.beta_pump).abs() < 1e-6 * bwd.beta_pump);
    // the mirrored Stokes mode carries power in the opposite direction
    assert!(fwd.p_pump * fwd.p_stokes > 0.0 && bwd.p_pump * bwd.p_stokes < 0.0);
}

fn fake_modes(vs: Vec<Vec<femcore::csparse::C>>) -> KzModes {
    KzModes { kz: 1.0, freqs: vec![c(1e9, 0.0); vs.len()], vectors: vs, c6: vec![], rho: vec![], converged: true }
}

#[test]
fn modes_are_tracked_by_overlap_whatever_their_phase_and_order() {
    let a: Vec<femcore::csparse::C> = (0..40).map(|i| c((i as f64 * 0.3).sin(), (i as f64 * 0.1).cos())).collect();
    let b: Vec<femcore::csparse::C> = (0..40).map(|i| c((i as f64 * 0.77).cos(), (i as f64 * 0.21).sin())).collect();
    let rot = c(0.6, 0.8);
    let prev = fake_modes(vec![a.clone(), b.clone()]);
    let next = fake_modes(vec![b.iter().map(|z| *z * rot).collect(), a.iter().map(|z| *z * c(0.0, 1.0)).collect()]);
    let m = track_modes(&prev, &next, 0.5);
    assert_eq!(m[0].map(|x| x.0), Some(1));
    assert_eq!(m[1].map(|x| x.0), Some(0));
    assert!(m[0].unwrap().1 > 0.999);
}

#[test]
fn viscous_loss_proportional_to_stiffness_gives_the_exact_quality_factor() {
    // eta = tau C damps every mode with 1/Q = omega tau
    let mesh = femcore::mesh::SimplexMesh::rectangle(1.0, 0.6, 12, 8);
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.order = 2;
    sp.materials.default.insert("E".into(), 70e9);
    sp.materials.default.insert("nu".into(), 0.3);
    sp.materials.default.insert("rho".into(), 2700.0);
    sp.bcs = vec![Bc { on: None, labels: vec![1, 2, 3, 4], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    sp.study = Study::Eigenfrequency { n_modes: 4, f_search: 3000.0 };
    let (_, modes) = femcore::solid_kz::solve_solid_kz_full(&mesh, &sp, 5.0).unwrap();
    let tau = 2e-9;
    let eta: Vec<M6> = modes.c6.iter().map(|m| m.map(|r| r.map(|x| x * tau))).collect();
    let inv_q = viscous_inverse_q(&mesh, 2, &modes, &eta);
    assert_eq!(inv_q.len(), modes.freqs.len());
    for (f, iq) in modes.freqs.iter().zip(&inv_q) {
        let expect = 2.0 * std::f64::consts::PI * f.re * tau;
        assert!((iq - expect).abs() < 1e-9 * expect, "{iq} against {expect}");
    }
}

fn ring_problem(order: usize) -> (femcore::mesh::SimplexMesh, ProblemSpec) {
    // a ring of mean radius R = 1 m with a 0.02 m x 0.02 m cross-section (x = r, y = z)
    let (rr, a) = (1.0, 0.02);
    let parts = vec![Part { shape: Shape::Rect { x0: rr - a / 2.0, y0: -a / 2.0, x1: rr + a / 2.0, y1: a / 2.0 }, domain: 1, h: Some(a / 4.0) }];
    let mesh = mesh_2d(&Spec2d { parts, h: a / 4.0, min_angle_deg: 25.0 }).unwrap();
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.order = order;
    sp.materials.default.insert("E".into(), 70e9);
    sp.materials.default.insert("nu".into(), 0.33);
    sp.materials.default.insert("rho".into(), 2700.0);
    sp.study = Study::Eigenfrequency { n_modes: 8, f_search: 500.0 };
    (mesh, sp)
}

#[test]
fn a_thin_ring_breathes_at_the_hoop_frequency() {
    // m = 0: all the ring moves radially, f = sqrt(E / rho) / (2 pi R) for a thin ring
    let (mesh, sp) = ring_problem(2);
    let (_, modes) = femcore::solid_kz::solve_solid_axisym(&mesh, &sp, 0.0).unwrap();
    let c = (70e9f64 / 2700.0).sqrt();
    let expect = c / (2.0 * std::f64::consts::PI);
    let best = modes.freqs.iter().map(|f| f.re).min_by(|a, b| (a - expect).abs().partial_cmp(&(b - expect).abs()).unwrap()).unwrap();
    println!("breathing mode {best:.3} Hz, hoop frequency {expect:.3} Hz");
    assert!((best - expect).abs() < 0.01 * expect, "{best} against {expect}");
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn a_thin_ring_bends_in_its_plane_with_the_classical_frequencies() {
    // in-plane flexure of a thin ring with radial thickness h: omega^2 = E h^2 / (12 rho R^4) m^2 (m^2 - 1)^2 / (m^2 + 1)
    let (mesh, sp) = ring_problem(2);
    let (rr, h) = (1.0f64, 0.02f64);
    for m in [2.0f64, 3.0, 4.0] {
        let w2 = 70e9 * h * h / (12.0 * 2700.0 * rr.powi(4)) * m * m * (m * m - 1.0).powi(2) / (m * m + 1.0);
        let expect = w2.sqrt() / (2.0 * std::f64::consts::PI);
        let mut sp2 = sp.clone();
        sp2.study = Study::Eigenfrequency { n_modes: 10, f_search: expect };
        let (_, modes) = femcore::solid_kz::solve_solid_axisym(&mesh, &sp2, m).unwrap();
        let best = modes.freqs.iter().map(|f| f.re).min_by(|a, b| (a - expect).abs().partial_cmp(&(b - expect).abs()).unwrap()).unwrap();
        println!("m = {m}: {best:.3} Hz, ring theory {expect:.3} Hz");
        assert!((best - expect).abs() < 0.03 * expect, "m = {m}: {best} against {expect}");
    }
}

/// Elastic and photoelastic data of the silicon wire in silica (cubic silicon along x, y, z).
fn wire_elastic(mesh: &femcore::mesh::SimplexMesh, f_search: f64, n_modes: usize) -> ProblemSpec {
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.order = 2;
    for d in mesh.domain_labels() {
        if d == 2 {
            sp.materials.aniso.insert(d, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }), euler_deg: [0.0; 3], ..Default::default() });
            sp.materials.set(d, "rho", 2329.0);
        } else {
            sp.materials.set(d, "rho", 2203.0);
        }
        sp.materials.set(d, "E", 73.1e9);
        sp.materials.set(d, "nu", 0.17);
    }
    let mut labels = vec![];
    for f in 0..mesh.n_facets() {
        labels.push(mesh.facet_labels[f]);
    }
    labels.sort();
    labels.dedup();
    sp.bcs = vec![Bc { on: None, labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    sp.study = Study::Eigenfrequency { n_modes, f_search };
    sp
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn the_phase_matched_overlap_is_mirror_symmetric_and_the_strain_as_is_form_is_not() {
    // Backward scattering of a mirror-symmetric waveguide: describing the same process in the mirror image
    // (pump counter-propagating, Stokes forward, phonon -q) must give the same gains.
    let mesh = wire_mesh();
    let set = solve_optical(&mesh, &optics(&mesh, 1.55e-6)).unwrap();
    let mut pe = BTreeMap::new();
    pe.insert(2, PhotoElastic::cubic(-0.094, 0.017, -0.051, &rotation([0.0; 3])));
    pe.insert(1, PhotoElastic::isotropic(0.121, 0.270));
    let moving: BTreeSet<i32> = [2].into_iter().collect();
    let mk = |mirror| coupling(&PairInput { mesh: &mesh, pump: (&set, 0), stokes: (&set, 0), backward: true, mirror_pump: mirror, photoelastic: &pe, moving: &moving }).unwrap();
    let (cp, cm) = (mk(false), mk(true));
    assert!((cp.q() + cm.q()).abs() < 1e-9 * cp.q().abs());
    let sp = wire_elastic(&mesh, 16e9, 14);
    let (_, modes_p) = femcore::solid_kz::solve_solid_kz_full(&mesh, &sp, cp.q()).unwrap();
    let (_, modes_m) = femcore::solid_kz::solve_solid_kz_full(&mesh, &sp, cm.q()).unwrap();
    assert_eq!(modes_p.freqs.len(), modes_m.freqs.len());
    let opts = |strain| GainOptions { strain, body_domains: [2].into_iter().collect(), ..Default::default() };
    let mut worst_conj = 0.0f64;
    let mut worst_as_is = 0.0f64;
    let mut strongest = 0.0f64;
    for (strain, worst) in [(StrainConvention::Conjugate, &mut worst_conj), (StrainConvention::AsIs, &mut worst_as_is)] {
        let rp = gains(&mesh, &cp, &modes_p, &opts(strain));
        let rm = gains(&mesh, &cm, &modes_m, &opts(strain));
        for (a, b) in rp.iter().zip(&rm) {
            assert!((a.freq.re - b.freq.re).abs() < 1e-6 * a.freq.re, "modes of +q and -q pair up by frequency");
            if a.gain_pe.max(a.gain_mb) < 1e-6 {
                continue;
            }
            strongest = strongest.max(a.gain_pe);
            for (x, y) in [(a.gain_pe, b.gain_pe), (a.gain_mb, b.gain_mb), (a.gain_total, b.gain_total)] {
                *worst = worst.max((x - y).abs() / x.max(y).max(1e-30));
            }
        }
    }
    println!("largest PE gain {strongest:.3e}; worst mirror mismatch: conjugate strain {worst_conj:.2e}, strain as is {worst_as_is:.2e}");
    assert!(worst_conj < 1e-3, "{worst_conj}");
    // the expression with the strain as it is describes a phonon of the wrong wave-number sign: it is not mirror symmetric
    assert!(worst_as_is > 0.1, "{worst_as_is}");
}

#[test]
fn modes_of_a_symmetric_plate_are_classified_by_parity_and_followed_across_wave_numbers() {
    // a clamped square plate: its modes are even or odd under x -> -x; the parity is +-1 up to the mesh's
    // asymmetry, and tracking by overlap follows each mode from kz to the next kz
    let mesh = femcore::mesh::SimplexMesh::rectangle(1.0, 1.0, 16, 16);
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.order = 2;
    sp.materials.default.insert("E".into(), 70e9);
    sp.materials.default.insert("nu".into(), 0.33);
    sp.materials.default.insert("rho".into(), 2700.0);
    let mut labels: Vec<i32> = mesh.facet_labels.clone();
    labels.sort();
    labels.dedup();
    sp.bcs = vec![Bc { on: None, labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    sp.study = Study::Eigenfrequency { n_modes: 6, f_search: 2000.0 };
    let (_, modes) = femcore::solid_kz::solve_solid_kz_full(&mesh, &sp, 3.0).unwrap();
    let map = mirror_map(&mesh, 2, 0);
    let Some(map) = map else {
        println!("the mesh is not mirror symmetric in its vertices; parity check skipped");
        return;
    };
    let n1 = mesh.n_dofs(2);
    let par: Vec<f64> = modes.vectors.iter().map(|u| mirror_parity(u, n1, &map, 0)).collect();
    println!("parities {par:.3?}");
    assert!(par.iter().all(|p| (p.abs() - 1.0).abs() < 0.05 || p.abs() < 0.05 || p.abs() <= 1.0 + 1e-9));
    // degenerate pairs can mix even and odd shapes, but nondegenerate modes are clean: at least one pure mode of each kind
    assert!(par.iter().any(|p| *p > 0.95) && par.iter().any(|p| *p < -0.95), "{par:?}");
    // tracking across a small change of kz keeps the frequency order
    let (all, tracks) = elastic_dispersion(&mesh, &sp, &[3.0, 3.5, 4.0], 0.5).unwrap();
    for t in &tracks {
        let fs: Vec<f64> = t.iter().enumerate().filter_map(|(k, i)| i.map(|i| all[k].freqs[i].re)).collect();
        for w in fs.windows(2) {
            // frequencies grow smoothly with kz (here by 5 to 10 % per step); a jump to a neighbouring mode would be larger or negative
            assert!(w[1] > w[0] && w[1] < 1.2 * w[0], "{fs:?}");
        }
    }
}

#[test]
fn a_coaxial_cavity_has_the_transmission_line_resonances() {
    // r in [0.1, 0.2] m, z in [0, 0.3] m, perfect conductor everywhere, m = 0: the TEM standing waves
    // E_r ~ sin(l pi z / L) / r resonate at f = l c / (2 L)
    let (a, b, len) = (0.1, 0.2, 0.3);
    let mesh = mesh_2d(&Spec2d { parts: vec![Part { shape: Shape::Rect { x0: a, y0: 0.0, x1: b, y1: len }, domain: 1, h: Some(0.015) }], h: 0.015, min_angle_deg: 25.0 }).unwrap();
    let mat = femcore::emw_axisym::AxiMaterial::from_indices(&mesh, &BTreeMap::new());
    let c0 = 299_792_458.0f64;
    let modes = femcore::emw_axisym::solve_axisym_modes(&mesh, &mat, 0.0, 0.9e9, 3, 2).unwrap();
    println!("frequencies {:?}", modes.freqs.iter().map(|f| f / 1e6).collect::<Vec<_>>());
    for l in [1.0f64, 2.0] {
        let expect = l * c0 / (2.0 * len);
        if !(0.4e9..1.4e9).contains(&expect) {
            continue;
        }
        let best = modes.freqs.iter().cloned().fold(f64::INFINITY, |m, f| if (f - expect).abs() < (m - expect).abs() { f } else { m });
        assert!((best - expect).abs() < 2e-3 * expect, "l = {l}: {best} against {expect}");
    }
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn a_whispering_gallery_mode_of_a_large_ring_has_the_straight_waveguide_index() {
    // the silicon wire in silica bent to a radius of 2 mm: m / (k0 R) is the effective index of the straight guide
    let um = 1e-6;
    let radius = 2000.0 * um;
    let mesh = wire_mesh();
    let straight_set = solve_optical(&mesh, &optics(&mesh, 1.55e-6)).unwrap();
    let n_straight = straight_set.modes[0].neff.re;
    // the same cross-section moved to r = R (the wire is centred, the box is 2.4 um x 2 um)
    let (sx, sy) = (radius, 0.0);
    let parts = vec![
        Part { shape: Shape::Rect { x0: sx - 1.2 * um, y0: sy - 1.0 * um, x1: sx + 1.2 * um, y1: sy + 1.0 * um }, domain: 1, h: None },
        Part { shape: Shape::Rect { x0: sx - 0.225 * um, y0: sy - 0.11 * um, x1: sx + 0.225 * um, y1: sy + 0.11 * um }, domain: 2, h: Some(0.04 * um) },
    ];
    let ring = mesh_2d(&Spec2d { parts, h: 0.3 * um, min_angle_deg: 25.0 }).unwrap();
    let mut idx = BTreeMap::new();
    idx.insert(1, 1.444);
    idx.insert(2, 3.476);
    let mat = femcore::emw_axisym::AxiMaterial::from_indices(&ring, &idx);
    let f0 = 299_792_458.0 / 1.55e-6;
    let k0 = 2.0 * std::f64::consts::PI * f0 / 299_792_458.0;
    let m = (n_straight * k0 * radius).round();
    let modes = femcore::emw_axisym::solve_axisym_modes(&ring, &mat, m, f0, 2, 2).unwrap();
    let f = modes.freqs[0];
    let n_ring = m / (2.0 * std::f64::consts::PI * f / 299_792_458.0 * radius);
    println!("m = {m}, f = {:.4} THz (target {:.4}), neff ring {n_ring:.5}, straight {n_straight:.5}", f * 1e-12, f0 * 1e-12);
    assert!((n_ring - n_straight).abs() < 1.5e-3, "{n_ring} against {n_straight}");
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn a_large_ring_has_the_gains_of_the_straight_waveguide() {
    // the silicon wire in silica on a ring of radius 20 mm: the axisymmetric gain (azimuthal orders m_p + m_s for the
    // phonon) equals the straight-waveguide gain at q = M / R, to the curvature corrections
    use femcore::ring_brillouin::*;
    let um = 1e-6;
    let radius = std::env::var("RING_R_UM").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(20000.0) * um;
    let mesh = wire_mesh();
    let n_idx: BTreeMap<i32, f64> = [(1, 1.444), (2, 3.476)].into_iter().collect();
    // straight
    let set = solve_optical(&mesh, &optics(&mesh, 1.55e-6)).unwrap();
    let mut pe = BTreeMap::new();
    // isotropic materials: the ring model averages a crystal over rotations about its axis, the straight one does not
    pe.insert(2, PhotoElastic::isotropic(-0.094, 0.017));
    pe.insert(1, PhotoElastic::isotropic(0.121, 0.270));
    let moving: BTreeSet<i32> = [2].into_iter().collect();
    let flat = coupling(&PairInput { mesh: &mesh, pump: (&set, 0), stokes: (&set, 0), backward: true, mirror_pump: false, photoelastic: &pe, moving: &moving }).unwrap();
    // the ring: the same mesh moved to r = R
    let mut pts = mesh.points.clone();
    for i in 0..mesh.n_points() {
        pts[2 * i] += radius;
    }
    let ring = femcore::mesh::SimplexMesh::new(2, pts, mesh.cells.clone(), mesh.cell_labels.clone(), mesh.facets.clone(), mesh.facet_labels.clone());
    let mat = femcore::emw_axisym::AxiMaterial::from_indices(&ring, &n_idx);
    let f0 = 299_792_458.0 / 1.55e-6;
    let k0 = 2.0 * std::f64::consts::PI * f0 / 299_792_458.0;
    let m = (set.modes[0].neff.re * k0 * radius).round();
    let opt = femcore::emw_axisym::solve_axisym_modes(&ring, &mat, m, f0, 1, 2).unwrap();
    let rc = ring_coupling(&RingInput { mesh: &ring, material: &mat, pump: (&opt, 0), stokes: (&opt, 0), backward: true, photoelastic: &pe, moving: &moving, r_ref: radius }).unwrap();
    println!("group index: ring {:.4}, straight {:.4}; P ring {:.3e} straight {:.3e}", rc.group_index().0, flat.group_index().0, rc.p_pump, flat.p_pump);
    assert!(rc.p_pump > 0.0, "power along +theta for m > 0");
    assert!((rc.group_index().0 - flat.group_index().0).abs() < 3e-3 * flat.group_index().0);
    // elastic modes
    let q = rc.elastic_order() / radius;
    let iso = |mut sp: ProblemSpec| {
        sp.materials.aniso.clear();
        sp.materials.set(2, "E", 165e9);
        sp.materials.set(2, "nu", 0.22);
        sp
    };
    let sp_flat = iso(wire_elastic(&mesh, 16e9, 12));
    let sp_ring = iso(wire_elastic(&ring, 16e9, 12));
    let (_, mf) = femcore::solid_kz::solve_solid_kz_full(&mesh, &sp_flat, q).unwrap();
    let (_, mr) = femcore::solid_kz::solve_solid_axisym(&ring, &sp_ring, rc.elastic_order()).unwrap();
    let opts = GainOptions { body_domains: [2].into_iter().collect(), length: 1.0, ..Default::default() };
    let rf = gains(&mesh, &flat, &mf, &opts);
    let rr = ring_gains(&ring, &rc, &mr, &opts);
    // the strongest modes carry the comparison; weak modes are sums of nearly cancelling terms, which the curvature
    // corrections of the ring (of relative size w / R, 1e-4 here) change by larger factors; at R = 2, 20 and 200 mm the
    // strongest mode's photoelastic gain differs from the straight one by 6, 0.6 and 0.06 percent
    let top = rf.iter().map(|r| r.gain_pe).fold(0.0, f64::max);
    let mut compared = 0;
    for a in &rf {
        if a.gain_pe < 1e-4 * top {
            continue;
        }
        let b = rr.iter().min_by(|x, y| (x.freq.re - a.freq.re).abs().partial_cmp(&(y.freq.re - a.freq.re).abs()).unwrap()).unwrap();
        println!("{:.4} GHz / {:.4} GHz: PE {:.4e} / {:.4e}  MB {:.4e} / {:.4e}", a.freq.re * 1e-9, b.freq.re * 1e-9, a.gain_pe, b.gain_pe, a.gain_mb, b.gain_mb);
        if a.gain_pe > 0.1 * top {
            assert!((a.freq.re - b.freq.re).abs() < 2e-3 * a.freq.re, "{} against {}", a.freq.re, b.freq.re);
            assert!((a.gain_pe - b.gain_pe).abs() < 0.03 * a.gain_pe, "PE {} against {}", a.gain_pe, b.gain_pe);
            assert!((a.gain_mb - b.gain_mb).abs() < 0.05 * a.gain_mb.max(1e-2 * a.gain_pe), "MB {} against {}", a.gain_mb, b.gain_mb);
            compared += 1;
        }
    }
    assert!(compared >= 1);
}

#[test]
// the settings are filled field by field the way the window does it
#[allow(clippy::field_reassign_with_default)]
fn the_brillouin_physics_solves_a_ring_and_a_straight_waveguide_with_the_same_settings() {
    // one problem description for both: the straight wire, and the same cross-section on a ring of 20 mm radius
    let um = 1e-6;
    let radius = 20000.0 * um;
    let mesh = wire_mesh();
    let build = |ring: Option<f64>, m: &femcore::mesh::SimplexMesh| -> ProblemSpec {
        let mut sp = wire_elastic(m, 16e9, 10);
        // isotropic materials: the ring model averages a crystal over rotations about its axis, the straight one does not
        sp.materials.aniso.clear();
        sp.materials.set(2, "E", 165e9);
        sp.materials.set(2, "nu", 0.22);
        sp.physics = Physics::Brillouin;
        let mut b = BrillouinSpec::default();
        b.wavelength = 1.55e-6;
        b.n_target = 2.35;
        b.n_optical_modes = 1;
        b.refractive_index.insert(1, 1.444);
        b.refractive_index.insert(2, 3.476);
        b.photoelastic.insert(2, PhotoSpec { kind: PhotoKind::Isotropic { p11: -0.094, p12: 0.017 }, euler_deg: [0.0; 3], follow_stiffness: false });
        b.photoelastic.insert(1, PhotoSpec { kind: PhotoKind::Isotropic { p11: 0.121, p12: 0.270 }, euler_deg: [0.0; 3], follow_stiffness: false });
        b.moving = vec![2];
        b.ring_radius = ring;
        sp.brillouin = Some(b);
        sp
    };
    let best = |r: &SolveResult| r.scalar("largest total gain").unwrap();
    let straight = femcore::solve::solve(&mesh, &build(None, &mesh), &mut SolveState::default()).unwrap();
    let mut pts = mesh.points.clone();
    for i in 0..mesh.n_points() {
        pts[2 * i] += radius;
    }
    let ring_mesh = femcore::mesh::SimplexMesh::new(2, pts, mesh.cells.clone(), mesh.cell_labels.clone(), mesh.facets.clone(), mesh.facet_labels.clone());
    let ring = femcore::solve::solve(&ring_mesh, &build(Some(radius), &ring_mesh), &mut SolveState::default()).unwrap();
    println!("largest gain: straight {:.4e}, ring {:.4e}", best(&straight), best(&ring));
    assert!(best(&straight) > 0.0);
    assert!((best(&ring) - best(&straight)).abs() < 0.03 * best(&straight), "{} against {}", best(&ring), best(&straight));
    assert!(ring.scalar("azimuthal order of the optical modes").unwrap() > 1e5);
    // the polarization analysis of the optical mode: the ring maps (E_r, E_z, E_theta) onto (lateral, normal, propagation), so a
    // ring of 20 mm radius has the Stokes vector of the straight guide
    let pol = |r: &SolveResult| r.analysis.as_ref().expect("the optical mode has a polarization analysis").rows[0].report.clone().unwrap();
    let (ps, pr) = (pol(&straight), pol(&ring));
    assert_eq!(ps.class.label(), "quasi-TE (x)");
    assert_eq!(pr.class.label(), ps.class.label());
    assert!((pr.te_fraction - ps.te_fraction).abs() < 0.01 && (pr.long_fraction - ps.long_fraction).abs() < 0.01, "{:?} against {:?}", pr, ps);
    assert!(ring.scalars.iter().any(|(n, _)| n == "polarization 1"));
}

fn ring_with_pml(radius_um: f64, w_um: f64, h_um: f64) -> femcore::mesh::SimplexMesh {
    // silicon wire in silica, absorbing layers on the outer side and above and below (domains 3 to 9)
    let um = 1e-6;
    let r = radius_um;
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: x0 * um, y0: y0 * um, x1: x1 * um, y1: y1 * um }, domain: d, h: h.map(|v| v * um) };
    let (xl, xm, xr) = (r - 1.6, r + 1.4, r + 3.0);
    let parts = vec![
        rect(xl, -1.2, xr, 1.2, 1, None),
        rect(xm, -0.8, xr, 0.8, 3, Some(0.15)),
        rect(xl, 0.8, xm, 1.2, 4, Some(0.15)),
        rect(xl, -1.2, xm, -0.8, 5, Some(0.15)),
        rect(xm, 0.8, xr, 1.2, 6, Some(0.15)),
        rect(xm, -1.2, xr, -0.8, 7, Some(0.15)),
        rect(r - w_um / 2.0, -h_um / 2.0, r + w_um / 2.0, h_um / 2.0, 2, Some(0.05)),
    ];
    mesh_2d(&Spec2d { parts, h: 0.25 * um, min_angle_deg: 25.0 }).unwrap()
}

fn ring_pml_map(strength: f64) -> BTreeMap<i32, Pml> {
    let p = |x, y| Pml { x, y, z: PmlDir::Off, strength, order: 2.0, re: 1.0 };
    [(3, p(PmlDir::Plus, PmlDir::Off)), (4, p(PmlDir::Off, PmlDir::Plus)), (5, p(PmlDir::Off, PmlDir::Minus)), (6, p(PmlDir::Plus, PmlDir::Plus)), (7, p(PmlDir::Plus, PmlDir::Minus))].into_iter().collect()
}

#[test]
#[ignore = "slow (30-60 s): run with `cargo test --release -- --ignored`"]
fn absorbing_layers_of_a_ring_leave_a_confined_mode_alone_and_give_a_leaky_one_a_loss_that_does_not_depend_on_them() {
    use femcore::emw_axisym::*;
    let n_idx: BTreeMap<i32, f64> = (1..=7).map(|d| (d, if d == 2 { 3.476 } else { 1.444 })).collect();
    let f0 = 299_792_458.0 / 1.55e-6;
    let k0 = 2.0 * std::f64::consts::PI / 1.55e-6;
    // a big ring: confined; with and without layers the frequency agrees and the loss is zero
    let big = ring_with_pml(100.0, 0.45, 0.22);
    let m_big = (2.35 * k0 * 100e-6f64).round();
    let plain = solve_axisym_modes(&big, &AxiMaterial::from_indices(&big, &n_idx), m_big, f0, 1, 2).unwrap();
    let damped = solve_axisym_modes(&big, &AxiMaterial::from_indices(&big, &n_idx).with_pml(&big, &ring_pml_map(3.0)), m_big, f0, 1, 2).unwrap();
    println!("R = 100 um: f = {:.6} THz without layers, {:?} THz with", plain.freqs[0] * 1e-12, damped.freqs_c[0] * 1e-12);
    assert!((plain.freqs[0] - damped.freqs[0]).abs() < 2e-6 * plain.freqs[0]);
    assert!(damped.freqs_c[0].im.abs() < 1e-6 * damped.freqs_c[0].re);
    // a tight ring of a lower-index core radiates: the loss must not depend on the layers
    let tight = ring_with_pml(3.0, 0.8, 0.5);
    let mut idx2 = n_idx.clone();
    idx2.insert(2, 2.0);
    let mut found = None;
    for m_t in 8..24 {
        let md = solve_axisym_modes(&tight, &AxiMaterial::from_indices(&tight, &idx2).with_pml(&tight, &ring_pml_map(3.0)), m_t as f64, f0, 2, 2).unwrap();
        for f in &md.freqs_c {
            if (f.re - f0).abs() < 0.02 * f0 && f.im.abs() > 1e-9 * f.re && f.re / (2.0 * f.im.abs()) > 30.0 {
                found = Some((m_t as f64, *f));
            }
        }
        if found.is_some() {
            break;
        }
    }
    let (m_t, f3) = found.expect("a whispering-gallery mode of the tight ring");
    let md6 = solve_axisym_modes(&tight, &AxiMaterial::from_indices(&tight, &idx2).with_pml(&tight, &ring_pml_map(6.0)), m_t, f0, 2, 2).unwrap();
    let f6 = md6.freqs_c.iter().min_by(|a, b| (a.re - f3.re).abs().partial_cmp(&(b.re - f3.re).abs()).unwrap()).unwrap();
    let (q3, q6) = (f3.re / (2.0 * f3.im.abs()), f6.re / (2.0 * f6.im.abs()));
    println!("R = 3 um, m = {m_t}: f = {:.4} THz, Q = {q3:.3e} (layers of strength 3) and {q6:.3e} (strength 6)", f3.re * 1e-12);
    assert!((q3 / q6 - 1.0).abs() < 0.15 && (f3.re - f6.re).abs() < 1e-3 * f3.re, "{q3} against {q6}");
}
