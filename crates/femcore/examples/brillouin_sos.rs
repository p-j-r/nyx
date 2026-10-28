//! Brillouin gain of the 300TE ring chip's rib waveguide (silicon rib on slab on r-plane sapphire), backward intra-mode scattering of the fundamental TE
//! mode (phonon near 18.35 GHz).
//!   cargo run --release -p femcore --example brillouin_sos -- [psi_deg] [f_search GHz] [n_modes]
use femcore::brillouin::*;
use femcore::material::*;
use femcore::mesher::*;
use femcore::problem::*;
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let psi = a.first().copied().unwrap_or(0.0).to_radians();
    let f_search = a.get(1).copied().unwrap_or(18.35) * 1e9;
    let n_modes = a.get(2).copied().unwrap_or(60.0) as usize;
    let envf = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(d);
    let um = 1e-6;
    let (w, htot, hslab) = (0.5 * um, 0.3 * um, 0.13 * um);
    let (xin, hin, lp, ha) = (2.0 * um, 2.0 * um, 1.0 * um, 1.2 * um);
    let x = xin + lp;
    let scale = envf("MESH_SCALE", 1.0);
    let (h_si, h_sap, h_pml, h_far) = (0.028 * scale * um, 0.05 * scale * um, 0.07 * scale * um, 0.12 * um);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: d, h };
    let ytop = htot + ha;
    let parts = vec![
        rect(-x, -hin - lp, x, ytop + lp, 4, None),
        // air PML pieces
        rect(-x, hslab, -xin, ytop, 41, Some(h_pml)),
        rect(xin, hslab, x, ytop, 42, Some(h_pml)),
        rect(-xin, ytop, xin, ytop + lp, 43, Some(h_pml)),
        rect(-x, ytop, -xin, ytop + lp, 44, Some(h_pml)),
        rect(xin, ytop, x, ytop + lp, 45, Some(h_pml)),
        // sapphire
        rect(-xin, -hin, xin, 0.0, 1, Some(h_sap)),
        rect(-x, -hin, -xin, 0.0, 11, Some(h_pml)),
        rect(xin, -hin, x, 0.0, 12, Some(h_pml)),
        rect(-xin, -hin - lp, xin, -hin, 13, Some(h_pml)),
        rect(-x, -hin - lp, -xin, -hin, 14, Some(h_pml)),
        rect(xin, -hin - lp, x, -hin, 15, Some(h_pml)),
        // silicon slab and rib
        rect(-xin, 0.0, xin, hslab, 2, Some(h_si)),
        rect(-x, 0.0, -xin, hslab, 21, Some(h_pml)),
        rect(xin, 0.0, x, hslab, 22, Some(h_pml)),
        rect(-w / 2.0, hslab, w / 2.0, htot, 3, Some(h_si)),
    ];
    let t0 = Instant::now();
    let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
    println!("mesh: {} cells, {} points ({:.1} s)", mesh.n_cells(), mesh.n_points(), t0.elapsed().as_secs_f64());
    let domains = mesh.domain_labels();
    let sap_all = [1, 11, 12, 13, 14, 15];
    let si_all = [2, 3, 21, 22];

    // ---- optics
    let lambda = 1.55e-6;
    let freq = 299_792_458.0 / lambda;
    let mut opt = ProblemSpec::new(Physics::Emw);
    opt.order = 2;
    for &d in &domains {
        let n: f64 = if sap_all.contains(&d) { 1.75 } else if si_all.contains(&d) { 3.476 } else { 1.0 };
        opt.materials.set(d, "eps_r", n * n);
        opt.materials.set(d, "mu_r", 1.0);
    }
    opt.study = Study::Modes { freq, n_target: 2.78, evanescent: false, n_modes: envf("OPT_MODES", 6.0) as usize };
    let t = Instant::now();
    let set = solve_optical(&mesh, &opt).unwrap();
    println!("optical modes ({:.1} s):", t.elapsed().as_secs_f64());
    for (j, m) in set.modes.iter().enumerate() {
        println!("  mode {j}: neff = {:.5}", m.neff.re);
    }
    let jp = (0..set.modes.len()).min_by(|&i, &j| (set.modes[i].neff.re - 2.7831).abs().partial_cmp(&(set.modes[j].neff.re - 2.7831).abs()).unwrap()).unwrap();
    println!("pump = Stokes = mode {jp} (neff {:.5})", set.modes[jp].neff.re);

    // ---- material tensors
    let b = 57.6f64.to_radians();
    let (s, cb) = (b.sin(), b.cos());
    let (sp, cp) = (psi.sin(), psi.cos());
    // crystal frame: x lateral, y normal, z propagation; p = c-axis projection, q = a axis
    let p = [sp, 0.0, cp];
    let q = [cp, 0.0, -sp];
    let y = [0.0, 1.0, 0.0];
    let cax = [s * p[0] + cb * y[0], s * p[1] + cb * y[1], s * p[2] + cb * y[2]];
    let aax = q;
    let cross = |u: [f64; 3], v: [f64; 3]| [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let a2 = cross(cax, aax);
    let cols = |c1: [f64; 3], c2: [f64; 3], c3: [f64; 3]| -> M3 { [[c1[0], c2[0], c3[0]], [c1[1], c2[1], c3[1]], [c1[2], c2[2], c3[2]]] };
    let mut r_sap = cols(aax, a2, cax);
    // diagnostic variants of the sapphire orientation (SAP_VARIANT): 1 tilt towards -z, 2 transposed rotation, 3 a axis along z
    match envf("SAP_VARIANT", 0.0) as i32 {
        1 => {
            let cf = [0.0, cb, -s];
            let a2f = cross(cf, aax);
            r_sap = cols(aax, a2f, cf);
        }
        2 => r_sap = [[r_sap[0][0], r_sap[1][0], r_sap[2][0]], [r_sap[0][1], r_sap[1][1], r_sap[2][1]], [r_sap[0][2], r_sap[1][2], r_sap[2][2]]],
        3 => {
            let a3 = [0.0, 0.0, 1.0];
            let c3 = [0.0, cb, s];
            let a2f = cross(c3, a3);
            r_sap = cols(a3, a2f, c3);
        }
        _ => {}
    }
    let r_si = cols(p, q, y);
    let eul_sap = euler_from_rotation(&r_sap);
    let eul_si = euler_from_rotation(&r_si);
    let sapphire = Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 };
    let silicon = Stiffness::Cubic { c11: envf("C11_SI", 165.7) * 1e9, c12: 63.9e9, c44: 79.6e9 };

    // ---- mechanics
    let mut mech = ProblemSpec::new(Physics::Solid);
    mech.order = 2;
    mech.inactive = domains.iter().copied().filter(|d| !sap_all.contains(d) && !si_all.contains(d)).collect();
    for &d in &domains {
        mech.materials.set(d, "E", 100e9);
        mech.materials.set(d, "nu", 0.25);
        mech.materials.set(d, "rho", 1.0);
    }
    for d in sap_all {
        mech.materials.aniso.insert(d, Anisotropy { stiffness: Some(sapphire.clone()), euler_deg: eul_sap, ..Default::default() });
        mech.materials.set(d, "rho", 3980.0);
    }
    for d in si_all {
        mech.materials.aniso.insert(d, Anisotropy { stiffness: Some(silicon.clone()), euler_deg: eul_si, ..Default::default() });
        mech.materials.set(d, "rho", 2329.0);
    }
    let mut labels = vec![];
    for f in 0..mesh.n_facets() {
        let (pa, pb) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
        let (mx, my) = (0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]));
        let outer = (mx.abs() - x).abs() < 1e-12 || (my - (-hin - lp)).abs() < 1e-12;
        if outer && my < hslab + 1e-12 {
            labels.push(mesh.facet_labels[f]);
        }
    }
    labels.sort();
    labels.dedup();
    mech.bcs = vec![Bc { labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let st = envf("PML_STRENGTH", 3.0);
    let pml = |xd, yd| Pml { x: xd, y: yd, z: PmlDir::Off, strength: st, order: 2.0, re: 1.0 };
    let mech_pml: BTreeMap<i32, Pml> = [
        (11, pml(PmlDir::Minus, PmlDir::Off)),
        (12, pml(PmlDir::Plus, PmlDir::Off)),
        (21, pml(PmlDir::Minus, PmlDir::Off)),
        (22, pml(PmlDir::Plus, PmlDir::Off)),
        (13, pml(PmlDir::Off, PmlDir::Minus)),
        (14, pml(PmlDir::Minus, PmlDir::Minus)),
        (15, pml(PmlDir::Plus, PmlDir::Minus)),
    ]
    .into_iter()
    .collect();
    mech.materials.pml = mech_pml.clone();
    mech.study = Study::Eigenfrequency { n_modes, f_search };

    // ---- coupling
    let beta = set.modes[jp].neff.re * 2.0 * std::f64::consts::PI / lambda;
    let kz = 2.0 * beta;
    println!("backward: q = {:.3} /um, searching {:.3} GHz", kz * 1e-6, f_search * 1e-9);
    let mut pe = BTreeMap::new();
    let p_si = PhotoElastic::cubic(-0.094, 0.017, -0.051, &r_si);
    let p_sap = PhotoElastic::trigonal(-0.23, 0.03, 0.02, 0.0, -0.2, -0.1, &r_sap);
    for d in [2, 3] {
        pe.insert(d, p_si.clone());
    }
    pe.insert(1, p_sap);
    let moving: BTreeSet<i32> = [2, 3].into_iter().collect();
    let cpl = coupling(&PairInput { mesh: &mesh, pump: (&set, jp), stokes: (&set, jp), backward: true, mirror_pump: false, photoelastic: &pe, moving: &moving }).unwrap();
    println!("optical: P' = {:.3e}, vg = {:.3e} m/s (ng = {:.3})", cpl.p_pump, cpl.vg_pump, cpl.group_index().0);
    let t = Instant::now();
    let (res, modes) = femcore::solid_kz::solve_solid_kz_full(&mesh, &mech, kz).unwrap();
    println!("{} elastic modes in {:.1} s ({})", modes.freqs.len(), t.elapsed().as_secs_f64(), res.log.join(" | "));
    let mk = |strain| GainOptions { body_domains: [2, 3].into_iter().collect(), pml_domains: mech_pml.keys().copied().collect(), strain, ..Default::default() };
    let rows = gains(&mesh, &cpl, &modes, &mk(StrainConvention::AsIs));
    let rows_c = gains(&mesh, &cpl, &modes, &mk(StrainConvention::Conjugate));
    let rows_n = gains(&mesh, &cpl, &modes, &GainOptions { mb_conjugate: false, ..mk(StrainConvention::AsIs) });
    let mut idx: Vec<usize> = (0..rows.len()).collect();
    idx.sort_by(|&a, &b| rows[b].gain_total.partial_cmp(&rows[a].gain_total).unwrap());
    println!("  f (GHz)      Q_eig      body  pml    gain_MB  gain_PE  total(as is)  total(conj)  |fPE|/|fMB| phase(deg)");
    for &i in idx.iter().take(10) {
        let (r, rc) = (&rows[i], &rows_c[i]);
        let ph = (r.f_pe.im.atan2(r.f_pe.re) - r.f_mb.im.atan2(r.f_mb.re)).to_degrees();
        let phc = (rc.f_pe.im.atan2(rc.f_pe.re) - rc.f_mb.im.atan2(rc.f_mb.re)).to_degrees();
        println!("  {:8.4}  {:10.3e}  {:5.2} {:5.2}  {:8.3}  {:8.3}  {:9.3}  {:9.3}   {:.2} / {:.2}  {:.0} / {:.0} (conj: PE {:.3})  MB without conj {:.3}", r.freq.re * 1e-9, r.q_eig, r.frac_body, r.frac_pml, r.gain_mb, r.gain_pe, r.gain_total, rc.gain_total, r.f_pe.norm() / r.f_mb.norm(), rc.f_pe.norm() / rc.f_mb.norm(), ph, phc, rc.gain_pe, rows_n[i].gain_mb);
    }
    if std::env::var_os("LIST").is_some() {
        let mut all: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].frac_body > 0.3).collect();
        all.sort_by(|&a, &b| rows[a].freq.re.partial_cmp(&rows[b].freq.re).unwrap());
        for i in all {
            println!("LIST {:8.3} GHz  Q {:10.3e}  body {:.2} pml {:.2}", rows[i].freq.re * 1e-9, rows[i].q_eig, rows[i].frac_body, rows[i].frac_pml);
        }
    }
}
