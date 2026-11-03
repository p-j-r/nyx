//! 800 nm wide silicon strip (300 nm thick, fully etched) on r-plane sapphire, air above: guided optical modes over 1300-1640 nm,
//! then the forward intermodal Brillouin gain of TE-launched pairs at 1550 nm (photoelastic PE and moving-boundary MB reported separately).
//! Wafer frame: x lateral, y propagation, z surface normal (the cross-section is x-z, meshed as x-y).
//!   cargo run --release -p femcore --example wg800 -- optics [out_dir]
//!   cargo run --release -p femcore --example wg800 -- acoustics [out_dir] [pump_mode_number] [f_search_GHz,...]
use femcore::sbs::*;
use femcore::chips::*;
use femcore::csparse::C;
use femcore::material::Anisotropy;
use femcore::mesh::SimplexMesh;
use femcore::geomdrive::MeshedGeom;
use femcore::problem::*;
use std::collections::BTreeSet;
use std::fmt::Write as _;

const UM: f64 = 1e-6;
const WIDTH: f64 = 0.8;
const HEIGHT: f64 = 0.3;
const N_CUT: f64 = N_SAPPHIRE_ORDINARY;

fn envf(k: &str, d: f64) -> f64 {
    std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(d)
}

/// the strip (`rib`) on sapphire in air with absorbing layers around: the chip of `chip_geometry` without the slab
fn model(psi: f64, scale: f64, n_modes: usize, f_search: f64) -> (MeshedGeom, ProblemSpec) {
    let c = ChipSpec { name: "WG800", htot: HEIGHT, etch: HEIGHT, width: WIDTH, xin: 2.0, n_eff: 2.6, polarisation: "TE" };
    let g = chip_mesh(&c, scale).unwrap();
    let mut spec = chip_problem(&c, psi, n_modes, f_search, &g);
    let b = spec.brillouin.as_mut().unwrap();
    b.moving = g.domains("rib");
    b.backward = false;
    (g, spec)
}

fn optical_spec(inner: &ProblemSpec, mesh: &SimplexMesh, lambda: f64, n_target: f64, n_modes: usize) -> ProblemSpec {
    let b = inner.brillouin.as_ref().unwrap();
    let mut opt = ProblemSpec::new(Physics::Emw);
    opt.order = inner.order;
    for d in mesh.domain_labels() {
        let n = b.refractive_index.get(&d).copied().unwrap_or(1.0);
        opt.materials.set(d, "eps_r", n * n);
        opt.materials.set(d, "mu_r", 1.0);
        if let Some(a) = inner.materials.aniso.get(&d) {
            if a.permittivity.is_some() {
                opt.materials.aniso.insert(d, Anisotropy { permittivity: a.permittivity.clone(), euler_deg: a.euler_deg, ..Default::default() });
            }
        }
    }
    opt.study = Study::Modes { freq: 299_792_458.0 / lambda, n_target, evanescent: false, n_modes };
    opt
}

struct Info {
    neff: f64,
    pol: &'static str,
    px: f64,
    py: f64,
    si_frac: f64,
    m: usize,
    v: usize,
    label: String,
}

/// polarisation (dominant lateral or normal field), fraction of |E|^2 in the silicon and the numbers of sign changes of the dominant field
/// along x at mid-height and along y through the silicon
fn classify(mesh: &SimplexMesh, si: &BTreeSet<i32>, set: &OpticalSet, j: usize) -> Info {
    let all = mode_energy_parts(mesh, set, j, None);
    let in_si = mode_energy_parts(mesh, set, j, Some(si));
    let tot: f64 = all.iter().sum();
    let (px, py) = (all[0] / tot, all[1] / tot);
    let te = px >= py;
    let comp = if te { 0 } else { 1 };
    let (w, h) = (WIDTH * UM, HEIGHT * UM);
    let nx = 161;
    let line_x: Vec<[f64; 2]> = (0..nx).map(|i| [-0.9 * w + 1.8 * w * i as f64 / (nx - 1) as f64, 0.5 * h]).collect();
    let ex = sample_mode_e(mesh, set, j, &line_x);
    let phase = {
        let mut best = (0.0, C::new(1.0, 0.0));
        for e in ex.iter().flatten() {
            if e[comp].norm() > best.0 {
                best = (e[comp].norm(), e[comp] / e[comp].norm());
            }
        }
        best.1.conj()
    };
    let signs = |v: Vec<f64>| -> usize {
        let m = v.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
        let mut last = 0.0;
        let mut n = 0;
        for x in v {
            if x.abs() > 0.1 * m {
                if last != 0.0 && x.signum() != last {
                    n += 1;
                }
                last = x.signum();
            }
        }
        n
    };
    let m = signs(ex.iter().map(|e| e.map_or(0.0, |e| (e[comp] * phase).re)).collect());
    // vertical: x slightly off-centre so that odd modes do not vanish
    let ny = 61;
    let line_y: Vec<[f64; 2]> = (0..ny).map(|i| [0.0 * w + 0.07 * w, 0.01 * h + 0.98 * h * i as f64 / (ny - 1) as f64]).collect();
    let ey = sample_mode_e(mesh, set, j, &line_y);
    let v = signs(ey.iter().map(|e| e.map_or(0.0, |e| (e[comp] * phase).re)).collect());
    let pol = if te { "TE" } else { "TM" };
    let label = if v > 0 { format!("{pol}{m}v{v}") } else { format!("{pol}{m}") };
    Info { neff: set.modes[j].neff.re, pol, px, py, si_frac: in_si.iter().sum::<f64>() / tot, m, v, label }
}

fn solve_modes(mesh: &SimplexMesh, inner: &ProblemSpec, lambda: f64) -> (OpticalSet, Vec<usize>) {
    let n_modes = envf("OPT_MODES", 24.0) as usize;
    let opt = optical_spec(inner, mesh, lambda, envf("N_TARGET", 3.4), n_modes);
    let set = solve_optical(mesh, &opt).unwrap();
    if std::env::var_os("DEBUG").is_some() {
        println!("raw: {}", set.modes.iter().map(|m| format!("{:.5}{:+.1e}i", m.neff.re, m.neff.im)).collect::<Vec<_>>().join(" "));
    }
    let mut idx: Vec<usize> = (0..set.modes.len()).filter(|&j| set.modes[j].neff.re > N_CUT + 0.002 && set.modes[j].neff.im.abs() < 1e-3).collect();
    idx.sort_by(|&a, &b| set.modes[b].neff.re.partial_cmp(&set.modes[a].neff.re).unwrap());
    // the Krylov solver can return the same mode several times: keep one of each effective index
    let mut kept: Vec<usize> = vec![];
    for j in idx {
        if kept.last().is_none_or(|&k| (set.modes[k].neff.re - set.modes[j].neff.re).abs() > 1e-5) {
            kept.push(j);
        }
    }
    let idx = kept;
    let lowest = set.modes.iter().map(|m| m.neff.re).fold(f64::MAX, f64::min);
    if lowest > N_CUT || !set.converged {
        eprintln!("warning at {:.0} nm: lowest neff found {lowest:.4}, converged {} - guided modes may be missing", lambda * 1e9, set.converged);
    }
    (set, idx)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let what = args.first().cloned().unwrap_or_else(|| "optics".into());
    let out = std::path::PathBuf::from(args.get(1).cloned().unwrap_or_else(|| ".".into()));
    std::fs::create_dir_all(&out).unwrap();
    let scale = envf("MESH_SCALE", 1.0);
    let psi = envf("PSI", 0.0);
    let f0 = envf("F_SEARCH", 8.0) * 1e9;
    let t0 = std::time::Instant::now();
    let (g, spec) = model(psi, scale, envf("N_ELASTIC", 40.0) as usize, f0);
    let si: BTreeSet<i32> = g.domains("rib").into_iter().collect();
    let mesh = g.mesh;
    println!("mesh: {} cells, {} points", mesh.n_cells(), mesh.n_points());
    let inner = spec.clone();
    if what == "optics" {
        let mut csv = String::from("lambda_nm,label,neff,pol,frac_Ex,frac_Ey,frac_in_Si,m,v\n");
        let (l0, l1, dl) = (envf("L0", 1300.0), envf("L1", 1640.0), envf("DL", 20.0));
        let mut lam = l0;
        while lam <= l1 + 1e-9 {
            let (set, idx) = solve_modes(&mesh, &inner, lam * 1e-9);
            let infos: Vec<Info> = idx.iter().map(|&j| classify(&mesh, &si, &set, j)).collect();
            print!("{lam:6.0} nm:");
            for i in &infos {
                print!("  {} {:.4}", i.label, i.neff);
                writeln!(csv, "{lam},{},{:.6},{},{:.4},{:.4},{:.4},{},{}", i.label, i.neff, i.pol, i.px, i.py, i.si_frac, i.m, i.v).unwrap();
            }
            println!();
            lam += dl;
        }
        std::fs::write(out.join("optical_modes.csv"), csv).unwrap();
        if std::env::var_os("NO_FIELDS").is_some() {
            return;
        }
        // field maps at the centre wavelength
        let lam = 1550e-9;
        let (set, idx) = solve_modes(&mesh, &inner, lam);
        let (nx, ny) = (96usize, 60usize);
        let (xa, xb, ya, yb) = (-1.2 * UM, 1.2 * UM, -0.7 * UM, 0.9 * UM);
        let pts: Vec<[f64; 2]> = (0..ny).flat_map(|jy| (0..nx).map(move |ix| [xa + (xb - xa) * ix as f64 / (nx - 1) as f64, yb - (yb - ya) * jy as f64 / (ny - 1) as f64])).collect();
        let mut fields = String::from("label,neff,nx,ny,xa_um,xb_um,ya_um,yb_um\n");
        for (k, &j) in idx.iter().enumerate() {
            let info = classify(&mesh, &si, &set, j);
            let e = sample_mode_e(&mesh, &set, j, &pts);
            let mut rows = String::new();
            for comp in 0..3 {
                for v in &e {
                    write!(rows, "{:.5e},", v.map_or(f64::NAN, |e| e[comp].norm())).unwrap();
                }
                rows.push('\n');
            }
            std::fs::write(out.join(format!("field_{k}.csv")), rows).unwrap();
            writeln!(fields, "{},{:.6},{nx},{ny},{},{},{},{}", info.label, info.neff, xa / UM, xb / UM, ya / UM, yb / UM).unwrap();
        }
        std::fs::write(out.join("fields_index.csv"), fields).unwrap();
    } else if what == "shapes" {
        // displacement fields of chosen elastic modes: shapes <out> <pump> <stokes> <f_GHz,...>
        let lambda = 1550e-9;
        let (set, idx) = solve_modes(&mesh, &inner, lambda);
        let pump_k = args.get(2).and_then(|s| s.parse::<usize>().ok()).unwrap_or(1) - 1;
        let stokes_k = args.get(3).and_then(|s| s.parse::<usize>().ok()).unwrap_or(2) - 1;
        let freqs: Vec<f64> = args.get(4).map(|s| s.split(',').filter_map(|x| x.parse::<f64>().ok()).collect()).unwrap_or_default();
        let b = inner.brillouin.as_ref().unwrap();
        let photoelastic: std::collections::BTreeMap<i32, PhotoElastic> = b.photoelastic.iter().map(|(d, p)| (*d, photo_from_spec(p, inner.materials.aniso.get(d).map(|a| a.euler_deg)))).collect();
        let moving: BTreeSet<i32> = si.clone();
        let cp = coupling(&PairInput { mesh: &mesh, pump: (&set, idx[pump_k]), stokes: (&set, idx[stokes_k]), backward: false, mirror_pump: false, photoelastic: &photoelastic, moving: &moving }).unwrap();
        let q = cp.q();
        let (nx, ny) = (96usize, 60usize);
        let (xa, xb, ya, yb) = (-1.2 * UM, 1.2 * UM, -0.7 * UM, 0.9 * UM);
        let pts: Vec<[f64; 2]> = (0..ny).flat_map(|jy| (0..nx).map(move |ix| [xa + (xb - xa) * ix as f64 / (nx - 1) as f64, yb - (yb - ya) * jy as f64 / (ny - 1) as f64])).collect();
        let tris: Vec<[[f64; 2]; 3]> = (0..mesh.n_cells()).map(|ci| { let v = mesh.cell(ci); let p = |i: usize| { let r = mesh.point(v[i] as usize); [r[0], r[1]] }; [p(0), p(1), p(2)] }).collect();
        let locate: Vec<Option<(usize, [f64; 3])>> = pts.iter().map(|&[x, y]| {
            for (ci, t) in tris.iter().enumerate() {
                let det = (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) - (t[2][0] - t[0][0]) * (t[1][1] - t[0][1]);
                let l1 = ((x - t[0][0]) * (t[2][1] - t[0][1]) - (t[2][0] - t[0][0]) * (y - t[0][1])) / det;
                let l2 = ((t[1][0] - t[0][0]) * (y - t[0][1]) - (x - t[0][0]) * (t[1][1] - t[0][1])) / det;
                let l0 = 1.0 - l1 - l2;
                if l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9 {
                    return Some((ci, [l0, l1, l2]));
                }
            }
            None
        }).collect();
        let n1 = mesh.n_dofs(2);
        println!("pair {} -> {}: q = {:.4} /um", infos_label(&set, &mesh, idx[pump_k]), infos_label(&set, &mesh, idx[stokes_k]), q * 1e-6);
        for &fg in &freqs {
            let mut mech = inner.clone();
            mech.physics = Physics::Solid;
            mech.brillouin = None;
            mech.study = Study::Eigenfrequency { n_modes: 10, f_search: fg * 1e9 };
            let (_r, modes) = femcore::sms_kz::solve_solid_kz_full(&mesh, &mech, q).unwrap();
            let j = (0..modes.freqs.len()).min_by(|&a, &b| (modes.freqs[a].re - fg * 1e9).abs().partial_cmp(&(modes.freqs[b].re - fg * 1e9).abs()).unwrap()).unwrap();
            let u = &modes.vectors[j];
            println!("  target {fg} GHz -> mode at {:.4} GHz", modes.freqs[j].re * 1e-9);
            // phase so that the largest component is real
            let mut big = (0.0, C::new(1.0, 0.0));
            for v in u.iter() {
                if v.norm() > big.0 {
                    big = (v.norm(), *v / v.norm());
                }
            }
            let ph = big.1.conj();
            let mut pd = Vec::new();
            let mut rows = String::new();
            for comp in 0..3 {
                for loc in &locate {
                    let val = match loc {
                        Some((ci, l)) => {
                            mesh.cell_dofs(*ci, 2, &mut pd);
                            let sh = femcore::fem::shape(l, 2).0;
                            let mut acc = C::new(0.0, 0.0);
                            for a in 0..sh.len() {
                                acc += u[comp * n1 + pd[a]] * sh[a];
                            }
                            (acc * ph).re / big.0
                        }
                        None => f64::NAN,
                    };
                    write!(rows, "{val:.5e},").unwrap();
                }
                rows.push('\n');
            }
            std::fs::write(out.join(format!("shape_p{}s{}_{:.2}GHz.csv", pump_k + 1, stokes_k + 1, modes.freqs[j].re * 1e-9)), rows).unwrap();
        }
    } else {
        let lambda = 1550e-9;
        let (set, idx) = solve_modes(&mesh, &inner, lambda);
        let infos: Vec<Info> = idx.iter().map(|&j| classify(&mesh, &si, &set, j)).collect();
        for (k, i) in infos.iter().enumerate() {
            println!("optical mode {}: {} neff {:.5}", k + 1, i.label, i.neff);
        }
        let pump_k = args.get(2).and_then(|s| s.parse::<usize>().ok()).unwrap_or(1) - 1;
        let searches: Vec<f64> = args.get(3).map(|s| s.split(',').filter_map(|x| x.parse::<f64>().ok()).collect()).unwrap_or_else(|| vec![f0 * 1e-9]);
        let b = inner.brillouin.as_ref().unwrap();
        let photoelastic: std::collections::BTreeMap<i32, PhotoElastic> = b.photoelastic.iter().map(|(d, p)| (*d, photo_from_spec(p, inner.materials.aniso.get(d).map(|a| a.euler_deg)))).collect();
        let moving: BTreeSet<i32> = si.clone();
        let pml_domains: BTreeSet<i32> = inner.materials.pml.keys().copied().collect();
        let mut report = String::new();
        for (ks, i_s) in infos.iter().enumerate() {
            if ks == pump_k || std::env::var("STOKES").ok().and_then(|v| v.parse::<usize>().ok()).is_some_and(|v| v != ks + 1) {
                continue;
            }
            let cp = coupling(&PairInput { mesh: &mesh, pump: (&set, idx[pump_k]), stokes: (&set, idx[ks]), backward: false, mirror_pump: false, photoelastic: &photoelastic, moving: &moving }).unwrap();
            let q = cp.q();
            println!("\n== pump {} -> Stokes {}: q = {:.4} /um (n_p - n_s = {:.4})", infos[pump_k].label, i_s.label, q * 1e-6, infos[pump_k].neff - i_s.neff);
            writeln!(report, "pair {} -> {} q {:.5e}", infos[pump_k].label, i_s.label, q).unwrap();
            for &fs in &searches {
                let mut mech = inner.clone();
                mech.physics = Physics::Solid;
                mech.brillouin = None;
                mech.study = Study::Eigenfrequency { n_modes: envf("N_ELASTIC", 40.0) as usize, f_search: fs * 1e9 };
                let t = std::time::Instant::now();
                let (_r, modes) = femcore::sms_kz::solve_solid_kz_full(&mesh, &mech, q).unwrap();
                let opts = GainOptions { body_domains: moving.clone(), pml_domains: pml_domains.clone(), strain: StrainConvention::Conjugate, q_cap: 1000.0, ..Default::default() };
                let rows = gains(&mesh, &cp, &modes, &opts);
                let mut order: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].frac_body > 0.15).collect();
                order.sort_by(|&a, &b| rows[b].gain_total.partial_cmp(&rows[a].gain_total).unwrap());
                println!("  search {fs} GHz ({:.0} s): {} modes, f from {:.3} to {:.3} GHz", t.elapsed().as_secs_f64(), rows.len(), rows.iter().map(|r| r.freq.re).fold(f64::MAX, f64::min) * 1e-9, rows.iter().map(|r| r.freq.re).fold(0.0, f64::max) * 1e-9);
                println!("  f (GHz)   Q         body/pml   MB        PE        total");
                for (n_shown, &i) in order.iter().enumerate() {
                    let r = &rows[i];
                    if n_shown >= 6 && r.gain_total < 1e-3 {
                        continue;
                    }
                    println!("  {:8.4}  {:9.2e} {:.2}/{:.2}  {:9.3e} {:9.3e} {:9.3e}", r.freq.re * 1e-9, r.q_eig, r.frac_body, r.frac_pml, r.gain_mb, r.gain_pe, r.gain_total);
                    writeln!(report, "{},{},{:.5},{:.4e},{:.4},{:.4},{:.5e},{:.5e},{:.5e}", infos[pump_k].label, i_s.label, r.freq.re * 1e-9, r.q_eig, r.frac_body, r.frac_pml, r.gain_mb, r.gain_pe, r.gain_total).unwrap();
                }
            }
        }
        std::fs::write(out.join(format!("acoustics_pump{}_stokes{}.txt", pump_k + 1, std::env::var("STOKES").unwrap_or_default())), report).unwrap();
    }
    println!("total {:.0} s", t0.elapsed().as_secs_f64());
}

fn infos_label(set: &OpticalSet, _mesh: &SimplexMesh, j: usize) -> String {
    format!("mode with neff {:.4}", set.modes[j].neff.re)
}
