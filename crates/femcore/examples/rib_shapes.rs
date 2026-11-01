//! (Mode shapes) Silicon rib waveguide (300 nm film, 130 nm slab) on r-plane sapphire: guided optical modes at one wavelength and the forward
//! intermodal Brillouin gains of the TE0 pump with every other guided mode (photoelastic PE and moving-boundary MB separately), for one rib
//! width. Driven by `run_sweep.sh`, which loops over widths and runs workers in parallel.
//!   WIDTH=0.8 cargo run --release -p femcore --example rib_sweep -- <out_dir>
//! Environment: WIDTH (um), LAMBDA (nm, 1550), OPT_MODES (12), N_ELASTIC (40), F_LIST (GHz, comma list), WORKER/NWORK (this process
//! takes the Stokes modes k with (k-1) % NWORK == WORKER), MESH_SCALE, OPTICS_ONLY.
use femcore::brillouin::*;
use femcore::chips::*;
use femcore::csparse::C;
use femcore::material::Anisotropy;
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use std::collections::BTreeSet;
use std::fmt::Write as _;

const UM: f64 = 1e-6;
const HTOT: f64 = 0.3;
const ETCH: f64 = 0.17;

fn envf(k: &str, d: f64) -> f64 {
    std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(d)
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

/// Effective index of the fundamental TE or TM mode of the unetched slab (air / silicon / sapphire n_o), by bisection.
fn slab_neff(lambda: f64, d: f64, te: bool) -> f64 {
    let (n1, ns, nc) = (N_SILICON, N_SAPPHIRE_ORDINARY, 1.0);
    let k0 = 2.0 * std::f64::consts::PI / lambda;
    let f = |n: f64| {
        let kx = k0 * (n1 * n1 - n * n).sqrt();
        let (gc, gs) = (k0 * (n * n - nc * nc).sqrt(), k0 * (n * n - ns * ns).sqrt());
        let (ac, as_) = if te { (gc / kx, gs / kx) } else { ((n1 / nc).powi(2) * gc / kx, (n1 / ns).powi(2) * gs / kx) };
        kx * d - ac.atan() - as_.atan()
    };
    let (mut lo, mut hi) = (ns + 1e-9, n1 - 1e-9);
    if f(lo) < 0.0 {
        return ns;
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

struct Info {
    neff: f64,
    pol: &'static str,
    px: f64,
    py: f64,
    si_frac: f64,
    tail: f64,
    label: String,
}

/// polarisation (dominant lateral or normal component), silicon fraction, lateral order (sign changes of the dominant field along the
/// rib) and the lateral tail: the largest |E| in the slab 0.4-1.4 um beside the rib relative to the peak (a bound mode decays there).
fn classify(mesh: &SimplexMesh, set: &OpticalSet, j: usize, w: f64) -> Info {
    let si: BTreeSet<i32> = [2, 3, 21, 22].into_iter().collect();
    let all = mode_energy_parts(mesh, set, j, None);
    let in_si = mode_energy_parts(mesh, set, j, Some(&si));
    let tot: f64 = all.iter().sum();
    let (px, py) = (all[0] / tot, all[1] / tot);
    let te = px >= py;
    let comp = if te { 0 } else { 1 };
    let hslab = (HTOT - ETCH) * UM;
    let ymid = hslab + 0.5 * (HTOT * UM - hslab);
    let nx = 201;
    let line_x: Vec<[f64; 2]> = (0..nx).map(|i| [-0.9 * w + 1.8 * w * i as f64 / (nx - 1) as f64, ymid]).collect();
    let ex = sample_mode_e(mesh, set, j, &line_x);
    let mut best = (0.0, C::new(1.0, 0.0));
    for e in ex.iter().flatten() {
        if e[comp].norm() > best.0 {
            best = (e[comp].norm(), e[comp] / e[comp].norm());
        }
    }
    let phase = best.1.conj();
    let m = {
        let v: Vec<f64> = ex.iter().map(|e| e.map_or(0.0, |e| (e[comp] * phase).re)).collect();
        let mx = v.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
        let (mut last, mut n) = (0.0, 0);
        for x in v {
            if x.abs() > 0.1 * mx {
                if last != 0.0 && x.signum() != last {
                    n += 1;
                }
                last = x.signum();
            }
        }
        n
    };
    let tail_pts: Vec<[f64; 2]> = (0..40).flat_map(|i| {
        let d = (0.4 + i as f64 * 1.0 / 39.0) * UM;
        [[w / 2.0 + d, 0.5 * hslab], [-w / 2.0 - d, 0.5 * hslab]]
    }).collect();
    let tail = sample_mode_e(mesh, set, j, &tail_pts).iter().flatten().map(|e| (e[0].norm_sqr() + e[1].norm_sqr() + e[2].norm_sqr()).sqrt()).fold(0.0f64, f64::max);
    let peak = {
        let pk: Vec<[f64; 2]> = (0..41).map(|i| [-0.5 * w + w * i as f64 / 40.0, ymid]).collect();
        sample_mode_e(mesh, set, j, &pk).iter().flatten().map(|e| (e[0].norm_sqr() + e[1].norm_sqr() + e[2].norm_sqr()).sqrt()).fold(1e-30, f64::max)
    };
    let pol = if te { "TE" } else { "TM" };
    Info { neff: set.modes[j].neff.re, pol, px, py, si_frac: in_si.iter().sum::<f64>() / tot, tail: tail / peak, label: format!("{pol}{m}") }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let out = std::path::PathBuf::from(args.first().cloned().unwrap_or_else(|| ".".into()));
    std::fs::create_dir_all(&out).unwrap();
    let w_um = envf("WIDTH", 0.8);
    let lambda = envf("LAMBDA", 1550.0) * 1e-9;
    let n_el = envf("N_ELASTIC", 40.0) as usize;
    let (worker, nwork) = (envf("WORKER", 0.0) as usize, envf("NWORK", 1.0) as usize);
    let t0 = std::time::Instant::now();
    let chip = ChipSpec { name: "RIB", htot: HTOT, etch: ETCH, width: w_um, xin: 1.75 + 0.5 * w_um, n_eff: 2.7, polarisation: "TE" };
    let (mesh, spec) = chip_model(&chip, envf("PSI", 0.0), envf("MESH_SCALE", 1.0), n_el, 8e9).unwrap();
    let inner = spec.to_model_frame(2);
    println!("width {w_um} um: mesh {} cells", mesh.n_cells());

    // ---- optical modes
    let n_modes = envf("OPT_MODES", 12.0) as usize;
    let set = solve_optical(&mesh, &optical_spec(&inner, &mesh, lambda, 3.4, n_modes)).unwrap();
    let mut idx: Vec<usize> = (0..set.modes.len()).filter(|&j| set.modes[j].neff.im.abs() < 1e-3).collect();
    idx.sort_by(|&a, &b| set.modes[b].neff.re.partial_cmp(&set.modes[a].neff.re).unwrap());
    let mut kept: Vec<usize> = vec![];
    for j in idx {
        if kept.last().is_none_or(|&k| (set.modes[k].neff.re - set.modes[j].neff.re).abs() > 1e-5) {
            kept.push(j);
        }
    }
    let (n_te, n_tm) = (slab_neff(lambda, (HTOT - ETCH) * UM, true), slab_neff(lambda, (HTOT - ETCH) * UM, false));
    let cut = n_te.max(N_SAPPHIRE_ORDINARY) + 0.005;
    let lowest = kept.last().map_or(0.0, |&j| set.modes[j].neff.re);
    println!("slab modes: TE0 {n_te:.4}, TM0 {n_tm:.4}; bound means neff above the slab TE0 index, {cut:.4}; lowest mode found {lowest:.4} (converged {})", set.converged);
    let mut infos: Vec<(usize, Info)> = vec![];
    let mut txt = String::from("# k  label  neff  pol  frac_Ex frac_Ey frac_in_Si tail  status\n");
    let mut k = 0;
    for &j in &kept {
        let info = classify(&mesh, &set, j, w_um * UM);
        let bound = info.neff > cut;
        let status = if bound { "bound" } else { "below the slab TE0 index: couples to slab modes, not a bound mode" };
        if bound {
            k += 1;
        }
        writeln!(txt, "{} {} {:.5} {} {:.3} {:.3} {:.3} {:.3} {status}", if bound { k } else { 0 }, info.label, info.neff, info.pol, info.px, info.py, info.si_frac, info.tail).unwrap();
        if bound {
            infos.push((j, info));
        }
    }
    if lowest > cut {
        txt.push_str("# WARNING: every mode found is above the cutoff, guided modes may be missing\n");
    }
    print!("{txt}");

    // ---- displacement shapes of chosen elastic modes: REQ = "k:f_GHz,k:f_GHz" (k = number of the Stokes mode)
    let b = inner.brillouin.as_ref().unwrap();
    let photoelastic: std::collections::BTreeMap<i32, PhotoElastic> = b.photoelastic.iter().map(|(d, p)| (*d, photo_from_spec(p, inner.materials.aniso.get(d).map(|a| a.euler_deg)))).collect();
    let moving: BTreeSet<i32> = [2, 3].into_iter().collect();
    let reqs: Vec<(usize, f64)> = std::env::var("REQ").unwrap_or_default().split(',').filter_map(|x| { let mut p = x.split(':'); Some((p.next()?.parse().ok()?, p.next()?.parse().ok()?)) }).collect();
    let (nx, ny) = (140usize, 64usize);
    let half = (0.5 * w_um + 1.0) * UM;
    let (xa, xb, ya, yb) = (-half, half, -0.8 * UM, 0.8 * UM);
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
    for (k, fg) in reqs {
        let cp = coupling(&PairInput { mesh: &mesh, pump: (&set, infos[0].0), stokes: (&set, infos[k - 1].0), backward: false, mirror_pump: false, photoelastic: &photoelastic, moving: &moving }).unwrap();
        let q = cp.q();
        let mut mech = inner.clone();
        mech.physics = Physics::Solid;
        mech.brillouin = None;
        mech.study = Study::Eigenfrequency { n_modes: 12, f_search: fg * 1e9 };
        let (_r, modes) = femcore::solid_kz::solve_solid_kz_full(&mesh, &mech, q).unwrap();
        let j = (0..modes.freqs.len()).min_by(|&a, &b| (modes.freqs[a].re - fg * 1e9).abs().partial_cmp(&(modes.freqs[b].re - fg * 1e9).abs()).unwrap()).unwrap();
        let u = &modes.vectors[j];
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
        let name = format!("shape_w{w_um}_{}to{}_{:.3}GHz.csv", infos[0].1.label, infos[k - 1].1.label, modes.freqs[j].re * 1e-9);
        println!("{name}: q {:.4} /um, Q {:.1}", q * 1e-6, modes.freqs[j].re / (2.0 * modes.freqs[j].im.abs()));
        std::fs::write(out.join(name), rows).unwrap();
    }
    println!("total {:.0} s", t0.elapsed().as_secs_f64());
}
