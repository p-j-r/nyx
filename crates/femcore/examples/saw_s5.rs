//! Stage S5 of the topological SAW study, the FE part: disorder in the 8-ridge SSH chains (lambda = 3 um, q0 = 2 pi / 3 um, 300 nm Si ridges,
//! psi = 0, mesh and PML as in S3). Disorder breaks the mirror symmetry, so every run is a full-domain solve.
//!  * S5a (`SINGLE=1`): one ridge at q0 for the widths in WIDTHS (um) -> `s5a_single.csv` (f, Q, dofs, time).
//!  * S5b: every row of `s5_runs.csv` (written by `design_s5.py`): per-ridge widths w1..w8 and edge-to-edge gaps g1..g7 (um); for each run all 8
//!    ridge-bound modes with Q, per-ridge energy weights, and the complex projection amplitudes on the nominal single-ridge mode -> `s5_fe_modes.csv`, `s5_fe_proj.csv`.
//!   cargo run --release -p femcore --example saw_s5 -- [out_dir]
//! Environment: SINGLE (run S5a only), WIDTHS (default 2.7,2.85,3.0,3.15,3.3), ONLY (comma list of run names), N_MODES (70), V_FILM (5790), SHIFT (0.98), SIDE (5), DEPTH (5).
mod saw_array;
mod saw_common;
use femcore::csparse::C;
use femcore::mesh::SimplexMesh;
use femcore::mesher::*;
use femcore::problem::*;
use saw_array::*;
use saw_common::*;
use std::fmt::Write as _;

/// Chain with ridge widths `w` and edge-to-edge gaps `g` (m), centred on the middle of its outer edges. Returns the mesh, ridge centres, the outer edges and the Array record.
fn build_chain(lam: f64, w: &[f64], g: &[f64], scale: f64, side: f64, depth_lam: f64) -> (SimplexMesh, Array, Vec<f64>, (f64, f64)) {
    let n = w.len();
    let mut xc = vec![0.5 * w[0]];
    for i in 1..n {
        xc.push(xc[i - 1] + 0.5 * w[i - 1] + g[i - 1] + 0.5 * w[i]);
    }
    let (l, r) = (xc[0] - 0.5 * w[0], xc[n - 1] + 0.5 * w[n - 1]);
    let mid = 0.5 * (l + r);
    let xc: Vec<f64> = xc.iter().map(|x| x - mid).collect();
    let (xl, xr) = (l - mid, r - mid);
    let xe = xl.abs().max(xr.abs());
    let (xin_l, xin_r) = (xl - side * lam, xr + side * lam);
    let (pt, d) = (lam, depth_lam * lam);
    let (xo_l, xo_r) = (xin_l - pt, xin_r + pt);
    let (h_far, h_surf, h_si, h_edge) = (lam / 8.0 * scale, lam / 16.0 * scale, 75e-9 * scale, 75e-9 * scale);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: dom, h };
    let mut parts = vec![rect(xin_l, -d, xin_r, 0.0, 1, None), rect(xin_l, -0.6 * lam, xin_r, 0.0, 1, Some(h_surf))];
    for i in 0..n {
        parts.push(rect(xc[i] - 0.5 * w[i] - 0.25 * UM, -0.25 * UM, xc[i] + 0.5 * w[i] + 0.25 * UM, 0.0, 1, Some(h_edge)));
    }
    parts.push(rect(xin_r, -d, xo_r, 0.0, 5, None));
    parts.push(rect(xin_l, -d - pt, xin_r, -d, 6, None));
    parts.push(rect(xin_r, -d - pt, xo_r, -d, 9, None));
    parts.push(rect(xo_l, -d, xin_l, 0.0, 4, None));
    parts.push(rect(xo_l, -d - pt, xin_l, -d, 8, None));
    for i in 0..n {
        parts.push(rect(xc[i] - 0.5 * w[i], 0.0, xc[i] + 0.5 * w[i], T_SI, 2, Some(h_si)));
    }
    let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
    // xo is only used for the outer-boundary label search: both outer sides are passed through `outer`
    (mesh, Array { lam, xc: xc.clone(), xe, xo: xo_r, d, pt }, vec![xo_l, xo_r], (xl, xr))
}

fn chain_problem(mesh: &SimplexMesh, arr: &Array, xo: &[f64], psi: f64) -> ProblemSpec {
    let outer = outer_labels(mesh, &|p| (p[0] - xo[0]).abs() < 1e-12 || (p[0] - xo[1]).abs() < 1e-12 || (p[1] + arr.d + arr.pt).abs() < 1e-12);
    let mut spec = base_spec(psi, false);
    spec.bcs = vec![Bc { on: None, labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    spec
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "runs/saw_topo".into()));
    let lam = 3.0 * UM;
    let (scale, nm, v_film, shift, side, depth) = (envf("MESH_SCALE", 1.0), envf("N_MODES", 70.0) as usize, envf("V_FILM", 5790.0), envf("SHIFT", 0.98), envf("SIDE", 5.0), envf("DEPTH", 5.0));
    let psi = envf("PSI", 0.0);
    let q0 = 2.0 * std::f64::consts::PI / lam;
    let mut log = String::new();
    macro_rules! say {
        ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); log.push_str(&s); log.push('\n'); }};
    }

    // ---- S5a: single ridge, several widths
    let widths: Vec<f64> = std::env::var("WIDTHS").ok().map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![2.7, 2.85, 3.0, 3.15, 3.3]);
    let mut csv_a = String::from("width_um,f_GHz,f_im_GHz,Q,si_frac,pml_frac,dofs,cells,solve_s\n");
    let mut f_nominal = 0.0;
    let do_single = std::env::var_os("SINGLE").is_some() || !out.join("s5a_single.csv").exists();
    if do_single {
        for &wu in &widths {
            let (mesh, arr, xo, _) = build_chain(lam, &[wu * UM], &[], scale, side, depth);
            let spec = chain_problem(&mesh, &arr, &xo, psi);
            let f_search = 0.985 * v_film * q0 / (2.0 * std::f64::consts::PI);
            let sol = solve_q(&mesh, &spec, q0, 40, f_search).unwrap();
            let (cx, labels) = (cell_x(&mesh), mesh.cell_labels.clone());
            let tot = elastic_energy_where(&mesh, &sol.modes, 2, &|_| true);
            let core = elastic_energy_where(&mesh, &sol.modes, 2, &|c| !PML_DOMAINS.contains(&labels[c]) && cx[c].abs() < 0.5 * wu * UM + 0.5 * lam);
            let pmlf = elastic_energy_where(&mesh, &sol.modes, 2, &|c| PML_DOMAINS.contains(&labels[c]));
            let sif = elastic_energy_where(&mesh, &sol.modes, 2, &|c| labels[c] == 2);
            let mut o: Vec<usize> = (0..sol.modes.freqs.len()).collect();
            o.sort_by(|&a, &b| sol.modes.freqs[a].re.partial_cmp(&sol.modes.freqs[b].re).unwrap());
            let j = *o.iter().find(|&&j| core[j] / tot[j] >= CORE_MIN && pmlf[j] / tot[j] <= PML_MAX).expect("ridge mode");
            let f = sol.modes.freqs[j];
            say!("single ridge w = {wu} um: f = {:.6} GHz, Q {:.2e}, Si {:.3}, {} dofs, {:.1} s", f.re * 1e-9, q_factor(f), sif[j] / tot[j], sol.dofs, sol.seconds);
            writeln!(csv_a, "{wu},{:.8},{:.4e},{:.4e},{:.4},{:.3e},{},{},{:.2}", f.re * 1e-9, f.im * 1e-9, q_factor(f), sif[j] / tot[j], pmlf[j] / tot[j], sol.dofs, mesh.n_cells(), sol.seconds).unwrap();
            if (wu - 3.0).abs() < 1e-9 {
                f_nominal = f.re;
            }
            std::fs::write(out.join("s5a_single.csv"), &csv_a).unwrap();
        }
        std::fs::write(out.join("s5a_log.txt"), &log).unwrap();
        if std::env::var_os("SINGLE").is_some() {
            return;
        }
    }

    // ---- reference single-ridge mode (nominal width) for the projections
    let (mesh1, arr1, xo1, _) = build_chain(lam, &[lam], &[], scale, side, depth);
    let spec1 = chain_problem(&mesh1, &arr1, &xo1, psi);
    let sol1 = solve_q(&mesh1, &spec1, q0, 40, 0.985 * v_film * q0 / (2.0 * std::f64::consts::PI)).unwrap();
    let (cx1, labels1) = (cell_x(&mesh1), mesh1.cell_labels.clone());
    let tot1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|_| true);
    let core1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|c| !PML_DOMAINS.contains(&labels1[c]) && cx1[c].abs() < lam);
    let pml1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|c| PML_DOMAINS.contains(&labels1[c]));
    let mut o1: Vec<usize> = (0..sol1.modes.freqs.len()).collect();
    o1.sort_by(|&a, &b| sol1.modes.freqs[a].re.partial_cmp(&sol1.modes.freqs[b].re).unwrap());
    let j1 = *o1.iter().find(|&&j| core1[j] / tot1[j] >= CORE_MIN && pml1[j] / tot1[j] <= PML_MAX).expect("single ridge mode");
    let f_single = sol1.modes.freqs[j1].re;
    if f_nominal == 0.0 {
        f_nominal = f_single;
    }
    say!("reference single ridge: f = {:.6} GHz", f_single * 1e-9);
    let phi: Vec<[C; 3]> = {
        let v = vertex_vectors(&mesh1, &sol1.modes, j1);
        let w = vertex_weights(&mesh1);
        let norm = v.iter().zip(&w).map(|(a, w)| w * (a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr())).sum::<f64>().sqrt();
        v.into_iter().map(|a| [a[0] / norm, a[1] / norm, a[2] / norm]).collect()
    };
    let loc1 = Locator::new(&mesh1);

    // ---- S5b: the chain runs
    let only: Option<Vec<String>> = std::env::var("ONLY").ok().map(|s| s.split(',').map(|x| x.to_string()).collect());
    let runs = std::fs::read_to_string(out.join("s5_runs.csv")).expect("s5_runs.csv (run design_s5.py)");
    let mut modes_csv = String::from("run,chain,kind,idx,rank,f_GHz,f_im_GHz,df_from_single_MHz,Q,si_frac,pml_frac,core_frac");
    let mut proj_csv = String::from("run,rank,f_GHz");
    for m in 1..=8 {
        write!(modes_csv, ",w{m}").unwrap();
        write!(proj_csv, ",cre{m},cim{m}").unwrap();
    }
    modes_csv.push('\n');
    proj_csv.push('\n');
    let mut timing = String::from("run,dofs,cells,attempt,solve_s\n");
    for line in runs.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        let (name, chain, kind, idx) = (f[0], f[1], f[2], f[3]);
        if let Some(o) = &only {
            if !o.iter().any(|x| x == name) {
                continue;
            }
        }
        let g: Vec<f64> = f[4..11].iter().map(|x| x.parse::<f64>().unwrap() * UM).collect();
        let w: Vec<f64> = f[11..19].iter().map(|x| x.parse::<f64>().unwrap() * UM).collect();
        let n = w.len();
        let (mesh, arr, xo, (xl, xr)) = build_chain(lam, &w, &g, scale, side, depth);
        let spec = chain_problem(&mesh, &arr, &xo, psi);
        let (cx, labels, wv) = (cell_x(&mesh), mesh.cell_labels.clone(), vertex_weights(&mesh));
        let bounds: Vec<(f64, f64)> = (0..n)
            .map(|j| (if j == 0 { f64::NEG_INFINITY } else { 0.5 * (arr.xc[j - 1] + 0.5 * w[j - 1] + arr.xc[j] - 0.5 * w[j]) }, if j + 1 == n { f64::INFINITY } else { 0.5 * (arr.xc[j] + 0.5 * w[j] + arr.xc[j + 1] - 0.5 * w[j + 1]) }))
            .collect();
        let proj_ref: Vec<(Vec<usize>, Vec<[C; 3]>, f64)> = (0..n)
            .map(|j| {
                let (lo, hi) = bounds[j];
                let idx: Vec<usize> = (0..mesh.n_points()).filter(|&i| mesh.point(i)[0] >= lo && mesh.point(i)[0] < hi).collect();
                let vals: Vec<[C; 3]> = idx.iter().map(|&i| loc1.eval(&phi, mesh.point(i)[0] - arr.xc[j], mesh.point(i)[1])).collect();
                let nphi: f64 = idx.iter().zip(&vals).map(|(&i, a)| wv[i] * (a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr())).sum();
                (idx, vals, nphi)
            })
            .collect();
        say!("\n== {name}: gaps {:?} um, width deviations {:?} nm; {} cells", g.iter().map(|x| (x / UM * 1000.0).round() / 1000.0).collect::<Vec<_>>(), w.iter().map(|x| ((x - lam) * 1e9).round()).collect::<Vec<_>>(), mesh.n_cells());
        let mut done = false;
        for (attempt, (sh, nmodes)) in [(shift, nm), (0.94, nm + 30), (1.0, nm + 30)].into_iter().enumerate() {
            let f_search = sh * v_film * q0 / (2.0 * std::f64::consts::PI);
            let sol = match solve_q(&mesh, &spec, q0, nmodes, f_search) {
                Ok(s) => s,
                Err(e) => {
                    say!("  solve failed: {e}");
                    continue;
                }
            };
            say!("  attempt {}: {} dofs, {:.1} s, {} modes", attempt + 1, sol.dofs, sol.seconds, sol.modes.freqs.len());
            writeln!(timing, "{name},{},{},{},{:.2}", sol.dofs, mesh.n_cells(), attempt + 1, sol.seconds).unwrap();
            let m = &sol.modes;
            let tot = elastic_energy_where(&mesh, m, 2, &|_| true);
            let pmlf = elastic_energy_where(&mesh, m, 2, &|c| PML_DOMAINS.contains(&labels[c]));
            let core = elastic_energy_where(&mesh, m, 2, &|c| !PML_DOMAINS.contains(&labels[c]) && cx[c] > xl - 0.5 * lam && cx[c] < xr + 0.5 * lam);
            let mut ord: Vec<usize> = (0..m.freqs.len()).collect();
            ord.sort_by(|&a, &b| m.freqs[a].re.partial_cmp(&m.freqs[b].re).unwrap());
            let cands: Vec<usize> = ord.iter().copied().filter(|&j| core[j] / tot[j] >= CORE_MIN && pmlf[j] / tot[j] <= PML_MAX).collect();
            if cands.len() < n {
                say!("    only {} ridge-bound modes, retrying", cands.len());
                continue;
            }
            let sel: Vec<usize> = cands.into_iter().take(n).collect();
            let sif = elastic_energy_where(&mesh, m, 2, &|c| labels[c] == 2);
            let e_reg: Vec<Vec<f64>> = (0..n).map(|r| elastic_energy_where(&mesh, m, 2, &|c| cx[c] >= bounds[r].0 && cx[c] < bounds[r].1)).collect();
            say!("  rank  f - f_s (MHz)   Q        end-ridge energy weight (ridges 1+8)");
            for (rank, &j) in sel.iter().enumerate() {
                let fz = m.freqs[j];
                let u = vertex_vectors(&mesh, m, j);
                let u2: f64 = (0..mesh.n_points()).map(|i| wv[i] * (u[i][0].norm_sqr() + u[i][1].norm_sqr() + u[i][2].norm_sqr())).sum();
                let wts: Vec<f64> = (0..n).map(|r| e_reg[r][j] / tot[j]).collect();
                say!("  {:>3}   {:+9.4}   {:.1e}   {:.3}", rank + 1, (fz.re - f_single) * 1e-6, q_factor(fz), wts[0] + wts[n - 1]);
                write!(modes_csv, "{name},{chain},{kind},{idx},{},{:.9},{:.4e},{:.6},{:.3e},{:.4},{:.3e},{:.4}", rank + 1, fz.re * 1e-9, fz.im * 1e-9, (fz.re - f_single) * 1e-6, q_factor(fz), sif[j] / tot[j], pmlf[j] / tot[j], core[j] / tot[j]).unwrap();
                for x in &wts {
                    write!(modes_csv, ",{:.5}", x).unwrap();
                }
                modes_csv.push('\n');
                write!(proj_csv, "{name},{},{:.9}", rank + 1, fz.re * 1e-9).unwrap();
                for r in 0..n {
                    let (idx, vals, nphi) = &proj_ref[r];
                    let mut num = C::new(0.0, 0.0);
                    for (&i, a) in idx.iter().zip(vals) {
                        num += (a[0].conj() * u[i][0] + a[1].conj() * u[i][1] + a[2].conj() * u[i][2]) * wv[i];
                    }
                    let c = num / (nphi.sqrt() * u2.sqrt());
                    write!(proj_csv, ",{:.6e},{:.6e}", c.re, c.im).unwrap();
                }
                proj_csv.push('\n');
            }
            done = true;
            break;
        }
        if !done {
            say!("  ridge-bound modes not found for {name}");
        }
        std::fs::write(out.join("s5_fe_modes.csv"), &modes_csv).unwrap();
        std::fs::write(out.join("s5_fe_proj.csv"), &proj_csv).unwrap();
        std::fs::write(out.join("s5_fe_timing.csv"), &timing).unwrap();
        std::fs::write(out.join("s5_fe_log.txt"), &log).unwrap();
        std::fs::write(out.join("s5_info.txt"), format!("f_single_GHz={:.8}\nf_nominal_w3_GHz={:.8}\n", f_single * 1e-9, f_nominal * 1e-9)).unwrap();
    }
}
