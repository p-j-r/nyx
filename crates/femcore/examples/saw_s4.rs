//! Stage S4 of the topological SAW study, the FE part: supermodes of the 10-ridge Thouless-pump chain at fixed pump phases phi, lambda = 3 um,
//! q0 = 2 pi / 3 um (300 nm Si ridges, width = lambda, psi = 0). The nine gaps at every phase come from `design_s4.py` (`s4_design.csv`:
//! kappa_m(phi) = kappa0 (1 + 0.7 cos(2 pi b m + phi)), b = 2/3, mapped to gaps through the S2 kappa(g) table). The chain has no mirror
//! symmetry, so every phase is one full-domain solve. For every ridge-bound mode: frequency, Q, elastic energy per ridge region, and the complex
//! projection amplitudes on the single-ridge mode of each ridge (as in S3); `postproc_s4.py` builds H_eff(phi) and propagates the pump.
//!   cargo run --release -p femcore --example saw_s4 -- [out_dir]
//! Environment: PHIS (phi/pi list; default: every row of s4_design.csv), N_MODES (80), V_FILM (5790), SHIFT (0.98), SIDE (5), DEPTH (5), MESH_SCALE.
mod saw_array;
mod saw_common;
use femcore::csparse::C;
use femcore::problem::*;
use saw_array::*;
use saw_common::*;
use std::fmt::Write as _;

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "runs/saw_topo".into()));
    let lam = 3.0 * UM;
    let (scale, nm, v_film, shift, side, depth) = (envf("MESH_SCALE", 1.0), envf("N_MODES", 80.0) as usize, envf("V_FILM", 5790.0), envf("SHIFT", 0.98), envf("SIDE", 5.0), envf("DEPTH", 5.0));
    let psi = envf("PSI", 0.0);
    let q0 = 2.0 * std::f64::consts::PI / lam;
    let want: Option<Vec<f64>> = std::env::var("PHIS").ok().map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect());
    let mut log = String::new();
    macro_rules! say {
        ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); log.push_str(&s); log.push('\n'); }};
    }
    // design table: phi/pi, g1..g9, ...
    let design: Vec<(f64, Vec<f64>)> = std::fs::read_to_string(out.join("s4_design.csv"))
        .expect("s4_design.csv (run design_s4.py first)")
        .lines()
        .skip(1)
        .map(|l| {
            let v: Vec<f64> = l.split(',').map(|x| x.parse().unwrap()).collect();
            (v[0], v[1..10].iter().map(|g| g * UM).collect())
        })
        .collect();
    let mut timing = String::from("phi_over_pi,dofs,cells,attempt,solve_s\n");

    // ---- single-ridge reference mode
    let (mesh1, arr1) = build_array(lam, &[], false, scale, side, depth);
    let outer1 = outer_labels(&mesh1, &|p| (p[0].abs() - arr1.xo).abs() < 1e-12 || (p[1] + arr1.d + arr1.pt).abs() < 1e-12);
    let mut spec1 = base_spec(psi, false);
    spec1.bcs = vec![Bc { on: None, labels: outer1, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let sol1 = solve_q(&mesh1, &spec1, q0, 40, 0.985 * v_film * q0 / (2.0 * std::f64::consts::PI)).unwrap();
    let (cx1, labels1) = (cell_x(&mesh1), mesh1.cell_labels.clone());
    let tot1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|_| true);
    let core1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|c| !PML_DOMAINS.contains(&labels1[c]) && cx1[c].abs() < lam);
    let pml1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|c| PML_DOMAINS.contains(&labels1[c]));
    let mut o1: Vec<usize> = (0..sol1.modes.freqs.len()).collect();
    o1.sort_by(|&a, &b| sol1.modes.freqs[a].re.partial_cmp(&sol1.modes.freqs[b].re).unwrap());
    let j1 = *o1.iter().find(|&&j| core1[j] / tot1[j] >= CORE_MIN && pml1[j] / tot1[j] <= PML_MAX).expect("single ridge mode");
    let f_single = sol1.modes.freqs[j1].re;
    say!("single ridge: f = {:.6} GHz ({} dofs, {:.1} s)", f_single * 1e-9, sol1.dofs, sol1.seconds);
    writeln!(timing, "single,{},{},1,{:.2}", sol1.dofs, mesh1.n_cells(), sol1.seconds).unwrap();
    let phi: Vec<[C; 3]> = {
        let v = vertex_vectors(&mesh1, &sol1.modes, j1);
        let w = vertex_weights(&mesh1);
        let norm = v.iter().zip(&w).map(|(a, w)| w * (a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr())).sum::<f64>().sqrt();
        v.into_iter().map(|a| [a[0] / norm, a[1] / norm, a[2] / norm]).collect()
    };
    let loc1 = Locator::new(&mesh1);

    let mut modes_csv = String::from("phi_over_pi,rank,f_GHz,f_im_GHz,df_from_single_MHz,Q,v_m_per_s,si_frac,pml_frac,core_frac,left2_energy,right2_energy,left2_proj,right2_proj");
    for m in 1..=10 {
        write!(modes_csv, ",w{m}").unwrap();
    }
    modes_csv.push('\n');
    let mut proj_csv = String::from("phi_over_pi,rank,f_GHz");
    for m in 1..=10 {
        write!(proj_csv, ",cre{m},cim{m}").unwrap();
    }
    proj_csv.push('\n');

    for (phi_pi, gaps) in &design {
        if let Some(w) = &want {
            if !w.iter().any(|x| (x - phi_pi).abs() < 1e-9) {
                continue;
            }
        }
        let n = gaps.len() + 1;
        let (mesh, arr) = build_array(lam, gaps, false, scale, side, depth);
        let core_half = arr.xe + 0.5 * lam;
        let outer = outer_labels(&mesh, &|p| (p[0].abs() - arr.xo).abs() < 1e-12 || (p[1] + arr.d + arr.pt).abs() < 1e-12);
        let mut spec = base_spec(psi, false);
        spec.bcs = vec![Bc { on: None, labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        let (cx, labels, wv) = (cell_x(&mesh), mesh.cell_labels.clone(), vertex_weights(&mesh));
        let bounds: Vec<(f64, f64)> = (0..n)
            .map(|j| (if j == 0 { f64::NEG_INFINITY } else { 0.5 * (arr.xc[j - 1] + arr.xc[j]) }, if j + 1 == n { f64::INFINITY } else { 0.5 * (arr.xc[j] + arr.xc[j + 1]) }))
            .collect();
        say!("\n== phi = {phi_pi} pi: gaps {:?} um; {} cells, box width {:.1} um", gaps.iter().map(|g| (g / UM * 100.0).round() / 100.0).collect::<Vec<_>>(), mesh.n_cells(), 2.0 * arr.xo / UM);
        // projection references: every ridge region
        let proj_ref: Vec<(Vec<usize>, Vec<[C; 3]>, f64)> = (0..n)
            .map(|j| {
                let (lo, hi) = bounds[j];
                let idx: Vec<usize> = (0..mesh.n_points()).filter(|&i| mesh.point(i)[0] >= lo && mesh.point(i)[0] < hi).collect();
                let vals: Vec<[C; 3]> = idx.iter().map(|&i| loc1.eval(&phi, mesh.point(i)[0] - arr.xc[j], mesh.point(i)[1])).collect();
                let nphi: f64 = idx.iter().zip(&vals).map(|(&i, a)| wv[i] * (a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr())).sum();
                (idx, vals, nphi)
            })
            .collect();
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
            say!("  attempt {}, shift {:.3} GHz: {} dofs, {:.1} s, {} modes", attempt + 1, f_search * 1e-9, sol.dofs, sol.seconds, sol.modes.freqs.len());
            writeln!(timing, "{phi_pi},{},{},{},{:.2}", sol.dofs, mesh.n_cells(), attempt + 1, sol.seconds).unwrap();
            let m = &sol.modes;
            let tot = elastic_energy_where(&mesh, m, 2, &|_| true);
            let pmlf = elastic_energy_where(&mesh, m, 2, &|c| PML_DOMAINS.contains(&labels[c]));
            let core = elastic_energy_where(&mesh, m, 2, &|c| !PML_DOMAINS.contains(&labels[c]) && cx[c].abs() < core_half);
            let mut ord: Vec<usize> = (0..m.freqs.len()).collect();
            ord.sort_by(|&a, &b| m.freqs[a].re.partial_cmp(&m.freqs[b].re).unwrap());
            let cands: Vec<usize> = ord.iter().copied().filter(|&j| core[j] / tot[j] >= CORE_MIN && pmlf[j] / tot[j] <= PML_MAX).collect();
            if cands.len() < n {
                say!("    only {} ridge-bound modes (need {n}), retrying", cands.len());
                continue;
            }
            if cands.len() > n {
                say!("    {} candidates, keeping the lowest {n}", cands.len());
            }
            let sel: Vec<usize> = cands.into_iter().take(n).collect();
            let sif = elastic_energy_where(&mesh, m, 2, &|c| labels[c] == 2);
            let e_reg: Vec<Vec<f64>> = (0..n).map(|r| elastic_energy_where(&mesh, m, 2, &|c| cx[c] >= bounds[r].0 && cx[c] < bounds[r].1)).collect();
            // mode shapes of the most left- and most right-localised modes (by projected weight on the two outer ridges of that end)
            let mut rows: Vec<(f64, usize, Vec<f64>, Vec<C>)> = vec![];
            for &j in &sel {
                let u = vertex_vectors(&mesh, m, j);
                let u2: f64 = (0..mesh.n_points()).map(|i| wv[i] * (u[i][0].norm_sqr() + u[i][1].norm_sqr() + u[i][2].norm_sqr())).sum();
                let c: Vec<C> = (0..n)
                    .map(|r| {
                        let (idx, vals, nphi) = &proj_ref[r];
                        let mut num = C::new(0.0, 0.0);
                        for (&i, a) in idx.iter().zip(vals) {
                            num += (a[0].conj() * u[i][0] + a[1].conj() * u[i][1] + a[2].conj() * u[i][2]) * wv[i];
                        }
                        num / (nphi.sqrt() * u2.sqrt())
                    })
                    .collect();
                let w: Vec<f64> = (0..n).map(|r| e_reg[r][j] / tot[j]).collect();
                rows.push((m.freqs[j].re, j, w, c));
            }
            rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            say!("  rank  f - f_s (MHz)   Q        Si     PML    left2 energy/proj  right2 energy/proj");
            let mut best_l = (0.0, 0usize);
            let mut best_r = (0.0, 0usize);
            for (rank, (f, j, w, c)) in rows.iter().enumerate() {
                let cs: f64 = c.iter().map(|z| z.norm_sqr()).sum();
                let (le, re) = (w[0] + w[1], w[n - 2] + w[n - 1]);
                let (lp, rp) = ((c[0].norm_sqr() + c[1].norm_sqr()) / cs, (c[n - 2].norm_sqr() + c[n - 1].norm_sqr()) / cs);
                if lp > best_l.0 {
                    best_l = (lp, *j);
                }
                if rp > best_r.0 {
                    best_r = (rp, *j);
                }
                let fz = m.freqs[*j];
                say!("  {:>3}   {:+9.3}   {:.1e}  {:.3}  {:.0e}   {:.3}/{:.3}   {:.3}/{:.3}", rank + 1, (f - f_single) * 1e-6, q_factor(fz), sif[*j] / tot[*j], pmlf[*j] / tot[*j], le, lp, re, rp);
                write!(modes_csv, "{phi_pi},{},{:.8},{:.4e},{:.5},{:.3e},{:.1},{:.4},{:.3e},{:.4},{:.5},{:.5},{:.5},{:.5}", rank + 1, f * 1e-9, fz.im * 1e-9, (f - f_single) * 1e-6, q_factor(fz), 2.0 * std::f64::consts::PI * f / q0, sif[*j] / tot[*j], pmlf[*j] / tot[*j], core[*j] / tot[*j], le, re, lp, rp).unwrap();
                for x in w {
                    write!(modes_csv, ",{:.5}", x).unwrap();
                }
                modes_csv.push('\n');
                write!(proj_csv, "{phi_pi},{},{:.8}", rank + 1, f * 1e-9).unwrap();
                for z in c {
                    write!(proj_csv, ",{:.6e},{:.6e}", z.re, z.im).unwrap();
                }
                proj_csv.push('\n');
            }
            // shapes at the end phases of the sweep
            if (phi_pi - 0.45).abs() < 1e-9 || (phi_pi - 2.21).abs() < 1e-9 {
                let n1 = mesh.n_dofs(2);
                for (label, (_, j)) in [("left", best_l), ("right", best_r)] {
                    let u = &m.vectors[j];
                    let big = (0..mesh.n_points()).max_by(|&a, &b| u[n1 + a].norm().partial_cmp(&u[n1 + b].norm()).unwrap()).unwrap();
                    let ph = u[n1 + big].conj() / u[n1 + big].norm();
                    let mut pts = String::from("x_um,z_um,uz_signed,uz_abs\n");
                    for i in 0..mesh.n_points() {
                        let p = mesh.point(i);
                        writeln!(pts, "{:.5},{:.5},{:.6e},{:.6e}", p[0] / UM, p[1] / UM, (u[n1 + i] * ph).re, u[n1 + i].norm()).unwrap();
                    }
                    std::fs::write(out.join(format!("s4_shape_phi{phi_pi}_{label}_pts.csv")), pts).unwrap();
                    say!("  mode shape {label}: {:.6} GHz", m.freqs[j].re * 1e-9);
                }
                let mut tri = String::from("a,b,c,domain\n");
                for c in 0..mesh.n_cells() {
                    let v = mesh.cell(c);
                    writeln!(tri, "{},{},{},{}", v[0], v[1], v[2], labels[c]).unwrap();
                }
                std::fs::write(out.join(format!("s4_shape_phi{phi_pi}_tri.csv")), tri).unwrap();
            }
            done = true;
            break;
        }
        if !done {
            say!("  ridge-bound modes not found at phi = {phi_pi} pi");
        }
        std::fs::write(out.join("s4_modes.csv"), &modes_csv).unwrap();
        std::fs::write(out.join("s4_proj.csv"), &proj_csv).unwrap();
        std::fs::write(out.join("s4_timing.csv"), &timing).unwrap();
        std::fs::write(out.join("s4_log.txt"), &log).unwrap();
        std::fs::write(out.join("s4_info.txt"), format!("f_single_GHz={:.8}\n", f_single * 1e-9)).unwrap();
    }
}
