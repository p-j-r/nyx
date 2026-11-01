//! Stage S1 of the topological SAW study: the ridge-bound SAW mode of one 300 nm silicon (001) ridge, fully etched to the
//! r-plane sapphire, width = lambda (lambda = 2, 3, 4 um), crystal azimuth psi = 0. Eigenfrequencies at q = q0 (0.8, 0.9, 1, 1.1, 1.2),
//! q0 = 2 pi / lambda. Sapphire 5 lambda deep, absorbing layers (1 lambda) on both sides and below, fixed outer boundary.
//!
//! The ridge-bound mode is the lowest-frequency mode with its elastic energy concentrated in the ridge: at least `CORE_MIN`
//! (0.5) of it within the ridge plus 0.5 lambda on each side, at most `PML_MAX` (0.2) in the absorbing layers. For it the program
//! writes f, Q (from the complex frequency), the elastic energy fractions in the silicon and in the layers, and the lateral decay
//! of |u| along the sapphire surface beside the ridge (fit of ln|u| against the distance from the ridge edge, per side).
//!   cargo run --release -p femcore --example saw_s1 -- [out_dir]
//! Environment: LAMBDAS (um), QRATIOS, MESH_SCALE, N_MODES (40), V_FILM (m/s, 5790: the shift is SHIFT * V_FILM q / 2 pi, SHIFT 0.985),
//! LSIDE (sapphire beside the ridge in wavelengths, 3), DEPTH (5), PML_STRENGTH (3), PSI (0), WIDTH (ridge width in wavelengths, 1), TAG.
mod saw_common;
use femcore::mesher::*;
use femcore::problem::*;
use saw_common::*;
use std::fmt::Write as _;

const CORE_MIN: f64 = 0.5;
const PML_MAX: f64 = 0.2;
const PML_DOMAINS: [i32; 5] = [4, 5, 6, 8, 9];

/// Decay of |u| along the sapphire surface (y = 0) beside the ridge on one side: (decay length m, |u| at one wavelength from the
/// edge relative to the edge, number of points in the fit). The fit runs from 0.1 um to where |u| falls below 3 % of the edge value
/// or the physical region ends.
fn lateral_decay(surface: &[(f64, f64)], lam: f64) -> (f64, f64, usize, f64) {
    // surface: (distance from the ridge edge, |u|) sorted by distance, distance > 0
    let edge = surface.iter().find(|p| p.0 >= 0.0).map_or(f64::NAN, |p| p.1);
    let pts: Vec<&(f64, f64)> = surface.iter().filter(|p| p.0 >= 0.1 * UM && p.1 > 0.03 * edge).collect();
    if pts.len() < 4 {
        return (f64::NAN, f64::NAN, pts.len(), f64::NAN);
    }
    let n = pts.len() as f64;
    let (sx, sy): (f64, f64) = pts.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0, a.1 + p.1.ln()));
    let (mx, my) = (sx / n, sy / n);
    let (sxx, sxy) = pts.iter().fold((0.0, 0.0), |a, p| (a.0 + (p.0 - mx).powi(2), a.1 + (p.0 - mx) * (p.1.ln() - my)));
    let slope = sxy / sxx;
    let at_lambda = {
        let target = lam;
        let mut v = f64::NAN;
        for w in surface.windows(2) {
            if w[0].0 <= target && w[1].0 >= target {
                v = (w[0].1 + (w[1].1 - w[0].1) * (target - w[0].0) / (w[1].0 - w[0].0)) / edge;
            }
        }
        v
    };
    // distance at which |u| first falls to 1/e of the edge value (no model assumed)
    let mut efold = f64::NAN;
    for w in surface.windows(2) {
        let (a, b) = (w[0].1 / edge, w[1].1 / edge);
        if efold.is_nan() && a > (-1.0f64).exp() && b <= (-1.0f64).exp() {
            efold = w[0].0 + (w[1].0 - w[0].0) * (a - (-1.0f64).exp()) / (a - b);
        }
    }
    (if slope < 0.0 { -1.0 / slope } else { f64::NAN }, at_lambda, pts.len(), efold)
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "runs/saw_topo".into()));
    std::fs::create_dir_all(&out).unwrap();
    let lambdas: Vec<f64> = std::env::var("LAMBDAS").ok().map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![2.0, 3.0, 4.0]);
    let qratios: Vec<f64> = std::env::var("QRATIOS").ok().map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![0.8, 0.9, 1.0, 1.1, 1.2]);
    let (scale, n_modes) = (envf("MESH_SCALE", 1.0), envf("N_MODES", 40.0) as usize);
    let (v_film, shift, lside, depth_lam, wid) = (envf("V_FILM", 5790.0), envf("SHIFT", 0.985), envf("LSIDE", 3.0), envf("DEPTH", 5.0), envf("WIDTH", 1.0));
    let psi = envf("PSI", 0.0);
    let tag = std::env::var("TAG").unwrap_or_default();

    let mut log = String::new();
    macro_rules! say {
        ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); log.push_str(&s); log.push('\n'); }};
    }
    let b = sapphire_bulk(psi);
    say!("bulk thresholds of the sapphire cut (m/s): SH {:.0}, sagittal shear {:.0}, longitudinal {:.0}; slowest bulk shear along y {:.0}", b.thr_sh, b.thr_sv, b.thr_l, b.v_sh.min(b.v_sv));

    let mut csv = String::from("lambda_um,width_um,q_rad_per_m,q_over_q0,f_GHz,f_im_GHz,Q,v_m_per_s,si_frac,pml_frac,core_frac,decay_len_left_um,decay_len_right_um,decay_len_fit_um,efold_um,u_at_1lambda_over_edge,fit_points,n_candidates,dofs,solve_s\n");
    let mut all = String::from("lambda_um,q_over_q0,rank,f_GHz,f_im_GHz,Q,v_m_per_s,si_frac,pml_frac,core_frac,candidate\n");
    for &lam_um in &lambdas {
        let lam = lam_um * UM;
        let q0 = 2.0 * std::f64::consts::PI / lam;
        let w = wid * lam;
        let (xin, pt, d) = (0.5 * w + lside * lam, lam, depth_lam * lam);
        let xo = xin + pt;
        let (h_far, h_surf, h_si, h_edge) = (lam / 8.0 * scale, lam / 16.0 * scale, 75e-9 * scale, 75e-9 * scale);
        let rect = |x0: f64, y0: f64, x1: f64, y1: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: dom, h };
        let parts = vec![
            rect(-xin, -d, xin, 0.0, 1, None),
            rect(-xin, -0.6 * lam, xin, 0.0, 1, Some(h_surf)),
            rect(-0.5 * w - 0.25 * UM, -0.25 * UM, 0.5 * w + 0.25 * UM, 0.0, 1, Some(h_edge)),
            rect(-xo, -d, -xin, 0.0, 4, None),
            rect(xin, -d, xo, 0.0, 5, None),
            rect(-xin, -d - pt, xin, -d, 6, None),
            rect(-xo, -d - pt, -xin, -d, 8, None),
            rect(xin, -d - pt, xo, -d, 9, None),
            rect(-0.5 * w, 0.0, 0.5 * w, T_SI, 2, Some(h_si)),
        ];
        let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
        let outer = outer_labels(&mesh, &|p| (p[0].abs() - xo).abs() < 1e-12 || (p[1] + d + pt).abs() < 1e-12);
        let mut base = ProblemSpec::new(Physics::Solid);
        base.order = 2;
        set_materials(&mut base, psi, &[1, 2, 4, 5, 6, 8, 9], &[1, 4, 5, 6, 8, 9], &[2]);
        for (dom, px, py) in [(4, PmlDir::Minus, PmlDir::Off), (5, PmlDir::Plus, PmlDir::Off), (6, PmlDir::Off, PmlDir::Minus), (8, PmlDir::Minus, PmlDir::Minus), (9, PmlDir::Plus, PmlDir::Minus)] {
            base.materials.pml.insert(dom, pml(px, py));
        }
        base.bcs = vec![Bc { on: None, labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        let labels = mesh.cell_labels.clone();
        let cx: Vec<f64> = (0..mesh.n_cells()).map(|c| mesh.cell(c).iter().map(|&v| mesh.point(v as usize)[0]).sum::<f64>() / 3.0).collect();
        let core_half = 0.5 * w + 0.5 * lam;
        say!("\n== lambda {lam_um} um, ridge {:.2} um x 0.3 um: mesh {} cells, physical half width {:.1} um, depth {:.1} um", w / UM, mesh.n_cells(), xin / UM, d / UM);
        for &qr in &qratios {
            let q = qr * q0;
            let mut sol = None;
            let mut cand: Vec<usize> = vec![];
            let mut metrics = (vec![], vec![], vec![], vec![]);
            let mut order: Vec<usize> = vec![];
            let mut total_s = 0.0;
            for (attempt, (sh, nm)) in [(shift, n_modes), (0.94, n_modes + 20), (1.02, n_modes + 20)].into_iter().enumerate() {
                let f_search = sh * v_film * q / (2.0 * std::f64::consts::PI);
                let s = match solve_q(&mesh, &base, q, nm, f_search) {
                    Ok(s) => s,
                    Err(e) => {
                        say!("  q/q0 {qr}: solve failed: {e}");
                        continue;
                    }
                };
                total_s += s.seconds;
                say!("  q/q0 {qr} (attempt {}, shift {:.3} GHz): {} dofs, {:.2} s, {} modes", attempt + 1, f_search * 1e-9, s.dofs, s.seconds, s.modes.freqs.len());
                let m = &s.modes;
                let tot = elastic_energy_where(&mesh, m, 2, &|_| true);
                let sif = elastic_energy_where(&mesh, m, 2, &|c| labels[c] == 2);
                let pmlf = elastic_energy_where(&mesh, m, 2, &|c| PML_DOMAINS.contains(&labels[c]));
                let core = elastic_energy_where(&mesh, m, 2, &|c| !PML_DOMAINS.contains(&labels[c]) && cx[c].abs() < core_half);
                order = (0..m.freqs.len()).collect();
                order.sort_by(|&a, &b| m.freqs[a].re.partial_cmp(&m.freqs[b].re).unwrap());
                cand = order.iter().copied().filter(|&j| core[j] / tot[j] >= CORE_MIN && pmlf[j] / tot[j] <= PML_MAX).collect();
                metrics = (tot.clone(), sif, pmlf, core);
                sol = Some(s);
                if !cand.is_empty() {
                    break;
                }
            }
            let Some(s) = sol else { continue };
            let m = &s.modes;
            let (tot, sif, pmlf, core) = &metrics;
            for (rank, &j) in order.iter().enumerate() {
                let f = m.freqs[j];
                writeln!(all, "{lam_um},{qr},{},{:.6},{:.4e},{:.4e},{:.1},{:.4},{:.4},{:.4},{}", rank + 1, f.re * 1e-9, f.im * 1e-9, q_factor(f), f.re * 2.0 * std::f64::consts::PI / q, sif[j] / tot[j], pmlf[j] / tot[j], core[j] / tot[j], cand.contains(&j) as u8).unwrap();
            }
            let Some(&j) = cand.first() else {
                say!("    no ridge-bound mode found among {} modes (lowest {:.4} GHz)", order.len(), m.freqs[order[0]].re * 1e-9);
                writeln!(csv, "{lam_um},{:.3},{q:.6e},{qr},NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,{},{total_s:.2}", w / UM, s.dofs).unwrap();
                continue;
            };
            let f = m.freqs[j];
            let vf = vertex_fields(&mesh, m, 2, j);
            // |u| along the sapphire surface beside the ridge, per side
            let side = |sign: f64| -> Vec<(f64, f64)> {
                let mut v: Vec<(f64, f64)> = (0..mesh.n_points())
                    .filter_map(|i| {
                        let p = mesh.point(i);
                        let dist = sign * p[0] - 0.5 * w;
                        (p[1].abs() < 1e-12 && dist > -1e-12 && dist < lside * lam + 1e-12).then_some((dist.max(0.0), vf[i][0]))
                    })
                    .collect();
                v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                v
            };
            let (dl, dr) = (lateral_decay(&side(-1.0), lam), lateral_decay(&side(1.0), lam));
            let avg = |a: f64, b: f64| 0.5 * (a + b);
            say!(
                "    ridge mode {:.5} GHz (Im {:+.2e}) Q {:.3e} v {:.0} m/s | Si {:.3} PML {:.2e} core {:.3} | log-fit decay length {:.2}/{:.2} um (L/R), 1/e distance {:.2}/{:.2} um, |u|(edge+lambda)/edge {:.3}/{:.3}",
                f.re * 1e-9, f.im * 1e-9, q_factor(f), 2.0 * std::f64::consts::PI * f.re / q, sif[j] / tot[j], pmlf[j] / tot[j], core[j] / tot[j], dl.0 * 1e6, dr.0 * 1e6, dl.3 * 1e6, dr.3 * 1e6, dl.1, dr.1
            );
            writeln!(
                csv,
                "{lam_um},{:.3},{q:.6e},{qr},{:.6},{:.4e},{:.4e},{:.1},{:.4},{:.4e},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{},{},{},{:.2}",
                w / UM,
                f.re * 1e-9,
                f.im * 1e-9,
                q_factor(f),
                2.0 * std::f64::consts::PI * f.re / q,
                sif[j] / tot[j],
                pmlf[j] / tot[j],
                core[j] / tot[j],
                dl.0 * 1e6,
                dr.0 * 1e6,
                avg(dl.0, dr.0) * 1e6,
                avg(dl.3, dr.3) * 1e6,
                avg(dl.1, dr.1),
                dl.2 + dr.2,
                cand.len(),
                s.dofs,
                total_s
            )
            .unwrap();
            if (qr - 1.0).abs() < 1e-9 {
                // mode shape at q0: vertices (x, z in um, |u|, |ux|, |uy| propagation, |uz| normal) and triangles
                let mut pts = String::from("x_um,z_um,u_abs,ux_abs,uy_abs,uz_abs\n");
                for i in 0..mesh.n_points() {
                    let p = mesh.point(i);
                    writeln!(pts, "{:.5},{:.5},{:.6e},{:.6e},{:.6e},{:.6e}", p[0] / UM, p[1] / UM, vf[i][0], vf[i][1], vf[i][2], vf[i][3]).unwrap();
                }
                let mut tri = String::from("a,b,c,domain\n");
                for c in 0..mesh.n_cells() {
                    let v = mesh.cell(c);
                    writeln!(tri, "{},{},{},{}", v[0], v[1], v[2], labels[c]).unwrap();
                }
                std::fs::write(out.join(format!("s1_mode_lam{lam_um}_pts{tag}.csv")), pts).unwrap();
                std::fs::write(out.join(format!("s1_mode_lam{lam_um}_tri{tag}.csv")), tri).unwrap();
                let side_csv = |name: &str, v: Vec<(f64, f64)>| {
                    let mut s = String::from("dist_from_edge_um,u_abs\n");
                    for p in v {
                        writeln!(s, "{:.5},{:.6e}", p.0 / UM, p.1).unwrap();
                    }
                    std::fs::write(out.join(format!("s1_surface_{name}_lam{lam_um}{tag}.csv")), s).unwrap();
                };
                side_csv("left", side(-1.0));
                side_csv("right", side(1.0));
            }
        }
    }
    std::fs::write(out.join(format!("s1_sweep{tag}.csv")), csv).unwrap();
    std::fs::write(out.join(format!("s1_all_modes{tag}.csv")), all).unwrap();
    std::fs::write(out.join(format!("s1_log{tag}.txt")), log).unwrap();
}
