//! Stage S0 of the topological SAW study: surface acoustic modes of r-plane sapphire (a) bare and (b) under a uniform 300 nm
//! silicon (001) film, as eigenfrequencies at the wave numbers q = 2 pi / lambda (lambda = 2, 3, 4 um) along the wafer y axis.
//!
//! The cross-section is a narrow strip (0.2 lambda wide: the first lateral standing wave, kx = pi/w, is then far above the window) of sapphire 5 lambda deep (+ a 1 lambda absorbing layer below), with
//! the silicon film (b) on top. The sapphire cut has a mirror plane normal to the lateral axis x, so at zero lateral wave
//! number the sagittal modes (u_y, u_z) and the shear-horizontal modes (u_x) decouple exactly, and the strip sides carry the
//! two exact conditions of the uniform (kx = 0) wave: "sag" fixes u_x on the sides (u_y, u_z free), "sh" fixes u_y and u_z
//! (u_x free). Each is solved separately.
//!   cargo run --release -p femcore --example saw_s0 -- [out_dir]
//! Environment: LAMBDAS (um, comma list), MESH_SCALE, N_MODES (30), PML_STRENGTH (3), SHIFT (shift as a fraction of the bulk
//! threshold velocity, 0.8), DEPTH (sapphire depth in wavelengths, 5), PSI (azimuth, 0).
mod saw_common;
use femcore::mesher::*;
use femcore::problem::*;
use saw_common::*;
use std::fmt::Write as _;

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "runs/saw_topo".into()));
    std::fs::create_dir_all(&out).unwrap();
    let lambdas: Vec<f64> = std::env::var("LAMBDAS").ok().map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![2.0, 3.0, 4.0]);
    let (scale, n_modes, shift) = (envf("MESH_SCALE", 1.0), envf("N_MODES", 30.0) as usize, envf("SHIFT", 0.8));
    let depth_lam = envf("DEPTH", 5.0);
    let psi = envf("PSI", 0.0);

    let b = sapphire_bulk(psi);
    let mut log = String::new();
    macro_rules! say {
        ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); log.push_str(&s); log.push('\n'); }};
    }
    say!("sapphire r-plane cut, psi = {psi} deg: bulk waves for propagation along y");
    say!("  at t = 0 (along the surface): SH (u_x) {:.0} m/s, quasi-shear sagittal {:.0} m/s, quasi-longitudinal {:.0} m/s", b.v_sh, b.v_sv, b.v_l);
    say!("  bulk thresholds (min over the sagittal angle of v/cos t): SH {:.0}, sagittal shear {:.0}, longitudinal {:.0} m/s", b.thr_sh, b.thr_sv, b.thr_l);
    say!("  mirror check: max Gamma_xy,xz / Gamma_xx = {:.1e}, max even-odd stiffness coupling = {:.1e} GPa", b.mirror_gamma, b.mirror_c);
    say!("  slowest bulk shear velocity along y: {:.0} m/s", b.v_sh.min(b.v_sv));

    let mut csv = String::from("lambda_um,film,bc,q_rad_per_m,rank,f_GHz,f_im_GHz,Q,v_m_per_s,v_over_slow_shear,surface_frac_top_lambda,pml_frac,si_frac,class\n");
    let mut summary = String::new();
    for &lam_um in &lambdas {
        let lam = lam_um * UM;
        let q = 2.0 * std::f64::consts::PI / lam;
        for film in [false, true] {
            // ---- mesh: x in [-w/2, w/2], y vertical (wafer z), sapphire below 0, film 0..300 nm, absorbing layer under the sapphire
            let (w, d, pt) = (envf("STRIP", 0.2) * lam, depth_lam * lam, lam);
            let (h_far, h_surf, h_si) = (lam / 8.0 * scale, lam / 16.0 * scale, 75e-9 * scale);
            let rect = |x0: f64, y0: f64, x1: f64, y1: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: dom, h };
            let mut parts = vec![
                rect(-w / 2.0, -d, w / 2.0, 0.0, 1, None),
                rect(-w / 2.0, -0.6 * lam, w / 2.0, 0.0, 1, Some(h_surf)),
                rect(-w / 2.0, -d - pt, w / 2.0, -d, 6, None),
            ];
            if film {
                parts.push(rect(-w / 2.0, 0.0, w / 2.0, T_SI, 2, Some(h_si)));
            }
            let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
            let (sides, bottom) = (outer_labels(&mesh, &|p| (p[0].abs() - w / 2.0).abs() < 1e-12), outer_labels(&mesh, &|p| (p[1] + d + pt).abs() < 1e-12));
            let doms: Vec<i32> = if film { vec![1, 2, 6] } else { vec![1, 6] };
            let mut base = ProblemSpec::new(Physics::Solid);
            base.order = 2;
            set_materials(&mut base, psi, &doms, &[1, 6], if film { &[2] } else { &[] });
            base.materials.pml.insert(6, pml(PmlDir::Off, PmlDir::Minus));
            let v_thr = b.thr_sh.min(b.thr_sv);
            let f_search = shift * v_thr * q / (2.0 * std::f64::consts::PI);
            say!("\n== lambda {lam_um} um, {}: mesh {} cells, shift {:.3} GHz", if film { "300 nm Si film" } else { "bare sapphire" }, mesh.n_cells(), f_search * 1e-9);
            for (bc_name, fixed_sides) in [("sag", [true, false, false]), ("sh", [false, true, true])] {
                let mut spec = base.clone();
                spec.bcs = vec![
                    Bc { on: None, labels: bottom.clone(), kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
                    Bc { on: None, labels: sides.clone(), kind: BcKind::Fixed { fixed: fixed_sides, values: [0.0; 3] } },
                ];
                let sol = match solve_q(&mesh, &spec, q, n_modes, f_search) {
                    Ok(s) => s,
                    Err(e) => {
                        say!("  {bc_name}: failed: {e}");
                        continue;
                    }
                };
                say!("  {bc_name}: {} dofs, {:.2} s, {} modes converged", sol.dofs, sol.seconds, sol.modes.freqs.len());
                let m = &sol.modes;
                let labels = mesh.cell_labels.clone();
                let tot = elastic_energy_where(&mesh, m, 2, &|_| true);
                let top = elastic_energy_where(&mesh, m, 2, &|c| labels[c] != 6 && centroid_y(&mesh, c) > -lam);
                let pmlf = elastic_energy_where(&mesh, m, 2, &|c| labels[c] == 6);
                let sif = elastic_energy_where(&mesh, m, 2, &|c| labels[c] == 2);
                let mut order: Vec<usize> = (0..m.freqs.len()).collect();
                order.sort_by(|&a, &b| m.freqs[a].re.partial_cmp(&m.freqs[b].re).unwrap());
                let slow = b.v_sh.min(b.v_sv);
                let mut rank = 0;
                let mut best = vec![];
                for &j in &order {
                    let f = m.freqs[j];
                    let v = f.re * lam;
                    let (sf, pf, sif_) = (top[j] / tot[j], pmlf[j] / tot[j], sif[j] / tot[j]);
                    let qq = q_factor(f);
                    let class = if sf >= 0.5 { if qq > 1e4 { "surface-bound" } else { "surface-leaky" } } else { "bulk/PML" };
                    rank += 1;
                    writeln!(csv, "{lam_um},{},{bc_name},{q:.6e},{rank},{:.6},{:.4e},{:.4e},{:.1},{:.4},{:.4},{:.4},{:.4},{class}", film as u8, f.re * 1e-9, f.im * 1e-9, qq, v, v / slow, sf, pf, sif_).unwrap();
                    if class != "bulk/PML" && best.len() < 3 {
                        best.push(format!("{:.4} GHz v={:.0} m/s (v/slow shear {:.3}) Q={:.3e} surf {:.2} pml {:.3} Si {:.2} {class}", f.re * 1e-9, v, v / slow, qq, sf, pf, sif_));
                    }
                }
                if best.is_empty() {
                    say!("    no surface-like mode among the {} lowest-distance modes (lowest f {:.4} GHz, v {:.0} m/s)", order.len(), m.freqs[order[0]].re * 1e-9, m.freqs[order[0]].re * lam);
                }
                for l in &best {
                    say!("    {l}");
                    writeln!(summary, "lambda {lam_um} um | {} | {bc_name} | {l} | {} dofs {:.2} s", if film { "film" } else { "bare" }, sol.dofs, sol.seconds).unwrap();
                }
            }
        }
    }
    let tag = std::env::var("TAG").unwrap_or_default();
    std::fs::write(out.join(format!("s0_modes{tag}.csv")), csv).unwrap();
    std::fs::write(out.join(format!("s0_summary{tag}.txt")), summary).unwrap();
    std::fs::write(out.join(format!("s0_log{tag}.txt")), log).unwrap();
}
