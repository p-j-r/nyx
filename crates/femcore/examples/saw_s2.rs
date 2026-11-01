//! Stage S2 of the topological SAW study: two identical 300 nm silicon ridges (width = lambda, fully etched, psi = 0) on r-plane
//! sapphire, centre-to-centre separation lambda + g. Even and odd supermodes at q0 = 2 pi / lambda versus the gap g.
//!
//! The structure and the elastic problem are mirror symmetric about the mid-plane x = 0 (sapphire cut and silicon have an exact
//! mirror plane normal to x). So the supermodes are found on the half domain x >= 0 with one condition on the mid-plane,
//! which gives them on the same mesh (no mesh asymmetry in the tiny splittings at large gaps) and exact parity:
//!   even sector: u_x = 0 on x = 0 (u_y, u_z even),    odd sector: u_y = u_z = 0 on x = 0 (u_z odd).
//! For the gaps in FULL_GAPS the two ridges are also solved on the full domain and the parity is read from u_z at probe points,
//! as a check of the half-domain approach and of the classification.
//!   cargo run --release -p femcore --example saw_s2 -- [out_dir]
//! Environment: LAMBDA (um, 3), GAPS (um, comma list), CHECK_GAPS (also solved at 0.95 q0 and 1.05 q0), FULL_GAPS, DUMP_GAP (mode
//! shapes), MESH_SCALE, N_MODES (50), V_FILM (5790), SHIFT (0.98), SIDE (sapphire beside the outer ridge edges, wavelengths, 5), DEPTH (5).
mod saw_common;
use femcore::csparse::C;
use femcore::mesh::SimplexMesh;
use femcore::mesher::*;
use femcore::problem::*;
use femcore::solid_kz::KzModes;
use saw_common::*;
use std::fmt::Write as _;

const CORE_MIN: f64 = 0.5;
const PML_MAX: f64 = 0.2;
const PML_DOMAINS: [i32; 5] = [4, 5, 6, 8, 9];

struct Geo {
    xc: f64,
    xe: f64,
    xo: f64,
    d: f64,
    pt: f64,
}

fn build(lam: f64, gap: f64, half: bool, scale: f64, side: f64, depth_lam: f64) -> (SimplexMesh, Geo) {
    let w = lam;
    let xc = 0.5 * (lam + gap);
    let xe = xc + 0.5 * w;
    let (xin, pt, d) = (xe + side * lam, lam, depth_lam * lam);
    let xo = xin + pt;
    let (h_far, h_surf, h_si, h_edge) = (lam / 8.0 * scale, lam / 16.0 * scale, 75e-9 * scale, 75e-9 * scale);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: dom, h };
    let xlo = if half { 0.0 } else { -xin };
    let mut parts = vec![rect(xlo, -d, xin, 0.0, 1, None), rect(xlo, -0.6 * lam, xin, 0.0, 1, Some(h_surf))];
    let signs: Vec<f64> = if half { vec![1.0] } else { vec![-1.0, 1.0] };
    // refinement zones around the ridges: one zone when the ridges are close, else one per ridge
    if gap < 1.0 * UM {
        parts.push(rect(xlo.max(-(xe + 0.25 * UM)), -0.25 * UM, xe + 0.25 * UM, 0.0, 1, Some(h_edge)));
    } else {
        for &s in &signs {
            let (a, b) = (s * (xc - 0.5 * w - 0.25 * UM), s * (xc + 0.5 * w + 0.25 * UM));
            parts.push(rect(a.min(b), -0.25 * UM, a.max(b), 0.0, 1, Some(h_edge)));
        }
    }
    parts.push(rect(xin, -d, xo, 0.0, 5, None));
    parts.push(rect(xlo, -d - pt, xin, -d, 6, None));
    parts.push(rect(xin, -d - pt, xo, -d, 9, None));
    if !half {
        parts.push(rect(-xo, -d, -xin, 0.0, 4, None));
        parts.push(rect(-xo, -d - pt, -xin, -d, 8, None));
    }
    for &s in &signs {
        let (a, b) = (s * (xc - 0.5 * w), s * (xc + 0.5 * w));
        parts.push(rect(a.min(b), 0.0, a.max(b), T_SI, 2, Some(h_si)));
    }
    let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
    (mesh, Geo { xc, xe, xo, d, pt })
}

fn base_spec(psi: f64, half: bool) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    let doms = [1, 2, 4, 5, 6, 8, 9];
    set_materials(&mut s, psi, &doms, &[1, 4, 5, 6, 8, 9], &[2]);
    let mut dirs = vec![(5, PmlDir::Plus, PmlDir::Off), (6, PmlDir::Off, PmlDir::Minus), (9, PmlDir::Plus, PmlDir::Minus)];
    if !half {
        dirs.extend([(4, PmlDir::Minus, PmlDir::Off), (8, PmlDir::Minus, PmlDir::Minus)]);
    }
    for (dom, px, py) in dirs {
        s.materials.pml.insert(dom, pml(px, py));
    }
    s
}

#[derive(Clone)]
struct Rec {
    f: C,
    q: f64,
    v: f64,
    si: f64,
    pml: f64,
    core: f64,
    ux_uz: f64,
    parity: f64,
    cand: bool,
}

fn nearest(mesh: &SimplexMesh, x: f64, z: f64) -> usize {
    (0..mesh.n_points())
        .min_by(|&a, &b| {
            let (pa, pb) = (mesh.point(a), mesh.point(b));
            ((pa[0] - x).powi(2) + (pa[1] - z).powi(2)).partial_cmp(&((pb[0] - x).powi(2) + (pb[1] - z).powi(2))).unwrap()
        })
        .unwrap()
}

/// Metrics of every mode, sorted by frequency. `probes`: vertex pairs (x, -x) for the parity of u_z (full domain only).
fn analyze(mesh: &SimplexMesh, m: &KzModes, q: f64, core_half: f64, probes: Option<&[(usize, usize)]>) -> Vec<Rec> {
    let labels = &mesh.cell_labels;
    let cx: Vec<f64> = (0..mesh.n_cells()).map(|c| mesh.cell(c).iter().map(|&v| mesh.point(v as usize)[0]).sum::<f64>() / 3.0).collect();
    let tot = elastic_energy_where(mesh, m, 2, &|_| true);
    let sif = elastic_energy_where(mesh, m, 2, &|c| labels[c] == 2);
    let pmlf = elastic_energy_where(mesh, m, 2, &|c| PML_DOMAINS.contains(&labels[c]));
    let core = elastic_energy_where(mesh, m, 2, &|c| !PML_DOMAINS.contains(&labels[c]) && cx[c].abs() < core_half);
    let n1 = mesh.n_dofs(2);
    let mut si_v: Vec<usize> = (0..mesh.n_cells()).filter(|&c| labels[c] == 2).flat_map(|c| mesh.cell(c).iter().map(|&v| v as usize)).collect();
    si_v.sort();
    si_v.dedup();
    let mut recs: Vec<Rec> = (0..m.freqs.len())
        .map(|j| {
            let u = &m.vectors[j];
            let (ux, uz) = (si_v.iter().map(|&i| u[i].norm_sqr()).sum::<f64>(), si_v.iter().map(|&i| u[n1 + i].norm_sqr()).sum::<f64>());
            let parity = probes.map_or(f64::NAN, |p| {
                let (mut num, mut a, mut b) = (0.0, 0.0, 0.0);
                for &(ip, im) in p {
                    num += (u[n1 + ip] * u[n1 + im].conj()).re;
                    a += u[n1 + ip].norm_sqr();
                    b += u[n1 + im].norm_sqr();
                }
                num / (a * b).sqrt().max(1e-300)
            });
            let f = m.freqs[j];
            let (sf, pf, cf) = (sif[j] / tot[j], pmlf[j] / tot[j], core[j] / tot[j]);
            Rec { f, q: q_factor(f), v: 2.0 * std::f64::consts::PI * f.re / q, si: sf, pml: pf, core: cf, ux_uz: ux / uz.max(1e-300), parity, cand: cf >= CORE_MIN && pf <= PML_MAX }
        })
        .collect();
    // sort by frequency, keeping the index order for the caller through the stable sort of (f, rec)
    recs.sort_by(|a, b| a.f.re.partial_cmp(&b.f.re).unwrap());
    recs
}

fn candidates(s: &Sector) -> Vec<&Rec> {
    s.recs.iter().filter(|r| r.cand).collect()
}

fn order_of(m: &KzModes) -> Vec<usize> {
    let mut o: Vec<usize> = (0..m.freqs.len()).collect();
    o.sort_by(|&a, &b| m.freqs[a].re.partial_cmp(&m.freqs[b].re).unwrap());
    o
}

struct Sector {
    recs: Vec<Rec>,
    sol: Solved,
    /// index into the sorted modes of the first candidate (j in `m.vectors`)
    first: Option<usize>,
    seconds: f64,
}

/// One solve (retrying with other shifts if no ridge-bound mode is found).
fn solve_sector(mesh: &SimplexMesh, spec: &ProblemSpec, q: f64, nm: usize, v_film: f64, shift: f64, core_half: f64, probes: Option<&[(usize, usize)]>) -> Option<Sector> {
    let mut total = 0.0;
    for (sh, n) in [(shift, nm), (0.94, nm + 20), (1.02, nm + 20)] {
        let f_search = sh * v_film * q / (2.0 * std::f64::consts::PI);
        let sol = match solve_q(mesh, spec, q, n, f_search) {
            Ok(s) => s,
            Err(e) => {
                println!("    solve failed: {e}");
                continue;
            }
        };
        total += sol.seconds;
        let recs = analyze(mesh, &sol.modes, q, core_half, probes);
        let ord = order_of(&sol.modes);
        let first = recs.iter().position(|r| r.cand).map(|k| ord[k]);
        if first.is_some() {
            return Some(Sector { recs, sol, first, seconds: total });
        }
    }
    None
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "runs/saw_topo".into()));
    std::fs::create_dir_all(&out).unwrap();
    let lam_um = envf("LAMBDA", 3.0);
    let lam = lam_um * UM;
    let list = |k: &str, d: &[f64]| -> Vec<f64> { std::env::var(k).ok().map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| d.to_vec()) };
    let gaps = list("GAPS", &[0.6, 0.9, 1.2, 1.8, 2.4, 3.6, 4.8, 6.0, 7.2, 9.0]);
    let check_gaps = list("CHECK_GAPS", &[1.2, 2.4, 7.2]);
    let full_gaps = list("FULL_GAPS", &[1.2, 7.2]);
    let dump_gap = envf("DUMP_GAP", 1.2);
    let (scale, nm, v_film, shift, side, depth) = (envf("MESH_SCALE", 1.0), envf("N_MODES", 50.0) as usize, envf("V_FILM", 5790.0), envf("SHIFT", 0.98), envf("SIDE", 5.0), envf("DEPTH", 5.0));
    let psi = envf("PSI", 0.0);
    let q0 = 2.0 * std::f64::consts::PI / lam;
    let tag = std::env::var("TAG").unwrap_or_default();
    let mut log = String::new();
    macro_rules! say {
        ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); log.push_str(&s); log.push('\n'); }};
    }
    let mut modes_csv = String::from("lambda_um,gap_um,domain,sector,q_over_q0,rank,f_GHz,f_im_GHz,Q,v_m_per_s,si_frac,pml_frac,core_frac,uxsq_over_uzsq_in_si,uz_parity,candidate\n");
    let mut super_csv = String::from("lambda_um,gap_um,q_over_q0,q_rad_per_m,f_even_GHz,f_odd_GHz,df_MHz,Q_even,Q_odd,si_even,si_odd,pml_even,pml_odd,core_even,core_odd,f2_even_GHz,f2_odd_GHz,dofs_even,dofs_odd,s_even,s_odd,cells\n");
    let mut full_csv = String::from("lambda_um,gap_um,rank,f_GHz,Q,v_m_per_s,si_frac,pml_frac,core_frac,uxsq_over_uzsq_in_si,uz_parity,candidate,dofs,solve_s\n");

    for &g_um in &gaps {
        let gap = g_um * UM;
        let (mesh, geo) = build(lam, gap, true, scale, side, depth);
        let core_half = geo.xe + 0.5 * lam;
        let sym = outer_labels(&mesh, &|p| p[0].abs() < 1e-12);
        let outer = outer_labels(&mesh, &|p| (p[0] - geo.xo).abs() < 1e-12 || (p[1] + geo.d + geo.pt).abs() < 1e-12);
        say!("\n== lambda {lam_um} um, gap {g_um} um (half domain: {} cells, outer ridge edge at {:.1} um, box half width {:.1} um)", mesh.n_cells(), geo.xe / UM, geo.xo / UM);
        let qratios: Vec<f64> = if check_gaps.iter().any(|&c| (c - g_um).abs() < 1e-9) { vec![1.0, 0.95, 1.05] } else { vec![1.0] };
        for &qr in &qratios {
            let q = qr * q0;
            let mut res: Vec<(&str, Sector)> = vec![];
            for (name, fixed) in [("even", [true, false, false]), ("odd", [false, true, true])] {
                let mut spec = base_spec(psi, true);
                spec.bcs = vec![
                    Bc { on: None, labels: outer.clone(), kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
                    Bc { on: None, labels: sym.clone(), kind: BcKind::Fixed { fixed, values: [0.0; 3] } },
                ];
                let Some(s) = solve_sector(&mesh, &spec, q, nm, v_film, shift, core_half, None) else {
                    say!("  q/q0 {qr} {name}: no ridge-bound mode found");
                    continue;
                };
                say!("  q/q0 {qr} {name}: {} dofs, {:.1} s, {} modes", s.sol.dofs, s.seconds, s.recs.len());
                for (rank, r) in s.recs.iter().enumerate() {
                    writeln!(modes_csv, "{lam_um},{g_um},half,{name},{qr},{},{:.6},{:.4e},{:.4e},{:.1},{:.4},{:.4e},{:.4},{:.4e},{},{}", rank + 1, r.f.re * 1e-9, r.f.im * 1e-9, r.q, r.v, r.si, r.pml, r.core, r.ux_uz, r.parity, r.cand as u8).unwrap();
                }
                if (qr - 1.0).abs() < 1e-9 && (g_um - dump_gap).abs() < 1e-9 {
                    let j = s.first.unwrap();
                    let n1 = mesh.n_dofs(2);
                    let u = &s.sol.modes.vectors[j];
                    // real, signed u_z: the phase of the largest u_z vertex
                    let big = (0..mesh.n_points()).max_by(|&a, &b| u[n1 + a].norm().partial_cmp(&u[n1 + b].norm()).unwrap()).unwrap();
                    let ph = u[n1 + big].conj() / u[n1 + big].norm();
                    let vf = vertex_fields(&mesh, &s.sol.modes, 2, j);
                    let mut pts = String::from("x_um,z_um,u_abs,ux_abs,uy_abs,uz_abs,uz_signed\n");
                    for i in 0..mesh.n_points() {
                        let p = mesh.point(i);
                        writeln!(pts, "{:.5},{:.5},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}", p[0] / UM, p[1] / UM, vf[i][0], vf[i][1], vf[i][2], vf[i][3], (u[n1 + i] * ph).re).unwrap();
                    }
                    let mut tri = String::from("a,b,c,domain\n");
                    for c in 0..mesh.n_cells() {
                        let v = mesh.cell(c);
                        writeln!(tri, "{},{},{},{}", v[0], v[1], v[2], mesh.cell_labels[c]).unwrap();
                    }
                    std::fs::write(out.join(format!("s2_mode_lam{lam_um}_g{g_um}_{name}_pts{tag}.csv")), pts).unwrap();
                    std::fs::write(out.join(format!("s2_mode_lam{lam_um}_g{g_um}_{name}_tri{tag}.csv")), tri).unwrap();
                }
                res.push((name, s));
            }
            if res.len() == 2 {
                let (ce, co) = (candidates(&res[0].1), candidates(&res[1].1));
                let (e, o) = (ce[0], co[0]);
                let f2 = |c: &Vec<&Rec>| c.get(1).map_or(f64::NAN, |r| r.f.re * 1e-9);
                say!("    f_even {:.6} GHz, f_odd {:.6} GHz, split {:.4} MHz (Q {:.2e}/{:.2e}, Si {:.3}/{:.3}, PML {:.1e}/{:.1e})", e.f.re * 1e-9, o.f.re * 1e-9, (o.f.re - e.f.re) * 1e-6, e.q, o.q, e.si, o.si, e.pml, o.pml);
                writeln!(
                    super_csv,
                    "{lam_um},{g_um},{qr},{q:.6e},{:.8},{:.8},{:.5},{:.3e},{:.3e},{:.4},{:.4},{:.3e},{:.3e},{:.4},{:.4},{:.6},{:.6},{},{},{:.1},{:.1},{}",
                    e.f.re * 1e-9,
                    o.f.re * 1e-9,
                    (o.f.re - e.f.re) * 1e-6,
                    e.q,
                    o.q,
                    e.si,
                    o.si,
                    e.pml,
                    o.pml,
                    e.core,
                    o.core,
                    f2(&ce),
                    f2(&co),
                    res[0].1.sol.dofs,
                    res[1].1.sol.dofs,
                    res[0].1.seconds,
                    res[1].1.seconds,
                    mesh.n_cells()
                )
                .unwrap();
            }
        }
        std::fs::write(out.join(format!("s2_modes_lam{lam_um}{tag}.csv")), &modes_csv).unwrap();
        std::fs::write(out.join(format!("s2_super_lam{lam_um}{tag}.csv")), &super_csv).unwrap();
        std::fs::write(out.join(format!("s2_log_lam{lam_um}{tag}.txt")), &log).unwrap();
    }

    // ---- full-domain check: both ridges, parity from u_z at probe points
    for &g_um in &full_gaps {
        let gap = g_um * UM;
        let (mesh, geo) = build(lam, gap, false, scale, side, depth);
        let core_half = geo.xe + 0.5 * lam;
        let outer = outer_labels(&mesh, &|p| (p[0].abs() - geo.xo).abs() < 1e-12 || (p[1] + geo.d + geo.pt).abs() < 1e-12);
        let mut spec = base_spec(psi, false);
        spec.bcs = vec![Bc { on: None, labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        let mut probes = vec![];
        for dx in [-0.35, 0.0, 0.35] {
            for z in [0.15 * UM, -0.1 * UM, -0.4 * UM] {
                let x = geo.xc + dx * lam;
                probes.push((nearest(&mesh, x, z), nearest(&mesh, -x, z)));
            }
        }
        say!("\n== full domain, lambda {lam_um} um, gap {g_um} um: {} cells", mesh.n_cells());
        let Some(s) = solve_sector(&mesh, &spec, q0, nm, v_film, shift, core_half, Some(&probes)) else {
            say!("  no ridge-bound mode found");
            continue;
        };
        say!("  {} dofs, {:.1} s, {} modes", s.sol.dofs, s.seconds, s.recs.len());
        for (rank, r) in s.recs.iter().enumerate() {
            writeln!(full_csv, "{lam_um},{g_um},{},{:.6},{:.4e},{:.1},{:.4},{:.4e},{:.4},{:.4e},{:.3},{},{},{:.1}", rank + 1, r.f.re * 1e-9, r.q, r.v, r.si, r.pml, r.core, r.ux_uz, r.parity, r.cand as u8, s.sol.dofs, s.seconds).unwrap();
            if r.cand && rank < 12 {
                say!("    {:.6} GHz Q {:.2e} Si {:.3} PML {:.1e} parity(u_z) {:+.3}", r.f.re * 1e-9, r.q, r.si, r.pml, r.parity);
            }
        }
        std::fs::write(out.join(format!("s2_full_lam{lam_um}{tag}.csv")), &full_csv).unwrap();
        std::fs::write(out.join(format!("s2_log_lam{lam_um}{tag}.txt")), &log).unwrap();
    }
}
