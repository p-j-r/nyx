//! Stage S3 of the topological SAW study: arrays of identical 300 nm Si ridges (width = lambda, full etch, psi = 0) on r-plane
//! sapphire: the 8-ridge SSH chains (topological and trivial) and a 13-ridge uniform array. Every array is mirror symmetric about
//! x = 0, so the supermodes come from two half-domain solves per array (even sector: u_x = 0 on the mid-plane, odd sector:
//! u_y = u_z = 0), as in S2. For every supermode the program writes the frequency, Q, the elastic energy per ridge region and the
//! complex projection amplitudes c_k,m of the mode onto the single-ridge mode (S1 mode shape, translated to ridge m); the
//! effective coupled-ridge Hamiltonian and the propagation maps are built from these by `postproc_s3.py`.
//!   cargo run --release -p femcore --example saw_s3 -- [out_dir]
//! Environment: LAMBDA (um, 3), ARRAYS (comma list of topo,trivial,uniform13), MESH_SCALE, N_MODES (70), V_FILM (5790),
//! SHIFT (0.98), SIDE (5), DEPTH (5), GAP_STRONG (1.2), GAP_WEAK (7.2), GAP_UNIFORM (1.8), TAG.
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

struct Array {
    lam: f64,
    /// centres of the N ridges (full array), left to right
    xc: Vec<f64>,
    xe: f64,
    xo: f64,
    d: f64,
    pt: f64,
}

fn centres(lam: f64, gaps: &[f64]) -> Vec<f64> {
    let mut x = vec![0.0];
    for g in gaps {
        x.push(x.last().unwrap() + lam + g);
    }
    let mid = 0.5 * x.last().unwrap();
    x.iter().map(|v| v - mid).collect()
}

fn build_array(lam: f64, gaps: &[f64], half: bool, scale: f64, side: f64, depth_lam: f64) -> (SimplexMesh, Array) {
    let xc = centres(lam, gaps);
    let xe = xc.last().unwrap() + 0.5 * lam;
    let (xin, pt, d) = (xe + side * lam, lam, depth_lam * lam);
    let xo = xin + pt;
    let (h_far, h_surf, h_si, h_edge) = (lam / 8.0 * scale, lam / 16.0 * scale, 75e-9 * scale, 75e-9 * scale);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, dom: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: dom, h };
    let xlo = if half { 0.0 } else { -xin };
    let mut parts = vec![rect(xlo, -d, xin, 0.0, 1, None), rect(xlo, -0.6 * lam, xin, 0.0, 1, Some(h_surf))];
    for &x in xc.iter().filter(|&&x| !half || x > -1e-12) {
        let (a, b) = ((x - 0.5 * lam - 0.25 * UM).max(xlo), x + 0.5 * lam + 0.25 * UM);
        parts.push(rect(a, -0.25 * UM, b, 0.0, 1, Some(h_edge)));
    }
    parts.push(rect(xin, -d, xo, 0.0, 5, None));
    parts.push(rect(xlo, -d - pt, xin, -d, 6, None));
    parts.push(rect(xin, -d - pt, xo, -d, 9, None));
    if !half {
        parts.push(rect(-xo, -d, -xin, 0.0, 4, None));
        parts.push(rect(-xo, -d - pt, -xin, -d, 8, None));
    }
    for &x in xc.iter().filter(|&&x| !half || x > -1e-12) {
        parts.push(rect((x - 0.5 * lam).max(xlo), 0.0, x + 0.5 * lam, T_SI, 2, Some(h_si)));
    }
    let mesh = mesh_2d(&Spec2d { parts, h: h_far, min_angle_deg: 25.0 }).unwrap();
    (mesh, Array { lam, xc, xe, xo, d, pt })
}

fn base_spec(psi: f64, half: bool) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    set_materials(&mut s, psi, &[1, 2, 4, 5, 6, 8, 9], &[1, 4, 5, 6, 8, 9], &[2]);
    let mut dirs = vec![(5, PmlDir::Plus, PmlDir::Off), (6, PmlDir::Off, PmlDir::Minus), (9, PmlDir::Plus, PmlDir::Minus)];
    if !half {
        dirs.extend([(4, PmlDir::Minus, PmlDir::Off), (8, PmlDir::Minus, PmlDir::Minus)]);
    }
    for (dom, px, py) in dirs {
        s.materials.pml.insert(dom, pml(px, py));
    }
    s
}

fn cell_x(mesh: &SimplexMesh) -> Vec<f64> {
    (0..mesh.n_cells()).map(|c| mesh.cell(c).iter().map(|&v| mesh.point(v as usize)[0]).sum::<f64>() / 3.0).collect()
}

/// Lumped mass weight of every vertex (rho x area / 3 of the cells around it).
fn vertex_weights(mesh: &SimplexMesh) -> Vec<f64> {
    let mut w = vec![0.0; mesh.n_points()];
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let (a, b, d) = (mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize));
        let area = 0.5 * ((b[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (b[1] - a[1])).abs();
        let rho = if mesh.cell_labels[c] == 2 { RHO_SI } else { RHO_SAPPHIRE };
        for &i in v {
            w[i as usize] += rho * area / 3.0;
        }
    }
    w
}

/// P1 interpolation of vertex data on a mesh, with a grid hash for point location.
struct Locator<'a> {
    mesh: &'a SimplexMesh,
    h: f64,
    bins: std::collections::HashMap<(i64, i64), Vec<u32>>,
}

impl<'a> Locator<'a> {
    fn new(mesh: &'a SimplexMesh) -> Self {
        let h = 0.5 * UM;
        let mut bins: std::collections::HashMap<(i64, i64), Vec<u32>> = Default::default();
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let p: Vec<&[f64]> = v.iter().map(|&i| mesh.point(i as usize)).collect();
            let (x0, x1) = (p.iter().map(|q| q[0]).fold(f64::MAX, f64::min), p.iter().map(|q| q[0]).fold(f64::MIN, f64::max));
            let (y0, y1) = (p.iter().map(|q| q[1]).fold(f64::MAX, f64::min), p.iter().map(|q| q[1]).fold(f64::MIN, f64::max));
            for i in (x0 / h).floor() as i64..=(x1 / h).floor() as i64 {
                for j in (y0 / h).floor() as i64..=(y1 / h).floor() as i64 {
                    bins.entry((i, j)).or_default().push(c as u32);
                }
            }
        }
        Locator { mesh, h, bins }
    }

    /// (cell, barycentric coordinates) of the point (x, y) if it lies in the mesh.
    fn locate(&self, x: f64, y: f64) -> Option<(usize, [f64; 3])> {
        let key = ((x / self.h).floor() as i64, (y / self.h).floor() as i64);
        for &c in self.bins.get(&key)? {
            let v = self.mesh.cell(c as usize);
            let (a, b, d) = (self.mesh.point(v[0] as usize), self.mesh.point(v[1] as usize), self.mesh.point(v[2] as usize));
            let det = (b[1] - d[1]) * (a[0] - d[0]) + (d[0] - b[0]) * (a[1] - d[1]);
            let l0 = ((b[1] - d[1]) * (x - d[0]) + (d[0] - b[0]) * (y - d[1])) / det;
            let l1 = ((d[1] - a[1]) * (x - d[0]) + (a[0] - d[0]) * (y - d[1])) / det;
            let l2 = 1.0 - l0 - l1;
            let eps = -1e-9;
            if l0 >= eps && l1 >= eps && l2 >= eps {
                return Some((c as usize, [l0, l1, l2]));
            }
        }
        None
    }

    fn eval(&self, data: &[[C; 3]], x: f64, y: f64) -> [C; 3] {
        match self.locate(x, y) {
            Some((c, l)) => {
                let v = self.mesh.cell(c);
                let mut out = [C::new(0.0, 0.0); 3];
                for k in 0..3 {
                    for a in 0..3 {
                        out[k] += data[v[a] as usize][k] * l[a];
                    }
                }
                out
            }
            None => [C::new(0.0, 0.0); 3],
        }
    }
}

fn vertex_vectors(mesh: &SimplexMesh, m: &KzModes, j: usize) -> Vec<[C; 3]> {
    let n1 = mesh.n_dofs(2);
    let u = &m.vectors[j];
    (0..mesh.n_points()).map(|i| [u[i], u[n1 + i], u[2 * n1 + i]]).collect()
}

struct Picked {
    sector: &'static str,
    /// index into the modes of the solution
    j: usize,
    f: C,
    q: f64,
    si: f64,
    pml: f64,
    core: f64,
    /// energy weight in every ridge of the full array (partition of the x axis at the middle of the gaps)
    w: Vec<f64>,
    /// complex projection amplitudes on the single-ridge mode of every ridge of the full array
    c: Vec<C>,
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "runs/saw_topo".into()));
    std::fs::create_dir_all(&out).unwrap();
    let lam_um = envf("LAMBDA", 3.0);
    let lam = lam_um * UM;
    let (scale, nm, v_film, shift, side, depth) = (envf("MESH_SCALE", 1.0), envf("N_MODES", 70.0) as usize, envf("V_FILM", 5790.0), envf("SHIFT", 0.98), envf("SIDE", 5.0), envf("DEPTH", 5.0));
    let (g_s, g_w, g_u) = (envf("GAP_STRONG", 1.2) * UM, envf("GAP_WEAK", 7.2) * UM, envf("GAP_UNIFORM", 1.8) * UM);
    let psi = envf("PSI", 0.0);
    let q0 = 2.0 * std::f64::consts::PI / lam;
    let tag = std::env::var("TAG").unwrap_or_default();
    let names: Vec<String> = std::env::var("ARRAYS").ok().map(|s| s.split(',').map(|x| x.to_string()).collect()).unwrap_or_else(|| vec!["topo".into(), "trivial".into(), "uniform13".into()]);
    let mut log = String::new();
    macro_rules! say {
        ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); log.push_str(&s); log.push('\n'); }};
    }
    let mut timing = String::from("array,sector,dofs,cells,attempt,solve_s\n");

    // ---- the single-ridge mode (reference shape for the projections)
    let (mesh1, arr1) = build_array(lam, &[], false, scale, side, depth);
    let outer1 = outer_labels(&mesh1, &|p| (p[0].abs() - arr1.xo).abs() < 1e-12 || (p[1] + arr1.d + arr1.pt).abs() < 1e-12);
    let mut spec1 = base_spec(psi, false);
    spec1.bcs = vec![Bc { on: None, labels: outer1, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let f1 = 0.985 * v_film * q0 / (2.0 * std::f64::consts::PI);
    let sol1 = solve_q(&mesh1, &spec1, q0, 40, f1).unwrap();
    let cx1 = cell_x(&mesh1);
    let labels1 = mesh1.cell_labels.clone();
    let tot1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|_| true);
    let core1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|c| !PML_DOMAINS.contains(&labels1[c]) && cx1[c].abs() < 0.5 * lam + 0.5 * lam);
    let pml1 = elastic_energy_where(&mesh1, &sol1.modes, 2, &|c| PML_DOMAINS.contains(&labels1[c]));
    let mut order1: Vec<usize> = (0..sol1.modes.freqs.len()).collect();
    order1.sort_by(|&a, &b| sol1.modes.freqs[a].re.partial_cmp(&sol1.modes.freqs[b].re).unwrap());
    let j1 = *order1.iter().find(|&&j| core1[j] / tot1[j] >= CORE_MIN && pml1[j] / tot1[j] <= PML_MAX).expect("single ridge mode");
    let f_single = sol1.modes.freqs[j1].re;
    say!("single ridge (full domain, {} cells): f = {:.6} GHz, Q {:.2e}, {} dofs, {:.1} s", mesh1.n_cells(), f_single * 1e-9, q_factor(sol1.modes.freqs[j1]), sol1.dofs, sol1.seconds);
    writeln!(timing, "single,full,{},{},1,{:.2}", sol1.dofs, mesh1.n_cells(), sol1.seconds).unwrap();
    let phi: Vec<[C; 3]> = {
        let v = vertex_vectors(&mesh1, &sol1.modes, j1);
        let w = vertex_weights(&mesh1);
        let norm = v.iter().zip(&w).map(|(a, w)| w * (a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr())).sum::<f64>().sqrt();
        v.into_iter().map(|a| [a[0] / norm, a[1] / norm, a[2] / norm]).collect()
    };
    let loc1 = Locator::new(&mesh1);

    for name in &names {
        let gaps: Vec<f64> = match name.as_str() {
            "topo" => (0..7).map(|i| if i % 2 == 0 { g_w } else { g_s }).collect(),
            "trivial" => (0..7).map(|i| if i % 2 == 0 { g_s } else { g_w }).collect(),
            "uniform13" => vec![g_u; 12],
            other => panic!("unknown array {other}"),
        };
        let n = gaps.len() + 1;
        let (mesh, arr) = build_array(lam, &gaps, true, scale, side, depth);
        let core_half = arr.xe + 0.5 * lam;
        let sym = outer_labels(&mesh, &|p| p[0].abs() < 1e-12);
        let outer = outer_labels(&mesh, &|p| (p[0] - arr.xo).abs() < 1e-12 || (p[1] + arr.d + arr.pt).abs() < 1e-12);
        // ridges of the half domain (centre >= 0), partition boundaries at the middle of the gaps
        let half_ridges: Vec<usize> = (0..n).filter(|&j| arr.xc[j] > -1e-12).collect();
        let bound = |j: usize| -> f64 { if j + 1 >= n { f64::INFINITY } else { 0.5 * (arr.xc[j] + arr.xc[j + 1]) } };
        let lower = |j: usize| -> f64 { if j == 0 || arr.xc[j] < 1e-12 { f64::NEG_INFINITY } else { 0.5 * (arr.xc[j - 1] + arr.xc[j]) } };
        let cx = cell_x(&mesh);
        let labels = mesh.cell_labels.clone();
        let wv = vertex_weights(&mesh);
        let n_even = (n + 1) / 2;
        let n_odd = n / 2;
        say!("\n== {name}: {n} ridges, gaps {:?} um; half domain {} cells, half width {:.1} um", gaps.iter().map(|g| g / UM).collect::<Vec<_>>(), mesh.n_cells(), arr.xo / UM);
        // projection reference: per half-ridge region, the vertices and the single-ridge mode there
        let mut proj_ref: Vec<(Vec<usize>, Vec<[C; 3]>, f64)> = vec![];
        for &j in &half_ridges {
            let (lo, hi) = (lower(j), bound(j));
            let idx: Vec<usize> = (0..mesh.n_points()).filter(|&i| mesh.point(i)[0] >= lo.max(0.0) - 1e-12 && mesh.point(i)[0] < hi).collect();
            let vals: Vec<[C; 3]> = idx.iter().map(|&i| loc1.eval(&phi, mesh.point(i)[0] - arr.xc[j], mesh.point(i)[1])).collect();
            let nphi: f64 = idx.iter().zip(&vals).map(|(&i, a)| wv[i] * (a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr())).sum();
            proj_ref.push((idx, vals, nphi));
        }
        let mut picked: Vec<Picked> = vec![];
        let mut dump: Vec<(String, KzModes, usize)> = vec![];
        for (sector, fixed, expected) in [("even", [true, false, false], n_even), ("odd", [false, true, true], n_odd)] {
            let mut spec = base_spec(psi, true);
            spec.bcs = vec![
                Bc { on: None, labels: outer.clone(), kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
                Bc { on: None, labels: sym.clone(), kind: BcKind::Fixed { fixed, values: [0.0; 3] } },
            ];
            let mut done = false;
            for (attempt, (sh, nmodes)) in [(shift, nm), (0.94, nm + 30), (1.0, nm + 30)].into_iter().enumerate() {
                let f_search = sh * v_film * q0 / (2.0 * std::f64::consts::PI);
                let sol = match solve_q(&mesh, &spec, q0, nmodes, f_search) {
                    Ok(s) => s,
                    Err(e) => {
                        say!("  {sector}: solve failed: {e}");
                        continue;
                    }
                };
                say!("  {sector} (attempt {}, shift {:.3} GHz): {} dofs, {:.1} s, {} modes", attempt + 1, f_search * 1e-9, sol.dofs, sol.seconds, sol.modes.freqs.len());
                writeln!(timing, "{name},{sector},{},{},{},{:.2}", sol.dofs, mesh.n_cells(), attempt + 1, sol.seconds).unwrap();
                let m = &sol.modes;
                let tot = elastic_energy_where(&mesh, m, 2, &|_| true);
                let pmlf = elastic_energy_where(&mesh, m, 2, &|c| PML_DOMAINS.contains(&labels[c]));
                let core = elastic_energy_where(&mesh, m, 2, &|c| !PML_DOMAINS.contains(&labels[c]) && cx[c].abs() < core_half);
                let mut ord: Vec<usize> = (0..m.freqs.len()).collect();
                ord.sort_by(|&a, &b| m.freqs[a].re.partial_cmp(&m.freqs[b].re).unwrap());
                let cands: Vec<usize> = ord.iter().copied().filter(|&j| core[j] / tot[j] >= CORE_MIN && pmlf[j] / tot[j] <= PML_MAX).collect();
                if cands.len() < expected {
                    say!("    only {} ridge-bound modes found (expected {expected}), retrying", cands.len());
                    continue;
                }
                if cands.len() > expected {
                    say!("    {} candidates, keeping the lowest {expected}", cands.len());
                }
                let sel: Vec<usize> = cands.into_iter().take(expected).collect();
                let sif = elastic_energy_where(&mesh, m, 2, &|c| labels[c] == 2);
                // energy per half-ridge region
                let e_reg: Vec<Vec<f64>> = half_ridges
                    .iter()
                    .map(|&j| {
                        let (lo, hi) = (lower(j), bound(j));
                        elastic_energy_where(&mesh, m, 2, &|c| cx[c] >= lo.max(0.0) && cx[c] < hi)
                    })
                    .collect();
                for &j in &sel {
                    let u = vertex_vectors(&mesh, m, j);
                    let u2: f64 = (0..mesh.n_points()).map(|i| wv[i] * (u[i][0].norm_sqr() + u[i][1].norm_sqr() + u[i][2].norm_sqr())).sum();
                    let mut w = vec![0.0; n];
                    let mut c = vec![C::new(0.0, 0.0); n];
                    for (h, &r) in half_ridges.iter().enumerate() {
                        let centre = arr.xc[r].abs() < 1e-9;
                        let mirror = n - 1 - r;
                        let e = e_reg[h][j] / tot[j];
                        let (idx, vals, nphi) = &proj_ref[h];
                        let mut num = C::new(0.0, 0.0);
                        for (&i, a) in idx.iter().zip(vals) {
                            num += (a[0].conj() * u[i][0] + a[1].conj() * u[i][1] + a[2].conj() * u[i][2]) * wv[i];
                        }
                        let ch = num / (nphi.sqrt() * u2.sqrt());
                        if centre {
                            w[r] = e;
                            c[r] = if sector == "even" { ch } else { C::new(0.0, 0.0) };
                        } else {
                            w[r] = 0.5 * e;
                            w[mirror] = 0.5 * e;
                            let s = if sector == "even" { 1.0 } else { -1.0 };
                            c[r] = ch / 2.0f64.sqrt();
                            c[mirror] = ch / 2.0f64.sqrt() * s;
                        }
                    }
                    picked.push(Picked { sector: if sector == "even" { "even" } else { "odd" }, j: 0, f: m.freqs[j], q: q_factor(m.freqs[j]), si: sif[j] / tot[j], pml: pmlf[j] / tot[j], core: core[j] / tot[j], w, c });
                    picked.last_mut().unwrap().j = dump.len();
                    dump.push((sector.to_string(), KzModes { kz: m.kz, freqs: vec![m.freqs[j]], vectors: vec![m.vectors[j].clone()], c6: vec![], rho: vec![], converged: true }, 0));
                }
                done = true;
                break;
            }
            if !done {
                say!("  {sector}: ridge-bound modes not found");
            }
        }
        // ---- write
        let mut order: Vec<usize> = (0..picked.len()).collect();
        order.sort_by(|&a, &b| picked[a].f.re.partial_cmp(&picked[b].f.re).unwrap());
        let mut modes_csv = String::from("rank,sector,f_GHz,f_im_GHz,df_from_single_MHz,Q,v_m_per_s,si_frac,pml_frac,core_frac,end_ridge_weight,proj_end_weight");
        for m in 1..=n {
            write!(modes_csv, ",w{m}").unwrap();
        }
        modes_csv.push('\n');
        let mut proj_csv = String::from("rank,sector,f_GHz");
        for m in 1..=n {
            write!(proj_csv, ",cre{m},cim{m}").unwrap();
        }
        proj_csv.push('\n');
        say!("  rank sector  f (GHz)    df (MHz)   Q        Si    PML     end-ridge energy weight  projected end weight");
        for (rank, &k) in order.iter().enumerate() {
            let p = &picked[k];
            let end = p.w[0] + p.w[n - 1];
            let csum: f64 = p.c.iter().map(|z| z.norm_sqr()).sum();
            let pend = (p.c[0].norm_sqr() + p.c[n - 1].norm_sqr()) / csum;
            say!("  {:>3}  {:<5} {:.6}  {:+8.3}  {:.1e}  {:.3}  {:.0e}   {:.3}   {:.3}", rank + 1, p.sector, p.f.re * 1e-9, (p.f.re - f_single) * 1e-6, p.q, p.si, p.pml, end, pend);
            write!(modes_csv, "{},{},{:.8},{:.4e},{:.5},{:.3e},{:.1},{:.4},{:.3e},{:.4},{:.5},{:.5}", rank + 1, p.sector, p.f.re * 1e-9, p.f.im * 1e-9, (p.f.re - f_single) * 1e-6, p.q, 2.0 * std::f64::consts::PI * p.f.re / q0, p.si, p.pml, p.core, end, pend).unwrap();
            for w in &p.w {
                write!(modes_csv, ",{:.5}", w).unwrap();
            }
            modes_csv.push('\n');
            write!(proj_csv, "{},{},{:.8}", rank + 1, p.sector, p.f.re * 1e-9).unwrap();
            for z in &p.c {
                write!(proj_csv, ",{:.6e},{:.6e}", z.re, z.im).unwrap();
            }
            proj_csv.push('\n');
        }
        std::fs::write(out.join(format!("s3_modes_{name}{tag}.csv")), modes_csv).unwrap();
        std::fs::write(out.join(format!("s3_proj_{name}{tag}.csv")), proj_csv).unwrap();
        // mode shapes of the modes with the largest end-ridge weight (one per sector) and of the mode closest to f_single
        let mut want: Vec<(String, usize)> = vec![];
        for sec in ["even", "odd"] {
            if let Some(&k) = order.iter().filter(|&&k| picked[k].sector == sec).max_by(|&&a, &&b| (picked[a].w[0] + picked[a].w[n - 1]).partial_cmp(&(picked[b].w[0] + picked[b].w[n - 1])).unwrap()) {
                want.push((format!("{sec}_maxend"), k));
            }
            if let Some(&k) = order.iter().filter(|&&k| picked[k].sector == sec).min_by(|&&a, &&b| (picked[a].f.re - f_single).abs().partial_cmp(&(picked[b].f.re - f_single).abs()).unwrap()) {
                want.push((format!("{sec}_nearest"), k));
            }
        }
        let n1 = mesh.n_dofs(2);
        for (label, k) in want {
            let u = &dump[picked[k].j].1.vectors[0];
            let big = (0..mesh.n_points()).max_by(|&a, &b| u[n1 + a].norm().partial_cmp(&u[n1 + b].norm()).unwrap()).unwrap();
            let ph = u[n1 + big].conj() / u[n1 + big].norm();
            let mut pts = String::from("x_um,z_um,uz_signed,uz_abs\n");
            for i in 0..mesh.n_points() {
                let p = mesh.point(i);
                writeln!(pts, "{:.5},{:.5},{:.6e},{:.6e}", p[0] / UM, p[1] / UM, (u[n1 + i] * ph).re, u[n1 + i].norm()).unwrap();
            }
            std::fs::write(out.join(format!("s3_shape_{name}_{label}_pts{tag}.csv")), pts).unwrap();
            say!("  mode shape {label}: {:.6} GHz ({})", picked[k].f.re * 1e-9, picked[k].sector);
        }
        let mut tri = String::from("a,b,c,domain\n");
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            writeln!(tri, "{},{},{},{}", v[0], v[1], v[2], labels[c]).unwrap();
        }
        std::fs::write(out.join(format!("s3_shape_{name}_tri{tag}.csv")), tri).unwrap();
        let mut info = String::new();
        writeln!(info, "n={n}\nf_single_GHz={:.8}\nridge_centres_um={}", f_single * 1e-9, arr.xc.iter().map(|x| format!("{:.4}", x / UM)).collect::<Vec<_>>().join(",")).unwrap();
        std::fs::write(out.join(format!("s3_info_{name}{tag}.txt")), info).unwrap();
        std::fs::write(out.join(format!("s3_timing{tag}.csv")), &timing).unwrap();
        std::fs::write(out.join(format!("s3_log{tag}.txt")), &log).unwrap();
    }
}
