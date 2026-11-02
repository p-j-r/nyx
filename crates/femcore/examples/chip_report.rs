//! Backward Brillouin gain of the four SOS ring chips computed in Rust with the full anisotropic tensors
//! (stiffness, photoelastic and permittivity of silicon and tilted r-plane sapphire), the phase-matched overlap integrals
//! and the exact moving-boundary form. Writes a markdown report with the photoelastic (PE) and moving-boundary (MB) gains
//! of every strong elastic mode, next to the as-is photoelastic form (strain without complex conjugate).
//!
//!   cargo run --release -p femcore --example chip_report -- [out.md] [chips, e.g. 300TE,500TM] [psi list, e.g. 0,45]
use femcore::brillouin::*;
use femcore::chips::*;
use femcore::problem::*;
use std::fmt::Write as _;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let out_path = args.first().cloned().unwrap_or_else(|| "chip_report.md".into());
    let chips: Vec<ChipSpec> = args.get(1).map(|s| s.split(',').filter_map(chip_by_name).collect()).unwrap_or_else(|| CHIPS.to_vec());
    let psis: Vec<f64> = args.get(2).map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![0.0, 45.0]);
    let f_search = |name: &str| match name {
        "300TE" => 18.35e9,
        "500TE" => 20.7e9,
        "300TM" => 15.6e9,
        _ => 19.0e9,
    };
    let mut md = String::new();
    writeln!(md, "# Brillouin gain of the four ring chips, rebuilt in Rust\n").unwrap();
    writeln!(md, "Backward scattering of the design mode with itself at 1550 nm; rib on slab on r-plane sapphire (c axis 57.6 degrees from the surface normal, tilted towards the propagation direction at azimuth psi = 0); wafer frame (x lateral, y propagation, z normal). Silicon c11 = 165.7 GPa, optical n_Si = {N_SILICON}, sapphire uniaxial n_o = {N_SAPPHIRE_ORDINARY}, n_e = {N_SAPPHIRE_EXTRAORDINARY} with the tilted crystal tensor (exact coupled mode solver). Gains are in 1/(W m) with the mechanical Q capped at 1000, power and group velocity from the mode fields.\n").unwrap();
    writeln!(md, "Columns: PE and MB are the photoelastic and moving-boundary gains alone (`|f_PE|`, `|f_MB|` in the gain formula); total combines both overlap integrals with their phase. *PE, as is* and *total, as is* use the strain without complex conjugate.\n").unwrap();
    for c in &chips {
        for &psi in &psis {
            let t = std::time::Instant::now();
            let (mesh, spec) = match chip_model(c, psi, 1.0, 40, f_search(c.name)) {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("{} psi {psi}: {e}", c.name);
                    continue;
                }
            };
            let inner = spec.clone();
            let study = match study_straight(&mesh, &inner, &SolveState::default()) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("{} psi {psi}: {e}", c.name);
                    continue;
                }
            };
            let neff = study.neff[0].re;
            let pr = &study.pairs[0];
            eprintln!("{} psi {psi}: neff {neff:.5}, q {:.3} /um, {} elastic modes, {:.0} s", c.name, pr.q * 1e-6, pr.rows.len(), t.elapsed().as_secs_f64());
            writeln!(md, "## {} at psi = {psi} degrees\n", c.name).unwrap();
            writeln!(md, "Design mode neff = {neff:.5} (group index {:.3}), q = {:.3} /um, {} cells, {} elastic modes, {:.0} s.\n", pr.group_index.0, pr.q * 1e-6, mesh.n_cells(), pr.rows.len(), t.elapsed().as_secs_f64()).unwrap();
            writeln!(md, "| f (GHz) | Q | body / PML | MB | PE | total | PE, as is | total, as is |").unwrap();
            writeln!(md, "|---|---|---|---|---|---|---|---|").unwrap();
            let mut idx: Vec<usize> = (0..pr.rows.len()).filter(|&i| pr.rows[i].frac_body > 0.3).collect();
            idx.sort_by(|&a, &b| pr.rows[b].gain_total.partial_cmp(&pr.rows[a].gain_total).unwrap().then(pr.rows[b].gain_pe.partial_cmp(&pr.rows[a].gain_pe).unwrap()));
            let mut shown = 0;
            for &i in &idx {
                let (r, ra) = (&pr.rows[i], &pr.rows_alt[i]);
                if r.gain_pe.max(r.gain_mb).max(ra.gain_pe) < 0.05 && shown >= 6 {
                    continue;
                }
                let f = r.freq.re * 1e-9;
                let q = if r.q_eig.is_finite() { format!("{:.2e}", r.q_eig) } else { "inf".into() };
                writeln!(md, "| {f:.3} | {q} | {:.2} / {:.2} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} |", r.frac_body, r.frac_pml, r.gain_mb, r.gain_pe, r.gain_total, ra.gain_pe, ra.gain_total).unwrap();
                shown += 1;
                if shown >= 16 {
                    break;
                }
            }
            writeln!(md).unwrap();
            let _ = std::fs::write(&out_path, &md);
        }
    }
    let _ = std::fs::write(&out_path, &md);
    println!("{md}");
}
