//! Use case: modes of a wire on 300 nm SOI at 1550 nm with the polarization analysis (Stokes parameters, Poincare sphere).
//!
//!   cargo run --release -p nyx-engine --example soi_modes                 # one wire (500 nm), table of the modes
//!   cargo run --release -p nyx-engine --example soi_modes -- sweep         # width sweep 300 ... 1000 nm
//!   cargo run --release -p nyx-engine --example soi_modes -- hybrid        # 300 nm wire, left wall sloped to 75 degrees: hybrid modes
//!
//! The Poincare sphere of the single wire is written to `soi_poincare.svg` and the table to `soi_modes.csv` in the current folder.
use femcore::problem::{ProblemSpec, SolveResult, SolveState};
use femcore::usecases::soi::*;

thread_local! {
    /// the polarization analysis of the last solve (the use cases hand out the result only)
    static LAST: std::cell::RefCell<Option<nyx_engine::analysis::ModeAnalysis>> = const { std::cell::RefCell::new(None) };
}

/// What the use cases solve with: the bridge to the physics crates, keeping the analysis of the solve in `LAST`.
fn solver(mesh: &femcore::mesh::SimplexMesh, spec: &ProblemSpec) -> Result<SolveResult, String> {
    let p = femcore::progress::current().unwrap_or_else(femcore::progress::Progress::new);
    let (r, a) = nyx_engine::legacy::solve_with_analysis(mesh, spec, &mut SolveState::default(), &p)?;
    LAST.with(|l| *l.borrow_mut() = a);
    Ok(r)
}

fn main() -> Result<(), String> {
    let arg = std::env::args().nth(1).unwrap_or_default();
    match arg.as_str() {
        "sweep" => {
            println!(
                "{:>8}  {:>9} {:>9} {:>9}  {:>7} {:>7}",
                "w (nm)", "TE0 neff", "TM0 neff", "TE1 neff", "TE0 TE", "TM0 TE"
            );
            for w in [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0] {
                let wire = SoiWire { width: w * 1e-9, n_modes: if w >= 700.0 { 4 } else { 2 }, ..Default::default() };
                let r = wire.solve(&solver)?;
                let mut rows: Vec<_> =
                    r.solutions.iter().map(|s| (s.neff.unwrap().0, s.te_fraction.unwrap_or(f64::NAN))).collect();
                rows.sort_by(|x, y| y.0.total_cmp(&x.0));
                let te: Vec<_> = rows.iter().filter(|x| x.1 > 0.5).collect();
                let tm: Vec<_> = rows.iter().filter(|x| x.1 <= 0.5).collect();
                let get = |v: &Vec<&(f64, f64)>, k: usize| v.get(k).map(|x| x.0).unwrap_or(f64::NAN);
                println!(
                    "{:>8.0}  {:>9.5} {:>9.5} {:>9.5}  {:>7.3} {:>7.3}",
                    w,
                    get(&te, 0),
                    get(&tm, 0),
                    get(&te, 1),
                    te.first().map(|x| x.1).unwrap_or(f64::NAN),
                    tm.first().map(|x| x.1).unwrap_or(f64::NAN)
                );
            }
        }
        other => {
            let n_modes = std::env::var("SOI_MODES").ok().and_then(|v| v.parse().ok()).unwrap_or(2);
            let wire = if other == "hybrid" {
                SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), n_modes, ..Default::default() }
            } else {
                SoiWire { n_modes, ..Default::default() }
            };
            let t0 = std::time::Instant::now();
            let r = wire.solve(&solver)?;
            println!("{} modes in {:.1} s ({} dofs)", r.solutions.len(), t0.elapsed().as_secs_f64(), r.dofs);
            for l in &r.log {
                println!("{l}");
            }
            let a = LAST.with(|l| l.borrow_mut().take()).ok_or("no analysis")?;
            std::fs::write("soi_poincare.svg", a.sphere_svg(360.0)).map_err(|e| e.to_string())?;
            std::fs::write("soi_modes.csv", a.csv()).map_err(|e| e.to_string())?;
            for pol in [SlabPol::Te, SlabPol::Tm] {
                println!("slab limit {:?}: {:?}", pol, wire.slab_reference(pol));
            }
        }
    }
    Ok(())
}
