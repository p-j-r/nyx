//! Use case: Bragg gratings in a 300 nm SOI wire (500 nm wide, 80 nm deep teeth, 300 nm period), symmetric and asymmetric.
//!
//!   cargo run --release -p nyx-engine --example bragg_gratings                 # symmetric, then one-sided
//!   cargo run --release -p nyx-engine --example bragg_gratings -- symmetric    # one kind: symmetric | asymmetric | shifted
//!
//! Each grating is one 3D Floquet cell at the band edge (phase pi); the two lowest modes bound the stop band. The table is the
//! Stokes/Poincare analysis of those two modes with respect to (x lateral, y propagation, z normal); the Poincare sphere of all gratings is
//! written to `bragg_poincare.svg`.
use femcore::problem::{ProblemSpec, SolveResult, SolveState};
use femcore::usecases::bragg::*;

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
    let kinds: Vec<(&str, GratingKind)> = match std::env::args().nth(1).as_deref() {
        Some("symmetric") => vec![("symmetric", GratingKind::Symmetric)],
        Some("asymmetric") => vec![("asymmetric", GratingKind::Asymmetric)],
        Some("shifted") => vec![("shifted", GratingKind::HalfPeriodShifted)],
        Some("uniform") => vec![("uniform", GratingKind::Uniform)],
        _ => vec![("symmetric", GratingKind::Symmetric), ("asymmetric", GratingKind::Asymmetric)],
    };
    for (name, kind) in kinds {
        let width = std::env::var("WIDTH_NM").ok().and_then(|v| v.parse::<f64>().ok()).map(|w| w * 1e-9);
        let g = if std::env::var("COARSE").is_ok() {
            BraggGrating { kind, ..BraggGrating::coarse() }
        } else {
            BraggGrating { kind, ..Default::default() }
        };
        let g = BraggGrating { width: width.unwrap_or(g.width), ..g };
        // "square": the 300 nm wide wire, where TE0 and TM0 are almost degenerate, at its own matched period
        let g = if std::env::var("SQUARE").is_ok() {
            let g = BraggGrating { width: 300e-9, ..g };
            BraggGrating { period: g.matched_period(&solver)?, ..g }
        } else {
            g
        };
        let g =
            BraggGrating { order: std::env::var("ORDER").ok().and_then(|v| v.parse().ok()).unwrap_or(g.order), ..g };
        println!("== {} ==", kind.name());
        println!(
            "matched period for 1550 nm: {:.1} nm (used: {:.1} nm)",
            1e9 * g.matched_period(&solver)?,
            1e9 * g.period
        );
        let t0 = std::time::Instant::now();
        let r = g.solve(&solver)?;
        let analysis = LAST.with(|l| l.borrow_mut().take());
        println!("{} modes in {:.1} s, {} dofs", r.solutions.len(), t0.elapsed().as_secs_f64(), r.dofs);
        for l in &r.log {
            if !l.starts_with("    ") {
                println!("{l}");
            }
        }
        let b = g.band_edges_of(&r, &solver)?;
        println!(
            "stop band: {:.3} to {:.3} THz, gap {:.3} % ({:.1} nm wide at {:.1} nm), n_g = {:.3}",
            b.f_low * 1e-12,
            b.f_high * 1e-12,
            100.0 * b.gap_fraction,
            b.stopband * 1e9,
            b.bragg_wavelength * 1e9,
            b.group_index
        );
        println!(
            "coupling kappa = {:.0} 1/m ({:.2} 1/mm), coupled-mode estimate from the section indices {:.0} 1/m",
            b.kappa,
            b.kappa * 1e-3,
            b.kappa_estimate
        );
        if let Some(a) = &analysis {
            std::fs::write(format!("bragg_{name}_poincare.svg"), a.sphere_svg(360.0)).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
