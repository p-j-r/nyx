//! Headless runs without a window, for workstations and clusters:
//!
//! ```text
//! nyx --batch model.nyx.xml [--physics es|heat|solid|emw|brillouin] [--output serial|parallel FILE.h5]
//! mpiexec -n 16 nyx --batch model.nyx.xml --output parallel results.h5     (build with --features phdf5)
//! ```
//!
//! The project (geometry or model file, materials, boundary conditions, study) is read, meshed, solved
//! with the same code as the window, and the results go to the file chosen in the project or on the
//! command line. In a parallel run every process solves the whole problem and writes its own block
//! of rows of the results file (the solvers themselves are not distributed).

use crate::state::*;
use femcore::problem::*;

pub fn run(args: &[String]) -> Result<(), String> {
    let mut path = None;
    let mut physics = None;
    let mut output: Option<(OutputMode, String)> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--physics" => {
                i += 1;
                physics = Some(match args.get(i).map(|s| s.as_str()) {
                    Some("es") => Physics::Electrostatics,
                    Some("heat") => Physics::Heat,
                    Some("solid") => Physics::Solid,
                    Some("emw") => Physics::Emw,
                    Some("brillouin") => Physics::Brillouin,
                    other => return Err(format!("--physics needs es, heat, solid, emw or brillouin, got {other:?}")),
                });
            }
            "--output" => {
                let mode = match args.get(i + 1).map(|s| s.as_str()) {
                    Some("serial") => OutputMode::Hdf5Serial,
                    Some("parallel") => OutputMode::Hdf5Parallel,
                    other => return Err(format!("--output needs serial or parallel and a file, got {other:?}")),
                };
                let file = args.get(i + 2).ok_or("--output needs a file name after the mode")?.clone();
                output = Some((mode, file));
                i += 2;
            }
            a if a.starts_with("--") => return Err(format!("unknown option {a}")),
            a => path = Some(a.to_string()),
        }
        i += 1;
    }
    let path = path.ok_or("usage: nyx --batch PROJECT.nyx.xml [--physics es|heat|solid|emw|brillouin] [--output serial|parallel FILE.h5]")?;
    let proj = Project::load(std::path::Path::new(&path))?;
    let phys = physics.unwrap_or(proj.active);
    let mut spec = proj.spec(phys).clone();
    if let Some((m, f)) = output {
        spec.output.mode = m;
        spec.output.path = f;
    }
    #[cfg(feature = "hdf5")]
    let (rank, size) = match spec.output.mode {
        OutputMode::Hdf5Parallel => femcore::h5io::start_runtime(femcore::h5io::H5Mode::Parallel)?,
        _ => (0, 1),
    };
    #[cfg(not(feature = "hdf5"))]
    let (rank, size) = (0usize, 1usize);
    let say = |s: String| {
        if rank == 0 {
            println!("{s}");
        }
    };
    let mesh = std::sync::Arc::new(proj.geometry.build_for(&spec)?);
    say(format!("{}: {} cells, {} points, {} process(es)", phys.name(), mesh.n_cells(), mesh.n_points(), size));
    let mut state = SolveState::default();
    let result = femcore::solve::solve(&mesh, &spec, &mut state)?;
    say(format!("solved in {:.2} s, {} dofs, {} solution(s)", result.seconds, result.dofs, result.solutions.len()));
    for (n, v) in &result.scalars {
        match v {
            Scalar::Num(x) => say(format!("  {n} = {x}")),
            Scalar::Text(t) => say(format!("  {n} = {t}")),
        }
    }
    if !result.log.is_empty() {
        say("what the solver did:".into());
        for l in &result.log {
            say(format!("  {l}"));
        }
    }
    match crate::worker::write_output(&mesh, &spec, &result) {
        Some(Ok(m)) => say(m),
        Some(Err(e)) => return Err(format!("results file: {e}")),
        None => say("no results file requested (use --output or set it in the project)".into()),
    }
    Ok(())
}
