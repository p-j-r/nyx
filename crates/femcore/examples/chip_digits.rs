//! Full-precision numbers of the Brillouin chip study, for before/after comparisons of solver changes:
//! the optical effective index and, for every elastic mode with a noticeable gain, frequency, Q and the
//! photoelastic, moving-boundary and total gain with 12 significant digits.
//!
//!   cargo run --release -p femcore --example chip_digits -- 300TE 0 > before.txt
//!   (chip name from `chips::CHIPS`, crystal azimuth psi in degrees)
use femcore::brillouin::*;
use femcore::chips::*;
use femcore::problem::*;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "300TE".into());
    let psi: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let chip = chip_by_name(&name).expect("unknown chip name");
    let f_search = match name.as_str() {
        "300TE" => 18.35e9,
        "500TE" => 20.7e9,
        "300TM" => 15.6e9,
        _ => 19.0e9,
    };
    let t = std::time::Instant::now();
    let scale: f64 = std::env::var("MESH_SCALE").ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
    let (mesh, spec) = chip_model(&chip, psi, scale, 40, f_search).expect("chip model");
    let inner = spec.clone();
    let study = study_straight(&mesh, &inner, &SolveState::default()).expect("study");
    eprintln!("{name} psi {psi}: {:.1} s", t.elapsed().as_secs_f64());
    println!("chip {name} psi {psi}");
    println!("neff {:.12e} {:.12e}", study.neff[0].re, study.neff[0].im);
    let pr = &study.pairs[0];
    println!("q {:.12e}  group_index {:.12e} {:.12e}", pr.q, pr.group_index.0, pr.group_index.1);
    let mut rows: Vec<usize> = (0..pr.rows.len()).filter(|&i| pr.rows[i].gain_pe.max(pr.rows[i].gain_mb) > 1e-3).collect();
    rows.sort_by(|&a, &b| pr.rows[a].freq.re.partial_cmp(&pr.rows[b].freq.re).unwrap());
    for i in rows {
        let r = &pr.rows[i];
        println!("f {:.12e} Q {:.6e} MB {:.12e} PE {:.12e} total {:.12e}", r.freq.re, r.q_eig, r.gain_mb, r.gain_pe, r.gain_total);
    }
}
