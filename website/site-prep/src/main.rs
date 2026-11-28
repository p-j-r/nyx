mod api;
mod benchmarks;
mod examples;
mod materials;
mod sanitize;
mod tokens;

use std::env;
use std::path::PathBuf;
use std::process;

struct CliArgs {
    repo: PathBuf,
    out: Option<PathBuf>,
    check: bool,
}

fn parse_args() -> Result<CliArgs, String> {
    let mut repo = PathBuf::from(".");
    let mut out = None;
    let mut check = false;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--repo" => {
                let val = args
                    .next()
                    .ok_or_else(|| "--repo requires a path argument".to_string())?;
                repo = PathBuf::from(val);
            }
            "--out" => {
                let val = args
                    .next()
                    .ok_or_else(|| "--out requires a path argument".to_string())?;
                out = Some(PathBuf::from(val));
            }
            "--check" => {
                check = true;
            }
            "-h" | "--help" => {
                println!(
                    "Usage: site-prep [--repo <path>] [--out <website_path>] [--check]\n\n\
                     Generates website/sass/_tokens.scss and website/data/*.json from repository\n\
                     sources, enforcing the public confidentiality deny-list."
                );
                process::exit(0);
            }
            other => {
                return Err(format!("unknown argument '{other}' (use --help)"));
            }
        }
    }

    Ok(CliArgs { repo, out, check })
}

fn run() -> Result<(), String> {
    let args = parse_args()?;
    let out_site = args.out.unwrap_or_else(|| args.repo.join("website"));

    let role_count = tokens::generate_tokens_scss(&args.repo, &out_site)?;
    let mat_count = materials::generate_materials_json(&args.repo, &out_site)?;
    let ex_count = examples::generate_examples_json(&args.repo, &out_site)?;
    let (suite_count, target_count) =
        benchmarks::generate_benchmarks_json(&args.repo, &out_site)?;
    let (cmd_count, catalog_hash) = api::generate_api_json(&args.repo, &out_site)?;

    // Always scan generated output files; when --check is passed, also scan all authored
    // site files under `out_site` (excluding `public` and `target`).
    let scanned = if args.check {
        sanitize::scan_website_tree(&out_site)?
    } else {
        let mut total = sanitize::scan_website_tree(&out_site.join("data"))?;
        total += sanitize::scan_website_tree(&out_site.join("sass"))?;
        total
    };

    println!(
        "site-prep: ok | roles={role_count}, materials={mat_count}, examples={ex_count}, \
         benchmark_suites={suite_count} (targets={target_count}), commands={cmd_count} \
         ({catalog_hash}), scanned_files={scanned}, confidentiality_hits=0"
    );
    Ok(())
}

fn main() {
    if let Err(err) = run() {
        eprintln!("site-prep error: {err}");
        process::exit(1);
    }
}
