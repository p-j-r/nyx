// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-Rust SQuADDS cQED benchmark runner, calibration driver, Richardson GCI & thread-determinism harness,
//! and anti-survivorship statistical gate evaluator (`Plan §4.4`, `§4.5`, `§5`, `§7.1`, `ADR 13`).
//!
//! Subcommands:
//!   squadds_bench calibrate --which <C0.1|C0.2|C0.3|all> --out <dir>
//!   squadds_bench single --layout <path.layout.json> [--mode <spec>] [--tier <t>] [--scale-h <f>] [--threads <N>] [--measure-peak-ram] [--append] --out <path>
//!   squadds_bench run --manifest <path> --inputs <dir> --tier <smoke|p1|p2|core|ablation|qcap|airbridge> [--suite <id>] [--mode <spec>] [--scale-h <f>] [--threads <N>] [--limit <K>] [--measure-peak-ram] [--append] --out <runs.jsonl>
//!   squadds_bench gci --layout <path.layout.json> [--mode <spec>] [--r <ratio>] [--levels <3|4>] [--threads <N>] [--append] --out <gci.jsonl>
//!   squadds_bench determinism --layout <path.layout.json> [--mode <spec>] [--threads-a 1] [--threads-b 8] [--append] --out <det.jsonl>
//!   squadds_bench gate --runs <runs.jsonl> --gate <G0|G1|G2|G3|all> --out <gate_report.json>
//!   squadds_bench bug-checks [--inputs <dir>] [--append] --out <bug_checks.jsonl>

#[path = "squadds_bench/mod.rs"]
mod squadds_bench;

use serde_json::json;
use squadds_bench::bug_checks::run_all_bug_checks;
use squadds_bench::calibrate::run_calibrations;
use squadds_bench::electro::run_electro_case;
use squadds_bench::em_hw::run_em_hw_case;
use squadds_bench::gate::evaluate_gates;
use squadds_bench::math::compute_roache_gci;
use squadds_bench::schema::{load_manifest_items, parse_json_lenient_nan, LayoutSpec, RunRecord};
use std::io::Write;
use std::path::{Path, PathBuf};

fn is_electro_dispatch(layout: &LayoutSpec, mode: &str) -> bool {
    let m = mode.trim().to_ascii_lowercase();
    if matches!(
        m.as_str(),
        "g2_dual_box"
            | "g2_dual_box_raw"
            | "far_field"
            | "open"
            | "amr_ps"
            | "qp_graded"
            | "thick_conformal"
            | "planar_bem"
            | "gds_import"
    ) {
        return true;
    }
    if matches!(
        m.as_str(),
        "cavity_qep"
            | "driven_s21"
            | "driven"
            | "driven_cap"
            | "driven_filter"
            | "qep_and_driven"
            | "hybrid"
            | "qc_all"
            | "m1"
            | "m2"
            | "m3"
            | "m4"
            | "w0"
            | "w1"
            | "w2"
            | "w3"
            | "muninn"
            | "dissipator_squid"
            | "airbridge_qep"
    ) {
        return false;
    }
    let s = layout.suite_id.to_ascii_lowercase();
    s.starts_with("s1.") || s.starts_with("ref") || s.contains("qcap") || s.contains("cap_matrix")
}

fn dispatch_single_in_process(
    layout_path: &Path,
    tier_override: Option<&str>,
    mode: &str,
    scale_h: f64,
    cluster_override: Option<&str>,
    weight_override: Option<f64>,
) -> Result<RunRecord, String> {
    let txt = std::fs::read_to_string(layout_path)
        .map_err(|e| format!("failed to read layout '{}': {e}", layout_path.display()))?;
    let mut layout: LayoutSpec =
        parse_json_lenient_nan(&txt).map_err(|e| format!("failed to parse layout '{}': {e}", layout_path.display()))?;
    if let Some(c) = cluster_override {
        layout.cluster_id = Some(c.to_string());
    }
    if let Some(w) = weight_override {
        layout.stratum_weight = Some(w);
    }
    let tier = tier_override.map(String::from).or_else(|| layout.tier.clone()).unwrap_or_else(|| "smoke".into());
    let layout_dir = layout_path.parent().unwrap_or(Path::new("."));

    let eff_mode =
        if mode.trim().is_empty() || mode.eq_ignore_ascii_case("auto") || mode.eq_ignore_ascii_case("default") {
            layout.mode.as_deref().or(layout.wm1_level.as_deref()).unwrap_or("auto")
        } else {
            mode
        };

    let rec = if is_electro_dispatch(&layout, eff_mode) {
        run_electro_case(&layout, layout_dir, &tier, eff_mode, scale_h)
    } else {
        run_em_hw_case(&layout, layout_dir, &tier, eff_mode, scale_h)
    };
    Ok(rec)
}

/// Executes `single` in a fresh child process so `nyx_engine::process::peak_memory()` measures exact per-case peak
/// working set & commit (`Plan §4.4`), and optionally enforces a custom `NYX_THREADS` count.
fn dispatch_single_child_process(
    layout_path: &Path,
    tier: &str,
    mode: &str,
    scale_h: f64,
    threads_opt: Option<usize>,
    cluster_override: Option<&str>,
    weight_override: Option<f64>,
) -> Result<RunRecord, String> {
    let exe = std::env::current_exe().map_err(|e| format!("current_exe failed: {e}"))?;
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let tmp_out = std::env::temp_dir().join(format!(
        "squadds_bench_child_{}_{}.json",
        std::process::id(),
        SEQ.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));

    let mut cmd = std::process::Command::new(&exe);
    cmd.arg("single")
        .arg("--layout")
        .arg(layout_path)
        .arg("--tier")
        .arg(tier)
        .arg("--mode")
        .arg(if mode.is_empty() { "auto" } else { mode })
        .arg("--scale-h")
        .arg(format!("{scale_h:.12}"))
        .arg("--out")
        .arg(&tmp_out);
    if let Some(nt) = threads_opt {
        cmd.arg("--threads").arg(nt.to_string());
        cmd.env("NYX_THREADS", nt.to_string());
    }
    if let Some(cid) = cluster_override {
        cmd.arg("--cluster-id").arg(cid);
    }
    if let Some(w) = weight_override {
        cmd.arg("--stratum-weight").arg(format!("{w:.12}"));
    }

    let out = cmd.output().map_err(|e| format!("failed to spawn squadds_bench child: {e}"))?;
    if !out.status.success() {
        let _ = std::fs::remove_file(&tmp_out);
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(format!("child squadds_bench exited with {:?}: {}", out.status.code(), stderr.trim()));
    }
    let rec_txt = std::fs::read_to_string(&tmp_out)
        .map_err(|e| format!("failed to read child output '{}': {e}", tmp_out.display()))?;
    let _ = std::fs::remove_file(&tmp_out);
    serde_json::from_str::<RunRecord>(rec_txt.trim()).map_err(|e| format!("failed to parse child RunRecord: {e}"))
}

fn open_jsonl_writer(path: &Path, append: bool) -> Result<std::fs::File, String> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
    }
    std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(!append)
        .append(append)
        .open(path)
        .map_err(|e| format!("cannot open output '{}': {e}", path.display()))
}

fn write_jsonl_record(file: &mut std::fs::File, rec: &RunRecord) -> Result<(), String> {
    let line = serde_json::to_string(rec).map_err(|e| e.to_string())?;
    writeln!(file, "{line}").map_err(|e| e.to_string())?;
    file.flush().map_err(|e| e.to_string())
}

fn get_opt(args: &[String], flag: &str) -> Option<String> {
    args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1)).cloned()
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

fn print_usage() {
    eprintln!(
        "Usage:\n  \
         squadds_bench calibrate --which <C0.1|C0.2|C0.3|all> --out <dir>\n  \
         squadds_bench single --layout <path.layout.json> [--mode <spec>] [--tier <t>] [--scale-h <f>] [--threads <N>] [--measure-peak-ram] [--append] --out <json_or_jsonl>\n  \
         squadds_bench run --manifest <path> --inputs <dir> --tier <smoke|p1|p2|core|ablation|qcap|airbridge> [--suite <id>] [--mode <spec>] [--scale-h <f>] [--threads <N>] [--limit <K>] [--measure-peak-ram] [--append] --out <runs.jsonl>\n  \
         squadds_bench gci --layout <path.layout.json> [--mode <spec>] [--r <ratio>] [--levels <3|4>] [--threads <N>] [--append] --out <gci.jsonl>\n  \
         squadds_bench determinism --layout <path.layout.json> [--mode <spec>] [--threads-a 1] [--threads-b 8] [--append] --out <det.jsonl>\n  \
         squadds_bench gate --runs <runs.jsonl> --gate <G0|G1|G2|G3|all> --out <gate_report.json>"
    );
}

fn main() {
    nyx_engine::process::opt_out_of_power_throttling();

    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || has_flag(&args, "--help") || has_flag(&args, "-h") {
        print_usage();
        return;
    }

    // Configure NYX_THREADS before the rayon pool in `femcore::par::threads()` is initialized!
    let threads_cli = get_opt(&args, "--threads").and_then(|s| s.parse::<usize>().ok()).filter(|&n| n >= 1);
    if let Some(nt) = threads_cli {
        std::env::set_var("NYX_THREADS", nt.to_string());
    }

    let sub = args[0].as_str();
    let rest = &args[1..];

    let result: Result<(), String> = match sub {
        "calibrate" => {
            let which = get_opt(rest, "--which").unwrap_or_else(|| "all".into());
            let out_dir = PathBuf::from(get_opt(rest, "--out").unwrap_or_else(|| "calibrations".into()));
            run_calibrations(&which, &out_dir).map(|summary| {
                println!("{}", serde_json::to_string_pretty(&summary).unwrap());
            })
        }

        "single" => {
            let Some(lay_str) = get_opt(rest, "--layout") else {
                eprintln!("error: single requires --layout <path.layout.json>");
                std::process::exit(2);
            };
            let layout_path = PathBuf::from(lay_str);
            let mode = get_opt(rest, "--mode").unwrap_or_else(|| "auto".into());
            let tier_opt = get_opt(rest, "--tier");
            let scale_h = get_opt(rest, "--scale-h")
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|&x| x.is_finite() && x > 0.0)
                .unwrap_or(1.0);
            let cluster_opt = get_opt(rest, "--cluster-id");
            let weight_opt = get_opt(rest, "--stratum-weight").and_then(|s| s.parse::<f64>().ok());
            let measure_ram = has_flag(rest, "--measure-peak-ram");
            let append = has_flag(rest, "--append");
            let out_opt = get_opt(rest, "--out").map(PathBuf::from);

            let rec_res = if measure_ram {
                dispatch_single_child_process(
                    &layout_path,
                    tier_opt.as_deref().unwrap_or("smoke"),
                    &mode,
                    scale_h,
                    threads_cli,
                    cluster_opt.as_deref(),
                    weight_opt,
                )
            } else {
                dispatch_single_in_process(
                    &layout_path,
                    tier_opt.as_deref(),
                    &mode,
                    scale_h,
                    cluster_opt.as_deref(),
                    weight_opt,
                )
            };

            rec_res.and_then(|rec| {
                if let Some(ref out_p) = out_opt {
                    let is_jsonl = out_p.extension().and_then(|s| s.to_str()) == Some("jsonl") || append;
                    if is_jsonl {
                        let mut f = open_jsonl_writer(out_p, append)?;
                        write_jsonl_record(&mut f, &rec)?;
                    } else {
                        if let Some(parent) = out_p.parent() {
                            if !parent.as_os_str().is_empty() {
                                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
                            }
                        }
                        std::fs::write(out_p, serde_json::to_string_pretty(&rec).unwrap())
                            .map_err(|e| e.to_string())?;
                    }
                } else {
                    println!("{}", serde_json::to_string_pretty(&rec).unwrap());
                }
                println!(
                    "[single] {} ({}) mode={} status={} cells={} dofs={} total={:.2}s peak_commit={:.2} GiB",
                    rec.case_id,
                    rec.suite_id,
                    rec.mode,
                    rec.status,
                    rec.cells,
                    rec.dofs,
                    rec.timings_s.total,
                    (rec.memory_bytes.peak_commit as f64) / (1024.0 * 1024.0 * 1024.0)
                );
                Ok(())
            })
        }

        "run" => {
            let manifest_path = PathBuf::from(get_opt(rest, "--manifest").unwrap_or_default());
            let inputs_dir = PathBuf::from(get_opt(rest, "--inputs").unwrap_or_else(|| {
                manifest_path.parent().map(|p| p.to_string_lossy().to_string()).unwrap_or_else(|| ".".into())
            }));
            let tier = get_opt(rest, "--tier").unwrap_or_else(|| "smoke".into());
            let suite_filter = get_opt(rest, "--suite");
            let mode_cli = get_opt(rest, "--mode");
            let scale_h = get_opt(rest, "--scale-h")
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|&x| x.is_finite() && x > 0.0)
                .unwrap_or(1.0);
            let limit = get_opt(rest, "--limit").and_then(|s| s.parse::<usize>().ok());
            let measure_ram = has_flag(rest, "--measure-peak-ram");
            let append = has_flag(rest, "--append");
            let Some(out_str) = get_opt(rest, "--out") else {
                eprintln!("error: run requires --out <runs.jsonl>");
                std::process::exit(2);
            };
            let out_path = PathBuf::from(out_str);

            load_manifest_items(&manifest_path, &inputs_dir, &tier, suite_filter.as_deref()).and_then(|mut items| {
                if let Some(lim) = limit {
                    items.truncate(lim);
                }
                if items.is_empty() {
                    return Err(format!(
                        "0 manifest items matched tier='{tier}' suite={suite_filter:?} in '{}'",
                        inputs_dir.display()
                    ));
                }

                let mut out_file = open_jsonl_writer(&out_path, append)?;

                // Determine mode expansion list (e.g., all 7 ablation modes when `--tier ablation` and `--mode` is omitted or `all`)
                let ablation_all = tier.eq_ignore_ascii_case("ablation")
                    && mode_cli
                        .as_deref()
                        .is_none_or(|m| m.eq_ignore_ascii_case("all") || m.eq_ignore_ascii_case("ablation7"));

                let mut n_done = 0usize;
                let mut n_ok = 0usize;
                for item in &items {
                    let modes_for_item: Vec<String> = if ablation_all {
                        vec![
                            "g2_dual_box_raw".into(),
                            "g2_dual_box".into(),
                            "far_field".into(),
                            "amr_ps".into(),
                            "qp_graded".into(),
                            "thick_conformal".into(),
                            "planar_bem".into(),
                        ]
                    } else if let Some(ref m_str) = mode_cli {
                        m_str.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect()
                    } else if let Some(ref mo) = item.mode_override {
                        vec![mo.clone()]
                    } else {
                        vec!["auto".into()]
                    };

                    for m in &modes_for_item {
                        let rec = if measure_ram {
                            dispatch_single_child_process(
                                &item.layout_path,
                                &tier,
                                m,
                                scale_h,
                                threads_cli,
                                item.cluster_id.as_deref(),
                                item.stratum_weight,
                            )
                            .unwrap_or_else(|err| RunRecord {
                                case_id: item.case_id.clone(),
                                suite_id: item.suite_id.clone(),
                                tier: tier.clone(),
                                mode: m.clone(),
                                cluster_id: item.cluster_id.clone().unwrap_or_else(|| item.case_id.clone()),
                                stratum_weight: item.stratum_weight.unwrap_or(1.0),
                                status: "solver_nonconv".into(),
                                failure_detail: Some(err),
                                ..Default::default()
                            })
                        } else {
                            dispatch_single_in_process(
                                &item.layout_path,
                                Some(&tier),
                                m,
                                scale_h,
                                item.cluster_id.as_deref(),
                                item.stratum_weight,
                            )?
                        };

                        if rec.status == "ok" {
                            n_ok += 1;
                        }
                        n_done += 1;
                        write_jsonl_record(&mut out_file, &rec)?;
                        println!(
                            "[{n_done}] {} ({}) mode={} -> {} (cells={}, dofs={}, {:.2}s)",
                            rec.case_id, rec.suite_id, rec.mode, rec.status, rec.cells, rec.dofs, rec.timings_s.total
                        );
                    }
                }
                println!(
                    "[run complete] {n_ok}/{n_done} ok ({:.1}%) -> {}",
                    100.0 * (n_ok as f64) / (n_done.max(1) as f64),
                    out_path.display()
                );
                Ok(())
            })
        }

        "gci" => (|| -> Result<(), String> {
            let lay_str = get_opt(rest, "--layout").ok_or("gci requires --layout <path.layout.json>")?;
            let layout_path = PathBuf::from(lay_str);
            let mode = get_opt(rest, "--mode").unwrap_or_else(|| "auto".into());
            let txt = std::fs::read_to_string(&layout_path).unwrap_or_default();
            let lay_probe: Option<LayoutSpec> = parse_json_lenient_nan(&txt).ok();
            let is_el = lay_probe.as_ref().is_some_and(|l| is_electro_dispatch(l, &mode));

            let levels: usize =
                get_opt(rest, "--levels").and_then(|s| s.parse().ok()).unwrap_or(if is_el { 4 } else { 3 });
            let r_ratio: f64 =
                get_opt(rest, "--r").and_then(|s| s.parse().ok()).unwrap_or(if is_el { 1.35 } else { 1.25 });
            let append = has_flag(rest, "--append");
            let out_str = get_opt(rest, "--out").ok_or("gci requires --out <gci.jsonl>")?;
            let out_path = PathBuf::from(out_str);

            let exponents: Vec<i32> = if levels >= 4 { vec![-1, 0, 1, 2] } else { vec![-1, 0, 1] };
            let mut recs = Vec::with_capacity(exponents.len());
            for &ell in &exponents {
                // Plan §5.2: h^(ell) = h^(0) / r^ell  => scale_h = r^(-ell)
                let scale_h = r_ratio.powi(-ell);
                let mut rec = dispatch_single_in_process(&layout_path, Some("gci"), &mode, scale_h, None, None)?;
                rec.mode = format!("{}_gci_l{ell:+}", rec.mode);
                recs.push(rec);
            }

            // Compute Roache GCI on the primary observable across the 3 finest levels
            if recs.len() >= 3 && recs.iter().all(|r| r.status == "ok") {
                let n = recs.len();
                let target_key = ["cross_to_cross", "top_to_top", "north_to_north", "f_r_hz", "f01_hz", "c_0_0_f"]
                    .into_iter()
                    .find(|&k| recs[n - 1].predictions_si.contains_key(k))
                    .map(String::from)
                    .or_else(|| recs[n - 1].predictions_si.keys().next().cloned());

                if let Some(tk) = target_key {
                    if let (Some(&yc), Some(&ym), Some(&yf)) = (
                        recs[n - 3].predictions_si.get(&tk),
                        recs[n - 2].predictions_si.get(&tk),
                        recs[n - 1].predictions_si.get(&tk),
                    ) {
                        let p_nom = if is_el { 2.0 } else { 1.5 };
                        let summary = compute_roache_gci(yc, ym, yf, r_ratio, p_nom);
                        println!(
                            "[gci] {} ({}) observable={tk}: p_obs={:.2} (clamped {:.2}), GCI_fine={:.4}%",
                            recs[n - 1].case_id,
                            recs[n - 1].suite_id,
                            summary.observed_order_raw,
                            summary.observed_order_clamped,
                            summary.gci_fine_pct
                        );
                        recs[n - 1].extra = json!({
                            "gci_observable": tk,
                            "gci_summary": summary,
                            "prev_extra": recs[n - 1].extra,
                        });
                    }
                }
            }

            let mut f = open_jsonl_writer(&out_path, append)?;
            for r in &recs {
                write_jsonl_record(&mut f, r)?;
            }
            Ok(())
        })(),

        "determinism" => {
            let Some(lay_str) = get_opt(rest, "--layout") else {
                eprintln!("error: determinism requires --layout <path.layout.json>");
                std::process::exit(2);
            };
            let layout_path = PathBuf::from(lay_str);
            let mode = get_opt(rest, "--mode").unwrap_or_else(|| "auto".into());
            let ta: usize = get_opt(rest, "--threads-a").and_then(|s| s.parse().ok()).unwrap_or(1);
            let tb: usize = get_opt(rest, "--threads-b").and_then(|s| s.parse().ok()).unwrap_or(8);
            let append = has_flag(rest, "--append");
            let Some(out_str) = get_opt(rest, "--out") else {
                eprintln!("error: determinism requires --out <det.jsonl>");
                std::process::exit(2);
            };
            let out_path = PathBuf::from(out_str);

            dispatch_single_child_process(&layout_path, "determinism", &mode, 1.0, Some(ta), None, None).and_then(
                |mut ra| {
                    dispatch_single_child_process(&layout_path, "determinism", &mode, 1.0, Some(tb), None, None)
                        .and_then(|mut rb| {
                            ra.mode = format!("{}_det_t{ta}", ra.mode);
                            rb.mode = format!("{}_det_t{tb}", rb.mode);
                            let mut max_rel = 0.0_f64;
                            for (k, &va) in &ra.predictions_si {
                                if let Some(&vb) = rb.predictions_si.get(k) {
                                    if va.abs() > 1.0e-30 {
                                        max_rel = max_rel.max((vb - va).abs() / va.abs());
                                    }
                                }
                            }
                            ra.predictions_si.insert("det_max_rel_diff".into(), max_rel);
                            rb.predictions_si.insert("det_max_rel_diff".into(), max_rel);
                            println!(
                        "[determinism] {} ({}): NYX_THREADS={ta} vs {tb} max_rel_diff={max_rel:.3e} (pass < 1e-13: {})",
                        ra.case_id,
                        ra.suite_id,
                        max_rel < 1.0e-13
                    );
                            let mut f = open_jsonl_writer(&out_path, append)?;
                            write_jsonl_record(&mut f, &ra)?;
                            write_jsonl_record(&mut f, &rb)?;
                            Ok(())
                        })
                },
            )
        }

        "gate" => {
            let Some(runs_str) = get_opt(rest, "--runs") else {
                eprintln!("error: gate requires --runs <runs.jsonl>");
                std::process::exit(2);
            };
            let gate = get_opt(rest, "--gate").unwrap_or_else(|| "all".into());
            let out_path = PathBuf::from(get_opt(rest, "--out").unwrap_or_else(|| "gate_report.json".into()));
            evaluate_gates(Path::new(&runs_str), &gate, &out_path).map(|rep| {
                println!("{}", serde_json::to_string_pretty(&rep["gates"]).unwrap());
                println!(
                    "[gate] overall completion: {}/{} ({:.1}%) -> {}",
                    rep["overall_completion"]["n_ok"],
                    rep["overall_completion"]["n_total"],
                    rep["overall_completion"]["completion_pct"].as_f64().unwrap_or(0.0),
                    out_path.display()
                );
            })
        }

        "bug-checks" | "bug_checks" => {
            let inputs_dir = get_opt(rest, "--inputs").map(PathBuf::from);
            let append = has_flag(rest, "--append");
            let out_path = PathBuf::from(get_opt(rest, "--out").unwrap_or_else(|| "bug_checks.jsonl".into()));
            let records = run_all_bug_checks(inputs_dir.as_deref());
            (|| -> Result<(), String> {
                let mut f = open_jsonl_writer(&out_path, append)?;
                for rec in &records {
                    write_jsonl_record(&mut f, rec)?;
                }
                let n_failed = records.iter().filter(|r| r.status != "ok").count();
                println!(
                    "[bug-checks] wrote {} records ({} failed / reproduced) -> {}",
                    records.len(),
                    n_failed,
                    out_path.display()
                );
                Ok(())
            })()
        }

        other => Err(format!("unknown subcommand '{other}'")),
    };

    if let Err(err) = result {
        eprintln!("squadds_bench error: {err}");
        std::process::exit(1);
    }
}
