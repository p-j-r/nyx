// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The script binding's contract (`script_api.v1`, ADR 93, ADR 19): what a script can rely on, written from the constants the code uses so the
//! exemplar cannot drift from them. The commands a script calls are `engine_commands`, not this file. `NYX_WRITE_CONTRACT=1` rewrites the exemplar;
//! do that only together with a line in `docs/architecture/contracts/VERSIONS.md`.

use nyx_engine::path_guard::PATH_ARGS;
use nyx_script::kernel::{MAX_TIMEOUT, POLLABLE};
use nyx_script::trace::{CONDS, FORMULA_ARGS, FUNCS_1, FUNCS_2};
use nyx_script::Limits;
use serde_json::{json, Value};
use std::path::PathBuf;

fn current() -> Value {
    let l = Limits::default();
    json!({
        "format": "script_api.v1",
        "extension": ".nys",
        "language": "Python syntax, run by Monty 1.6.0 in a worker process (no CPython)",
        "object": "nyx",
        "calls": ["nyx.<command>(name=value, ...)", "nyx.call(name, **args)", "nyx.call(name, {...})", "nyx.<command>({...})"],
        "value": "the last expression of the script, as JSON",
        "output": "print() lines, truncated at a UTF-8 char boundary so total bytes stay within output_bytes",
        "engine_error": "RuntimeError('<command>: <message>'), catchable",
        "static_check": "names and required-ness of commands and arguments, against the stub of nyx api stubs; value types are not checked",
        "limits": {
            "seconds": l.seconds,
            "max_timeout_seconds": MAX_TIMEOUT.as_secs(),
            "memory_bytes": l.memory_bytes,
            "depth": l.depth,
            "calls": l.calls,
            "output_bytes": l.output_bytes,
            "seconds_counts": "script time only; the time the engine takes to answer a command is not counted",
            "timeout_validation": "non-finite, negative, or Duration-overflowing seconds/timeout_s values return a structured validation error; finite values are clamped to max_timeout_seconds",
        },
        "path_guard": { "arguments": PATH_ARGS, "root": "the workspace given to the run, else the folder of the script; relative paths are taken from it" },
        "network": {
            "refused_commands": ["fetch_material", "search_remote_materials(live=true)"],
            "rule": "outbound network access is denied during sandboxed script execution (RuntimeError)",
        },
        "determinism": { "clock": "refused (RuntimeError)", "random": "seed 0", "sleep": "returns at once" },
        "files": "open, pathlib and os raise PermissionError",
        "run_record": { "artifact_kind": "script", "file": "script-<blake3, 32 hex digits>.nys" },
        "worker": { "program": "nyx-script-worker", "environment": "NYX_SCRIPT_WORKER", "argument": "subprocess" },
        "traced_function": {
            "form": "nyx.<command>(expr=f, ...) with f a def or a lambda of one argument, the point p; also trace(f) for its text",
            "formula_arguments": FORMULA_ARGS,
            "point": "p[\"E\"] (a vector, with .x .y .z), p[\"eps_r\"], p[\"Ex\"], p[\"nx\"]: a name of a formula; p.x p.y p.z the coordinates",
            "functions_of_one_value": FUNCS_1,
            "functions_of_two_values": FUNCS_2.iter().map(|(n, _)| *n).collect::<Vec<_>>(),
            "conditions": CONDS.iter().map(|(n, _)| *n).chain(["lnot"]).collect::<Vec<_>>(),
            "other": ["where(cond, a, b)", "material(name)", "region(name)", "abs(x)", "+ - * / ** and unary minus"],
            "called": "once per command, with a symbolic point; the text it builds is the engine formula of the same words (parsed by femcore::formula)",
            "branch_rule": "if, while, a if c else b, and, or, not, assert, a comparison, math.f and min/max on a value that depends on the point are refused before the script runs, with the line; write where(cond, a, b)",
            "guard_limit": "syntactic: a value the guard cannot follow (a function chosen at run time) is not seen, and then `if sym:` takes the true branch",
            "numbers": "written with repr: finite int or float; an integer stays an integer",
        },
        "kernel": {
            "session": "one Monty session outlives a script: variables, functions and imports stay from one script to the next; the static check knows them",
            "params": "the dict given to the run, bound in the script as `params` (an empty dict when none)",
            "waiting_commands": POLLABLE,
            "waiting_rule": "started without waiting, then polled by the host; when the caller's time is spent the script stays paused at the call and the run answers `running`; the value is the value of the blocking call; the polling is not a script command",
            "blocked_job": "a job that cannot start (held, queue paused) is an error, not an endless wait",
            "futures": { "namespace": "nyx.aio", "commands": POLLABLE, "use": "await asyncio.gather(nyx.aio.solve_run(...), ...)", "pending_at_end": "the job of a future left pending when a script ends is cancelled", "missing_in_monty_1_6_0": ["asyncio.gather(return_exceptions=True)", "asyncio.create_task"] },
            "destructive_commands": "project_new, project_load, project_open, store_open and remove_*: RuntimeError unless the run allows them; `nyx run` allows them",
            "worker_loss": "a worker that dies is replaced and restored from the dump taken after the last script that finished",
            "check": "type-checks in a copy of the session; what the part of the script before its first command raises at run time (it has no params) is not an error of the check",
            "guard_limit_kernel": "the traced-function guard reads one script: a function defined by an earlier script and given as a formula in a later one is not seen",
        },
    })
}

#[test]
fn the_script_api_matches_its_exemplar() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts/script_api.v1.json");
    let now = current();
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::write(&path, serde_json::to_string_pretty(&now).unwrap() + "\n").unwrap();
        return;
    }
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e} (NYX_WRITE_CONTRACT=1 writes it)", path.display()));
    let exemplar: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(exemplar, now, "the script API changed: bump the version in VERSIONS.md and rewrite the exemplar with NYX_WRITE_CONTRACT=1");
}
