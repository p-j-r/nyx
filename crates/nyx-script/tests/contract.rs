//! The script binding's contract (`script_api.v1`, ADR 93, ADR 19): what a script can rely on, written from the constants the code uses so the
//! exemplar cannot drift from them. The commands a script calls are `engine_commands`, not this file. `NYX_WRITE_CONTRACT=1` rewrites the exemplar;
//! do that only together with a line in `docs/architecture/contracts/VERSIONS.md`.

use nyx_engine::path_guard::PATH_ARGS;
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
        "output": "print() lines, up to the output limit",
        "engine_error": "RuntimeError('<command>: <message>'), catchable",
        "static_check": "names and required-ness of commands and arguments, against the stub of nyx api stubs; value types are not checked",
        "limits": {
            "seconds": l.seconds,
            "memory_bytes": l.memory_bytes,
            "depth": l.depth,
            "calls": l.calls,
            "output_bytes": l.output_bytes,
            "seconds_counts": "script time only; the time the engine takes to answer a command is not counted",
        },
        "path_guard": { "arguments": PATH_ARGS, "root": "the workspace given to the run, else the folder of the script; relative paths are taken from it" },
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
