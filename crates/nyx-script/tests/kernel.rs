//! The kernel (productivity plan S7, ADR [94]): one session that outlives a script, a command that waits for a job that cannot block for ever,
//! `params`, the futures of `nyx.aio`, the refusal of the commands that discard the project, and a worker that dies replaced from the last dump.

use std::path::PathBuf;
use std::time::Duration;

use nyx_engine::Engine;
use nyx_script::kernel::{Kernel, Step};
use nyx_script::{Limits, Options};
use serde_json::{json, Value};

fn root(name: &str) -> PathBuf {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!("nyx-kernel-test-{name}-{}-{}", std::process::id(), N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn kernel(e: &Engine, limits: Limits) -> Kernel {
    let opts = Options { root: root("k"), worker: Some(PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))) };
    Kernel::new(e, limits, opts).unwrap()
}

/// Runs to the end (no time limit) and returns the value.
fn go(k: &mut Kernel, e: &mut Engine, src: &str) -> Result<Value, String> {
    match k.run(e, "t.nys", src, None, None, true)? {
        Step::Done(o) => Ok(o.value),
        Step::Running(p) => panic!("still running, waiting for {}", p.waiting),
    }
}

/// A 2D silicon wire that solves in a few seconds, not yet solved.
const MODEL: &str = r#"
nyx.project_new(dim=2, title="kernel")
nyx.set_param(name="w", value=5e-7)
nyx.add_rect(name="box", min=[-2.0e-6, -2.0e-6], max=[2.0e-6, 0.0])
nyx.add_rect(name="clad", min=[-2.0e-6, 0.0], max=[2.0e-6, 2.0e-6])
nyx.add_rect(name="core", min=["-w/2", 0.0], max=["w/2", 2.2e-7])
nyx.add_overlay(regions=[{"name": "BOX", "input": "box"}, {"name": "Clad", "input": "clad"}, {"name": "Si", "input": "core"}])
nyx.add_physics(physics="ems", study="modes", study_params={"freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 2})
for reg, mat in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")]:
    nyx.set_material(physics="ems", region=reg, material=mat)
nyx.set_mesh(h=5e-7, adapt_scale=0.4)
"#;

#[test]
fn variables_and_functions_live_from_one_script_to_the_next_and_check_leaves_nothing() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    go(&mut k, &mut e, "x = 41\ndef f(a):\n    return a + x\n").unwrap();
    assert_eq!(go(&mut k, &mut e, "f(1)").unwrap(), json!(42));
    // the static check knows the variables of the session, and still refuses what is not there
    k.check(&mut e, "c.nys", "f(x) + 1").unwrap();
    let err = k.check(&mut e, "c.nys", "zz + 1").unwrap_err();
    assert!(err.contains("zz"), "{err}");
    let err = k.check(&mut e, "c.nys", "nyx.set_param(nam='w', value=1)").unwrap_err();
    assert!(err.contains("nam") || err.contains("name"), "{err}");
    // a checked script is not run: its assignment is not in the session
    k.check(&mut e, "c.nys", "q = 5").unwrap();
    assert!(go(&mut k, &mut e, "q").unwrap_err().contains("q"), "q must not exist");
    // and a failed script does not end the session
    assert_eq!(go(&mut k, &mut e, "x").unwrap(), json!(41));
}

#[test]
fn params_reach_the_script_and_a_default_stands_in_for_a_missing_one() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    let src = "params.get('w', 7.0) * 2";
    let Step::Done(o) = k.run(&mut e, "t.nys", src, Some(&json!({ "w": 3.0 })), None, true).unwrap() else { panic!() };
    assert_eq!(o.value, json!(6.0));
    let Step::Done(o) = k.run(&mut e, "t.nys", src, None, None, true).unwrap() else { panic!() };
    assert_eq!(o.value, json!(14.0));
    // a script that reads params type-checks
    k.check(&mut e, "c.nys", "params['w'] + 1").unwrap();
}

#[test]
fn a_command_that_waits_for_a_job_yields_then_continues_and_gives_the_value_of_the_blocking_call() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    go(&mut k, &mut e, MODEL).unwrap();
    let src = "print('before')\nr = nyx.solve_run(physics='ems')\nprint('after')\nr";
    // 1 ms: the solve takes longer, so the script stays paused at the call
    let first = k.run(&mut e, "t.nys", src, None, Some(Duration::from_millis(1)), true).unwrap();
    let Step::Running(p) = first else { panic!("a solve cannot finish in 1 ms") };
    assert_eq!(p.printed, vec!["before".to_string()]);
    assert!(p.waiting.contains("solve_run"), "{}", p.waiting);
    assert!(k.waiting().is_some());
    // another script cannot start while this one waits
    let err = k.run(&mut e, "t2.nys", "1", None, None, true).unwrap_err();
    assert!(err.contains("still running") && err.contains("script_wait"), "{err}");
    let done = k.wait(&mut e, Some(Duration::from_secs(120))).unwrap();
    let Step::Done(o) = done else { panic!("the solve takes less than two minutes") };
    assert_eq!(o.printed, vec!["after".to_string()]);
    assert!(o.value.get("job").is_some() || o.value.get("results").is_some() || o.value.is_object(), "{}", o.value);
    // the commands the host sent to poll are not the script's: one solve_run is one call
    assert_eq!(o.calls, 1);
    assert!(k.waiting().is_none());
    // the engine has the result, as after a blocking solve
    let info = e.call("result_info", &Value::Null).unwrap();
    assert!(info.to_string().contains("ems"), "{info}");
}

#[test]
fn cancel_abandons_a_waiting_script_and_the_variables_of_the_finished_ones_stay() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    go(&mut k, &mut e, MODEL).unwrap();
    go(&mut k, &mut e, "kept = 123").unwrap();
    let first = k.run(&mut e, "t.nys", "lost = 1\nnyx.solve_run(physics='ems')", None, Some(Duration::from_millis(1)), true).unwrap();
    assert!(matches!(first, Step::Running(_)));
    k.cancel(&mut e).unwrap();
    assert!(k.waiting().is_none());
    assert!(k.cancel(&mut e).is_err(), "nothing is left to cancel");
    assert_eq!(go(&mut k, &mut e, "kept").unwrap(), json!(123));
    assert!(go(&mut k, &mut e, "lost").is_err(), "the abandoned script's variables are gone");
}

#[test]
fn a_job_that_cannot_start_is_an_error_not_an_endless_wait() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    go(&mut k, &mut e, MODEL).unwrap();
    go(&mut k, &mut e, "nyx.queue_pause()").unwrap();
    let t = std::time::Instant::now();
    let err = k.run(&mut e, "t.nys", "nyx.solve_run(physics='ems')", None, None, true).unwrap_err();
    assert!(err.contains("cannot start") && err.contains("paused"), "{err}");
    assert!(t.elapsed().as_secs() < 20, "took {:?}", t.elapsed());
    // a wait with its own timeout on such a job returns its entry, as the engine's does
    go(&mut k, &mut e, "nyx.queue_resume()").unwrap();
}

#[test]
fn commands_that_discard_the_project_need_the_flag() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    for src in ["nyx.project_new(dim=2, title='x')", "nyx.call('project_new', dim=2, title='x')"] {
        let err = k.run(&mut e, "t.nys", src, None, None, false).unwrap_err();
        assert!(err.contains("allow_destructive"), "{err}");
    }
    // a script can still catch the refusal and the session goes on
    let Step::Done(o) = k.run(&mut e, "t.nys", "try:\n    nyx.project_new(dim=2)\n    r = 'ran'\nexcept RuntimeError:\n    r = 'refused'\nr", None, None, false).unwrap() else { panic!() };
    assert_eq!(o.value, json!("refused"));
    assert!(matches!(k.run(&mut e, "t.nys", "nyx.project_new(dim=2, title='x')", None, None, true).unwrap(), Step::Done(_)));
    // a command that does not discard anything needs nothing
    assert!(matches!(k.run(&mut e, "t.nys", "nyx.project_info()", None, None, false).unwrap(), Step::Done(_)));
}

#[test]
fn futures_gather_the_waits_and_the_value_is_the_value_of_the_blocking_call() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    go(&mut k, &mut e, MODEL).unwrap();
    let r = go(&mut k, &mut e, "import asyncio\nrs = await asyncio.gather(nyx.aio.solve_run(physics='ems'), nyx.aio.solve_run(physics='ems'))\nlen(rs)").unwrap();
    assert_eq!(r, json!(2));
    // a single awaited future gives what the blocking call gives
    let a = go(&mut k, &mut e, "r = await nyx.aio.solve_run(physics='ems')\nsorted(r.keys())").unwrap();
    let b = go(&mut k, &mut e, "sorted(nyx.solve_run(physics='ems').keys())").unwrap();
    assert_eq!(a, b);
    // a failure raises at the await and the script can catch it
    let r = go(&mut k, &mut e, "try:\n    await nyx.aio.job_wait(id=987654)\n    r = 'no error'\nexcept RuntimeError as ex:\n    r = 'caught'\nr").unwrap();
    assert_eq!(r, json!("caught"));
    // an await in a script that is not a coroutine of its own also type-checks
    k.check(&mut e, "c.nys", "import asyncio\nawait asyncio.gather(nyx.aio.solve_run(physics='ems'))").unwrap();
}

#[test]
fn a_future_left_pending_has_its_job_cancelled() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    go(&mut k, &mut e, MODEL).unwrap();
    // started and never awaited: the script ends, the job must not be left running for nobody
    go(&mut k, &mut e, "f = nyx.aio.solve_run(physics='ems')\n1").unwrap();
    let t = std::time::Instant::now();
    loop {
        let jobs = e.call("job_list", &Value::Null).unwrap().to_string();
        if !jobs.contains("\"running\"") && !jobs.contains("\"queued\"") {
            break;
        }
        assert!(t.elapsed().as_secs() < 60, "the job is still going: {jobs}");
        std::thread::sleep(Duration::from_millis(100));
    }
}

#[test]
fn a_worker_that_dies_is_replaced_and_the_variables_of_the_finished_scripts_are_back() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits { memory_bytes: 64 << 20, seconds: 20.0, ..Limits::default() });
    go(&mut k, &mut e, "kept = 5").unwrap();
    let err = go(&mut k, &mut e, "big = [0]\nwhile True:\n    big = big + big\n").unwrap_err();
    assert!(err.to_lowercase().contains("memory"), "{err}");
    assert_eq!(go(&mut k, &mut e, "kept + 1").unwrap(), json!(6));
    // the engine is the same one throughout
    e.call("project_new", &json!({ "dim": 2, "title": "alive" })).unwrap();
}

#[test]
fn a_script_error_names_the_line_and_the_session_goes_on() {
    let mut e = Engine::new();
    let mut k = kernel(&e, Limits::default());
    let err = go(&mut k, &mut e, "a = 1\nnyx.set_param(name='w', value=None)\n").unwrap_err();
    assert!(err.contains("line 2") && err.contains("set_param"), "{err}");
    assert_eq!(go(&mut k, &mut e, "a = 2\na + 1").unwrap(), json!(3));
}
