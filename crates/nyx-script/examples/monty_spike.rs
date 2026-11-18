//! S4 spike runner for the Monty binding: `monty_spike sample [file]` runs the strip sweep, `monty_spike bench` measures the cost of a call,
//! `monty_spike limits` runs the runaway scripts.

use std::time::Instant;

use nyx_engine::Engine;
use nyx_script::{monty_host, Limits, Outcome};
use serde_json::json;

fn show(label: &str, r: Result<Outcome, String>, t: Instant) {
    match r {
        Ok(o) => {
            for l in &o.printed {
                println!("  | {l}");
            }
            println!("{label}: ok in {:.2} s, {} engine calls, value {}", t.elapsed().as_secs_f64(), o.calls, o.value);
        }
        Err(e) => println!("{label}: FAILED after {:.2} s\n{e}", t.elapsed().as_secs_f64()),
    }
}

fn fresh() -> Engine {
    let mut e = Engine::new();
    e.call("project_new", &json!({ "dim": 2, "title": "spike" })).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": 5e-7 })).unwrap();
    e
}

fn bench() {
    let n: u64 = 20_000;
    let limits = Limits { seconds: 600.0, calls: 1_000_000, ..Limits::default() };
    println!("{:<34} {:>12}", "per call, microseconds", "monty");
    for (name, call) in [
        ("set_param (edit)", "nyx.set_param(name=\"w\", value=5e-7)"),
        ("project_info (read, reply of a map)", "nyx.project_info()"),
        ("list_materials (read, reply of a list)", "nyx.list_materials()"),
    ] {
        let mut e = fresh();
        let src = format!("for i in range({n}):\n    {call}\n");
        let t = Instant::now();
        monty_host::run(&mut e, "bench.py", &src, &limits).unwrap();
        println!("{name:<34} {:>12.2}", t.elapsed().as_secs_f64() / n as f64 * 1e6);
    }
    let mut e = fresh();
    let m: u64 = 2_000_000;
    let t = Instant::now();
    monty_host::run(&mut e, "loop.py", &format!("s = 0.0\nfor i in range({m}):\n    s += i * 0.5\ns\n"), &limits).unwrap();
    println!("{:<34} {:>12.3}   (loop of arithmetic, no engine call)", "per iteration, microseconds", t.elapsed().as_secs_f64() / m as f64 * 1e6);
    let t = Instant::now();
    for _ in 0..50 {
        monty_host::run(&mut e, "empty.py", "1", &limits).unwrap();
    }
    println!("{:<34} {:>12.1}   (an empty script, bindings of {} commands)", "start-up, microseconds", t.elapsed().as_secs_f64() / 50.0 * 1e6, e.catalog().len());
}

/// What a person writing Python would try. Prints the value, or the last line of the error.
fn probe() {
    let lim = Limits { seconds: 5.0, ..Limits::default() };
    let mut e = fresh();
    for (label, src) in [
        ("import math", "import math\nmath.sqrt(2) + math.pi"),
        ("from math import pi", "from math import pi\npi"),
        ("import json", "import json\njson.dumps({'a': 1})"),
        ("import time", "import time\ntime.time()"),
        ("import random", "import random\nrandom.random()"),
        ("import numpy", "import numpy as np\nnp.linspace(0, 1, 3)"),
        ("class", "class A:\n    x = 1\nA().x"),
        ("dataclass-free tuple", "a, b = (1, 2)\nb, a"),
        ("% format", "'%5.2f' % 3.14159"),
        ("str.format", "'{:.3f}'.format(1.2345)"),
        ("f-string spec", "x = 12345.6789\nf'{x:>12.3f} {x:e} {x:,.1f}'"),
        ("sorted key lambda", "sorted([3, 1, 2], key=lambda v: -v)"),
        ("zip/enumerate", "[(i, a * b) for i, (a, b) in enumerate(zip([1, 2], [3, 4]))]"),
        ("sum/min/max/round/divmod", "[sum([1, 2]), min(3, 4), max([5, 6]), round(2.567, 2), divmod(7, 2), 2 ** 10, 7 // 2]"),
        ("generator expr", "sum(x * x for x in range(5))"),
        ("slicing", "[1, 2, 3, 4][::-1]"),
        ("dict comprehension", "{k: k * 2 for k in range(3)}"),
        ("global in function", "c = 0\ndef f():\n    global c\n    c += 1\nf()\nc"),
        ("closure", "def mk(n):\n    def g(x):\n        return x + n\n    return g\nmk(2)(3)"),
        ("try/finally", "try:\n    1 / 0\nexcept ZeroDivisionError as e:\n    r = 'caught'\nfinally:\n    r = r + '!'\nr"),
        ("assert", "assert 1 == 2, 'one is not two'"),
        ("with statement", "with open('x') as f:\n    pass"),
        ("print sep", "print(1, 2, sep=',')\n3"),
        ("last stmt is assignment", "result = 5"),
        ("last expr", "result = 5\nresult"),
        ("float repr", "[0.1 + 0.2, 1 / 3, 1e-9, 5e-7, 10 ** 20, 2 ** 70]"),
        ("int dict key", "x = {1: 'a'}\nnyx.call('project_info', **{'a': 1})"),
        ("unknown command", "nyx.nope()"),
        ("misspelled argument", "nyx.set_param(nme='w', value=1)"),
        ("positional arguments", "nyx.set_param('w', 1)"),
        ("engine error is catchable", "try:\n    nyx.set_param(name='w')\nexcept RuntimeError as e:\n    msg = str(e)\nmsg"),
        ("nan to engine", "nyx.set_param(name='w', value=float('nan'))"),
        ("undefined name", "print(undefined_thing)"),
        ("syntax error", "x = (1,\n"),
        ("tab/indent error", "if True:\nprint(1)\n"),
        ("error in a later line", "a = 1\nb = 2\nnyx.set_param(name='w')\nc = 3\n"),
        ("input()", "input('x')"),
        ("exec/eval", "eval('1 + 1')"),
        ("match statement", "match 1:\n    case 1:\n        r = 'one'\nr"),
        ("type hints", "def f(x: float) -> float:\n    return x * 2\nf(2.0)"),
        ("walrus", "(n := 5) + n"),
        ("list methods", "x = [3, 1]\nx.append(2)\nx.sort()\nx.extend([9])\nx"),
        ("str methods", "' a,b '.strip().split(',')"),
        ("shadow nyx", "nyx = 5\nnyx"),
    ] {
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| monty_host::run(&mut e, "x.py", src, &lim))).unwrap_or_else(|_| Err("PANIC inside the interpreter".into()));
        let msg: String = match r {
            Ok(o) => format!("ok: {} {}", o.value, o.printed.join("|")),
            Err(m) => format!("ERR {}", m.lines().rfind(|l| !l.trim().is_empty() && !l.contains("Traceback") && !l.starts_with(' ')).unwrap_or("")),
        };
        println!("{label:<28} {}", msg.chars().take(150).collect::<String>());
    }
    // one full error message, as a person sees it
    let src = "a = 1\nb = 2\nnyx.set_param(name='w')\nc = 3\n";
    if let Err(m) = monty_host::run(&mut e, "model.py", src, &lim) {
        println!("--- full error of a bad call, line 3 of model.py:\n{m}");
    }
}

fn limits_demo(alloc_case: Option<usize>) {
    let tight = Limits { seconds: 1.0, calls: 200, depth: 50, heap_bytes: 64 << 20, ..Limits::default() };
    let mut e = fresh();
    // the cases marked true allocate without bound: they abort the process unless a memory limit is armed, so each runs alone (`limits N`)
    let cases = [
        ("while True", "while True:\n    pass\n", false),
        ("huge loop", "s = 0\nfor i in range(10**12):\n    s += i\ns\n", false),
        ("unbounded list", "x = [0]\nwhile True:\n    x = x + x\n", true),
        ("one huge list", "len([0] * 10**12)", true),
        ("unbounded string", "s = 'x'\nwhile True:\n    s = s + s\n", true),
        ("deep recursion", "def f(n):\n    return f(n + 1)\nf(0)\n", false),
        ("too many engine calls", "for i in range(1000):\n    nyx.project_info()\n", false),
        ("import os", "import os\nos.listdir('.')\n", false),
        ("open a file", "open('C:/Windows/win.ini').read()", false),
        ("print flood", "while True:\n    print('x' * 1000)\n", false),
    ];
    for (i, (label, src, alloc)) in cases.into_iter().enumerate() {
        if alloc != alloc_case.is_some() || alloc_case.is_some_and(|n| n != i) {
            continue;
        }
        let t = Instant::now();
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| monty_host::run(&mut e, "x.py", src, &tight))).unwrap_or_else(|_| Err("PANIC inside the interpreter".into()));
        let msg = match r {
            Ok(o) => format!("NOT STOPPED, value {}", o.value),
            Err(m) => m.lines().rfind(|l| !l.trim().is_empty() && !l.contains("Traceback") && !l.starts_with(' ')).unwrap_or("").chars().take(110).collect(),
        };
        println!("{label:<24} stopped in {:>6.2} s: {msg}", t.elapsed().as_secs_f64());
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("bench") => bench(),
        Some("limits") => limits_demo(args.get(1).and_then(|a| a.parse().ok())),
        Some("probe") => probe(),
        Some("sample") => {
            let file = args.get(1).cloned().unwrap_or_else(|| "docs/scripting-samples/strip_sweep.monty".to_string());
            let src = std::fs::read_to_string(&file).unwrap_or_else(|e| panic!("{file}: {e}"));
            let mut e = Engine::new();
            let t = Instant::now();
            show("monty", monty_host::run(&mut e, &file, &src, &Limits::default()), t);
        }
        _ => eprintln!("usage: monty_spike bench | limits | sample [file]"),
    }
}
