// strip_sweep.rhai -- SAMPLE ONLY: nyx has no Rhai binding yet. Written to compare the language with strip_sweep.star.
//
// Assumed binding (what the host would register; generated from Engine::catalog()):
//   nyx::call(name, map)        -> map      any engine command, the same JSON as an MCP tool call
//   nyx::add_rect(map), nyx::set_material(map), ...   one function per catalog command, taking one map
// Errors from the engine are thrown as Rhai errors, so try/catch works.
// The host sets limits (max operations, call depth, string/array sizes) and gives no file or network access.
//   fmt(x, digits)              fixed-point text of a number (Rhai has no printf; the same host function as in the Starlark sample)
//
// Workflow: a silicon strip on oxide at 1550 nm.
//   1. build the model with a parameter `w` (strip width)
//   2. mesh convergence on TE0: refine until neff changes by less than 2e-3
//   3. sweep the width at the converged mesh size, print a small table

const C0 = 299792458.0;
const LAMBDA = 1.55e-6;
const TOL = 2e-3;

// ---------- 1. model ----------
nyx::project_new(#{ dim: 2, title: "SOI strip width sweep" });

nyx::set_param(#{ name: "w", value: 500e-9 });
nyx::set_param(#{ name: "t", value: 220e-9 });

nyx::add_rect(#{ name: "box",  min: [-2.0e-6, -2.0e-6], max: [2.0e-6, 0.0] });
nyx::add_rect(#{ name: "clad", min: [-2.0e-6,  0.0],    max: [2.0e-6, 2.0e-6] });
nyx::add_rect(#{ name: "core", min: ["-w/2", 0.0],      max: ["w/2", "t"] });   // fields may be parameter expressions
nyx::add_overlay(#{ regions: [
    #{ name: "BOX",  input: "box" },
    #{ name: "Clad", input: "clad" },
    #{ name: "Si",   input: "core" },        // painter order: the core wins
]});

nyx::add_physics(#{
    physics: "ems",
    study: "modes",
    study_params: #{ freq: C0 / LAMBDA, n_target: 2.9, n_modes: 4 },
});
nyx::set_material(#{ physics: "ems", region: "BOX",  material: "silica" });
nyx::set_material(#{ physics: "ems", region: "Clad", material: "silica" });
nyx::set_material(#{ physics: "ems", region: "Si",   material: "silicon" });

// ---------- helpers ----------
// Rhai script functions do NOT see the caller's variables (only constants and their arguments).
// That is a feature here: a function's result depends only on what is passed in.
fn solve_modes() {
    nyx::solve_run(#{ physics: "ems", wait: true });
    let res = nyx::result_info(#{ physics: "ems" })[0];
    if res.stale { throw "result is stale"; }
    let pol = nyx::analyse(#{ analysis: "polarization", physics: "ems" });
    // shape of the polarization rows: { label, value:[re,im] (the effective index), report:{ class: "quasi-TE (x)", te_fraction, ... } }
    let te = ();
    let tm = ();
    for row in pol.rows {
        let cls = row.report.class;
        if cls.starts_with("quasi-TE") && te == () { te = row.value[0]; }
        if cls.starts_with("quasi-TM") && tm == () { tm = row.value[0]; }
    }
    #{ te0: te, tm0: tm, dofs: res.dofs, seconds: res.seconds }
}

// ---------- 2. mesh convergence ----------
let h = 200e-9;
nyx::set_mesh(#{ h: h, adapt_scale: 0.4 });
let prev = solve_modes();
let converged = false;

while h > 25e-9 {                                  // Rhai has while; the host's operation limit stops a runaway loop
    h /= 2.0;
    nyx::set_mesh(#{ h: h, adapt_scale: 0.4 });
    let cur = solve_modes();
    let d = (cur.te0 - prev.te0).abs();
    print(`h = ${fmt(h * 1e9, 1)} nm  TE0 = ${fmt(cur.te0, 6)}  change = ${fmt(d, 6)}  dofs = ${cur.dofs}`);
    prev = cur;
    if d < TOL { converged = true; break; }
}
if !converged { throw `TE0 not converged to ${TOL} above h = 25 nm`; }
let h_conv = h * 2.0;                              // the coarser of the last pair already met the tolerance
nyx::set_mesh(#{ h: h_conv, adapt_scale: 0.4 });

// ---------- 3. width sweep ----------
let rows = [];
for w_nm in range(400, 701, 100) {                  // integers; range(from, to_exclusive, step)
    nyx::set_param(#{ name: "w", value: w_nm.to_float() * 1e-9 });
    try {
        let r = solve_modes();
        rows.push(#{ w_nm: w_nm, te0: r.te0, tm0: r.tm0 });
    } catch (err) {
        rows.push(#{ w_nm: w_nm, error: `${err}` });  // keep the other points
    }
}

// Only what is printed or returned goes back to the agent; fields stay in the engine.
print(`converged h = ${fmt(h_conv * 1e9, 1)} nm`);
print("w_nm  TE0       TM0");
for r in rows {
    if "error" in r { print(`${r.w_nm}  failed: ${r.error}`); }
    else { print(`${r.w_nm}   ${fmt(r.te0, 5)}   ${fmt(r.tm0, 5)}`); }   // no printf in Rhai: the host adds fmt(x, digits), as for Starlark
}
rows                                                // the script's value: returned to the caller as JSON
