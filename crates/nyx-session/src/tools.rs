//! The tools: each one is a single step (geometry, mesh, problem, solve, mode selection, checks) or the recipe that chains them.
//!
//! Every tool takes handles and plain values and returns a small JSON object. Nothing here knows about MCP; `nyx-mcp` only carries the
//! JSON. The physics is `femcore::usecases::soi` (a silicon strip on buried oxide, solved by the full-vector mode solver of `femcore`),
//! used as it is.

use crate::guide;
use crate::params::*;
use crate::workspace::*;
use femcore::meshing::geomdrive;
use femcore::mesh::SimplexMesh;
use femcore::polarization::PolClass;
use femcore::problem::{ProblemSpec, SolveResult, SolveState};
use femcore::progress::Progress;
use femcore::usecases::soi::{slab_neff, Cladding, SlabPol, SoiWire};
use serde_json::{json, Value};

/// Largest mesh this server solves inside its own process; bigger problems belong in a separate solver process (not built yet).
pub const MAX_CELLS: usize = 200_000;
/// Mesh convergence criterion: the repository's own test accepts 1.5e-3 between 50 nm and 30 nm elements.
pub const CONVERGED: f64 = 1.5e-3;
/// Factor on the element sizes for the refined solve of `mesh_converge`.
pub const REFINE: f64 = 0.6;

const TEMPLATES: &[&str] = &["soi_strip"];
const POLS: &[&str] = &["TE", "TM"];
const CLADS: &[&str] = &["air", "oxide"];
const CHECKS: &[&str] = &["slab_bound", "mesh_converge"];

const THICKNESS: Param = Param { name: "thickness", ty: Ty::Num { min: 1e-8, max: 2e-6 }, unit: "m", doc: "silicon device layer thickness (300 nm SOI = 3e-7)", policy: Policy::Required, matters: Matters::High };
const WIDTH: Param = Param {
    name: "width",
    ty: Ty::Num { min: 5e-8, max: 1e-5 },
    unit: "m",
    doc: "strip width at its foot",
    policy: Policy::Ask { json: "5e-7", why: "a common single-mode strip width and the repository's documented default; TE0 neff changes by about 0.9 between 300 nm and 1000 nm wide wires", suggest: "[4.5e-7, 5e-7, 7e-7]" },
    matters: Matters::High,
};
const BOX: Param = Param { name: "box", ty: Ty::Num { min: 1e-7, max: 2e-5 }, unit: "m", doc: "buried oxide thickness", policy: Policy::Default { json: "3e-6", why: "3 um buried oxide as assumed in the request (the repository's code default is 2 um)" }, matters: Matters::Low };
const CLADDING: Param = Param {
    name: "cladding",
    ty: Ty::Enum(CLADS),
    unit: "",
    doc: "material above the strip",
    policy: Policy::Default { json: "\"air\"", why: "no cladding was named; air for a strip 'on silica' (the repository's use case defaults to oxide)" },
    matters: Matters::Medium,
};
const WAVELENGTH: Param = Param { name: "wavelength", ty: Ty::Num { min: 4e-7, max: 1e-5 }, unit: "m", doc: "free-space wavelength", policy: Policy::Default { json: "1.55e-6", why: "telecom C band; silicon 3.4757 and SiO2 1.4440 are the 1550 nm values" }, matters: Matters::Medium };
const POL: Param = Param { name: "pol", ty: Ty::Enum(POLS), unit: "", doc: "mode family wanted", policy: Policy::Default { json: "\"TE\"", why: "TE0 is the usual meaning of 'the fundamental mode'" }, matters: Matters::Medium };

const GEOM_TEMPLATE_P: &[Param] = &[
    Param { name: "template", ty: Ty::Enum(TEMPLATES), unit: "", doc: "soi_strip = rectangular silicon strip on buried oxide", policy: Policy::Required, matters: Matters::High },
    THICKNESS,
    WIDTH,
    BOX,
    CLADDING,
    Param { name: "margin", ty: Ty::Num { min: 1e-7, max: 1e-5 }, unit: "m", doc: "distance of the perfectly conducting window from the core", policy: Policy::Default { json: "1e-6", why: "the repository's setting; the guided modes do not reach it" }, matters: Matters::Low },
];
const MESH_MAKE_P: &[Param] = &[
    Param { name: "geom", ty: Ty::Handle('g'), unit: "", doc: "geometry handle from geom_template", policy: Policy::Required, matters: Matters::High },
    Param { name: "h_core", ty: Ty::Num { min: 2e-9, max: 5e-7 }, unit: "m", doc: "element size in the silicon", policy: Policy::Default { json: "4e-8", why: "the repository's converged setting: TE neff changes by under 1.5e-3 between 50 nm and 30 nm" }, matters: Matters::Medium },
    Param { name: "h_near", ty: Ty::Num { min: 5e-9, max: 1e-6 }, unit: "m", doc: "element size in the oxide next to the core", policy: Policy::Default { json: "9e-8", why: "about twice the core size, as in the repository" }, matters: Matters::Low },
    Param { name: "h_far", ty: Ty::Num { min: 2e-8, max: 2e-6 }, unit: "m", doc: "element size far from the core", policy: Policy::Default { json: "3e-7", why: "the repository's setting" }, matters: Matters::Low },
];
const PROBLEM_DEFINE_P: &[Param] = &[
    Param { name: "mesh", ty: Ty::Handle('m'), unit: "", doc: "mesh handle from mesh_make", policy: Policy::Required, matters: Matters::High },
    WAVELENGTH,
    POL,
    Param { name: "n_modes", ty: Ty::Int { min: 1, max: 6 }, unit: "", doc: "number of modes to compute", policy: Policy::Default { json: "2", why: "the wanted mode and the next one; more modes near the cladding index make the solve slow" }, matters: Matters::Low },
    Param { name: "n_target", ty: Ty::Num { min: 1.0, max: 4.0 }, unit: "", doc: "search index: the solver returns the modes nearest to it. Omit for automatic: just above the analytic slab index of the silicon layer", policy: Policy::Optional, matters: Matters::Low },
];
const SOLVE_RUN_P: &[Param] = &[Param { name: "problem", ty: Ty::Handle('p'), unit: "", doc: "problem handle from problem_define", policy: Policy::Required, matters: Matters::High }];
const RESULT_MODES_P: &[Param] = &[
    Param { name: "result", ty: Ty::Handle('r'), unit: "", doc: "result handle from solve_run", policy: Policy::Required, matters: Matters::High },
    Param { name: "pol", ty: Ty::Enum(POLS), unit: "", doc: "select this family; omit to list all modes", policy: Policy::Optional, matters: Matters::Medium },
    Param { name: "order", ty: Ty::Int { min: 0, max: 9 }, unit: "", doc: "0 = the mode of that family with the highest effective index (TE0, TM0), 1 = the next", policy: Policy::Optional, matters: Matters::Low },
];
const CHECK_RUN_P: &[Param] = &[
    Param { name: "result", ty: Ty::Handle('r'), unit: "", doc: "result handle from solve_run", policy: Policy::Required, matters: Matters::High },
    Param { name: "checks", ty: Ty::List(CHECKS), unit: "", doc: "slab_bound: the wire must lie below the infinite slab of the same thickness. mesh_converge: solve again on a 40 percent finer mesh. Omit for both", policy: Policy::Optional, matters: Matters::Medium },
];
const RECIPE_RUN_P: &[Param] = &[
    Param { name: "recipe", ty: Ty::Enum(&["wg_modes_2d"]), unit: "", doc: "wg_modes_2d = effective index of a mode of a silicon strip on oxide", policy: Policy::Required, matters: Matters::High },
    THICKNESS,
    WIDTH,
    BOX,
    CLADDING,
    WAVELENGTH,
    POL,
    Param { name: "order", ty: Ty::Int { min: 0, max: 9 }, unit: "", doc: "0 = fundamental (TE0 or TM0)", policy: Policy::Default { json: "0", why: "the fundamental mode" }, matters: Matters::Medium },
];
const GET_GUIDE_P: &[Param] = &[Param { name: "topic", ty: Ty::Enum(guide::TOPICS), unit: "", doc: "which card; omit for nyx-instructions", policy: Policy::Optional, matters: Matters::Low }];

pub const TOOLS: &[ToolDef] = &[
    ToolDef { name: "get_guide", title: "Read a nyx guide", description: "Read a short guide. Call with topic nyx-instructions before using other nyx tools. Other topics: optical-modes (reading mode tables), troubleshooting (errors).", params: GET_GUIDE_P, access: Access::Read },
    ToolDef { name: "session_status", title: "Session status", description: "Orientation: session, revision, the handles that exist and what they are. Cheap; call it when a result's nyx.project_hash is not the one you last saw.", params: &[], access: Access::Read },
    ToolDef {
        name: "recipe_run",
        title: "Effective index of a waveguide mode",
        description: "One call from question to answer: effective index of the TE or TM mode of a silicon strip on buried oxide (2D cross-section, full-vector modes, polarization analysis, slab bound and mesh-convergence checks). Takes about a second to a few seconds. If a slot that matters (width) is missing it returns status needs_input with a question for the user. The result has an answer block, the assumptions, the checks and the handles of every step.",
        params: RECIPE_RUN_P,
        access: Access::Expensive,
    },
    ToolDef { name: "geom_template", title: "Build a geometry from a template", description: "Step 1. Build the named cross-section (soi_strip) and return a geometry handle, its region and boundary names and extent. Instant. Unstated values are assumed and listed in 'assumed'.", params: GEOM_TEMPLATE_P, access: Access::SessionWrite },
    ToolDef { name: "mesh_make", title: "Mesh a geometry", description: "Step 2. Triangulate the geometry with element sizes chosen by rule (reasons in 'assumed'). Returns a mesh handle and cell count. Under a second for the strip.", params: MESH_MAKE_P, access: Access::SessionWrite },
    ToolDef { name: "problem_define", title: "Define the optical mode problem", description: "Step 3. Set materials, wavelength and how many modes to find, with the search index placed automatically above the analytic slab index so the fundamental mode is found first. Returns a problem handle. Instant.", params: PROBLEM_DEFINE_P, access: Access::SessionWrite },
    ToolDef { name: "solve_run", title: "Solve the mode problem", description: "Step 4. Run the full-vector mode solver and the polarization analysis. About a second for the strip. Returns a result handle and a compact mode table sorted by effective index.", params: SOLVE_RUN_P, access: Access::Expensive },
    ToolDef { name: "result_modes", title: "Pick a mode from a result", description: "Step 5. Select the TE0, TM0 (or higher order) mode by polarization class and effective index, or list all modes. Cheap. Warns when a higher-index mode of another class exists.", params: RESULT_MODES_P, access: Access::Read },
    ToolDef { name: "check_run", title: "Check a result", description: "Step 6. Numeric checks: slab_bound (cheap) and mesh_converge (one more solve on a finer mesh). Returns pass/fail with the numbers.", params: CHECK_RUN_P, access: Access::Expensive },
];

pub fn find(name: &str) -> Option<&'static ToolDef> {
    TOOLS.iter().find(|t| t.name == name)
}

type Imp = fn(&mut Workspace, &Bound) -> ToolResult;

fn imp(name: &str) -> Option<Imp> {
    Some(match name {
        "get_guide" => get_guide,
        "session_status" => session_status,
        "recipe_run" => recipe_run,
        "geom_template" => geom_template,
        "mesh_make" => mesh_make,
        "problem_define" => problem_define,
        "solve_run" => solve_run,
        "result_modes" => result_modes,
        "check_run" => check_run,
        _ => return None,
    })
}

/// What a call came to.
#[derive(Debug)]
pub enum Outcome {
    Ok(Value),
    NeedsInput(Value),
    Err(ToolError),
    UnknownTool,
}

/// Runs a tool and returns its JSON without the session header (recipes use this for their steps).
pub fn invoke(ws: &mut Workspace, name: &str, args: &Value) -> Result<Value, Stop> {
    let def = find(name).ok_or_else(|| ToolError::new("unknown_tool", format!("no tool '{name}'")).with("valid", TOOLS.iter().map(|t| t.name).collect::<Vec<_>>()))?;
    let b = def.bind(args)?;
    let f = imp(def.name).expect("every tool has an implementation");
    let mut v = f(ws, &b)?;
    if !b.assumed.is_empty() {
        v["assumed"] = Value::Array(b.assumed.iter().map(Assumed::to_json).collect());
    }
    Ok(v)
}

/// The entry point a server calls: runs the tool, adds the session header, writes the journal.
pub fn call(ws: &mut Workspace, name: &str, args: &Value) -> Outcome {
    if find(name).is_none() {
        return Outcome::UnknownTool;
    }
    match invoke(ws, name, args) {
        Ok(mut v) => {
            v["nyx"] = ws.header();
            let handles: Vec<&str> = v.as_object().into_iter().flatten().filter_map(|(_, x)| x.as_str()).filter(|s| s.len() > 3 && s.as_bytes()[1] == b'_').collect();
            let summary = json!({ "handles": handles });
            ws.journal(name, args, "ok", summary);
            Outcome::Ok(v)
        }
        Err(Stop::NeedsInput(mut v)) => {
            v["nyx"] = ws.header();
            ws.journal(name, args, "needs_input", json!({ "slot": v["slot"] }));
            Outcome::NeedsInput(v)
        }
        Err(Stop::Err(e)) => {
            ws.journal(name, args, "error", json!({ "error": e.code }));
            Outcome::Err(e)
        }
    }
}

/// Text for the host: a guide is plain text, everything else compact JSON.
pub fn render(tool: &str, v: &Value) -> String {
    if tool == "get_guide" {
        if let Some(t) = v.get("text").and_then(Value::as_str) {
            return t.to_string();
        }
    }
    v.to_string()
}

fn r4(x: f64) -> f64 {
    (x * 1e4).round() / 1e4
}
fn r6(x: f64) -> f64 {
    (x * 1e6).round() / 1e6
}

fn nm(x: f64) -> f64 {
    r4(x * 1e9)
}

// ------------------------------------------------------------------ steps

fn get_guide(_ws: &mut Workspace, b: &Bound) -> ToolResult {
    let topic = b.opt_str("topic").unwrap_or("nyx-instructions");
    let text = guide::card(topic).ok_or_else(|| ToolError::new("bad_value", format!("no guide '{topic}'")).with("valid", guide::TOPICS))?;
    Ok(json!({ "topic": topic, "text": text }))
}

fn session_status(ws: &mut Workspace, _b: &Bound) -> ToolResult {
    let handles: Vec<Value> = ws.handles().into_iter().map(|(h, k)| json!({ "handle": h, "kind": k })).collect();
    Ok(json!({
        "session": ws.session,
        "folder": ws.dir.display().to_string(),
        "revision": ws.revision(),
        "handles": handles,
        "tools": TOOLS.iter().map(|t| t.name).collect::<Vec<_>>(),
        "journal_written": !ws.journal_failed()
    }))
}

fn make_wire(g: &GeomParams, m: &MeshParams, p: Option<&ProblemParams>) -> SoiWire {
    let mut w = SoiWire {
        width: g.width,
        thickness: g.thickness,
        box_thickness: g.box_thickness,
        cladding: if g.air_cladding { Cladding::Air } else { Cladding::Oxide },
        margin: g.margin,
        h_core: m.h_core,
        h_near: m.h_near,
        h_far: m.h_far,
        ..Default::default()
    };
    if let Some(p) = p {
        w.wavelength = p.wavelength;
        w.n_modes = p.n_modes;
        w.n_target = p.n_target;
    }
    w
}

/// Effective index of the slab of the silicon layer between claddings of the lowest index above and below: an upper bound for any wire of
/// that thickness (for an air-clad wire on oxide it is a loose bound).
fn slab_index(w: &SoiWire, air: bool, fam: Family) -> Option<f64> {
    let n_clad = if air { 1.0 } else { w.n_sio2 };
    slab_neff(w.n_si, n_clad, w.thickness, w.wavelength, if fam == Family::Tm { SlabPol::Tm } else { SlabPol::Te }, 0)
}

fn geom_template(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let gp = GeomParams { thickness: b.num("thickness"), width: b.num("width"), box_thickness: b.num("box"), air_cladding: b.str("cladding") == "air", margin: b.num("margin") };
    let wire = make_wire(&gp, &MeshParams { h_core: 4e-8, h_near: 9e-8, h_far: 3e-7 }, None);
    let ev = geomdrive::evaluate(&wire.model().list).map_err(|e| ToolError::new("geometry_failed", e))?;
    let labels = ev.label_map().map_err(|e| ToolError::new("geometry_failed", e.to_string()))?;
    let (helpers, regions): (Vec<_>, Vec<_>) = labels.domains.iter().map(|d| d.name.clone()).partition(|n| n.ends_with("_near"));
    let mut boundaries: Vec<String> = labels.boundaries.iter().map(|b| b.name.clone()).collect();
    boundaries.sort();
    boundaries.dedup();
    let bb = ev.bounds();
    let key = serde_json::to_string(&gp).unwrap_or_default();
    let (handle, existed) = ws.insert('g', &key, Item::Geom(GeomItem { params: gp.clone() }));
    if !existed {
        ws.write_file(&format!("{handle}.geom.json"), &serde_json::to_string_pretty(&wire.model().list).unwrap_or_default());
    }
    Ok(json!({
        "geom": handle,
        "cached": existed,
        "template": "soi_strip",
        "regions": regions,
        "boundaries": boundaries,
        "mesh_helper_zones_hidden": helpers.len(),
        "extent_nm": { "x": [nm(bb.lo[0]), nm(bb.hi[0])], "z": [nm(bb.lo[2]), nm(bb.hi[2])] },
        "silicon": { "thickness_nm": nm(gp.thickness), "width_nm": nm(gp.width) },
        "model_file": format!("{handle}.geom.json"),
        "next": "mesh_make"
    }))
}

fn mesh_make(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let gh = b.str("geom").to_string();
    let gp = ws.geom(&gh)?.params.clone();
    let mp = MeshParams { h_core: b.num("h_core"), h_near: b.num("h_near"), h_far: b.num("h_far") };
    if mp.h_core > mp.h_near || mp.h_near > mp.h_far {
        return Err(ToolError::new("bad_sizes", "element sizes must satisfy h_core <= h_near <= h_far").with("got", json!({"h_core": mp.h_core, "h_near": mp.h_near, "h_far": mp.h_far})).into());
    }
    let key = format!("{gh}|{}", serde_json::to_string(&mp).unwrap_or_default());
    if let Some(h) = existing(ws, 'm', &key) {
        let m = ws.mesh(&h)?;
        return Ok(mesh_json(&h, true, m));
    }
    let wire = make_wire(&gp, &mp, None);
    let meshed = wire.mesh().map_err(|e| ToolError::new("mesh_failed", e))?;
    if meshed.mesh.n_cells() > MAX_CELLS {
        return Err(ToolError::new("budget_exceeded", format!("{} cells exceed the {} this server solves in-process", meshed.mesh.n_cells(), MAX_CELLS)).with("cells", meshed.mesh.n_cells()).into());
    }
    let (h, _) = ws.insert('m', &key, Item::Mesh(MeshItem { geom: gh, params: mp, meshed }));
    let m = ws.mesh(&h)?;
    Ok(mesh_json(&h, false, m))
}

/// The handle `insert` would give for `key`, if that artifact exists already.
fn existing(ws: &Workspace, prefix: char, key: &str) -> Option<String> {
    let h = fnv64(&format!("{prefix}|{key}"));
    for digits in [8usize, 12, 16] {
        let handle = format!("{prefix}_{:0width$x}", h >> (64 - 4 * digits as u32), width = digits);
        if ws.handles().iter().any(|(k, _)| *k == handle) {
            return Some(handle);
        }
    }
    None
}

fn mesh_json(h: &str, cached: bool, m: &MeshItem) -> Value {
    json!({
        "mesh": h,
        "cached": cached,
        "geom": m.geom,
        "cells": m.meshed.mesh.n_cells(),
        "vertices": m.meshed.mesh.n_points(),
        "sizes_nm": { "core": nm(m.params.h_core), "near": nm(m.params.h_near), "far": nm(m.params.h_far) },
        "next": "problem_define"
    })
}

fn problem_define(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let mh = b.str("mesh").to_string();
    let (gp, mp, meshed_ok) = {
        let m = ws.mesh(&mh)?;
        (ws.geom(&m.geom)?.params.clone(), m.params.clone(), ())
    };
    let _ = meshed_ok;
    let fam = if b.str("pol") == "TM" { Family::Tm } else { Family::Te };
    let mut pp = ProblemParams { wavelength: b.num("wavelength"), n_modes: b.int("n_modes") as usize, n_target: 0.0, auto_target: b.opt_num("n_target").is_none(), pol: b.str("pol").to_string() };
    let probe = make_wire(&gp, &mp, Some(&pp));
    let mut decisions = vec![];
    match b.opt_num("n_target") {
        Some(t) => pp.n_target = t,
        None => {
            let slab = slab_index(&probe, gp.air_cladding, fam).ok_or_else(|| {
                ToolError::new("not_guided", format!("a {} nm silicon layer does not guide a {} slab mode at {} nm", nm(gp.thickness), fam.name(), nm(pp.wavelength)))
                    .with("hint", "increase the thickness or shorten the wavelength")
            })?;
            pp.n_target = slab + 0.02;
            decisions.push(json!({ "what": "search index", "value": r4(pp.n_target), "why": format!("{} slab index {:.4} (upper bound for a wire of this thickness) plus 0.02, so the fundamental mode is the nearest", fam.name(), slab) }));
        }
    }
    decisions.push(json!({ "what": "solver", "value": "full-vector Nedelec(2)+P2 eigenmodes", "why": "lossless real materials, perfectly conducting window, no absorbing layers" }));
    let key = format!("{mh}|{}", serde_json::to_string(&pp).unwrap_or_default());
    if let Some(h) = existing(ws, 'p', &key) {
        return Ok(json!({ "problem": h, "cached": true, "mesh": mh, "next": "solve_run" }));
    }
    let wire = make_wire(&gp, &mp, Some(&pp));
    let spec: ProblemSpec = wire.spec(&ws.mesh(&mh)?.meshed);
    let (h, _) = ws.insert('p', &key, Item::Problem(ProblemItem { mesh: mh.clone(), params: pp.clone(), spec }));
    Ok(json!({
        "problem": h,
        "cached": false,
        "mesh": mh,
        "wavelength_nm": nm(pp.wavelength),
        "n_modes": pp.n_modes,
        "n_target": r4(pp.n_target),
        "materials": { "Si": wire.n_si, "SiO2": wire.n_sio2, "clad": if gp.air_cladding { 1.0 } else { wire.n_sio2 } },
        "decisions": decisions,
        "next": "solve_run"
    }))
}

fn family_of(c: &PolClass) -> Family {
    match c {
        PolClass::QuasiTe => Family::Te,
        PolClass::QuasiTm => Family::Tm,
        _ => Family::Other,
    }
}

fn mode_recs(r: &SolveResult) -> Result<Vec<ModeRec>, ToolError> {
    let a = r.analysis.as_ref().ok_or_else(|| ToolError::new("solver_error", "the mode solver returned no polarization analysis"))?;
    let mut v = vec![];
    for (sol, row) in r.solutions.iter().zip(&a.rows) {
        let (re, im) = sol.neff.ok_or_else(|| ToolError::new("solver_error", "a mode has no effective index"))?;
        let m = match &row.report {
            Some(m) => ModeRec { neff_re: re, neff_im: im, class: m.class.label().to_string(), family: family_of(&m.class), te_fraction: m.te_fraction, long_fraction: m.long_fraction, dop: m.dop, degenerate: row.group.is_some() },
            None => ModeRec { neff_re: re, neff_im: im, class: "no intensity".into(), family: Family::Other, te_fraction: 0.0, long_fraction: 0.0, dop: 0.0, degenerate: row.group.is_some() },
        };
        v.push(m);
    }
    v.sort_by(|x, y| y.neff_re.total_cmp(&x.neff_re));
    Ok(v)
}

fn mode_json(rank: usize, m: &ModeRec) -> Value {
    json!({
        "rank": rank,
        "neff": r6(m.neff_re),
        "neff_im": if m.neff_im.abs() < 1e-9 { 0.0 } else { m.neff_im },
        "class": m.class,
        "family": m.family.name(),
        "te_share": r4(m.te_fraction),
        "longitudinal_share": r4(m.long_fraction),
        "degenerate": m.degenerate
    })
}

struct Solved {
    modes: Vec<ModeRec>,
    seconds: f64,
    dofs: usize,
}

fn run_solve(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<Solved, ToolError> {
    let progress = Progress::new();
    let r = nyx_engine::legacy::solve_with(mesh, spec, &mut SolveState::default(), &progress).map_err(|e| {
        if e == femcore::progress::CANCELLED {
            ToolError::new("cancelled", "the solve was cancelled")
        } else {
            ToolError::new("solver_error", e)
        }
    })?;
    Ok(Solved { modes: mode_recs(&r)?, seconds: r.seconds, dofs: r.dofs })
}

fn solve_run(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let ph = b.str("problem").to_string();
    let key = format!("{ph}|solve");
    if let Some(h) = existing(ws, 'r', &key) {
        let r = ws.result(&h)?;
        return Ok(result_json(&h, true, r));
    }
    let solved = {
        let p = ws.problem(&ph)?;
        let m = ws.mesh(&p.mesh)?;
        run_solve(&m.meshed.mesh, &p.spec)?
    };
    let item = ResultItem { problem: ph.clone(), modes: solved.modes, seconds: solved.seconds, dofs: solved.dofs };
    let digest = json!({ "problem": ph, "seconds": item.seconds, "unknowns": item.dofs, "modes": item.modes });
    let (h, _) = ws.insert('r', &key, Item::Result(item));
    ws.write_file(&format!("{h}.result.json"), &serde_json::to_string_pretty(&digest).unwrap_or_default());
    let r = ws.result(&h)?;
    Ok(result_json(&h, false, r))
}

fn result_json(h: &str, cached: bool, r: &ResultItem) -> Value {
    json!({
        "result": h,
        "cached": cached,
        "problem": r.problem,
        "seconds": r4(r.seconds),
        "unknowns": r.dofs,
        "modes": r.modes.iter().enumerate().map(|(i, m)| mode_json(i, m)).collect::<Vec<_>>(),
        "next": "result_modes"
    })
}

fn select_mode<'a>(modes: &'a [ModeRec], fam: Family, order: usize) -> Result<&'a ModeRec, ToolError> {
    let same: Vec<&ModeRec> = modes.iter().filter(|m| m.family == fam).collect();
    same.get(order).copied().ok_or_else(|| {
        ToolError::new("no_such_mode", format!("the solver returned {} {} mode(s); order {order} does not exist", same.len(), fam.name()))
            .with("found", modes.iter().enumerate().map(|(i, m)| mode_json(i, m)).collect::<Vec<_>>())
            .with("hint", "ask for more modes in problem_define (n_modes), or move n_target")
    })
}

fn result_modes(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let r = ws.result(b.str("result"))?;
    let all: Vec<Value> = r.modes.iter().enumerate().map(|(i, m)| mode_json(i, m)).collect();
    let Some(pol) = b.opt_str("pol") else {
        return Ok(json!({ "result": b.str("result"), "modes": all }));
    };
    let fam = if pol == "TM" { Family::Tm } else { Family::Te };
    let order = b.opt_int("order").unwrap_or(0) as usize;
    let m = select_mode(&r.modes, fam, order)?;
    let mut warnings = vec![];
    if let Some(top) = r.modes.first() {
        if top.neff_re > m.neff_re + 1e-9 && top.family != fam {
            warnings.push(format!("a mode of class '{}' with a higher effective index ({:.5}) exists; check that '{}{}' is what the user means", top.class, top.neff_re, fam.name(), order));
        }
    }
    if m.degenerate {
        warnings.push("this mode is part of a degenerate group; its class and field are an arbitrary mixture".to_string());
    }
    Ok(json!({
        "result": b.str("result"),
        "selected": { "name": format!("{}{}", fam.name(), order), "mode": mode_json(0, m) },
        "warnings": warnings,
        "modes": all
    }))
}

fn check_run(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let res = ws.result(b.str("result"))?;
    let prob = ws.problem(&res.problem)?;
    let mesh = ws.mesh(&prob.mesh)?;
    let gp = ws.geom(&mesh.geom)?.params.clone();
    let wire = make_wire(&gp, &mesh.params, Some(&prob.params));
    let wanted: Vec<String> = if b.get("checks").is_some() { b.list("checks") } else { CHECKS.iter().map(|s| s.to_string()).collect() };
    let mut out = vec![];
    for c in &wanted {
        match c.as_str() {
            "slab_bound" => {
                let mut rows = vec![];
                let mut pass = true;
                for fam in [Family::Te, Family::Tm] {
                    let Some(m) = res.modes.iter().find(|m| m.family == fam) else { continue };
                    let Some(slab) = slab_index(&wire, gp.air_cladding, fam) else { continue };
                    let ok = m.neff_re < slab;
                    pass &= ok;
                    rows.push(json!({ "family": fam.name(), "neff": r6(m.neff_re), "slab": r6(slab), "below_slab": ok }));
                }
                out.push(json!({ "check": "slab_bound", "pass": pass, "rows": rows, "note": "a wire cannot exceed the infinite slab of the same thickness; for an air-clad wire the slab is taken with air on both sides, a loose bound" }));
            }
            "mesh_converge" => {
                let mp2 = MeshParams { h_core: mesh.params.h_core * REFINE, h_near: mesh.params.h_near * REFINE, h_far: mesh.params.h_far };
                let wire2 = make_wire(&gp, &mp2, Some(&prob.params));
                let g2 = wire2.mesh().map_err(|e| ToolError::new("mesh_failed", e))?;
                if g2.mesh.n_cells() > MAX_CELLS {
                    return Err(ToolError::new("budget_exceeded", "the refined mesh is too large to solve in-process").with("cells", g2.mesh.n_cells()).into());
                }
                let spec2 = wire2.spec(&g2);
                let s2 = run_solve(&g2.mesh, &spec2)?;
                let mut rows = vec![];
                let mut pass = true;
                for fam in [Family::Te, Family::Tm] {
                    let (Some(a), Some(c)) = (res.modes.iter().find(|m| m.family == fam), s2.modes.iter().find(|m| m.family == fam)) else { continue };
                    let delta = (a.neff_re - c.neff_re).abs();
                    let ok = delta < CONVERGED;
                    pass &= ok;
                    rows.push(json!({ "family": fam.name(), "neff": r6(a.neff_re), "neff_refined": r6(c.neff_re), "delta": delta, "converged": ok }));
                }
                out.push(json!({
                    "check": "mesh_converge",
                    "pass": pass && !rows.is_empty(),
                    "rows": rows,
                    "core_size_nm": [nm(mesh.params.h_core), nm(mp2.h_core)],
                    "criterion": format!("delta < {CONVERGED} (the repository's own mesh test accepts this between 50 nm and 30 nm elements)"),
                    "refined_cells": g2.mesh.n_cells(),
                    "seconds": r4(s2.seconds)
                }));
            }
            _ => {}
        }
    }
    Ok(json!({ "result": b.str("result"), "checks": out, "all_passed": out.iter().all(|c| c["pass"] == json!(true)) }))
}

fn recipe_run(ws: &mut Workspace, b: &Bound) -> ToolResult {
    let mut ledger: Vec<Value> = b.assumed.iter().map(Assumed::to_json).collect();
    let note = |v: &Value, ledger: &mut Vec<Value>| {
        if let Some(a) = v.get("assumed").and_then(Value::as_array) {
            ledger.extend(a.iter().cloned());
        }
    };
    let (t, w) = (b.num("thickness"), b.num("width"));
    let clad = b.str("cladding").to_string();
    let g = invoke(ws, "geom_template", &json!({ "template": "soi_strip", "thickness": t, "width": w, "box": b.num("box"), "cladding": clad }))?;
    note(&g, &mut ledger);
    let m = invoke(ws, "mesh_make", &json!({ "geom": g["geom"] }))?;
    note(&m, &mut ledger);
    let p = invoke(ws, "problem_define", &json!({ "mesh": m["mesh"], "wavelength": b.num("wavelength"), "pol": b.str("pol") }))?;
    note(&p, &mut ledger);
    let r = invoke(ws, "solve_run", &json!({ "problem": p["problem"] }))?;
    let rh = r["result"].as_str().unwrap_or_default().to_string();
    let sel = invoke(ws, "result_modes", &json!({ "result": rh, "pol": b.str("pol"), "order": b.int("order") }))?;
    let chk = invoke(ws, "check_run", &json!({ "result": rh }))?;

    let mode = &sel["selected"]["mode"];
    let name = sel["selected"]["name"].as_str().unwrap_or("mode");
    let neff = mode["neff"].as_f64().unwrap_or(f64::NAN);
    let clad_txt = if clad == "air" { "air" } else { "SiO2" };
    let mut warnings: Vec<String> = sel["warnings"].as_array().into_iter().flatten().filter_map(|x| x.as_str().map(String::from)).collect();
    let mut checked = vec![];
    for c in chk["checks"].as_array().into_iter().flatten() {
        let pass = c["pass"] == json!(true);
        let line = match c["check"].as_str() {
            Some("slab_bound") => {
                let row = c["rows"].as_array().and_then(|r| r.iter().find(|x| x["family"] == mode["family"]).cloned());
                match row {
                    Some(r) => format!("below the {} infinite-slab index {:.4}: {}", r["family"].as_str().unwrap_or(""), r["slab"].as_f64().unwrap_or(f64::NAN), if pass { "yes" } else { "NO" }),
                    None => "slab bound: nothing to compare".to_string(),
                }
            }
            Some("mesh_converge") => {
                let row = c["rows"].as_array().and_then(|r| r.iter().find(|x| x["family"] == mode["family"]).cloned());
                match row {
                    Some(r) => format!("changed by {:.1e} when the core elements went from {} to {} nm: {}", r["delta"].as_f64().unwrap_or(f64::NAN), c["core_size_nm"][0], c["core_size_nm"][1], if pass { "converged" } else { "NOT converged" }),
                    None => "mesh convergence: no mode to compare".to_string(),
                }
            }
            _ => String::new(),
        };
        if !pass {
            warnings.push(format!("check failed: {line}"));
        }
        checked.push(line);
    }
    let assumed_txt: Vec<String> = ledger.iter().map(|a| format!("{} = {} ({})", a["name"].as_str().unwrap_or(""), a["value"], a["why"].as_str().unwrap_or(""))).collect();
    let failed = chk["all_passed"] != json!(true);
    let answer = format!(
        "{}{name} effective index = {neff:.5} ({}, TE share {:.3}) for a {:.0} nm x {:.0} nm silicon strip on {:.2} um SiO2, {:.0} nm wavelength, {clad_txt} cladding.\nAssumed: {}.\nChecked: {}.\nSolver: {} unknowns, {} cells, {:.2} s. Perfectly conducting window {:.1} um from the core; no silicon substrate under the oxide.",
        if failed { "UNVERIFIED (a check failed, see 'Checked'; do not quote this number without saying so): " } else { "" },
        mode["class"].as_str().unwrap_or(""),
        mode["te_share"].as_f64().unwrap_or(f64::NAN),
        nm(t),
        nm(w),
        b.num("box") * 1e6,
        nm(b.num("wavelength")),
        if assumed_txt.is_empty() { "nothing (every value was given)".to_string() } else { assumed_txt.join("; ") },
        checked.join("; "),
        r["unknowns"],
        m["cells"],
        r["seconds"].as_f64().unwrap_or(f64::NAN),
        1.0
    );
    Ok(json!({
        "status": if failed { "unverified" } else { "done" },
        "answer": answer,
        "result": { "name": name, "neff": mode["neff"], "neff_im": mode["neff_im"], "class": mode["class"], "te_share": mode["te_share"], "longitudinal_share": mode["longitudinal_share"] },
        "assumed": ledger,
        "checks": chk["checks"],
        "warnings": warnings,
        "handles": { "geom": g["geom"], "mesh": m["mesh"], "problem": p["problem"], "result": rh },
        "all_modes": sel["modes"]
    }))
}
