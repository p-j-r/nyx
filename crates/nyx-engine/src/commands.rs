// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The commands: every action a front-end can take, by name, with JSON arguments, a JSON schema and a JSON reply. The window calls the typed
//! methods of `Engine` directly or these commands; the MCP server turns `catalog()` into its tool list and `call()` into its tool calls,
//! so a tool is a command and nothing else. The geometry tools are generated from `geomdesc` (gate G5: every operation, none picked by
//! hand), the physics-dependent schemas come from the registry.

use crate::engine::Engine;
use crate::geomdesc as gd;
use crate::geometry;
use crate::project::{MeshSettings, Project};
use crate::queue::QueueOpts;
use nyx_physics::{contract, BcSpec, Error, MaterialAssign, ProblemSpec, StudySpec};
use serde::Serialize;
use serde_json::{json, Map, Value};

/// How much a command changes or costs; a host may ask the person before an expensive or editing command.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Access {
    /// reads the project or the results
    Read,
    /// changes the project
    Edit,
    /// starts a mesh build or a solve, or writes files
    Expensive,
}

#[derive(Clone, Debug, Serialize)]
pub struct CommandDesc {
    pub name: String,
    pub title: String,
    pub description: String,
    pub access: Access,
    /// JSON Schema of the arguments
    pub schema: Value,
}

type R = Result<Value, Error>;

// ------------------------------------------------------------------------------------------------------------- schema helpers

fn text(doc: &str) -> Value {
    json!({ "type": "string", "description": doc })
}
fn number(doc: &str) -> Value {
    json!({ "type": "number", "description": doc })
}
fn int(doc: &str) -> Value {
    json!({ "type": "integer", "description": doc })
}
fn boolean(doc: &str) -> Value {
    json!({ "type": "boolean", "description": doc })
}
/// An expression of fields (ADR 90): a field name, a number, a handle or an operation object; `overlap` and `confinement` also take a list of them.
fn expr_schema(doc: &str) -> Value {
    json!({ "type": ["string", "number", "object", "array"], "description": doc })
}
fn object(doc: &str) -> Value {
    json!({ "type": "object", "description": doc })
}
fn post_analyses_schema() -> Value {
    json!({
        "type": "array",
        "description": "read-only post-solve analyses or field reductions run on full in-memory fields before fields are dropped (either { \"analysis\": \"...\", \"physics\"?: \"...\", \"params\"?: {...} } shorthand or { \"cmd\": \"analyse\" | \"integrate\" | \"confinement\" | \"extremum\" | \"probe\" | \"line_cut\" | \"overlap\", \"args\": {...} })",
        "items": {
            "type": "object",
            "description": "one read-only post-solve analysis or field-reduction request",
            "properties": {
                "analysis": text("shorthand analysis id to run via `analyse` (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator, uq_cqed)"),
                "physics": text("three-letter physics code for shorthand `analysis`; defaults to the solved physics that supports the analysis"),
                "params": object("optional parameters object for shorthand `analysis`"),
                "cmd": text("explicit read-only post-solve command name (analyse, integrate, confinement, extremum, probe, line_cut, overlap, result_info, result_field, list_analyses, field, field_op)"),
                "args": object("arguments object passed to `cmd`"),
                "label": text("optional label attached to this post-analysis entry in the reply")
            },
            "additionalProperties": false
        }
    })
}
fn schema(props: Vec<(&str, Value)>, required: &[&str]) -> Value {
    let p: Map<String, Value> = props.into_iter().map(|(k, v)| (k.to_string(), v)).collect();
    json!({ "type": "object", "properties": Value::Object(p), "required": required, "additionalProperties": false })
}

fn cmd(name: &str, title: &str, description: &str, access: Access, schema: Value) -> CommandDesc {
    CommandDesc { name: name.into(), title: title.into(), description: description.into(), access, schema }
}

// ------------------------------------------------------------------------------------------------------------- argument helpers

fn args_obj(a: &Value) -> Result<&Map<String, Value>, Error> {
    static EMPTY: std::sync::OnceLock<Map<String, Value>> = std::sync::OnceLock::new();
    match a {
        Value::Null => Ok(EMPTY.get_or_init(Map::new)),
        Value::Object(m) => Ok(m),
        other => Err(Error::invalid(format!("the arguments must be an object, not {other}"))),
    }
}

fn get_str<'a>(a: &'a Value, k: &str) -> Result<&'a str, Error> {
    args_obj(a)?
        .get(k)
        .and_then(Value::as_str)
        .ok_or_else(|| Error::invalid(format!("`{k}` is required and must be text")))
}

fn opt_str<'a>(a: &'a Value, k: &str) -> Option<&'a str> {
    a.get(k).and_then(Value::as_str)
}

fn opt_bool(a: &Value, k: &str, default: bool) -> bool {
    a.get(k).and_then(Value::as_bool).unwrap_or(default)
}

fn get_u64(a: &Value, k: &str) -> Result<u64, Error> {
    a.get(k)
        .and_then(Value::as_u64)
        .ok_or_else(|| Error::invalid(format!("`{k}` is required and must be a whole number")))
}

fn ok() -> R {
    Ok(json!({ "ok": true }))
}

// ------------------------------------------------------------------------------------------------------------- the catalog

impl Engine {
    /// Every command of this build, with its schema. The geometry `add_*` commands come from `geomdesc`; the choices of physics from
    /// the registry.
    pub fn catalog(&self) -> Vec<CommandDesc> {
        use Access::*;
        let ids: Vec<&str> = self.registry.ids();
        let phys = || json!({ "type": "string", "enum": ids, "description": "the three-letter code of the physics" });
        let wl = || json!({ "type": "number", "description": "metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve" });
        let mut c = vec![
            cmd("list_commands", "List the commands", "Every command with its title, access and schema.", Read, schema(vec![], &[])),
            cmd("project_new", "New project", "Starts an empty project (dropping the current one, its mesh and results).", Edit,
                schema(vec![("title", text("project title")), ("dim", int("2 (x-z cross-section, the default) or 3")), ("h", number("global element size in metres"))], &[])),
            cmd("project_info", "Project summary", "Title, geometry size, the problem of each physics, the mesh, the results (and whether they are stale), the store binding and the jobs.", Read, schema(vec![], &[])),
            cmd("project_json", "Project as JSON", "The whole project as JSON text (the file that is saved).", Read, schema(vec![], &[])),
            cmd("project_script", "Project as a script", "The project as a script (Python, `.nys`) that rebuilds it: project_new, set_param (each parameter reads `params`), one add_* per feature, set_mesh, add_material, set_problem, set_active. It does not solve. It is replayed on a scratch engine and offered only when it gives the same project; `form` says how it is written: `commands`, `geometry_json` (the geometry as one set_geometry, for ids and pins that must be kept) or `project_json` (one project_load of the whole project: exact, not readable). Starts with project_new, so over MCP it needs allow_destructive.", Read, schema(vec![], &[])),
            cmd("project_load", "Load a project from JSON", "Replaces the project by JSON text, after checking it.", Edit, schema(vec![("json", text("project JSON text"))], &["json"])),
            cmd("project_open", "Open a project file", "Reads a project JSON file.", Edit, schema(vec![("path", text("file path"))], &["path"])),
            cmd("project_save", "Save the project to a file", "Writes the project as a JSON file.", Expensive, schema(vec![("path", text("file path"))], &["path"])),
            cmd("list_operations", "List the geometry operations", "Every geometry operation with its fields, types and defaults (what the add_* commands take).", Read, schema(vec![], &[])),
            cmd("list_features", "List the features", "The feature list: parameters, and each feature with its kind and arguments.", Read, schema(vec![], &[])),
            cmd("set_geometry", "Replace the geometry", "Replaces the whole feature list by one given as JSON (what the window's geometry editor writes); it is checked first.", Edit,
                schema(vec![("geometry", object("a feature list: { dim, params, features }"))], &["geometry"])),
            cmd("set_problem", "Replace the problem of a physics", "Replaces the whole problem (study, order, materials, conditions, options) of a physics by one given as JSON, after checking it against the descriptors. Adds the physics if the project has no problem for it yet.", Edit,
                schema(vec![("physics", phys()), ("problem", object("a problem: { physics, order, study, materials, bcs, params }"))], &["physics", "problem"])),
            cmd("update_feature", "Change a feature", "Changes some fields of a feature; give only the fields to change, in the form the add command takes.", Edit,
                schema(vec![("name", text("feature name")), ("fields", object("fields to change"))], &["name", "fields"])),
            cmd("remove_feature", "Remove a feature", "Removes a feature that no later feature uses.", Edit, schema(vec![("name", text("feature name"))], &["name"])),
            cmd("rename_feature", "Rename a feature", "", Edit, schema(vec![("name", text("feature name")), ("new_name", text("new name"))], &["name", "new_name"])),
            cmd("suppress_feature", "Suppress or restore a feature", "A suppressed feature stays in the list and is skipped.", Edit,
                schema(vec![("name", text("feature name")), ("suppressed", boolean("true (default) to suppress, false to restore"))], &["name"])),
            cmd("tag_feature", "Tag a feature", "Adds (or removes) a tag of a feature.", Edit,
                schema(vec![("name", text("feature name")), ("tag", text("tag")), ("remove", boolean("true removes the tag"))], &["name", "tag"])),
            cmd("set_param", "Set a model parameter", "Sets (or creates) a named parameter that expressions in the features can use. A number, or an expression over other parameters.", Edit,
                schema(vec![("name", text("parameter name")), ("value", json!({ "oneOf": [{ "type": "number" }, { "type": "string" }], "description": "a number in metres (or the unit the features give it), or an expression over other parameters" }))], &["name", "value"])),
            cmd("remove_param", "Remove a model parameter", "", Edit, schema(vec![("name", text("parameter name"))], &["name"])),
            cmd("geometry_names", "Evaluate the geometry and list its names", "Evaluates the model and returns the region and boundary names a problem can use, and the bounding box. Instant; it reports a geometry that does not evaluate.", Read, schema(vec![], &[])),
            cmd("list_physics", "List the physics of this build", "The physics compiled in, their studies, and which the project loads.", Read, schema(vec![], &[])),
            cmd("describe_physics", "Describe a physics", "Studies, conditions, properties, fields and couplings of one physics, and the JSON Schema of its problem.", Read, schema(vec![("physics", phys())], &["physics"])),
            cmd("add_physics", "Add a physics to the project", "Adds a problem for a physics (a project loads one to three). The study defaults to the first the physics lists.", Edit,
                schema(vec![("physics", phys()), ("study", text("study kind (see describe_physics)")), ("study_params", object("study parameters")), ("order", int("element order, default 1"))], &["physics"])),
            cmd("set_study", "Set the study of a physics", "", Edit,
                schema(vec![("physics", phys()), ("kind", text("study kind")), ("params", object("study parameters")), ("order", int("element order"))], &["physics", "kind"])),
            cmd("set_material", "Assign a material to a region", "A material by name (library or project: list_materials), a named material placed at an angle (`{\"material\": \"sapphire\", \"euler_deg\": [0, 57.6, 0]}`, Bunge angles in degrees), or inline scalar/tensor properties by descriptor name (`{\"E\": 1.7e11, \"nu\": 0.28}`, diagonal `{\"eps_r\": [9.4, 9.4, 11.6]}`, symmetric 6-vector `{\"eps_r\": [xx, yy, zz, xy, xz, yz]}`, or symmetric 3x3 matrix `{\"eps_r\": [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]}`).", Edit,
                schema(vec![("physics", phys()), ("region", text("region name")), ("material", json!({ "oneOf": [{ "type": "string" }, { "type": "object" }], "description": "material name, { material, euler_deg }, or inline scalar/tensor property map (scalar, [xx, yy, zz], [xx, yy, zz, xy, xz, yz], or symmetric 3x3 matrix)" }))], &["physics", "region", "material"])),
            cmd("add_condition", "Add a boundary or domain condition", "`on` is a selector over the geometry (`outer & plane(x = 0)`), a boundary name or a region name.", Edit,
                schema(vec![("physics", phys()), ("on", text("selector, boundary or region")), ("kind", text("condition kind (see describe_physics)")), ("params", object("condition parameters"))], &["physics", "on", "kind"])),
            cmd("remove_condition", "Remove a condition", "By its position (from 0) in the physics' list, as project_info shows.", Edit, schema(vec![("physics", phys()), ("index", int("position from 0"))], &["physics", "index"])),
            cmd("set_physics_params", "Set the options of a physics", "", Edit, schema(vec![("physics", phys()), ("params", object("physics options"))], &["physics", "params"])),
            cmd("remove_physics", "Remove a physics from the project", "", Edit, schema(vec![("physics", phys())], &["physics"])),
            cmd("set_active", "Set the active physics", "The one solved by default.", Edit, schema(vec![("physics", phys())], &["physics"])),
            cmd("check_project", "Check the project", "Checks the geometry, the element size and every problem against the descriptors of its physics, without meshing.", Read, schema(vec![], &[])),
            cmd("import_gds", "Import GDSII conductive sheets", "Reads a binary GDSII (.gds) stream in a 3D project, flattens SREF/AREF hierarchy and PATH/BOUNDARY/BOX geometry into SI metres, unions each requested [layer, datatype] pair, optionally splits disconnected conductive islands and renames islands from TEXT labels on `pin_layer`, fits circular arcs within `arc_tol`, and adds `Polygons` + `Sheet` features at height `z`. Returns `{ sheets, layers, bbox, unpadded_bbox, padded_bbox, lateral_pad, pins }`.", Edit,
                schema(vec![
                    ("path", text("file path to the binary GDSII (.gds) stream")),
                    ("layers", object("map from base sheet name to [layer, datatype] GDSII integer pair, e.g. {\"gnd\": [1, 0], \"sig\": [1, 10]}")),
                    ("z", number("wafer-normal z height of the imported conductive sheets in metres, default 0")),
                    ("split_islands", boolean("split disconnected polygons on each layer into separate sheets (<name>_0, <name>_1, ... or matched TEXT pin names), default true")),
                    ("pin_layer", json!({ "type": "array", "items": { "type": "integer" }, "minItems": 2, "maxItems": 2, "description": "[layer, texttype] GDSII pair whose TEXT labels rename matching conductive islands" })),
                    ("simplify_tol", number("collinear/chord vertex simplification tolerance in metres (shared cross-layer contact vertices stay pinned), default 0")),
                    ("arc_tol", number("1-DOF endpoint-anchored circular arc fitting tolerance in metres (0 disables arc fitting), default 0")),
                    ("cell", text("top structure cell name to flatten; omit to auto-select the unreferenced root cell")),
                    ("layer_bbox", boolean("when true (default), return bbox over imported layers only; when false, return full top-cell bbox")),
                    ("lateral_pad", number("non-negative lateral buffer margin in metres added around the returned bbox, default 0")),
                ], &["path", "layers"])),
            cmd("add_chip_stack", "Add a 3D chip dielectric stack", "Constructs a 3D chip dielectric stack (`VacBot`, `Si`, `VacTop` boxes and `Overlay`) over a 2D `[xmin, ymin, xmax, ymax]` bounding box (explicitly given or inferred from existing 2D/sheet features in the project), expanded laterally by `lateral_pad >= 0`.", Edit,
                schema(vec![
                    ("bbox", json!({ "type": "array", "items": { "type": "number" }, "minItems": 4, "maxItems": 4, "description": "optional unpadded 2D bounding box [xmin, ymin, xmax, ymax] in metres; omit to infer from existing 2D/sheet features" })),
                    ("lateral_pad", number("non-negative lateral buffer margin in metres added on all 4 horizontal sides, default 0")),
                    ("z", number("wafer-normal z coordinate of the top substrate interface in metres, default 0")),
                    ("t_sub", number("substrate thickness in metres (alias: substrate_thickness), default 350e-6")),
                    ("substrate_thickness", number("alias for t_sub: substrate thickness in metres")),
                    ("h_vac_top", number("top vacuum thickness above z in metres (aliases: t_vac_top, vac_top_thickness), default 3 * t_sub")),
                    ("t_vac_top", number("alias for h_vac_top: top vacuum thickness in metres")),
                    ("vac_top_thickness", number("alias for h_vac_top: top vacuum thickness in metres")),
                    ("h_vac_bot", number("bottom vacuum thickness below the substrate in metres (0 omits bottom vacuum; aliases: t_vac_bot, vac_bot_thickness), default h_vac_top")),
                    ("t_vac_bot", number("alias for h_vac_bot: bottom vacuum thickness in metres")),
                    ("vac_bot_thickness", number("alias for h_vac_bot: bottom vacuum thickness in metres")),
                    ("vac_bot_box", text("feature name of the bottom vacuum Box, default \"vac_bot\"")),
                    ("substrate_box", text("feature name of the substrate Box, default \"si_sub\"")),
                    ("vac_top_box", text("feature name of the top vacuum Box, default \"vac_top\"")),
                    ("vac_bot_region", text("region name of the bottom vacuum in Overlay, default \"VacBot\"")),
                    ("substrate_region", text("region name of the substrate in Overlay, default \"Si\"")),
                    ("vac_top_region", text("region name of the top vacuum in Overlay, default \"VacTop\"")),
                    ("overlay", boolean("create or update the project's Overlay feature with the stack regions, default true")),
                ], &[])),
            cmd("set_mesh", "Set the mesh settings", "Give only what changes. `adapt_scale`, `adapt_grade`, `adapt_max_cells`, `gap_cells`, `max_cells`, `dz_min`, `max_z_slabs`, `z_growth_ratio`, `max_steiner_multiplier`, `singular_edges`, `quarter_point_xi`, `pml_mode`, `structured_pml`, or `sweep` set to null return to their defaults. `sweep` extrudes a 2D cross-section or a 3D geometry model along y through layered intervals with optional `fill` or per-layer `regions` remap.", Edit,
                schema(vec![("h", number("global element size in metres")), ("adapt_scale", json!({ "type": ["number", "null"], "description": "adaptive meshing: the factor on the sizes the model asks for (smaller is finer); null is the simple mesher" })), ("adapt_grade", json!({ "type": ["number", "null"], "description": "adaptive / AMR metric gradation slope (> 0); null uses 0.4" })), ("adapt_max_cells", json!({ "type": ["integer", "null"], "minimum": 1, "description": "adaptive / AMR maximum cell count (>= 1); null uses 2000000" })), ("physics_sizing", boolean("size the elements from the waves of the physics (points per wavelength)")), ("points_per_wavelength", json!({ "type": ["number", "null"], "description": "points per shortest wavelength when physics_sizing is on; null: 12 for order 1, 6 for order 2" })), ("pml_mode", json!({ "type": ["string", "null"], "enum": ["auto", "structured", "inline", null], "description": "how absorbing layers are meshed: \"auto\" (default, structured native Quad/Prism/Hex layers when PML regions exist), \"structured\", or \"inline\" (legacy unstructured simplices)" })), ("structured_pml", json!({ "type": ["boolean", "null"], "description": "legacy boolean alias for pml_mode (true -> \"structured\", false -> \"inline\")" })), ("max_dev", json!({ "type": ["number", "null"], "description": "largest chord deviation from a curved edge or face, metres; null: h / 50 for the adaptive mesher and no limit for the simple 3D mesher" })), ("curved", boolean("3D: curved (quadratic) isoparametric geometry for cells on curved surfaces and fitted sheet arc rims (ems, dfs, and sheet3d), default false")),
                    ("gap_cells", json!({ "type": ["number", "null"], "description": "3D sheet stack: target minimum number of 2D triangles across narrow coplanar gaps between conductor sheets (positive number, e.g. 2.5); null disables automatic gap clamping" })),
                    ("z_coarsen", boolean("3D sheet stack: progressively coarsen the 2D sheet triangulation away from the conductor plane along z via conforming pure-tetrahedral transition slabs, default false")),
                    ("max_cells", json!({ "type": ["integer", "null"], "minimum": 1, "description": "maximum number of mesh cells allowed before aborting meshing (overrides the default 2,000,000 cap on adaptive meshing and enforces pre-allocation guard on 3D sheet stacks); null restores default" })),
                    ("dz_min", json!({ "type": ["number", "null"], "exclusiveMinimum": 0, "description": "3D sheet stack: floor on initial vertical z-slab thickness near sheet planes in metres; null disables the floor" })),
                    ("max_z_slabs", json!({ "type": ["integer", "null"], "minimum": 1, "description": "3D sheet stack: maximum total number of horizontal z-slabs across active stack layers (must be >= active stack layers); null disables the cap" })),
                    ("z_growth_ratio", json!({ "type": ["number", "null"], "exclusiveMinimum": 1.0, "description": "3D sheet stack: geometric growth ratio (> 1.0) for vertical z-slab thickness away from sheet planes; null uses default 1.35" })),
                    ("max_steiner_multiplier", json!({ "type": ["integer", "null"], "minimum": 0, "description": "2D Ruppert refinement: cap on additional Steiner vertices as a multiple of initial PSLG vertices; null uses default max(20 * n_init, 100000)" })),
                    ("strict_mesh", boolean("fail meshing if 2D Ruppert refinement terminates early before completion or drops unhealed collinear sliver triangles, default true")),
                    ("singular_edges", json!({ "type": ["string", "null"], "enum": ["none", "quarter_point", null], "description": "3D sheet stack: edge singularity enrichment on zero-thickness conductor rims (\"none\" or \"quarter_point\"); null uses study defaults" })),
                    ("quarter_point_xi", json!({ "type": ["number", "null"], "minimum": 0.25, "maximum": 0.5, "description": "3D sheet stack: Barsoum quarter-point mid-edge parameter xi_QP in [0.25, 0.50]; null uses 0.25" })),
                    ("corner_grading", boolean("3D sheet stack: radical Morrison-Lewis corner grading toward sharp 2D sheet corners, default false")),
                    ("target_tol", json!({ "type": ["number", "null"], "exclusiveMinimum": 0, "description": "target relative discretization tolerance used to certify error bounds (e.g. Prager-Synge flux_bounds) before a convergence check; null uses default 1e-2" })),
                    ("adapt_passes", json!({ "type": ["integer", "null"], "minimum": 1, "description": "number of outer goal-oriented AMR passes for ems::cavity (Maxwell residual) and sms (Hessian estimator); null or 1 disables AMR" })),
                    ("adapt_rel_tol", json!({ "type": ["number", "null"], "exclusiveMinimum": 0, "description": "relative quantity-of-interest tolerance between consecutive ems::cavity / sms AMR passes; null uses 1e-3" })),
                    ("adapt_growth", json!({ "type": ["number", "null"], "exclusiveMinimum": 1.0, "description": "target mesh cell growth multiplier per ems::cavity / sms AMR pass (> 1.0); null uses 1.5" })),
                    ("retain_amr_data", boolean("retain transient cellwise ElectroAmrData indicators on SolveResult instead of stripping them after solve, default false")),
                    ("sweep", json!({ "type": ["array", "null"], "items": { "type": "object", "properties": { "thickness": { "type": "number", "description": "thickness of the layer along y, metres" }, "n": { "type": "integer", "description": "cells across the thickness" }, "domain_offset": { "type": "integer", "description": "added to the domain numbers of the layer's cells" }, "fill": { "type": "string", "description": "region name: every cell of the layer belongs to it" }, "regions": { "type": "object", "description": "optional per-layer region remap { source_region: target_region } in this layer" } }, "required": ["thickness", "n"] }, "description": "layers to sweep a 2D cross-section or 3D geometry model along y" }))], &[])),
            cmd("mesh_build", "Mesh the geometry", "Meshes the geometry with the current settings. With wait false it starts a job and returns its id.", Expensive, schema(vec![("wait", boolean("wait for the mesh, default true"))], &[])),
            cmd("solve_run", "Solve", "Solves the active physics (or the named ones, in coupling order), meshing first when needed. With wait false it starts a job. When the project is bound to a workflow the run is recorded.", Expensive,
                schema(vec![("physics", json!({ "oneOf": [{ "type": "string" }, { "type": "array", "items": { "type": "string" } }], "description": "physics code(s); default the active one" })), ("wait", boolean("wait for the result, default true")), ("isolate", boolean("run the solve in a worker process: running out of memory or crashing fails this job and not the program, and a cancel stops it at once; the cost is a results round trip. Default: what the program was set up to do")), ("shared", boolean("hand the solve to the queue of the workspace, which every program that uses the workspace shares (the window, the MCP server, `nyx queue`): a queue worker process runs the items one at a time in a worker, and keeps running when the program that asked is closed. Default: what the program was set up to do; implies a worker process")), ("force", boolean("start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits")), ("save_fields", boolean("when false, omit mesh and nodal/cell field arrays from worker HDF5 I/O and stored results after running any post_analyses, keeping scalars and eigenvalues; default true")), ("retain_amr_data", boolean("retain transient cellwise ElectroAmrData indicators on the stored result instead of stripping them after solve; default false")), ("post_analyses", post_analyses_schema()), ("memory_wait_timeout_s", number("maximum seconds (>= 0) to wait when free RAM is below the solve estimate (overrides NYX_MEMORY_WAIT_TIMEOUT_S and engine default)"))], &[])),
            cmd("job_status", "Status of a job", "State, current stage, fraction done and elapsed time.", Read, schema(vec![("id", int("job id"))], &["id"])),
            cmd("job_list", "List the jobs", "", Read, schema(vec![], &[])),
            cmd("job_cancel", "Cancel a job", "The solver notices at its next check; a factorisation in progress finishes first.", Expensive, schema(vec![("id", int("job id"))], &["id"])),
            cmd("job_wait", "Wait for a job", "Blocks until the job has finished and applies its outcome. With `timeout_s` it gives up after that many seconds and returns the job's entry with `waiting: true` (the job goes on; wait again).", Expensive,
                schema(vec![("id", int("job id")), ("timeout_s", number("longest wait in seconds; omit to wait until the job has finished"))], &["id"])),
            cmd("queue_submit", "Submit work to the queue", "Queues a `solve` (default; `physics` is a code or a list), a `mesh` build or a mesh `convergence` check (`physics`, `factor`, `grids`, `tolerance`) and returns its entry at once. The heavy lane runs one job at a time in order of priority (1 a person's request, 2 an agent's, the default, 3 sweeps and batches; a waiting job rises one level after a minute), then in turns between owners, then in order of submission. The job carries a frozen copy of the project: later edits do not change it (`project_changed` says they happened). `supersede` replaces a queued request of the same kind and targets; `hold` keeps it from starting until `queue_release`. A job that comes after another of the same project reuses the mesh the first built.", Expensive,
                schema(vec![("kind", json!({ "type": "string", "enum": ["solve", "mesh", "convergence"], "description": "default solve" })), ("physics", json!({ "oneOf": [{ "type": "string" }, { "type": "array", "items": { "type": "string" } }], "description": "physics code(s); default the active one" })), ("factor", number("convergence: comparison/present element size, default 0.6")), ("grids", json!({ "type": "integer", "enum": [2, 3], "description": "convergence: number of mesh levels (2 default, or 3 for Roache GCI)" })), ("coarsen", boolean("convergence: invert factor < 1 into a coarsening ladder relative to present mesh")), ("tolerance", number("convergence: largest accepted relative difference")), ("match", json!({ "type": "string", "enum": ["index", "nearest"], "description": "convergence: solution pairing strategy (index or nearest)" })), ("extrapolate", json!({ "type": "string", "enum": ["richardson_log_spd", "richardson_linear"], "description": "convergence: optional Richardson extrapolation mode" })), ("apply_extrapolated", boolean("convergence: promote extrapolated capacitance matrix onto stored result")), ("quantities", json!({ "type": "array", "items": { "type": "string", "description": "quantity filter" }, "minItems": 1, "description": "convergence: quantity filter list" })),
                    ("priority", int("1, 2 or 3; lower runs first")), ("supersede", boolean("replace a queued request of the same kind and targets")), ("hold", boolean("start held")),
                    ("shared", boolean("hand the solve to the queue of the workspace, which every program that uses the workspace shares (the window, the MCP server, `nyx queue`): a queue worker process runs the items one at a time in a worker, and keeps running when the program that asked is closed. Default: what the program was set up to do; implies a worker process")), ("detached", boolean("solve only, with `shared`: write the item into the queue and return at once, with nothing of this session waiting for it; the answer is the item (`shared_id`, and `results`, the file that `open_results` reads once its state is `done`). Closing the program leaves it running")),
                    ("isolate", boolean("solve in a worker process (only for kind solve): running out of memory or crashing fails this job and not the program, and a cancel stops it at once; the cost is a results round trip. Default: what the program was set up to do")), ("force", boolean("start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits")),
                    ("save_fields", boolean("solve only: when false, omit mesh and nodal/cell field arrays from worker HDF5 I/O and stored results after running any post_analyses; default true")), ("retain_amr_data", boolean("solve only: retain transient cellwise ElectroAmrData indicators on the stored result instead of stripping them after solve; default false")), ("post_analyses", post_analyses_schema()), ("memory_wait_timeout_s", number("solve only: maximum seconds (>= 0) to wait when free RAM is below the solve estimate (overrides NYX_MEMORY_WAIT_TIMEOUT_S and engine default)"))], &[])),
            cmd("queue_list", "List the queue", "Every job of this session and every item that another program wrote into the workspace queue (`foreign`, with a `shared_id`), oldest first, with its state (`queued`, `running`, `done`, `cancelled`, `failed`), lane, owner, priority, position in the queue, the time it waited and ran, its stage, what it was asked (`kind`, `targets`), `project_changed` and `applied`.", Read,
                schema(vec![("owner", text("only this owner: gui, ai or cli")), ("state", text("only this state"))], &[])),
            cmd("queue_get", "One job of the queue", "The entry of `queue_list` for one job and the last lines of its stage log. Give `id` (a job of this session) or `shared_id` (an item of the workspace queue).", Read, schema(vec![("id", int("job id")), ("shared_id", text("an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote"))], &[])),
            cmd("queue_cancel", "Cancel a queued or running job", "A queued job never runs. A running one stops at its next check (a factorisation in progress finishes first; a solve in a worker process is stopped at once); `cancel_requested` says so until it stops. Give `id` or `shared_id`.", Expensive, schema(vec![("id", int("job id")), ("shared_id", text("an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote"))], &[])),
            cmd("queue_hold", "Hold a queued job", "Keeps a queued job from starting until `queue_release`. Give `id` or `shared_id`.", Expensive, schema(vec![("id", int("job id")), ("shared_id", text("an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote"))], &[])),
            cmd("queue_release", "Release a held job", "Give `id` or `shared_id`.", Expensive, schema(vec![("id", int("job id")), ("shared_id", text("an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote"))], &[])),
            cmd("queue_set_priority", "Change the priority of a queued job", "1 runs first, 3 last. Give `id` or `shared_id`.", Expensive, schema(vec![("id", int("job id")), ("shared_id", text("an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote")), ("priority", int("1, 2 or 3"))], &["priority"])),
            cmd("queue_pause", "Pause the queue", "Nothing new starts; what runs finishes.", Expensive, schema(vec![], &[])),
            cmd("queue_resume", "Resume the queue", "", Expensive, schema(vec![], &[])),
            cmd("queue_wait", "Wait for a job, at most a while", "Waits for the job up to `timeout_s` seconds (required: a tool call must not block for ever) and applies its outcome. At the timeout it returns the job's entry with `waiting: true`, state and stage, and the job goes on.", Expensive,
                schema(vec![("id", int("job id")), ("timeout_s", number("longest wait in seconds"))], &["id", "timeout_s"])),
            cmd("sweep_run", "Run a sweep: one recorded run per point", "Queues one solve per explicit point (`points` or `series`) or streams a deterministic LHS / Sobol QMC ensemble (`sampling`, checking geometry and `SimplexMesh::check()` per sample and recording invalid samples as `skipped_invalid_geometry` without aborting the ensemble). An explicit point is a list of edits (`{ \"cmd\": \"set_param\", \"args\": { \"name\": \"w\", \"value\": 4e-7 } }`) applied to a copy of the project on demand without upfront `Vec<Engine>` cloning. `physics` is the code or list solved at every point (default: what a solve would solve). Points run at priority 3 on the engine's mesh/solve stage caches; optionally pass `checkpoint_path` (`.h5`) to stream points and skipped metadata incrementally to HDF5 `SweepCheckpoint`. Read the table with `sweep_table`.", Expensive,
                schema(vec![("name", text("a name for the table")), ("physics", json!({ "oneOf": [{ "type": "string" }, { "type": "array", "items": { "type": "string" } }], "description": "code or list of codes to solve at every point" })),
                    ("points", json!({ "type": "array", "items": { "type": "object", "properties": { "label": { "type": "string", "description": "a name for the point in the table; default `point N`" }, "edits": { "type": "array", "description": "the edits that make this point, applied in order to a copy of the project", "items": { "type": "object", "properties": { "cmd": { "type": "string", "description": "a command that changes the project" }, "args": { "type": "object", "description": "its arguments" } }, "required": ["cmd"] } } }, "required": ["edits"] }, "description": "the points, each a list of edits" })),
                    ("series", json!({ "type": "object", "properties": { "cmd": { "type": "string", "description": "the edit command repeated at every point" }, "args": { "type": "object", "description": "its arguments that stay the same" }, "field": { "type": "string", "description": "the argument that takes each value" }, "values": { "type": "array", "description": "the values, one point each" } }, "required": ["cmd", "field", "values"], "description": "one command repeated with `field` set to each value" })),
                    ("sampling", json!({ "type": "object", "description": "deterministic LHS / Sobol QMC streaming design-of-experiments (alternative to `points` or `series`)", "properties": {
                        "method": { "type": "string", "enum": ["lhs", "sobol"], "description": "QMC generator (`lhs` SplitMix64 Latin Hypercube Sampling, default; or low-discrepancy `sobol`)" },
                        "n_samples": { "type": "integer", "minimum": 1, "description": "number of QMC samples to stream (1..=4096)" },
                        "seed": { "type": "integer", "minimum": 0, "description": "deterministic 64-bit SplitMix64 seed" },
                        "scramble": { "type": "boolean", "description": "apply deterministic SplitMix64 XOR digital scramble to Sobol coordinates; default false" },
                        "offset": { "type": "integer", "minimum": 0, "description": "starting index in Sobol sequence (default 1 to avoid the origin corner)" },
                        "skip_invalid_geometry": { "type": "boolean", "description": "record geometry, meshing, or SimplexMesh::check failures as `skipped_invalid_geometry` without aborting the ensemble; default true" },
                        "parameters": { "type": "array", "description": "sampled parameter dimensions", "items": { "type": "object", "description": "one parametric margin (uniform `[min, max]` or normal `{ mean, sigma | rel_tol }`)", "properties": {
                            "name": { "type": "string", "description": "column label in checkpoint and table" },
                            "cmd": { "type": "string", "description": "edit command (default `set_param`)" },
                            "args": { "type": "object", "description": "base arguments kept constant (e.g. `{ \"name\": \"w\" }`)" },
                            "field": { "type": "string", "description": "argument key receiving each sampled value (default `value`)" },
                            "distribution": { "type": "string", "enum": ["uniform", "normal"], "description": "marginal distribution (default `uniform`, or `normal` when mean/sigma given)" },
                            "min": { "type": "number", "description": "lower bound (uniform range min, or optional normal clamp min)" },
                            "max": { "type": "number", "description": "upper bound (uniform range max, or optional normal clamp max)" },
                            "bounds": { "type": "array", "items": { "type": "number" }, "minItems": 2, "maxItems": 2, "description": "alternative `[min, max]` uniform range" },
                            "mean": { "type": "number", "description": "nominal mean for `normal` marginal" },
                            "sigma": { "type": "number", "description": "1-sigma standard deviation for `normal` marginal" },
                            "rel_tol": { "type": "number", "description": "relative 1-sigma tolerance `|mean| * rel_tol` for `normal` marginal" }
                        } } }
                    }, "required": ["n_samples", "parameters"] })),
                    ("checkpoint_path", text("optional .h5 file path for incremental HDF5 SweepCheckpoint streaming of solved points and skipped invalid geometries")),
                    ("checkpoint_lzf", boolean("compress checkpoint columns with Shuffle+LZF; default false")),
                    ("priority", int("1, 2 or 3; default 3")), ("hold", boolean("queue the points held; release each with queue_release")), ("reuse", json!({ "type": "string", "enum": ["none", "warm"], "description": "`warm`: each point starts its eigensolves from the modes of the previous point instead of a random block (the first point only records them), which shortens the iterations of a sweep that keeps its mesh (a material, a study value, a boundary condition, a frequency with the mesh not sized by the physics; a point whose mesh or number of unknowns differs starts cold). The numbers agree to the solver's tolerance but not bit for bit with a cold solve, so a warm point is marked `exploratory` in its row and run (`reuse`) until it is solved cold. Default `none`: every point is a cold solve" })), ("shared", boolean("hand the solve to the queue of the workspace, which every program that uses the workspace shares (the window, the MCP server, `nyx queue`): a queue worker process runs the items one at a time in a worker, and keeps running when the program that asked is closed. Default: what the program was set up to do; implies a worker process")), ("isolate", boolean("solve every point in a worker process, as for solve_run")), ("force", boolean("start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits")), ("save_fields", boolean("when false, omit mesh and nodal/cell field arrays from worker HDF5 I/O and stored results at every sweep point after running any post_analyses; default true")), ("post_analyses", post_analyses_schema()), ("memory_wait_timeout_s", number("maximum seconds (>= 0) to wait per sweep point when free RAM is below the solve estimate (overrides NYX_MEMORY_WAIT_TIMEOUT_S and engine default)"))], &[])),
            cmd("sweep_table", "The table of a sweep", "Applies the points that have finished and returns one row per point: its label and edits, its state (`queued`, `running`, `done`, `skipped_invalid_geometry`, `cancelled`, `failed`), the numbers of each physics solved (eigenvalues, values, numeric scalars), the unknowns and cells, the seconds and the run it was recorded as, plus `skipped_invalid_geometry` counts and optional HDF5 `checkpoint_path` when running QMC/checkpointed sweeps. With `wait_s` it waits up to that many seconds for the points to finish (`complete` says whether they all have). Without `sweep` it lists the sweeps of this session.", Read,
                schema(vec![("sweep", int("the sweep's number")), ("wait_s", number("longest wait in seconds, default 0"))], &[])),
            cmd("sweep_cancel", "Cancel the points of a sweep that have not finished", "Queued points never run; a running one stops at its next check.", Expensive, schema(vec![("sweep", int("the sweep's number"))], &["sweep"])),
            cmd("poll_jobs", "Apply finished jobs", "Applies every job that has finished since the last call and returns what each came to.", Read, schema(vec![], &[])),
            cmd("result_info", "Describe the results", "For each physics solved: the study, size, time, solutions with their values, fields (node-averaged, real) and `cell_fields` (the per-cell fields the reductions read: name, unit, complex or not), scalars, unconditional solver `diagnostics` and numerical `error_budget`, `stale` (whether what this result was made from has changed since: its own problem, the materials, the mesh or a physics it takes fields from) and `age_s` (seconds since it was stored), and `validated`: `{state, tier, reason, check?}` where `state` is `exploratory` (no mesh-convergence check with a tolerance has passed on this problem), `checked` (one did or flux bounds certified it) or `check_not_affordable` (not checked, and a check is expected not to fit this machine). A stale result is never checked; every answer that carries `stale` carries `validated`.", Read, schema(vec![("physics", phys())], &[])),
            cmd("result_field", "Read a field", "The values of one field of one solution (sampled evenly when there are more than `max_values`), with `stale` and `age_s` of the result it was read from.", Read,
                schema(vec![("physics", phys()), ("solution", int("solution number from 0, default 0")), ("field", text("field name")), ("max_values", int("default 2000"))], &["physics", "field"])),
            cmd("get_guide", "Read a nyx guide", "Read a short guide. Call with topic nyx-instructions before using other nyx tools. Other topics: geometry, optical-modes (reading mode tables), cqed (superconducting qubit capacitance, inductance, EPR, Hamiltonian and TLS participation), troubleshooting (errors).", Read,
                schema(vec![("topic", json!({ "type": "string", "enum": crate::guide::TOPICS, "description": "which card; omit for nyx-instructions" }))], &[])),
            cmd("field", "A handle to a field", "A handle to a field of a solution, to compute with (`integrate`, `extremum`, `probe`, `line_cut`, `overlap`, `confinement`, `field_op`). A handle is a small value, not the field: it names the physics, the solution and the field, and is read from the stored result each time it is used. A name is a per-cell field of the solution (`Ex`, `uz`, ...: complex when the solve is) or else a flat field of `result_info` (`T`, `|E|`: real). A vector name (`E` when `Ex` and `Ez` exist) answers one handle per component.", Read,
                schema(vec![("physics", phys()), ("solution", int("solution number from 0, default 0")), ("name", text("field name, e.g. Ex, or a vector name such as E"))], &["name"])),
            cmd("field_op", "Arithmetic on handles", "Combines handles: `add`, `mul` (two or more operands), `sub`, `div`, `neg`, `abs`, `abs2` (|z|^2 with no square root), `re`, `im`, `conj`, `sqrt` and `pow` (with `n`, a whole number 1 to 16). An operand is a handle, a field name of `physics`/`solution`, a number, or another operation `{op, args}`. The operation is checked against the results and returned as one handle.", Read,
                schema(vec![("op", text("add, sub, mul, div, neg, abs, abs2, re, im, conj, sqrt or pow")), ("args", json!({ "type": "array", "description": "the operands" })), ("n", int("the power, for pow")), ("wavelength", wl()), ("physics", phys()), ("solution", int("default 0"))], &["op", "args"])),
            cmd("integrate", "Integrate an expression over a region", "The integral of an expression of fields over a region (a region name, comma list or `name*`; default the whole mesh), with the degree-4 rule of the finite elements (exact for the product of two P2 fields). A real expression answers a number, a complex one [re, im]. `measure` is the area or volume of the region. `revolution: true` weights by 2 pi r for a body of revolution (a 2D r-z mesh). BATCH REDUCTIONS: pass `solutions` (`\"all\"` or a list of 0-based indices), `by_region: true`, or a batch `expr` list `[\"abs2(Ex)\", \"abs2(Ey)\"]` / named map `{\"energy\": \"...\", \"ex2\": \"...\"}` to evaluate every `(solution, region, expr)` combination in a single bit-exact quadrature pass. Answers `stale` and `age_s` of what it read. FORMULAS (every `expr` of the reductions takes one): a text such as `abs(Ex)^2 * eps_r` or `where(x > 0, T, 0)`, parsed once and evaluated at each point in Rust. Numbers (`2j` is imaginary), `pi`, `e`, `+ - * / ^`, comparisons, `and or not`; functions `abs abs2 re im conj sqrt exp sin cos arg log log10 tan sinh cosh tanh atan step min max atan2` and `where(cond, a, b)` (only the branch that applies is evaluated, so `where(x != 0, 1/x, 0)` is safe). Names: the coordinates `x y z`; the fields of the result (`Ex`, `T`; in backticks when the name is not an identifier, `` `|E|` ``; `` `ems:Ex` `` or `` `ems:1:Ex` `` in backticks for another physics or solution); material values by their descriptor name (`eps_r`, `rho`; `eps` means `eps_r`; `material(\"rho\")` to be explicit); `region(\"Si\")` is 1 inside the region and 0 outside. A name like `E` stands for the vector (Ex, Ey, Ez) when the solution has those components: `dot(a, b)`, `inner(a, b)` (conjugates a), `cross`, `norm(v)`, `norm2(v)`, `E.x`. Orders, comparisons, `log` and `tan` read real values: wrap a complex one in `abs`, `re` or `im`. A material value is the one the solver saw in each cell. A tensor property (`eps_r`, `k`, `alpha`) has its entries as `eps_r_xx`, `eps_r_yy`, `eps_r_zz`, `eps_r_xy`, `eps_r_xz`, `eps_r_yz` in wafer axes with the crystal orientation applied (an isotropic material is the scalar times the identity); the plain `eps_r` is refused on a mesh where some material is a tensor. BOUNDARY INTEGRALS: with `boundary` (an outer boundary such as \"window.top\", an interface \"oxide|silicon\", or an interior 3D sheet) the sum is over the facets of the boundary instead of the cells: a line integral on a 2D mesh, a surface integral on a 3D one, with `revolution` weighting by 2 pi r. The formula may use `nx`, `ny`, `nz`, the unit normal in wafer axes, so the flux of a vector F is `dot(F, [nx, ny, nz])`. At an interface or interior sheet give `side` (a region name, or \"+z\"/\"-z\" for a sheet inside one region); the normal points out of it. On an outer boundary no `side` is needed. `wavelength` (m) reads the material values at that wavelength, not the study's; the fields stay those of the solve.", Read,
                schema(vec![("expr", expr_schema("a formula string such as \"abs(Ex)^2 * eps_r\", a handle/AST object ({op, args}, {field, ...}), a batch list of expressions [\"abs2(Ex)\", \"abs2(Ey)\"], or a batch named map {\"energy\": \"...\", \"ex2\": \"...\"}")), ("region", text("region name(s); default all")), ("by_region", boolean("accumulate per-region integrals across every domain (or within `region`) in a single pass, default false")), ("boundary", text("integrate over the facets of this boundary (an outer boundary name, an interface \"a|b\", or a 3D sheet) instead of over cells: a line integral on a 2D mesh, a surface integral on 3D. The formula may use nx, ny, nz, the unit normal. Not with region or by_region")), ("side", text("with boundary: the region (or \"+z\"/\"-z\" for a sheet inside one region) an interface is read in; the normal points out of it. Not needed on an outer boundary")), ("revolution", boolean("weight by 2 pi r, default false")), ("wavelength", wl()), ("physics", phys()), ("solution", int("0-based solution index, default 0 (mutually exclusive with solutions)")), ("solutions", json!({ "oneOf": [{ "type": "string", "enum": ["all"] }, { "type": "array", "items": { "type": "integer", "minimum": 0 } }], "description": "\"all\" or a list of 0-based solution indices to integrate in a single pass" }))], &["expr"])),
            cmd("extremum", "Largest or smallest value of an expression", "The largest (`op: max`, default) or smallest (`min`) real value of an expression in a region. By default reports the extremum at the nodes while checking P2 candidate cells for larger interior peaks (`interior_extremum`); pass `interior: true` to return the P2 sub-sampled/Newton-refined interior extremum directly as `value`. A complex expression has no order: wrap it in `abs`, `abs2`, `re` or `im`.", Read,
                schema(vec![("expr", expr_schema("a formula text such as \"abs(Ex)^2 * eps_r\" (see integrate), a handle, a field name, a number or {op, args}")), ("op", json!({ "type": "string", "enum": ["max", "min"], "description": "default max" })), ("interior", boolean("for P2 fields, return the sub-sampled and Newton-refined interior extremum directly in `value` instead of the nodal value; default false still checks peak candidates and surfaces `interior_extremum` when higher")), ("region", text("region name(s); default all")), ("wavelength", wl()), ("physics", phys()), ("solution", int("default 0"))], &["expr"])),
            cmd("probe", "Value of an expression at points", "The value of an expression at points, interpolated in the cell that holds each point. Points are in mesh axes: [x, z] on a 2D mesh, [x, y, z] on a 3D one. A point outside the mesh answers `inside: false` and a null value.", Read,
                schema(vec![("expr", expr_schema("a formula text such as \"abs(Ex)^2 * eps_r\" (see integrate), a handle, a field name, a number or {op, args}")), ("points", json!({ "type": "array", "description": "a list of points, each a list of coordinates", "items": { "type": "array", "items": { "type": "number" } } })), ("wavelength", wl()), ("physics", phys()), ("solution", int("default 0"))], &["expr", "points"])),
            cmd("line_cut", "Values of an expression along a line", "The value of an expression at `n` equally spaced points from `from` to `to` (inclusive), with the distance `s` along the line. Points outside the mesh answer `inside: false`.", Read,
                schema(vec![("expr", expr_schema("a formula text such as \"abs(Ex)^2 * eps_r\" (see integrate), a handle, a field name, a number or {op, args}")), ("from", json!({ "type": "array", "items": { "type": "number" }, "description": "start point: [x, z] on a 2D mesh, [x, y, z] on a 3D one" })), ("to", json!({ "type": "array", "items": { "type": "number" }, "description": "end point" })), ("n", int("number of points, default 101")), ("wavelength", wl()), ("physics", phys()), ("solution", int("default 0"))], &["expr", "from", "to"])),
            cmd("overlap", "Overlap of two fields", "The overlap integral of two fields (or two lists of components, e.g. [Ex, Ey, Ez]) over a region: the complex integral of sum_k conj(a_k) b_k, the two norms (integrals of sum_k |c_k|^2) and the coefficient |overlap|^2 / (norm_a norm_b), a number from 0 to 1.", Read,
                schema(vec![("a", expr_schema("a field, a formula, an expression, or a list of components")), ("b", expr_schema("a field, a formula, an expression, or a list of components")), ("region", text("region name(s); default all")), ("revolution", boolean("weight by 2 pi r, default false")), ("wavelength", wl()), ("physics", phys()), ("solution", int("default 0"))], &["a", "b"])),
            cmd("confinement", "Fraction of a field in a region", "The integral of sum_k |c_k|^2 over one or more regions divided by the same over the whole mesh, for the components `fields` (e.g. [Ex, Ey, Ez]). Without `weight` it has no material weight; give `weight` a real formula such as \"eps_r\" for the energy fraction (the integrand is weight * sum_k |c_k|^2, in each region and over the whole mesh). Pass `solutions` (`\"all\"` or a list of 0-based indices), `by_region: true`, or `regions` (`[str]`) to accumulate all requested solutions and regions in a single bit-exact quadrature pass.", Read,
                schema(vec![("fields", expr_schema("a field, a formula, an expression, or a list of components")), ("expr", expr_schema("optional alias for `fields`: a field, formula, expression, or list of components")), ("weight", text("a real formula that multiplies the integrand at each point, e.g. eps_r; default none")), ("region", json!({ "oneOf": [{ "type": "string" }, { "type": "array", "items": { "type": "string" } }], "description": "the region (or list of regions) the fraction is of; omit when `regions` or `by_region: true` is given" })), ("regions", json!({ "type": "array", "items": { "type": "string" }, "description": "batch list of region names to compute confinement fractions over in a single pass" })), ("by_region", boolean("compute confinement fractions across every region of the mesh (or within `region`) in a single pass, default false")), ("revolution", boolean("weight by 2 pi r, default false")), ("wavelength", wl()), ("physics", phys()), ("solution", int("0-based solution index, default 0 (mutually exclusive with solutions)")), ("solutions", json!({ "oneOf": [{ "type": "string", "enum": ["all"] }, { "type": "array", "items": { "type": "integer", "minimum": 0 } }], "description": "\"all\" or a list of 0-based solution indices to evaluate in a single pass" }))], &["fields"])),
            cmd("add_absorbing_layer", "Add an outward absorbing layer", "Adds an outward zero-padding absorbing layer (`Op::Pml` with `mode = \"Outward\"`, `window = null`, `coord = \"Auto\"`, `layer_grading = 2.0`, `thickness = \"auto\"`, `layers = \"auto\"`) around the exterior of the physical domain across Cartesian, Cylindrical, or Spherical topologies.", Edit,
                schema(vec![
                    ("name", text("feature name; omit for pml1, pml2, ...")),
                    ("sides", json!({ "type": "array", "items": { "enum": ["XMin", "XMax", "YMin", "YMax", "ZMin", "ZMax"] }, "description": "active sides (empty: all)" })),
                    ("thickness", json!({ "oneOf": [{ "type": "number" }, { "type": "string" }], "description": "layer thickness in metres, or \"auto\" (default \"auto\")" })),
                    ("layers", json!({ "oneOf": [{ "type": "number" }, { "type": "string" }], "description": "number of structured element layers across the shell, or \"auto\" (default \"auto\")" })),
                    ("grading", json!({ "oneOf": [{ "type": "number" }, { "type": "string" }], "description": "polynomial stretch profile order m, default 2.0" })),
                    ("layer_grading", json!({ "oneOf": [{ "type": "number" }, { "type": "string" }, { "type": "null" }], "description": "geometric mesh layer grading ratio (>= 1; default 2.0)" })),
                    ("coord", json!({ "type": "string", "enum": ["Auto", "Cartesian", "Cylindrical", "Spherical", "auto", "cartesian", "cylindrical", "spherical"], "description": "coordinate topology, default \"Auto\"" })),
                    ("center", json!({ "type": ["array", "null"], "items": { "oneOf": [{ "type": "number" }, { "type": "string" }] }, "minItems": 2, "maxItems": 3, "description": "optional centre [x, y, z] (or [x, z] in 2D) for Cylindrical or Spherical shells" })),
                    ("axis", json!({ "type": ["string", "null"], "enum": ["X", "Y", "Z", null], "description": "optional cylinder axis (X, Y, or Z)" })),
                    ("regions", json!({ "type": "array", "items": { "type": "string" }, "description": "regions to extrude (empty: all)" })),
                    ("mode", json!({ "type": "string", "enum": ["Inward", "Outward", "inward", "outward"], "description": "growth mode, default \"Outward\"" })),
                    ("window", json!({ "oneOf": [{ "type": "string" }, { "type": "integer" }, { "type": "null" }], "description": "optional window feature when mode is Inward" })),
                ], &[])),
            cmd("add_pml_conditions", "Add the absorbing-layer conditions", "Adds a `pml` (or `infinite_element` for static studies) condition for every absorbing piece the geometry made with `add_pml` or `add_absorbing_layer` (replacing earlier ones on those pieces), so the physics absorbs where the geometry says. When `strength <= 0` or `auto_scale: true`, parameters are jointly auto-calibrated across materials sharing each PML slab/shell from the study's wavelengths.", Edit,
                schema(vec![
                    ("physics", phys()),
                    ("strength", number("imaginary stretch at the outer edge (default 3; <= 0 auto-calibrates from study wavelengths)")),
                    ("profile", json!({ "type": ["string", "null"], "enum": ["auto", "kappa_polynomial", "rational", "exponential", null], "description": "stretching profile (\"auto\" selects kappa_polynomial for eigenproblems/waveguides and rational for direct driven/static problems)" })),
                    ("auto_scale", boolean("jointly auto-calibrate strength, profile, and kappa_max/curvature from study wavelengths across materials sharing each PML slab/shell")),
                    ("target_reflection", number("target normal-incidence amplitude reflection R_0 when auto_scale is enabled, default 1e-6")),
                ], &["physics"])),
            cmd("pml_diagnostics", "Absorbing-layer mesh and stability diagnostics", "Returns absorbing-layer diagnostics across all PML domains: uncut `element_counts`, `min_pml_scaled_jacobian`, `pml_warnings`, normalized penetration depth `s_dist_bounds` `[0.0, 1.0]`, and pre-solve Christoffel backward-wave check (`has_backward_wave`, `min_projection`, `recommended_mpml_ratio`).", Read, schema(vec![], &[])),
            cmd("list_materials", "List the materials", "The materials a problem can name: those of the library (material files) and those the project carries.", Read, schema(vec![], &[])),
            cmd("describe_material", "Describe a material", "What a material gives (elastic class, optical model, thermal, electrical) and its TOML text.", Read, schema(vec![("name", text("material name"))], &["name"])),
            cmd("add_material", "Add a material to the project", "Adds (or replaces) a material the project carries, as TOML text (units are written in the file; `extends = \"silicon\"` starts from a library material). Read at once, so a mistake is reported here.", Edit,
                schema(vec![("name", text("material name")), ("toml", text("the material file text"))], &["name", "toml"])),
            cmd("remove_material", "Remove a project material", "", Edit, schema(vec![("name", text("material name"))], &["name"])),
            cmd("list_remote_materials", "List downloadable remote materials", "Curated presets from open material databases (`refractiveindex.info`, `AFLOW`, `NumBAT`) that `fetch_material` can download and convert on demand.", Read,
                schema(vec![("domain", json!({ "type": "string", "enum": ["all", "photonics", "phononics", "superconducting"], "description": "physics domain filter, default all" }))], &[])),
            cmd("import_material", "Import an external material record", "Converts raw payload text from an open database (`rii` YAML, `aflow` JSON, `mp` JSON or `numbat` JSONC) into validated TOML, optionally extending a library material and adding it to the project or writing `<name>.toml` into `folder`.", Edit,
                schema(vec![
                    ("name", text("material name without slashes or commas")),
                    ("format", json!({ "type": "string", "enum": ["rii", "aflow", "mp", "numbat"], "description": "source database format: rii, aflow, mp or numbat" })),
                    ("content", text("raw YAML, JSON or JSONC payload text from the external database")),
                    ("description", text("optional one-line comparative summary written into the material description field")),
                    ("extends", text("optional base material name in the library to inherit unprovided properties from")),
                    ("add_to_project", boolean("add the converted material to the active project (default true)")),
                    ("folder", text("optional destination folder to write <name>.toml into")),
                ], &["name", "format", "content"])),
            cmd("fetch_material", "Download a material from an open database", "Downloads a preset or database entry over HTTPS (no redirects) from an allowlisted host (`api.materialsproject.org`, `aflowlib.org`, `aflow.org`, `aflowlib.duke.edu`, `refractiveindex.info`, or curated `raw.githubusercontent.com` repo paths), converts it into validated TOML, and optionally adds it to the project or writes `<name>.toml` into `folder`. Refused in sandboxed script mode (`path_root`).", Expensive,
                schema(vec![
                    ("source", text("curated preset id (from `list_remote_materials`), `catalog`, `rii`, `aflow`, `mp` or `numbat`")),
                    ("target", text("preset id, database entry key, Materials Project id (e.g. mp-149), AFLOW query, NumBAT stem or allowlisted https URL (refused for `mp` or when an API key is present); omit when `source` is a preset id")),
                    ("name", text("material name in TOML (`[A-Za-z0-9_-]+`, 1..=64, no Windows reserved device names); defaults to the preset id or derived entry stem")),
                    ("description", text("optional one-line comparative summary written into the material description field")),
                    ("extends", text("optional base material name in the library to inherit unprovided properties from")),
                    ("add_to_project", boolean("add the downloaded material to the active project (default true)")),
                    ("folder", text("optional destination folder to write <name>.toml into")),
                    ("api_key", text("optional Materials Project API key (falls back to MP_API_KEY environment variable; sent only to api.materialsproject.org via curl stdin, never stored)")),
                ], &["source"])),
            cmd("search_remote_materials", "Search curated and live open material databases", "Searches the curated remote material presets by keyword, formula, domain or source database, and optionally queries live AFLOW and Materials Project HTTPS endpoints (live queries are refused in sandboxed script mode).", Read,
                schema(vec![
                    ("query", text("search keyword, chemical formula or material identifier (e.g. Si, AlN, LiNbO3, mp-149)")),
                    ("domain", json!({ "type": "string", "enum": ["all", "photonics", "phononics", "superconducting"], "description": "physics domain filter, default all" })),
                    ("source", json!({ "type": "string", "enum": ["all", "rii", "aflow", "mp", "numbat"], "description": "database filter, default all" })),
                    ("live", boolean("also query live AFLOW and Materials Project HTTPS endpoints (refused in sandboxed script mode), default false")),
                    ("api_key", text("optional Materials Project API key for live mp queries (falls back to MP_API_KEY environment variable)")),
                ], &["query"])),
            cmd("evaluate_material", "Evaluate a material at temperature, frequency and wavelength", "Evaluates a material's temperature-, frequency-, wavelength- and electric-field-dependent properties: thermal curves at T, optical (n, k) and complex permittivity at (wavelength, T), microwave complex permittivity and TLS loss tangent at (frequency, T, e_field), piezoelectric Legendre duals (eps_clamped, eps_free, C^E, C^D), photoelastic, electro-optic, viscoelastic, and BCS superconducting surface impedance.", Read,
                schema(vec![
                    ("name", text("material name in the library or active project")),
                    ("temperature", number("evaluation temperature in K (default: material reference_temperature)")),
                    ("frequency", number("evaluation frequency in Hz for microwave dielectric and BCS superconducting models")),
                    ("wavelength", number("evaluation free-space wavelength in m for optical dispersion models")),
                    ("e_field", number("electric-field amplitude in V/m for saturable TLS dielectric loss")),
                    ("thickness", number("film thickness in m for BCS superconducting sheet impedance")),
                    ("film_model", json!({ "type": "string", "enum": ["kerr", "two_sided", "thin_limit"], "description": "thin-film boundary formula for BCS sheet impedance, default kerr" })),
                    ("euler_deg", json!({ "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3, "description": "Bunge Z-X-Z crystal orientation angles in degrees [phi1, Phi, phi2]" })),
                ], &["name"])),
            cmd("bcs_surface_impedance", "BCS thin-film surface impedance", "Evaluates Zimmermann arbitrary-purity BCS complex conductivity and thin-film surface impedance (Z_s = R_s + i omega L_k_sq) for a superconducting material or parameter set.", Read,
                schema(vec![
                    ("material", text("name of a library or project material with a [superconducting] block")),
                    ("freq", json!({ "type": "number", "minimum": 1.0, "description": "frequency in Hz" })),
                    ("thickness", json!({ "type": "number", "minimum": 1e-12, "description": "film thickness in metres" })),
                    ("temperature", number("temperature in K, default 0.01")),
                    ("film_model", json!({ "type": "string", "enum": ["two_sided", "kerr", "thin_limit"], "default": "kerr", "description": "thin-film boundary convention: two_sided (interior sheet jump), kerr (one-sided backed wall) or thin_limit (Pearl d << lambda_eff)" })),
                    ("tc", number("critical temperature T_c in K (when material is omitted or overridden)")),
                    ("gap_0", number("zero-temperature energy gap Delta(0) in J")),
                    ("gap_0_eV", number("zero-temperature energy gap Delta(0) in eV")),
                    ("lambda_L", number("clean zero-temperature London penetration depth in metres")),
                    ("xi_0", number("BCS coherence length in metres")),
                    ("mean_free_path", number("electron mean free path in metres")),
                    ("sigma_n", number("normal-state conductivity in S/m")),
                    ("dynes_gamma", number("Dynes pair-breaking broadening Gamma in J")),
                    ("r_res", number("residual sheet resistance floor in Ohm/sq")),
                    ("x_qp_floor", number("non-equilibrium normalized quasiparticle density floor")),
                ], &["freq", "thickness"])),
            cmd("convergence_check", "Check mesh convergence", "Solves the physics on 1 comparison mesh (2-grid default) or 2 comparison meshes (3-grid Roache GCI when `grids: 3`) with element sizes scaled by `factor` (`< 1` refines, `> 1` coarsens, or inverted coarsening ladder `[1.0, 1/factor, (1/factor)^2]` when `coarsen: true`; the active project is never mutated) and compares every eigenvalue, value and scalar across grids. By default `scale_sizes` is true so non-suppressed `add_size` rules (`h`, `max_deviation`, `min_layers`) and `gap_cells` scale proportionally on comparison copies. With `tolerance` each row says pass or fail, `converged` summarises, and the verdict is kept (`validated` becomes `checked` with tier `gci_verified` / `two_grid_checked` or stays `exploratory` / `non_convergent`). Sub-solves are isolated through `nyx-worker` and admitted against machine RAM when available. When mode ordering changes between meshes, `match: \"nearest\"` pairs mutual nearest eigenvalues and `track` filters candidate modes by target (`f`, `n_b` or `n_eff`) and/or physical mode signatures. With `extrapolate` (`\"richardson_log_spd\"` or `\"richardson_linear\"`), computes Richardson extrapolation (`extrapolate_order` defaults to 1.0 on 3D/sheet capacitance matrices, else `2 * order`) and promotes extrapolated capacitance matrices when `apply_extrapolated` is true (tagging `SpatialErrorEstimate.extrapolated_scalars_only = true`).", Expensive,
                schema(vec![("physics", phys()), ("factor", number("comparison/present element-size ratio, positive and not 1 (< 1 refines, > 1 coarsens), default 0.6")), ("grids", json!({ "type": "integer", "enum": [2, 3], "description": "number of mesh levels to solve: 2 (default fast check) or 3 (Roache GCI with observed order and safety factor)" })), ("coarsen", boolean("when true and factor < 1, invert factor into a coarsening ladder [1.0, 1/factor, (1/factor)^2] relative to the present mesh when memory headroom is tight, default false")), ("coarsen_factors", json!({ "type": "array", "items": { "type": "number" }, "description": "explicit coarsening element-size multipliers (> 1.0) relative to the present mesh (1 entry for grids=2, 2 entries for grids=3)" })), ("scale_sizes", boolean("scale non-suppressed add_size features (h, max_deviation, min_layers) proportionally with factor on comparison copies, default true")), ("tolerance", number("largest accepted relative difference (or 3-grid GCI), e.g. 1e-3; omit to only report")), ("force", boolean("start convergence sub-solves even when machine memory estimate flags low headroom, default false")), ("isolate", boolean("override worker process isolation for convergence sub-solves (default uses nyx-worker when configured)")),
                    ("match", json!({ "type": "string", "enum": ["index", "nearest"], "description": "how the solutions of the meshes are paired: index (default) compares solution i with solution i; nearest pairs each mode with the nearest eigenvalue of the other mesh in the complex plane and lists the modes without a partner" })),
                    ("track", json!({
                        "type": "object",
                        "description": "follow mode(s) matching a target eigenvalue (f, n_b or n_eff) and/or physical mode-signature filters across meshes",
                        "properties": {
                            "f": { "type": "number", "description": "target eigenvalue real part (e.g. frequency in Hz)" },
                            "n_b": { "type": "number", "description": "target Bloch index real part" },
                            "n_eff": { "type": "number", "description": "target effective index real part" },
                            "region": { "type": "string", "description": "region name whose energy share (`share in <region> N`) must be >= min_share" },
                            "min_share": { "type": "number", "description": "minimum energy share in `region`, in [0, 1]" },
                            "component": { "type": "string", "description": "field component name (e.g. Ex, Ey, Ez) whose share must be >= min_component" },
                            "min_component": { "type": "number", "description": "minimum component share in [0, 1]" },
                            "max_layer_fraction": { "type": "number", "description": "maximum absorbing-layer energy fraction (`layer fraction N`), in [0, 1]" }
                        },
                        "additionalProperties": false
                    })),
                    ("extrapolate", json!({ "type": "string", "enum": ["richardson_log_spd", "richardson_linear"], "description": "optional two-grid Richardson extrapolation mode: richardson_log_spd (SPD-preserving log-space Maxwell capacitance extrapolation) or richardson_linear" })),
                    ("extrapolate_order", number("asymptotic convergence order p > 0 for Richardson extrapolation and GCI theoretical order; default 1.0 for 3D/sheet capacitance matrices, else 2 * order")),
                    ("apply_extrapolated", boolean("when true, promote the extrapolated capacitance matrix onto the stored latest result of physics and mark extrapolated_scalars_only in error_budget, default false")),
                    ("quantities", json!({
                        "type": "array",
                        "items": {
                            "type": "string",
                            "description": "quantity name or scalar name to include in the verdict (e.g. \"eigenvalue[0].re\", \"scalar C_0_0 (F)\", or \"C_0_0 (F)\")"
                        },
                        "minItems": 1,
                        "description": "non-empty list of quantity names to check; only matching rows count toward converged, worst_relative, and worst_quantity"
                    })),
                    ("max_cells", number("refuse, before any solve, a check whose comparison mesh has more cells than this (an adaptive 2D mesh has 6 to 8 times as many cells at factor 0.6)")),
                    ("max_seconds", number("refuse, before any solve, a check estimated (from the solves this machine has run; no refusal without an estimate) to take longer than this for the comparison mesh")),
                    ("dry_run", boolean("only build the comparison mesh(es) and answer their cells and what a solve of each is expected to cost; nothing is solved")),
                    ("wait", boolean("wait for the answer, default true; false starts a job that job_list, job_status and job_cancel see and job_wait answers"))], &["physics"])),
            cmd("export_results", "Write the results to an HDF5 file", "Writes the project, the mesh and the latest result of every solved physics (all fields, compressed, exact) into an HDF5 file and an XDMF index of the stored datasets next to it, for nyx and its tools only (not for ParaView or VisIt). Needs a build with HDF5 support; `parallel` also needs `phdf5`.", Expensive, schema(vec![("path", text("file path, .h5")), ("parallel", boolean("write with the MPI-IO driver; every process of an MPI run calls this with the same arguments and writes its own rows (needs a build with the `phdf5` feature)")), ("wait", boolean("wait for the file, default true; false starts a job that writes what the engine holds now, under its name only when complete"))], &["path"])),
            cmd("results_file_info", "Describe an HDF5 results file", "The project title, the mesh size and for each physics the study, size, time, number of solutions, scalars, solver diagnostics and numerical error budget of a results file.", Read, schema(vec![("path", text("file path"))], &["path"])),
            cmd("open_results", "Open an HDF5 results file", "Makes the project of a results file the project, builds its mesh again (meshing is deterministic) and checks it is the stored mesh, then takes the stored results as the latest results, so every analysis, plot and export works on them as after a solve. Replaces the present project. Needs a build with HDF5 support. With `wait` false it starts a job (reading and meshing happen on it; the project is replaced when job_wait or poll_jobs applies it).", Edit, schema(vec![("path", text("file path, .h5")), ("wait", boolean("wait for it, default true"))], &["path"])),
            cmd("list_analyses", "List the analyses", "The analyses of results (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator, uq_cqed), which solved physics each applies to and, per physics, whether its result is `stale`.", Read, schema(vec![], &[])),
            cmd("analyse", "Analyse a result", "Runs an analysis on the latest result of a physics: polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator, or uq_cqed. For `capacitance_matrix`, optional `params` perform ground projection or floating Schur elimination (`ground`, `eliminate: \"drop\" | \"floating_schur\"`), grounded/differential lumped qubit reduction (`qubits: [{ name, pads, E_J_Hz? }]`), resonator ZPF / dispersive coupling extraction (`resonators: [{ name, terminal, f_r_Hz, Z0_ohm?, resonator_type?, C_r_F? }]`), and matched-asymptotic finite-metal-thickness correction (`metal_thickness_m`, `edge_cutoff_m`). For `hamiltonian`, outputs include `f_cav_Hz`, `g_matrix_Hz`, `kappa_cav_Hz`, and `t1_purcell_s`, and `params` accept `sweep_port` (0-based index or port name), `squid_asymmetry`, `squid_flux_sweep`, and `ej_hz_sweep`. For `tls_participation`, interface thicknesses/permittivities/loss tangents auto-load from the material TOML `[interfaces.{ma,ms,sa}]`, `params` accept `trench_depth_m` (substrate over-etch depth in metres), and passing `electrostatic_physics` (e.g. `\"dfs\"`, with optional `port_terminals`, `modal_voltages`, `modal_freqs_hz`) evaluates hybrid quasi-electrostatic modal TLS participation and `T1_TLS_s` on the fine electrostatic mesh. For `network`, extracts damped companion-QEP modes, ROM-synthesized Z/Y/S matrices, Schur-complement Black-Box Quantization (`bbq_modes`), Taubin algebraic circle-fit resonator parameters (`s21_fit`), and Touchstone v1.1 (`.sNp`) text from either `CircuitRom` or `DrivenPortData`. For `piezo_resonator`, extracts series/parallel resonance frequencies (`f_s_Hz`, `f_p_Hz`), IEEE coupling factors (`k_t2`, `k_eff2`), unwrapped phase-slope quality factors (`Q_s`, `Q_p`), and bounded 6-parameter mBVD equivalent circuit (`C_0`, `R_0`, `R_s`, `C_m`, `L_m`, `R_m`) from `DrivenPortData`. For `uq_cqed`, propagates fabrication and constitutive 1-sigma `uncertainties` (plus `seed`, `n_samples`, `electrostatic_physics` companion bounds) through zero-solve EPR/Hamiltonian/Capacitance/TLS post-processing while keeping ASME V&V 20 epistemic brackets strictly separated from aleatoric variances. The reply carries `stale`, `age_s`, and `validated` of the result it was run on.", Read,
                schema(vec![("analysis", text("analysis id (list_analyses)")), ("physics", phys()), ("params", object("optional analysis parameters (capacitance_matrix: ground, eliminate, qubits, resonators, metal_thickness_m, edge_cutoff_m; hamiltonian: basis, fock_trunc, sweep_port, flux_sweep, squid_asymmetry, squid_flux_sweep, ej_hz_sweep; tls_participation: trench_depth_m, t_ma/t_ms/t_sa, electrostatic_physics, port_terminals, modal_voltages, modal_freqs_hz; network: freqs_hz, terminations_ohm, active_ports, s21_ports, resonator_branch, touchstone_format; piezo_resonator: port, terminations, f_window_hz, z_ref_ohm, fix_r_s_ohm; uq_cqed: seed, n_samples, basis, uncertainties, electrostatic_physics, qubits, resonators, metal_thickness_m; epr options)"))], &["analysis", "physics"])),
            cmd("figure", "Draw the modes or 3D fields of a result", "Writes a 2D or 3D figure of the latest result of a physics as SVG or PNG (the extension of `path` decides). 2D kinds (`gallery`, `components`) align each mode's phase and draw on the mesh section with region outlines and in-plane arrows. 3D kinds (`slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d`) evaluate `expr` (e.g. `re(Ex)`, `abs2(E)`) directly on un-averaged 8-subtet P2 `CellFields` so normal-field discontinuities across material interfaces stay sharp, render with Z-buffered/translucent orthographic projection (`camera`, `show_edges`), and support zero-copy Bloch unit-cell tiling (`repeat_periodic` with `exp(-i m periodic_phase)` applied before `expr`).", Read,
                schema(vec![("physics", phys()), ("path", text("file to write, .svg or .png")),
                    ("kind", json!({ "type": "string", "enum": ["gallery", "components", "slice3d", "cutaway3d", "isosurface3d", "mesh3d"], "description": "figure layout: gallery (default), components, slice3d (3D planar slice), cutaway3d (3D clipped cutaway), isosurface3d (3D translucent marching-tets isosurfaces) or mesh3d (3D region/element surface view)" })),
                    ("modes", json!({ "type": "array", "items": { "type": "integer", "minimum": 1 }, "description": "the modes to draw, numbered from 1; default the first max_modes" })),
                    ("max_modes", json!({ "type": "integer", "minimum": 1, "description": "how many modes when `modes` is not given, default 12" })),
                    ("columns", json!({ "type": "integer", "minimum": 1, "maximum": 8, "description": "panels per row of a gallery, default 4" })),
                    ("cmap", json!({ "type": "string", "enum": ["inferno", "viridis", "cividis", "coolwarm"], "description": "colour map: inferno (default for unsigned magnitudes), viridis, cividis, coolwarm (default symmetric map for signed 3D formulas)" })),
                    ("warp", json!({ "type": "number", "minimum": 0, "description": "largest motion of the deformed shape as a fraction of the window width; default 0.05 for a displacement (vector field `u`) and 0 for any other field; 0 draws the undeformed section" })),
                    ("arrows", boolean("arrows of the in-plane motion, default true")),
                    ("window", json!({ "type": "array", "items": { "type": "number" }, "description": "[x0, x1, z0, z1] in metres; default the extent of the modes" })),
                    ("at", json!({ "type": "number", "description": "the y of the plane a 3D mesh is cut on (or slice3d axis offset), in metres; default its first plane or midpoint" })),
                    ("core", json!({ "type": "array", "items": { "type": "integer" }, "description": "domain numbers that carry the mode (the guide, not the substrate): the phase and the scale are taken over them; default all" })),
                    ("family", text("the vector field to draw, by its prefix (u, E); default the first")),
                    ("expr", text("3D scalar formula string evaluated per-cell without vertex averaging, e.g. \"re(Ex)\", \"abs2(E)\", \"im(Hz)\"; default |u| or |E|")),
                    ("plane", json!({ "oneOf": [{ "type": "string", "enum": ["x", "y", "z"] }, { "type": "object" }], "description": "slice plane for slice3d: \"x\"/\"y\"/\"z\", { \"axis\": \"x\"|\"y\"|\"z\", \"at\": metres }, or { \"origin\": [x, y, z], \"normal\": [nx, ny, nz] }" })),
                    ("clip", json!({ "oneOf": [{ "type": "object" }, { "type": "array", "items": { "type": ["number", "null"] }, "minItems": 6, "maxItems": 6 }], "description": "3D cutaway clip bounds in metres: { xmin, xmax, ymin, ymax, zmin, zmax } (or x/y/z [min, max] pairs) or [xmin, xmax, ymin, ymax, zmin, zmax]" })),
                    ("isosurfaces", json!({ "type": "array", "items": { "type": "number" }, "description": "threshold levels for isosurface3d (default [0.25, 0.6] of peak, or [-0.6, -0.25, 0.25, 0.6] on signed formulas)" })),
                    ("repeat_periodic", json!({ "type": "integer", "minimum": 1, "description": "number of unit-cell copies to tile along periodic_shift with Bloch phase exp(-i m periodic_phase) applied before expr, default 1" })),
                    ("periodic_shift", json!({ "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3, "description": "explicit [dx, dy, dz] lattice translation in metres for repeat_periodic; auto-detected from periodic faces or y-extent when omitted" })),
                    ("periodic_phase", number("explicit Bloch phase advance per period in radians for repeat_periodic; auto-detected from result scalars or periodic BC when omitted")),
                    ("show_edges", json!({ "oneOf": [{ "type": "string", "enum": ["none", "regions", "mesh"] }, { "type": "boolean" }], "description": "3D wireframe overlay: \"regions\" (default for field views), \"mesh\" (default for mesh3d, or true) or \"none\" (or false)" })),
                    ("camera", json!({ "oneOf": [{ "type": "object", "properties": { "azimuth_deg": { "type": "number" }, "elevation_deg": { "type": "number" }, "azimuth": { "type": "number" }, "elevation": { "type": "number" } } }, { "type": "array", "items": { "type": "number" }, "minItems": 2, "maxItems": 2 }], "description": "3D orthographic camera angles in degrees: { azimuth_deg, elevation_deg } or [azimuth_deg, elevation_deg]" })),
                    ("panel_width", json!({ "type": "number", "description": "width of a panel in drawing units, default 250" })),
                    ("scale", json!({ "type": "number", "description": "pixels per drawing unit of a PNG, default 2" }))], &["physics", "path"])),
            cmd("store_save", "Save into a workspace", "Binds the project to a workflow folder (created in `workspace` the first time) or stores a new version of it. Later solves are recorded in the run table.", Expensive,
                schema(vec![("workspace", text("workspace folder, needed the first time"))], &[])),
            cmd("store_open", "Open a workflow", "Opens a workflow folder and makes its experiment the project.", Edit, schema(vec![("folder", text("workflow folder")), ("experiment", text("experiment id, default the first"))], &["folder"])),
            cmd("store_runs", "Run table", "The runs of the bound experiment: status, time, params and metrics.", Read, schema(vec![], &[])),
        ];
        for d in gd::all() {
            let mut s = gd::schema(&d);
            if let Some(p) = s.get_mut("properties").and_then(Value::as_object_mut) {
                p.insert("name".into(), text("a name for the feature; default the kind and a number (`rect1`)"));
            }
            c.push(cmd(
                &gd::tool_name(d.kind),
                d.title,
                &format!("{} [{}; dimensions {:?}]", d.doc, d.group.name(), d.dims),
                Edit,
                s,
            ));
        }
        c
    }

    // --------------------------------------------------------------------------------------------------------------------- dispatch

    /// Runs the command `name` with JSON arguments and returns its JSON reply.
    pub fn call(&mut self, name: &str, args: &Value) -> R {
        let mut a = args_obj(args)?.clone();
        if let Some(root) = &self.path_root {
            crate::path_guard::confine_args(root, name, &mut a)?;
        }
        let args = &Value::Object(a.clone());
        if let Some(kind) = geom_kind(name) {
            let mut fields = a.clone();
            let given = fields.remove("name");
            let n = geometry::add_feature(
                &mut self.project.geometry,
                kind,
                given.as_ref().and_then(Value::as_str),
                &Value::Object(fields),
            )?;
            return Ok(json!({ "name": n, "features": self.project.geometry.features.len() }));
        }
        match name {
            "list_commands" => Ok(json!(self.catalog().iter().map(|c| json!({ "name": c.name, "title": c.title, "access": c.access })).collect::<Vec<_>>())),
            "project_new" => {
                let dim = args.get("dim").and_then(Value::as_u64).unwrap_or(2);
                if dim != 2 && dim != 3 {
                    return Err(Error::invalid("`dim` is 2 or 3"));
                }
                let mut p = Project::new(opt_str(args, "title").unwrap_or("Untitled"), dim as u8);
                if let Some(h) = args.get("h").and_then(Value::as_f64) {
                    p.mesh.h = h;
                }
                self.set_project(p);
                ok()
            }
            "project_info" => Ok(self.info()),
            "project_json" => Ok(json!({ "json": self.project.to_json() })),
            "project_script" => {
                let e = crate::scriptout::export(&self.project)?;
                Ok(json!({ "script": e.script, "form": e.form, "commands": e.commands }))
            }
            "project_load" => self.load(Project::from_json(get_str(args, "json")?)?),
            "project_open" => {
                let path = get_str(args, "path")?;
                let t = std::fs::read_to_string(path).map_err(|e| Error::invalid(format!("{path}: {e}")))?;
                self.load(Project::from_json(&t)?)
            }
            "project_save" => {
                let path = get_str(args, "path")?;
                std::fs::write(path, self.project.to_json()).map_err(|e| Error::invalid(format!("{path}: {e}")))?;
                ok()
            }
            "list_operations" => Ok(json!(gd::all().iter().map(|d| json!({
                "kind": d.kind, "tool": gd::tool_name(d.kind), "title": d.title, "doc": d.doc, "group": d.group.name(), "dims": d.dims,
                "fields": d.fields.iter().map(|f| json!({ "name": f.name, "type": format!("{:?}", f.ty), "required": f.required, "default": f.default, "doc": f.doc })).collect::<Vec<_>>(),
            })).collect::<Vec<_>>())),
            "list_features" => Ok(geometry::describe(&self.project.geometry)),
            "set_geometry" => {
                let g: geom::model::FeatureList = serde_json::from_value(args.get("geometry").cloned().unwrap_or(Value::Null))
                    .map_err(|e| Error::invalid(format!("`geometry` is not a feature list: {e}")))?;
                geometry::set_all(&mut self.project.geometry, g)?;
                ok()
            }
            "set_problem" => {
                let id = get_str(args, "physics")?.to_string();
                let p = crate::physics_for(&self.registry, &id)?;
                let spec: ProblemSpec = serde_json::from_value(args.get("problem").cloned().unwrap_or(Value::Null))
                    .map_err(|e| Error::invalid(format!("`problem` is not a problem: {e}")))?;
                if spec.physics != id {
                    return Err(Error::invalid(format!("the problem is for '{}', not '{id}'", spec.physics)));
                }
                nyx_physics::validate_spec(p, &spec)?;
                match self.project.spec_mut(&id) {
                    Some(s) => *s = spec,
                    None => {
                        self.project.physics.push(spec);
                        if self.project.active.is_none() {
                            self.project.active = Some(id);
                        }
                    }
                }
                ok()
            }
            "update_feature" => {
                geometry::update_feature(&mut self.project.geometry, get_str(args, "name")?, args.get("fields").unwrap_or(&Value::Null))?;
                ok()
            }
            "remove_feature" => {
                geometry::remove_feature(&mut self.project.geometry, get_str(args, "name")?)?;
                ok()
            }
            "rename_feature" => {
                geometry::rename_feature(&mut self.project.geometry, get_str(args, "name")?, get_str(args, "new_name")?)?;
                ok()
            }
            "suppress_feature" => {
                geometry::suppress_feature(&mut self.project.geometry, get_str(args, "name")?, opt_bool(args, "suppressed", true))?;
                ok()
            }
            "tag_feature" => {
                geometry::tag_feature(&mut self.project.geometry, get_str(args, "name")?, get_str(args, "tag")?, opt_bool(args, "remove", false))?;
                ok()
            }
            "set_param" => {
                geometry::set_param(&mut self.project.geometry, get_str(args, "name")?, args.get("value").unwrap_or(&Value::Null))?;
                ok()
            }
            "remove_param" => {
                geometry::remove_param(&mut self.project.geometry, get_str(args, "name")?)?;
                ok()
            }
            "geometry_names" => geometry::names(&self.project.geometry),
            "list_physics" => Ok(json!(self.registry.ids().iter().map(|id| {
                let p = self.registry.get(id).expect("listed");
                json!({ "physics": id, "name": p.name(), "studies": p.studies().iter().map(|s| s.kind.clone()).collect::<Vec<_>>(), "in_project": self.project.spec(id).is_some() })
            }).collect::<Vec<_>>())),
            "describe_physics" => {
                let p = crate::physics_for(&self.registry, get_str(args, "physics")?)?;
                Ok(contract::descriptor_snapshot(p))
            }
            "add_physics" => self.add_physics(args),
            "set_study" => {
                let (kind, order) = (get_str(args, "kind")?.to_string(), args.get("order").and_then(Value::as_u64));
                let params = args.get("params").cloned().unwrap_or(Value::Null);
                let spec = self.project.need_spec(get_str(args, "physics")?)?;
                spec.study = StudySpec { kind, params };
                if let Some(o) = order {
                    spec.order = o as u8;
                }
                ok()
            }
            "set_material" => {
                let region = get_str(args, "region")?.to_string();
                self.check_region(&region)?;
                let m: MaterialAssign = serde_json::from_value(args.get("material").cloned().unwrap_or(Value::Null))
                    .map_err(|_| Error::invalid("`material` is a material name (text), {\"material\", \"euler_deg\"}, or an object of scalar or tensor properties"))?;
                self.project.need_spec(get_str(args, "physics")?)?.materials.insert(region, m);
                ok()
            }
            "add_condition" => {
                let bc = BcSpec {
                    on: get_str(args, "on")?.into(),
                    kind: get_str(args, "kind")?.into(),
                    params: args.get("params").cloned().unwrap_or(Value::Null),
                };
                let spec = self.project.need_spec(get_str(args, "physics")?)?;
                spec.bcs.push(bc);
                Ok(json!({ "index": spec.bcs.len() - 1 }))
            }
            "remove_condition" => {
                let i = get_u64(args, "index")? as usize;
                let spec = self.project.need_spec(get_str(args, "physics")?)?;
                if i >= spec.bcs.len() {
                    return Err(Error::invalid(format!("there is no condition {i} (there are {})", spec.bcs.len())));
                }
                spec.bcs.remove(i);
                ok()
            }
            "set_physics_params" => {
                let params = args.get("params").cloned().unwrap_or(Value::Null);
                self.project.need_spec(get_str(args, "physics")?)?.params = params;
                ok()
            }
            "remove_physics" => {
                let id = get_str(args, "physics")?.to_string();
                self.project.need_spec(&id)?;
                self.project.physics.retain(|s| s.physics != id);
                if self.project.active.as_deref() == Some(&id) {
                    self.project.active = self.project.physics.first().map(|s| s.physics.clone());
                }
                ok()
            }
            "set_active" => {
                let id = get_str(args, "physics")?.to_string();
                self.project.need_spec(&id)?;
                self.project.active = Some(id);
                ok()
            }
            "check_project" => {
                self.project.check(&self.registry)?;
                ok()
            }
            "import_gds" => self.import_gds(args),
            "add_chip_stack" => self.add_chip_stack(args),
            "set_mesh" => self.set_mesh(args),
            "mesh_build" => {
                if opt_bool(args, "wait", true) {
                    self.mesh_now()
                } else {
                    Ok(json!({ "job": self.start_mesh()? }))
                }
            }
            "solve_run" => {
                let targets = match args.get("physics") {
                    None | Some(Value::Null) => None,
                    Some(Value::String(s)) => Some(vec![s.clone()]),
                    Some(Value::Array(a)) => Some(a.iter().filter_map(|x| x.as_str().map(String::from)).collect()),
                    Some(_) => return Err(Error::invalid("`physics` is a code or a list of codes")),
                };
                let opts = QueueOpts {
                    isolate: args.get("isolate").and_then(Value::as_bool),
                    force: opt_bool(args, "force", false),
                    shared: args.get("shared").and_then(Value::as_bool),
                    save_fields: crate::queue::save_fields_of(args)?,
                    post_analyses: crate::queue::post_analyses_of(args)?,
                    memory_wait_timeout_s: crate::queue::memory_wait_timeout_of(args)?,
                    retain_amr_data: crate::queue::retain_amr_data_of(args)?,
                    ..Default::default()
                };
                let id = self.start_solve_with(targets, &opts)?;
                if opt_bool(args, "wait", true) {
                    self.wait(id)
                } else {
                    Ok(json!({ "job": id }))
                }
            }
            "job_status" => {
                let id = get_u64(args, "id")?;
                self.jobs.info(id).map(|i| json!(i)).ok_or_else(|| Error::invalid(format!("no job {id}")))
            }
            "job_list" => Ok(json!(self.jobs.list())),
            "job_cancel" => Ok(json!({ "cancelled": self.jobs.cancel(get_u64(args, "id")?) })),
            "job_wait" => match args.get("timeout_s") {
                None | Some(Value::Null) => self.wait(get_u64(args, "id")?),
                Some(v) => {
                    let t = v
                        .as_f64()
                        .filter(|&x| x >= 0.0 && x.is_finite())
                        .ok_or_else(|| Error::invalid("`timeout_s` is a finite number of seconds, zero or more"))?;
                    let dur = std::time::Duration::try_from_secs_f64(t)
                        .map(|d| d.min(std::time::Duration::from_secs(86_400)))
                        .map_err(|_| Error::invalid("`timeout_s` overflows Duration"))?;
                    self.wait_timeout(get_u64(args, "id")?, dur)
                }
            },
            "queue_submit" | "queue_list" | "queue_get" | "queue_cancel" | "queue_hold" | "queue_release" | "queue_set_priority" | "queue_pause"
            | "queue_resume" | "queue_wait" => self.queue_command(name, args),
            "sweep_run" | "sweep_table" | "sweep_cancel" => self.sweep_command(name, args),
            "poll_jobs" => Ok(json!(self
                .poll()
                .into_iter()
                .map(|(id, r)| match r {
                    Ok(v) => json!({ "job": id, "ok": v }),
                    Err(e) => json!({ "job": id, "error": { "code": e.code(), "message": e.to_string() } }),
                })
                .collect::<Vec<_>>())),
            "result_info" => self.result_info(opt_str(args, "physics")),
            "result_field" => self.result_field(args),
            "field" => self.cmd_field(args),
            "field_op" => self.cmd_field_op(args),
            "integrate" => self.cmd_integrate(args),
            "extremum" => self.cmd_extremum(args),
            "probe" => self.cmd_probe(args),
            "line_cut" => self.cmd_line_cut(args),
            "overlap" => self.cmd_overlap(args),
            "confinement" => self.cmd_confinement(args),
            "get_guide" => {
                let topic = opt_str(args, "topic").unwrap_or("nyx-instructions");
                let text = crate::guide::card(topic).ok_or_else(|| {
                    Error::invalid(format!("there is no guide '{topic}' (they are: {})", crate::guide::TOPICS.join(", ")))
                })?;
                Ok(json!({ "topic": topic, "text": text }))
            }
            "add_absorbing_layer" => {
                let mut fields = a.clone();
                let given = fields.remove("name");
                fields.entry("mode").or_insert_with(|| json!("Outward"));
                fields.entry("coord").or_insert_with(|| json!("Auto"));
                fields.entry("thickness").or_insert_with(|| json!("auto"));
                fields.entry("layers").or_insert_with(|| json!("auto"));
                fields.entry("grading").or_insert_with(|| json!(2.0));
                fields.entry("layer_grading").or_insert_with(|| json!(2.0));
                let n = geometry::add_feature(
                    &mut self.project.geometry,
                    "Pml",
                    given.as_ref().and_then(Value::as_str),
                    &Value::Object(fields),
                )?;
                Ok(json!({ "name": n, "features": self.project.geometry.features.len() }))
            }
            "add_pml_conditions" => {
                let physics = get_str(args, "physics")?;
                let profile = opt_str(args, "profile");
                let auto_scale = args.get("auto_scale").and_then(Value::as_bool);
                let target_r0 = args.get("target_reflection").and_then(Value::as_f64);
                let strength_arg = args.get("strength").and_then(Value::as_f64);
                if profile.is_none() && auto_scale.is_none() && target_r0.is_none() {
                    self.add_pml_conditions(physics, strength_arg.unwrap_or(3.0))
                } else {
                    self.add_pml_conditions_with(physics, strength_arg, profile, auto_scale, target_r0)
                }
            }
            "pml_diagnostics" => self.pml_diagnostics(),
            "list_materials" => self.list_materials(),
            "describe_material" => self.describe_material(get_str(args, "name")?),
            "add_material" => self.add_material(get_str(args, "name")?, get_str(args, "toml")?),
            "remove_material" => {
                self.remove_material(get_str(args, "name")?)?;
                ok()
            }
            "list_remote_materials" => self.list_remote_materials(opt_str(args, "domain")),
            "import_material" => self.import_material(
                get_str(args, "name")?,
                get_str(args, "format")?,
                get_str(args, "content")?,
                opt_str(args, "description"),
                opt_str(args, "extends"),
                opt_bool(args, "add_to_project", true),
                opt_str(args, "folder").map(std::path::Path::new),
            ),
            "fetch_material" => self.fetch_material_with_key(
                get_str(args, "source")?,
                opt_str(args, "target"),
                opt_str(args, "name"),
                opt_str(args, "description"),
                opt_str(args, "extends"),
                opt_bool(args, "add_to_project", true),
                opt_str(args, "folder").map(std::path::Path::new),
                opt_str(args, "api_key"),
            ),
            "search_remote_materials" => self.search_remote_materials(
                get_str(args, "query")?,
                opt_str(args, "domain"),
                opt_str(args, "source"),
                opt_bool(args, "live", false),
                opt_str(args, "api_key"),
            ),
            "evaluate_material" => {
                let euler = match args.get("euler_deg").and_then(Value::as_array) {
                    Some(a) if a.len() == 3 => Some([
                        a[0].as_f64().unwrap_or(0.0),
                        a[1].as_f64().unwrap_or(0.0),
                        a[2].as_f64().unwrap_or(0.0),
                    ]),
                    _ => None,
                };
                self.evaluate_material(
                    get_str(args, "name")?,
                    args.get("temperature").and_then(Value::as_f64),
                    args.get("frequency").and_then(Value::as_f64),
                    args.get("wavelength").and_then(Value::as_f64),
                    args.get("e_field").and_then(Value::as_f64),
                    args.get("thickness").and_then(Value::as_f64),
                    opt_str(args, "film_model"),
                    euler,
                )
            }
            "bcs_surface_impedance" => self.bcs_surface_impedance_cmd(args),
            "convergence_check" => {
                let (physics, factor, tolerance) = (
                    get_str(args, "physics")?,
                    args.get("factor").and_then(Value::as_f64).unwrap_or(0.6),
                    args.get("tolerance").and_then(Value::as_f64),
                );
                let matching = crate::convergence::Matching::from_args(args)?;
                let limits = crate::convergence::Limits::from_args(args)?;
                if opt_bool(args, "dry_run", false) {
                    // what the check would cost: both meshes are built, nothing is solved
                    return self.convergence_job().with_limits(limits).plan(physics, factor);
                }
                let start = |e: &mut Engine| {
                    e.start_convergence_limited(physics, factor, tolerance, &matching, &crate::queue::QueueOpts::default(), limits.clone())
                };
                if opt_bool(args, "wait", true) {
                    // through the queue, like a solve: it waits for the heavy lane and reuses the coarse solve of an earlier request
                    let id = start(self)?;
                    self.wait(id)
                } else {
                    Ok(json!({ "job": start(self)? }))
                }
            }
            "export_results" => {
                let parallel = args.get("parallel").and_then(Value::as_bool).unwrap_or(false);
                if opt_bool(args, "wait", true) {
                    self.export_results(get_str(args, "path")?, parallel)
                } else {
                    Ok(json!({ "job": self.start_export_results(get_str(args, "path")?, parallel)? }))
                }
            }
            "results_file_info" => self.results_file_info(get_str(args, "path")?),
            "open_results" => {
                if opt_bool(args, "wait", true) {
                    self.open_results(get_str(args, "path")?)
                } else {
                    Ok(json!({ "job": self.start_open_results(get_str(args, "path")?)? }))
                }
            }
            "figure" => self.figure(args),
            "list_analyses" => Ok(self.list_analyses()),
            "analyse" => self.analyse_with(
                get_str(args, "analysis")?,
                get_str(args, "physics")?,
                args.get("params").unwrap_or(&Value::Null),
            ),
            "store_save" => self.save_to_workspace(opt_str(args, "workspace").map(Into::into)),
            "store_open" => {
                let folder = get_str(args, "folder")?.to_string();
                self.open_workflow(std::path::Path::new(&folder), opt_str(args, "experiment"))
            }
            "store_runs" => match &self.link {
                None => Err(Error::invalid("the project is not bound to a workflow; use store_save first")),
                Some(l) => Ok(json!(l.runs())),
            },
            other => Err(Error::Unsupported(format!("there is no command '{other}'; list_commands lists them"))),
        }
    }

    // --------------------------------------------------------------------------------------------------------------------- parts

    fn load(&mut self, p: Project) -> R {
        p.check(&self.registry)?;
        self.set_project(p);
        ok()
    }

    fn info(&self) -> Value {
        let p = &self.project;
        json!({
            "title": p.title,
            "dim": p.geometry.dim,
            "features": p.geometry.features.len(),
            "params": p.geometry.params.iter().map(|q| json!({ "name": q.name, "value": q.value })).collect::<Vec<_>>(),
            "mesh_settings": p.mesh,
            "mesh": self.model().map(|m| m.info()),
            "active": p.active,
            "physics": p.physics.iter().map(|s| json!({
                "physics": s.physics, "study": s.study.kind, "order": s.order,
                "materials": s.materials.keys().collect::<Vec<_>>(), "conditions": s.bcs,
            })).collect::<Vec<_>>(),
            "results": self.results.iter().map(|(k, s)| json!({
                "physics": k, "stale": self.is_stale(k), "age_s": self.result_age(k), "validated": self.validated_json(k),
                "diagnostics": s.result.diagnostics, "error_budget": s.result.error_budget,
            })).collect::<Vec<_>>(),
            "store": self.link.as_ref().map(|l| json!({ "workflow": l.dir.wf.id, "folder": l.dir.root, "experiment": l.experiment, "version": l.dir.wf.version })),
            "jobs_running": self.jobs.running(),
        })
    }

    fn add_physics(&mut self, args: &Value) -> R {
        let id = get_str(args, "physics")?.to_string();
        let p = crate::physics_for(&self.registry, &id)?;
        if self.project.spec(&id).is_some() {
            return Err(Error::invalid(format!(
                "the project already has '{id}'; change it with set_study and the other commands"
            )));
        }
        if self.project.physics.len() >= 3 {
            return Err(Error::invalid("a project loads one to three physics"));
        }
        let kind = match opt_str(args, "study") {
            Some(k) => k.to_string(),
            None => p
                .studies()
                .first()
                .map(|s| s.kind.clone())
                .ok_or_else(|| Error::invalid(format!("{id} lists no study")))?,
        };
        let order = args.get("order").and_then(Value::as_u64).unwrap_or(1) as u8;
        let mut spec = ProblemSpec::new(&id, order, &kind);
        spec.study.params = args.get("study_params").cloned().unwrap_or(Value::Null);
        self.project.physics.push(spec);
        if self.project.active.is_none() {
            self.project.active = Some(id.clone());
        }
        Ok(json!({ "physics": id, "study": kind, "order": order }))
    }

    /// Refuses a region name that the overlay of the geometry does not have (names of absorbing pieces, with a dot, are not checked).
    fn check_region(&self, region: &str) -> Result<(), Error> {
        let names: Vec<&str> = self
            .project
            .geometry
            .features
            .iter()
            .filter_map(|f| match &f.op {
                geom::model::Op::Overlay { regions, .. } => Some(regions.iter().map(|r| r.name.as_str())),
                _ => None,
            })
            .flatten()
            .collect();
        let base = region.trim_end_matches('*');
        if names.is_empty() || region.contains('.') || region.contains(',') || names.contains(&base) {
            return Ok(());
        }
        Err(Error::invalid(format!(
            "the geometry has no region '{region}' (the regions are: {}); add the overlay first",
            names.join(", ")
        )))
    }

    fn set_mesh(&mut self, args: &Value) -> R {
        let a = args_obj(args)?;
        let mut s: MeshSettings = self.project.mesh.clone();
        for (k, v) in a {
            match k.as_str() {
                "h" => s.h = v.as_f64().ok_or_else(|| Error::invalid("`h` is a number in metres"))?,
                "adapt_scale" => {
                    s.adapt_scale = if v.is_null() {
                        None
                    } else {
                        Some(v.as_f64().ok_or_else(|| Error::invalid("`adapt_scale` is a number or null"))?)
                    }
                }
                "physics_sizing" => {
                    s.physics_sizing = v.as_bool().ok_or_else(|| Error::invalid("`physics_sizing` is true or false"))?
                }
                "points_per_wavelength" => {
                    s.points_per_wavelength = if v.is_null() {
                        None
                    } else {
                        Some(v.as_f64().ok_or_else(|| Error::invalid("`points_per_wavelength` is a number or null"))?)
                    }
                }
                "pml_mode" => {
                    s.pml_mode = if v.is_null() {
                        crate::project::PmlModeRequest::Auto
                    } else {
                        serde_json::from_value(v.clone()).map_err(|_| {
                            Error::invalid("`pml_mode` must be \"auto\", \"structured\", \"inline\", or null")
                        })?
                    };
                    if !a.contains_key("structured_pml") {
                        s.structured_pml = None;
                    }
                }
                "structured_pml" => {
                    if v.is_null() {
                        s.structured_pml = None;
                    } else {
                        let b = v
                            .as_bool()
                            .ok_or_else(|| Error::invalid("`structured_pml` is true, false, or null"))?;
                        s.structured_pml = Some(b);
                        if !a.contains_key("pml_mode") {
                            s.pml_mode = crate::project::PmlModeRequest::Auto;
                        }
                    }
                }
                "max_dev" => {
                    s.max_dev = if v.is_null() {
                        None
                    } else {
                        Some(v.as_f64().ok_or_else(|| Error::invalid("`max_dev` is a number in metres or null"))?)
                    }
                }
                "curved" => s.curved = v.as_bool().ok_or_else(|| Error::invalid("`curved` is true or false"))?,
                "gap_cells" => {
                    s.gap_cells = if v.is_null() {
                        None
                    } else {
                        let gc = v.as_f64().ok_or_else(|| Error::invalid("`gap_cells` is a positive number or null"))?;
                        if !(gc > 0.0 && gc.is_finite()) {
                            return Err(Error::invalid("`gap_cells` must be a positive finite number"));
                        }
                        Some(gc)
                    }
                }
                "z_coarsen" => {
                    s.z_coarsen = v.as_bool().ok_or_else(|| Error::invalid("`z_coarsen` is true or false"))?
                }
                "adapt_grade" => {
                    s.adapt_grade = if v.is_null() {
                        None
                    } else {
                        let gr = v
                            .as_f64()
                            .ok_or_else(|| Error::invalid("`adapt_grade` is a positive number or null"))?;
                        if !(gr > 0.0 && gr.is_finite()) {
                            return Err(Error::invalid("`adapt_grade` must be a positive finite number"));
                        }
                        Some(gr)
                    }
                }
                "adapt_max_cells" => {
                    s.adapt_max_cells = if v.is_null() {
                        None
                    } else {
                        let mc = v
                            .as_u64()
                            .filter(|&n| n >= 1)
                            .ok_or_else(|| Error::invalid("`adapt_max_cells` is a positive integer (>= 1) or null"))?;
                        let mc = usize::try_from(mc).map_err(|_| Error::invalid("`adapt_max_cells` is too large"))?;
                        Some(mc)
                    }
                }
                "max_cells" => {
                    s.max_cells = if v.is_null() {
                        None
                    } else {
                        let mc = v.as_u64().ok_or_else(|| Error::invalid("`max_cells` is an integer >= 1 or null"))?;
                        let mc = usize::try_from(mc).map_err(|_| Error::invalid("`max_cells` is too large"))?;
                        if mc < 1 {
                            return Err(Error::invalid("`max_cells` must be an integer >= 1"));
                        }
                        Some(mc)
                    }
                }
                "dz_min" => {
                    s.dz_min = if v.is_null() {
                        None
                    } else {
                        let dz = v.as_f64().ok_or_else(|| Error::invalid("`dz_min` is a positive number in metres or null"))?;
                        if !(dz > 0.0 && dz.is_finite()) {
                            return Err(Error::invalid("`dz_min` must be a positive finite number"));
                        }
                        Some(dz)
                    }
                }
                "max_z_slabs" => {
                    s.max_z_slabs = if v.is_null() {
                        None
                    } else {
                        let ms = v.as_u64().ok_or_else(|| Error::invalid("`max_z_slabs` is an integer >= 1 or null"))?;
                        let ms = usize::try_from(ms).map_err(|_| Error::invalid("`max_z_slabs` is too large"))?;
                        if ms < 1 {
                            return Err(Error::invalid("`max_z_slabs` must be an integer >= 1"));
                        }
                        Some(ms)
                    }
                }
                "z_growth_ratio" => {
                    s.z_growth_ratio = if v.is_null() {
                        None
                    } else {
                        let zr = v.as_f64().ok_or_else(|| Error::invalid("`z_growth_ratio` is a number > 1.0 or null"))?;
                        if !(zr > 1.0 && zr.is_finite()) {
                            return Err(Error::invalid("`z_growth_ratio` must be a finite number strictly greater than 1.0"));
                        }
                        Some(zr)
                    }
                }
                "max_steiner_multiplier" => {
                    s.max_steiner_multiplier = if v.is_null() {
                        None
                    } else {
                        let mm = v.as_u64().ok_or_else(|| Error::invalid("`max_steiner_multiplier` is a non-negative integer or null"))?;
                        let mm = usize::try_from(mm).map_err(|_| Error::invalid("`max_steiner_multiplier` is too large"))?;
                        Some(mm)
                    }
                }
                "strict_mesh" => {
                    s.strict_mesh = v.as_bool().ok_or_else(|| Error::invalid("`strict_mesh` is true or false"))?
                }
                "singular_edges" => {
                    s.singular_edges = if v.is_null() {
                        None
                    } else {
                        let se = v
                            .as_str()
                            .ok_or_else(|| Error::invalid("`singular_edges` is \"none\", \"quarter_point\", or null"))?;
                        if !matches!(se, "none" | "quarter_point") {
                            return Err(Error::invalid(format!(
                                "`singular_edges` must be \"none\" or \"quarter_point\", got \"{se}\""
                            )));
                        }
                        Some(se.to_string())
                    }
                }
                "quarter_point_xi" => {
                    s.quarter_point_xi = if v.is_null() {
                        None
                    } else {
                        let xi = v
                            .as_f64()
                            .filter(|&x| x.is_finite() && (0.25..=0.50).contains(&x))
                            .ok_or_else(|| Error::invalid("`quarter_point_xi` is a number in [0.25, 0.50] or null"))?;
                        Some(xi)
                    }
                }
                "corner_grading" => {
                    s.corner_grading = v.as_bool().ok_or_else(|| Error::invalid("`corner_grading` is true or false"))?
                }
                "target_tol" => {
                    s.target_tol = if v.is_null() {
                        None
                    } else {
                        let t = v
                            .as_f64()
                            .filter(|&x| x > 0.0 && x.is_finite())
                            .ok_or_else(|| Error::invalid("`target_tol` must be a finite number > 0 or null"))?;
                        Some(t)
                    }
                }
                "adapt_passes" => {
                    s.adapt_passes = if v.is_null() {
                        None
                    } else {
                        let p = v
                            .as_u64()
                            .and_then(|x| usize::try_from(x).ok())
                            .filter(|&x| x >= 1)
                            .ok_or_else(|| Error::invalid("`adapt_passes` must be an integer >= 1 or null"))?;
                        Some(p)
                    }
                }
                "adapt_rel_tol" => {
                    s.adapt_rel_tol = if v.is_null() {
                        None
                    } else {
                        let t = v
                            .as_f64()
                            .filter(|&x| x > 0.0 && x.is_finite())
                            .ok_or_else(|| Error::invalid("`adapt_rel_tol` must be a finite number > 0 or null"))?;
                        Some(t)
                    }
                }
                "adapt_growth" => {
                    s.adapt_growth = if v.is_null() {
                        None
                    } else {
                        let g = v
                            .as_f64()
                            .filter(|&x| x > 1.0 && x.is_finite())
                            .ok_or_else(|| Error::invalid("`adapt_growth` must be a finite number > 1.0 or null"))?;
                        Some(g)
                    }
                }
                "retain_amr_data" => {
                    s.retain_amr_data =
                        v.as_bool().ok_or_else(|| Error::invalid("`retain_amr_data` is true or false"))?
                }
                "sweep" => {
                    s.sweep = if v.is_null() {
                        None
                    } else {
                        Some(serde_json::from_value(v.clone()).map_err(|e| Error::invalid(format!("`sweep`: {e}")))?)
                    }
                }
                other => {
                    return Err(Error::invalid(format!(
                        "set_mesh has no field '{other}' (h, adapt_scale, adapt_grade, adapt_max_cells, physics_sizing, points_per_wavelength, pml_mode, structured_pml, max_dev, curved, gap_cells, z_coarsen, max_cells, dz_min, max_z_slabs, z_growth_ratio, max_steiner_multiplier, strict_mesh, singular_edges, quarter_point_xi, corner_grading, target_tol, adapt_passes, adapt_rel_tol, adapt_growth, retain_amr_data, sweep)"
                    )))
                }
            }
        }
        if !(s.h > 0.0 && s.h.is_finite()) {
            return Err(Error::invalid("the element size must be a positive finite number"));
        }
        self.project.mesh = s;
        ok()
    }

    fn result_info(&self, physics: Option<&str>) -> R {
        let rows: Vec<Value> = self
            .results
            .iter()
            .filter(|(k, _)| physics.is_none_or(|p| p == k.as_str()))
            .map(|(k, v)| {
                let r = &v.result;
                json!({
                    "physics": k, "study": r.study, "order": r.order, "dim": r.dim, "dofs": r.dofs, "seconds": r.seconds,
                    "stale": self.is_stale(k), "age_s": self.result_age(k), "validated": self.validated_json(k),
                    "diagnostics": r.diagnostics, "error_budget": r.error_budget,
                    "solutions": r.solutions.iter().enumerate().map(|(i, s)| json!({
                        "index": i, "label": s.label, "value": s.value, "eigenvalue": s.eigenvalue,
                        "fields": s.fields.iter().map(|f| json!({ "name": f.name, "unit": f.unit, "components": f.components, "n": f.values.len() })).collect::<Vec<_>>(),
                        "cell_fields": s.cell_fields.iter().map(|f| json!({ "name": f.name, "unit": f.unit, "complex": !f.im.is_empty() })).collect::<Vec<_>>(),
                    })).collect::<Vec<_>>(),
                    "scalars": r.scalars,
                })
            })
            .collect();
        if rows.is_empty() {
            return Err(Error::invalid(match physics {
                Some(p) => format!("there is no result for '{p}' yet; solve first"),
                None => "there are no results yet; solve first".into(),
            }));
        }
        Ok(json!(rows))
    }

    fn result_field(&self, args: &Value) -> R {
        let p = get_str(args, "physics")?;
        let stored = self
            .results
            .get(p)
            .ok_or_else(|| Error::invalid(format!("there is no result for '{p}' yet; solve first")))?;
        let i = args.get("solution").and_then(Value::as_u64).unwrap_or(0) as usize;
        let sol = stored.result.solutions.get(i).ok_or_else(|| {
            Error::invalid(format!("the result has {} solutions; there is no {i}", stored.result.solutions.len()))
        })?;
        let name = get_str(args, "field")?;
        let f = sol.field(name).ok_or_else(|| {
            Error::invalid(format!(
                "no field '{name}' (the fields are: {})",
                sol.fields.iter().map(|f| f.name.as_str()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        let max = args.get("max_values").and_then(Value::as_u64).unwrap_or(2000).max(1) as usize;
        let stride = f.values.len().div_ceil(max).max(1);
        let values: Vec<f64> = f.values.iter().step_by(stride).copied().collect();
        Ok(
            json!({ "name": f.name, "unit": f.unit, "components": f.components, "n": f.values.len(), "stride": stride, "values": values, "stale": self.is_stale(p), "age_s": self.result_age(p), "validated": self.validated_json(p) }),
        )
    }
}

/// The geometry operation an `add_*` command stands for.
fn geom_kind(name: &str) -> Option<&'static str> {
    gd::all().into_iter().find(|d| gd::tool_name(d.kind) == name).map(|d| d.kind)
}
