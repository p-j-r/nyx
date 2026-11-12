//! The commands: every action a front-end can take, by name, with JSON arguments, a JSON schema and a JSON reply. The window calls the typed
//! methods of `Engine` directly or these commands; the MCP server turns `catalog()` into its tool list and `call()` into its tool calls,
//! so a tool is a command and nothing else. The geometry tools are generated from `geomdesc` (gate G5: every operation, none picked by
//! hand), the physics-dependent schemas come from the registry.

use crate::engine::Engine;
use crate::geomdesc as gd;
use crate::geometry;
use crate::project::{MeshSettings, Project};
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
fn object(doc: &str) -> Value {
    json!({ "type": "object", "description": doc })
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
        let mut c = vec![
            cmd("list_commands", "List the commands", "Every command with its title, access and schema.", Read, schema(vec![], &[])),
            cmd("project_new", "New project", "Starts an empty project (dropping the current one, its mesh and results).", Edit,
                schema(vec![("title", text("project title")), ("dim", int("2 (x-z cross-section, the default) or 3")), ("h", number("global element size in metres"))], &[])),
            cmd("project_info", "Project summary", "Title, geometry size, the problem of each physics, the mesh, the results (and whether they are stale), the store binding and the jobs.", Read, schema(vec![], &[])),
            cmd("project_json", "Project as JSON", "The whole project as JSON text (the file that is saved).", Read, schema(vec![], &[])),
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
                schema(vec![("name", text("parameter name")), ("value", json!({ "oneOf": [{ "type": "number" }, { "type": "string" }] }))], &["name", "value"])),
            cmd("remove_param", "Remove a model parameter", "", Edit, schema(vec![("name", text("parameter name"))], &["name"])),
            cmd("geometry_names", "Evaluate the geometry and list its names", "Evaluates the model and returns the region and boundary names a problem can use, and the bounding box. Instant; it reports a geometry that does not evaluate.", Read, schema(vec![], &[])),
            cmd("list_physics", "List the physics of this build", "The physics compiled in, their studies, and which the project loads.", Read, schema(vec![], &[])),
            cmd("describe_physics", "Describe a physics", "Studies, conditions, properties, fields and couplings of one physics, and the JSON Schema of its problem.", Read, schema(vec![("physics", phys())], &["physics"])),
            cmd("add_physics", "Add a physics to the project", "Adds a problem for a physics (a project loads one to three). The study defaults to the first the physics lists.", Edit,
                schema(vec![("physics", phys()), ("study", text("study kind (see describe_physics)")), ("study_params", object("study parameters")), ("order", int("element order, default 1"))], &["physics"])),
            cmd("set_study", "Set the study of a physics", "", Edit,
                schema(vec![("physics", phys()), ("kind", text("study kind")), ("params", object("study parameters")), ("order", int("element order"))], &["physics", "kind"])),
            cmd("set_material", "Assign a material to a region", "A material by name (library or project: list_materials), a named material placed at an angle (`{\"material\": \"sapphire\", \"euler_deg\": [0, 57.6, 0]}`, Bunge angles in degrees), or scalar properties by descriptor name (`{\"E\": 1.7e11, \"nu\": 0.28}`).", Edit,
                schema(vec![("physics", phys()), ("region", text("region name")), ("material", json!({ "oneOf": [{ "type": "string" }, { "type": "object" }], "description": "material name or properties" }))], &["physics", "region", "material"])),
            cmd("add_condition", "Add a boundary or domain condition", "`on` is a selector over the geometry (`outer & plane(x = 0)`), a boundary name or a region name.", Edit,
                schema(vec![("physics", phys()), ("on", text("selector, boundary or region")), ("kind", text("condition kind (see describe_physics)")), ("params", object("condition parameters"))], &["physics", "on", "kind"])),
            cmd("remove_condition", "Remove a condition", "By its position (from 0) in the physics' list, as project_info shows.", Edit, schema(vec![("physics", phys()), ("index", int("position from 0"))], &["physics", "index"])),
            cmd("set_physics_params", "Set the options of a physics", "", Edit, schema(vec![("physics", phys()), ("params", object("physics options"))], &["physics", "params"])),
            cmd("remove_physics", "Remove a physics from the project", "", Edit, schema(vec![("physics", phys())], &["physics"])),
            cmd("set_active", "Set the active physics", "The one solved by default.", Edit, schema(vec![("physics", phys())], &["physics"])),
            cmd("check_project", "Check the project", "Checks the geometry, the element size and every problem against the descriptors of its physics, without meshing.", Read, schema(vec![], &[])),
            cmd("set_mesh", "Set the mesh settings", "Give only what changes. `adapt_scale` or `sweep` set to null return to the simple mesher / a 2D mesh.", Edit,
                schema(vec![("h", number("global element size in metres")), ("adapt_scale", json!({ "type": ["number", "null"] })), ("physics_sizing", boolean("size the elements from the waves of the physics (points per wavelength)")), ("points_per_wavelength", json!({ "type": ["number", "null"] })), ("structured_pml", boolean("rebuild the absorbing layers as structured layers")), ("max_dev", json!({ "type": ["number", "null"], "description": "largest chord deviation from a curved edge, metres" })),
                    ("sweep", json!({ "type": ["array", "null"], "items": { "type": "object", "properties": { "thickness": { "type": "number" }, "n": { "type": "integer" }, "domain_offset": { "type": "integer" } }, "required": ["thickness", "n"] }, "description": "layers to sweep the 2D mesh along y" }))], &[])),
            cmd("mesh_build", "Mesh the geometry", "Meshes the geometry with the current settings. With wait false it starts a job and returns its id.", Expensive, schema(vec![("wait", boolean("wait for the mesh, default true"))], &[])),
            cmd("solve_run", "Solve", "Solves the active physics (or the named ones, in coupling order), meshing first when needed. With wait false it starts a job. When the project is bound to a workflow the run is recorded.", Expensive,
                schema(vec![("physics", json!({ "oneOf": [{ "type": "string" }, { "type": "array", "items": { "type": "string" } }], "description": "physics code(s); default the active one" })), ("wait", boolean("wait for the result, default true"))], &[])),
            cmd("job_status", "Status of a job", "State, current stage, fraction done and elapsed time.", Read, schema(vec![("id", int("job id"))], &["id"])),
            cmd("job_list", "List the jobs", "", Read, schema(vec![], &[])),
            cmd("job_cancel", "Cancel a job", "The solver notices at its next check; a factorisation in progress finishes first.", Expensive, schema(vec![("id", int("job id"))], &["id"])),
            cmd("job_wait", "Wait for a job", "Blocks until the job has finished and applies its outcome.", Expensive, schema(vec![("id", int("job id"))], &["id"])),
            cmd("poll_jobs", "Apply finished jobs", "Applies every job that has finished since the last call and returns what each came to.", Read, schema(vec![], &[])),
            cmd("result_info", "Describe the results", "For each physics solved: the study, size, time, solutions with their values, fields and scalars, and whether the project changed since.", Read, schema(vec![("physics", phys())], &[])),
            cmd("result_field", "Read a field", "The values of one field of one solution (sampled evenly when there are more than `max_values`).", Read,
                schema(vec![("physics", phys()), ("solution", int("solution number from 0, default 0")), ("field", text("field name")), ("max_values", int("default 2000"))], &["physics", "field"])),
            cmd("get_guide", "Read a nyx guide", "Read a short guide. Call with topic nyx-instructions before using other nyx tools. Other topics: geometry, optical-modes (reading mode tables), troubleshooting (errors).", Read,
                schema(vec![("topic", json!({ "type": "string", "enum": crate::guide::TOPICS, "description": "which card; omit for nyx-instructions" }))], &[])),
            cmd("add_pml_conditions", "Add the absorbing-layer conditions", "Adds a `pml` condition for every absorbing piece the geometry made with add_pml (replacing earlier ones on those pieces), so the physics absorbs where the geometry says. `strength` is the imaginary stretch at the outer edge (default 3).", Edit,
                schema(vec![("physics", phys()), ("strength", number("imaginary stretch at the outer edge, default 3"))], &["physics"])),
            cmd("list_materials", "List the materials", "The materials a problem can name: those of the library (material files) and those the project carries.", Read, schema(vec![], &[])),
            cmd("describe_material", "Describe a material", "What a material gives (elastic class, optical model, thermal, electrical) and its TOML text.", Read, schema(vec![("name", text("material name"))], &["name"])),
            cmd("add_material", "Add a material to the project", "Adds (or replaces) a material the project carries, as TOML text (units are written in the file; `extends = \"silicon\"` starts from a library material). Read at once, so a mistake is reported here.", Edit,
                schema(vec![("name", text("material name")), ("toml", text("the material file text"))], &["name", "toml"])),
            cmd("remove_material", "Remove a project material", "", Edit, schema(vec![("name", text("material name"))], &["name"])),
            cmd("convergence_check", "Check mesh convergence", "Solves the physics again on a finer mesh (element size times `factor`; the project is not changed) and compares every eigenvalue, value and scalar with the present result. With `tolerance` each row says pass or fail and `converged` summarises. Use it before quoting a number: a mode index that moves by more than your tolerance is not converged.", Expensive,
                schema(vec![("physics", phys()), ("factor", number("finer/present element size, between 0 and 1, default 0.6")), ("tolerance", number("largest accepted relative difference, e.g. 1e-3; omit to only report"))], &["physics"])),
            cmd("export_results", "Write the results to an HDF5 file", "Writes the project, the mesh and the latest result of every solved physics (all fields, compressed, exact) into an HDF5 file and an XDMF description next to it for ParaView and VisIt. Needs a build with HDF5 support.", Expensive, schema(vec![("path", text("file path, .h5"))], &["path"])),
            cmd("results_file_info", "Describe an HDF5 results file", "The project title, the mesh size and for each physics the study, size, time, number of solutions and scalars of a results file.", Read, schema(vec![("path", text("file path"))], &["path"])),
            cmd("open_results", "Open an HDF5 results file", "Makes the project of a results file the project, builds its mesh again (meshing is deterministic) and checks it is the stored mesh, then takes the stored results as the latest results, so every analysis, plot and export works on them as after a solve. Replaces the present project. Needs a build with HDF5 support.", Edit, schema(vec![("path", text("file path, .h5"))], &["path"])),
            cmd("list_analyses", "List the analyses", "The analyses of results (polarization of the modes) and which solved physics each applies to.", Read, schema(vec![], &[])),
            cmd("analyse", "Analyse a result", "Runs an analysis on the latest result of a physics: the polarization of the modes gives each mode its class (quasi-TE, quasi-TM, hybrid), Stokes parameters, degree of polarization and the degenerate groups.", Read,
                schema(vec![("analysis", text("analysis id (list_analyses)")), ("physics", phys())], &["analysis", "physics"])),
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
        let a = args_obj(args)?;
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
                    .map_err(|_| Error::invalid("`material` is a material name (text) or an object of scalar properties"))?;
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
                if opt_bool(args, "wait", true) {
                    self.solve_now(targets)
                } else {
                    Ok(json!({ "job": self.start_solve(targets)? }))
                }
            }
            "job_status" => {
                let id = get_u64(args, "id")?;
                self.jobs.info(id).map(|i| json!(i)).ok_or_else(|| Error::invalid(format!("no job {id}")))
            }
            "job_list" => Ok(json!(self.jobs.list())),
            "job_cancel" => Ok(json!({ "cancelled": self.jobs.cancel(get_u64(args, "id")?) })),
            "job_wait" => self.wait(get_u64(args, "id")?),
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
            "get_guide" => {
                let topic = opt_str(args, "topic").unwrap_or("nyx-instructions");
                let text = crate::guide::card(topic).ok_or_else(|| {
                    Error::invalid(format!("there is no guide '{topic}' (they are: {})", crate::guide::TOPICS.join(", ")))
                })?;
                Ok(json!({ "topic": topic, "text": text }))
            }
            "add_pml_conditions" => self.add_pml_conditions(get_str(args, "physics")?, args.get("strength").and_then(Value::as_f64).unwrap_or(3.0)),
            "list_materials" => self.list_materials(),
            "describe_material" => self.describe_material(get_str(args, "name")?),
            "add_material" => self.add_material(get_str(args, "name")?, get_str(args, "toml")?),
            "remove_material" => {
                self.remove_material(get_str(args, "name")?)?;
                ok()
            }
            "convergence_check" => self.convergence_check(
                get_str(args, "physics")?,
                args.get("factor").and_then(Value::as_f64).unwrap_or(0.6),
                args.get("tolerance").and_then(Value::as_f64),
            ),
            "export_results" => self.export_results(get_str(args, "path")?),
            "results_file_info" => self.results_file_info(get_str(args, "path")?),
            "open_results" => self.open_results(get_str(args, "path")?),
            "list_analyses" => Ok(self.list_analyses()),
            "analyse" => self.analyse(get_str(args, "analysis")?, get_str(args, "physics")?),
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
            "results": self.results.iter().map(|(k, v)| json!({ "physics": k, "stale": v.stamp != self.stamp() })).collect::<Vec<_>>(),
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
                "sweep" => {
                    s.sweep = if v.is_null() {
                        None
                    } else {
                        Some(serde_json::from_value(v.clone()).map_err(|e| Error::invalid(format!("`sweep`: {e}")))?)
                    }
                }
                other => {
                    return Err(Error::invalid(format!(
                        "set_mesh has no field '{other}' (h, adapt_scale, physics_sizing, sweep)"
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
        let stamp = self.stamp();
        let rows: Vec<Value> = self
            .results
            .iter()
            .filter(|(k, _)| physics.is_none_or(|p| p == k.as_str()))
            .map(|(k, v)| {
                let r = &v.result;
                json!({
                    "physics": k, "study": r.study, "order": r.order, "dim": r.dim, "dofs": r.dofs, "seconds": r.seconds,
                    "stale": v.stamp != stamp,
                    "solutions": r.solutions.iter().enumerate().map(|(i, s)| json!({
                        "index": i, "label": s.label, "value": s.value, "eigenvalue": s.eigenvalue,
                        "fields": s.fields.iter().map(|f| json!({ "name": f.name, "unit": f.unit, "components": f.components, "n": f.values.len() })).collect::<Vec<_>>(),
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
            json!({ "name": f.name, "unit": f.unit, "components": f.components, "n": f.values.len(), "stride": stride, "values": values }),
        )
    }
}

/// The geometry operation an `add_*` command stands for.
fn geom_kind(name: &str) -> Option<&'static str> {
    gd::all().into_iter().find(|d| gd::tool_name(d.kind) == name).map(|d| d.kind)
}
