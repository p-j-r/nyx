# nyx workflows: dashboard, experiments, `.nyx` files and run tracking

Status: design only, 2026-11-01. Nothing in `app`, `femcore` or `geom` is changed by this document; every item that would touch them is marked **[needs app edit]** so you can decide when.

## 0. Summary

You described how nyx comes together: a **dashboard** that manages independent **workflows**; each workflow holds related **experiments** (simulations); each experiment has three sections, **Concept, Build, Analysis**; every workflow leaves a small searchable `.nyx` file; tracking of many runs should feel like MLflow but without Excel; everything Rust, fast, uncluttered, and ready for an AI to drive later.

My proposal in six lines:

1. **Four levels:** workspace (the folder behind the dashboard) > workflow > experiment > run. A *run* is one Compute (one point of a sweep is one run). This is MLflow's experiment/run split with one extra level on top, because you want related experiments grouped.
2. **A new pure-Rust crate `nyx-store`** (no egui, serde only) owns the manifest, run records, search and retention. The GUI, `nyx --batch`, scripts and the future MCP server all call it, so none of them can do something the others cannot. Building it needs **no edit** to `app`, `femcore` or `geom`.
3. **`workflow.nyx`** is a small JSON manifest (a few KB, a few hundred KB at the very most) holding, per experiment, an auto-generated **concept card** (geometry family and dimensions, materials, physics, study, wavelength) plus a **digest of its runs** (counts, parameter ranges, metric ranges). Nobody types metadata for it to be findable. The existing `.nyx.xml` stays as the full spec of an experiment.
4. **Runs live in append-only `runs.jsonl`**, one line each: parameters, metrics, status, hashes, links to heavy files. Heavy files (fields, HDF5) are deletable; metrics are kept forever.
5. **Dashboard rules:** pinned plus recently touched, hard cap, nothing is ever deleted by the rules, only moved to "Earlier" where search finds it. Empty drafts vanish on their own.
6. **Three sections replace six tabs.** The existing editors are reused; what changes is where they sit.

Section 10 is the feature-by-feature gap list with a suggested build order; section 11 lists the decisions I took that you may want to reverse.

## 1. Vocabulary

| Term | Meaning | Today's nearest thing |
|---|---|---|
| Workspace | folder `Documents\nyx\` (setting) holding all workflows | `Documents\nyx_runs` (flat autosave folder) |
| Workflow | an independent line of work with a title, tags, notes, status. One tab | none: the app holds exactly one `Project` |
| Experiment | one simulation definition inside a workflow: geometry, materials, physics, study. Can be forked from another | `state::Project` (`.nyx.xml`) |
| Run | one execution of an experiment's frozen spec, with parameters, metrics, log, artifacts | a `temp_untitled_<stamp>` autosave triplet |
| Concept card | auto-generated summary of an experiment's inputs; what search matches against | none |

"Related or connected" experiments are expressed two ways: **fork** (a copy remembering its parent, shown as a small tree) and **link** (any experiment to any experiment, with a one-word label such as "same geometry", "validates").

## 2. What exists today (checked in the source)

* `crates/app/src/app.rs`: `App` holds **one** `proj: Project`, one `doc: Doc` (mesh), one `results: Vec<ResultView>`, one viewport, one job slot (`job`, "only one job runs at a time"). `new_project` **replaces** the project ("nothing is kept of the old geometry"); File > New opens `NewDialog` (physics, 2D/3D, structured).
* `crates/app/src/state.rs`: `Project { version, geometry: GeometrySpec, source, active: Physics, specs: Vec<ProblemSpec>, colormap }`, written as `.nyx.xml` through `femcore::xmlio` (a `serde_json::Value` to XML bridge) or `.nyx.json`. `Tab` has six values: Geometry, Materials, Physics, Study, Results, Analysis.
* `crates/app/src/worker.rs`: every solve autosaves `temp_untitled_<stamp>.nyx.xml` + `.log.txt` (+ `.hdf5` with the feature) into `nyx_runs`, forever. This is the only run history and it is just a pile of files. `runs/rib` in the repo shows the cost: `run_sweep*.sh`, `all_gains.csv`, `summarize.py`, `plot_phonons.py` are a hand-made experiment tracker.
* `crates/app/src/tabs.rs`: `tab_analysis` is the polarization analysis only (Stokes parameters, Poincare sphere, table from `femcore::analysis`), shown when `result.analysis` exists. `tab_results` shows fields. The menu has "Python (optional)" which runs a script on project JSON plus result CSV, and `export_csv` for one result.
* `femcore::problem::SolveResult` has `scalars: Vec<(String, Scalar)>` (`Num` or `Text`), typed `Solution.neff`, `analysis`, `log`, `seconds`, `dofs`. Results are deterministic (README, `tests/determinism.rs`), so a hash of the inputs identifies a result.
* `geom::model::FeatureList` has named parameters (`with(&[("w", 450e-9)])`) and `validate()`; `Geometry::Model` carries one. This is where "300 nm" can be read from without parsing text.
* `nyx --batch project.nyx.xml` solves one project and writes an HDF5 file; it records nothing.
* `docs/SWEEP_CACHING_PROPOSAL.md` already proposes a project **parameter table**, a Sweep section and a results table (its section 5). `docs/MCP_ASSISTANT_PLAN.md` already proposes a session folder, `session.nyx.xml`, a content-hashed artifact cache and a `journal.jsonl` (its sections 8 and 21). This design is meant to absorb both rather than compete (section 8).
* Not present: any workflow or run concept, any index or search, any table widget (`egui_extras` is not a dependency; `egui_plot` is), any analysis-function registry, any handling of a `.nyx` extension (`open_path` accepts only xml, json, h5).

## 3. On disk

```text
Documents\nyx\                          workspace (setting; migrate nyx_runs into it once)
  workspace.json                        tiny: pinned order, open tabs, retention settings
  index.cache                           optional, rebuildable, safe to delete (section 7)
  0001_soi-te-modes_k3f9a2\             <slug>_<short id>; the folder never moves when the title changes
    workflow.nyx                        the manifest (JSON, KB)
    notes.md                            free text, optional
    e1_5xq\
      experiment.nyx.xml                the existing Project, unchanged
      runs.jsonl                        append-only run records (the MLflow table)
      specs\<hash>.nyx.xml              frozen specs, content addressed, de-duplicated
      runs\r0007_7hb\                   heavy and deletable
        log.txt  results.h5  analysis\...
```

Rules: ids are short, time-ordered, base32 (sorting by id is sorting by creation, no file reads). Meshes are never stored (rebuilt deterministically, as the MCP plan also concludes). Everything is written by write-temp-then-rename; `runs.jsonl` is written by single `O_APPEND` lines so a detached worker (MCP plan section 20) and the window can both add runs without a lock. The manifest's run digest is **derived** from `runs.jsonl` and rebuilt when missing or stale, so a crash never loses runs.

### 3.1 `workflow.nyx` (illustrative; field names are proposals)

```json
{
  "nyx": 1, "kind": "workflow", "id": "k3f9a2",
  "title": "Plasmon scattering on 300 nm SOS frame",
  "status": "active", "pinned": false, "origin": "gui",
  "created": "2026-11-01T00:08:15Z", "touched": "2026-11-02T09:41:03Z",
  "tags": ["plasmon", "scattering"],
  "experiments": [{
    "id": "e1", "title": "Frame, TE modes vs width", "parent": null, "links": [],
    "spec_hash": "b3:9d1c...",
    "card": {
      "physics": "emw", "study": "modes", "dim": 2, "wavelength_m": 1.55e-6,
      "geometry": {
        "family": "silicon-on-sapphire rib",
        "summary": "Si 300 nm x 500 nm on sapphire, air clad",
        "dims": [{"name": "thickness", "of": "Si", "m": 3.0e-7}, {"name": "w", "of": "Si", "m": 5.0e-7}]
      },
      "materials": ["Si", "Al2O3", "air"],
      "keywords": ["rib", "waveguide", "sapphire", "SOS"]
    },
    "runs": {
      "count": 21, "ok": 21, "failed": 0, "last": "2026-11-02T09:41:03Z",
      "params":  {"w": {"min": 4e-7, "max": 6e-7, "n": 21}},
      "metrics": {"neff_re": {"min": 2.31, "max": 2.62}, "te_fraction": {"min": 0.81, "max": 0.97}},
      "starred": ["r0012_x2k"]
    }
  }]
}
```

Size: card about 0.5 KB, run digest about 0.3 KB, so a 50-experiment workflow is well under 100 KB. The per-run detail never goes in here, which is how it stays small even for a 10,000-point sweep.

### 3.2 A run record (one line of `runs.jsonl`)

```json
{"id":"r0007_7hb","t":"2026-11-02T09:41:03Z","sec":0.11,"status":"ok","spec":"b3:9d1c...","femcore":"0.1.0",
 "params":{"w":5.0e-7,"lambda":1.55e-6,"h":4.0e-8,"dofs":18432},
 "metrics":{"neff_re":2.4412,"neff_im":1.2e-7,"te_fraction":0.93},
 "tables":{"modes":"runs/r0007_7hb/modes.json"},
 "tags":[], "star":false, "note":"",
 "artifacts":[{"kind":"log","path":"runs/r0007_7hb/log.txt"},{"kind":"fields","path":"runs/r0007_7hb/results.h5","bytes":4128000}]}
```

`params` are numbers with their SI value (units are in the name table, not repeated per run); `metrics` are scalars only; anything with rows (mode tables, spectra) goes to a small side file referenced under `tables`. The scalars come from typed fields (`Solution.neff`, analysis rows), not from the formatted `Scalar::Text` entries, which the MCP plan already flags as unsafe to parse.

### 3.3 The concept card is generated, never typed

From the experiment's own spec: the `FeatureList` parameter table and region bounding boxes (dims), region and material names, `Physics`, `Study` kind, wavelength or frequency, dimension, plus a fixed synonym table (sapphire = Al2O3, SOI = silicon on insulator, ...). When a geometry came from a template (the 16 examples, `usecases::soi::SoiWire`), the template name is the `family`; otherwise the family is derived from the material stack. The user's title, tags and notes are added on top, and the future AI can add tags and a one-line summary. Re-generated whenever the spec changes.

## 4. Dashboard, tabs and lifecycle

**Top tab strip:** `[Dashboard] [workflow A] [workflow B] ... [+]`. `+`, Ctrl+N and File > New all do the same thing: make a workflow. Ctrl+W closes a tab (the workflow is already saved, closing only hides it from the strip), Ctrl+Tab cycles, Ctrl+Shift+T reopens the last closed. Open tabs are restored on restart.

**Autosave replaces Save.** A workflow is saved on every meaningful change (debounced, atomic). "Save" and "Save as" disappear from the menu; "Export workflow" (zip of everything except heavy artifacts) and "Import" take their place for sharing.

**New creates a draft in memory.** It touches the disk only on the first meaningful change (a geometry edit, a Compute, a title, a note). A draft that was never changed leaves nothing behind. This alone removes the `temp_untitled_*` clutter.

**Dashboard shows, in this order:** *Running/queued* (live status), *Pinned*, *Recent* (touched within the window), then a collapsed *Earlier (N)* with a search box that is always visible. Each card: title, tags, status dot, one-line concept summary, last touched, run count, a sparkline of the main metric if the last experiment has a sweep. Cards are the same size and show nothing else.

| Rule | Default | Why |
|---|---|---|
| "Touched" = edited, ran, or noted. Merely viewing does not count | | looking at an old one must not revive it |
| Not touched for **3 days** and not pinned: moves to *Earlier* | 3 days, setting | your number; automatic, silent |
| Marked **Done**: moves to *Earlier* at once | | explicit completion |
| Dashboard cap: **12** cards; the oldest unpinned falls to *Earlier* | 12 | keeps one screen |
| Nothing is deleted by these rules | | "disappears" means hidden from the dashboard, still found by search |
| Drafts with no change vanish at once; abandoned empty workflows (no run, no edit beyond the title) are dropped after 1 day | | no clutter |
| Failed or cancelled runs are hidden in the table by default and their heavy files are removed after 7 days | 7 days | |
| Heavy artifacts of unstarred runs: keep the newest 5 per experiment, delete older; starred forever; metrics rows forever | 5 | disk stays bounded, history stays complete |
| Explicit **Delete workflow** moves the folder to the OS recycle bin, never a silent erase | | |

## 5. Inside a workflow

**Experiment strip** on the left: experiments as a short list with fork lineage shown as indentation; buttons *New*, *Fork this*, *Link to...*. Selecting one shows its three sections as a segmented control at the top: **Concept | Build | Analysis**.

| Section | What it holds | Reuses (today's tab / code) |
|---|---|---|
| **Concept** | geometry (shapes, model, solids), materials per region, the physics kind and the few essential physical inputs (wavelength or frequency range, symmetry, polarization). The concept card is shown here as a read-only summary strip, so you see what search will see | Geometry tab (`tab_geometry`, `model_editor`, `shapes_editor`, `solids_editor`), Materials tab, the physics choice from `NewDialog` |
| **Build** | mesh rules (element sizes, extrusion layers, PML), boundary conditions, the study (modes, static, eigenfrequency, transient), solver options, parameter table and sweep set-up, the **Compute** button with progress and Cancel, the queue | Physics tab (conditions), Study tab, `Geometry::Model { h, extrude }`, `worker::spawn_solve`, `Progress` |
| **Analysis** | the runs table, plots of metric against parameter, comparison of selected runs, field viewer for one run, analysis functions (section 6.3), the polarization view | Results tab (`tab_results`, `ResultView`, viewport), Analysis tab (`tab_analysis`), `export_csv`, the "Python (optional)" menu |

The split follows the data: Concept is what you would write in a paper's set-up paragraph, Build is how the computer does it, Analysis is what comes out. Anything that changes the spec (Concept, Build) changes the spec hash, and the next Compute makes a new run against the new frozen spec; an old run never changes when you edit afterwards.

## 6. Experiment tracking without Excel (the MLflow take)

**Take from MLflow:** the experiment > run > params/metrics/tags/artifacts model; the runs table with sortable param and metric columns; compare-selected-runs; a tiny filter language (`params.w between 4e-7 and 6e-7 and metrics.neff_re > 2.4`); star/tag/notes on runs. **Leave out:** the server, a SQL backend, Python and autologging, the model registry (for now), heavy dashboards.

### 6.1 Runs table

One table per experiment, fast for 10,000 rows (virtualised rows; `egui_extras::TableBuilder` or a small hand-made grid). Columns: status, time, then every param and metric as a column, choosable and remembered. Click to sort, type to filter with the tiny language, shift-click to select several rows, star, tag, note. *Export CSV* stays for people who want it elsewhere, but is never required to see or compare anything.

### 6.2 Plots and compare

`egui_plot` (already a dependency): metric vs parameter (line/scatter), one series per selected group, error-free and exportable as PNG (the `image` crate with `png` is already in the app, used for `shot`). *Compare* shows selected runs side by side: parameter diff (only the differing keys), overlaid curves of the same metric, and the two field views in the existing viewport. Parallel-coordinates and 2D heat maps over two parameters come later and need nothing new in the data.

### 6.3 Analysis functions

"Run functions on results to plot data" needs a small registry: a function has an id, a description, the inputs it accepts (one run, a set of runs, a column), typed parameters, and an output (scalars, a table, a plot). Built-ins: *metric vs parameter*, *group by*, *min/max/argbest*, *difference of two runs*, *dispersion* (neff vs wavelength, group index), the existing *polarization* view as the first entry. Outputs are saved under the experiment as analysis artifacts, so they are searchable and an AI can call them by id. The "Python (optional)" menu becomes one more kind of function ("external script") instead of a top-level menu.

## 7. Finding "that plasmon scattering on a 300 nm SOS frame"

Three layers, all in memory, no search engine dependency:

1. **Facets** from the cards: physics, study, dimension, materials, family, wavelength range, parameters present.
2. **Text** over title, tags, notes, keywords and the summary: lower-cased tokens with the synonym table (SOS, silicon on sapphire, sapphire, Al2O3 ...).
3. **Numbers with units:** the query parser turns "300 nm" into 3e-7 m and matches any card dimension or run parameter within a tolerance (default 2 %), also inside a sweep's min..max range.

The query box on the dashboard accepts plain words ("plasmon 300 nm sapphire frame") and the filter syntax of 6.1. Scan time for a few thousand manifests is milliseconds; `index.cache` is only a start-up accelerator. Reading all manifests at start-up is parallel (rayon) and reads headers only: a few hundred workflows should open the dashboard in well under 100 ms (to be measured, not promised). `runs.jsonl` is only parsed when an experiment is opened.

Two honest caveats. **"SOS" is ambiguous:** the repo contains both SOI (`usecases::soi`) and silicon-on-sapphire chips (`chips.rs`, the rib on sapphire); the card records the real material stack, so both are found by materials either way, and the synonym table should list SOS as silicon-on-sapphire unless you meant SOI. And **"plasmon scattering" is not a study nyx has today** (physics: electrostatics, heat, solid, EM waves with mode and eigenfrequency studies, Brillouin); it is findable only as a user or AI tag until a driven scattering study exists, and that is a solver feature outside this design.

## 8. Ready for AI, MCP and surrogate training

**One command layer.** Every operation is a function in `nyx-store` (plus the session functions of the MCP plan), not code inside a button handler:

```text
workflow.create / open / close / pin / done          experiment.create / fork / link / set_spec
run.start (blocking or detached) / cancel / star / tag / note
runs.query(filter, columns) / runs.compare(ids) / runs.dataset(experiment, inputs, outputs)
search(query) / analysis.run(id, selection, params) / card.get / manifest.get
```

The window, `nyx --batch`, Rust scripts and an MCP server are four callers of the same functions. This is what makes "just write a prompt" possible: an agent needs only `search`, `runs.query`, `analysis.run` and `run.start`, all of which return small, typed, bounded results (counts and tables, never fields), consistent with the MCP plan's principle that fields stay out of the model's context.

**Merging with the MCP plan.** Its `Workspace` session folder and `session.nyx.xml` become a workflow with `origin: "ai"`; its `session.nyx.xml` is the experiment's `experiment.nyx.xml`; its content-hashed handle cache is the spec hash plus the run table (the same inputs and femcore version means the same run, so `run.start` can answer from `runs.jsonl` instead of solving). Its `journal.jsonl` stays a separate file because it records *operations* (who did what), whereas `runs.jsonl` records *results*. Result: a user opens in the window exactly what the AI built, and the AI sees exactly what the user built.

**Surrogates.** What a prompt like "train a surrogate of neff over width and thickness" needs from this design, all cheap now:
* typed numeric params and metrics with bounds and units, in one query (`runs.dataset` returns inputs/outputs as columns; CSV or NPY first, HDF5 when built with it);
* deterministic hashes, so a dataset is reproducible and a surrogate names exactly which runs trained it;
* run status and starred/outlier flags, so failed or non-converged points are excluded explicitly;
* a design-of-experiments helper (Latin hypercube or Sobol points) that feeds the sweep machinery, and `SweepSession` from the sweep proposal so 200 points are not 200 cold solves;
* a **surrogate experiment kind** (`kind: "surrogate"`): an artifact (weights or GP data), its input bounds, the training run ids, a held-out validation error, and a hard warning on queries outside the training domain. First model in Rust can be small (RBF or Gaussian process on scalar metrics; `linfa`/`candle`/`burn` are options for later); the data model does not care.

## 9. Speed and "things should be removed"

Speed: manifests and run lines are plain serde JSON (no database); appends are O(1); the table is virtualised; the dashboard reads headers only; solves keep their own thread (and later a detached worker), so the UI never waits for I/O. Scheduling: solvers already use all cores, so a **global queue with one running solve** (shown on the dashboard) is the default; per-workflow queues are a setting for small 2D jobs.

To remove or fold away (each frees menu or screen space):

* File > Save, Save as, the `autosave` toggle and folder: replaced by autosave of the workspace. **[needs app edit]**
* The Examples menu (16 entries): become templates in the New-experiment dialog.
* The "Python (optional)" menu: becomes an analysis function.
* Six inspector tabs: three sections.
* `temp_untitled_*` files: gone; runs are records.
* `NewDialog` as a modal: replaced by a small inline template picker.
* Duplicated global state in `App` (one project, one mesh, one results list): moves into a per-workflow view struct.

## 10. Gaps, feature by feature

"Edit" means a change to `app` (`femcore` and `geom` need none for any item here). Sizes: S under a day, M a few days, L a week or more of agent time.

| # | Feature | Exists | Gap | Edit | Size |
|---|---|---|---|---|---|
| F0 | `nyx-store` crate: ids, manifest and run schemas, atomic I/O, JSONL append, repair, unit tests | nothing | new crate, depends on `femcore` (xmlio, problem) and `geom` read-only; reads `.nyx.xml` as `serde_json::Value` and deserialises only `FeatureList` and `ProblemSpec` | no | M |
| F1 | Concept card generator | `FeatureList` params, `Selector`, `ProblemSpec`, `Study` | summariser plus synonym table; family names for templates | no | M |
| F2 | Run recording | `worker::autosave_run`, `SolveResult`, `seconds`, `dofs` | write a run record on every solve (`App::poll` on `Msg::Solved`), typed metrics extraction, frozen spec by hash, replace `temp_untitled_*` | yes | M |
| F3 | Retention and housekeeping | none (autosave grows forever) | rules of section 4, run on open and after each solve | no (logic), yes (hook) | S |
| F4 | Multi-workflow tabs | single `App.proj/doc/results/vp`, `new_project` replaces | per-workflow view struct, tab strip, state restore, `with_context` test harness adapted | yes | L |
| F5 | Job queue across workflows | one `job`, "refuse while busy" | global queue, per-tab status, cancel per job | yes | M |
| F6 | Dashboard and lifecycle rules | File > New dialog | cards, Pinned/Recent/Earlier, rules table, drafts, autosave of workflows | yes | M |
| F7 | Concept / Build / Analysis sections | `Tab` enum (6), `inspector()`, `script.rs` `tab` command | regroup existing editors; no new editors | yes | S-M |
| F8 | Experiments in a workflow, fork and link | none | experiment strip, fork copies spec and records parent, link labels | yes (UI), no (store) | M |
| F9 | Runs table, filter language, star/tag/note | `ResultView` shows one result | table widget (`egui_extras` is not yet a dependency), filter parser in `nyx-store` | yes (UI) | M |
| F10 | Plots and compare | `egui_plot`, viewport, `shot` PNG | metric-vs-parameter plots, run diff, side-by-side fields | yes | M |
| F11 | Parameter table and sweeps | `FeatureList::with` only; proposal in `SWEEP_CACHING_PROPOSAL.md` section 5 and order of work steps 0 to 2 | shared parameter table in the project, Sweep section in Build, one run per point, resumable | yes (and femcore for caching later) | L |
| F12 | Search | none | facets, text with synonyms, numeric-with-units parser, optional `index.cache` | no (store), yes (dashboard box) | M |
| F13 | Analysis function registry | polarization view; Python menu; `export_csv` | registry, built-ins of 6.3, outputs saved as artifacts | yes | M-L |
| F14 | `.nyx` file handling | `open_path` knows xml, json, h5 | associate `.nyx`, open a workflow into a tab, Export/Import zip | yes | S |
| F15 | Headless and batch | `nyx --batch` solves one project, records nothing | `--workflow DIR --experiment ID` records runs; `--sweep` from the sweep proposal | yes | S |
| F16 | Command layer and MCP hookup | proposal only (`nyx-session`, journal) | expose F0-F15 as the same functions the MCP plan calls; workflow with `origin: "ai"` | no (separate crate) | M |
| F17 | Dataset export and surrogate kind | nothing | `runs.dataset`, DOE helper, surrogate artifact with provenance and domain bounds | no | M, later |
| F18 | Detached solve worker | proposal (MCP plan section 20) | `nyx-solve` process; runs append to `runs.jsonl` | no (new binary) | L, later |

### 10.1 Suggested order (value first, riskiest refactor deliberately not first)

1. **F0 + F1 + F3 (store, card, retention):** testable alone with the existing examples; no UI; no edits to read-only code. Proves the schema and the size claims on real `.nyx.xml` files (the 16 examples and `runs/rib`).
2. **F2:** every Compute becomes a recorded run in the current single-window app. The first visible gain and the data every later feature needs.
3. **F9 + F10 (table, plots, compare)** in the current window, for one experiment. This is the MLflow value, delivered before any restructuring.
4. **F4 + F5 + F6 + F7 + F8 + F14 (tabs, queue, dashboard, sections, experiments, `.nyx` files):** the big refactor of `App`, done once the store is proven.
5. **F12 (search)** as soon as there are a few workflows to search.
6. **F11 (parameters and sweeps)**, following the sweep proposal order; **F13 (analysis functions).**
7. **F15, F16, F17, F18** for MCP and surrogates.

Items 1 and a first cut of 2 need no go-ahead on the read-only rule except for the small hook in F2; everything from 2 on needs your say-so to edit `app`.

## 11. Decisions I took (reverse any) and open questions

1. **`.nyx` is one manifest per workflow, not per experiment.** One file answers "which work was this", and a workflow is the unit you will want to find and share. The alternative (one per experiment) makes search read many more files for the same information.
2. **JSON, not XML, for `.nyx` and runs.** Faster, appendable, and trivial for MCP; the experiment spec stays `.nyx.xml` so every existing project and its round-trip tests keep working. XML for the manifest is possible but makes appending awkward.
3. **3 days is the default for leaving the dashboard, "Done" leaves at once, nothing is deleted by rules.** Also a 12-card cap and drafts that never touch disk.
4. **`Documents\nyx\` as the workspace**, with a one-time import of existing `nyx_runs` autosaves as a single "Imported runs" workflow, or left alone if you prefer. Tell me if the workspace should live elsewhere.
5. **"SOS" read as silicon-on-sapphire in the synonym table** (the repo has sapphire chips); SOI is separate. Tell me if you meant SOI.
6. **Open question:** should an experiment's runs be allowed to mix different specs (editing the geometry between runs and keeping all rows in one table) or should a changed spec always start a new experiment? I assumed the former, with the spec hash shown as a column and a "group by spec" toggle, because it matches how people actually iterate. The latter is cleaner but produces many tiny experiments.
7. **Risk:** the `App` refactor (F4) touches about 2,800 lines of `app.rs` and `tabs.rs` and every test that uses `with_context`; keeping the per-workflow view struct as a plain move of existing fields (no behaviour change) is what keeps it safe.
