/*
 * nyx architecture model (Structurizr DSL, C4).
 *
 * Source of truth for the picture; the code is the source of truth for the facts. Every element names the code it
 * stands for (property "code") and the ADRs that justify it (property "adr"). docs/architecture/CODE_TO_ADR.md is the
 * same mapping file by file. Known problems the model shows are listed in docs/architecture/ROUGH_SPOTS.md and are
 * marked here with the tag "rough" on the relationship or element they concern.
 *
 * Written 2026-11-03 against main a047506. View: `structurizr-lite` on this folder, or
 * `structurizr-cli export -workspace workspace.dsl -format mermaid|plantuml`.
 */
workspace "nyx" "Desktop finite-element application for photonics, phononics and stimulated Brillouin scattering, in pure Rust." {

    !identifiers hierarchical
    !adrs adr

    model {

        researcher = person "Researcher" "Designs waveguides, rings and phononic cells; runs mode, Brillouin and eigenfrequency studies; compares runs." "Person"
        aiHost = softwareSystem "AI host" "Claude Code or Claude Desktop (any MCP host). Drives nyx through typed tools." "External"
        viewer = softwareSystem "Results viewer" "ParaView or any XDMF/HDF5 reader. Optional; only with the hdf5 feature." "External"
        wsFolder = softwareSystem "Workspace folder" "Documents/nyx on disk: workflows, runs, .nyx archives, project files." "External,Store"

        nyx = softwareSystem "nyx" "Geometry, meshing, finite-element solvers and result analysis for wave problems in wafer-based devices." {
            properties {
                "adr" "0001, 0002, 0004, 0013"
                "code" "Cargo.toml (workspace: geom, femcore, app, nyx-store), ai/Cargo.toml (second workspace)"
            }

            app = container "app" "The egui window (binary nyx): geometry editor, mesh view, physics tabs, results, analysis, workflow dashboard; also the headless batch runner." "Rust, eframe/egui" "Executable" {
                properties {
                    "adr" "0023, 0014, 0015"
                    "code" "crates/app/src"
                }
                project = component "Project and state" "Project (saved .nyx.json / .nyx.xml), GeometrySpec, ResultView, tab state. Note: Project is defined here, not in a library." "state.rs, examples.rs" "rough"
                worker = component "Solve worker" "Background jobs: mesh, solve, cancel and progress; writes results." "worker.rs, batch.rs" {
                    properties {
                        "adr" "0009, 0011"
                    }
                }
                tabs = component "Tabs and viewport" "Model, Mesh, Physics, Results, Analysis tabs; 2D/3D drawing; colour maps." "tabs.rs, geomui.rs, viewport.rs, colormap.rs"
                shell = component "Workflow shell" "Dashboard, workflow view, sheets, share; talks to nyx-store." "shell.rs, wfview.rs, sheets.rs" {
                    properties {
                        "adr" "0014"
                    }
                }
                script = component "Script commands" "Text commands that drive the same Project." "script.rs"
                connectMenu = component "Connect menu" "Registers nyx-mcp with Claude Code and Claude Desktop." "connect.rs" {
                    properties {
                        "adr" "0015"
                    }
                }
            }

            store = container "nyx-store" "Workspace, workflow, experiment and run storage: manifests, append-only run log, versions, search, .nyx archives." "Rust library" "Library" {
                properties {
                    "adr" "0014, 0019"
                    "code" "crates/nyx-store/src"
                }
                wsComp = component "Workspace" "workspace.json, list of workflows, settings." "workspace.rs"
                flow = component "Workflow" "workflow.json, history.jsonl, notes.md, specs/, one folder per experiment." "flow.rs, model.rs, versions.rs"
                runs = component "Runs" "runs.jsonl records and their digest." "runs.rs"
                archive = component "Archive" ".nyx = deterministic ZIP with a manifest; safe extraction." "archive.rs"
                search = component "Search and card" "Generated concept card, synonym search, number-with-unit search." "card.rs, search.rs"
            }

            geom = container "geom" "Geometry kernel: serialisable feature list, parameters, arrangement into a region complex, selectors, LabelMap, exact queries, layer stacks. No dependency on femcore." "Rust library (serde only)" "Library" {
                properties {
                    "adr" "0002, 0003, 0001"
                    "code" "crates/geom/src, crates/geom/golden"
                }
                featureList = component "Feature list" "FeatureList, parameters, expressions: the single source of truth for a geometry." "model.rs, expr.rs, print.rs, diff.rs"
                eval = component "Evaluator" "Tiers, cache, content signatures; turns a FeatureList into Evaluated." "eval.rs"
                arrangement = component "Arrangement and complex" "Planar arrangement, curves, regions with shared interfaces." "arrange.rs, complex.rs, curve.rs, prim.rs, planar.rs, scene.rs"
                stack = component "Layer stack (2.5D/3D)" "Prismatic models as layer stacks, 3D topology." "stack.rs, topo3.rs"
                selectors = component "Selectors and labels" "Rule language that picks entities; LabelMap compiles names to the integers solvers read." "select.rs, labels.rs"
                sim = component "Simulation data" "PML regions, size rules, periodic pairs as data." "sim.rs"
                queries = component "Queries and tessellation" "Exact measures and closest points; display-only tessellation." "query.rs, tess.rs"
            }

            femcore = container "femcore" "Solver library: meshes, mesh generation and adaptation, sparse and eigen solvers, physics kernels, analysis, project files. Also the home of use cases and the four ring chips." "Rust library" "Library" {
                properties {
                    "adr" "0001, 0004, 0005, 0006, 0007, 0008, 0009, 0010, 0011, 0012, 0013, 0021, 0022"
                    "code" "crates/femcore/src"
                }

                problem = component "Problem and result types" "ProblemSpec, Materials, Bc, Study, SolveResult, Solution, Field; serde. The interface every front-end shares." "problem.rs, material.rs, expr.rs, validate.rs, progress.rs, par.rs" {
                    properties {
                        "adr" "0001, 0019, 0011"
                        "contract" "contracts/problem_spec, contracts/solve_result"
                    }
                }
                solveEntry = component "solve()" "The one entry point: validates, dispatches by physics, refuses NaN, appends the stage log." "solve.rs" {
                    properties {
                        "adr" "0009"
                    }
                }
                section = component "Wafer frame" "The one place that knows axis numbers: plane/local tensors, Voigt order, field names, propagation sign." "section.rs" {
                    properties {
                        "adr" "0001"
                        "contract" "contracts/frames"
                    }
                }
                dfs = component "dfs: diffusion" "Electrostatics and heat transfer on P1/P2 simplices." "dfs.rs" "Kernel"
                sms = component "sms: solid mechanics" "Static, eigenfrequency, kz, PML, low-reflecting, Floquet, axisymmetric." "sms*.rs, periodic.rs" "Kernel"
                ems = component "ems: EM waves" "2D modes (real, complex/PML/tilted crystal), axisymmetric, 3D cavities." "ems*.rs" "Kernel"
                sbs = component "sbs: Brillouin" "Optical modes, photoelastic overlap, phase-matched elastic modes, gain; ring variant." "sbs.rs, sbs_ring.rs" "Kernel"
                analysis = component "Analysis" "Polarization (Stokes, Poincare) and mode classification, in the solvers local axes." "analysis.rs, polarization.rs" {
                    properties {
                        "adr" "0021"
                    }
                }
                useCases = component "Use cases and chips" "Ready-made SOI, Bragg and ring-chip problems. Application content living in the solver library." "usecases/, chips.rs" "rough" {
                    properties {
                        "adr" "0022"
                    }
                }
                meshTypes = component "Mesh types" "SimplexMesh (dim at run time), nyx CellMesh (named element types), numbering of unknowns, draw-only view." "mesh.rs, cellmesh.rs, view.rs, fields.rs" {
                    properties {
                        "adr" "0004, 0005"
                        "contract" "contracts/mesh"
                    }
                }
                geomdrive = component "Geom to mesh" "Evaluates a geom model, meshes it, names domains and boundaries through LabelMap, builds PML materials." "geomdrive.rs, mesher.rs" {
                    properties {
                        "adr" "0002, 0003"
                    }
                }
                meshing = component "Meshing interface" "MeshPlan: Simple or Adaptive mesher, structured PML, mesh cache, physics sizing, estimators and the solve-estimate-adapt loop." "meshing.rs, meshing_pml.rs, physics_size.rs, estimate.rs, estimate3.rs" {
                    properties {
                        "adr" "0006, 0007, 0017"
                    }
                }
                mesher3d = component "3D legacy meshers" "Isosurface-stuffing tets from primitives, Rivara refinement, structured 3D PML." "mesher3d.rs, refine3d.rs, pml3.rs" "rough"
                meshersPort = component "meshers (ported)" "Own mesh types, metric adaptation 2D and 3D, structured PML generator, .fmsh format, validation harness; interop converts at the boundary." "meshers/" {
                    properties {
                        "adr" "0005, 0007, 0016"
                    }
                }
                solversPort = component "solvers (ported)" "Krylov, eigen, preconditioners, time stepping, optimisation, for real and complex systems. Used by examples only; not wired into the physics." "solvers/" "rough" {
                    properties {
                        "adr" "0008"
                    }
                }
                linalg = component "Linear algebra (original)" "Sparse matrices, faer wrappers, AMG, shift-invert Arnoldi, LOBPCG, element basics. The code the physics really calls." "sparse.rs, csparse.rs, linalg.rs, amg.rs, eigen.rs, lobpcg.rs, fem.rs" "rough"
                files = component "Files" "XML encoding of any serde value; HDF5/XDMF writer behind a feature." "xmlio.rs, h5io.rs" {
                    properties {
                        "adr" "0013, 0019"
                    }
                }
            }

            session = container "nyx-session" "Typed tools an AI host can call: build geometry, mesh, define a problem, solve, select modes, check. Content-hashed handles, journal, assumption ledger." "Rust library (separate workspace ai/)" "Library" {
                properties {
                    "adr" "0015, 0019"
                    "code" "ai/nyx-session/src"
                }
                tools = component "Tools" "geom_template, mesh_make, problem_define, solve_run, result_modes, check_run, recipe_run, session_status." "tools.rs"
                params = component "Parameter tables" "One table per tool: JSON Schema, validation, defaults, assumption ledger." "params.rs" {
                    properties {
                        "contract" "contracts/mcp_tools"
                    }
                }
                handles = component "Session and handles" "g_, m_, p_, r_ handles by FNV-1a content hash; journal.jsonl; geom JSON per geometry." "workspace.rs"
                guide = component "Guide cards" "Short texts served by get_guide." "guide.rs"
            }

            mcp = container "nyx-mcp" "Model Context Protocol server over stdio (newline-delimited JSON-RPC 2.0), hand-written, no async runtime." "Rust binary" "Executable" {
                properties {
                    "adr" "0015"
                    "code" "ai/nyx-mcp/src/main.rs"
                }
            }
            connect = container "nyx-connect" "Finds Claude Code and Claude Desktop installs and registers nyx-mcp with them; command-line twin of the Connect menu." "Rust library and binary" "Executable" {
                properties {
                    "adr" "0015"
                    "code" "ai/nyx-connect/src"
                }
            }

            projectFile = container "Project file" "Everything saved for one model: geometry, materials, boundary conditions, studies. JSON or XML." "File .nyx.json / .nyx.xml" "Store" {
                properties {
                    "adr" "0019"
                    "contract" "contracts/problem_spec (the ProblemSpec part); Project itself has none yet"
                }
            }
            workspaceStore = container "Workspace data" "workspace.json, workflow.json, history.jsonl, runs.jsonl, specs/, .nyx archives." "Files and folders" "Store" {
                properties {
                    "adr" "0014"
                    "contract" "contracts/store"
                }
            }
            sessionStore = container "Session folder" "journal.jsonl, g_*.geom.json, result digests of an AI session." "Files" "Store" {
                properties {
                    "adr" "0015"
                }
            }
        }

        # context
        researcher -> nyx.app "Edits models, runs studies, reads results"
        researcher -> aiHost "Asks in words"
        aiHost -> nyx.mcp "Calls tools" "MCP over stdio (JSON-RPC 2.0)"
        nyx.app -> viewer "Writes HDF5/XDMF results" "optional hdf5 feature"
        nyx.app -> wsFolder "Reads and writes"

        # containers
        nyx.app -> nyx.femcore "Meshes, solves, analyses" "Rust API"
        nyx.app -> nyx.geom "Edits and evaluates feature lists" "Rust API"
        nyx.app -> nyx.store "Saves workflows, runs, archives" "Rust API"
        nyx.app -> nyx.projectFile "Saves and loads" "serde JSON / XML"
        nyx.store -> nyx.femcore "Uses ProblemSpec, SolveResult and xmlio" "Rust API" "rough"
        nyx.store -> nyx.workspaceStore "Reads and writes"
        nyx.femcore -> nyx.geom "Evaluates models, reads selectors and labels" "Rust API"
        nyx.session -> nyx.femcore "Meshes, solves, analyses" "Rust API"
        nyx.session -> nyx.sessionStore "Journals every call"
        nyx.mcp -> nyx.session "Dispatches tools" "Rust API"
        nyx.connect -> nyx.mcp "Registers and probes" "spawns the binary"
        nyx.app -> nyx.connect "Connect menu duplicates this logic in connect.rs" "" "rough"

        # components of femcore
        nyx.femcore.solveEntry -> nyx.femcore.problem "Reads ProblemSpec, returns SolveResult"
        nyx.femcore.solveEntry -> nyx.femcore.dfs "Dispatch (Electrostatics, Heat)"
        nyx.femcore.solveEntry -> nyx.femcore.sms "Dispatch (Solid)"
        nyx.femcore.solveEntry -> nyx.femcore.ems "Dispatch (Emw 2D / 3D)"
        nyx.femcore.solveEntry -> nyx.femcore.sbs "Dispatch (Brillouin)"
        nyx.femcore.solveEntry -> nyx.femcore.meshTypes "Takes SimplexMesh"
        nyx.femcore.dfs -> nyx.femcore.section "Wafer names and tensors"
        nyx.femcore.sms -> nyx.femcore.section "Wafer names, Voigt order"
        nyx.femcore.ems -> nyx.femcore.section "Wafer names, propagation sign"
        nyx.femcore.sbs -> nyx.femcore.section "Wafer names"
        nyx.femcore.sbs -> nyx.femcore.ems "Optical modes"
        nyx.femcore.sbs -> nyx.femcore.sms "Phase-matched elastic modes (kz)"
        nyx.femcore.dfs -> nyx.femcore.linalg "Assembles and factorises" "" "rough"
        nyx.femcore.sms -> nyx.femcore.linalg "Assembles, eigensolves" "" "rough"
        nyx.femcore.ems -> nyx.femcore.linalg "Assembles, eigensolves" "" "rough"
        nyx.femcore.analysis -> nyx.femcore.problem "Reads Solution fields"
        nyx.femcore.useCases -> nyx.femcore.geomdrive "Builds problems from geom models"
        nyx.femcore.useCases -> nyx.femcore.problem "Fills ProblemSpec"
        nyx.femcore.geomdrive -> nyx.geom "Evaluates models"
        nyx.femcore.geomdrive -> nyx.femcore.meshTypes "Returns MeshedGeom with SimplexMesh"
        nyx.femcore.meshing -> nyx.femcore.geomdrive "Plans the mesh of a model"
        nyx.femcore.meshing -> nyx.femcore.meshersPort "Adapts through interop" "Rust API"
        nyx.femcore.meshing -> nyx.femcore.problem "Physics sizing reads ProblemSpec"
        nyx.femcore.meshersPort -> nyx.femcore.meshTypes "interop converts both ways" "" "rough"
        nyx.femcore.mesher3d -> nyx.femcore.meshTypes "Builds SimplexMesh"
        nyx.femcore.mesher3d -> nyx.femcore.meshersPort "pml3 and physics_size::adapt_tets overlap with meshers::adapt::pml" "" "rough"
        nyx.femcore.solversPort -> nyx.femcore.linalg "Intended replacement; no caller yet" "" "rough"
        nyx.femcore.files -> nyx.femcore.problem "Serialises specs"

        # components of the app
        nyx.app.tabs -> nyx.app.project "Edits"
        nyx.app.worker -> nyx.app.project "Reads a frozen spec"
        nyx.app.worker -> nyx.femcore.solveEntry "Runs"
        nyx.app -> nyx.femcore.solveEntry "Runs (worker thread)"
        nyx.app.worker -> nyx.femcore.meshing "Meshes"
        nyx.app.shell -> nyx.store "Saves runs"
        nyx.app.script -> nyx.app.project "Edits"
        nyx.app.connectMenu -> nyx.mcp "Registers"

        # components of the store
        nyx.store.flow -> nyx.store.wsComp "Lives in"
        nyx.store.runs -> nyx.store.flow "Appends to"
        nyx.store.archive -> nyx.store.flow "Packs and unpacks"
        nyx.store.search -> nyx.store.flow "Reads manifests"
        nyx.store.search -> nyx.femcore "Card is generated from a ProblemSpec" "xmlio, problem" "rough"

        # components of geom
        nyx.geom.eval -> nyx.geom.featureList "Evaluates"
        nyx.geom.eval -> nyx.geom.arrangement "2D"
        nyx.geom.eval -> nyx.geom.stack "Prismatic 3D"
        nyx.geom.selectors -> nyx.geom.eval "Resolves against Evaluated"
        nyx.geom.sim -> nyx.geom.selectors "Uses"
        nyx.geom.queries -> nyx.geom.eval "Reads Evaluated"

        # components of the session
        nyx.session.tools -> nyx.session.params "Validates"
        nyx.session.tools -> nyx.session.handles "Records"
        nyx.session.tools -> nyx.femcore "Calls" "Rust API"
        nyx.session.tools -> nyx.session.guide "Serves"
        nyx.mcp -> nyx.session.tools "tools/call"

        # data flow a researcher sees
        nyx.femcore.useCases -> nyx.femcore.meshing "mesh_model, mesh_evaluated"
    }

    views {

        systemContext nyx "Context" "Who and what nyx talks to." {
            include *
            autolayout lr
        }

        container nyx "Containers" "Crates, binaries and files of nyx. Dashed red = rough spot (see ROUGH_SPOTS.md)." {
            include *
            autolayout lr
        }

        component nyx.femcore "femcore" "Inside the solver library." {
            include *
            autolayout tb
        }

        component nyx.geom "geom" "Inside the geometry kernel." {
            include *
            autolayout lr
        }

        component nyx.app "app" "Inside the egui application." {
            include *
            autolayout lr
        }

        component nyx.store "store" "Inside the workflow store." {
            include *
            autolayout lr
        }

        component nyx.session "session" "Inside the AI session layer." {
            include *
            autolayout lr
        }

        dynamic nyx.femcore "SolveACrossSection" "What happens on solve(mesh, spec, state)." {
            nyx.app -> nyx.femcore.solveEntry "solve_with(mesh, spec, state, progress)"
            nyx.femcore.solveEntry -> nyx.femcore.problem "validate::check_inputs, order is 1 or 2"
            nyx.femcore.solveEntry -> nyx.femcore.ems "dispatch (Emw 2D)"
            nyx.femcore.ems -> nyx.femcore.section "wafer axes and names"
            nyx.femcore.ems -> nyx.femcore.linalg "assemble, shift-invert eigensolve"
            nyx.femcore.solveEntry -> nyx.femcore.problem "check_result, append stage log"
            autolayout lr
        }

        dynamic nyx "AIRecipe" "An AI host asks for the effective index of a silicon wire." {
            aiHost -> nyx.mcp "recipe_run"
            nyx.mcp -> nyx.session "dispatch"
            nyx.session -> nyx.femcore "geom -> mesh -> ProblemSpec -> solve"
            nyx.session -> nyx.sessionStore "journal the call"
            autolayout lr
        }

        styles {
            element "Person" {
                shape Person
                background #1d4e89
                color #ffffff
            }
            element "External" {
                background #8a8a8a
                color #ffffff
            }
            element "Software System" {
                background #1168bd
                color #ffffff
            }
            element "Container" {
                background #438dd5
                color #ffffff
            }
            element "Component" {
                background #85bbf0
                color #000000
            }
            element "Executable" {
                shape RoundedBox
            }
            element "Library" {
                shape Component
            }
            element "Store" {
                shape Cylinder
            }
            element "Kernel" {
                background #2e7d32
                color #ffffff
            }
            element "rough" {
                stroke #b02a37
                strokeWidth 5
            }
            relationship "Relationship" {
                routing Orthogonal
            }
            relationship "rough" {
                color #b02a37
                dashed true
            }
        }
    }
}
