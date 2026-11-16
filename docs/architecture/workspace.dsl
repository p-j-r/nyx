/*
 * nyx architecture model (Structurizr DSL, C4): context, containers (one per crate), components (every container) and two
 * dynamic views. Code level (C4 level 4) is not a Structurizr view: it is a class-diagram page per component, in c4/4-code/.
 * As of the end of the restructure (phase 9, 2026-11-13).
 *
 * Hand-written here: people, systems, containers, the container relationships (exactly the layering table of
 * crates/nyx-engine/tests/layering.rs) and the dynamic views. Generated, by c4/tools/gen.pl from the Rust source and
 * c4/spec/components.spec: the components of every container (c4/model/<container>.components.dsl), the relationships between
 * components and from components to containers (c4/model/relationships.dsl: every `use` and path reference of the source) and
 * the component views (c4/model/views.dsl). Do not edit the generated files; change the spec or the source and run
 * `bash c4/tools/build.sh`.
 *
 * Every element names the code it stands for (property "code") and the ADRs that justify it (property "adr").
 * docs/architecture/CODE_TO_ADR.md is the same mapping file by file.
 *
 * View: `structurizr-lite` on this folder, `structurizr-cli export -workspace workspace.dsl -format mermaid|plantuml`, or the
 * pages and the browser viewer in c4/ (see c4/README.md).
 */
workspace "nyx" "Desktop finite-element application for photonics and phononics, in pure Rust." {

    !identifiers hierarchical
    !adrs adr

    model {

        researcher = person "Researcher" "Designs waveguides, rings and phononic cells; runs mode and eigenfrequency studies; compares runs." "Person"
        aiHost = softwareSystem "AI host" "Claude Code or Claude Desktop (any MCP host). Drives nyx through the engine commands as MCP tools." "External"
        viewer = softwareSystem "Results viewer" "ParaView or VisIt. Reads the XDMF description and the HDF5 results file." "External"
        wsFolder = softwareSystem "Workspace folder" "Documents/nyx on disk: workflows, runs, .nyx archives, project files." "External,Store"

        nyx = softwareSystem "nyx" "Geometry, meshing, finite-element solvers and result analysis for wave problems in wafer-based devices." {
            properties {
                "adr" "0001, 0002, 0013, 0026"
                "code" "Cargo.toml (one workspace: crates/*)"
            }

            app = container "app" "The egui window (binary nyx): geometry editor, mesh view, Study, Physics and Materials tabs drawn from the physics descriptors, results, analysis, workflow dashboard. A front-end of the engine: it holds an Engine and calls commands." "Rust, eframe/egui" "Executable" {
                properties {
                    "adr" "0057, 0014"
                    "code" "crates/app/src"
                }
                !include c4/model/app.components.dsl
            }

            mcp = container "nyx-mcp" "The MCP server (stdio, JSON-RPC 2.0): every tool is one engine command, with the JSON Schema the engine generates." "Rust binary" "Executable" {
                properties {
                    "adr" "0057, 0015"
                    "code" "crates/nyx-mcp/src"
                    "contract" "contracts/mcp_tools.v2.json"
                }
                !include c4/model/mcp.components.dsl
            }

            connect = container "nyx-connect" "Finds Claude Code and Claude Desktop installs and registers nyx-mcp with them; the Connect menu and its command-line twin." "Rust library and binary" "Executable" {
                properties {
                    "adr" "0028"
                    "code" "crates/nyx-connect/src"
                    "contract" "contracts/connect_report.v1.json"
                }
                !include c4/model/connect.components.dsl
            }

            engine = container "nyx-engine" "The Engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them, and the binding to a workflow. Every action is a named command with a JSON Schema." "Rust library" "Library" {
                properties {
                    "adr" "0057, 0058, 0060, 0061, 0062, 0063"
                    "code" "crates/nyx-engine/src"
                    "contract" "contracts/engine_commands.v1.json, engine_project.v2.json, geom_operations.v1.json"
                }
                !include c4/model/engine.components.dsl
            }

            physics = container "nyx-physics" "The Physics trait, the problem envelope (ProblemSpec), descriptors, SolveResult, the registry type, validation, JSON Schemas." "Rust library" "Library" {
                properties {
                    "adr" "0025, 0063"
                    "code" "crates/nyx-physics/src"
                    "contract" "contracts/problem_spec.v3.json, solve_result.v2.json, physics_descriptors.v1.json"
                }
                !include c4/model/physics.components.dsl
            }
            dfs = container "nyx-dfs" "Diffusion physics: heat and electrostatics. Typed problem, three-stage kernels, descriptors, goldens." "Rust library" "Library,Kernel" {
                properties {
                    "adr" "0041, 0009"
                    "code" "crates/nyx-dfs/src"
                    "contract" "contracts/physics_dfs.v1.json"
                }
                !include c4/model/dfs.components.dsl
            }
            sms = container "nyx-sms" "Solid mechanics: static, eigenfrequency, PML, Floquet, low-reflecting boundaries, axisymmetric eigenfrequency study (bodies of revolution, axis conditions, cylindrical layers)." "Rust library" "Library,Kernel" {
                properties {
                    "adr" "0042, 0009"
                    "code" "crates/nyx-sms/src"
                    "contract" "contracts/physics_sms.v1.json"
                }
                !include c4/model/sms.components.dsl
            }
            ems = container "nyx-ems" "EM waves: waveguide modes (real and complex), axisymmetric, cavity in 3D." "Rust library" "Library,Kernel" {
                properties {
                    "adr" "0045, 0009"
                    "code" "crates/nyx-ems/src"
                    "contract" "contracts/physics_ems.v1.json"
                }
                !include c4/model/ems.components.dsl
            }
            analysis = container "nyx-analysis" "Analyses of results: the polarization (Stokes parameters, Poincare sphere, classes, degenerate groups) as a function of a SolveResult and its mesh." "Rust library" "Library" {
                properties {
                    "adr" "0056, 0021"
                    "code" "crates/nyx-analysis/src"
                }
                !include c4/model/analysis.components.dsl
            }
            store = container "nyx-store" "Workspace, workflow, experiment and run storage: manifests, append-only run log, versions, search, .nyx archives. Knows nothing of a project." "Rust library" "Library" {
                properties {
                    "adr" "0014, 0039, 0059, 0064"
                    "code" "crates/nyx-store/src"
                    "contract" "contracts/store_workflow.v1.json, store_run.v1.json"
                }
                !include c4/model/store.components.dsl
            }

            femcore = container "femcore" "The FEM core: meshes, meshing and adaptation, linear algebra behind one driver, finite-element basics, the material specification and library, expressions." "Rust library" "Library" {
                properties {
                    "adr" "0032, 0033, 0029, 0006"
                    "code" "crates/femcore/src"
                }
                !include c4/model/femcore.components.dsl
            }
            geom = container "geom" "Lean CSG-over-B-rep geometry kernel: feature lists, named regions and boundaries, selectors, exact projection oracles. Depends on serde and optional rayon only." "Rust library" "Library" {
                properties {
                    "adr" "0002, 0003, 0063"
                    "code" "crates/geom/src"
                    "contract" "crates/geom/golden (feature lists and expected regions)"
                }
                !include c4/model/geom.components.dsl
            }

            projectFile = container "Project file" "Everything saved for one model: geometry, materials, boundary conditions, studies. JSON with a format key." "File .nyx.json" "Store" {
                properties {
                    "adr" "0063"
                    "contract" "contracts/engine_project.v2.json"
                }
            }
            resultsFile = container "Results file" "Project, mesh and every solution of every solved physics, compressed, exact; with an XDMF description next to it." "File .h5 and .xdmf" "Store" {
                properties {
                    "adr" "0062"
                    "contract" "contracts/VERSIONS.md (results_file, layout 3)"
                }
            }
            workspaceStore = container "Workspace data" "workspace.json, workflow.json, history.jsonl, runs.jsonl, specs/ (experiment.nyx.json, specs/<hash>.nyx.json), run folders with results.h5, .nyx archives." "Files and folders" "Store" {
                properties {
                    "adr" "0014, 0064"
                    "contract" "contracts/store_workflow.v1.json, store_run.v1.json"
                }
            }
            materialLibrary = container "Material library" "materials/*.toml next to the program: silicon, silica, sapphire, metals, nitrides, niobate and more, with sources." "Files" "Store" {
                properties {
                    "adr" "0029, 0034"
                }
            }
        }

        # context
        researcher -> nyx.app "Edits models, runs studies, reads results"
        researcher -> aiHost "Asks in words"
        aiHost -> nyx.mcp "Calls tools" "MCP over stdio (JSON-RPC 2.0)"
        viewer -> nyx.resultsFile "Reads" "XDMF / HDF5"
        nyx.app -> wsFolder "Reads and writes"

        # containers: exactly the layering table of crates/nyx-engine/tests/layering.rs
        nyx.app -> nyx.engine "Calls commands, polls jobs, draws results" "Rust API"
        nyx.app -> nyx.connect "Connect menu" "Rust API"
        nyx.mcp -> nyx.engine "One tool per command" "Rust API"
        nyx.engine -> nyx.physics "Registry, ProblemSpec, SolveResult" "Rust API"
        nyx.engine -> nyx.dfs "Compiled-in physics (feature dfs)" "Rust API"
        nyx.engine -> nyx.sms "Compiled-in physics (feature sms)" "Rust API"
        nyx.engine -> nyx.ems "Compiled-in physics (feature ems)" "Rust API"
        nyx.engine -> nyx.analysis "Analyses of results" "Rust API"
        nyx.engine -> nyx.store "Workflows and runs; passes spec text and a card" "Rust API"
        nyx.engine -> nyx.femcore "Meshing, materials, plots" "Rust API"
        nyx.engine -> nyx.geom "Feature lists, labels, selectors" "Rust API"
        nyx.dfs -> nyx.physics "Implements Physics" "Rust API"
        nyx.sms -> nyx.physics "Implements Physics" "Rust API"
        nyx.ems -> nyx.physics "Implements Physics" "Rust API"
        nyx.dfs -> nyx.femcore "Mesh, assembly, linear algebra" "Rust API"
        nyx.sms -> nyx.femcore "Mesh, assembly, linear algebra" "Rust API"
        nyx.ems -> nyx.femcore "Mesh, assembly, linear algebra" "Rust API"
        nyx.analysis -> nyx.physics "Reads SolveResult" "Rust API"
        nyx.analysis -> nyx.femcore "Mesh geometry, quadrature" "Rust API"
        nyx.physics -> nyx.femcore "Mesh, materials" "Rust API"
        nyx.femcore -> nyx.geom "Evaluates models, reads selectors and labels" "Rust API"

        # containers and data
        nyx.engine -> nyx.projectFile "Saves and loads" "serde JSON"
        nyx.engine -> nyx.resultsFile "Writes and opens" "HDF5 (serial or MPI-IO), XDMF"
        nyx.engine -> nyx.materialLibrary "Loads" "TOML"
        nyx.store -> nyx.workspaceStore "Reads and writes"

        # components: read from the use and path references of the Rust source by c4/tools/gen.pl (from c4/spec/components.spec)
        !include c4/model/relationships.dsl
    }

    views {

        systemContext nyx "Context" "Who uses nyx and what it talks to." {
            include *
            autolayout lr
        }

        container nyx "Containers" "The crates of the workspace and the files they write." {
            include *
            autolayout lr
        }

        # one component view per container, generated
        !include c4/model/views.dsl

        dynamic nyx "SolveFromTheWindow" "What happens when the window solves the active physics." {
            nyx.app -> nyx.engine "start_solve(targets)"
            nyx.engine -> nyx.femcore "mesh the model when none is current"
            nyx.engine -> nyx.ems "Physics::solve(mesh, spec, ctx)"
            nyx.ems -> nyx.femcore "assemble, shift-invert eigensolve"
            nyx.engine -> nyx.store "record the run, freeze the spec"
            nyx.engine -> nyx.resultsFile "write results.h5 and results.xdmf"
            nyx.app -> nyx.engine "poll: take the mesh and the results"
            autolayout lr
        }

        dynamic nyx "AICommand" "An AI host builds a model and solves it through the same commands (an instance: the anchor-loss run of a disk on a post, set_study axisymmetric_eigenfrequency then add_pml_conditions, solve_run and result_info, which returns Q)." {
            aiHost -> nyx.mcp "tools/call add_rect, add_overlay, add_physics, solve_run"
            nyx.mcp -> nyx.engine "call(command, arguments)"
            nyx.engine -> nyx.geom "feature list edit and evaluation"
            nyx.engine -> nyx.femcore "mesh"
            nyx.engine -> nyx.ems "solve"
            nyx.engine -> nyx.store "auto_store: the workflow shows in the dashboard"
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
            relationship "Relationship" {
                routing Orthogonal
            }
        }
    }
}
