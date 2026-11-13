/*
 * nyx architecture model (Structurizr DSL, C4), as of the end of the restructure (phase 9, 2026-11-13).
 *
 * Source of truth for the picture; the code is the source of truth for the facts. Every element names the code it
 * stands for (property "code") and the ADRs that justify it (property "adr"). docs/architecture/CODE_TO_ADR.md is the
 * same mapping file by file. The crate layering is enforced by crates/nyx-engine/tests/layering.rs: the relationships
 * between containers below are exactly its table.
 *
 * View: `structurizr-lite` on this folder, or
 * `structurizr-cli export -workspace workspace.dsl -format mermaid|plantuml`.
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
            }

            mcp = container "nyx-mcp" "The MCP server (stdio, JSON-RPC 2.0): every tool is one engine command, with the JSON Schema the engine generates." "Rust binary" "Executable" {
                properties {
                    "adr" "0057, 0015"
                    "code" "crates/nyx-mcp/src"
                    "contract" "contracts/mcp_tools.v2.json"
                }
            }

            connect = container "nyx-connect" "Finds Claude Code and Claude Desktop installs and registers nyx-mcp with them; the Connect menu and its command-line twin." "Rust library and binary" "Executable" {
                properties {
                    "adr" "0028"
                    "code" "crates/nyx-connect/src"
                    "contract" "contracts/connect_report.v1.json"
                }
            }

            engine = container "nyx-engine" "The Engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them, and the binding to a workflow. Every action is a named command with a JSON Schema." "Rust library" "Library" {
                properties {
                    "adr" "0057, 0058, 0060, 0061, 0062, 0063"
                    "code" "crates/nyx-engine/src"
                    "contract" "contracts/engine_commands.v1.json, engine_project.v2.json, geom_operations.v1.json"
                }
                project = component "Project" "The saved state: feature list, mesh settings, one ProblemSpec per physics, project materials. The format key (1) is checked at the file boundary." "project.rs" {
                    properties {
                        "adr" "0063"
                    }
                }
                commands = component "Commands" "The catalog (name, access, JSON Schema) and the dispatcher; the window and the MCP server use the same entries." "commands.rs, engine.rs, geometry.rs, geomdesc.rs" {
                    properties {
                        "adr" "0057, 0038"
                    }
                }
                jobs = component "Job runner" "Background mesh and solve jobs with progress and cancel; results are applied on poll." "jobs.rs, solve.rs, convergence.rs" {
                    properties {
                        "adr" "0061"
                    }
                }
                meshing = component "Meshing seam" "MeshRequest to a MeshedModel (mesh, names, label resolver); the only caller of the femcore mesher." "meshing.rs, pml.rs" {
                    properties {
                        "adr" "0054, 0050"
                    }
                }
                materials = component "Materials" "Named, oriented and project materials resolved from the library and the project at the study wavelength." "materials.rs" {
                    properties {
                        "adr" "0058"
                    }
                }
                analyse = component "Analyses" "Runs the analyses of nyx-analysis on the latest result of a physics." "analyse.rs" {
                    properties {
                        "adr" "0056"
                    }
                }
                resultsIO = component "Results files" "HDF5 results (serial or MPI-IO), the XDMF description, and opening a results file as a project with its results." "h5io.rs, xdmf.rs, results_file.rs" {
                    properties {
                        "adr" "0062, 0063"
                    }
                }
                storeLink = component "Store link" "Binds the project to a workflow of the workspace, commits versions and records runs (card and metrics made here)." "store_link.rs" {
                    properties {
                        "adr" "0039, 0059"
                    }
                }
                registry = component "Registry" "The physics compiled in (cargo features dfs, sms, ems)." "lib.rs" {
                    properties {
                        "adr" "0025"
                    }
                }
            }

            physics = container "nyx-physics" "The Physics trait, the problem envelope (ProblemSpec), descriptors, SolveResult, the registry type, validation, JSON Schemas." "Rust library" "Library" {
                properties {
                    "adr" "0025, 0063"
                    "code" "crates/nyx-physics/src"
                    "contract" "contracts/problem_spec.v3.json, solve_result.v2.json, physics_descriptors.v1.json"
                }
            }
            dfs = container "nyx-dfs" "Diffusion physics: heat and electrostatics. Typed problem, three-stage kernels, descriptors, goldens." "Rust library" "Library,Kernel" {
                properties {
                    "adr" "0041, 0009"
                    "code" "crates/nyx-dfs/src"
                    "contract" "contracts/physics_dfs.v1.json"
                }
            }
            sms = container "nyx-sms" "Solid mechanics: static, eigenfrequency, PML, Floquet, low-reflecting boundaries, axisymmetric." "Rust library" "Library,Kernel" {
                properties {
                    "adr" "0042, 0009"
                    "code" "crates/nyx-sms/src"
                    "contract" "contracts/physics_sms.v1.json"
                }
            }
            ems = container "nyx-ems" "EM waves: waveguide modes (real and complex), axisymmetric, cavity in 3D." "Rust library" "Library,Kernel" {
                properties {
                    "adr" "0045, 0009"
                    "code" "crates/nyx-ems/src"
                    "contract" "contracts/physics_ems.v1.json"
                }
            }
            analysis = container "nyx-analysis" "Analyses of results: the polarization (Stokes parameters, Poincare sphere, classes, degenerate groups) as a function of a SolveResult and its mesh." "Rust library" "Library" {
                properties {
                    "adr" "0056, 0021"
                    "code" "crates/nyx-analysis/src"
                }
            }
            store = container "nyx-store" "Workspace, workflow, experiment and run storage: manifests, append-only run log, versions, search, .nyx archives. Knows nothing of a project." "Rust library" "Library" {
                properties {
                    "adr" "0014, 0039, 0059, 0064"
                    "code" "crates/nyx-store/src"
                    "contract" "contracts/store_workflow.v1.json, store_run.v1.json"
                }
            }

            femcore = container "femcore" "The FEM core: meshes, meshing and adaptation, linear algebra behind one driver, finite-element basics, the material specification and library, expressions." "Rust library" "Library" {
                properties {
                    "adr" "0032, 0033, 0029, 0006"
                    "code" "crates/femcore/src"
                }
                mesh = component "Mesh" "The native mesh (cells of one kind each), the exchange mesh Mesh the physics take, and the discretisation (edges, P2 numbering)." "mesh/" {
                    properties {
                        "adr" "0005, 0032"
                        "contract" "contracts/mesh_cells.v1.json"
                    }
                }
                meshing = component "Meshing" "MeshRequest to mesh: the lattice mesher, structured absorbing layers, wavelength sizing, drive by geom." "meshing/" {
                    properties {
                        "adr" "0006, 0017, 0052, 0054, 0055"
                        "contract" "contracts/mesh_request.v1.json, mesh_report.v1.json"
                    }
                }
                adapt = component "Adaptation" "Metric adaptation, residual error estimators, the solve-estimate-adapt loop in 2D and 3D." "adapt/" {
                    properties {
                        "adr" "0007, 0053"
                    }
                }
                linalg = component "Linear algebra" "Csr, direct and iterative solvers, shift-invert Krylov-Schur, LOBPCG, AMG, behind one driver." "linalg/" {
                    properties {
                        "adr" "0033, 0011"
                        "contract" "contracts/linalg_options.v1.json"
                    }
                }
                fem = component "Finite-element basics" "Shape functions and quadrature, edge elements, PML coefficients, periodic maps, material properties per domain." "fem.rs, edge.rs, pml.rs, periodic.rs, props.rs" {
                    properties {
                        "adr" "0009"
                    }
                }
                material = component "Material specification" "TOML files with explicit units, tensors, dispersion models evaluated at the study wavelength, the library loader." "material/, optical.rs" {
                    properties {
                        "adr" "0029, 0030, 0034, 0035"
                        "contract" "contracts/material_schema.v1.json"
                    }
                }
                section = component "Frame and sections" "The one place that knows the wafer axes (x lateral, z normal, waves along +y)." "section.rs" {
                    properties {
                        "adr" "0001"
                    }
                }
                expr = component "Expressions" "The safe expression grammar for initial fields and boundary values; a copy of geom's, kept equal by a corpus test." "expr.rs" {
                    properties {
                        "adr" "0065"
                    }
                }
            }
            geom = container "geom" "Lean CSG-over-B-rep geometry kernel: feature lists, named regions and boundaries, selectors, exact projection oracles. Depends on serde and optional rayon only." "Rust library" "Library" {
                properties {
                    "adr" "0002, 0003, 0063"
                    "code" "crates/geom/src"
                    "contract" "crates/geom/golden (feature lists and expected regions)"
                }
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

        # components of the engine
        nyx.engine.commands -> nyx.engine.project "Edits"
        nyx.engine.commands -> nyx.engine.jobs "Starts mesh and solve jobs"
        nyx.engine.commands -> nyx.engine.analyse "analyse"
        nyx.engine.commands -> nyx.engine.resultsIO "export_results, open_results"
        nyx.engine.commands -> nyx.engine.storeLink "store_save, store_open, runs"
        nyx.engine.jobs -> nyx.engine.meshing "Builds the mesh when none is current"
        nyx.engine.jobs -> nyx.engine.materials "Resolves materials at the study wavelength"
        nyx.engine.jobs -> nyx.engine.registry "Finds the physics of a problem"
        nyx.engine.jobs -> nyx.engine.storeLink "Records a run"
        nyx.engine.storeLink -> nyx.engine.resultsIO "Writes results.h5 in the run folder"
        nyx.engine.meshing -> nyx.femcore.meshing "MeshRequest"
        nyx.engine.meshing -> nyx.geom "Evaluates the feature list"
        nyx.engine.materials -> nyx.femcore.material "Library and spec"
        nyx.engine.registry -> nyx.dfs "Registers"
        nyx.engine.registry -> nyx.sms "Registers"
        nyx.engine.registry -> nyx.ems "Registers"
        nyx.engine.analyse -> nyx.analysis "Calls"
        nyx.engine.storeLink -> nyx.store "Commits specs, records runs"

        # components of femcore
        nyx.femcore.meshing -> nyx.femcore.mesh "Builds"
        nyx.femcore.meshing -> nyx.femcore.adapt "Adapts to sizes"
        nyx.femcore.meshing -> nyx.geom "Drives the geometry"
        nyx.femcore.adapt -> nyx.femcore.mesh "Reads and rebuilds"
        nyx.femcore.fem -> nyx.femcore.mesh "Element numbering"
        nyx.femcore.fem -> nyx.femcore.section "Axis names"
        nyx.femcore.fem -> nyx.femcore.material "Properties per domain"
        nyx.femcore.fem -> nyx.femcore.linalg "Assembled matrices"
        nyx.femcore.material -> nyx.femcore.expr "Expressions in material files"
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

        component nyx.engine "Engine" "Inside the engine." {
            include *
            autolayout tb
        }

        component nyx.femcore "femcore" "Inside the FEM core." {
            include *
            autolayout tb
        }

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

        dynamic nyx "AICommand" "An AI host builds a model and solves it through the same commands." {
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
