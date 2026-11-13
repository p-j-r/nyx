# nyx C4 set

The full C4 model of nyx at the end of the restructure (phase 9, 2026-11-13): **context, containers (one per crate), components (every container) and code (every component)**, plus two dynamic views. Structurizr is the source for levels 1 to 3; level 4 is a class-diagram page per component generated from the Rust source. Short summary and findings: [MEMO.md](MEMO.md).

| Level | Folder | What | Made from |
|---|---|---|---|
| 1 Context | [`1-context/`](1-context/context.md) | nyx, the researcher, the AI host, the results viewer, the workspace folder | `workspace.dsl`, view `Context` |
| 2 Containers | [`2-containers/`](2-containers/containers.md) | the 12 crates and the 4 kinds of file they write; the arrows are the layering table of `crates/nyx-engine/tests/layering.rs` | view `Containers` |
| 3 Components | [`3-components/`](3-components/README.md) | 91 components, one view per container (femcore also in three parts); a component is a group of modules with one job | generated components, views `Components_<container>` |
| 4 Code | [`4-code/`](4-code/README.md) | one page per component: files, what it uses and what uses it, a class diagram of its public types, a complete table of its public items | generated from the source |
| Dynamic | [`5-dynamic/`](5-dynamic/README.md) | a solve from the window, an AI host building and solving | dynamic views |

Other folders: `export/` (every Structurizr view as Mermaid, PlantUML and workspace JSON), `viewer/index.html` (all 110 diagrams in one page with a menu; opens in a browser, loads Mermaid from a CDN), `model/` (generated Structurizr fragments; in the standalone folder they are in `structurizr/model/`), `spec/` (the component spec), `tools/` (generator and build scripts), `reports/` (what the source scan found).

## Where things are

* **In the repository**: this folder is `docs/architecture/c4/`. The master workspace is [`../workspace.dsl`](../workspace.dsl) (it includes `model/`); the ADRs it links to are [`../adr/`](../adr/).
* **Standalone copy**: `Documents\opomech\nyx-c4\` (same folders; master workspace `structurizr/workspace.dsl`, ADRs copied to `structurizr/adr/`). Made by `bash tools/make_standalone.sh`. It is a build output of this folder, not a second source: change the spec or the source, rebuild both.

## Viewing

* **GitHub or VS Code** render the `.md` pages (Mermaid blocks).
* **`viewer/index.html`** in a browser (needs internet for the Mermaid script); a fit-to-width switch is at the top.
* **Structurizr Lite** on `docs/architecture` (or on `structurizr/` in the standalone folder): `docker run -it --rm -p 8080:8080 -v <folder>:/usr/local/structurizr structurizr/lite`. It shows the same views with its own layout and imports the ADRs as decisions.
* **Structurizr CLI**: `structurizr-cli validate -workspace workspace.dsl`; `export -format mermaid|plantuml|json`.

## How it is made

```bash
bash tools/build.sh                 # in the repository: gen.pl, validate, export (mermaid, plantuml, json), pages, viewer
bash tools/make_standalone.sh       # the standalone folder, same steps
```

Needs Perl, Java and the Structurizr CLI release (`STRUCTURIZR_CLI=<folder with lib/>`; found automatically when it sits next to the repository folders as `structurizr-cli/unpacked`). No Python, no network, no cargo.

1. [`spec/components.spec`](spec/components.spec) says which files form which component (every `.rs` file of every crate's `src/` must be in exactly one, or the build stops), with a title, the ADRs and a description of each.
2. [`tools/gen.pl`](tools/gen.pl) reads the Rust source: public items, fields, variants, signatures, `impl` blocks, and every `use` and `crate::` / `femcore::` / `nyx_*::` path (re-exports followed inside a crate; comments, strings and `#[cfg(test)]` modules ignored). From that it writes the component definitions, the **relationships between components** (and from a component to another crate's container, with the modules named in the label), the component views and all code pages.
3. The Structurizr CLI validates and exports; `gen.pl --pages` writes the pages and the viewer from the exports.

What is hand-written: the people, systems, containers, container relationships and dynamic views in `../workspace.dsl`, and the component list and descriptions in the spec. What is derived: everything else. `reports/source-scan.txt` compares the crate-to-crate references found in the source with the container relationships of the DSL: every reference is declared, and every declared relationship but `mcp -> engine` has a reference (that one is declared in the spec because `nyx-mcp` is one file shared by two components, which the scan does not read).

## Limits

* `gen.pl` is a reader of Rust text, not a compiler. Macros that define items, `pub use` of glob re-exports across crates and trait-method calls are not followed; a reference it cannot resolve is attributed to the crate root. Test code (`tests/`, `examples/`, `#[cfg(test)]`) is outside the model.
* A component edge means "a module of this component names an item of that one", not that a call happens at run time. Rust allows module cycles inside a crate and there are some (e.g. most of the engine's components name each other, because `impl Engine` blocks are spread over its files); the pages show them as they are.
* A class diagram shows at most 36 types and a few members each. The table under it is complete (400 rows at most). A file that exports nothing (a binary's `main.rs`) shows its private items, marked `-`.
* Dependencies of `femcore` on `faer`, `rayon`, `hdf5` and other external crates are not modelled.
* Mermaid exported by Structurizr uses HTML labels. GitHub renders them; some viewers do not.

## Keeping it true

Rebuild after a change of crate layout or a new file (the build fails on a file that belongs to no component). Commit the generated files with the change, as `CODE_TO_ADR.md` is. A decision that changes the picture still gets an ADR (see `../README.md`).
