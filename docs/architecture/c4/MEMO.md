# Memo: the full C4 set for nyx

2026-11-13, on top of `main` (phase 9 plus the Structurizr validation). For paul.

**What you have.** C4 at all four levels, for all of nyx. Level 1: context. Level 2: the 12 crates as containers, with the four kinds of file they write; the arrows are the layering table. Level 3: 91 components, a view for each container (femcore, with 22 components and about 80 internal arrows, also in three smaller parts). Level 4: a page for each of the 91 components with a class diagram of its public types, what it uses, what uses it, and a table of all 2,500 public items. Two dynamic views (a solve from the window, an AI host building a model). Everything is a Structurizr model validated with the CLI (exit 0), exported to Mermaid, PlantUML and JSON, and all 110 Mermaid diagrams were checked to render in Chrome.

**Where.** In the repository: `docs/architecture/c4/` (start at its `README.md`), master model `docs/architecture/workspace.dsl`. Standalone: `Documents\opomech\nyx-c4\` with subfolders `1-context`, `2-containers`, `3-components`, `4-code`, `5-dynamic`, `structurizr`, `export`, `viewer` (open `viewer\index.html` for all diagrams in one page).

**How it was made, and why you can trust the lower levels.** Levels 1 and 2 and the dynamic views are the hand-written model from phase 9. Components and everything below are read from the code: a script lists which files belong to which component (a spec I wrote, 91 entries; the build fails if a file belongs to none), then reads the Rust source for the public items and for every `use` and path reference between modules. The relationships between components are therefore what the code says, not what I remember. As a check, the crate-to-crate references found this way are all in the layering table, and every table entry but one has a reference in the code (`nyx-mcp` to `nyx-engine` is declared by hand because the scan skips that one file). The scan covers `src/` of all crates, 78,868 lines; tests and examples are not in the model.

**What the code view showed.** Nothing that contradicts the layering. Two things worth knowing. First, module cycles inside crates are common (they are legal in Rust): most of the engine's components name each other partly because `impl Engine` blocks are spread over its files, and six of the seven components of `geom` are one cycle. The cross-crate picture is clean, the inside of the engine is not layered. Second, in `femcore`, `adapt::work` (the simplex working mesh) and `mesh` name each other, and so do the 2D/3D meshers and the structured PML; these are the places to look if femcore is ever split. I did not add any of this to `ROUGH_SPOTS.md` since none is a defect found by a test, but I can.

**Limits.** The class diagrams are capped (36 types, a few members each) and the tables are complete; the reader is a text reader, not the compiler, so macro-generated items and glob re-exports across crates are not seen. Details in the README.

**Not done, on purpose.** The worked 1550 nm optical-mode use case with a C4 dynamic view (a later step). `diagrams.md` (hand-written Mermaid) is kept and now points to this set.

**For you to run** (I did not move `main` and did not push):

```bash
jj bookmark set main -r restructure-9b-c4
jj new main
```

The standalone folder is outside the repository and needs no merge. To rebuild it after a change: `bash docs/architecture/c4/tools/make_standalone.sh`.
