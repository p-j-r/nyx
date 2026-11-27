# nyx architecture

Written 2026-11-03; brought up to date at the end of the restructure (phase 9, 2026-11-13). Since then the ADRs (now up to 133) and `CODE_TO_ADR.md` have been kept current; `workspace.dsl`, `diagrams.md` and the C4 spec predate the installer crate `nyx-install` (ADR 122), which they do not model yet. The ADRs of the early days describe the structure before the restructure; those that no longer hold say so in their Status line and point to the one that replaced them (ADR 65).

**Validated** on 2026-11-13 with the Structurizr CLI (the release zip, run with Java: `java -cp "<cli>/lib/*" com.structurizr.cli.StructurizrCliApplication validate -workspace workspace.dsl`; it exits 0 and prints nothing when the DSL is valid, and names the line of an error otherwise). The export to Mermaid (`-format mermaid`) writes the six views. `diagrams.md` stays hand-written, because the export has no colour classes or notes; keep the two in step.

| What | Where |
|---|---|
| **Full C4 set**: context, containers, components and code pages per component, exports and a browser viewer, read from the Rust source. The spec, model fragments and tools are tracked; the generated pages are not (see its README) | [`c4/`](c4/README.md) |
| Master Structurizr model (people, systems, containers, container relationships, two dynamic views; includes the generated components of `c4/`) as Structurizr DSL | [`workspace.dsl`](workspace.dsl) |
| Rough spots the model showed on 2026-11-03 (a snapshot against the old layout, partly overtaken; ids cited in code comments as `R<n>`) | [`ROUGH_SPOTS.md`](ROUGH_SPOTS.md) |
| Decisions found in the history, docs and code (119 ADRs, Nygard format, numbered up to 133 with some numbers reserved; the first 23 were written retrospectively) | [`adr/`](adr/), index in [`ADR_INDEX.md`](ADR_INDEX.md) |
| Every file of the repository mapped to the ADRs that cover it | [`CODE_TO_ADR.md`](CODE_TO_ADR.md), checked by `check_code_to_adr.sh` |
| Interface contracts: exemplar files and the tests that pin them | [`contracts/`](contracts/README.md) |
| Mermaid diagrams of the final structure (hand-written, render on GitHub or mermaid.live) | [`diagrams.md`](diagrams.md) |
| The restructure plan (executed, phases 0 to 9), its proposals, decision rounds and handoff log | [`../history/`](../history/README.md) |

## Viewing the model

Structurizr Lite on this folder (`docker run -it --rm -p 8080:8080 -v <this folder>:/usr/local/structurizr structurizr/lite`) opens `workspace.dsl` and imports `adr/` as decisions. `structurizr-cli export -workspace workspace.dsl -format mermaid` produces diagrams without a server. The master `workspace.dsl` includes the generated component model of `c4/model/`; see `c4/README.md` to rebuild, view and export it.

Elements carry three properties: `code` (where it lives), `adr` (which decisions explain it) and, where there is one, `contract`. The tag `rough` (red outline, red dashed lines) marks an element or relationship that appeared in the rough-spots snapshot.

## Running the checks

```bash
cargo test --release -p femcore --test contracts
cargo test --release -p geom --test contracts
cargo test --release -p nyx-store --test contracts
cargo test --release -p nyx-mcp --test contracts
bash docs/architecture/check_code_to_adr.sh
```

`cargo test -p nyx-store --test contracts -- --include-ignored` shows the one known-broken contract (R25).

## Keeping it true

* A new file needs a row in `CODE_TO_ADR.md` (the script fails otherwise).
* A change to an interface changes an exemplar and `contracts/VERSIONS.md` in the same commit (ADR 19).
* A decision that changes the picture changes `workspace.dsl` and gets a new ADR; an ADR is never edited to say something different, it is superseded by a later one.
