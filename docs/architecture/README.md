# nyx architecture

Written 2026-11-03 against main `a047506`; brought up to date at the end of the restructure (phase 9, 2026-11-13: `workspace.dsl`, `diagrams.md`, ADRs 24 to 65, `CODE_TO_ADR.md`, `ROUGH_SPOTS.md`). The ADRs of the early days describe the structure before the restructure; those that no longer hold say so in their Status line and point to the one that replaced them (ADR 65).

**Validated** on 2026-11-13 with the Structurizr CLI (the release zip, run with Java: `java -cp "<cli>/lib/*" com.structurizr.cli.StructurizrCliApplication validate -workspace workspace.dsl`; it exits 0 and prints nothing when the DSL is valid, and names the line of an error otherwise). The export to Mermaid (`-format mermaid`) writes the six views. `diagrams.md` stays hand-written, because the export has no colour classes or notes; keep the two in step.

| What | Where |
|---|---|
| C4 model (context, containers, components, two dynamic views) as Structurizr DSL | [`workspace.dsl`](workspace.dsl) |
| Rough spots the model showed (28 items, with evidence, severity and a fix direction; the resolutions are at the end) | [`ROUGH_SPOTS.md`](ROUGH_SPOTS.md) |
| Decisions found in the history, docs and code (65 ADRs, Nygard format; the first 23 written retrospectively) | [`adr/`](adr/), index in [`ADR_INDEX.md`](ADR_INDEX.md) |
| Every file of the repository mapped to the ADRs that cover it | [`CODE_TO_ADR.md`](CODE_TO_ADR.md), checked by `check_code_to_adr.sh` |
| Interface contracts: exemplar files and the tests that pin them | [`contracts/`](contracts/README.md) |
| Mermaid diagrams of the final structure (hand-written, render on GitHub or mermaid.live) | [`diagrams.md`](diagrams.md) |
| **The restructure plan** (authoritative, self-contained, for the next working session); status in `docs/RESTRUCTURE_STATUS.md` | [`RESTRUCTURE_PLAN.md`](RESTRUCTURE_PLAN.md) |
| Rationale behind it (earlier proposals and decision rounds, superseded as plans) | [`PROPOSED_STRUCTURE.md`](PROPOSED_STRUCTURE.md), [`PROPOSED_STRUCTURE_2.md`](PROPOSED_STRUCTURE_2.md) |
| Decision rounds 2 and 3 (explanations: sbs removal, multiphysics, materials, plug-ins, goldens, crate cost, dispersion) | [`DECISIONS_ROUND_2.md`](DECISIONS_ROUND_2.md), [`DECISIONS_ROUND_3.md`](DECISIONS_ROUND_3.md) |

## Viewing the model

Structurizr Lite on this folder (`docker run -it --rm -p 8080:8080 -v <this folder>:/usr/local/structurizr structurizr/lite`) opens `workspace.dsl` and imports `adr/` as decisions. `structurizr-cli export -workspace workspace.dsl -format mermaid` produces diagrams without a server. **The DSL was written by hand and has not been run through Structurizr**: no Structurizr tool is installed on this machine and I did not download one. If it reports a syntax error, that is a defect in the file, not in the model.

Elements carry three properties: `code` (where it lives), `adr` (which decisions explain it) and, where there is one, `contract`. The tag `rough` (red outline, red dashed lines) marks an element or relationship that appears in `ROUGH_SPOTS.md`.

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
