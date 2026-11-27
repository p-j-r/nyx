# Design history

Plans, proposals, studies, handoff notes and status logs from the work that built nyx. They record **how and why** the code got to its present shape. They are **not current documentation**: most were written before the work landed, some describe modules that no longer exist, and their links and file paths were not kept up to date.

For the current state read the [README](../../README.md), [How it works](../HOW_IT_WORKS.md), the [architecture overview](../architecture/README.md) and the [ADRs](../architecture/adr/); an ADR wins over anything in this folder.

| Folder or file | What it was | What replaced it |
|---|---|---|
| `RESTRUCTURE_STATUS.md`, `architecture/RESTRUCTURE_PLAN.md`, `architecture/PROPOSED_STRUCTURE*.md`, `architecture/DECISIONS_ROUND_*.md` | the restructure plan, its drafts and its session handoff log (phases 0 to 9, closed 2026-11-13) | `CLAUDE.md`, `docs/architecture/README.md`, ADR 20 and the ADRs of each phase |
| `UX_STATUS.md`, `ux/UX_*.md` | the UX review trail of 2026-11-17 | ADR 123 and 124 (`docs/ux/tokens.json` stays) |
| `LATENCY_STATUS.md`, `FRAME_UNIFICATION_STATUS.md`, `KERNEL_PORT_STATUS.md` | handoff logs of finished work | ADR 1, 10, 98 |
| `GEOM_HANDOFF*.md`, `GEOMETRY_ARCHITECTURE.md`, `GEOMETRY_GUI_PLAN.md` | briefs for the session that wrote `crates/geom` | `crates/geom` and ADR 2, 3, 18 |
| `MCP_ASSISTANT_PLAN.md`, `WORKFLOWS_DESIGN.md`, `SWEEP_CACHING_PROPOSAL.md`, `PLOTTING_PLAN.md`, `RESULTS_FORMAT_STUDY.md` | designs for work that has since been built differently | ADR 94, 14, 87 and 96, 99, 107 |
| `PML_FEATURE_REQUEST.md`, `PML_IMPLEMENTATION_PLAN.md`, `SQUADDS_BUG_FEATURE_REQUESTS.md` | specifications handed to an implementing session | ADR 128, 125 to 127 |
| `GDSFACTORY_INTEGRATION.md` | design note, nothing built | `import_gds` (ADR 121) |
| `chip_report.md`, `photonic_use_cases.md` | reports on the Brillouin ring chips and use cases, removed in ADR 24 | none |
| `ports/` | staging notes of the `femsolve` and `femmesh` ports, written without access to nyx. The licence and provenance files stay in `docs/ports/` | `femcore::linalg`, `femcore::meshing` |

Removed from the repository altogether (still in the git history of commit `bdf0c78`): the market note, the UI mockups and screenshots of the pre-reset window, the `geom_handoff.zip` snapshot and a one-off results-format measurement script.
