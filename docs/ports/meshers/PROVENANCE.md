# Provenance and licence notice for `meshers/`

Plain statement first, details below.

* **No source code was copied, translated or adapted from any other project.** Not from Mmg, ParMmg, TetGen, Gmsh, CGAL, libMesh, Netgen,
  rapidmesh (forbidden by the project's licence rules), and not from Omega_h, PRAgMaTIc, NASA refine, deal.II, MFEM, PETSc/SLEPc or
  `rustc-hash` either. **No third-party source file was opened** while writing `femmesh`; no such source was available to the writing
  session and the licence rules would have required a LICENSE check and a provenance row first.
* Everything is original Rust written from (a) published mathematical formulas and algorithms, recalled from memory and
  **not re-read** during the work, and (b) first-principles derivations noted in the module docs. Where a citation is "recalled" it is
  marked so in [`femmesh/docs/REFERENCES.md`](femmesh/docs/REFERENCES.md); treat recalled citations as unconfirmed until someone has the paper.
* `femmesh` has **no dependencies** (`[dependencies]` is empty). Nothing is linked, called or shelled out to.
* Licence of this folder: **BSD-2-Clause** (`femmesh/LICENSE`). The licence text was copied from this repository's solver port
  (branch `claude/stoic-curie-u5aij3`, `LICENSE`, same copyright holder line) so that both ports carry the same terms. That is the only file copied from elsewhere in this repository.
  nyx's own licence is unknown to the author of this folder; check it before moving code into `femcore`.

## Where ideas come from

"Standard method" means a textbook or widely published technique that I implemented from understanding. "Inspired by" means the
overall design follows a known system, taken from general public knowledge of how it works, with no source read.

| Area (file) | Relationship | What exactly |
|---|---|---|
| Metric algebra, intersection by simultaneous reduction, log-Euclidean mean (`adapt/metric.rs`) | standard method | Anisotropic-adaptation literature (Frey and George; Alauzet; Arsigny et al. 2006). Verified by property tests, not by comparison with any code. |
| Edge length with geometric size variation, `l = (la-lb)/ln(la/lb)` (`adapt/metric.rs`) | standard method | Re-derived as the integral of `|e|/h(t)` with log-linear `h`; checked numerically in a test. |
| L^p-optimal Hessian-to-metric and complexity normalisation (`lp_metrics`) | standard method | Loseille and Alauzet, continuous mesh framework (recalled). Fused kernel, global eigenvalue floor, log-space determinant and Illinois-type root finder for the complexity scale are my own design. |
| Gradation (`grade`) | own parameterisation | Linear size-growth law `h' = h + g*dist`. Not verified to equal Alauzet's formulation. |
| Operator set: split / collapse / flip / smooth driven by metric edge length, thresholds `1/sqrt(2)` and `sqrt(2)`, "fraction of edges in band" statistic (`adapt/ops.rs`, `drive.rs`) | **inspired by** Omega_h, PRAgMaTIc, Mmg style adaptation | The *idea* and the conventional thresholds are common knowledge in the field. The implementations, the quality gates, the link-condition check, the ring-centroid smoothing candidate, the `converge_frac` stopping rule, the split/flip interleaving and the frozen-chain handling are my own. I have not seen how any of those codes implements them. |
| Boundary/interface/corner classification from labels (`adapt/topo.rs`) | own design | Written from the handoff's description of nyx's labels. |
| Geometry oracle with analytic projectors (`adapt/geom.rs`) | own design | |
| Structured PML meshes with conforming prism/hex layers; simplex split as a cone from the smallest vertex id and quad diagonal through the smallest id (`adapt/pml.rs`, `cellmesh.rs`) | own design | The "diagonal through the smallest global id" device is a standard way to make independent splits conform; the specific construction is mine. Region/label encoding is my own and is *not* nyx's `chips.rs` scheme. |
| PML complex coordinate stretching, polynomial profile, design reflection (`validate/helmholtz.rs`) | standard method | Berenger 1994 / textbook derivation (recalled). |
| Scalar Helmholtz FE harness (`validate/`) | own, standard FE | P1/Q1/P1xlinear/trilinear elements, 3-point triangle, 4-point tet and Gauss rules (standard tabulated rules, rebuilt from their defining points), banded LU, inverse iteration. **Not** derived from PETSc/SLEPc. |
| Reverse Cuthill-McKee (`validate/band.rs`) | standard algorithm | Cuthill and McKee 1969 / George; implemented from the description. |
| Jacobi symmetric eigen-solver (`adapt/metric.rs`) | standard algorithm | |
| Fx-style hasher (`adapt/topo.rs`) | standard scheme, reimplemented | The multiply/rotate/xor step and its 64-bit constant are the well-known Fx scheme used in Firefox and `rustc-hash` (MIT/Apache-2.0). Written from memory as a 10-line function; `rustc-hash` source was not opened. If this matters to a licence audit, replace it or add the upstream notice. |
| Deterministic fixed-chunk parallel reductions (`par.rs`) | own, standard technique | |
| `.fmsh`: varint (LEB128), zigzag, run-length, delta coding (`io/fmsh.rs`) | standard encodings | |
| `.fmsh`: Gorilla-style XOR float coding | **inspired by** Pelkonen et al., "Gorilla", VLDB 2015 (recalled) | Byte-aligned variant of my own; not their bit stream. |
| `.fmsh`: CRC-32 (IEEE) | standard algorithm | Bitwise implementation of the public polynomial; test vector `123456789 -> 0xCBF43926`. |
| `.fmsh`: adaptive binary range coder (`io/entropy.rs`) | **inspired by** the LZMA range coder | The carry-propagating low/range/cache scheme with adaptive bit probabilities is the public LZMA scheme. Re-implemented from memory (LZMA uses 11-bit probabilities; this uses 15-bit, and the order-0/1 byte models are my choice); no xz, 7-zip or LZMA SDK source was opened. A reader comparing it with LZMA's `rc` will see the same structure; this is by design, and the algorithm is a published one. |
| `.fmsh` layout, byte planes, dictionary coding, result steps/fields | own design | |
| Test fixtures, analytic fields (Gaussian, layers) | own | |
| Degree-5 7-point triangle rule in `examples/adapt_convergence.rs` | standard tabulated rule | Radon/Dunavant-type points, formula `(6 +- sqrt 15)/21`. |

## What the next session must do about provenance

1. Keep this file true. If any third-party source is read, check its LICENSE first (rule 2 of the brief), and add a row here
   (file, derived from, upstream file, licence), reproduce the upstream notice in `femmesh/THIRD_PARTY_NOTICES.md` and put an SPDX header on the derived files.
2. Never read Mmg, ParMmg, TetGen, Gmsh, CGAL, libMesh, Netgen or rapidmesh. Omega_h, PRAgMaTIc, refine, deal.II, MFEM and PETSc only after the licence check.
3. When code moves into nyx `femcore`, carry this file's table (or a pointer to it) with it.
4. The original brief is in [`docs/nyx_meshing_adaptivity_handoff.md`](docs/nyx_meshing_adaptivity_handoff.md); section 2 there is the authoritative rule set.
