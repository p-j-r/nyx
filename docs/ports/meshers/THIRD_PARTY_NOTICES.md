# Third-party notices

`femmesh` has **no third-party dependencies** (`[dependencies]` is empty) and **no code copied, translated or derived from
any third-party project**. The full statement, including which algorithms are standard methods and which designs are only
*inspired by* known systems (Omega_h / PRAgMaTIc / Mmg-style operator sets, the LZMA range-coder scheme, the Fx hash scheme,
Gorilla float coding), is in [`../PROVENANCE.md`](../PROVENANCE.md). Citations: [`docs/REFERENCES.md`](docs/REFERENCES.md).

No source of Omega_h, PRAgMaTIc, NASA refine, deal.II, MFEM, PETSc, `rustc-hash` or the LZMA SDK has been opened. If that
happens (see the licence rules in `../docs/nyx_meshing_adaptivity_handoff.md`, section 2), add a provenance row below, reproduce
the upstream notice here, and put an SPDX header on the derived files.

Never read, implement from papers only: Mmg, ParMmg, TetGen, Gmsh, CGAL, libMesh, Netgen, rapidmesh.

| femmesh file | derived from | upstream file | licence |
|---|---|---|---|
| (none) | | | |

The only file copied from elsewhere is `LICENSE` (BSD-2-Clause text, from this repository's solver port, branch `claude/stoic-curie-u5aij3`).
