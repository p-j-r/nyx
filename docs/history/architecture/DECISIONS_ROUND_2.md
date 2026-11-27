# Decisions, round 2: answers and explanations

> **Merged into `RESTRUCTURE_PLAN.md`, which is authoritative.** This file is kept for the reasoning behind the decisions.

2026-11-03. Proposal only; nothing implemented. Follows `PROPOSED_STRUCTURE_2.md`.

## 1. Remove sbs altogether

Accepted. The physics set becomes `dfs` (electrostatics, heat), `sms` (solid mechanics), `ems` (EM waves).

**What goes**
* `sbs.rs` (1,480 lines), `sbs_ring.rs`, `chips.rs`, `Physics::Brillouin`, `BrillouinSpec` and everything only it used: photoelastic data, moving-boundary term, viscosity, `ring_radius`, `mirror_axis`, the pump/Stokes mode pairing.
* Tests and goldens: 3 Brillouin and 2 ring-Brillouin golden tests, `tests/brillouin.rs`, `ring_frame.rs` parts that exercise Brillouin, examples `wg800`, `chip_*`, `rib_*`.
* In the contracts: the `sbs` branch of the problem-spec exemplar.
* ADR 22 and the chips part of ADR 21 are marked superseded.

**What stays, because it is not Brillouin**
* `sms_kz` (solid modes with an out-of-plane wave number: phononic waveguides) and the axisymmetric variants `ems_axisym`, `sms_axisym` (rings and disks). They were built largely for the Brillouin work but are ordinary physics studies. If you want to go smaller still, say which.
* The analysis (Stokes, polarization) stays; it only ever needed the optical mode fields.

**Side effects to know about**
* 5 of the 7 failing golden tests are Brillouin or ring-Brillouin. They disappear; the two left are `ems` (`golden_soi_wire_modes`, `golden_soi_wire_with_sloped_wall_hybrid_modes`), so decision 6 shrinks to two tests.
* The 800 nm intermodal Brillouin study (`runs/wg800`, awaiting your validation per my notes) uses this code. I suggest tagging the commit before removal (`sbs-last`) so it can be rebuilt from history; it would not be in the main tree.
* The 1,900 lines of `sbs` also depended on `ems` and `sms` internals (the optical modes, the kz solver). Removing it should let those two be simplified; that is a clean-up for later, not a requirement.

## 2. A few physics loaded, not twenty-five

Agreed, and here are real examples of what "a few" means. Couplings people actually do:

| Coupling | Physics involved | Direction |
|---|---|---|
| Thermal expansion | heat then solid | one way (this exists today as `SolveState`) |
| Joule heating | electrostatics/current then heat | one way |
| Thermo-optic shift | heat then EM waves (dn/dT) | one way |
| Electro-optic (Pockels) | electrostatics then EM waves | one way |
| Piezoelectric / electrostriction | solid and electrostatics | two way |
| Opto-mechanics (forces from light) | EM waves and solid | two way |

So a project usually loads **one to three** physics. Design consequences:
* The registry may know many physics (they are compiled in), but **a project names the few it uses**; the window and MCP show only those.
* Couplings are **explicit links** between a field one physics exports and a quantity another imports (`dfs.T` into `sms.thermal_load`, `dfs.T` into `ems.dn_dT`). Each physics only declares what it exports and imports. Adding the 25th physics adds one set of declarations, not 24 pairwise adapters.
* Version 1 of the engine supports **one-way chains** (solve A, feed B). That covers the first four rows. Two-way coupling (rows 5 and 6) needs iteration or a joint system and is a later step.
* No hard cap on the number loaded; a soft warning above four is enough.

## 3. Retire XML

Accepted. Projects and specs become JSON. The XML reader stays in exactly one place (a `legacy` module in the engine), used only by the migration that turns old `.nyx.xml` files into the new format; it is deleted once your old files are migrated. `xmlio`, the XML option of the file dialogs and the XML round-trip tests go away.

## 4. Materials: easy to create, complete when needed

Today's `material.rs` can only describe four crystal classes by named constants (cubic, hexagonal, trigonal, orthotropic), and `sapphire()`, `silicon()` are written in code. You want to be able to write the **entire elastic matrix in reduced (Voigt) form**. Proposal: a material is a small TOML file where every property can be given in whichever form is natural, and nyx checks it.

Rules
* **Units are always written** (`unit = "GPa"`), never assumed.
* **Voigt order is the textbook one**: 1 = xx, 2 = yy, 3 = zz, 4 = yz, 5 = xz, 6 = xy, with engineering shear strain, in the **crystal frame**. (This is already what `material.rs` uses internally before it converts to the solver's order, so nothing changes numerically.) The orientation (Bunge Euler angles) places the crystal in the wafer frame.
* **Three ways to give stiffness**, all producing the same 6x6 matrix:
  1. by class and independent constants (`class = "cubic"`, `c11`, `c12`, `c44`); classes: isotropic (`E`, `nu` or `lambda`, `mu`), cubic, hexagonal, trigonal, tetragonal, orthotropic, monoclinic, triclinic;
  2. the reduced matrix, **upper triangle only** (the rest is symmetric), or the full 6x6;
  3. the compliance matrix `S` instead of `C`.
* **Checks with plain messages**: symmetric, positive definite (stable), declared class consistent with the numbers (a "cubic" matrix with c11 different on the diagonal is reported), unknown keys reported rather than ignored.
* The same idea for every tensor property: permittivity, conductivity and thermal expansion as one number (isotropic), three (diagonal) or the six independent entries; later piezoelectric `d` as 3x6.
* Optional `notes`, `source` (a reference), and a `temperature` table later. A `.nyx` file embeds the materials it uses, so it opens anywhere.

```toml
# silicon.toml  (illustrative values)
name = "Si"
source = "Hopcroft et al., J. Microelectromech. Syst. 2010 (illustrative: check before use)"
notes = "single crystal, crystal axes along the cube edges"

[optical]
eps_r = 12.08                 # scalar means isotropic

[elastic]
unit = "GPa"
voigt = [                     # upper triangle, rows 1..6; symmetric completion is implied
  [165.7, 63.9, 63.9,  0.0,  0.0,  0.0],
  [       165.7, 63.9, 0.0,  0.0,  0.0],
  [              165.7, 0.0,  0.0,  0.0],
  [                     79.6, 0.0,  0.0],
  [                           79.6, 0.0],
  [                                 79.6],
]
density = { value = 2329, unit = "kg/m3" }

[thermal]
expansion = { value = 2.6e-6, unit = "1/K" }
```

The same silicon by class: `[elastic] class = "cubic"`, `c11 = 165.7`, `c12 = 63.9`, `c44 = 79.6`, `unit = "GPa"`.

## 5. Compile-time versus run-time plug-ins, in detail

The question: **can someone add a new physics to nyx without rebuilding nyx?**

### Compile-time (what I proposed)
A physics is a Rust crate. The nyx program is built with a chosen list of them (cargo features). A user who wants a physics that is not in their copy installs a different build of nyx. A developer who writes a new physics adds a crate and one line to the registry, then rebuilds.
* Good: one executable, no loading code, no security question, the compiler checks everything, debuggable, and the determinism rules (same input, same numbers) are enforced by tests in the same repository.
* Limit: a third party cannot add a physics to an installed nyx. "Downloadable extras" are only **data**: `.nyx` files, material files, geometry libraries. Data cannot run code, so opening one is safe.

### Run-time plug-ins
A physics ships as a separate file that nyx loads while it is running. There are three ways, each with a different price:

| Way | What it is | Price |
|---|---|---|
| Native library (`.dll`) | a Rust or C library loaded with `libloading` | Rust has **no stable binary interface**: the plug-in must expose a C-style interface (flat arrays for the mesh and fields, a versioned header), be built with the same compiler assumptions, and a crash in the plug-in takes nyx down. Running a downloaded native library means running arbitrary code with the user's rights. |
| WebAssembly module | a sandboxed `.wasm` file run by an embedded runtime | Safe and portable, but slower, awkward for multithreading and big memory, and the solver library cannot just be linked in: the plug-in would need its own linear algebra or a host call interface for it. |
| Separate process | the physics is its own program that nyx starts and talks to over a pipe, as the MCP server does | Safest (a crash does not kill nyx, can run on another machine), any language; meshes and fields are copied across the pipe; you design and version a wire protocol. |

What changes in the code if you want plug-ins:
1. The `Physics` trait must cross a boundary, so everything passed through it has to be **plain data**: descriptors and specs as JSON, meshes and fields as flat numeric arrays, results as arrays plus names. (My trait is already close; the cost of keeping this option open is only discipline.)
2. A **loader and registry**: where plug-ins live, how nyx discovers them, a version check between nyx and the plug-in, and what the window shows when one is missing or incompatible.
3. A **trust model**: signing, a user prompt, or only allowing sandboxed formats.
4. The tests (contracts, goldens) must run against the plug-in interface too, doubling the interface to keep stable.
5. Distribution and support: builds per operating system, and "which plug-in version works with which nyx" becomes a real question.

My recommendation: **compile-time now, with the trait written so that only plain data crosses it.** That keeps the cheap path today and leaves the "separate process" way open later, because that is the way that fits nyx best (it already runs solves as detached workers and already has an MCP server). Decide only if you expect people outside the core developers to add physics soon; if not, nothing needs to change.

## 6. The golden platform or tolerance, in detail

### What a golden is
A golden file stores the numbers a kernel produced for a fixed test problem (for example the effective index and field samples of a silicon wire), taken at a moment when the kernel was known to be right. The test solves the same problem again and checks that the new numbers equal the stored ones to a tolerance (1e-9 for diffusion, 1e-7 for EM waves, 1e-6 for the Brillouin cases). It exists so that a refactor which should change nothing (like the frame port) can prove it changed nothing.

### Why 7 fail on this machine
The failures are tiny: about 1e-6 relative in fields and up to 2e-5 in one effective index. They are the same with one thread, and they are already present at the commit that wrote the goldens. So they are not caused by a code change; the stored numbers came from a different machine. Floating-point results depend on the processor and compiler:
* fused multiply-add, vector width and the way the matrix library (faer) blocks a matrix multiply change the last bits of every sum;
* an eigensolver stops when its residual is below a tolerance, so tiny differences in the last bits change where it stops, and the answer differs at about the solver's own tolerance (around 1e-6);
* the compiler version can reorder operations.

All of these are within the accuracy of the method; none is a physics error. But a test tolerance tighter than the solver's own tolerance fails on any machine other than the one that made the file.

### The options and what each changes
| Option | What you do | Effect | Downside |
|---|---|---|---|
| A. Reference platform | declare one environment the oracle (a named machine or CI image) and enforce goldens strictly only there | strict regression detection (1e-7) | this Windows PC can never pass in strict mode unless it is the reference; other machines only warn |
| B. Looser tolerance | raise the tolerance per case to what the solver can actually reproduce (about 1e-4 relative for n_eff) | passes on every machine | a regression smaller than 1e-4 goes unseen |
| C. Both | strict on the reference platform, loose elsewhere | best detection and everyone can run tests | a little more setup |
| D. Closed-form benchmarks | test against exact solutions (slab n_eff, cavity frequencies, analytic plate) instead of stored numbers | independent of any machine | not every case has a closed form; keep a few goldens for regression |

### What changes after removing sbs
Five of the seven failing tests are Brillouin and disappear. Two remain, both `ems` silicon-wire cases. For those, option B at about 1e-4, with a closed-form check of the slab index alongside, is enough. The question for you reduces to: **is this Windows machine the reference?** If yes, the simplest honest action is to regenerate those two goldens here, say so in the commit, and keep strict tolerances; other machines would then fail by the same 1e-6. If not, loosen the tolerance. Either way no physics changes; only the test configuration does.
