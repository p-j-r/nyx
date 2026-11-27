# Decisions, round 3

> **Merged into `RESTRUCTURE_PLAN.md`, which is authoritative.** This file is kept for the reasoning behind the decisions.

2026-11-03. Proposal only. Follows `DECISIONS_ROUND_2.md`.

## 2. Separate crates do not slow the software

Correct, with two caveats that I would guard by measurement.

* A Rust workspace of many crates compiles into **one executable**. Crate boundaries are a source-organisation and compile-time matter; at run time there is no call between programs.
* The only run-time cost of the physics interface is the registry lookup and one virtual call to `Physics::solve` per solve. A solve takes seconds to minutes; the call costs nanoseconds.
* Where a split *could* hurt: a small function in one crate called millions of times from another (for example an element-matrix helper in femcore called per element from a physics crate). Generic functions are compiled into the caller and inline normally; non-generic ones only inline across crates if marked `#[inline]` or if link-time optimisation is on. The release profile today sets only `opt-level = 3`.
* Guard: before moving anything, record timings of the existing benchmarks (`examples/bench`, `profile3d`, `mem_profile`, `cavity_bench`, `solve_scaling`); after each move compare. If something slows, the fixes are `#[inline]` on the helper or `lto = "thin"` in the release profile, both one-line changes. Build time usually *improves* with more crates, since only the changed one recompiles.
* The hot loops stay inside one crate by design: assembly and kernels in the physics crate, linear algebra in femcore.

## 3. No migration of old projects

Accepted. Consequences, so they are explicit:

* The project and spec format is **not versioned during the refactor**. The `format` key and its reader are introduced **once, at the end**, when the engine, the physics interface and the materials are integrated. The first new files are `format: 1`.
* Old `.nyx.xml` projects and existing store workflows stop opening in the new build. Old builds remain available from git (I would tag the last commit before the switch, `last-old-format`).
* The XML code (`xmlio`, the legacy reader) is deleted with the rest, not kept as a migration module.
* **The physics is still proven unchanged**: the golden tests build their problems in test code, not from project files. They are rewritten against the new problem API and must give the same numbers as the old goldens (that is the proof that the refactor changed no physics), and only then are the goldens retired or kept.
* Contract exemplars (`problem_spec.v1.json`, `store_workflow.v1.json`) are rewritten once at the switch, with a `VERSIONS.md` entry.

## 4. Wavelength-dependent optical properties

### What exists now
Optical materials are one number per domain: `eps_r` (or a refractive index in the Brillouin settings), optionally a tensor. Nothing depends on wavelength; the study's wavelength only sets the frequency.

### Proposal: dispersion lives in the material, evaluated at the study's wavelength
The material file (TOML, round 2) gets an `[optical]` section whose entries are **models**, not numbers. The physics asks the material library for one value, `material.optical(lambda, T)`, and gets the permittivity tensor at that wavelength. The mesh, the solvers and the study never see a formula.

Models (complex values allowed; time convention `e^{j omega t}`, so loss is `n - j k`, matching the wafer-frame convention of ADR 1):

| Model | Form | Typical use |
|---|---|---|
| constant | `n = 3.4757` or `eps_r = 12.08` | quick runs |
| Sellmeier (sum of terms) | `n^2 = 1 + sum_i B_i lambda^2 / (lambda^2 - C_i)` | silica, silicon, sapphire, most dielectrics |
| Cauchy, Schott, Herzberger | polynomial or series in lambda | glasses, old data |
| Lorentz / Drude / Lorentz-Drude | poles in frequency | metals, doped semiconductors, absorption lines |
| tabulated `lambda, n, k` | a table with interpolation | measured data, ellipsometry |
| thermo-optic | `dn/dT` (linear), later a temperature-dependent Sellmeier | heated waveguides, a later coupling with the heat physics |

Details that matter
* **Units are explicit**: wavelength unit declared in the file (`lambda_unit = "um"`, the convention of most published coefficients), so coefficients can be copied from papers and from the refractiveindex.info database without conversion errors. A small importer for that database's file format is cheap and gives hundreds of materials.
* **Interpolation of tables**: monotone cubic (PCHIP) by default, so there is no ringing between points; linear on request. Interpolate `n` and `k` against wavelength, not permittivity. **No silent extrapolation**: outside the tabulated range nyx stops with a message that names the range, or extrapolates only if the file says `extrapolate = "linear"`.
* **Anisotropy**: each principal axis has its own model (ordinary and extraordinary Sellmeier for sapphire, lithium niobate), in the crystal frame, rotated into the wafer frame with the same Euler angles as the elastic tensor.
* **Group index and dispersion** need `dn/dlambda`. Formula models give it analytically; tables via the spline derivative. The analysis crate can then report group index and GVD without finite differences on the solver.
* **Validation**: Sellmeier poles must lie outside the declared range, `n` must be real and positive where the file says "transparent", and a table must be sorted; each failure is a plain message.
* **Temperature** is an optional argument of the evaluation, so the heat-to-optics coupling of round 2 needs no new mechanism.

### Where it runs in a study
* A mode study at one wavelength: evaluate every material once, build the matrices, solve. No nonlinearity.
* A **wavelength sweep** (`neff(lambda)`, dispersion curves): the engine loops over wavelengths, re-evaluating the materials; each solve is independent, so sweeps parallelise and can reuse the mesh. This is the common case.
* **Cavity eigenfrequencies with dispersive materials** are a genuinely nonlinear eigenproblem (the matrix depends on the frequency you are solving for). The ported solver library already contains the tools for it (`nep`, `nleigs`, contour methods); I would only start this after the linear case works, and use a fixed-point iteration (evaluate at the current frequency, solve, repeat) as the first version.

### Where it lives in the structure
* The model types and their evaluation: `femcore::material` (next to the elastic and tensor specs), since every physics may need them.
* The material files: TOML, optional library, embedded in a `.nyx` file for the materials it uses.
* Which optical property a study needs: declared by the `ems` physics in its property descriptors (`optical.eps` at the study wavelength); the engine does the evaluation, so `ems` itself stays free of dispersion code.

## 6. Why the numbers depend on the computer, and what to do about it

First the honest part: **I have not established the cause of the seven failures.** What I know is that they are about 1e-6 relative, identical with one thread, and present at the commit that wrote the goldens. The rest is the standard explanation, which I would test rather than assume.

### What "depends on the computer" means
Computers add numbers in floating point, and the last bits depend on how a sum is ordered and whether a multiply and an add are fused into one instruction. Different processors and different builds of a library choose differently. This is true of every numerical program (MATLAB, COMSOL, NumPy): none of them promises identical bits on different machines unless it is specially built for it. It does **not** mean the answer is wrong. A change of 1e-6 in an effective index is far below the error of the method: the mesh discretisation error is typically 1e-3 to 1e-4 for these problems, so the physical answer is the same on every machine.

### Where nyx is already careful
Same input, same machine, same build gives the same bits for **any thread count** (fixed-chunk reductions, ADR 11). That is the reproducibility that matters for a research workflow. What is not promised today is the same bits **across machines**.

### The likely causes here, in order of suspicion
1. **Iterative eigensolver stopping.** An eigensolver stops when its residual is below a tolerance. If that tolerance is around 1e-6, two machines that differ in the last bit stop at different iterations and the eigenvalue differs by about the tolerance. This would explain 1e-6 exactly. It is also fixable: converge tighter than anyone compares, for example residual 1e-10, and the machine-to-machine difference falls with it, to round-off size (around 1e-12), at the cost of a few more iterations.
2. **Matrix-multiply library dispatch.** faer picks instructions (FMA, vector width, blocking) by processor at run time. The repository already has a switch for it: the feature `portable-kernels` in `femcore`. It affects the **ported** solver library; the failing kernels may use the older stack, so it might not help them.
3. **Compiler and maths-library differences**, for example `exp`, `sin` and `sqrt` rounding in the last bit across platforms.

### Proposal (a test, then a rule)
* **Experiment first** (read-only, in my workspace): run the seven tests with `portable-kernels` and, separately, with a tighter eigensolver tolerance, and see which one closes the gap. This costs a build and two test runs and changes no code in the repository. I will not run it unless you say so, since you asked for proposals only.
* **The rule nyx states**: same machine and build, any thread count: identical bits. Different machine: agreement to the solver tolerance, which the software sets to **at least 100 times tighter than the mesh error it is meant to resolve**, and which is reported in the log.
* **Tests follow the rule**: stored-number tests compare at a tolerance justified by the solver tolerance (for example 1e-8 once the solvers converge to 1e-10), and every physics that has a closed form (slab effective index, cavity frequencies, analytic plate modes) is tested against that form, which depends on no machine at all.
* If the experiment shows that the old stack cannot be made machine-independent to 1e-8, the fallback is a **reference platform** for the strict tests (round 2, option C), with the closed-form tests running everywhere.
