# Feature Request: Physics-Controlled Structured PML & Open-Boundary Geometry, Native Mixed-Cell Meshing/Solvers, and Curvilinear Coordinate Stretching (2D, 2D Axisymmetric, 3D)

- **Target ADR:** ADR 128 (next available after ADRs 123–127 on `main`; numbers `70..79` and `81` are reserved)
- **Status:** Proposed (Feature Request & First-Principles Physics Specification)
- **Reference Baseline:** COMSOL Learning Center #52791 (*Geometry and Mesh Setup for Modeling Regions of Infinite Extent*) and COMSOL 6.3 *Perfectly Matched Layers* Reference Manual
- **Companion Document:** [`docs/PML_IMPLEMENTATION_PLAN.md`](PML_IMPLEMENTATION_PLAN.md)

---

## Executive Summary

`nyx` currently models open and radiating domains using unstructured inline meshes inside Perfectly Matched Layers (`structured_pml = false` by default), coarsens PML elements by `pml_relax = 1.5`, evaluates coordinate stretching once per simplex centroid against per-domain axis-aligned bounding boxes (AABBs), forces solvers to split all cells into simplices (`Discretisation = SimplexMesh`), and supports only inward Cartesian box slicing via legacy `geom::Op::Pml`. Meanwhile, COMSOL Multiphysics' reference workflow prescribes mapped/swept meshes aligned with the stretching direction, separate corner/edge domains, at least 5–8 elements through the layer thickness, and Cartesian, Cylindrical, and Spherical topologies across 2D, 2D axisymmetric, and 3D.

Rather than blindly copying COMSOL's legacy CAD-slicing workflow—and rejecting compromises driven solely by implementation effort—this feature request derives the underlying wave physics, elastodynamics, and finite-element mathematics ("the why") and adopts the **highest-efficiency, fastest-solving, and most physically accurate design for `nyx`'s real workloads** (integrated photonics waveguides and WGM resonators, anisotropic crystal phononics/optomechanics, layered-wafer superconducting cQED chips, and electro-/magnetostatics):

1. **Physics-controlled structured meshing ON by default (`PmlMode::Auto`)**, supporting all **7 canonical geometries** (2D Cartesian, 2D Cylindrical, 2D Axisymmetric Spherical, 2D Axisymmetric Cylindrical, 3D Spherical, 3D Cartesian, and 3D Cylindrical) as conformal **`Tet` $\leftrightarrow$ `Prism` $\leftrightarrow$ `Hex`** (3D) and **`Tri` $\leftrightarrow$ `Quad`** (2D) meshes with **zero transition pyramids**.
2. **Native Mixed-Cell (`Quad` / `Prism` / `Hex`) FEM Kernels (`Option B`, fulfilling `ADR 5` & `ADR 32`):** Instead of splitting structured PML cells into 2× (`Quad` $\to$ 2 tris), 3× (`Prism` $\to$ 3 tets), or 6× (`Hex` $\to$ 6 tets) chiral simplices, physics kernels assemble natively on `Quad`, `Prism`, and `Hex` elements. This cuts PML elements by **2×–6×**, slashes PML DOFs by **37%–65%** (speeding up sparse `LBL^H`/`LU` factorization, which dominates `nyx` runtime), eliminates oblique Nédélec edges and split-diagonal doublet symmetry breaking, and **unlocks locking-free Rational PMLs and Rational Infinite Elements** across both $P_1$ and $P_2$.
3. **New Outward Boundary-Extruded Absorbing Layer (Deprecating Legacy Inward `Op::Pml`):** Users model **only the physical interior $\Omega_{\text{phys}}$** (`Box`, `Cylinder`, `Sphere`, `Rect`, `Circle`, or `add_chip_stack`) without manually padding outer boxes or creating dummy `window` primitives. The mesher sweeps structured `Quad`/`Prism`/`Hex` layers outward along analytical coordinate rays, automatically inheriting each touched base-facet material ($\partial_u\varepsilon \equiv 0$ guaranteed across multi-layer wafers and piercing waveguides) and extruding boundary edges of interior conductive `Op::Sheet`s into lateral PML sheets. Legacy inward `Op::Pml` is retained strictly for backward compatibility of existing `.nyx.json` projects and golden fixtures.
4. **Ray-invariance enforcement** forbidding physically invalid topologies (e.g., rejecting Spherical PML across layered horizontal wafers or piercing waveguides, and preferring wafer-normal Cylindrical PML for stratified chips and WGM disks).
5. **Quadrature-point coordinate-free curvilinear Jacobians $\mathbf{J}_s(\mathbf{x}_q)$**, supporting **both** bounded $\kappa$-Polynomial ABCS stretching (optimal for linear eigenfrequency studies and iterative `AMS`/`GMRES` solves) **and** locking-free Rational stretching (optimal for direct-solved open radiation/scattering in 5–6 layers), plus **dimensionless physics-calibrated scaling** (eliminating magic `strength = 3.0`) and lockstep curved concentric shells ($R_0 \dots R_{N_L}$) that prevent thin-layer sagitta inversion.
6. **Anisotropic elastic backward-wave protection in `nyx-sms`** (pre-solve Hellmann–Feynman Christoffel check, M-PML / adiabatic stiffness damping, signed-$Q$ instability detection, and clamped outer PML termination) plus a **3-tier open-boundary hierarchy in `nyx-dfs`** (`planar_bem`, `far_field` Robin / curl-Robin, and locking-free native `Quad`/`Prism`/`Hex` Rational & Exponential Infinite Elements for layered 3D chips).

---

## 1. First-Principles Physics & FEM Mathematics: Why Structured, Aligned Meshing Matters

### 1.1 The 7 Canonical Topologies for Regions of Infinite Extent

In COMSOL #52791, unbounded domains—whether **Perfectly Matched Layers** (frequency-domain waves), **Infinite Elements** (real-stretched diffusion/statics), or **Absorbing Layers** (time-domain)—share the exact same **7 canonical spatial decompositions**:

| # | Dimension | Coordinate System | Required Domain Decomposition | Separate Corners / Edges? | Reference Geometry Parameters |
|---|---|---|---|---|---|
| 1 | **2D Planar** (`x, z`) | **Cartesian** | 4 side bands (`±x`, `±z`) + 4 corner rectangles | **Yes** (4 corners) | Axis-aligned outer frame |
| 2 | **2D Planar** (`x, z`) | **Cylindrical** | 1 concentric annular ring ($R_{\text{in}} \le r \le R_{\text{out}}$) | **No** (0) | Center $(x_0, z_0)$ |
| 3 | **2D Axisymmetric** (`r >= 0, z`) | **Spherical** | 1 concentric annular semicircle meeting $r = 0$ ($R_{\text{in}} \le R \le R_{\text{out}}$) | **No** (0) | Center $(0, z_0)$ on axis $r = 0$ |
| 4 | **2D Axisymmetric** (`r >= 0, z`) | **Cylindrical** | 1 outer radial band (`+r`) + top/bottom axial bands (`±z`) + 2 corner rectangles (never on $r = 0$) | **Yes** (2 outer corners) | Symmetry axis $r = 0$ |
| 5 | **3D** (`x, y, z`) | **Spherical** | 1 concentric spherical shell ($R_{\text{in}} \le R \le R_{\text{out}}$) | **No** (0) | Center $(x_0, y_0, z_0)$ |
| 6 | **3D** (`x, y, z`) | **Cartesian** | 6 face slabs + 12 edge bars + 8 corner blocks (26 sub-domains) | **Yes** (12 edges + 8 corners) | Axis-aligned box |
| 7 | **3D** (`x, y, z`) | **Cylindrical** | 1 radial cylindrical shell + 2 axial endcaps + 2 corner annuli | **Yes** (2 corner annuli) | Center $\mathbf{x}_0$ + unit axis $\hat{\mathbf{a}}$ |

Physical sizing guidelines from COMSOL #52791 and the COMSOL 6.3 PML manual:
- Physical thickness $d$ is typically $\sim 1/10$ of the modeling space (kept regular to avoid thin distorted cells; in frequency-domain PMLs, coordinate scaling decouples mathematical absorption length from physical thickness).
- Scatterers and radiating sources should stay $\gtrsim \lambda/8$ away from the PML interface when evanescent components are strong (or real stretching $\kappa > 1$ / rational scaling must be active).
- Normalized PML depth `sDist<i>` must lie cleanly in $[0, 1]$ per active axis before solving, and radiated power into the PML ($P_{\text{pml}}$) should remain strictly positive for passive outgoing radiation.

### 1.2 Why Structured (Mapped 2D / Swept 3D) Meshes Beat Unstructured Simplices (5 Mechanisms)

Under complex coordinate stretching $u \mapsto \tilde{u}(u) = u_0 + \int_{u_0}^u s_u(\xi)\,d\xi$ with $s_u(u) = \kappa(u) - i\sigma_u(u)$ (under `nyx`'s $e^{+i\omega t}$ / $e^{-i\beta y}$ sign convention in `crates/femcore/src/pml.rs`), the continuous PDE is reflectionless across $u = u_0$ for all incidence angles because tangential field traces match identically across the interface. Unstructured triangles or tetrahedra (`PmlMode::Inline`) break exact analytic continuation through five distinct numerical mechanisms:

1. **Broken Discrete Bloch–Floquet Translation Invariance ($R_{\text{if}}$):**
   When PML cells are extruded node-for-node along the stretching ray $\hat{\mathbf{u}}$, every layer sheet $u = u_k$ ($k = 0, \dots, N_L$) has the **exact same tangential mesh topology and mass/curl/stiffness stencil** as the interior boundary $u = u_0$. An outgoing discrete Bloch mode $(k_u^h, \mathbf{k}_\perp)$ couples **only** to the identical tangential wavevector $\mathbf{k}_\perp$ inside the PML, reducing reflection to the 1D discrete WKB error $R_{\text{if}} = \mathcal{O}\!\left((k_u h_u)^{2p}\,\sigma_{\max}(h_u/d)^m\right)$. On an unstructured mesh, cell edges cut obliquely across iso-$u$ contours at varying angles, breaking tangential translation symmetry and scattering specular energy into parasitic high-$|\mathbf{k}'_\perp|$ grazing and evanescent modes where attenuation $\exp(-2\cos\theta\int\sigma\,du)$ collapses as $\cos\theta \to 0$.

2. **Eliminating Centroid Staircase Corrugation:**
   `nyx` currently evaluates $s_u$ once per simplex at its centroid (`crates/femcore/src/pml.rs:73-81`). On unstructured simplices, adjacent centroids sit at staggered normal depths $\Delta u \sim h_u$, turning smooth coordinate stretching into a **corrugated impedance sawtooth** of jump $\Delta s \sim |s'_u| h_u$ and an $\mathcal{O}(h_u)$ first-order reflection floor. **Both** aligned layer meshes *and* pointwise Gauss-quadrature evaluation of $\mathbf{J}_s(\mathbf{x}_q)$ are required to eliminate this floor.

3. **Anisotropic Metric Decoupling & Matrix Conditioning:**
   For Jacobian $\mathbf{J}_s = \partial\tilde{\mathbf{x}}/\partial\mathbf{x}$, the pulled-back stiffness/Maxwell metric tensor is:
   $$\mathbf{G}_s = (\det\mathbf{J}_s)\,\mathbf{J}_s^{-1}\mathbf{J}_s^{-T}.$$
   Under uniaxial stretching ($s_\perp = 1$), $G_{uu} = 1/s_u$ and $G_{\perp\perp} = s_u$, giving metric anisotropy $|G_{\perp\perp}/G_{uu}| = |s_u|^2$. On an aligned tensor-product layer element ($\mathbf{J}_K = \mathrm{diag}(h_u, \mathbf{J}_\perp)$), normal and tangential operators separate cleanly:
   $$\frac{1}{s_u h_u^2}\partial_{\hat{\xi}_u}\phi_i\,\partial_{\hat{\xi}_u}\phi_j + s_u\,\nabla_\perp\phi_i\cdot\nabla_\perp\phi_j.$$
   Thin normal layers ($h_u \ll h_\perp$) balance the stretched complex step ($|\Delta\tilde{u}| = |s_u| h_u \approx h_\perp$), whereas rotated unstructured simplices leak $\mathcal{O}(|s_u|)$ stiffness across axes, driving $\kappa(\mathbf{K}_e) \propto |s_u|^2 (h_{\max}/h_{\min})^2$ and stalling `AMS`, `AMG`, Krylov–Schur, and `LBL^H`.

4. **Preventing Discrete Complex-Dispersion Aliasing ($\ge 5\text{--}8$ Layers):**
   A PML terminated by PEC (`nyx-ems`) or Dirichlet (`nyx-sms` / `nyx-dfs`) suffers two competing errors: continuous outer-wall reflection
   $$R_{\text{wall}}(\theta) = \exp\!\left(-2 k_0 n \cos\theta \int_0^d \sigma_u(\xi)\,d\xi\right),$$
   and discrete complex-dispersion reflection per layer $h_u = d/N_L$ with complex phase-decay step $\Delta\tilde{\phi}_k = \frac{k_u d}{N_L}s_u(t_k)$. When $N_L \le 3$, $|\Delta\tilde{\phi}| > 2\text{--}4$, causing $P_1$/$P_2$ basis functions to alias $e^{-\alpha u}$ with Padé sign oscillations that reflect energy *before* reaching the outer wall. Coarsening PMLs via `WavesRequest::pml_relax = 1.5` (`crates/femcore/src/meshing/request.rs:40-56`) is physically backwards because $|\tilde{k}_u| = k_u |s_u| > k_u$ *shortens* the effective wavelength inside the layer.

5. **Separate Corner/Edge Domains & Why "Zero Transition Pyramids" Is a Major Advantage:**
   Across internal PML seams (e.g., $z = z_{\max}$ for $x > x_{\max}$), continuity of $\tilde{\mathbf{x}}(x, z) = (\tilde{x}(x), \tilde{z}(z))$ requires corner domains to stretch *all* intersecting axes simultaneously ($\mathbf{J}_s = \mathrm{diag}(s_x, s_z)$ in 2D; 2-axis in 3D box edges and cylindrical rim annuli; 3-axis in 3D box corners). An undivided L-shaped frame creates $45^\circ$ miter discontinuities ($\nabla\times\tilde{\mathbf{x}} \neq \mathbf{0}$) and cannot be mapped/swept without skewed cells.

   **Why zero transition pyramids is an architectural and numerical advantage (not a limitation):**
   - In 3D meshing, a **pyramid** (5 nodes: 1 quad base + 4 triangular sides meeting at an apex) is only ever needed as an awkward glue cell when someone tries to attach a **hexahedron's quad face** directly onto a **tetrahedron's triangle face**.
   - Pyramids are notorious trouble-makers in FEM: their basis functions cannot be pure polynomials on the reference pyramid—they require rational functions ($\propto \frac{\xi\eta}{1 - \zeta}$) with a singular Jacobian derivative at the apex, demand specialized Duffy/collapsed quadrature, complicate conforming $H(\mathrm{curl})$ Nédélec basis construction and inter-element orientation across `Tet`–`Pyramid`–`Hex` boundaries, and degrade conditioning when stretched inside a PML.
   - Because `nyx`'s physical interior mesh ends on **boundary triangles** (`Tet` faces), sweeping those triangles outward across 3D faces and spherical/cylindrical shells produces **triangular prisms (`Prism`)**:
     - The **triangle base** of each `Prism` matches the interior `Tet` triangle face **node-for-node**;
     - The **quad side walls** of adjacent `Prism`s match **quad-for-quad** with the **hexahedra (`Hex`)** swept across 3D box edges, 3D box corners, and 3D cylindrical rim annuli!
   - Thus **`Tet` $\leftrightarrow$ `Prism` $\leftrightarrow$ `Hex`** is **100% conforming everywhere with zero pyramids**—the exact same topology COMSOL generates when sweeping prisms on faces and hexes on edges/corners, with strictly polynomial/tensor-product finite elements throughout!

---

## 2. Critical First-Principles Audit: The Most Physical Design for `nyx`

Below are seven first-principles engineering verdicts establishing the **highest-accuracy, highest-efficiency, fastest-solving** architecture for `nyx`.

### Verdict 1 — Analytic-Continuation Ray-Invariance Theorem (Topology Selection & Guard)

By Teixeira–Chew / Collino–Monk analytic continuation, stretching $u \mapsto \tilde{u}(u)$ at fixed transverse coordinates $\mathbf{v}$ is reflectionless if and only if:
1. Material tensors are invariant along every stretching ray: $\partial_u \boldsymbol{\varepsilon}(u, \mathbf{v}) = 0$, $\partial_u \boldsymbol{\mu}(u, \mathbf{v}) = 0$, $\partial_u \mathbf{C}(u, \mathbf{v}) = 0$, $\partial_u \rho(u, \mathbf{v}) = 0$; and
2. Every material interface or waveguide boundary crossing $\partial\Omega_{\text{phys}}$ is a level set $\mathbf{v} = \text{const}$.

| Practice / Naive Assumption | Physical Failure Mode | Optimal `nyx` Design |
|---|---|---|
| Using **3D Spherical PML** everywhere (as in COMSOL antenna tutorials) | **Physically illegal on layered chips** (`Si`/`SiO2`/`Sapphire`/`TFLN`/`Air`) and **piercing waveguides**: horizontal plane $z = z_k \neq z_0$ satisfies $R\cos\vartheta = z_k$, so $\partial_R\varepsilon \neq 0$. Analytic continuation maps the plane to complex $\tilde{R}\cos\vartheta = z_k$ while the mesh keeps the jump at real $z = z_k$, causing an **$\mathcal{O}(|\Delta\varepsilon/\varepsilon|)$ reflection floor that never vanishes as $h \to 0$** (Oskooi et al., *Opt. Express* 16, 2008). Likewise, radial rays flare piercing waveguides by $\tilde{R}/R$. | **Use Cylindrical PML ($\hat{\mathbf{a}} = \hat{\mathbf{z}}$)** as default for layered wafers, WGM rings/disks, and anchor-loss resonators ($\partial_r\varepsilon(z) \equiv 0$, $\partial_r\mathbf{C}(z) \equiv 0$ identically across arbitrary stacks; endcaps & corner annuli lie in homogeneous air/substrate). Reserve **Spherical PML** for compact radiators/cavities/scatterers in **homogeneous** backgrounds ($6\times$ smaller volume, 0 corners). Require **Cartesian PML** for piercing waveguides/CPW lines, 2D `(x, z)` layered cross-sections, and Bloch cells. Add a pre-solve **ray-invariance guard** rejecting illegal Spherical/Cylindrical setups with a clear physics diagnostic. |

### Verdict 2 — Dual Stretching Laws Unlocked by Native Tensor Elements (Eigen vs. Driven/Radiation)

For an evanescent wave tail $e^{-\gamma u}$ ($\gamma = \sqrt{|\mathbf{k}_\perp|^2 - k^2} > 0$) under $s(t) = \kappa(t) - i\psi(t)$:
$$\exp(-\gamma\tilde{u}) = \exp\!\left(-\gamma\int_0^u \kappa(\xi)\,d\xi\right)\cdot\exp\!\left(+i\gamma\int_0^u \psi(\xi)\,d\xi\right).$$
Notice that **imaginary stretching $-i\psi$ produces zero evanescent attenuation**; only real coordinate stretching $\kappa(t) > 1$ (or rational infinite mapping) damps evanescent tails and suppresses reactive frequency pulling $\Delta f/f \propto e^{-2\gamma\bar{\kappa}d}$.

With native `Quad`/`Prism`/`Hex` elements (**Verdict 6**, which eliminate split-simplex tangential locking against singular poles), `nyx` can provide **both** physically optimal frequency-domain stretching laws—matched cleanly to each study class—while rejecting $\omega$-dependent CFS-PML in eigenproblems:

| Study Class in `nyx` | Why a Single Law Fails Everywhere | Optimal `nyx` Stretching Law |
|---|---|---|
| **Eigenfrequency studies** (`ems::cavity`, `ems::axisym`, `sms::eigenfrequency`, `sms::axisym`) & **Iterative `AMS`/`GMRES` solves** | Singular Rational mapping $[0, d) \to \infty$ collapses the continuous radiation spectrum onto a dense arc of spurious branch-cut modes crowding `f_search` in Krylov–Schur, and $\frac{1}{(1-t)^2}$ inflates $\kappa(\mathbf{K}_e)$ for iterative preconditioners. Explicit $\omega$-dependent CFS-PML ($\frac{\tilde{\sigma}}{\alpha+i\omega}$) breaks linear eigenproblem structure for zero benefit. | **Bounded $\kappa$-Polynomial ABCS Complex Coordinate Dilation** ($\omega$-independent):<br>$$s(t) = \big(1 + (\kappa_{\max}-1)t^m\big) - i\,\psi_{\max} t^m \quad (m = 2).$$<br>By Aguilar–Balslev–Combes–Simon (ABCS) spectral deformation theory, this keeps $K(s)u = \omega^2 M(s)u$ **strictly linear** with exact $\omega$-invariant complex resonance poles $\omega_n \in \mathbb{C}$, clean `AMS` conditioning, and strong evanescent damping via $\kappa > 1$ (8 layers at $P_2$, or 6 with geometric grading). |
| **Direct-solved Driven / Open-Radiation / Scattering** (`ems::driven` single-$f$ or direct `LBL^H` solves) | Grazing-angle and strong evanescent fields near subwavelength scatterers require high attenuation across thin physical layers ($5\text{--}6$ layers). | **COMSOL-style Rational Stretching** (now **locking-free** on native `Quad`/`Prism`/`Hex`!):<br>$$s(t) = \frac{1}{(1 - c_{\text{rat}} t)^2} - i\,\psi_{\text{rat}}\frac{t^m}{(1 - c_{\text{rat}} t)^2}$$<br>absorbs broadband propagating + grazing + evanescent waves to $R < 10^{-8}$ in just **5–6 native layers**, cutting driven-solve DOFs further. |
| **Broadband `ems::driven` PROM sweeps** (`ADR 106`) & **Waveguide / Bloch modes** (`ems::modes`, `ems::bloch`, `sms::modes`) | `driven.rs:309-424` relies on $\omega$-independent $(K, C, M)$ matrices for Hermite Tangential Galerkin PROM. | $\omega$-independent **$\kappa$-Polynomial** or fixed-band **Rational** calibrated across $[f_{\min}, f_{\max}]$ when $f_{\max}/f_{\min} \le 3$ (preserving 100% exact affine PROM; split into octave sub-bands when $f_{\max}/f_{\min} > 3$); waveguide/Bloch modes calibrate at exact study $\omega$. |

### Verdict 3 — Anisotropic Elastic Backward Waves & Outer Clamping in `nyx-sms`

In `nyx`'s $e^{+i\omega t}$ convention, an outgoing wave ($v_{g,u} = \langle\mathcal{P}_u\rangle/\langle\mathcal{E}\rangle > 0$) propagates as $\exp(-i k_u \tilde{u}) = e^{-i k_u \operatorname{Re}\tilde{u}} \exp(-k_u \int\sigma\,du)$. Whenever **$k_u v_{g,u} < 0$** (antiparallel phase and group projection along $\hat{\mathbf{u}}$), standard PML ($s_u = 1 - i\sigma$) **exponentially amplifies** the wave instead of damping it. This occurs via two mechanisms in `nyx`'s core materials:
1. **Bulk crystal slowness skewness** in rotated `Si [110]`, tilted R-plane `Sapphire` (`ADR 101`), `LiNbO3`, and quartz (since $k_u v_{g,u} + \mathbf{k}_\perp\cdot\mathbf{v}_{g,\perp} = \omega > 0$, normal incidence is stable while oblique qSV/qSH waves near slowness concavities can have $k_u v_{g,u} < 0$).
2. **Guided $S_1$ Lamb backward waves** in suspended membranes (even isotropic `Si3N4`/`SiO2`) where $\mathbf{k}\cdot\mathbf{v}_g < 0$.

- **Why CFS-PML fails here:** For $\alpha, \tilde{\sigma}, \omega > 0$, $\operatorname{Im}(s_{\text{CFS}}) < 0$ has the exact same sign and still amplifies backward waves.
- **Why COMSOL-like traction-free outer `nyx-sms` walls are suboptimal:** `nyx-sms` currently leaves outer PML walls free unless tagged `fixed`, creating zero-frequency rigid-body modes (`crates/nyx-sms/src/pml.rs:261`) and outer-wall Rayleigh surface waves.
- **Optimal `nyx-sms` Design:**
  1. Reuse `femcore::meshing::physics_size::wave_speed_range`'s 400-point Fibonacci Christoffel sweep to evaluate exact Hellmann–Feynman group velocity $v_{g,i} = \frac{1}{2\rho v_p} p_a p_b C_{aibj}\hat{q}_j$ at zero extra meshing cost and compute $p_{\min}(\hat{\mathbf{u}}) = \max_{\mathbf{v}_g\cdot\hat{\mathbf{u}} > 0}\frac{-(\hat{\mathbf{q}}\cdot\hat{\mathbf{u}})(\mathbf{v}_g\cdot\hat{\mathbf{u}})}{\mathbf{q}_\perp\cdot\mathbf{v}_{g,\perp}}$.
  2. Support `stabilization = "auto" | "off" | "mpml" | "adiabatic"`: `"mpml"` sets transverse damping $s_\perp = 1 - i p\sigma(t)$ with $p = 1.25\,p_{\min}$ when $p_{\min} > 0$ (curing Mechanism 1); `"adiabatic"` applies $\tilde{\mathbf{C}}(t) = (1 + i\eta_{\max}t^3)\mathbf{C}$ (Oskooi & Johnson 2008), which damps strictly by $\operatorname{sign}(v_{g,u})$ for *both* bulk and Lamb backward waves while keeping $K(\eta)u = \omega^2 M u$ linear.
  3. **Fix sign-masking bug (`crates/nyx-sms/src/pml.rs:232` and `lrb.rs:292`):** Replace `f.im.abs()` in $Q$ computation with signed handling that flags `f.im < -IM_FLOOR * f.re` as an unstable amplifying PML mode.
  4. **Clamp outer PML facets (`u = 0` Dirichlet) by default** for PML exterior boundaries in `nyx-sms` (opt-out available).

### Verdict 4 — Multi-Material Manifold Continuity & Dimensionless Physics-Calibrated Scaling

When multiple materials (e.g. `Air` $n=1.0$, `SiO2` $n=1.44$, `Si` $n=3.48$, or elastic P- and S-waves) cross the same PML slab/shell interface, **they MUST share the exact same coordinate map $\tilde{u}(u)$**. Giving adjacent sub-domains different `strength` or `order` tears $\tilde{\mathbf{x}}$ discontinuously across the material interface, violating tangential Nédélec/displacement continuity and creating artificial internal sources.

Therefore, calibrate every coupled PML group jointly using closed-form dimensionless formulas targeting **$R_0 = 10^{-6}$ ($-120\ \text{dB}$ power reflection)**, polynomial degree **$m = 2$**, and decoupled geometric layer growth ratio **$g = w_{N_L}/w_1 = 2.0$**:

- **Thickness & $\kappa$-Polynomial stretching (calibrated to fastest wave / lowest index $\lambda_{\max}$):**
  $$d_{\text{auto}} = 0.5\,\lambda_{\max}, \qquad \psi_{\max} = \operatorname{clamp}\!\left(\frac{(m+1)\ln(1/R_0)}{4\pi\,(d/\lambda_{\max})},\; 1.5,\; 25.0\right),$$
  $$\bar{\kappa} = \max\!\left(1.0,\; \frac{\ln(1/R_0)}{4\pi\,(d/\lambda_{\max})}\right), \qquad \kappa_{\max} = \operatorname{clamp}\!\big(1 + (m+1)(\bar{\kappa}-1),\; 1.0,\; \kappa_{\text{cap}}\big)$$
  *(for eigenfrequency cavity studies, cap $\kappa_{\max} \le 2.0$ to minimize branch-cut rotation into the spectral search disk; when $d \ge 1.1\lambda_{\max}$, natural evanescent decay beats $R_0$ and $\kappa_{\max} \to 1.0$ automatically).*
- **Layer count (calibrated to slowest wave / highest index $\lambda_{\min}$ and element order $p \in \{1, 2\}$):**
  - For $\kappa$-Polynomial:
    $$N_L(p) = \max\!\left(N_{L,\min}(p),\; \left\lceil G_p \cdot \frac{d}{\lambda_{\min}}\sqrt{\bar{\kappa}^2 + \left(\frac{\psi_{\max}}{m+1}\right)^2} \right\rceil\right), \qquad \begin{cases} p = 2: & N_{L,\min} = 6,\;\; G_2 = 5.0, \\ p = 1: & N_{L,\min} = 10,\; G_1 = 10.0. \end{cases}$$
  - For Rational profile (`PmlProfile::Rational` on direct-solved open problems): default $N_L(2) = 5\text{--}6$ ($N_L(1) = 8$) with geometric layer grading $g = 2.0$.
- **Per-physics bounds:**
  - `nyx-ems`: $\lambda_{\max} = c_0 / (f_{\text{low}}\,n_{\min,\text{PML}})$, $\lambda_{\min} = c_0 / (f_{\text{high}}\,n_{\max,\text{PML}})$.
  - `nyx-sms`: $\lambda_{\max} = c_{\max}/f_{\text{search}}$ (fastest quasi-P wave sets $d, \psi_{\max}$), $\lambda_{\min} = c_{\min}/(\text{headroom}\cdot f_{\text{search}})$ (slowest quasi-S wave sets $N_L(p) \approx 8\text{--}10$ at $P_2$ since $c_P/c_S \approx 1.7\text{--}2.2$).

### Verdict 5 — New Outward Boundary-Extruded Absorbing Layer vs. Deprecated Legacy `Op::Pml`

**Why legacy `Op::Pml` was awkward and why we replace its workflow (keeping it only for backward compatibility):**

- **What legacy `Op::Pml` did (inward CAD box slicer):**
  Legacy `Op::Pml` (`crates/geom/src/scene.rs:1229-1338`) could **not** add or extrude any new geometry outside a model. Instead, it forced users to:
  1. Manually pad every outer background region (`air_box`, `sapphire_box`, `substrate_box`) outward by `+thickness`;
  2. Draw an extra dummy `Rect`/`Box` named `"window"` solely to mark the unpadded boundary;
  3. Run `Op::Pml`, which boolean-sliced the already-oversized CAD regions into up to 26 separate sub-domains (`air.pml_xmin`, ..., `air.pml_xmax_ymax_zmax`); and
  4. Manually manage wildcard `"region*"` material assignments and duplicate `"pml"` conditions across all 26 sub-domains—while completely failing for spherical hosts (`stack.rs:261-263`), cylindrical shells, and outward `Op::Sheet` continuation.
  In short, it coupled physics wavelengths directly to CAD dimensions and pushed COMSOL's B-Rep CAD-slicing tax onto users.

- **What the new primary outward absorbing layer does (`mode = "outward"`, or `Op::AbsorbingLayer` / outward `Op::Pml`):**
  1. **Zero CAD padding or dummy `window`:** Users model **only the physical interior $\Omega_{\text{phys}}$** (`Box`, `Cylinder`, `Sphere`, `Rect`, `Circle`, or `add_chip_stack`) and declare an absorbing boundary/envelope on its exterior boundary (with `coord = "auto" | "cartesian" | "cylindrical" | "spherical"`, `thickness = "auto"` or $d$, `layers = "auto"` or $N_L$). Changing wavelength or frequency never edits CAD geometry.
  2. **Clean outward structured extrusion:** `femcore` meshes $\Omega_{\text{phys}}$ alone, then extrudes native structured `Quad` (2D) or `Prism` + `Hex` (3D) layers outward along analytical coordinate rays—naturally supporting **all 7 topologies**, including 3D Spheres and 3D Cylinders, without requiring 3D non-prismatic B-Rep booleans in `geom` (`ADR 3`).
  3. **Automatic material & sheet continuation:** Every extruded PML column automatically inherits the physical material ID of its touched base interior facet (guaranteeing $\partial_u\varepsilon \equiv 0$ across multi-layer chips and piercing waveguides), and boundary edges of interior conductive `Op::Sheet`s hitting $\partial\Omega_{\text{phys}}$ are extruded radially/normally into lateral PML sheets so CPW lines terminate without modal reflection.
  4. **Why keep legacy inward `Op::Pml` at all?** **Strictly for backward compatibility**: existing `.nyx.json` files, `.nys` scripts, and `crates/geom/golden/*.json` (`13_silicon_on_sapphire`, `17_disk_anchor_loss`, `chip_500te`, `chip_500tm`) must still parse and run without breaking user projects. Internally, even legacy inward `Op::Pml` skips its inline CAD slices during meshing and routes straight into the exact same outward `Quad`/`Prism`/`Hex` structured extruder!

### Verdict 6 — Option B (Native `Quad` / `Prism` / `Hex` FEM Kernels) Wins on Efficiency, Solver Speed & Quality

Earlier prototypes in `nyx` (`ADR 5`, `ADR 32`, `ADR 52`) split structured `Quad`/`Prism`/`Hex` layers into triangles/tetrahedra (`SimplexMesh::from_cells_split`) strictly to avoid adding `Quad`/`Prism`/`Hex` element routines to the physics crates (`Option A`). Evaluating Option A vs. Option B strictly on **efficiency, solver speed, memory, and physical/numerical quality**—without shrinking from kernel work—shows **Option B is decisively superior**:

| Criterion | Option A: Split `Quad`/`Prism`/`Hex` $\to$ `SimplexMesh` | **Option B (Chosen): Native `Tri`/`Quad` (2D) & `Tet`/`Prism`/`Hex` (3D) FEM Kernels** |
|---|---|---|
| **PML Element Count** | Inflated by **2×** (`Quad` $\to$ 2 tris), **3×** (`Prism` $\to$ 3 tets), **6×** (`Hex` $\to$ 6 tets). | **1× (2×–6× fewer PML elements)**, plus separable tensor-product Gauss quadrature. |
| **PML DOF Count & Sparse Factorization Speed (`LBL^H` / `LU`, 85–95% of `nyx` solve time)** | Split-tet $P_1$ Nédélec uses 7 DOFs/vertex; $P_2$ Nédélec uses **38 DOFs/vertex** in prism/hex stacks; $P_2$ Lagrange uses 8 DOFs/vertex; needs 8–10 layers for polynomial accuracy. | **Massive DOF & runtime reduction:**<br>• **$P_1$ Nédélec (`ems` ord 1, `dfs::magneto`):** **43% fewer** in `Prism` (4 vs 7 DOFs/v), **57% fewer** in `Hex` (3 vs 7 DOFs/v).<br>• **$P_2$ Lagrange (`sms`, `dfs` serendipity `8/15/20`-node):** **25% fewer** in `Quad`, **37.5% fewer** in `Prism`, **50% fewer** in `Hex`.<br>• **$P_2$ Nédélec (`ems` ord 2):** **37% fewer** in `Hex` (24 vs 38 DOFs/v) & fewer in `Prism` (36 vs 38).<br>• Combined with needing only **5–6 layers** instead of 8–10 for Rational/tensor accuracy $\implies$ **40%–65% fewer total PML DOFs**, accelerating `LBL^H` ($\mathcal{O}(N^{1.5\text{--}2})$) and reducing peak RAM substantially. |
| **Oblique Nédélec Edges** | Every split diagonal cuts obliquely across iso-$u$ contours, mixing exponential normal decay $e^{-\int k_u\sigma\,du}$ with tangential phase $e^{-i\mathbf{k}_\perp\cdot\mathbf{x}_\perp}$. | **Zero oblique edges.** Every PML edge is strictly tangential (constant $s_u$) or strictly normal ($\parallel \hat{\mathbf{u}}$) $\implies$ zero normal/tangential polarization cross-talk. |
| **Exact Symmetry & Degenerate Modes** | Dompierre lowest-vertex-ID coning + 2D diagonal cuts break $C_{4v}$, mirror, and rotational symmetry ($\Delta f/f_{\text{PML}} \neq 0$ for WGM/degenerate doublets). | **Exact geometric symmetry**—zero diagonal artifacts, **exact zero ($\Delta f/f \sim 10^{-14}$) PML-induced doublet splitting** in symmetric geometries. |
| **Anisotropic Locking Under High Stretch $\|s_u\|$ (Rational PML & Infinite Elements)** | Split simplices lack $(1-\xi)\eta$ ($P_1$) and $(1-\xi)\eta^2$ ($P_2$), causing **logarithmic tangential locking** and $\kappa(\mathbf{K}_e) \sim 10^6\text{--}10^{10}$ under singular $(1-t)^{-2}$ poles for any dipole/multipole or oblique wave. | **Zero tangential locking at both $P_1$ and $P_2$!** Tensor-product `Quad`/`Prism`/`Hex` spaces contain $(1-\xi)\eta^k$ identically, canceling $(1-\xi)^{-2}$ poles analytically and unlocking **true COMSOL Rational PML & Rational Infinite Elements** in 5–6 layers! |
| **Curved Concentric Shells (`ADR 82`)** | Diagonal splits across curved cylindrical/spherical layers need synthetic diagonal mid-nodes and distorted sub-tet Jacobians. | **Clean tensor-product `CurvedQuad`, `CurvedPrism`, `CurvedHex`** with orthogonal product Jacobian $\mathbf{J}_K = \mathbf{J}_\perp \otimes h_R$, keeping $\det\mathbf{J}_K/\det\mathbf{J}_{\text{straight}} \approx 1$ for arbitrarily thin curved shells. |
| **Fit with `nyx` Interfaces** | Avoids updating local element assembly. | `Physics::solve` (`crates/nyx-physics/src/physics.rs:101`) **already** passes `Mesh::D2(CellMesh<2>)` / `Mesh::D3(CellMesh<3>)` (`CellType::{Tri, Quad, Tet, Prism, Hex}`); `.fmsh` v2 **already** stores codes `0..4`; **and `Ams` / `AmsP2` (`crates/femcore/src/linalg/ams.rs:76, 429`) are already 100% element-shape agnostic** (taking pure algebraic `G`, `p1_to_p2`, `edge_vectors`)! Only local element assembly + Stage-3 `SimplexMesh` display projection for `view3d`/`h5io` change. |

**Verdict:** **Adopt Option B unconditionally.** Native mixed-cell (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) FEM assembly delivers unmatched speed, memory efficiency, symmetry preservation, and locking-free Rational PML/IE accuracy while fulfilling the original vision of `ADR 5` and `ADR 32`.

### Verdict 7 — 3-Tier Open-Boundary Hierarchy for `nyx-dfs` (& `nyx-sms::static`)

With native `Quad`/`Prism`/`Hex` unlocking locking-free tensor-product mappings, `nyx-dfs` and `nyx-sms::static` get a complete, physically optimal 3-tier hierarchy:

1. **Tier 1 — 2D planar zero-thickness conductors on layered wafers (`electrostatic`):** `planar_bem` (`ADR 115`) remains unconditionally optimal (analytical image-series Green's function to $r = \infty$, zero exterior volume mesh).
2. **Tier 2 — Homogeneous 2D/3D exterior (`electrostatic`, extend to `heat_stationary` and `magnetostatic`):** `ADR 120` `far_field` Robin / Order-2 Bayliss–Turkel BCs use **0 extra DOFs**, achieve $\mathcal{O}((a/R)^6)$ accuracy, keep SPD conditioning, and preserve Prager–Synge $H(\mathrm{div})$ lower bounds. Extend `femcore::farfield::assemble_farfield_robin` to `heat_stationary` ($\varepsilon_0\varepsilon_r \to \kappa$ around $T_\infty$) and add 1st-order curl-Robin $(\mu_0\mu_r)^{-1}(\nabla\times\mathbf{A})\times\hat{\mathbf{n}} + \frac{1}{\mu_0 r}\mathbf{A}_t = \mathbf{0}$ to `nyx-dfs::magneto`.
3. **Tier 3 — Layered/anisotropic 3D wafers with thick conductors/trenches/airbridges & `nyx-sms::static`:** Reuse the outward structured extruder with **5–6 native `Quad`/`Prism`/`Hex` layers** supporting:
   - **True Rational Infinite Elements** ($s_u(t) = \frac{d_{\text{pole}}(d_{\text{pole}}+d)}{(d_{\text{pole}}+d - t d)^2}$ mapping $[R_{\text{in}}, R_{\text{out}}) \to [R_{\text{in}}, \infty)$—now **locking-free** on native `Quad`/`Prism`/`Hex` for direct solves and `P2`!), and
   - **Bounded Exponential Real Stretching + Order-1 `far_field` Robin** ($\tilde{u}(t) = u_0 + L_{\text{char}}\frac{e^{\beta t}-1}{\beta}$, $\chi = e^\beta \in [20, 50]$—optional gentler conditioning mode for iterative `AMG`+`PCG`).
   - **Transient Heat Guard (`heat_transient`):** Static Infinite Elements remain **physically invalid** for parabolic diffusion ($\rho c_p \det\mathbf{J}_{\text{IE}}$ distorts thermal heat capacity across timescales). Reject `infinite_element` on `heat_transient` and guide users to unstretched ($s_u \equiv 1$) geometrically graded shells of thickness $d = \sqrt{4 D_{\max} t_{\text{end}}\ln(1/R_0)}$.

---

## 3. Code-Verified Gap Audit of Current `nyx`

| Area | Current Behavior in `nyx` | Code Reference | Gap to Close |
|---|---|---|---|
| **Solver mesh contract** | `Discretisation = SimplexMesh`; `Mesh::discretise()` immediately splits all `Quad`/`Prism`/`Hex` cells into simplices (`Option A`). | `crates/nyx-physics/src/lib.rs:47-51`, `crates/femcore/src/mesh/exchange.rs:95-100`, `crates/femcore/src/mesh/mod.rs:520-703` | Introduce native `MixedDiscretisation` (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) + shared `fem_mixed`/`edge_mixed` kernels (`Option B`), keeping Stage-3 projection for `h5io`/`view3d`. |
| **Default PML meshing mode** | `structured_pml` defaults to `false`; all 17 shipped examples run `PmlMode::Inline` and coarsen PMLs via `pml_relax = 1.5`. | `crates/femcore/src/meshing/request.rs:40,74`, `crates/nyx-engine/src/project.rs:49`, `crates/nyx-engine/src/meshing.rs:355-363` | Default to `PmlMode::Auto` (structured ON whenever PML/IE exists); remove `pml_relax = 1.5` coarsening. |
| **Supported topologies** | Legacy `Op::Pml` computes only the Cartesian AABB of `window`; `pml_box::build_pml` rejects non-box faces; `geom::stack` rejects standalone 3D `Sphere` hosts. | `crates/geom/src/scene.rs:1229-1338`, `crates/femcore/src/meshing/pml_box.rs:276-278`, `crates/geom/src/stack.rs:261-263` | Add 2D Cylindrical, 2D-axi Spherical, 3D Cylindrical, and 3D Spherical structured PMLs + ray-invariance guard. |
| **Geometry workflow & UI** | `Op::Pml` cuts inward only after manual domain padding + dummy `window`; GUI Physics tab has no button for `add_pml_conditions`. | `crates/geom/src/scene.rs:1319-1335`, `crates/app/src/tabs.rs:918-1004` | Add outward boundary-extruded PML (no padding/`window`, inheriting materials & `Op::Sheet` edges); keep legacy `Op::Pml` back-compat only; add GUI Physics button & envelope preview. |
| **Stretching evaluation & profiles** | `cell_stretches` evaluates $s$ at simplex centroids against per-domain AABBs; polynomial only with hardcoded `strength = 3.0`. | `crates/femcore/src/pml.rs:13-97`, `crates/femcore/src/props.rs:29-62`, `crates/nyx-engine/src/pml.rs:14-49` | Evaluate coordinate-free $\mathbf{J}_s(\mathbf{x}_q)$ at Gauss points from explicit region specs; add $\kappa$-Polynomial + Rational stretching + dimensionless auto-scaling. |
| **Curvilinear EM/SMS kernels & $r=0$** | `nyx-ems` stores `mu_inv` as diagonal `[C; 3]`; `nyx-sms::pml` uses `grad / s[k]`; `nyx-ems::axisym` rejects $r \le 0$ and lacks $r = 0$ axis constraints. | `crates/nyx-ems/src/cavity.rs:255`, `crates/nyx-sms/src/pml.rs:126-129`, `crates/nyx-ems/src/axisym.rs:273-275,316-340` | Widen `mu_inv` to symmetric $3\times 3$ `[[C; 3]; 3]`; apply $\mathbf{J}_s^{-T}\nabla\phi$ in `nyx-sms`; allow $r = 0$ with $m$-dependent axis conditions in `nyx-ems::axisym`. |
| **Backward waves & outer `sms` BC** | `nyx-sms` lacks backward-wave detection/stabilization, hides `f.im < 0` via `f.im.abs()`, and defaults outer PML walls to traction-free. | `crates/nyx-sms/src/pml.rs:232,261`, `crates/nyx-sms/src/lrb.rs:292` | Add pre-solve Christoffel check, `mpml`/`adiabatic` stabilization, signed-$Q$ instability alerts, and default outer clamping (`u = 0`). |
| **Layer count, grading bug & convergence** | `Op::Pml::grading` is overloaded as stretch `order` while `PmlSpec::grading` is hardcoded to `1.0`; `convergence::refined()` skips `Op::Pml::layers`. | `crates/femcore/src/meshing/geomdrive.rs:150`, `crates/femcore/src/meshing/pml_geom.rs:99`, `crates/femcore/src/meshing/pml_geom3.rs:106`, `crates/nyx-engine/src/convergence.rs:58-118` | Separate `layer_grading` (default `2.0`) from stretch `order` (default `2.0`); warn on `< 5` layers; scale PML `layers` in `refined()`. |
| **2D wasted work & `sweep_model_3d`** | 2D `structured_pml` meshes inline PML triangles then discards them via `physical_part`; 3D `sweep_model_3d` ignores `structured_pml`. | `crates/femcore/src/meshing/pml_geom.rs:143-149`, `crates/nyx-engine/src/meshing.rs:403-441` | Skip PML regions in 2D physical CDT/Ruppert (`mesh_geom_2d`) and wire structured PML into `sweep_model_3d`. |
| **DFS open boundaries** | `nyx-dfs` supports `open`/`far_field` only in `electrostatic`; nothing for `heat_*`, `magnetostatic`, or `sms::static`. | `crates/femcore/src/farfield.rs:404-581`, `crates/nyx-dfs/src/magneto.rs:90`, `crates/nyx-dfs/src/lib.rs:177-315` | Extend `far_field` to `heat_stationary`, add curl-Robin to `magneto`, and add native Rational/Exponential IE for layered 3D chips / `sms::static`. |

---

## 4. Formal Functional Requirements

### FR-1 — Complete 7-Topology Coverage & Ray-Invariance Guard
- Support structured geometry specification, native mixed-cell meshing (`Tri`/`Quad` and `Tet`/`Prism`/`Hex` with **zero pyramids**), and coordinate-free stretching for all 7 canonical infinite-extent configurations (2D Cartesian, 2D Cylindrical, 2D Axisymmetric Spherical, 2D Axisymmetric Cylindrical, 3D Spherical, 3D Cartesian, 3D Cylindrical).
- Auto-detect topology from the outer physical boundary primitive (`Rect`/`Box` $\to$ Cartesian, `Circle` $\to$ Cylindrical or 2D-axi Spherical semicircle at $x = 0$, `Cylinder` $\to$ 3D Cylindrical, `Sphere` $\to$ 3D Spherical) unless explicitly overridden.
- Enforce a pre-solve **ray-invariance validation guard** checking material constancy along stretching rays on $\partial\Omega_{\text{phys}}$ (e.g., rejecting Spherical PML when off-center horizontal layers or piercing waveguides intersect the boundary, recommending Cylindrical or Cartesian PML).

### FR-2 — Primary Outward Absorbing Layer & Backward-Compatible Legacy `Op::Pml`
- Introduce the first-class **outward boundary-extruded absorbing layer** (`mode = "outward"` / `AbsorbingLayer`, default for new workflows): users define only $\Omega_{\text{phys}}$ (no outer region padding or dummy `window`), and `femcore` sweeps `Quad`/`Prism`/`Hex` layers outward, **automatically inheriting the touched base facet's material** and **continuing intersecting `Op::Sheet` edges** into lateral PML sheets.
- Retain legacy inward `Op::Pml` (`mode = "inward"`) strictly for backward compatibility with existing `.nyx.json` files and goldens, routing it internally into the same outward structured builder after skipping its CAD-sliced PML pieces.

### FR-3 — Physics-Controlled Automatic Structured Meshing & Dimensionless Calibration
- Replace `structured_pml: bool` default `false` with `pml_mode: Auto | Structured | Inline` (default `Auto`, backward-compatible with legacy `structured_pml`), automatically building native mapped/swept `Quad`/`Prism`/`Hex` layers whenever PML or IE domains exist.
- Decouple `layer_grading` (geometric layer growth ratio, default `2.0`) from stretch polynomial `order` (default `2.0`).
- Remove `pml_relax = 1.5` coarsening inside structured PMLs; implement multi-material continuous dimensionless auto-scaling for $d_{\text{auto}}$, $\psi_{\max}$, $\kappa_{\max}$, and $N_L(p)$ ($\ge 5\text{--}6$ Rational / $\ge 6\text{--}8$ $\kappa$-Polynomial at $P_2$), emitting a `MeshReport` warning whenever any PML specifies `< 5` layers.
- Skip inline PML domains during 2D CDT/Ruppert physical meshing, add structured PML support to `sweep_model_3d`, and automatically scale absorbing-layer `layers` inside `nyx_engine::convergence::refined`.

### FR-4 — Native Mixed-Cell (`Quad` / `Prism` / `Hex`) FEM Architecture (`Option B`) & Lockstep Curving
- Amend `ADR 4` while fulfilling `ADR 5` and `ADR 32`: introduce `MixedDiscretisation` (`Tri`+`Quad` in 2D, `Tet`+`Prism`+`Hex` in 3D) with canonical $D_3$/$D_4$ face-orientation matching, shared tensor-product Lagrange ($Q_1$, serendipity $Q_2'$, full $Q_2$) and conforming Nédélec ($N_1$, $N_2$) element spaces, and algebraic `gradient_p1` / `gradient_p2` / `p1_to_p2` operators wired directly into existing `Ams` and `AmsP2`.
- Keep `CellMesh<D>` uncut through `MeshedModel` into `Physics::solve(&Mesh, ...)`; implement Stage-3 projection onto `SimplexMesh`-indexed `CellField`s so `h5io.rs`, `xdmf.rs`, `view3d.rs`, and `reduce.rs` work seamlessly without visualization breakage.
- Curve every concentric cylindrical/spherical shell $R_0, \dots, R_{N_L}$ coherently in lockstep inside `pml_cyl`/`pml_sphere` via `CurvedQuad` / `CurvedPrism` / `CurvedHex` (eliminating thin-layer sagitta inversion).
- Preserve fixed `SimplexMesh::from_cells_split_with` (positive tet orientation + `edge_nodes` retention) for legacy `Inline` / fallback paths.

### FR-5 — Quadrature-Point Curvilinear Stretching, Dual Profiles & Kernel Support
- Replace centroid + AABB stretching with pointwise Gauss-quadrature evaluation of `pml_jacobian(spec, x_q, dim, axisym)` using explicit region geometry metadata.
- Support both $\kappa$-Polynomial ABCS (default for eigenfrequency, iterative `AMS`, waveguides, broadband PROM) and locking-free Rational stretching (for direct-solved driven/scattering problems in 5–6 layers).
- Upgrade `mu_inv` from `[C; 3]` to symmetric $3\times 3$ `[[C; 3]; 3]` across `nyx-ems`, apply full $\mathbf{J}_s(\mathbf{x}_q)^{-T}\nabla\phi$ in `nyx-sms`, and allow meshes touching $r = 0$ in `nyx-ems::axisym` with $m$-dependent axis regularity ($m=0$: $E_\theta=0$; $|m|=1$: $E_z=0$; $|m|\ge 2$: both $0$).

### FR-6 — `nyx-sms` Elastic Backward-Wave Protection & Outer Clamping
- Run a zero-cost pre-solve Hellmann–Feynman Christoffel check across PML materials, logging $p_{\min}(\hat{\mathbf{u}})$ whenever antiparallel waves ($k_u v_{g,u} < 0$) exist.
- Support `stabilization = "auto" | "off" | "mpml" | "adiabatic"`, fix the `f.im.abs()` sign-masking bug so amplifying `f.im < 0` modes are flagged as unstable, and default outer PML exterior walls in `nyx-sms` to clamped (`u = 0` Dirichlet).

### FR-7 — 3-Tier Open-Boundary Completion for `nyx-dfs` (& `nyx-sms::static`)
- Extend `femcore::farfield::assemble_farfield_robin` to `nyx-dfs::heat_stationary`.
- Add 1st-order asymptotic curl-Robin BC to `nyx-dfs::magneto`.
- Add outward-extruded native `Quad`/`Prism`/`Hex` **Rational Infinite Elements** (locking-free in 5–6 layers for direct solves) and optional **Bounded Exponential + Order-1 Robin** (`infinite_element`) for layered/anisotropic 3D `nyx-dfs::{electro, magneto, heat_stationary}` and `nyx-sms::static`; explicitly reject `infinite_element` on `heat_transient` with guidance to use unstretched diffusion-length layers.

### FR-8 — Pre-Solve & Post-Solve Diagnostics
- Expose normalized distance `sDist ∈ [0, 1]` as a pre-solve verification scalar/field, report `pml_layers`, `min_pml_scaled_jacobian`, `element_counts_by_type`, and `< 5` layer warnings in `MeshReport`, and expose radiated PML power $P_{\text{pml}}$ alongside `layer fraction N`.

### FR-9 — Front-End Parity (Engine, Script, MCP, GUI, Examples)
- Update `add_absorbing_layer` / `add_pml`, `add_pml_conditions` (supporting auto-scaling and curvilinear params), `set_mesh`, `.nys` prelude/stubs, `crates/nyx-engine/cards/geometry.md` (preserving the `< 3600` byte budget test in `guide.rs:31`), and GUI controls (coordinate/mode/profile pickers, `Auto / Structured / Inline` mesh selector, translucent outward preview envelope, and one-click Physics-tab **"Add PML Conditions from Geometry"** button).
- Keep existing examples working cleanly via legacy compatibility, modernize `09_pml_slab`, `11_anchor_loss`, `13_silicon_on_sapphire`, and `17_disk_anchor_loss` to outward structured `Quad`/`Prism`/`Hex` setups, and ship canonical 2D Cylindrical, 2D Axisymmetric Spherical, and 3D Cylindrical/Spherical PML examples.

---

## 5. Non-Functional Requirements & Quantitative Acceptance Criteria

1. **Determinism:** Bit-identical mesh topology, coordinates, and solve results across `NYX_THREADS=1` and `NYX_THREADS=4`.
2. **Architecture, Layering & Contracts:** Pass `crates/nyx-engine/tests/layering.rs`, `contracts` tests (`NYX_WRITE_CONTRACT=1` + `VERSIONS.md` bump in the same commit), selective formatting (`ADR 65`), and `bash docs/architecture/check_code_to_adr.sh`. Existing golden tests remain untouched (pinning legacy unstructured golden fixtures to `pml_mode: "inline"` only where required).
3. **Numerical & Efficiency Benchmarks (all `--release`):**
   - **De Rham & `AMS` Exactness on Mixed Meshes:** $\|K G_1\|_\infty < 10^{-13}$ and $\|K G_2\|_\infty < 10^{-13}$ across `Tet + Prism + Hex` interfaces; `Ams` and `AmsP2` iteration counts match or beat pure-tet baselines.
   - **DOF & Solver Speedup vs. Split-Tet:** $\ge 35\%$ reduction in PML DOFs and measurable `LBL^H` factorization speedup on 3D cavity/driven PML problems relative to `Option A` split tets.
   - **2D Planar Cartesian & Cylindrical PML:** Reflection $R < 10^{-6}$ ($\kappa$-Polynomial) and $R < 10^{-8}$ (Rational, 5–6 layers) with zero degenerate doublet splitting ($\Delta f / f < 10^{-12}$) on $C_{4v}$ / circular symmetric tests.
   - **2D Axisymmetric Spherical PML ($r = 0$ touching) & Cylindrical PML:** Mie / dipole radiation $Q$ and frequency for $m = 0, 1, 2$ within $< 0.5\%$ of analytical values.
   - **3D Cartesian, Cylindrical (layered wafer), & Spherical Curved-Shell PML:** Zero inverted elements (`min scaled Jacobian > 0.3` even at $h_R / R_{\text{in}} = 0.01$); $P_2$ cavity/driven $R < 10^{-6}\text{--}10^{-8}$; ray-invariance guard rejects Spherical PML on multi-layer wafers while Cylindrical PML converges monotonically under `convergence_check`.
   - **`nyx-sms` Backward Waves & Clamping:** Tilted `Sapphire`/`Si [110]` and isotropic $S_1$ Lamb membrane exhibit zero spurious `f.im < 0` amplifying modes under `"auto"`/`"mpml"`/`"adiabatic"` and zero rigid-body modes when outer-clamped.
   - **`nyx-dfs` Open Boundaries:** Isolated 3D sphere capacitance error $< 0.1\%$ at $R_{\text{in}} = 2R$; `magneto` curl-Robin reduces box-wall inductance bias by $> 10\times$ vs. PEC walls; 5-layer native `Prism`/`Hex` Rational Infinite Elements solve monopole, dipole, and quadrupole layered-chip test cases to $< 0.1\%$ without tangential locking.
