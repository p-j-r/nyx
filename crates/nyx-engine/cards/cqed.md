Superconducting qubits & cQED: GDS/sheet geometry, capacitance & inductance (`dfs`), cavity/driven EPR, network & Hamiltonian (`ems`), and TLS loss.

Geometry & meshing (`project_new(dim=3)` chips, `dim=2` CPW)
- GDSII: `nyx.import_gds(path="chip.gds", layers={"gnd":[1,0],"sig":[1,10]}, z=0.0, split_islands=True, pin_layer=[1,20], simplify_tol=2e-7, layer_bbox=True)` creates `Polygons` + `Sheet` features at `z` (`sheets`, `bbox`, `pins`).
- Stack (`add_sheet` + `add_overlay`): build `VacBot`, `Si`, `VacTop` (`add_box`), place sheets via `add_polygons`/`add_rect` + `add_sheet(name="c2", input="pad2", z=0.0)`, and list volumes in `add_overlay` (`c0`, `c2` boundaries, `c2:edge` rims).
- Mesh: `nyx.set_mesh(h=2e-4, gap_cells=2.5, z_coarsen=True, max_cells=500000, dz_min=5e-6, max_z_slabs=16, z_growth_ratio=1.35, max_steiner_multiplier=20, strict_mesh=False)` controls gap clamping, z-slab growth/caps, pre-allocation aborts, and Ruppert diagnostics.
- Execution: `solve_run`/`queue_submit`/`sweep_run` accept `save_fields=False`, `post_analyses=[{"analysis":"hamiltonian"}]` (runs in-worker before stripping fields), and `memory_wait_timeout_s`.

Electrostatics, magnetostatics & convergence (`dfs`)
- `terminal`/`ground`/`open` BCs: `add_condition(physics="dfs", on="c2", kind="terminal", params={"name":"c2"})` solves `N` terminals in 1 factorization. Use `kind="open"` (optional `origin`) on outer box faces for $1/r$ Robin far-field termination.
- `electrostatic`: `study_params={"formulation":"volumetric_fem", "singular_edges":"quarter_point"}` or `{"formulation":"planar_bem", "richardson":True}`; `set_physics_params(physics="dfs", params={"solver":"auto", "amg_line_smooth":True, "tol":1e-10})`.
- `magnetostatic`: `study_params={"formulation":"auto"}` with `london_sheet`, `current_terminal`, `flux_loop`, and `analyse(analysis="inductance_matrix", physics="dfs")`.
- Convergence: `nyx.convergence_check(physics="dfs", factor=1/1.5, scale_sizes=True, quantities=[...], extrapolate="richardson_log_spd", tolerance=0.05)` scales `add_size` and `dz_min` linearly without squaring `adapt_scale`.

Cavity & driven EM (`ems`)
- `cavity`: supports `lumped_jj` (`squid_asymmetry`), `lumped_port` (`z0_ohm`, `mode="single"|"cpw_even"`, `complex_damping`), and `surface_impedance` (`lossless=True` uses real eigensolve; omitted `thickness`/`film_model` auto-load from material TOML).
- `driven`: supports `freqs_hz`, `seed_freqs_hz`, `max_hdm_solves`, `strict_prom_tol`, `waveport` (`n_eff_guess`, `z0_ohm`), `lumped_port`, and complex dielectric `loss_tangent`/`eps_r_im`.

Post-solve analyses (`nyx.analyse`)
- `capacitance_matrix`: omit `params` for raw Maxwell/mutual/$E_C$ matrices; pass `ground`, `eliminate`, `qubits`, `resonators`, and `metal_thickness_m`/`edge_cutoff_m` for 3D finite-thickness correction.
- `epr` & `hamiltonian`: `hamiltonian` returns `f01_Hz`, `f12_Hz`, `alpha_Hz`, `chi_matrix_Hz`, `f_cav_Hz`, `g_matrix_Hz`, `kappa_cav_Hz`, `t1_purcell_s`; `params` accept `basis`, `fock_trunc`, `sweep_port`, `flux_sweep`, `squid_asymmetry`, `squid_flux_sweep`, `ej_hz_sweep`.
- `network` (`CircuitRom` or `DrivenPortData`): returns damped `qep_modes`, ROM `Z`/`Y`/`S`, Schur `bbq_modes`, Taubin `s21_fit` (`"hanger"|"transmission"|"reflection"`), and `touchstone_snp`.
- `tls_participation` (`dfs`/`ems`): auto-loads `[interfaces.{ma,ms,sa}]` from material TOML (`silicon_4k`, `aluminium_sc_emft`); `params` override interfaces and accept `trench_depth_m`.
