Superconducting qubits & cQED: GDS/sheet geometry, capacitance & inductance (`dfs`), cavity EPR & Hamiltonian (`ems`), and TLS loss.

Geometry & meshing (`project_new(dim=3)` chips, `dim=2` CPW)
- GDSII: `nyx.import_gds(path="chip.gds", layers={"gnd":[1,0],"sig":[1,10]}, z=0.0, split_islands=True, pin_layer=[1,20], simplify_tol=2e-7, layer_bbox=True)` creates `Polygons` + `Sheet` features at `z` (`sheets`, `bbox`, `pins`).
- Stack (`add_sheet` + `add_overlay`): build `VacBot`, `Si`, `VacTop` (`add_box`), place sheets via `add_polygons`/`add_rect` + `add_sheet(name="c2", input="pad2", z=0.0)`, and list volumes in `add_overlay` (`c0`, `c2` boundaries, `c2:edge` rims).
- Mesh: `nyx.set_mesh(h=2e-4, gap_cells=2.5, z_coarsen=True)` refines coplanar gaps symmetrically and coarsens pure-tet slabs away from `z_sheet`.

Electrostatics, magnetostatics & convergence (`dfs`)
- `terminal`/`ground` vs box BCs: `add_condition(physics="dfs", on="c2", kind="terminal", params={"name":"c2"})` solves `N` terminals in 1 factorization. Ground `c0` with `kind="ground"` (or tag `c0` `"terminal"` and eliminate in `capacitance_matrix`). Box walls are Neumann unless grounded (`ground` on `zmin`/`zmax`).
- `electrostatic`: `add_physics(physics="dfs", study="electrostatic", order=2, study_params={"formulation":"volumetric_fem", "singular_edges":"quarter_point"})` shifts P2 mid-edge nodes at `sheet:edge` rims to `0.25`; `study_params={"formulation":"planar_bem", "richardson":True}` runs fast 2D layered MoM. `set_physics_params(physics="dfs", params={"solver":"auto", "amg_line_smooth":True, "tol":1e-10})` uses Block-CG + line AMG.
- `magnetostatic`: `study_params={"formulation":"auto"}` (`"cross_section"`, `"planar_stream"`, `"mixed_3d"`) with `london_sheet` (`L_k_sq` or `material`+`thickness`+`temperature`+`film_model`; `bcs_surface_impedance`), `current_terminal`, `flux_loop`, and `analyse(analysis="inductance_matrix", physics="dfs")`.
- Convergence: `nyx.convergence_check(physics="dfs", factor=1/1.5, scale_sizes=True, extrapolate="richardson_log_spd", tolerance=0.05)` scales `add_size` rules (`validated:"checked"`) and extrapolates `CapacitanceMatrixData`.

Post-solve analyses (`nyx.analyse`)
- `capacitance_matrix`: omit `params` for raw `maxwell_matrix_F`/`mutual_matrix_F`/`ec_matrix_Hz`; pass `params={"ground":"c0", "eliminate":"drop", "qubits":[{"name":"Q1","pads":["c2","c3"],"E_J_Hz":18e9}], "resonators":[{"name":"R1","terminal":"c1","f_r_Hz":6.8e9,"Z0_ohm":50.0,"resonator_type":"quarter_wave"}]}` (`"floating_schur"`, `Xmon`/`TransmonPocket`, `"half_wave"`/`"lumped"` `C_r_F`) for `maxwell_reduced_F`, `condition_number`, `C_Sigma_F`, `E_C_Hz`, `C_eff_F`, `beta`, `f01_Hz`, `anharmonicity_Hz`, `couplings` (`g_Hz`, `chi_Hz`).
- `epr` & `hamiltonian` (`ems` `cavity` + `lumped_jj` + `surface_impedance`): `epr` (`lj_override`) -> `p_matrix`, `phi_zpf_matrix`, `alpha_Hz`, `chi_matrix_Hz`, `kinetic_inductance_fraction`; `hamiltonian` (`basis="normal_fock"`, `fock_trunc=[7,5]`, `flux_sweep` or `"craig_bampton"`) -> `f01_Hz`, `f12_Hz`, `alpha_Hz`, `chi_matrix_Hz`.
- `tls_participation` (`dfs`/`ems`): `params={"solution":0, "t_ma":2e-9, "t_ms":2e-9, "t_sa":2e-9, "eps_ma":10.0, "eps_ms":11.45, "eps_sa":10.0, "tan_delta_ma":1e-3, "tan_delta_ms":1e-3, "tan_delta_sa":1e-3, "d":100e-9}` -> `p_ma`, `p_ms`, `p_sa`, `loss_tan_total`, `q_tls`.
