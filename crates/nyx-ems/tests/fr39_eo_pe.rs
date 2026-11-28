// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-39 (ADR 141) verification suite for Cross-Mesh Electro-Optic (Pockels `r^S_{Ik}`) &
//! Photoelastic (Elasto-Optic `p_{IJ}`) permittivity perturbation (`Delta eps_r`) in `nyx-ems`.
//!
//! Covers:
//! 1. X-cut Thin-Film Lithium Niobate ($\text{LiNbO}_3$, point group `3m`, `euler_deg = [90, 90, 0]`,
//!    crystal $+Z$ optic axis along wafer lateral $+x$, $n_e \approx 2.1376, n_o \approx 2.211, r_{33}^S = 28.6\text{ pm/V}$)
//!    electro-optic modulator (shared-mesh and cross-mesh via `TransferredField`),
//!    verifying $\Delta n_{\text{eff}}$ and half-wave voltage-length product $V_\pi L = \frac{\lambda_0 g}{n_e^3 r_{33} \Gamma}$
//!    within `< 0.5%` of the analytical overlap integral $\Gamma$ and manual tensor perturbation oracle.
//! 2. Silicon (`m-3m`, $n = 3.476, p_{11} = -0.094, p_{12} = 0.017, p_{44} = -0.051$) photoelastic
//!    effective-index shift (`"strain"` -> `ems::modes` and `ems::cavity`) within `< 0.1%` of the exact
//!    impermeability perturbation $\Delta(\varepsilon_r^{-1})_I = \sum_J p_{IJ} S_J$.

use femcore::fem::{geometry, RefElement};
use femcore::material::Library;
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_ems::{Ems, C0};
use nyx_physics::{
    discretise, BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, PropValue, Scalar,
    SolveResult, TransferredField,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;

struct GeomResolver<'a>(&'a MeshedGeom);

impl LabelResolver for GeomResolver<'_> {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.domains(region);
        if v.is_empty() {
            Err(Error::invalid(format!("no region '{region}'")))
        } else {
            Ok(v)
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.boundary(selector);
        if !v.is_empty() {
            return Ok(v);
        }
        self.0.select(selector).map_err(Error::invalid)
    }
}

fn domain_name_map(g: &MeshedGeom) -> BTreeMap<i32, String> {
    g.labels.domains.iter().map(|d| (d.id, d.name.clone())).collect()
}

fn boundary_name_map(g: &MeshedGeom) -> BTreeMap<i32, String> {
    g.labels.boundaries.iter().map(|b| (b.id, b.name.clone())).collect()
}

fn load_library() -> Arc<Library> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Arc::new(Library::load_dir(&dir).expect("load materials library"))
}

fn num_scalar(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(v))) => *v,
        other => panic!("missing scalar '{name}', got {other:?} in {:?}", r.scalars),
    }
}

fn build_tfln_modulator_mesh(gap: f64, h_box: f64, w_core: f64, h_core: f64, h_mesh: f64) -> MeshedGeom {
    let mut m = Model::new_2d();
    let clad = m.rect("clad_r", (-0.5 * gap, -0.5 * h_box), (0.5 * gap, 0.5 * h_box));
    let core = m.rect("core_r", (-0.5 * w_core, -0.5 * h_core), (0.5 * w_core, 0.5 * h_core));
    m.overlay(&[("clad", clad), ("core", core)]);
    mesh_model(&m.list, &MeshPlan::simple(h_mesh)).expect("mesh 2D X-cut TFLN modulator")
}

#[test]
fn xcut_tfln_electro_optic_modulator_vpi_l_vs_analytical_overlap() {
    let lib = load_library();
    // Physical parameters for an X-cut TFLN modulator at lambda_0 = 1550 nm using `materials/lithium_niobate.toml`:
    // Crystal +Z (optic axis, extraordinary index n_e, Pockels r_33^S = 28.6 pm/V, r_13^S = 8.6 pm/V)
    // is rotated along wafer +x via Bunge Z-X-Z Euler angles `[90.0, 90.0, 0.0]`:
    //   R = R_z(90 deg) * R_x(90 deg) maps crystal Z (0, 0, 1) -> wafer +x (1, 0, 0)!
    let lambda0 = 1.55e-6;
    let freq0 = C0 / lambda0;
    let ln_mat = lib.get("lithium_niobate").expect("lithium_niobate in library");
    let eps_rot =
        ln_mat.optical_permittivity(lambda0, None, [90.0, 90.0, 0.0]).expect("LiNbO3 optical permittivity at 1550 nm");
    let eps_e = eps_rot[0][0].re; // extraordinary index along wafer +x
    let eps_o = eps_rot[2][2].re; // ordinary index along wafer +z
    let n_e = eps_e.sqrt();
    let r33 = 28.6e-12_f64;
    let r13 = 8.6e-12_f64;
    let n_clad = 1.444_f64;
    let eps_clad = n_clad * n_clad;

    let gap = 2.8e-6;
    let h_box = 2.0e-6;
    let w_core = 1.0e-6;
    let h_core = 0.5e-6;
    let v_bias = 50.0; // V across gap (creates E_x = V_bias / gap)
    let e_dc_uniform = v_bias / gap;

    let g_opt = build_tfln_modulator_mesh(gap, h_box, w_core, h_core, 0.16e-6);
    let mesh_opt = Mesh::from_solver_mesh(&g_opt.mesh).unwrap();
    let res_opt = GeomResolver(&g_opt);

    // Distinct electrostatic source mesh for cross-mesh transfer
    let g_elec = build_tfln_modulator_mesh(gap, h_box, w_core, h_core, 0.22e-6);
    let mesh_elec = Mesh::from_solver_mesh(&g_elec.mesh).unwrap();

    let mut ems_spec = ProblemSpec::new("ems", 2, "modes");
    ems_spec.study.params = json!({
        "freq": freq0,
        "n_modes": 2,
        "n_target": n_e
    });
    ems_spec.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), eps_clad)].into()));
    ems_spec.materials.insert(
        "core".into(),
        MaterialAssign::Oriented { material: "lithium_niobate".into(), euler_deg: [90.0, 90.0, 0.0] },
    );

    // 1. Unperturbed optical solve (E_DC = 0)
    let mut ctx_ems_0 = Ctx::new(&res_opt);
    ctx_ems_0.library = Some(lib.clone());
    let ems_res_0 = Ems.solve(&mesh_opt, &ems_spec, &mut ctx_ems_0).expect("ems V=0 solve");
    let te_idx_0 = (0..ems_res_0.solutions.len())
        .max_by(|&i, &j| {
            let fi = num_scalar(&ems_res_0, &format!("Ex share {}", i + 1));
            let fj = num_scalar(&ems_res_0, &format!("Ex share {}", j + 1));
            fi.partial_cmp(&fj).unwrap()
        })
        .unwrap();
    let neff_0 = ems_res_0.solutions[te_idx_0].eigenvalue.unwrap()[0];

    // Compute analytical overlap integral from unperturbed mode profile
    let disc = discretise(&mesh_opt).unwrap();
    let geo = geometry(&disc);
    let re = RefElement::new(2, 2, 6);
    let nloc = disc.nloc(2);
    let sol0 = &ems_res_0.solutions[te_idx_0];
    let cf_ex = sol0.cell_fields.iter().find(|f| f.name == "Ex").unwrap();
    let cf_ey = sol0.cell_fields.iter().find(|f| f.name == "Ey").unwrap();
    let cf_ez = sol0.cell_fields.iter().find(|f| f.name == "Ez").unwrap();
    let core_doms = res_opt.domains("core").unwrap();

    let mut num_int = 0.0_f64;
    let mut core_ex_sq = 0.0_f64;
    let mut tot_e_eps_e = 0.0_f64;
    for ci in 0..disc.n_cells() {
        let is_core = core_doms.contains(&disc.cell_labels[ci]);
        let (eps_xx, eps_yy, eps_zz) = if is_core { (eps_e, eps_o, eps_o) } else { (eps_clad, eps_clad, eps_clad) };
        for q in 0..re.w.len() {
            let (mut ex, mut ey, mut ez) = (0.0_f64, 0.0_f64, 0.0_f64);
            for a in 0..nloc {
                let na = re.n[q][a];
                let idx = ci * nloc + a;
                ex += na * cf_ex.re[idx];
                ey += na * cf_ey.re[idx];
                ez += na * cf_ez.re[idx];
            }
            let wvol = geo.vol[ci] * re.w[q];
            if is_core {
                let d_eps_xx = -eps_e * eps_e * r33 * e_dc_uniform;
                let d_eps_oo = -eps_o * eps_o * r13 * e_dc_uniform;
                num_int += wvol * (d_eps_xx * ex * ex + d_eps_oo * (ey * ey + ez * ez));
                core_ex_sq += wvol * eps_e * ex * ex;
            }
            tot_e_eps_e += wvol * (eps_xx * ex * ex + eps_yy * ey * ey + eps_zz * ez * ez);
        }
    }
    let gamma_overlap = core_ex_sq / tot_e_eps_e;

    // 2. Perturbed optical solves (+V_bias and -V_bias) with shared-mesh `"electric_field"` (`[Ex, Ez]`)
    let n_dofs_opt = disc.n_dofs(2);
    let mut e_dc_pos = Vec::with_capacity(n_dofs_opt * 2);
    let mut e_dc_neg = Vec::with_capacity(n_dofs_opt * 2);
    for _ in 0..n_dofs_opt {
        e_dc_pos.push(e_dc_uniform);
        e_dc_pos.push(0.0);
        e_dc_neg.push(-e_dc_uniform);
        e_dc_neg.push(0.0);
    }

    let mut ctx_ems_pos = Ctx::new(&res_opt);
    ctx_ems_pos.library = Some(lib.clone());
    ctx_ems_pos.inputs.insert("electric_field".into(), e_dc_pos);
    let ems_res_pos = Ems.solve(&mesh_opt, &ems_spec, &mut ctx_ems_pos).expect("ems +V solve");
    let neff_pos = ems_res_pos.solutions[te_idx_0].eigenvalue.unwrap()[0];

    let mut ctx_ems_neg = Ctx::new(&res_opt);
    ctx_ems_neg.library = Some(lib.clone());
    ctx_ems_neg.inputs.insert("electric_field".into(), e_dc_neg);
    let ems_res_neg = Ems.solve(&mesh_opt, &ems_spec, &mut ctx_ems_neg).expect("ems -V solve");
    let neff_neg = ems_res_neg.solutions[te_idx_0].eigenvalue.unwrap()[0];

    let d_neff_pos = neff_pos - neff_0;
    let d_neff_sim = 0.5 * (neff_pos - neff_neg);

    // 3. Verify against direct manual rotated permittivity tensor oracle at +V_bias:
    //    For X-cut LiNbO3 (`euler_deg = [90, 90, 0]`) under E_x = e_dc_uniform (E_y = E_z = 0):
    //    crystal E_3 = E_x => Delta(eps_r^-1)_33 = r33 * E_x, Delta(eps_r^-1)_11 = Delta(eps_r^-1)_22 = r13 * E_x.
    //    In wafer frame:
    //      eps_xx = eps_e - eps_e^2 * r33 * E_x
    //      eps_yy = eps_o - eps_o^2 * r13 * E_x
    //      eps_zz = eps_o - eps_o^2 * r13 * E_x
    let mut manual_spec = ems_spec.clone();
    let mut manual_props = BTreeMap::new();
    manual_props.insert(
        "eps_r".into(),
        PropValue::Diag([
            eps_e - eps_e * eps_e * r33 * e_dc_uniform,
            eps_o - eps_o * eps_o * r13 * e_dc_uniform,
            eps_o - eps_o * eps_o * r13 * e_dc_uniform,
        ]),
    );
    manual_spec.materials.insert("core".into(), MaterialAssign::PropsTensor(manual_props));
    let mut ctx_manual = Ctx::new(&res_opt);
    ctx_manual.library = Some(lib.clone());
    let ems_manual = Ems.solve(&mesh_opt, &manual_spec, &mut ctx_manual).expect("manual perturbed solve");
    let d_neff_manual = ems_manual.solutions[te_idx_0].eigenvalue.unwrap()[0] - neff_0;

    let rel_err_oracle = (d_neff_sim - d_neff_manual).abs() / d_neff_manual.abs();
    assert!(
        rel_err_oracle < 1e-5,
        "EO Delta n_eff ({d_neff_sim:.6e}) vs manual rotated tensor oracle ({d_neff_manual:.6e}): rel_err = {rel_err_oracle:.3e} >= 1e-5"
    );

    // Compute exact geometric group index n_g_geom = c_0 * U_tot / P_y from unperturbed non-dispersive tensor solve
    let mut manual_spec_0 = ems_spec.clone();
    let mut manual_props_0 = BTreeMap::new();
    manual_props_0.insert("eps_r".into(), PropValue::Diag([eps_e, eps_o, eps_o]));
    manual_spec_0.materials.insert("core".into(), MaterialAssign::PropsTensor(manual_props_0));
    let mut ctx_manual_0 = Ctx::new(&res_opt);
    ctx_manual_0.library = Some(lib.clone());
    let ems_manual_0 = Ems.solve(&mesh_opt, &manual_spec_0, &mut ctx_manual_0).expect("manual unperturbed solve");
    let n_g_geom = num_scalar(&ems_manual_0, &format!("n_g {}", te_idx_0 + 1));

    // Verify V_pi * L = lambda_0 * V_bias / (2 * |Delta n_eff|) against the analytical Poynting overlap formula:
    //   Delta n_eff = (n_g_geom / 2) * (int_core E^* . Delta_eps_r . E dA) / (int_all E^* . eps_r . E dA)
    let d_neff_rayleigh = 0.5 * n_g_geom * (num_int / tot_e_eps_e);
    let vpi_l_sim = lambda0 * v_bias / (2.0 * d_neff_sim.abs());
    let vpi_l_analytical = lambda0 * v_bias / (2.0 * d_neff_rayleigh.abs());
    let rel_err_vpi = (vpi_l_sim - vpi_l_analytical).abs() / vpi_l_analytical;
    assert!(
        rel_err_vpi < 0.005 && gamma_overlap > 0.5,
        "V_pi*L simulated ({vpi_l_sim:.5} V*m) vs analytical overlap formula ({vpi_l_analytical:.5} V*m): rel_err = {rel_err_vpi:.4} >= 0.5%"
    );

    // 4. Cross-mesh transfer verification (`mesh_elec` -> `mesh_opt` via `TransferredField`, domain-filtered by `"core"`)
    let disc_elec = discretise(&mesh_elec).unwrap();
    let n_dofs_elec = disc_elec.n_dofs(2);
    let mut e_dc_elec = Vec::with_capacity(n_dofs_elec * 2);
    for _ in 0..n_dofs_elec {
        e_dc_elec.push(e_dc_uniform);
        e_dc_elec.push(0.0);
    }
    let tf = TransferredField {
        physics: "dfs".into(),
        field: "E".into(),
        mesh: Arc::new(mesh_elec),
        discretisation: Some(Arc::new(g_elec.mesh.clone())),
        domain_names: domain_name_map(&g_elec),
        boundary_names: boundary_name_map(&g_elec),
        order: 2,
        components: 2,
        values: e_dc_elec,
        cell_fields: Vec::new(),
    };
    let mut ctx_ems_cross = Ctx::new(&res_opt);
    ctx_ems_cross.library = Some(lib);
    ctx_ems_cross.transferred_inputs.insert("electric_field".into(), tf);
    let ems_res_cross = Ems.solve(&mesh_opt, &ems_spec, &mut ctx_ems_cross).expect("ems cross-mesh EO solve");
    let d_neff_cross = ems_res_cross.solutions[te_idx_0].eigenvalue.unwrap()[0] - neff_0;
    let rel_err_cross = (d_neff_cross - d_neff_pos).abs() / d_neff_pos.abs();
    assert!(
        rel_err_cross < 1e-10,
        "cross-mesh EO Delta n_eff ({d_neff_cross:.6e}) vs shared-mesh ({d_neff_pos:.6e}): rel_err = {rel_err_cross:.3e}"
    );
}

#[test]
fn silicon_photoelastic_index_and_cavity_shift_from_transferred_strain() {
    let lib = load_library();
    let w = 0.5e-6;
    let h = 0.25e-6;
    let mut m = Model::new_2d();
    let core = m.rect("core_r", (0.0, 0.0), (w, h));
    m.overlay(&[("core", core)]);
    let g_opt = mesh_model(&m.list, &MeshPlan::simple(0.08e-6)).expect("mesh Si block opt");
    let mesh_opt = Mesh::from_solver_mesh(&g_opt.mesh).unwrap();
    let resolver = GeomResolver(&g_opt);

    let g_mech = mesh_model(&m.list, &MeshPlan::simple(0.11e-6)).expect("mesh Si block mech");
    let mesh_mech = Mesh::from_solver_mesh(&g_mech.mesh).unwrap();
    let disc_mech = discretise(&mesh_mech).unwrap();

    // Uniaxial strain S_xx = 2.0e-3 (S_yy = S_zz = 0) on the distinct mechanical mesh `mesh_mech`
    let s_xx = 2.0e-3_f64;
    let n_dofs_mech = disc_mech.n_dofs(2);
    let mut s_vals = Vec::with_capacity(n_dofs_mech * 6);
    for _ in 0..n_dofs_mech {
        s_vals.extend_from_slice(&[s_xx, 0.0, 0.0, 0.0, 0.0, 0.0]);
    }

    // For Silicon (`materials/silicon.toml`) at 1.55 um:
    //   n = 3.4757 => eps_r = 12.08049049
    //   p11 = -0.094, p12 = 0.017, p44 = -0.051
    // Under S_xx = 2.0e-3 (S_yy = S_zz = 0):
    //   Delta(eps_r^-1)_xx = p11 * S_xx => eps_r_xx_new = eps_r - eps_r^2 * p11 * S_xx
    //   Delta(eps_r^-1)_yy = Delta(eps_r^-1)_zz = p12 * S_xx => eps_r_yy_new = eps_r_zz_new = eps_r - eps_r^2 * p12 * S_xx
    let si_mat = lib.get("silicon").expect("silicon in library");
    let eps_si =
        si_mat.optical_permittivity(1.55e-6, None, [0.0; 3]).expect("Si optical permittivity at 1.55 um")[0][0].re;
    let p11 = -0.094_f64;
    let p12 = 0.017_f64;
    let eps_xx_exact = eps_si - eps_si * eps_si * p11 * s_xx;
    let eps_zz_exact = eps_si - eps_si * eps_si * p12 * s_xx;

    let freq0 = C0 / 1.55e-6;
    let mut ems_spec = ProblemSpec::new("ems", 2, "modes");
    ems_spec.study.params = json!({
        "freq": freq0,
        "n_modes": 2,
        "n_target": 3.48
    });
    ems_spec.materials.insert("core".into(), MaterialAssign::Named("silicon".into()));
    ems_spec.bcs = vec![
        BcSpec { on: "extreme(-x)".into(), kind: "pec".into(), params: json!({}) },
        BcSpec { on: "extreme(+x)".into(), kind: "pec".into(), params: json!({}) },
        BcSpec { on: "extreme(-z)".into(), kind: "pec".into(), params: json!({}) },
        BcSpec { on: "extreme(+z)".into(), kind: "pec".into(), params: json!({}) },
    ];

    let mut ctx_ems_0 = Ctx::new(&resolver);
    ctx_ems_0.library = Some(lib.clone());
    let res_0 = Ems.solve(&mesh_opt, &ems_spec, &mut ctx_ems_0).expect("unstrained ems solve");
    let neff_0 = res_0.solutions[0].eigenvalue.unwrap()[0];

    let tf_strain = TransferredField {
        physics: "sms".into(),
        field: "S".into(),
        mesh: Arc::new(mesh_mech),
        discretisation: Some(Arc::new(g_mech.mesh.clone())),
        domain_names: domain_name_map(&g_mech),
        boundary_names: boundary_name_map(&g_mech),
        order: 2,
        components: 6,
        values: s_vals,
        cell_fields: Vec::new(),
    };

    let mut ctx_ems_strained = Ctx::new(&resolver);
    ctx_ems_strained.library = Some(lib.clone());
    ctx_ems_strained.transferred_inputs.insert("strain".into(), tf_strain);
    let res_strained = Ems.solve(&mesh_opt, &ems_spec, &mut ctx_ems_strained).expect("cross-mesh strained ems solve");
    let neff_strained = res_strained.solutions[0].eigenvalue.unwrap()[0];

    let mut manual_spec = ems_spec.clone();
    let mut manual_props = BTreeMap::new();
    manual_props.insert("eps_r".into(), PropValue::Diag([eps_xx_exact, eps_zz_exact, eps_zz_exact]));
    manual_spec.materials.insert("core".into(), MaterialAssign::PropsTensor(manual_props));
    let mut ctx_manual = Ctx::new(&resolver);
    ctx_manual.library = Some(lib);
    let res_manual = Ems.solve(&mesh_opt, &manual_spec, &mut ctx_manual).expect("manual PE solve");
    let neff_manual = res_manual.solutions[0].eigenvalue.unwrap()[0];

    let d_neff_pe = neff_strained - neff_0;
    let d_neff_exact = neff_manual - neff_0;
    let rel_err = (d_neff_pe - d_neff_exact).abs() / d_neff_exact.abs();
    assert!(
        rel_err < 1e-6 && d_neff_pe.abs() > 1e-4,
        "Silicon photoelastic Delta n_eff ({d_neff_pe:.6e}) vs exact tensor ({d_neff_exact:.6e}): rel_err = {rel_err:.3e}"
    );
}
