// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! AF-6 integration tests for AC quasi-static dielectric & conductive loss (`tan_delta`, `sigma`, `freq`),
//! anisotropic loss tensors, library `silicon.toml` microwave loss at `13 GHz`, and the positive 2nd-order
//! capacitance shift (`Delta C^(2) > 0`) on a strong-loss (`tan_delta = 0.5`) two-layer series capacitor.

use femcore::fem::{geometry, quad_rule, shape};
use femcore::material::Library;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use geom::model::{Model, RegionKind};
use nyx_dfs::ac_loss::solve_electrostatics_with_ac_loss;
use nyx_dfs::{Bc, BcKind, Dfs, Problem, Study, EPS0};
use nyx_physics::{
    decode::materials_by_label, CapacitanceMatrixData, Ctx, Error, LabelResolver, MaterialAssign, Physics, ProblemSpec,
    PropValue, StudySpec,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::path::Path;
use std::sync::Arc;

#[test]
fn weak_loss_13ghz_conductance_matches_exact_silicon_energy_participation_and_lossless_is_none() {
    let mut m = Model::new_2d();
    let sub = m.rect("sub", (-20.0e-6, -20.0e-6), (20.0e-6, 0.0));
    let vac = m.rect("vac", (-20.0e-6, 0.0), (20.0e-6, 20.0e-6));
    let p1 = m.rect("p1", (-10.0e-6, 0.0), (-2.0e-6, 1.0e-6));
    let p2 = m.rect("p2", (2.0e-6, 0.0), (10.0e-6, 1.0e-6));
    m.overlay_full(
        &[
            ("sub", sub, RegionKind::Material, None),
            ("vac", vac, RegionKind::Material, None),
            ("p1", p1, RegionKind::Hole, None),
            ("p2", p2, RegionKind::Hole, None),
        ],
        &[],
    );
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(3.0e-6)).unwrap();
    let sub_dom = g.domain("sub").unwrap();
    let vac_dom = g.domain("vac").unwrap();

    let mut p1_labs = Vec::new();
    let mut p2_labs = Vec::new();
    let mut outer_labs = Vec::new();
    let tol = 1e-9;
    for f in 0..g.mesh.n_facets() {
        let fv = g.mesh.facet(f);
        let a = g.mesh.point(fv[0] as usize);
        let b = g.mesh.point(fv[1] as usize);
        let xm = 0.5 * (a[0] + b[0]);
        let zm = 0.5 * (a[1] + b[1]);
        let lab = g.mesh.facet_labels[f];
        if zm >= -tol && zm <= 1.0e-6 + tol && (-10.0e-6 - tol..=-2.0e-6 + tol).contains(&xm) {
            if !p1_labs.contains(&lab) {
                p1_labs.push(lab);
            }
        } else if zm >= -tol && zm <= 1.0e-6 + tol && (2.0e-6 - tol..=10.0e-6 + tol).contains(&xm) {
            if !p2_labs.contains(&lab) {
                p2_labs.push(lab);
            }
        } else if ((xm.abs() - 20.0e-6).abs() < tol || (zm.abs() - 20.0e-6).abs() < tol) && !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }

    let eps_si = 11.45;
    let tan_delta_si = 1.0e-7;
    let freq = 13.0e9;
    let omega = 2.0 * PI * freq;

    // 1. Lossless run (tan_delta == 0, sigma == 0) -> conductance_matrix_s and q_dielectric_matrix must be None
    let mut prob_lossless = Problem::new(Study::Electrostatic, 2);
    prob_lossless.materials.set(sub_dom, "eps_r", eps_si);
    prob_lossless.materials.set(vac_dom, "eps_r", 1.0);
    prob_lossless.bcs = vec![
        Bc { labels: p1_labs.clone(), kind: BcKind::Terminal("t1".into()) },
        Bc { labels: p2_labs.clone(), kind: BcKind::Terminal("t2".into()) },
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
    ];
    let (_, ac_lossless) = solve_electrostatics_with_ac_loss(&g.mesh, &prob_lossless, Some(freq)).unwrap();
    assert!(ac_lossless.conductance_matrix_s.is_none());
    assert!(ac_lossless.q_dielectric_matrix.is_none());

    // 2. Lossy silicon substrate (tan_delta = 1e-7 at freq = 13 GHz)
    let mut prob = prob_lossless.clone();
    prob.materials.set(sub_dom, "tan_delta", tan_delta_si);
    let (res, ac) = solve_electrostatics_with_ac_loss(&g.mesh, &prob, Some(freq)).unwrap();
    let cap = CapacitanceMatrixData::of(&res).unwrap();
    let g_mat = ac.conductance_matrix_s.as_ref().expect("conductance matrix must be present");
    let q_mat = ac.q_dielectric_matrix.as_ref().expect("Q matrix must be present");

    // Symmetry to < 1e-14
    let g_scale = g_mat[0][0].abs().max(g_mat[1][1].abs()).max(1e-30);
    assert!(((g_mat[0][1] - g_mat[1][0]).abs() / g_scale) < 1e-14);

    // Direct cell-by-cell bulk electric energy participation ratio p_Si = U_Si / U_total for terminal 0
    let geo = geometry(&g.mesh);
    let quad2 = quad_rule(2, 4);
    let sol0 = &res.solutions[0];
    let ex = sol0.cell_fields.iter().find(|f| f.name == "Ex").unwrap();
    let ez = sol0.cell_fields.iter().find(|f| f.name == "Ez").unwrap();
    let nloc = ex.nloc;
    let (mut u_si, mut u_tot) = (0.0, 0.0);
    for ci in 0..g.mesh.n_cells() {
        let is_si = g.mesh.cell_labels[ci] == sub_dom;
        let eps = if is_si { eps_si } else { 1.0 };
        let mut e2_int = 0.0;
        for (q, l) in quad2.lam.iter().enumerate() {
            let (sh, _) = shape(&[l[0], l[1], l[2]], 2);
            let base = ci * nloc;
            let (mut vx, mut vz) = (0.0, 0.0);
            for (a, &sh_a) in sh.iter().enumerate().take(nloc) {
                vx += sh_a * ex.re[base + a];
                vz += sh_a * ez.re[base + a];
            }
            e2_int += quad2.w[q] * (vx * vx + vz * vz);
        }
        let ue = 0.5 * EPS0 * eps * geo.vol[ci] * e2_int;
        if is_si {
            u_si += ue;
        }
        u_tot += ue;
    }
    let p_si_exact = u_si / u_tot;
    let p_si_from_g = g_mat[0][0] / (omega * cap.maxwell_matrix_f[0][0] * tan_delta_si);
    assert!(
        (p_si_from_g - p_si_exact).abs() < 1e-13,
        "G_00 / (omega * C_00 * tan_delta_Si) = {p_si_from_g} must equal p_Si = {p_si_exact} to < 1e-13"
    );

    let q00 = q_mat[0][0].expect("Q_00 must be Some");
    assert!(q00.is_finite());
    let expected_q00 = 1.0 / (p_si_exact * tan_delta_si);
    assert!(((q00 - expected_q00) / expected_q00).abs() < 1e-13);
}

#[test]
fn dc_ohmic_sigma_conductance_with_zero_frequency_has_no_division_by_zero() {
    let w = 10.0e-6;
    let d = 5.0e-6;
    let sigma = 2.5e-3; // S/m
    let mut m = Model::new_2d();
    let blk = m.rect("blk", (0.0, 0.0), (w, d));
    m.overlay(&[("blk", blk)]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(2.5e-6)).unwrap();
    let dom = g.domain("blk").unwrap();
    let top = g.boundary("blk.top");
    let bot = g.boundary("blk.bottom");

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(dom, "eps_r", 4.0);
    prob.materials.set(dom, "sigma", sigma);
    prob.bcs = vec![Bc { labels: top, kind: BcKind::Terminal("top".into()) }, Bc { labels: bot, kind: BcKind::Ground }];

    let (_, ac) = solve_electrostatics_with_ac_loss(&g.mesh, &prob, Some(0.0)).unwrap();
    let g_mat = ac.conductance_matrix_s.expect("DC sigma > 0 must produce conductance_matrix_s");
    assert!(ac.q_dielectric_matrix.is_none(), "q_dielectric_matrix must be None at freq == 0");
    let g_exact = sigma * w / d;
    assert!(((g_mat[0][0] - g_exact) / g_exact).abs() < 1e-13);
}

#[test]
fn strong_loss_two_layer_series_capacitor_proves_positive_second_order_capacitance_shift() {
    // Two-layer series parallel-plate capacitor of width W = 10 um and total thickness d = 10 um:
    // Layer 1 (top, z in [0, 5 um]): eps_1 = 4.0, tan_delta_1 = 0.0
    // Layer 2 (bottom, z in [-5 um, 0]): eps_2 = 4.0, tan_delta_2 = 0.5
    let w = 10.0e-6;
    let d_half = 5.0e-6;
    let eps1 = 4.0;
    let eps2 = 4.0;
    let td2 = 0.5;
    let freq = 10.0e9;
    let omega = 2.0 * PI * freq;

    let mut m = Model::new_2d();
    let lay1 = m.rect("lay1", (0.0, 0.0), (w, d_half));
    let lay2 = m.rect("lay2", (0.0, -d_half), (w, 0.0));
    m.overlay(&[("lay1", lay1), ("lay2", lay2)]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(2.5e-6)).unwrap();
    let dom1 = g.domain("lay1").unwrap();
    let dom2 = g.domain("lay2").unwrap();
    let top = g.boundary("lay1.top");
    let bot = g.boundary("lay2.bottom");

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(dom1, "eps_r", eps1);
    prob.materials.set(dom2, "eps_r", eps2);
    prob.materials.set(dom2, "tan_delta", td2);
    prob.bcs = vec![Bc { labels: top, kind: BcKind::Terminal("top".into()) }, Bc { labels: bot, kind: BcKind::Ground }];

    let (res, ac) = solve_electrostatics_with_ac_loss(&g.mesh, &prob, Some(freq)).unwrap();
    let c0 = CapacitanceMatrixData::of(&res).unwrap().maxwell_matrix_f[0][0];
    let c0_exact = EPS0 * 2.0 * w / d_half; // C1 * C2 / (C1 + C2) with C1 = C2 = 4 * eps0 * w / d_half -> C0 = 2 * eps0 * w / d_half
    assert!(((c0 - c0_exact) / c0_exact).abs() < 1e-13);

    // 1. Pure 2nd-order shift Delta C^(2) = + (K_loss v)_f^T K_ff^-1 (K_loss v)_f = +1/4 * C^(0) * tan^2(delta_2) > 0
    let dc2 = ac.delta_c_second_order_f.as_ref().expect("Delta C^(2) must be computed when max tan_delta > 1e-2")[0][0];
    let dc2_exact = 0.25 * c0_exact * td2 * td2;
    assert!(dc2 > 0.0, "2nd-order capacitance shift must be strictly POSITIVE (+), got {dc2}");
    assert!(
        ((dc2 - dc2_exact) / dc2_exact).abs() < 1e-12,
        "Delta C^(2) = {dc2} must match +1/4 C0 tan^2(delta) = {dc2_exact}"
    );

    // 2. Full complex series capacitor admittance Y(omega) = G(omega) + i * omega * C(omega):
    //    tilde_C = C0 * ((4 + 2 * td2^2) - i * 2 * td2) / (4 + td2^2)
    let c_ac_exact = c0_exact * (4.0 + 2.0 * td2 * td2) / (4.0 + td2 * td2);
    let g_ac_exact = omega * c0_exact * (2.0 * td2) / (4.0 + td2 * td2);

    let c_ac_num = ac.corrected_maxwell_matrix_f[0][0];
    let g_ac_num = ac.conductance_matrix_s.as_ref().unwrap()[0][0];

    let rel_err_c = ((c_ac_num - c_ac_exact) / c_ac_exact).abs();
    let rel_err_g = ((g_ac_num - g_ac_exact) / g_ac_exact).abs();
    assert!(
        rel_err_c < 1.0e-3,
        "AC series capacitance C(omega) = {c_ac_num} must match exact {c_ac_exact} to < 0.1% (rel_err = {rel_err_c})"
    );
    assert!(
        rel_err_g < 1.0e-3,
        "AC series conductance G(omega) = {g_ac_num} must match exact {g_ac_exact} to < 0.1% (rel_err = {rel_err_g})"
    );
}

#[test]
fn silicon_library_at_13ghz_and_anisotropic_tensor_loss_decode_without_optical_error() {
    struct MapResolver;
    impl LabelResolver for MapResolver {
        fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
            match region {
                "si" => Ok(vec![1]),
                "aniso" => Ok(vec![2]),
                _ => Err(Error::invalid(region)),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            Err(Error::invalid(s))
        }
    }

    let mat_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    let lib = Arc::new(Library::load_dir(&mat_dir).expect("materials library must load"));

    let mut materials = BTreeMap::new();
    materials.insert("si".into(), MaterialAssign::Named("silicon".into()));
    let mut aniso_map = BTreeMap::new();
    aniso_map.insert("eps_r".into(), PropValue::Diag([9.4, 9.4, 11.6]));
    aniso_map.insert("tan_delta".into(), PropValue::Diag([1.0e-6, 1.0e-6, 5.0e-6]));
    aniso_map.insert("sigma".into(), PropValue::Diag([1.0e-5, 1.0e-5, 3.0e-5]));
    materials.insert("aniso".into(), MaterialAssign::PropsTensor(aniso_map));

    let spec = ProblemSpec {
        format: 1,
        physics: "dfs".into(),
        order: 2,
        study: StudySpec { kind: "electrostatic".into(), params: json!({ "freq": 13.0e9 }) },
        params: json!({}),
        materials,
        bcs: vec![],
    };

    let resolver = MapResolver;
    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib);
    ctx.frequency = Some(13.0e9);

    let decoded = materials_by_label(&spec, &ctx, &Dfs.properties())
        .expect("silicon.toml at 13 GHz in dfs electrostatic must not trigger [optical] wavelength error");
    assert!((decoded.get(1, "eps_r").unwrap_or(0.0) - 11.7).abs() < 1e-12);
    assert!((decoded.get(1, "tan_delta").unwrap_or(0.0) - 5.0e-6).abs() < 1e-15);

    let an2 = decoded.aniso.get(&2).expect("anisotropic loss tensors must be stored on domain 2");
    assert!(an2.permittivity.is_some());
    assert!(an2.permittivity_im.is_some());
    assert!(an2.conductivity.is_some());
}
