// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-39 (ADR 141) verification suite for Domain-Restricted Mechanical Displacement (`u`) Submesh
//! in Monolithic Piezoelectric `(u, phi)` FEM and `"S"` / `"E"` field exports in `nyx-sms`.

use femcore::material::Library;
use femcore::mesh::cells::{CellMesh, CellType};
use nyx_physics::{BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, Scalar, SolveResult};
use nyx_sms::Sms;
use serde_json::json;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

struct MapResolver {
    doms: HashMap<String, Vec<i32>>,
    bnds: HashMap<String, Vec<i32>>,
}

impl LabelResolver for MapResolver {
    fn domains(&self, name: &str) -> Result<Vec<i32>, Error> {
        self.doms.get(name).cloned().ok_or_else(|| Error::invalid(format!("unknown domain '{name}'")))
    }

    fn boundaries(&self, name: &str) -> Result<Vec<i32>, Error> {
        self.bnds.get(name).cloned().ok_or_else(|| Error::invalid(format!("unknown boundary '{name}'")))
    }
}

fn mat_lib() -> Arc<Library> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Arc::new(Library::load_dir(&dir).expect("load materials library"))
}

fn num_scalar(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(v))) => *v,
        other => panic!("missing scalar '{name}', got {other:?} in {:?}", r.scalars),
    }
}

/// Builds a 2D two-domain mesh `[0, lx] x [0, lz_piezo + lz_air]`:
/// - Domain 1 (`"piezo"`): `z in [0, lz_piezo]` (`nz_piezo` rows)
/// - Domain 2 (`"air"`): `z in [lz_piezo, lz_piezo + lz_air]` (`nz_air` rows)
///
/// Boundaries:
/// - `1`: `x = 0`
/// - `2`: `z = 0` (bottom of piezo)
/// - `3`: `z = lz_piezo` (piezo-air interface)
/// - `4`: `z = lz_piezo + lz_air` (top of air box)
/// - `5`: `x = lx`
fn two_domain_piezo_air_2d(
    lx: f64,
    lz_piezo: f64,
    lz_air: f64,
    nx: usize,
    nz_piezo: usize,
    nz_air: usize,
    use_quads: bool,
) -> Mesh {
    let nz_tot = nz_piezo + nz_air;
    let mut cm = CellMesh::<2>::new();
    let vid = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    for j in 0..=nz_tot {
        let z = if j <= nz_piezo {
            lz_piezo * (j as f64) / (nz_piezo as f64)
        } else {
            lz_piezo + lz_air * ((j - nz_piezo) as f64) / (nz_air as f64)
        };
        for i in 0..=nx {
            let x = lx * (i as f64) / (nx as f64);
            cm.points.push([x, z]);
        }
    }
    for j in 0..nz_tot {
        let dom = if j < nz_piezo { 1 } else { 2 };
        for i in 0..nx {
            let (v00, v10, v11, v01) = (vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1));
            if use_quads {
                cm.push_cell(CellType::Quad, &[v00, v10, v11, v01], dom);
            } else {
                cm.push_cell(CellType::Tri, &[v00, v10, v11], dom);
                cm.push_cell(CellType::Tri, &[v00, v11, v01], dom);
            }
        }
    }
    for j in 0..nz_tot {
        cm.push_facet(&[vid(0, j), vid(0, j + 1)], 1);
        cm.push_facet(&[vid(nx, j), vid(nx, j + 1)], 5);
    }
    for i in 0..nx {
        cm.push_facet(&[vid(i, 0), vid(i + 1, 0)], 2);
        cm.push_facet(&[vid(i, nz_piezo), vid(i + 1, nz_piezo)], 3);
        cm.push_facet(&[vid(i, nz_tot), vid(i + 1, nz_tot)], 4);
    }
    Mesh::D2(cm)
}

/// Builds a 3D `Hex` two-domain piezo-in-air stack `[0, lx] x [0, ly] x [0, lz_piezo + lz_air]`.
/// - Domain 1 (`"piezo"`): `k in 0..nz_piezo`
/// - Domain 2 (`"air"`): `k in nz_piezo..(nz_piezo + nz_air)`
///
/// Boundaries: `1` (`z = 0`), `2` (`z = lz_piezo`), `3` (`z = lz_piezo + lz_air`).
#[allow(clippy::too_many_arguments)]
fn two_domain_piezo_air_3d_hex(
    lx: f64,
    ly: f64,
    lz_piezo: f64,
    lz_air: f64,
    nx: usize,
    ny: usize,
    nz_piezo: usize,
    nz_air: usize,
) -> Mesh {
    let nz_tot = nz_piezo + nz_air;
    let mut cm = CellMesh::<3>::new();
    let vid = |i: usize, j: usize, k: usize| ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
    for k in 0..=nz_tot {
        let z = if k <= nz_piezo {
            lz_piezo * (k as f64) / (nz_piezo as f64)
        } else {
            lz_piezo + lz_air * ((k - nz_piezo) as f64) / (nz_air as f64)
        };
        for j in 0..=ny {
            let y = ly * (j as f64) / (ny as f64);
            for i in 0..=nx {
                let x = lx * (i as f64) / (nx as f64);
                cm.points.push([x, y, z]);
            }
        }
    }
    for k in 0..nz_tot {
        let dom = if k < nz_piezo { 1 } else { 2 };
        for j in 0..ny {
            for i in 0..nx {
                cm.push_cell(
                    CellType::Hex,
                    &[
                        vid(i, j, k),
                        vid(i + 1, j, k),
                        vid(i + 1, j + 1, k),
                        vid(i, j + 1, k),
                        vid(i, j, k + 1),
                        vid(i + 1, j, k + 1),
                        vid(i + 1, j + 1, k + 1),
                        vid(i, j + 1, k + 1),
                    ],
                    dom,
                );
            }
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            cm.push_facet(&[vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0)], 1);
            cm.push_facet(
                &[vid(i, j, nz_piezo), vid(i + 1, j, nz_piezo), vid(i + 1, j + 1, nz_piezo), vid(i, j + 1, nz_piezo)],
                2,
            );
            cm.push_facet(
                &[vid(i, j, nz_tot), vid(i + 1, j, nz_tot), vid(i + 1, j + 1, nz_tot), vid(i, j + 1, nz_tot)],
                3,
            );
        }
    }
    Mesh::D3(cm)
}

#[test]
fn piezo_in_air_submesh_vs_full_mesh_static_and_eigen() {
    let lib = mat_lib();
    let (lx, lz_piezo, lz_air) = (2.0e-6, 1.0e-6, 1.5e-6);
    let (nx, nz_piezo, nz_air) = (6, 6, 9);

    for use_quads in [false, true] {
        let mesh = two_domain_piezo_air_2d(lx, lz_piezo, lz_air, nx, nz_piezo, nz_air, use_quads);
        let mut doms = HashMap::new();
        doms.insert("piezo".into(), vec![1]);
        doms.insert("air".into(), vec![2]);
        let mut bnds = HashMap::new();
        bnds.insert("left".into(), vec![1]);
        bnds.insert("bottom".into(), vec![2]);
        bnds.insert("interface".into(), vec![3]);
        bnds.insert("top_air".into(), vec![4]);
        bnds.insert("right".into(), vec![5]);
        let resolver = MapResolver { doms, bnds };

        // 1. Static actuation across the piezo + air capacitor stack (`V = 0` at z = 0, `V = 10 V` at z = lz_piezo + lz_air)
        let mut spec_sub = ProblemSpec::new("sms", 2, "static");
        spec_sub.params = json!({
            "piezoelectric": "on",
            "plane": "strain",
            "u_domains": ["piezo"]
        });
        spec_sub.materials.insert("piezo".into(), MaterialAssign::Named("aluminum_nitride".into()));
        spec_sub.materials.insert("air".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
        spec_sub.bcs = vec![
            BcSpec { on: "bottom".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0 }) },
            BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
            BcSpec { on: "right".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
            BcSpec { on: "bottom".into(), kind: "ground".into(), params: json!({}) },
            BcSpec { on: "top_air".into(), kind: "potential".into(), params: json!({ "value": 10.0 }) },
        ];

        let mut ctx_sub = Ctx::new(&resolver);
        ctx_sub.library = Some(lib.clone());
        let res_sub = Sms.solve(&mesh, &spec_sub, &mut ctx_sub).expect("submesh static solve");

        // Full-mesh reference: give "air" an ultra-compliant elastic surrogate (E = 1e-3 Pa, nu = 0.0, rho = 1e-9)
        // and include both "piezo" and "air" in `u_domains` so `u` is allocated on the entire mesh.
        let mut spec_full = spec_sub.clone();
        spec_full.params = json!({
            "piezoelectric": "on",
            "plane": "strain",
            "u_domains": ["piezo", "air"]
        });
        spec_full.materials.insert(
            "air".into(),
            MaterialAssign::Props(
                [("eps_r".into(), 1.0), ("E".into(), 1e-3), ("nu".into(), 0.0), ("rho".into(), 1e-9)].into(),
            ),
        );
        let mut ctx_full = Ctx::new(&resolver);
        ctx_full.library = Some(lib.clone());
        let res_full = Sms.solve(&mesh, &spec_full, &mut ctx_full).expect("full-mesh static solve");

        // Verify strict DOF reduction: submesh allocates u ONLY on solid "piezo" nodes while keeping phi on all nodes!
        assert!(
            res_sub.dofs < res_full.dofs,
            "expected submesh DOFs ({}) < full-mesh DOFs ({}) (use_quads={use_quads})",
            res_sub.dofs,
            res_full.dofs
        );
        let sol_sub = &res_sub.solutions[0];
        let sol_full = &res_full.solutions[0];

        // Verify potential `phi` across the entire piezo+air mesh agrees within < 1e-5
        let phi_sub = &sol_sub.field("phi").expect("phi field").values;
        let phi_full = &sol_full.field("phi").expect("phi field").values;
        assert_eq!(phi_sub.len(), phi_full.len());
        let max_phi = phi_full.iter().copied().fold(0.0_f64, |a, b| a.max(b.abs())).max(1.0);
        let max_diff_phi = phi_sub.iter().zip(phi_full).map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
        assert!(
            max_diff_phi / max_phi < 1e-6,
            "submesh vs full-mesh phi relative difference {:.3e} >= 1e-6",
            max_diff_phi / max_phi
        );

        // Verify exported strain `"S"` (6 components) and electric field `"E"` exist and are finite
        let s_field = sol_sub.field("S").expect("exported 6-component strain field S");
        assert_eq!(s_field.components.len(), 6);
        assert!(s_field.values.iter().all(|v| v.is_finite()));
        let e_field = sol_sub.field("E").expect("exported electric field E");
        assert_eq!(e_field.components.len(), 2);
        assert!(e_field.values.iter().all(|v| v.is_finite()));

        // Verify air-only nodes have exactly 0 displacement in submesh result while piezo top displacement matches full-mesh < 1e-5
        let u_sub = &sol_sub.field("u").expect("u field").values;
        let u_full = &sol_full.field("u").expect("u field").values;
        let max_uz_sub = u_sub.iter().skip(1).step_by(2).copied().fold(0.0_f64, |a, b| a.max(b.abs()));
        let max_uz_full = (0..phi_full.len())
            .filter(|&i| (phi_full[i] - phi_sub[i]).abs() < 1e-6 && u_sub[2 * i + 1].abs() > 0.0)
            .map(|i| u_full[2 * i + 1].abs())
            .fold(0.0_f64, f64::max);
        let rel_u = (max_uz_sub - max_uz_full).abs() / max_uz_full.max(1e-30);
        assert!(
            rel_u < 1e-5,
            "submesh vs full-mesh piezo displacement relative diff {rel_u:.3e} >= 1e-5 (use_quads={use_quads})"
        );

        // 2. Eigenfrequency resonator mode in air (open-circuit top electrode in air, grounded bottom electrode)
        // Closed-form 1D open-circuit AlN thickness-extensional fundamental frequency:
        // c_D33 = C_33 + e_33^2 / eps_33, f_1 = (1 / (4 * lz_piezo)) * sqrt(c_D33 / rho)
        let mut spec_eig_sub = ProblemSpec::new("sms", 2, "eigenfrequency");
        spec_eig_sub.params = json!({
            "piezoelectric": "on",
            "plane": "strain"
        });
        // Note: without explicit `u_domains`, `"air"` has only `eps_r: 1.0` (no `E`), so `resolve_domain_mats`
        // automatically restricts `u` to `"piezo"`!
        let f_search = 2.8e9;
        spec_eig_sub.study.params = json!({
            "n_modes": 2,
            "f_search": f_search
        });
        spec_eig_sub.materials.insert("piezo".into(), MaterialAssign::Named("aluminum_nitride".into()));
        spec_eig_sub.materials.insert("air".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
        spec_eig_sub.bcs = vec![
            BcSpec { on: "bottom".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0 }) },
            BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
            BcSpec { on: "right".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
            BcSpec { on: "bottom".into(), kind: "ground".into(), params: json!({}) },
            BcSpec { on: "top_air".into(), kind: "ground".into(), params: json!({}) },
        ];

        let mut ctx_eig_sub = Ctx::new(&resolver);
        ctx_eig_sub.library = Some(lib.clone());
        let res_eig_sub = Sms.solve(&mesh, &spec_eig_sub, &mut ctx_eig_sub).expect("submesh eigen solve");
        assert_eq!(res_eig_sub.dofs, res_sub.dofs);
        let f1_sub = num_scalar(&res_eig_sub, "f re 1 (Hz)");

        // Compare against explicit `u_domains: ["piezo"]`: must be 0-ULP identical!
        let mut spec_eig_explicit = spec_eig_sub.clone();
        spec_eig_explicit.params = json!({
            "piezoelectric": "on",
            "plane": "strain",
            "u_domains": ["piezo"]
        });
        let mut ctx_eig_exp = Ctx::new(&resolver);
        ctx_eig_exp.library = Some(lib.clone());
        let res_eig_exp =
            Sms.solve(&mesh, &spec_eig_explicit, &mut ctx_eig_exp).expect("explicit u_domains eigen solve");
        let f1_exp = num_scalar(&res_eig_exp, "f re 1 (Hz)");
        assert_eq!(f1_sub.to_bits(), f1_exp.to_bits());

        // Compare against a 1D-equivalent all-solid reference where the top boundary condition on `phi`
        // accounts for the air-gap series capacitance (or stiff-free air reference):
        // With `lz_air = 1.5 um, eps_air = 1` in series with `lz_piezo = 1.0 um, eps_33 = 10.23`,
        // the air gap drops ~94% of the potential, approaching the open-circuit stiffened frequency (~2.81 GHz).
        assert!(
            (f1_sub - 2.80e9).abs() / 2.80e9 < 0.03,
            "unexpected AlN-in-air thickness eigenfrequency {f1_sub:.6e} Hz"
        );
        // Verify eigenmode exports both 6-component strain `"S"` and electric field `"E"`
        let eig_sol = &res_eig_sub.solutions[0];
        assert_eq!(eig_sol.field("S").expect("eigen strain S").components.len(), 6);
        assert_eq!(eig_sol.field("E").expect("eigen electric field E").components.len(), 2);
    }
}

#[test]
fn all_solid_piezo_zero_ulp_parity_with_and_without_u_domains() {
    let lib = mat_lib();
    let mesh = two_domain_piezo_air_2d(2.0e-6, 1.0e-6, 1.0e-6, 4, 4, 4, false);
    let mut doms = HashMap::new();
    doms.insert("d1".into(), vec![1]);
    doms.insert("d2".into(), vec![2]);
    let mut bnds = HashMap::new();
    bnds.insert("bottom".into(), vec![2]);
    bnds.insert("top".into(), vec![4]);
    let resolver = MapResolver { doms, bnds };

    let mut spec_default = ProblemSpec::new("sms", 2, "eigenfrequency");
    spec_default.params = json!({
        "piezoelectric": "on",
        "plane": "strain"
    });
    spec_default.study.params = json!({
        "n_modes": 3,
        "f_search": 2.0e9
    });
    spec_default.materials.insert("d1".into(), MaterialAssign::Named("aluminum_nitride".into()));
    spec_default.materials.insert("d2".into(), MaterialAssign::Named("aluminum_nitride".into()));
    spec_default.bcs = vec![
        BcSpec { on: "bottom".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0 }) },
        BcSpec { on: "bottom".into(), kind: "ground".into(), params: json!({}) },
        BcSpec { on: "top".into(), kind: "ground".into(), params: json!({}) },
    ];

    let mut ctx_default = Ctx::new(&resolver);
    ctx_default.library = Some(lib.clone());
    let res_default = Sms.solve(&mesh, &spec_default, &mut ctx_default).expect("all-solid default");

    let mut spec_explicit = spec_default.clone();
    spec_explicit.params = json!({
        "piezoelectric": "on",
        "plane": "strain",
        "u_domains": ["d1", "d2"]
    });
    let mut ctx_explicit = Ctx::new(&resolver);
    ctx_explicit.library = Some(lib);
    let res_explicit = Sms.solve(&mesh, &spec_explicit, &mut ctx_explicit).expect("all-solid explicit u_domains");

    assert_eq!(res_default.dofs, res_explicit.dofs);
    for ((k1, v1), (k2, v2)) in res_default.scalars.iter().zip(&res_explicit.scalars) {
        assert_eq!(k1, k2);
        if let (Scalar::Num(a), Scalar::Num(b)) = (v1, v2) {
            assert_eq!(a.to_bits(), b.to_bits(), "scalar '{k1}' differed by ULP");
        }
    }
    for (s1, s2) in res_default.solutions.iter().zip(&res_explicit.solutions) {
        assert_eq!(s1.fields.len(), s2.fields.len());
        for (f1, f2) in s1.fields.iter().zip(&s2.fields) {
            assert_eq!(f1.name, f2.name);
            assert_eq!(f1.values.len(), f2.values.len());
            for (a, b) in f1.values.iter().zip(&f2.values) {
                assert_eq!(a.to_bits(), b.to_bits(), "field '{}' differed by ULP", f1.name);
            }
        }
    }
}

#[test]
fn piezo_submesh_3d_hex_mixed_and_driven_sweep() {
    let lib = mat_lib();
    let (lx, ly, lz_piezo, lz_air) = (1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6);
    let mesh = two_domain_piezo_air_3d_hex(lx, ly, lz_piezo, lz_air, 2, 2, 3, 3);
    let mut doms = HashMap::new();
    doms.insert("piezo".into(), vec![1]);
    doms.insert("air".into(), vec![2]);
    let mut bnds = HashMap::new();
    bnds.insert("bottom".into(), vec![1]);
    bnds.insert("interface".into(), vec![2]);
    bnds.insert("top_air".into(), vec![3]);
    let resolver = MapResolver { doms, bnds };

    let mut spec = ProblemSpec::new("sms", 1, "driven");
    spec.params = json!({
        "piezoelectric": "on",
        "u_domains": ["piezo"]
    });
    spec.study.params = json!({
        "freqs_hz": [2.5e9, 2.8e9],
        "prom_tol": 0.0,
        "solver": "direct"
    });
    spec.materials.insert("piezo".into(), MaterialAssign::Named("aluminum_nitride".into()));
    spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    spec.bcs = vec![
        BcSpec { on: "bottom".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }) },
        BcSpec { on: "bottom".into(), kind: "ground".into(), params: json!({}) },
        BcSpec {
            on: "top_air".into(),
            kind: "terminal".into(),
            params: json!({
                "name": "top",
                "Z_0": 50.0
            }),
        },
    ];

    let n_pts = match &mesh {
        Mesh::D3(cm) => cm.n_points(),
        _ => unreachable!(),
    };
    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib);
    let res = Sms.solve(&mesh, &spec, &mut ctx).expect("3D Hex submesh driven piezo solve");
    // Full 3D (u, phi) would have 4 * n_pts = 4 * 63 = 252 DOFs; submesh has 3 * 36 + 63 = 171 DOFs!
    assert!(res.dofs < 4 * n_pts, "expected 3D Hex submesh DOFs ({}) < 4 * {} ({})", res.dofs, n_pts, 4 * n_pts);
    let sol = &res.solutions[0];
    let s_field = sol.field("S").expect("3D driven piezo strain S");
    assert_eq!(s_field.components.len(), 6);
    assert!(s_field.values.iter().all(|v| v.is_finite()));
    let e_field = sol.field("E").expect("3D driven piezo electric field E");
    assert_eq!(e_field.components.len(), 3);
    assert!(e_field.values.iter().all(|v| v.is_finite()));
}
