// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Finite-difference vs discrete adjoint sensitivity verification for `nyx-sms` (ADR 142):
//! static structural compliance (`-u^T dK/dp u + 2 df/dp^T u` via `opt::adjoint_gradient`) and
//! 3D solid / 2D elastic waveguide modal eigenfrequency sensitivities (`sens::eigenvalue_derivative`
//! and `sens::match_modes`).

use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::geomdrive3::mesh_geom_3d;
use femcore::topology::{compute_shape_velocity, displace_mesh};
use geom::model::Model;
use nyx_physics::{
    BcSpec, Ctx, DesignParam, Error, LabelResolver, MaterialAssign, Physics, ProblemSpec, SensitivityObjective,
    SensitivityRequest,
};
use nyx_sms::Sms;
use serde_json::json;
use std::collections::BTreeMap;

struct TestResolver<'a>(&'a MeshedGeom);

impl LabelResolver for TestResolver<'_> {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.domains(region);
        if v.is_empty() {
            Err(Error::invalid(format!("unknown domain '{region}'")))
        } else {
            Ok(v)
        }
    }

    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        let b = self.0.boundary(selector);
        if !b.is_empty() {
            return Ok(b);
        }
        self.0.select(selector).map_err(Error::invalid)
    }
}

fn inline_props(pairs: &[(&str, f64)]) -> MaterialAssign {
    let mut m = BTreeMap::new();
    for &(k, v) in pairs {
        m.insert(k.to_string(), v);
    }
    MaterialAssign::Props(m)
}

#[test]
fn solid_static_compliance_adjoint_matches_fd_below_1e_6() {
    // 2D cantilever beam [0, 1.2] x [0, t_beam], clamped at x = 0, vertical traction on x = 1.2
    let mut m = Model::new_2d();
    m.param("t_beam", 0.30);
    let beam = m.rect("beam", (0.0, 0.0), (1.2, "t_beam"));
    m.overlay(&[("beam", beam)]);

    let ev0 = evaluate(&m.list).unwrap();
    let meshed = mesh_geom_2d(ev0, &GeomMeshOpts::new(0.12)).unwrap();
    let mesh = meshed.into_exchange_mesh().unwrap();
    let vel_t = compute_shape_velocity(&m.list, &meshed.eval, &mesh, "t_beam", None).unwrap();
    let resolver = TestResolver(&meshed);

    let make_spec = |ym: f64, nu: f64| -> ProblemSpec {
        let mut spec = ProblemSpec::new("sms", 2, "static");
        spec.materials.insert("beam".into(), inline_props(&[("E", ym), ("nu", nu), ("rho", 2330.0)]));
        spec.bcs.push(BcSpec {
            on: "plane(x = 0)".into(),
            kind: "fixed".into(),
            params: json!({ "ux": 0.0, "uz": 0.0 }),
        });
        spec.bcs.push(BcSpec {
            on: "plane(x = 1.2)".into(),
            kind: "load".into(),
            params: json!({ "traction": [0.0, 0.0, -1.0e6] }),
        });
        spec
    };

    let ym0 = 1.6e11;
    let nu0 = 0.27;
    let spec0 = make_spec(ym0, nu0);
    let mut ctx = Ctx::new(&resolver);

    let mut req = SensitivityRequest::new(
        SensitivityObjective::Compliance { target: None },
        vec![
            DesignParam::ParamVar { name: "t_beam".into(), min: None, max: None },
            DesignParam::MaterialScalar { region: "beam".into(), property: "E".into(), min: None, max: None },
            DesignParam::MaterialScalar { region: "beam".into(), property: "nu".into(), min: None, max: None },
        ],
    );
    req.shape_velocities.insert("t_beam".into(), vel_t.clone());
    let rep = Sms.sensitivity(&mesh, &spec0, &mut ctx, &req).unwrap();

    let req_e_only = SensitivityRequest::new(
        SensitivityObjective::Compliance { target: None },
        vec![DesignParam::MaterialScalar { region: "beam".into(), property: "E".into(), min: None, max: None }],
    );

    // 1. Shape derivative wrt `t_beam` (both stiffness dK/dt and tip boundary traction area df/dt vary!)
    let ht = 1e-6;
    let mp = displace_mesh(&mesh, ht, &vel_t).unwrap();
    let mm = displace_mesh(&mesh, -ht, &vel_t).unwrap();
    let cp = Sms.sensitivity(&mp, &spec0, &mut ctx, &req_e_only).unwrap().objective;
    let cm = Sms.sensitivity(&mm, &spec0, &mut ctx, &req_e_only).unwrap().objective;
    let fd_t = (cp - cm) / (2.0 * ht);
    let adj_t = rep.gradients[0][0];
    let rel_t = (adj_t - fd_t).abs() / fd_t.abs();
    assert!(rel_t < 1e-6, "cantilever compliance dJ/d(t_beam) adjoint ({adj_t:e}) vs FD ({fd_t:e}) rel = {rel_t:e}");

    // 2. Material derivatives wrt `E` and `nu`
    let he = 1e-4 * ym0;
    let ce_p = Sms.sensitivity(&mesh, &make_spec(ym0 + he, nu0), &mut ctx, &req_e_only).unwrap().objective;
    let ce_m = Sms.sensitivity(&mesh, &make_spec(ym0 - he, nu0), &mut ctx, &req_e_only).unwrap().objective;
    let fd_e = (ce_p - ce_m) / (2.0 * he);
    let rel_e = (rep.gradients[1][0] - fd_e).abs() / fd_e.abs();
    assert!(
        rel_e < 1e-6,
        "cantilever compliance dJ/dE adjoint ({:e}) vs FD ({fd_e:e}) rel = {rel_e:e}",
        rep.gradients[1][0]
    );

    let hnu = 1e-5;
    let cnu_p = Sms.sensitivity(&mesh, &make_spec(ym0, nu0 + hnu), &mut ctx, &req_e_only).unwrap().objective;
    let cnu_m = Sms.sensitivity(&mesh, &make_spec(ym0, nu0 - hnu), &mut ctx, &req_e_only).unwrap().objective;
    let fd_nu = (cnu_p - cnu_m) / (2.0 * hnu);
    let rel_nu = (rep.gradients[2][0] - fd_nu).abs() / fd_nu.abs();
    assert!(
        rel_nu < 1e-6,
        "cantilever compliance dJ/d(nu) adjoint ({:e}) vs FD ({fd_nu:e}) rel = {rel_nu:e}",
        rep.gradients[2][0]
    );

    // 3. SIMP + Helmholtz-filtered TopologyDensity compliance derivative
    let simplex = mesh.discretise().unwrap();
    let nc = simplex.n_cells();
    let rho0: Vec<f64> = (0..nc).map(|i| 0.55 + 0.30 * ((i as f64 + 1.0) * 0.83).sin()).collect();
    let dir: Vec<f64> = (0..nc).map(|i| 0.4 + 0.5 * ((i as f64 + 2.0) * 0.67).cos().abs()).collect();
    let eval_topo = |rho_vec: Vec<f64>, ctx_ref: &mut Ctx| {
        let r = SensitivityRequest::new(
            SensitivityObjective::Compliance { target: None },
            vec![DesignParam::TopologyDensity {
                domain: "beam".into(),
                filter_radius: 0.14,
                simp_penalty: 3.0,
                ramp_q: None,
                beta: 4.0,
                eta: 0.5,
                val_min: 1e-3,
                val_max: 1.0,
                val_min_im: 0.0,
                val_max_im: 0.0,
                plasmon_damping: 0.0,
                property: None,
                initial_density: None,
                densities: rho_vec,
            }],
        );
        Sms.sensitivity(&mesh, &spec0, ctx_ref, &r).unwrap()
    };
    let rep_topo = eval_topo(rho0.clone(), &mut ctx);
    let hr = 1e-6;
    let rp: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r + hr * d).collect();
    let rm: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r - hr * d).collect();
    let fd_topo = (eval_topo(rp, &mut ctx).objective - eval_topo(rm, &mut ctx).objective) / (2.0 * hr);
    let adj_topo: f64 = rep_topo.gradients[0].iter().zip(&dir).map(|(g, d)| g * d).sum();
    let rel_topo = (adj_topo - fd_topo).abs() / fd_topo.abs();
    assert!(rel_topo < 1e-6, "SIMP compliance topology adjoint ({adj_topo:e}) vs FD ({fd_topo:e}) rel = {rel_topo:e}");
}

#[test]
fn solid_eigenfrequency_and_waveguide_mode_sensitivities_match_fd() {
    // 1. 3D clamped elastic bar eigenfrequency sensitivity wrt length `L` and modulus `E`
    let mut m3 = Model::new_3d();
    m3.param("L", 1.0e-5);
    let rect = m3.rect("sec", (0.0, 0.0), ("L", 2.5e-6));
    let bar = m3.extrude("bar", rect, 0.0, 2.0e-6);
    m3.overlay(&[("bar", bar)]);

    let ev3 = evaluate(&m3.list).unwrap();
    let meshed3 = mesh_geom_3d(ev3, &GeomMeshOpts::new(1.5e-6)).unwrap();
    let mesh3 = meshed3.into_exchange_mesh().unwrap();
    let vel_l = compute_shape_velocity(&m3.list, &meshed3.eval, &mesh3, "L", None).unwrap();
    let resolver3 = TestResolver(&meshed3);

    let make_spec3 = |ym: f64| -> ProblemSpec {
        let mut spec = ProblemSpec::new("sms", 1, "eigenfrequency");
        spec.study.params = json!({ "n_modes": 3, "f_search": 0.0 });
        spec.materials.insert("bar".into(), inline_props(&[("E", ym), ("nu", 0.28), ("rho", 2330.0)]));
        spec.bcs.push(BcSpec {
            on: "plane(x = 0)".into(),
            kind: "fixed".into(),
            params: json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }),
        });
        spec
    };

    let ym0 = 1.65e11;
    let spec3 = make_spec3(ym0);
    let mut ctx3 = Ctx::new(&resolver3);

    let mut req3 = SensitivityRequest::new(
        SensitivityObjective::Eigenfrequency { mode: 0, target_hz: None, scale_hz: None },
        vec![
            DesignParam::ParamVar { name: "L".into(), min: None, max: None },
            DesignParam::MaterialScalar { region: "bar".into(), property: "E".into(), min: None, max: None },
        ],
    );
    req3.shape_velocities.insert("L".into(), vel_l.clone());
    let rep3 = Sms.sensitivity(&mesh3, &spec3, &mut ctx3, &req3).unwrap();

    let req3_e = SensitivityRequest::new(
        SensitivityObjective::Eigenfrequency { mode: 0, target_hz: None, scale_hz: None },
        vec![DesignParam::MaterialScalar { region: "bar".into(), property: "E".into(), min: None, max: None }],
    );

    let hl = 1e-11;
    let m3p = displace_mesh(&mesh3, hl, &vel_l).unwrap();
    let m3m = displace_mesh(&mesh3, -hl, &vel_l).unwrap();
    let mut r_track = req3_e.clone();
    r_track.reference_modes = rep3.mode_basis.clone();
    let fp = Sms.sensitivity(&m3p, &spec3, &mut ctx3, &r_track).unwrap().objective;
    let fm = Sms.sensitivity(&m3m, &spec3, &mut ctx3, &r_track).unwrap().objective;
    let fd_l = (fp - fm) / (2.0 * hl);
    let rel_l = (rep3.gradients[0][0] - fd_l).abs() / fd_l.abs();
    assert!(
        rel_l < 1e-6,
        "3D eigenfrequency df0/dL adjoint ({:e}) vs FD ({fd_l:e}) rel = {rel_l:e}",
        rep3.gradients[0][0]
    );

    let he = 1e6;
    let fe_p = Sms.sensitivity(&mesh3, &make_spec3(ym0 + he), &mut ctx3, &r_track).unwrap().objective;
    let fe_m = Sms.sensitivity(&mesh3, &make_spec3(ym0 - he), &mut ctx3, &r_track).unwrap().objective;
    let fd_e = (fe_p - fe_m) / (2.0 * he);
    let rel_e = (rep3.gradients[1][0] - fd_e).abs() / fd_e.abs();
    assert!(rel_e < 1e-7, "3D eigenfrequency df0/dE rel = {rel_e:e}");

    // 2. 2D acoustic waveguide `modes` sensitivity wrt core width `w`
    let mut mw = Model::new_2d();
    mw.param("w", 1.2e-6);
    let clad = mw.rect("clad", (-2.5e-6, -1.5e-6), (2.5e-6, 1.5e-6));
    let core = mw.rect("core", ("-w/2", -0.3e-6), ("w/2", 0.3e-6));
    mw.overlay(&[("clad", clad), ("core", core)]);
    let evw = evaluate(&mw.list).unwrap();
    let meshedw = mesh_geom_2d(evw, &GeomMeshOpts::new(0.45e-6)).unwrap();
    let meshw = meshedw.into_exchange_mesh().unwrap();
    let vel_w = compute_shape_velocity(&mw.list, &meshedw.eval, &meshw, "w", None).unwrap();
    let resw = TestResolver(&meshedw);

    let mut specw = ProblemSpec::new("sms", 1, "modes");
    specw.study.params = json!({ "n_modes": 3, "f_search": 0.0, "wavevector": 1.5e6 });
    specw.materials.insert("clad".into(), inline_props(&[("E", 7.2e10), ("nu", 0.17), ("rho", 2200.0)]));
    specw.materials.insert("core".into(), inline_props(&[("E", 4.5e10), ("nu", 0.25), ("rho", 3100.0)]));
    specw.bcs.push(BcSpec {
        on: "outer".into(),
        kind: "fixed".into(),
        params: json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }),
    });

    let mut ctxw = Ctx::new(&resw);
    let mut reqw = SensitivityRequest::new(
        SensitivityObjective::Eigenfrequency { mode: 0, target_hz: None, scale_hz: None },
        vec![DesignParam::ParamVar { name: "w".into(), min: None, max: None }],
    );
    reqw.shape_velocities.insert("w".into(), vel_w.clone());
    let repw = Sms.sensitivity(&meshw, &specw, &mut ctxw, &reqw).unwrap();

    let hw = 1e-12;
    let mwp = displace_mesh(&meshw, hw, &vel_w).unwrap();
    let mwm = displace_mesh(&meshw, -hw, &vel_w).unwrap();
    let mut r_mat = SensitivityRequest::new(
        SensitivityObjective::Eigenfrequency { mode: 0, target_hz: None, scale_hz: None },
        vec![DesignParam::MaterialScalar { region: "core".into(), property: "rho".into(), min: None, max: None }],
    );
    r_mat.reference_modes = repw.mode_basis.clone();
    let fwp = Sms.sensitivity(&mwp, &specw, &mut ctxw, &r_mat).unwrap().objective;
    let fwm = Sms.sensitivity(&mwm, &specw, &mut ctxw, &r_mat).unwrap().objective;
    let fd_w = (fwp - fwm) / (2.0 * hw);
    let rel_w = (repw.gradients[0][0] - fd_w).abs() / fd_w.abs();
    assert!(
        rel_w < 1e-6,
        "2D elastic waveguide df/dw adjoint ({:e}) vs FD ({fd_w:e}) rel = {rel_w:e}",
        repw.gradients[0][0]
    );
}
