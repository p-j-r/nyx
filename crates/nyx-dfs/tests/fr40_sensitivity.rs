// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Finite-difference vs discrete adjoint sensitivity verification for `nyx-dfs` (ADR 142):
//! electrostatics (`Capacitance` and `Compliance`) and stationary heat transfer (`Compliance`
//! and `TargetTemperature` mean / smooth `L_p` peak) across CAD shape parameters, scalar materials,
//! and Helmholtz-filtered SIMP/RAMP topology densities.

use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::topology::{compute_shape_velocity, displace_mesh};
use geom::model::{Model, RegionKind};
use nyx_dfs::Dfs;
use nyx_physics::{
    BcSpec, Ctx, DesignParam, Error, LabelResolver, MaterialAssign, Physics, ProblemSpec, SensitivityObjective,
    SensitivityRequest,
};
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
        if matches!(selector, "t1_hole" | "t2_hole" | "box_wall") {
            let m = &self.0.mesh;
            let ext = m.exterior_facets();
            let mut labs = Vec::new();
            for (f, &is_ext) in ext.iter().enumerate() {
                if !is_ext {
                    continue;
                }
                let fv = m.facet(f);
                let p0 = m.point(fv[0] as usize);
                let p1 = m.point(fv[1] as usize);
                let xm = 0.5 * (p0[0] + p1[0]);
                let zm = 0.5 * (p0[1] + p1[1]);
                let hit = match selector {
                    "t1_hole" => xm > -1.1 && xm < 0.0 && zm.abs() < 0.5,
                    "t2_hole" => xm > 0.0 && xm < 1.1 && zm.abs() < 0.5,
                    "box_wall" => (xm.abs() - 1.5).abs() < 1e-6 || (zm.abs() - 1.0).abs() < 1e-6,
                    _ => false,
                };
                if hit {
                    let lab = m.facet_labels[f];
                    if !labs.contains(&lab) {
                        labs.push(lab);
                    }
                }
            }
            labs.sort_unstable();
            if labs.is_empty() {
                return Err(Error::invalid(format!("empty hole/wall boundary '{selector}'")));
            }
            return Ok(labs);
        }
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
fn electrostatic_capacitance_adjoint_matches_fd_below_1e_7() {
    // 2D coupled strip capacitor in layered dielectric box
    let mut m = Model::new_2d();
    m.param("w", 0.4);
    let box_dom = m.rect("box", (-1.5, -1.0), (1.5, 1.0));
    let sub = m.rect("sub", (-1.5, -1.0), (1.5, 0.0));
    let sig1 = m.rect("sig1", ("-0.5 - w/2", 0.0), ("-0.5 + w/2", 0.15));
    let sig2 = m.rect("sig2", ("0.5 - w/2", 0.0), ("0.5 + w/2", 0.15));
    m.overlay_full(
        &[
            ("air", box_dom, RegionKind::Material, None),
            ("sub", sub, RegionKind::Material, None),
            ("sig1", sig1, RegionKind::Hole, None),
            ("sig2", sig2, RegionKind::Hole, None),
        ],
        &[],
    );

    let ev0 = evaluate(&m.list).unwrap();
    let meshed = mesh_geom_2d(ev0, &GeomMeshOpts::new(0.18)).unwrap();
    let mesh = meshed.into_exchange_mesh().unwrap();
    let vel_w = compute_shape_velocity(&m.list, &meshed.eval, &mesh, "w", None).unwrap();
    let resolver = TestResolver(&meshed);

    let make_spec = |eps_sub: f64| -> ProblemSpec {
        let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
        spec.materials.insert("air".into(), inline_props(&[("eps_r", 1.0)]));
        spec.materials.insert("sub".into(), inline_props(&[("eps_r", eps_sub)]));
        spec.bcs.push(BcSpec { on: "t1_hole".into(), kind: "terminal".into(), params: json!({ "name": "t1" }) });
        spec.bcs.push(BcSpec { on: "t2_hole".into(), kind: "terminal".into(), params: json!({ "name": "t2" }) });
        spec.bcs.push(BcSpec { on: "box_wall".into(), kind: "ground".into(), params: serde_json::Value::Null });
        spec
    };

    let spec0 = make_spec(9.8);
    let mut ctx = Ctx::new(&resolver);

    // 1. Self-capacitance C_00 sensitivity wrt CAD parameter `w` and substrate `eps_r`
    let mut req = SensitivityRequest::new(
        SensitivityObjective::Capacitance { term_i: 0, term_j: 0, target_f: None },
        vec![
            DesignParam::ParamVar { name: "w".into(), min: None, max: None },
            DesignParam::MaterialScalar { region: "sub".into(), property: "eps_r".into(), min: None, max: None },
        ],
    );
    req.shape_velocities.insert("w".into(), vel_w.clone());
    let rep = Dfs.sensitivity(&mesh, &spec0, &mut ctx, &req).unwrap();

    // Central FD wrt `w` on fixed-topology displaced mesh
    let hw = 1e-6;
    let mp = displace_mesh(&mesh, hw, &vel_w).unwrap();
    let mm = displace_mesh(&mesh, -hw, &vel_w).unwrap();
    let req_eps_only = SensitivityRequest::new(
        SensitivityObjective::Capacitance { term_i: 0, term_j: 0, target_f: None },
        vec![DesignParam::MaterialScalar { region: "sub".into(), property: "eps_r".into(), min: None, max: None }],
    );
    let cp = Dfs.sensitivity(&mp, &spec0, &mut ctx, &req_eps_only).unwrap().objective;
    let cm = Dfs.sensitivity(&mm, &spec0, &mut ctx, &req_eps_only).unwrap().objective;
    let fd_w = (cp - cm) / (2.0 * hw);
    let adj_w = rep.gradients[0][0];
    let rel_w = (adj_w - fd_w).abs() / fd_w.abs();
    assert!(rel_w < 1e-7, "d(C_00)/dw adjoint ({adj_w:e}) vs FD ({fd_w:e}) relative error = {rel_w:e}");

    // Central FD wrt `eps_r`
    let he = 1e-5;
    let ce_p = Dfs.sensitivity(&mesh, &make_spec(9.8 + he), &mut ctx, &req_eps_only).unwrap().objective;
    let ce_m = Dfs.sensitivity(&mesh, &make_spec(9.8 - he), &mut ctx, &req_eps_only).unwrap().objective;
    let fd_eps = (ce_p - ce_m) / (2.0 * he);
    let adj_eps = rep.gradients[1][0];
    let rel_eps = (adj_eps - fd_eps).abs() / fd_eps.abs();
    assert!(rel_eps < 1e-8, "d(C_00)/d(eps_r) adjoint ({adj_eps:e}) vs FD ({fd_eps:e}) relative error = {rel_eps:e}");

    // 2. Mutual capacitance C_01 wrt Helmholtz-filtered TopologyDensity in `sub`
    let sub_doms = resolver.domains("sub").unwrap();
    let simplex = mesh.discretise().unwrap();
    let nc_sub = simplex.cell_labels.iter().filter(|l| sub_doms.contains(l)).count();
    let rho0: Vec<f64> = (0..nc_sub).map(|i| 0.35 + 0.45 * ((i as f64 + 1.0) * 0.71).sin().abs()).collect();
    let dir: Vec<f64> = (0..nc_sub).map(|i| 0.4 + 0.6 * ((i as f64 + 2.0) * 1.13).cos().abs()).collect();

    let eval_topo = |rho_vec: Vec<f64>, ctx_ref: &mut Ctx| {
        let r = SensitivityRequest::new(
            SensitivityObjective::Capacitance { term_i: 0, term_j: 1, target_f: None },
            vec![DesignParam::TopologyDensity {
                domain: "sub".into(),
                filter_radius: 0.2,
                simp_penalty: 3.0,
                ramp_q: None,
                beta: 4.0,
                eta: 0.5,
                val_min: 1.0,
                val_max: 11.7,
                val_min_im: 0.0,
                val_max_im: 0.0,
                plasmon_damping: 0.0,
                property: None,
                initial_density: None,
                densities: rho_vec,
            }],
        );
        Dfs.sensitivity(&mesh, &spec0, ctx_ref, &r).unwrap()
    };

    let rep_topo = eval_topo(rho0.clone(), &mut ctx);
    let ht = 1e-6;
    let rp: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r + ht * d).collect();
    let rm: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r - ht * d).collect();
    let jp = eval_topo(rp, &mut ctx).objective;
    let jm = eval_topo(rm, &mut ctx).objective;
    let fd_topo = (jp - jm) / (2.0 * ht);
    let adj_topo: f64 = rep_topo.gradients[0].iter().zip(&dir).map(|(g, d)| g * d).sum();
    let rel_topo = (adj_topo - fd_topo).abs() / fd_topo.abs();
    assert!(rel_topo < 1e-7, "d(C_01)/d(rho) adjoint ({adj_topo:e}) vs FD ({fd_topo:e}) relative error = {rel_topo:e}");
}

#[test]
fn heat_stationary_compliance_and_temperature_adjoints_match_fd() {
    let mut m = Model::new_2d();
    m.param("h_pad", 0.25);
    let sub = m.rect("sub", (-1.0, -0.8), (1.0, 0.0));
    let pad = m.rect("pad", (-0.3, 0.0), (0.3, "h_pad"));
    m.overlay(&[("sub", sub), ("pad", pad)]);

    let ev0 = evaluate(&m.list).unwrap();
    let meshed = mesh_geom_2d(ev0, &GeomMeshOpts::new(0.16)).unwrap();
    let mesh = meshed.into_exchange_mesh().unwrap();
    let vel_h = compute_shape_velocity(&m.list, &meshed.eval, &mesh, "h_pad", None).unwrap();
    let resolver = TestResolver(&meshed);

    let make_spec = |k_sub: f64, q_pad: f64| -> ProblemSpec {
        let mut spec = ProblemSpec::new("dfs", 2, "heat_stationary");
        spec.materials.insert("sub".into(), inline_props(&[("k", k_sub), ("Q", 0.0)]));
        spec.materials.insert("pad".into(), inline_props(&[("k", 40.0), ("Q", q_pad)]));
        spec.bcs.push(BcSpec {
            on: "plane(z = -0.8)".into(),
            kind: "temperature".into(),
            params: json!({ "value": 300.0 }),
        });
        spec.bcs.push(BcSpec {
            on: "pad.top".into(),
            kind: "convection".into(),
            params: json!({ "h": 500.0, "t_ext": 300.0 }),
        });
        spec
    };

    let spec0 = make_spec(130.0, 2.5e6);
    let mut ctx = Ctx::new(&resolver);

    // 1. Thermal Compliance wrt shape `h_pad` and `sub.k`
    let mut req_comp = SensitivityRequest::new(
        SensitivityObjective::Compliance { target: None },
        vec![
            DesignParam::ParamVar { name: "h_pad".into(), min: None, max: None },
            DesignParam::MaterialScalar { region: "sub".into(), property: "k".into(), min: None, max: None },
        ],
    );
    req_comp.shape_velocities.insert("h_pad".into(), vel_h.clone());
    let rep_comp = Dfs.sensitivity(&mesh, &spec0, &mut ctx, &req_comp).unwrap();

    let hh = 1e-6;
    let mp = displace_mesh(&mesh, hh, &vel_h).unwrap();
    let mm = displace_mesh(&mesh, -hh, &vel_h).unwrap();
    let req_k_only = SensitivityRequest::new(
        SensitivityObjective::Compliance { target: None },
        vec![DesignParam::MaterialScalar { region: "sub".into(), property: "k".into(), min: None, max: None }],
    );
    let cp = Dfs.sensitivity(&mp, &spec0, &mut ctx, &req_k_only).unwrap().objective;
    let cm = Dfs.sensitivity(&mm, &spec0, &mut ctx, &req_k_only).unwrap().objective;
    let fd_h = (cp - cm) / (2.0 * hh);
    let rel_h = (rep_comp.gradients[0][0] - fd_h).abs() / fd_h.abs();
    assert!(rel_h < 1e-6, "heat compliance dJ/d(h_pad) rel = {rel_h:e}");

    let hk = 1e-3;
    let ck_p = Dfs.sensitivity(&mesh, &make_spec(130.0 + hk, 2.5e6), &mut ctx, &req_k_only).unwrap().objective;
    let ck_m = Dfs.sensitivity(&mesh, &make_spec(130.0 - hk, 2.5e6), &mut ctx, &req_k_only).unwrap().objective;
    let fd_k = (ck_p - ck_m) / (2.0 * hk);
    let rel_k = (rep_comp.gradients[1][0] - fd_k).abs() / fd_k.abs();
    assert!(rel_k < 1e-8, "heat compliance dJ/dk rel = {rel_k:e}");

    // 2. Smooth L_p peak temperature (`p_norm = 6.0`) in `pad` wrt `TopologyDensity` in `sub`
    let sub_doms = resolver.domains("sub").unwrap();
    let simplex = mesh.discretise().unwrap();
    let nc_sub = simplex.cell_labels.iter().filter(|l| sub_doms.contains(l)).count();
    let rho0: Vec<f64> = (0..nc_sub).map(|i| 0.4 + 0.3 * ((i as f64) * 0.9).sin()).collect();
    let dir: Vec<f64> = (0..nc_sub).map(|i| 0.5 + 0.4 * ((i as f64) * 0.5).cos()).collect();

    let eval_peak = |rho_vec: Vec<f64>, ctx_ref: &mut Ctx| {
        let r = SensitivityRequest::new(
            SensitivityObjective::TargetTemperature { domain: Some("pad".into()), target_k: None, p_norm: Some(6.0) },
            vec![DesignParam::TopologyDensity {
                domain: "sub".into(),
                filter_radius: 0.15,
                simp_penalty: 3.0,
                ramp_q: Some(3.0),
                beta: 3.0,
                eta: 0.5,
                val_min: 1.5,
                val_max: 150.0,
                val_min_im: 0.0,
                val_max_im: 0.0,
                plasmon_damping: 0.0,
                property: Some("k".into()),
                initial_density: None,
                densities: rho_vec,
            }],
        );
        Dfs.sensitivity(&mesh, &spec0, ctx_ref, &r).unwrap()
    };

    let rep_peak = eval_peak(rho0.clone(), &mut ctx);
    let ht = 1e-6;
    let rp: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r + ht * d).collect();
    let rm: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r - ht * d).collect();
    let fd_peak = (eval_peak(rp, &mut ctx).objective - eval_peak(rm, &mut ctx).objective) / (2.0 * ht);
    let adj_peak: f64 = rep_peak.gradients[0].iter().zip(&dir).map(|(g, d)| g * d).sum();
    let rel_peak = (adj_peak - fd_peak).abs() / fd_peak.abs();
    assert!(
        rel_peak < 1e-6,
        "L_p peak temperature topology adjoint ({adj_peak:e}) vs FD ({fd_peak:e}) rel = {rel_peak:e}"
    );
}
