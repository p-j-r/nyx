// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Finite-difference vs discrete adjoint and modal sensitivity verification for `nyx-ems` (ADR 142):
//! * 2D waveguide `d(n_eff)/dw` (`< 1e-7` rel), `d(n_eff)/d(eps_r)`, SIMP/optical `TopologyDensity`,
//!   and non-symmetric magneto-optic Hall `sheet_2d` (`left_eigenvector` verification).
//! * 3D cavity `d(f_0)/dR` (`< 1e-7` rel) and 3D Floquet-Bloch `nonlinear_eigenvalue_derivative`.
//! * 3D multi-port `driven` `d(|S_21|^2)/d(rho_e)` (`< 1e-6` rel) across BOTH full-order HDM
//!   (`prom_tol = 0`) and reduced Hermite Tangential Galerkin PROM (`prom_tol = 1e-6`).

use femcore::meshing::geomdrive::{block, evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::geomdrive3::mesh_geom_3d;
use femcore::topology::{compute_shape_velocity, displace_mesh};
use geom::model::Model;
use nyx_ems::Ems;
use nyx_physics::{
    BcSpec, Ctx, DesignParam, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, SensitivityObjective,
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
        let b = self.0.boundary(selector);
        if !b.is_empty() {
            return Ok(b);
        }
        if let Some((lhs, rhs)) = selector.split_once('|') {
            let swapped = format!("{}|{}", rhs.trim(), lhs.trim());
            let bs = self.0.boundary(&swapped);
            if !bs.is_empty() {
                return Ok(bs);
            }
        }
        self.0.select(selector).map_err(Error::invalid)
    }
}

struct BlockResolver {
    domain: i32,
    xmin: Vec<i32>,
    xmax: Vec<i32>,
    ymin: Vec<i32>,
    ymax: Vec<i32>,
    zmin: Vec<i32>,
    zmax: Vec<i32>,
}

impl LabelResolver for BlockResolver {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        match region {
            "wg" | "cell" => Ok(vec![self.domain]),
            other => Err(Error::invalid(format!("unknown domain '{other}'"))),
        }
    }

    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        match selector {
            "xmin" => Ok(self.xmin.clone()),
            "xmax" => Ok(self.xmax.clone()),
            "ymin" => Ok(self.ymin.clone()),
            "ymax" => Ok(self.ymax.clone()),
            "zmin" => Ok(self.zmin.clone()),
            "zmax" => Ok(self.zmax.clone()),
            other => Err(Error::invalid(format!("unknown boundary '{other}'"))),
        }
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
fn waveguide_neff_and_hall_sheet_sensitivities_match_fd_below_1e_7() {
    // 2D SOI strip waveguide: d(n_eff)/dw, d(n_eff)/d(eps_r), and TopologyDensity
    let mut m = Model::new_2d();
    m.param("w", 500e-9);
    let clad = m.rect("clad", (-1.5e-6, -1.0e-6), (1.5e-6, 1.0e-6));
    let core = m.rect("core", ("-w/2", -110e-9), ("w/2", 110e-9));
    m.overlay(&[("oxide", clad), ("silicon", core)]);

    let ev0 = evaluate(&m.list).unwrap();
    let meshed = mesh_geom_2d(ev0, &GeomMeshOpts::new(160e-9)).unwrap();
    let mesh = meshed.into_exchange_mesh().unwrap();
    let vel_w = compute_shape_velocity(&m.list, &meshed.eval, &mesh, "w", None).unwrap();
    let resolver = TestResolver(&meshed);

    let make_spec = |eps_si: f64| -> ProblemSpec {
        let mut spec = ProblemSpec::new("ems", 1, "modes");
        spec.study.params = json!({ "freq": 1.934e14, "n_modes": 2 });
        spec.materials.insert("oxide".into(), inline_props(&[("eps_r", 2.085)]));
        spec.materials.insert("silicon".into(), inline_props(&[("eps_r", eps_si)]));
        spec
    };

    let spec0 = make_spec(12.08);
    let mut ctx = Ctx::new(&resolver);

    let mut req = SensitivityRequest::new(
        SensitivityObjective::EffectiveIndex { mode: 0, target_neff: None },
        vec![
            DesignParam::ParamVar { name: "w".into(), min: None, max: None },
            DesignParam::MaterialScalar { region: "silicon".into(), property: "eps_r".into(), min: None, max: None },
        ],
    );
    req.shape_velocities.insert("w".into(), vel_w.clone());
    let rep = Ems.sensitivity(&mesh, &spec0, &mut ctx, &req).unwrap();

    let mut req_eps = SensitivityRequest::new(
        SensitivityObjective::EffectiveIndex { mode: 0, target_neff: None },
        vec![DesignParam::MaterialScalar { region: "silicon".into(), property: "eps_r".into(), min: None, max: None }],
    );
    req_eps.reference_modes = rep.mode_basis.clone();

    // Central FD wrt `w` on fixed-topology displaced mesh
    let hw = 5e-11;
    let mp = displace_mesh(&mesh, hw, &vel_w).unwrap();
    let mm = displace_mesh(&mesh, -hw, &vel_w).unwrap();
    let np = Ems.sensitivity(&mp, &spec0, &mut ctx, &req_eps).unwrap().objective;
    let nm = Ems.sensitivity(&mm, &spec0, &mut ctx, &req_eps).unwrap().objective;
    let fd_w = (np - nm) / (2.0 * hw);
    let adj_w = rep.gradients[0][0];
    let rel_w = (adj_w - fd_w).abs() / fd_w.abs();
    assert!(rel_w < 1e-7, "2D waveguide d(neff)/dw adjoint ({adj_w:e}) vs FD ({fd_w:e}) rel = {rel_w:e}");

    // Central FD wrt `silicon.eps_r`
    let he = 1e-5;
    let ne_p = Ems.sensitivity(&mesh, &make_spec(12.08 + he), &mut ctx, &req_eps).unwrap().objective;
    let ne_m = Ems.sensitivity(&mesh, &make_spec(12.08 - he), &mut ctx, &req_eps).unwrap().objective;
    let fd_eps = (ne_p - ne_m) / (2.0 * he);
    let adj_eps = rep.gradients[1][0];
    let rel_eps = (adj_eps - fd_eps).abs() / fd_eps.abs();
    assert!(rel_eps < 1e-8, "2D waveguide d(neff)/d(eps_r) adjoint ({adj_eps:e}) vs FD ({fd_eps:e}) rel = {rel_eps:e}");

    // 2D waveguide `TopologyDensity` directional-derivative FD check (optical interpolation + Helmholtz pullback)
    let core_doms = resolver.domains("silicon").unwrap();
    let simplex2 = mesh.discretise().unwrap();
    let n_core = simplex2.cell_labels.iter().filter(|l| core_doms.contains(l)).count();
    let rho0: Vec<f64> = (0..n_core).map(|i| 0.45 + 0.35 * ((i as f64 + 1.0) * 0.77).sin()).collect();
    let dir: Vec<f64> = (0..n_core).map(|i| 0.3 + 0.5 * ((i as f64 + 2.0) * 0.61).cos().abs()).collect();
    let make_topo_req = |rv: Vec<f64>, ref_m| {
        let mut rq = SensitivityRequest::new(
            SensitivityObjective::EffectiveIndex { mode: 0, target_neff: None },
            vec![DesignParam::TopologyDensity {
                domain: "silicon".into(),
                property: None,
                filter_radius: 120e-9,
                beta: 4.0,
                eta: 0.5,
                simp_penalty: 1.0,
                ramp_q: None,
                val_min: 2.085,
                val_max: 12.08,
                val_min_im: 0.0,
                val_max_im: 0.0,
                plasmon_damping: 0.0,
                initial_density: None,
                densities: rv,
            }],
        );
        rq.reference_modes = ref_m;
        rq
    };
    let rep_wg_topo = Ems.sensitivity(&mesh, &spec0, &mut ctx, &make_topo_req(rho0.clone(), None)).unwrap();
    let hr = 1e-6;
    let rp: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r + hr * d).collect();
    let rm: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r - hr * d).collect();
    let ntp =
        Ems.sensitivity(&mesh, &spec0, &mut ctx, &make_topo_req(rp, rep_wg_topo.mode_basis.clone())).unwrap().objective;
    let ntm =
        Ems.sensitivity(&mesh, &spec0, &mut ctx, &make_topo_req(rm, rep_wg_topo.mode_basis.clone())).unwrap().objective;
    let fd_wg_topo = (ntp - ntm) / (2.0 * hr);
    let adj_wg_topo: f64 = rep_wg_topo.gradients[0].iter().zip(&dir).map(|(g, d)| g * d).sum();
    let rel_wg_topo = (adj_wg_topo - fd_wg_topo).abs() / fd_wg_topo.abs();
    assert!(
        rel_wg_topo < 1e-6,
        "2D waveguide TopologyDensity directional d(neff)/d(rho) adjoint ({adj_wg_topo:e}) vs FD ({fd_wg_topo:e}) rel = {rel_wg_topo:e}"
    );

    // Non-symmetric Hall-biased 2D graphene sheet (`b_0 = 6.0` T -> `A^T != A`, left != right eigenvector!)
    let mut spec_hall = make_spec(12.08);
    spec_hall.bcs.push(BcSpec {
        on: "oxide|silicon".into(),
        kind: "sheet_2d".into(),
        params: json!({
            "model": "kubo_graphene",
            "mu_c_eV": 0.45,
            "tau_s": 0.4e-12,
            "temperature": 300.0,
            "b_0": 6.0
        }),
    });
    let req_hall = SensitivityRequest::new(
        SensitivityObjective::EffectiveIndex { mode: 0, target_neff: None },
        vec![DesignParam::MaterialScalar { region: "silicon".into(), property: "eps_r".into(), min: None, max: None }],
    );
    let rep_hall = Ems.sensitivity(&mesh, &spec_hall, &mut ctx, &req_hall).unwrap();
    let mut req_hall_tr = req_hall.clone();
    req_hall_tr.reference_modes = rep_hall.mode_basis.clone();
    let mut sp_p = spec_hall.clone();
    sp_p.materials.insert("silicon".into(), inline_props(&[("eps_r", 12.08 + he)]));
    let mut sp_m = spec_hall.clone();
    sp_m.materials.insert("silicon".into(), inline_props(&[("eps_r", 12.08 - he)]));
    let nh_p = Ems.sensitivity(&mesh, &sp_p, &mut ctx, &req_hall_tr).unwrap().objective;
    let nh_m = Ems.sensitivity(&mesh, &sp_m, &mut ctx, &req_hall_tr).unwrap().objective;
    let fd_hall = (nh_p - nh_m) / (2.0 * he);
    let rel_hall = (rep_hall.gradients[0][0] - fd_hall).abs() / fd_hall.abs();
    assert!(
        rel_hall < 1e-6,
        "non-symmetric Hall sheet waveguide d(neff)/d(eps_r) left-eigenvector adjoint ({:e}) vs FD ({fd_hall:e}) rel = {rel_hall:e}",
        rep_hall.gradients[0][0]
    );
}

#[test]
fn cavity_3d_f0_vs_radius_and_bloch_sensitivities_match_fd_below_1e_7() {
    // 1. 3D cylindrical PEC cavity: d(f0)/dR wrt CAD radius `R`
    let mut m3 = Model::new_3d();
    m3.param("R", 0.010);
    let circ = m3.circle("c", (0.0, 0.0), "R");
    let cyl = m3.extrude("cyl", circ, 0.0, 0.012);
    m3.overlay(&[("cav", cyl)]);

    let ev3 = evaluate(&m3.list).unwrap();
    let meshed3 = mesh_geom_3d(ev3, &GeomMeshOpts::new(0.0045)).unwrap();
    let mesh3 = meshed3.into_exchange_mesh().unwrap();
    let vel_r = compute_shape_velocity(&m3.list, &meshed3.eval, &mesh3, "R", None).unwrap();
    let res3 = TestResolver(&meshed3);

    let mut spec_cav = ProblemSpec::new("ems", 1, "cavity");
    spec_cav.study.params = json!({ "n_modes": 3, "f_search": 10.0e9 });
    spec_cav.materials.insert("cav".into(), inline_props(&[("eps_r", 1.0)]));

    let mut ctx3 = Ctx::new(&res3);
    let mut req_cav = SensitivityRequest::new(
        SensitivityObjective::Eigenfrequency { mode: 0, target_hz: None, scale_hz: None },
        vec![DesignParam::ParamVar { name: "R".into(), min: None, max: None }],
    );
    req_cav.shape_velocities.insert("R".into(), vel_r.clone());
    let rep_cav = Ems.sensitivity(&mesh3, &spec_cav, &mut ctx3, &req_cav).unwrap();

    let mut req_cav_eps = SensitivityRequest::new(
        SensitivityObjective::Eigenfrequency { mode: 0, target_hz: None, scale_hz: None },
        vec![DesignParam::MaterialScalar { region: "cav".into(), property: "eps_r".into(), min: None, max: None }],
    );
    req_cav_eps.reference_modes = rep_cav.mode_basis.clone();

    let hr = 1e-8;
    let m3p = displace_mesh(&mesh3, hr, &vel_r).unwrap();
    let m3m = displace_mesh(&mesh3, -hr, &vel_r).unwrap();
    let fp = Ems.sensitivity(&m3p, &spec_cav, &mut ctx3, &req_cav_eps).unwrap().objective;
    let fm = Ems.sensitivity(&m3m, &spec_cav, &mut ctx3, &req_cav_eps).unwrap().objective;
    let fd_r = (fp - fm) / (2.0 * hr);
    let adj_r = rep_cav.gradients[0][0];
    let rel_r = (adj_r - fd_r).abs() / fd_r.abs();
    assert!(rel_r < 1e-7, "3D cavity df0/dR adjoint ({adj_r:e}) vs FD ({fd_r:e}) rel = {rel_r:e}");

    // 2. 3D Floquet-Bloch unit cell: d(n_B)/d(eps_r) via `nonlinear_eigenvalue_derivative`
    let blk = block(0.4e-6, 0.4e-6, 0.4e-6, 0.2e-6, 2).unwrap();
    let mesh_bloch = Mesh::from_solver_mesh(&blk.mesh).unwrap();
    let res_bloch = BlockResolver {
        domain: blk.domain,
        xmin: blk.xmin,
        xmax: blk.xmax,
        ymin: blk.ymin,
        ymax: blk.ymax,
        zmin: blk.zmin,
        zmax: blk.zmax,
    };
    let make_bloch_spec = |eps_val: f64| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "bloch");
        s.study.params = json!({ "freq": 1.934e14, "n_target": 2.0, "n_modes": 2 });
        s.materials.insert("cell".into(), inline_props(&[("eps_r", eps_val)]));
        s.bcs.push(BcSpec {
            on: "ymin".into(),
            kind: "periodic".into(),
            params: json!({ "dest": "ymax", "phase": 0.0 }),
        });
        s
    };
    let mut ctx_b = Ctx::new(&res_bloch);
    let req_b = SensitivityRequest::new(
        SensitivityObjective::EffectiveIndex { mode: 0, target_neff: None },
        vec![DesignParam::MaterialScalar { region: "cell".into(), property: "eps_r".into(), min: None, max: None }],
    );
    let rep_b = Ems.sensitivity(&mesh_bloch, &make_bloch_spec(4.0), &mut ctx_b, &req_b).unwrap();
    let mut req_b_tr = req_b.clone();
    req_b_tr.reference_modes = rep_b.mode_basis.clone();
    let hb = 1e-5;
    let nb_p = Ems.sensitivity(&mesh_bloch, &make_bloch_spec(4.0 + hb), &mut ctx_b, &req_b_tr).unwrap().objective;
    let nb_m = Ems.sensitivity(&mesh_bloch, &make_bloch_spec(4.0 - hb), &mut ctx_b, &req_b_tr).unwrap().objective;
    let fd_b = (nb_p - nb_m) / (2.0 * hb);
    let adj_b = rep_b.gradients[0][0];
    let rel_b = (adj_b - fd_b).abs() / fd_b.abs();
    assert!(rel_b < 1e-7, "3D Bloch d(n_B)/d(eps_r) adjoint ({adj_b:e}) vs FD ({fd_b:e}) rel = {rel_b:e}");
}

#[test]
fn driven_3d_s21_adjoint_matches_fd_for_both_hdm_and_hermite_prom() {
    // 3D WR-style rectangular waveguide [0, 0.02] x [0, 0.03] x [0, 0.01] with waveports at ymin and ymax
    let blk = block(0.02, 0.03, 0.01, 0.0075, 4).unwrap();
    let nc = blk.mesh.n_cells();
    let mesh_wg = Mesh::from_solver_mesh(&blk.mesh).unwrap();
    let res_wg = BlockResolver {
        domain: blk.domain,
        xmin: blk.xmin,
        xmax: blk.xmax,
        ymin: blk.ymin,
        ymax: blk.ymax,
        zmin: blk.zmin,
        zmax: blk.zmax,
    };

    let make_driven_spec = |prom_tol: f64| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({
            "f_min": 9.5e9,
            "f_max": 11.5e9,
            "n_points": 5,
            "prom_tol": prom_tol,
            "max_hdm_solves": 8
        });
        s.materials.insert("wg".into(), inline_props(&[("eps_r", 1.0)]));
        s.bcs.push(BcSpec { on: "ymin".into(), kind: "waveport".into(), params: json!({ "name": "p1", "mode": 1 }) });
        s.bcs.push(BcSpec { on: "ymax".into(), kind: "waveport".into(), params: json!({ "name": "p2", "mode": 1 }) });
        s
    };

    // Spatially varying interior cellwise density (keeping cells touching ymin/ymax ports uniform at 0.0 so port eigenmodes stay invariant)
    let mut port_verts = vec![false; blk.mesh.n_points()];
    for f in 0..blk.mesh.n_facets() {
        let fl = blk.mesh.facet_labels[f];
        if res_wg.ymin.contains(&fl) || res_wg.ymax.contains(&fl) {
            for &vid in blk.mesh.facet(f) {
                port_verts[vid as usize] = true;
            }
        }
    }
    let mut rho0 = vec![0.0f64; nc];
    let mut dir = vec![0.0f64; nc];
    for c in 0..nc {
        let touches_port = blk.mesh.cell(c).iter().any(|&v| port_verts[v as usize]);
        if !touches_port {
            rho0[c] = 0.25 + 0.55 * ((c as f64 + 1.0) * 0.91).sin().abs();
            dir[c] = 0.3 + 0.7 * ((c as f64 + 2.0) * 0.57).cos().abs();
        }
    }
    assert!(dir.iter().any(|&d| d > 0.0), "interior design cells must be non-empty");

    let eval_driven = |prom_tol: f64, rho_vec: Vec<f64>| {
        let mut ctx = Ctx::new(&res_wg);
        let spec = make_driven_spec(prom_tol);
        let req = SensitivityRequest::new(
            SensitivityObjective::Transmission { out_port: 1, in_port: 0, freq_index: Some(2), target: None },
            vec![DesignParam::TopologyDensity {
                domain: "wg".into(),
                filter_radius: 0.0,
                simp_penalty: 2.0,
                ramp_q: None,
                beta: 1.0,
                eta: 0.5,
                val_min: 1.0,
                val_max: 3.5,
                val_min_im: 0.0,
                val_max_im: 0.0,
                plasmon_damping: 0.0,
                property: None,
                initial_density: None,
                densities: rho_vec,
            }],
        );
        Ems.sensitivity(&mesh_wg, &spec, &mut ctx, &req).unwrap()
    };

    // 1. Full-Order HDM (`prom_tol = 0.0`): Adjoint vs Central FD
    let rep_hdm = eval_driven(0.0, rho0.clone());
    let hd = 1e-6;
    let rp: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r + hd * d).collect();
    let rm: Vec<f64> = rho0.iter().zip(&dir).map(|(r, d)| r - hd * d).collect();
    let jp = eval_driven(0.0, rp.clone()).objective;
    let jm = eval_driven(0.0, rm.clone()).objective;
    let fd_hdm = (jp - jm) / (2.0 * hd);
    let adj_hdm: f64 = rep_hdm.gradients[0].iter().zip(&dir).map(|(g, d)| g * d).sum();
    let rel_hdm = (adj_hdm - fd_hdm).abs() / fd_hdm.abs();
    assert!(
        rel_hdm < 1e-6,
        "driven HDM d(|S_21|^2)/d(rho_e) adjoint ({adj_hdm:e}) vs FD ({fd_hdm:e}) rel = {rel_hdm:e}"
    );

    // 2. Reduced Hermite Tangential Galerkin PROM (`prom_tol = 1e-6`): Adjoint vs FD and vs HDM
    let rep_prom = eval_driven(1e-6, rho0.clone());
    let adj_prom: f64 = rep_prom.gradients[0].iter().zip(&dir).map(|(g, d)| g * d).sum();
    let rel_prom_vs_fd = (adj_prom - fd_hdm).abs() / fd_hdm.abs();
    assert!(
        rel_prom_vs_fd < 1e-6,
        "driven Hermite PROM d(|S_21|^2)/d(rho_e) adjoint ({adj_prom:e}) vs FD ({fd_hdm:e}) rel = {rel_prom_vs_fd:e}"
    );
}
