// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification benchmarks for Phase 1D (`nyx-sms` native mixed-cell assembly, curvilinear PMLs,
//! pre-solve Christoffel backward-wave stabilization, signed quality factor `Q = Re f / (2 Im f)`
//! without `.abs()` masking, outer PML clamping, and locking-free Static Infinite Elements).

mod common;

use femcore::adapt::work::simplex_to_adapt;
use femcore::material::{stiffness_global, Anisotropy, Stiffness};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::SimplexMesh;
use femcore::meshing::physics_size::backward_wave_mpml_ratio;
use femcore::meshing::pml_box::{build_pml, PmlSpec};
use femcore::progress::Progress;
use femcore::props::{Periodic, Pml, PmlCoord, PmlDir, PmlStabilization};
use nyx_physics::{BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, Scalar};
use nyx_sms::{run, solve_solid, Bc, BcKind, IeProfile, InfiniteElement, Plane, Problem, Sms, Study};
use serde_json::json;
use std::collections::HashMap;

fn num_scalar(r: &nyx_physics::SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(x))) => *x,
        other => panic!("expected numeric scalar '{name}', got {other:?} in {:?}", r.scalars),
    }
}

/// Builds a conforming 2D `Quad` strip `[0, lx_phys + lx_pml] x [0, lz]` with domain `1` in
/// `[0, lx_phys]` and domain `2` in `[lx_phys, lx_phys + lx_pml]`.
/// Boundary labels: `1` (`x = 0`), `2` (`z = 0`), `3` (`z = lz`), `4` (`x = lx_phys + lx_pml`).
fn quad_strip_with_pml(lx_phys: f64, lx_pml: f64, lz: f64, nx_phys: usize, nx_pml: usize, nz: usize) -> Mesh {
    let nx = nx_phys + nx_pml;
    let mut cm = CellMesh::<2>::new();
    let vid = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    for j in 0..=nz {
        let z = lz * (j as f64) / (nz as f64);
        for i in 0..=nx {
            let x = if i <= nx_phys {
                lx_phys * (i as f64) / (nx_phys as f64)
            } else {
                lx_phys + lx_pml * ((i - nx_phys) as f64) / (nx_pml as f64)
            };
            cm.points.push([x, z]);
        }
    }
    for j in 0..nz {
        for i in 0..nx {
            let dom = if i < nx_phys { 1 } else { 2 };
            cm.push_cell(CellType::Quad, &[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)], dom);
        }
    }
    for j in 0..nz {
        cm.push_facet(&[vid(0, j), vid(0, j + 1)], 1);
        cm.push_facet(&[vid(nx, j), vid(nx, j + 1)], 4);
    }
    for i in 0..nx {
        cm.push_facet(&[vid(i, 0), vid(i + 1, 0)], 2);
        cm.push_facet(&[vid(i, nz), vid(i + 1, nz)], 3);
    }
    Mesh::D2(cm)
}

/// Builds a conforming 1st-octant (`x >= 0, y >= 0, z >= 0`) cubed-sphere shell mesh from inner
/// radius `r0` to `r_in` (`n_phys` radial layers, domain `1`) and from `r_in` to `r_out`
/// (`n_ie` radial layers, domain `2`), using either native `Hex` (`use_prisms = false`) or native
/// `Prism` (`use_prisms = true`) elements, with optional quadratic spherical arc `edge_nodes`.
///
/// Boundary facet labels:
/// - `1`: symmetry plane `x = 0`
/// - `2`: symmetry plane `y = 0`
/// - `3`: symmetry plane `z = 0`
/// - `4`: inner spherical cavity `R = r0`
/// - `5`: outer spherical boundary `R = r_out`
fn octant_spherical_shell_mesh(
    n_ang: usize,
    radii_spec: [f64; 3],
    n_phys: usize,
    n_ie: usize,
    use_prisms: bool,
    curved: bool,
) -> Mesh {
    let [r0, r_in, r_out] = radii_spec;
    let nr = n_phys + n_ie;
    let mut radii = Vec::with_capacity(nr + 1);
    for k in 0..=n_phys {
        radii.push(r0 + (r_in - r0) * (k as f64) / (n_phys as f64));
    }
    for k in 1..=n_ie {
        radii.push(r_in + (r_out - r_in) * (k as f64) / (n_ie as f64));
    }

    // Deduplicate surface unit-sphere directions across the 3 cubed-sphere octant patches
    // using exact quantized integer cube coordinates in `0..=n_ang`.
    let mut surf_map: HashMap<[usize; 3], usize> = HashMap::new();
    let mut unit_dirs: Vec<[f64; 3]> = Vec::new();
    let mut patch_surf_vid = vec![vec![vec![0usize; n_ang + 1]; n_ang + 1]; 3];

    for (patch, p_grid) in patch_surf_vid.iter_mut().enumerate() {
        for (ia, row) in p_grid.iter_mut().enumerate() {
            for (ib, slot) in row.iter_mut().enumerate() {
                let cube_ijk = match patch {
                    0 => [n_ang, ia, ib],
                    1 => [ib, n_ang, ia],
                    _ => [ia, ib, n_ang],
                };
                let sid = *surf_map.entry(cube_ijk).or_insert_with(|| {
                    // Equiangular tan-projection reduces cubed-sphere cell distortion.
                    let map_c = |idx: usize| -> f64 {
                        if idx == n_ang {
                            1.0
                        } else {
                            (std::f64::consts::FRAC_PI_4 * (idx as f64) / (n_ang as f64)).tan()
                        }
                    };
                    let raw = [map_c(cube_ijk[0]), map_c(cube_ijk[1]), map_c(cube_ijk[2])];
                    let norm = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
                    let id = unit_dirs.len();
                    unit_dirs.push([raw[0] / norm, raw[1] / norm, raw[2] / norm]);
                    id
                });
                *slot = sid;
            }
        }
    }

    let n_surf = unit_dirs.len();
    let vid = |sid: usize, kr: usize| (kr * n_surf + sid) as u32;

    let mut cm = CellMesh::<3>::new();
    for &r in &radii {
        for &u in &unit_dirs {
            cm.points.push([r * u[0], r * u[1], r * u[2]]);
        }
    }

    // Volume cells (`Hex` or `Prism`)
    for kr in 0..nr {
        let dom = if kr < n_phys { 1 } else { 2 };
        for p_grid in &patch_surf_vid {
            for ia in 0..n_ang {
                for ib in 0..n_ang {
                    let s00 = p_grid[ia][ib];
                    let s10 = p_grid[ia + 1][ib];
                    let s11 = p_grid[ia + 1][ib + 1];
                    let s01 = p_grid[ia][ib + 1];
                    if use_prisms {
                        for tri in [[s00, s10, s11], [s00, s11, s01]] {
                            cm.push_cell(
                                CellType::Prism,
                                &[
                                    vid(tri[0], kr),
                                    vid(tri[1], kr),
                                    vid(tri[2], kr),
                                    vid(tri[0], kr + 1),
                                    vid(tri[1], kr + 1),
                                    vid(tri[2], kr + 1),
                                ],
                                dom,
                            );
                        }
                    } else {
                        cm.push_cell(
                            CellType::Hex,
                            &[
                                vid(s00, kr),
                                vid(s10, kr),
                                vid(s11, kr),
                                vid(s01, kr),
                                vid(s00, kr + 1),
                                vid(s10, kr + 1),
                                vid(s11, kr + 1),
                                vid(s01, kr + 1),
                            ],
                            dom,
                        );
                    }
                }
            }
        }
    }

    // Inner (`R = r0`, label 4) and outer (`R = r_out`, label 5) spherical facets
    for p_grid in &patch_surf_vid {
        for ia in 0..n_ang {
            for ib in 0..n_ang {
                let s00 = p_grid[ia][ib];
                let s10 = p_grid[ia + 1][ib];
                let s11 = p_grid[ia + 1][ib + 1];
                let s01 = p_grid[ia][ib + 1];
                if use_prisms {
                    for tri in [[s00, s10, s11], [s00, s11, s01]] {
                        cm.push_facet(&[vid(tri[0], 0), vid(tri[1], 0), vid(tri[2], 0)], 4);
                        cm.push_facet(&[vid(tri[0], nr), vid(tri[1], nr), vid(tri[2], nr)], 5);
                    }
                } else {
                    cm.push_facet(&[vid(s00, 0), vid(s10, 0), vid(s11, 0), vid(s01, 0)], 4);
                    cm.push_facet(&[vid(s00, nr), vid(s10, nr), vid(s11, nr), vid(s01, nr)], 5);
                }
            }
        }
    }

    // Symmetry-plane quad facets (`x = 0` -> 1, `y = 0` -> 2, `z = 0` -> 3)
    for kr in 0..nr {
        for idx in 0..n_ang {
            // x = 0: patch 1 (ib = 0, ia = idx..idx+1) and patch 2 (ia = 0, ib = idx..idx+1)
            for &(s0, s1) in &[
                (patch_surf_vid[1][idx][0], patch_surf_vid[1][idx + 1][0]),
                (patch_surf_vid[2][0][idx], patch_surf_vid[2][0][idx + 1]),
            ] {
                cm.push_facet(&[vid(s0, kr), vid(s1, kr), vid(s1, kr + 1), vid(s0, kr + 1)], 1);
            }
            // y = 0: patch 2 (ib = 0, ia = idx..idx+1) and patch 0 (ia = 0, ib = idx..idx+1)
            for &(s0, s1) in &[
                (patch_surf_vid[2][idx][0], patch_surf_vid[2][idx + 1][0]),
                (patch_surf_vid[0][0][idx], patch_surf_vid[0][0][idx + 1]),
            ] {
                cm.push_facet(&[vid(s0, kr), vid(s1, kr), vid(s1, kr + 1), vid(s0, kr + 1)], 2);
            }
            // z = 0: patch 0 (ib = 0, ia = idx..idx+1) and patch 1 (ia = 0, ib = idx..idx+1)
            for &(s0, s1) in &[
                (patch_surf_vid[0][idx][0], patch_surf_vid[0][idx + 1][0]),
                (patch_surf_vid[1][0][idx], patch_surf_vid[1][0][idx + 1]),
            ] {
                cm.push_facet(&[vid(s0, kr), vid(s1, kr), vid(s1, kr + 1), vid(s0, kr + 1)], 3);
            }
        }
    }

    if curved {
        // Project every spherical-shell edge midpoint (endpoints at the same radial index `kr`) onto the sphere `radii[kr]`.
        let disc = Mesh::D3(cm.clone()).discretise_mixed().unwrap();
        let mut curved_edges = Vec::new();
        for &(a, b) in &disc.edges {
            let kr_a = (a as usize) / n_surf;
            let kr_b = (b as usize) / n_surf;
            if kr_a == kr_b {
                let r = radii[kr_a];
                let ua = unit_dirs[(a as usize) % n_surf];
                let ub = unit_dirs[(b as usize) % n_surf];
                let m = [ua[0] + ub[0], ua[1] + ub[1], ua[2] + ub[2]];
                let mn = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
                curved_edges.push(([a, b], [r * m[0] / mn, r * m[1] / mn, r * m[2] / mn]));
            }
        }
        cm.set_edge_nodes(curved_edges);
    }

    Mesh::D3(cm)
}

#[test]
fn christoffel_backward_wave_detection_flagging_and_stabilization() {
    // 1. Verify pre-solve Christoffel backward-wave detector on both rotated cubic Si (`[22.5, 0, 0]`)
    //    and tilted R-plane sapphire (`[0, 57.6, 0]` / `[90, 57.6, -90]`).
    let si = Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 };
    let c6_si_2d = stiffness_global(&si, [90.0, 22.5, -90.0]);
    let chk_si = backward_wave_mpml_ratio(&c6_si_2d, 2329.0, [1.0, 0.0, 0.0]);
    assert!(
        chk_si.has_backward_wave && chk_si.min_projection < -0.05 && chk_si.recommended_mpml_ratio > 0.03,
        "rotated Si must exhibit Christoffel backward-wave slowness cusp: {chk_si:?}"
    );

    let sap = Stiffness::Trigonal {
        c11: 497.35e9,
        c12: 163.97e9,
        c13: 112.20e9,
        c14: -23.58e9,
        c33: 499.11e9,
        c44: 147.39e9,
    };
    let c6_sap = stiffness_global(&sap, [90.0, 57.6, -90.0]);
    let chk_sap = backward_wave_mpml_ratio(&c6_sap, 3980.0, [1.0, 0.0, 0.0]);
    assert!(
        chk_sap.has_backward_wave && chk_sap.recommended_mpml_ratio > 0.0,
        "tilted R-plane sapphire must exhibit backward-wave branch: {chk_sap:?}"
    );

    // 2. Solve an anisotropic rotated-crystal 2D native `Quad` strip (`x` fixed at 0, `+x` PML,
    //    oblique Bloch wavevector `k_z` via periodic top/bottom boundaries so bulk quasi-shear waves
    //    hit the `+x` PML right at the backward-wave slowness cusp!).
    //    Also test a strongly cusped cubic crystal (`c11 = 165.7e9, c12 = 135e9, c44 = 79.6e9`, Zener A = 5.18,
    //    rotated by 22.5 deg in `(x, z)`) so the amplifying mode appears cleanly among the lowest modes!
    let cusped_si = Stiffness::Cubic { c11: 165.7e9, c12: 135.0e9, c44: 79.6e9 };
    let c6_cusp = stiffness_global(&cusped_si, [90.0, 22.5, -90.0]);
    let chk_cusp = backward_wave_mpml_ratio(&c6_cusp, 2329.0, [1.0, 0.0, 0.0]);
    assert!(chk_cusp.has_backward_wave && chk_cusp.min_projection < -0.2, "{chk_cusp:?}");

    let lz = 1.0e-6;
    let mesh = quad_strip_with_pml(2.0e-6, 1.2e-6, lz, 12, 8, 4);
    assert!(!mesh.is_pure_simplex(), "expected native Quad mesh");

    let run_with_stab = |stab: PmlStabilization| {
        let mut p = Problem::new(Study::Eigenfrequency { n_modes: 10, f_search: 2.2e9 }, 2);
        p.plane = Plane::Strain;
        for dom in [1, 2] {
            p.materials.set(dom, "rho", 2329.0);
            p.materials.aniso.insert(
                dom,
                Anisotropy { stiffness: Some(cusped_si.clone()), euler_deg: [90.0, 22.5, -90.0], ..Default::default() },
            );
        }
        p.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } }];
        p.periodic = vec![Periodic { source: vec![2], dest: vec![3], phase: 2.2 }];
        p.materials.pml.insert(
            2,
            Pml {
                x: PmlDir::Plus,
                strength: 2.0,
                order: 2.0,
                re: 1.0,
                stabilization: stab,
                clamp_outer: true,
                ..Default::default()
            },
        );
        run(&mesh, &p, &Progress::new()).unwrap()
    };

    // Stabilization = Off: must produce amplifying `f.im < -1e-9 * f.re` backward-wave mode(s),
    // negative signed `Q < 0`, `[AMPLIFYING backward-wave mode]` label, `unstable_modes > 0`, and warning log!
    let r_off = run_with_stab(PmlStabilization::Off);
    let unstable_off = num_scalar(&r_off, "unstable_modes");
    assert!(
        unstable_off >= 1.0,
        "expected unstabilized backward-wave PML to report unstable_modes >= 1, got {unstable_off}; scalars={:?}",
        r_off.scalars
    );
    assert!(
        r_off.solutions.iter().any(|s| s.label.contains("[AMPLIFYING backward-wave mode]")),
        "expected [AMPLIFYING backward-wave mode] in solution labels: {:?}",
        r_off.solutions.iter().map(|s| &s.label).collect::<Vec<_>>()
    );
    assert!(
        r_off.log.iter().any(|l| l.contains("WARNING:") && l.contains("amplifying")),
        "expected amplifying mode warning in solver log: {:?}",
        r_off.log
    );
    // Verify negative signed Q (NOT masked by `.abs()`) and negative radiated power `P_pml < 0` on the amplifying mode!
    let mut saw_neg_q = false;
    for j in 1..=r_off.solutions.len() {
        let fim = num_scalar(&r_off, &format!("f im {j} (Hz)"));
        let fre = num_scalar(&r_off, &format!("f re {j} (Hz)"));
        if fim < -1e-9 * fre {
            let qj = num_scalar(&r_off, &format!("Q{j}"));
            let pj = num_scalar(&r_off, &format!("P_pml {j} (W)"));
            assert!(qj < 0.0, "signed Q must be negative when f.im < 0, got Q{j} = {qj}");
            assert!(pj <= 0.0, "P_pml must be non-positive when f.im < 0, got P_pml {j} = {pj}");
            saw_neg_q = true;
        }
    }
    assert!(saw_neg_q);

    // Stabilization = Auto, Mpml, and Adiabatic: ALL must eliminate `f.im < -1e-9 * f.re` (`unstable_modes == 0`)!
    for (name, stab) in [
        ("auto", PmlStabilization::Auto),
        ("mpml", PmlStabilization::Mpml { ratio: chk_cusp.recommended_mpml_ratio.max(0.35) }),
        ("adiabatic", PmlStabilization::Adiabatic { eta_max: 1.5 }),
    ] {
        let r_stab = run_with_stab(stab);
        let unstable = num_scalar(&r_stab, "unstable_modes");
        assert_eq!(
            unstable,
            0.0,
            "{name} stabilization must eliminate all amplifying modes (got unstable_modes={unstable}): {:?}",
            r_stab.solutions.iter().map(|s| &s.label).collect::<Vec<_>>()
        );
        assert!(
            !r_stab.solutions.iter().any(|s| s.label.contains("AMPLIFYING")),
            "{name} must have zero AMPLIFYING labels"
        );
        for j in 1..=r_stab.solutions.len() {
            let fim = num_scalar(&r_stab, &format!("f im {j} (Hz)"));
            let fre = num_scalar(&r_stab, &format!("f re {j} (Hz)"));
            assert!(fim >= -1e-9 * fre, "{name} mode {j}: fre={fre}, fim={fim}");
        }
        if name == "auto" {
            assert!(
                r_stab.log.iter().any(|l| l.contains("Christoffel check") && l.contains("p_min")),
                "Auto stabilization must log Christoffel diagnostic with p_min: {:?}",
                r_stab.log
            );
        }
    }
}

#[test]
fn outer_clamped_pml_eliminates_rigid_body_zero_modes_in_2d_and_3d() {
    // 3D floating elastic block surrounded on all 6 sides by structured `Prism` + `Hex` PML (`build_pml`),
    // with ZERO `BcKind::Fixed` conditions!
    let core_3d = femcore::meshing::geomdrive::block(1.0e-3, 0.8e-3, 0.6e-3, 0.45e-3, 3).unwrap().mesh;
    let phys_3d = simplex_to_adapt::<3>(&core_3d).unwrap();
    let pml_spec_3d = PmlSpec { outer_label: Some(10), ..PmlSpec::<3>::uniform(0.3e-3, 2) };
    let built_3d = build_pml(&phys_3d, &pml_spec_3d).unwrap();
    let mesh_3d = Mesh::D3(built_3d.mesh);
    assert!(!mesh_3d.is_pure_simplex(), "expected Tet + Prism + Hex mixed mesh");

    let make_3d_problem = |clamp_outer: bool| {
        let mut p = Problem::new(Study::Eigenfrequency { n_modes: 7, f_search: 0.0 }, 1);
        for (k, v) in [("E", 169e9), ("nu", 0.22), ("rho", 2329.0)] {
            p.materials.default.insert(k.into(), v);
        }
        for reg in &built_3d.regions {
            let dir = |c: i8| match c {
                1 => PmlDir::Plus,
                -1 => PmlDir::Minus,
                _ => PmlDir::Off,
            };
            let pml = Pml {
                x: dir(reg.code[0]),
                y: dir(reg.code[1]),
                z: dir(reg.code[2]),
                strength: 2.5,
                order: 2.0,
                re: 1.0,
                clamp_outer,
                ..Default::default()
            };
            for &(_, pml_label) in &reg.labels {
                p.materials.pml.insert(pml_label, pml.clone());
            }
        }
        p
    };

    let r_clamped = run(&mesh_3d, &make_3d_problem(true), &Progress::new()).unwrap();
    let f_min_clamped = r_clamped.solutions.iter().map(|s| s.value.unwrap()).fold(f64::INFINITY, f64::min);
    assert!(
        f_min_clamped > 1.0e5,
        "outer-clamped floating 3D PML must have ZERO rigid-body ~0 Hz modes, got lowest f = {f_min_clamped} Hz"
    );
    assert!(!r_clamped.log.iter().any(|l| l.contains("rigid-body")));

    let r_unclamped = run(&mesh_3d, &make_3d_problem(false), &Progress::new()).unwrap();
    let n_zero = r_unclamped.solutions.iter().filter(|s| s.value.unwrap().abs() < 100.0).count();
    assert!(
        n_zero >= 3,
        "unclamped floating 3D PML (`clamp_outer = false`) should retain ~0 Hz rigid-body modes, got {:?}",
        r_unclamped.solutions.iter().map(|s| s.value.unwrap()).collect::<Vec<_>>()
    );
}

#[test]
fn static_infinite_elements_match_analytical_3d_kelvin_1_over_r_stiffness_at_p1_and_p2() {
    // Exact 3D Kelvin / Stokes-Boussinesq rigid-sphere translation benchmark in an infinite isotropic
    // elastic space (`R in [R_0, infinity)`):
    // Applying uniform traction `t = [0, 0, t_z]` on `R = R_0` produces exact uniform displacement
    // `u(R_0) = [0, 0, U_0]` with analytical far-field stiffness (`u ~ 1/R` + `1/R^3`):
    //
    //   K_exact = F_z / U_0 = 24 * pi * G * (1 - nu) * R_0 / (5 - 6 * nu),   F_z = 4 * pi * R_0^2 * t_z.
    //
    // Tested on 5-layer native `Prism` AND 5-layer native `Hex` Rational Infinite Elements at BOTH
    // `order = 1` (`P_1` / `Q_1`) and `order = 2` (`Order2Serendipity`: 15-node `Prism`, 20-node `Hex`),
    // plus `Exponential` Infinite Elements with outer 1st-order spring Robin BC.
    let (e, nu) = (100.0e9_f64, 0.25_f64);
    let g = e / (2.0 * (1.0 + nu));
    let r0 = 1.0_f64;
    let r_in = 2.0_f64;
    let r_out = 3.5_f64;
    let tz = 1.0e6_f64;
    let f_total = 4.0 * std::f64::consts::PI * r0 * r0 * tz;
    let k_exact = 24.0 * std::f64::consts::PI * g * (1.0 - nu) * r0 / (5.0 - 6.0 * nu);
    let u0_exact = f_total / k_exact;

    for (cell_name, use_prisms) in [("Hex", false), ("Prism", true)] {
        for order in [1usize, 2usize] {
            let n_ang = if order == 1 { 12 } else { 7 };
            let n_phys = if order == 1 { 10 } else { 5 };
            let mesh = octant_spherical_shell_mesh(
                n_ang,
                [r0, r_in, r_out],
                n_phys,
                5, // 5-layer native Prism / Hex Infinite Element shell!
                use_prisms,
                order == 2,
            );
            assert!(!mesh.is_pure_simplex());

            let mut p = Problem::new(Study::Static, order);
            p.materials.default.insert("E".into(), e);
            p.materials.default.insert("nu".into(), nu);
            // Octant symmetry BCs for translation along z:
            // x = 0 (label 1): ux = 0; y = 0 (label 2): uy = 0; z = 0 (label 3): ux = 0, uy = 0
            // Inner sphere R = R_0 (label 4): uniform traction [0, 0, t_z]
            // Outer sphere R = R_out (label 5): NO explicit BC — clamped automatically by Rational IE!
            p.bcs = vec![
                Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
                Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
                Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, false], values: [0.0; 3] } },
                Bc { labels: vec![4], kind: BcKind::Load([0.0, 0.0, tz]) },
            ];
            p.infinite_elements.insert(
                2,
                InfiniteElement {
                    pml: Pml {
                        coord: PmlCoord::Spherical { center: [0.0; 3], r_in, thickness: r_out - r_in },
                        clamp_outer: true,
                        ..Default::default()
                    },
                    pole_dist: None, // auto-calibrates `d_pole = r_in` -> exact `1 - t = r_in / R_tilde`!
                    profile: IeProfile::Rational,
                    decay: 30.0,
                    clamp_outer: true,
                },
            );

            let res = run(&mesh, &p, &Progress::new()).unwrap();
            // Average `uz` across the inner cavity vertices `R = R_0` (first `n_surf` vertices of `cm`)
            let uz = res.solutions[0].flat("uz").unwrap();
            let Mesh::D3(cm) = &mesh else { unreachable!() };
            let mut u_sum = 0.0;
            let mut u_cnt = 0usize;
            for (vi, pt) in cm.points.iter().enumerate() {
                let rad = (pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]).sqrt();
                if (rad - r0).abs() < 1e-9 {
                    u_sum += uz[vi];
                    u_cnt += 1;
                }
            }
            let u0_num = u_sum / (u_cnt as f64);
            let k_num = f_total / u0_num;
            let rel_err = ((k_num - k_exact) / k_exact).abs();
            println!(
                "Rational IE 5-layer {cell_name} P{order} (dofs={}): K_num={k_num:.6e} N/m vs K_exact={k_exact:.6e} N/m (rel_err={:.4}%, u0_num={u0_num:.6e}, u0_exact={u0_exact:.6e})",
                res.dofs,
                100.0 * rel_err
            );
            assert!(
                rel_err < 0.005,
                "5-layer {cell_name} Rational IE at order {order} exceeded 0.5% stiffness error: rel_err = {:.4}% (u0_num={u0_num:.6e}, u0_exact={u0_exact:.6e})",
                100.0 * rel_err
            );
        }
    }

    // Also verify `IeProfile::Exponential` + outer 1st-order asymptotic spring Robin (`clamp_outer = false`)
    // via the full JSON `ProblemSpec` envelope (`Sms::solve`), AND verify dynamic studies reject `"infinite_element"`.
    let exp_mesh = octant_spherical_shell_mesh(7, [r0, r_in, r_out], 5, 5, false, true);
    struct OctantResolver;
    impl LabelResolver for OctantResolver {
        fn domains(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "core" => Ok(vec![1]),
                "inf" => Ok(vec![2]),
                _ => Err(Error::invalid(format!("unknown domain {name}"))),
            }
        }
        fn boundaries(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "sym_x" => Ok(vec![1]),
                "sym_y" => Ok(vec![2]),
                "sym_z" => Ok(vec![3]),
                "cavity" => Ok(vec![4]),
                _ => Err(Error::invalid(format!("unknown boundary {name}"))),
            }
        }
    }

    let mut spec = ProblemSpec::new("sms", 2, "static");
    for reg in ["core", "inf"] {
        spec.materials.insert(
            reg.into(),
            MaterialAssign::Props([("E".to_string(), e), ("nu".to_string(), nu), ("rho".to_string(), 2500.0)].into()),
        );
    }
    spec.bcs = vec![
        BcSpec { on: "sym_x".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
        BcSpec { on: "sym_y".into(), kind: "fixed".into(), params: json!({ "uy": 0.0 }) },
        BcSpec { on: "sym_z".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uy": 0.0 }) },
        BcSpec { on: "cavity".into(), kind: "load".into(), params: json!({ "traction": [0.0, 0.0, tz] }) },
        BcSpec {
            on: "inf".into(),
            kind: "infinite_element".into(),
            params: json!({
                "coord": "spherical",
                "profile": "exponential",
                "decay": 80.0
            }),
        },
    ];
    let r_exp = Sms.solve(&exp_mesh, &spec, &mut Ctx::new(&OctantResolver)).unwrap();
    let uz_exp = r_exp.solutions[0].flat("uz").unwrap();
    let Mesh::D3(cm_exp) = &exp_mesh else { unreachable!() };
    let u0_exp: f64 = cm_exp
        .points
        .iter()
        .enumerate()
        .filter(|(_, pt)| ((pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]).sqrt() - r0).abs() < 1e-9)
        .map(|(vi, _)| uz_exp[vi])
        .sum::<f64>()
        / (cm_exp
            .points
            .iter()
            .filter(|pt| ((pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]).sqrt() - r0).abs() < 1e-9)
            .count() as f64);
    let k_exp = f_total / u0_exp;
    assert!(
        ((k_exp - k_exact) / k_exact).abs() < 0.015,
        "Exponential IE + outer Robin stiffness error too large: K_exp={k_exp:.6e} vs K_exact={k_exact:.6e}"
    );

    // Dynamic studies must reject `"infinite_element"` cleanly with guidance to use `"pml"`.
    let mut bad_dyn = spec.clone();
    bad_dyn.study.kind = "eigenfrequency".into();
    bad_dyn.study.params = json!({ "n_modes": 2, "f_search": 1.0e5 });
    bad_dyn.bcs.retain(|b| b.kind != "load");
    let err = Sms.solve(&exp_mesh, &bad_dyn, &mut Ctx::new(&OctantResolver)).unwrap_err().to_string();
    assert!(err.contains("pml"), "expected error instructing user to use 'pml', got: {err}");
}

#[test]
fn mixed_cell_lamb_wave_pml_matches_or_beats_split_simplex_with_fewer_dofs_and_axisym_spherical_r0_works() {
    // 1. 2D Lamb-wave resonator strip (`Tri` physical core + structured `Quad` lateral PMLs via `build_pml::<2>`):
    //    Compare native `Tri + Quad` `Order2Serendipity` solve against split-triangle `SimplexMesh` solve
    //    and against a wider high-resolution reference PML.
    let plate_tri = femcore::meshing::geomdrive::plate(2.0e-3, 0.4e-3, 0.12e-3).unwrap().mesh;
    let phys_2d = simplex_to_adapt::<2>(&plate_tri).unwrap();
    let mut pspec = PmlSpec::<2>::uniform(0.0, 0);
    pspec.outer_label = Some(10);
    pspec.thickness[0] = [0.6e-3, 0.6e-3];
    pspec.layers[0] = [4, 4];
    let built_2d = build_pml(&phys_2d, &pspec).unwrap();
    let mixed_mesh = Mesh::D2(built_2d.mesh.clone());
    assert!(!mixed_mesh.is_pure_simplex());
    let (split_simplex, rep) = SimplexMesh::from_cells_split(&built_2d.mesh).unwrap();
    assert!(rep.quads > 0);

    let mut prob = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 6.5e6 }, 2);
    prob.plane = Plane::Strain;
    for (k, v) in [("E", 169e9), ("nu", 0.22), ("rho", 2329.0)] {
        prob.materials.default.insert(k.into(), v);
    }
    for reg in &built_2d.regions {
        let dir = |c: i8| match c {
            1 => PmlDir::Plus,
            -1 => PmlDir::Minus,
            _ => PmlDir::Off,
        };
        let pml = Pml {
            x: dir(reg.code[0]),
            z: dir(reg.code[1]),
            strength: 3.5,
            order: 2.0,
            re: 1.0,
            clamp_outer: true,
            ..Default::default()
        };
        for &(_, lab) in &reg.labels {
            prob.materials.pml.insert(lab, pml.clone());
        }
    }

    let r_mixed = run(&mixed_mesh, &prob, &Progress::new()).unwrap();
    let r_split = solve_solid(&split_simplex, &prob).unwrap();
    assert!(
        r_mixed.dofs < r_split.dofs,
        "native Tri+Quad Order2Serendipity must use fewer DOFs than split-Tri P2: {} vs {}",
        r_mixed.dofs,
        r_split.dofs
    );

    // Find the lowest-loss Lamb mode confined to the core (`layer fraction < 0.35`) in both solves:
    let pick_core_mode = |r: &nyx_physics::SolveResult| -> (f64, f64) {
        (1..=r.solutions.len())
            .filter_map(|j| {
                let lf = num_scalar(r, &format!("layer fraction {j}"));
                let fre = num_scalar(r, &format!("f re {j} (Hz)"));
                let fim = num_scalar(r, &format!("f im {j} (Hz)"));
                (lf < 0.45 && fim > 0.0).then_some((lf, fre, fim))
            })
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .map(|(_, fre, fim)| (fre, fre / (2.0 * fim)))
            .expect("expected core-confined Lamb mode")
    };
    let (f_mix, q_mix) = pick_core_mode(&r_mixed);
    let (f_spl, q_spl) = pick_core_mode(&r_split);
    assert!(((f_mix - f_spl) / f_spl).abs() < 2e-3, "native mixed Lamb frequency {f_mix} Hz vs split {f_spl} Hz");
    assert!(((q_mix - q_spl) / q_spl).abs() < 0.05, "native mixed Lamb Q {q_mix} vs split {q_spl}");

    // 2. 2D Axisymmetric Spherical PML touching `r = 0` (`Quad` + `Tri` half-disk in `(r >= 0, z)`):
    //    Inner elastic sphere `R in [0, 1.0e-3]` (domain 1) surrounded by Spherical PML `R in [1.0e-3, 1.5e-3]`
    //    (domain 2) touching `r = 0` at both poles `theta = 0` and `theta = pi`!
    let (n_th, nr_core, nr_pml) = (16usize, 6usize, 4usize);
    let nr_tot = nr_core + nr_pml;
    let (r_core, r_pml_out) = (1.0e-3_f64, 1.5e-3_f64);
    let mut cm_ax = CellMesh::<2>::new();
    cm_ax.points.push([0.0, 0.0]); // center vertex 0
    let ring_vid = |kr: usize, it: usize| -> u32 { (1 + (kr - 1) * (n_th + 1) + it) as u32 };
    for kr in 1..=nr_tot {
        let rad = if kr <= nr_core {
            r_core * (kr as f64) / (nr_core as f64)
        } else {
            r_core + (r_pml_out - r_core) * ((kr - nr_core) as f64) / (nr_pml as f64)
        };
        for it in 0..=n_th {
            let th = std::f64::consts::PI * (it as f64) / (n_th as f64);
            let r = (rad * th.sin()).max(0.0);
            let z = rad * th.cos();
            cm_ax.points.push([r, z]);
        }
    }
    for it in 0..n_th {
        cm_ax.push_cell(CellType::Tri, &[0, ring_vid(1, it + 1), ring_vid(1, it)], 1);
    }
    for kr in 1..nr_tot {
        let dom = if kr < nr_core { 1 } else { 2 };
        for it in 0..n_th {
            cm_ax.push_cell(
                CellType::Quad,
                &[ring_vid(kr, it), ring_vid(kr + 1, it), ring_vid(kr + 1, it + 1), ring_vid(kr, it + 1)],
                dom,
            );
        }
    }
    // Boundary facets along outer sphere `R = r_pml_out` (label 2) and symmetry axis `r = 0` (label 1)
    for it in 0..n_th {
        cm_ax.push_facet(&[ring_vid(nr_tot, it), ring_vid(nr_tot, it + 1)], 2);
    }
    cm_ax.push_facet(&[0, ring_vid(1, 0)], 1);
    cm_ax.push_facet(&[0, ring_vid(1, n_th)], 1);
    for kr in 1..nr_tot {
        cm_ax.push_facet(&[ring_vid(kr, 0), ring_vid(kr + 1, 0)], 1);
        cm_ax.push_facet(&[ring_vid(kr, n_th), ring_vid(kr + 1, n_th)], 1);
    }
    let ax_mesh = Mesh::D2(cm_ax);
    assert!(!ax_mesh.is_pure_simplex());

    let mut ax_prob =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 4, f_search: 2.5e6, azimuthal_order: 0 }, 2);
    // Slightly softer surrounding medium (`domain 2`) so breathing sphere modes radiate cleanly with finite Q!
    for (k, v) in [("E", 180e9), ("nu", 0.25), ("rho", 3200.0)] {
        ax_prob.materials.set(1, k, v);
    }
    for (k, v) in [("E", 40e9), ("nu", 0.25), ("rho", 2200.0)] {
        ax_prob.materials.set(2, k, v);
    }
    ax_prob.materials.pml.insert(
        2,
        Pml {
            coord: PmlCoord::Spherical { center: [0.0; 3], r_in: r_core, thickness: r_pml_out - r_core },
            strength: 3.5,
            order: 2.0,
            clamp_outer: true,
            ..Default::default()
        },
    );

    let r_ax = run(&ax_mesh, &ax_prob, &Progress::new()).unwrap();
    assert!(!r_ax.solutions.is_empty());
    assert_eq!(num_scalar(&r_ax, "unstable_modes"), 0.0);
    // Verify `u_r == 0` and `u_theta == 0` regularity on the `r = 0` symmetry axis (`m = 0`), and `P_pml > 0`!
    let simplex_ax = nyx_physics::discretise(&ax_mesh).unwrap();
    let coords_ax = simplex_ax.dof_coordinates(2);
    let re_ur = r_ax.solutions[0].flat("Re ux").unwrap();
    let abs_ut = r_ax.solutions[0].flat("|uy|").unwrap();
    for sd in 0..simplex_ax.n_dofs(2) {
        if coords_ax[2 * sd].abs() < 1e-12 {
            assert!(re_ur[sd].abs() < 1e-12, "u_r must vanish on r = 0, got {}", re_ur[sd]);
            assert!(abs_ut[sd].abs() < 1e-12, "|u_theta| must vanish on r = 0 at m = 0, got {}", abs_ut[sd]);
        }
    }
    let has_radiating = (1..=r_ax.solutions.len()).any(|j| {
        num_scalar(&r_ax, &format!("f im {j} (Hz)")) > 0.0 && num_scalar(&r_ax, &format!("P_pml {j} (W)")) > 0.0
    });
    assert!(has_radiating, "expected positive Im f and P_pml on axisymmetric spherical PML: {:?}", r_ax.scalars);
}
