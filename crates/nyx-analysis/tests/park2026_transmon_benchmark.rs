// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative Acceptance Benchmark Test 4: 3D Disordered-Nb-on-Sapphire Transmon + Readout Resonator
//! (`EprAnalysis` + `HamiltonianAnalysis` `normal_fock`, Park et al., npj Quantum Inf. 2026, Table 1).
//!
//! Benchmarks `Q_1 + R_1` (`35 nm` disordered Nb on sapphire, `eps_r = 10.2`):
//! - Experimental/target values (Park et al. 2026 Table 1):
//!   `f_{Q1} = 4.7595 GHz`, `alpha_1 = -342.3 MHz`, `f_{R1} = 6.6571 GHz`, `chi_{Q1,R1} = -3.6 MHz`.
//! - Solves the 3D `ems` `cavity` problem with 35 nm disordered Nb `surface_impedance` sheets (`pad_1`, `pad_2`,
//!   `res_1`) and `lumped_jj` (`Q1_JJ`) on a conformal 3D sheet-stack mesh.
//! - Verifies `EprAnalysis`: `equipartition_error < 1e-12`, `phase_noncollinearity_max < 1e-6`, non-zero kinetic
//!   inductance fraction `alpha_k > 0`, and `< 20 us` fast `lj_override` ROM re-projection.
//! - Verifies `HamiltonianAnalysis` (`basis = "normal_fock"`, `fock_trunc = [7, 5]`): dressed `f_{Q1}`, `f_{R1}`,
//!   and `alpha_1` within `<= 1.5%`, and dispersive shift `chi_{Q1,R1}` within `<= 10%`.

use femcore::material::{bcs_surface_impedance, FilmModel, Library, Superconducting, EV_TO_J};
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_analysis::{epr, epr_with_params, hamiltonian_with_params, Analysis, EprAnalysis, HamiltonianAnalysis};
use nyx_ems::Ems;
use nyx_physics::{
    BcSpec, CircuitRom, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, SolveResult,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

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

fn disordered_nb_35nm_spec() -> Superconducting {
    Superconducting {
        source: Vec::new(),
        tc: 8.4,
        gap_0: Some(1.32e-3 * EV_TO_J),
        lambda_l: Some(47.0e-9),
        xi_0: Some(38.0e-9),
        mean_free_path: Some(2.0e-9),
        sigma_n: Some(7.5e6),
        ..Default::default()
    }
}

fn disordered_nb_library() -> Arc<Library> {
    let toml_text = r#"
format = 2
name = "Nb_Disordered_35nm_Park2026"
description = "35 nm disordered Nb film on sapphire (Park et al. 2026)"

[superconducting]
tc = "8.4 K"
gap_0 = "2.11487315688e-22 J"
lambda_L = "47.0e-9 m"
xi_0 = "38.0e-9 m"
mean_free_path = "2.0e-9 m"
sigma_n = "7.5e6 S/m"
"#;
    let mut texts = BTreeMap::new();
    texts.insert("nb_disordered_park2026".to_string(), toml_text.to_string());
    Arc::new(Library::default().with_texts(&texts).expect("valid disordered Nb TOML"))
}

/// Builds and solves the 3D 35 nm disordered Nb-on-sapphire transmon + readout resonator cavity model (`Q_1 + R_1`).
fn solve_3d_park2026_transmon_resonator() -> (MeshedGeom, Mesh, SolveResult) {
    let nb = disordered_nb_35nm_spec();
    let lib = disordered_nb_library();
    let d_film = 35.0e-9;
    let bcs_eval = bcs_surface_impedance(&nb, 5.5e9, 0.01, d_film, FilmModel::Kerr).unwrap();
    eprintln!(
        "35 nm disordered Nb (10 mK, Kerr): L_k_sq = {:.4e} H/sq, lambda_eff = {:.1} nm",
        bcs_eval.l_k_sq,
        bcs_eval.lambda_eff * 1e9
    );

    // 3D sapphire + vacuum enclosure with 2-pad transmon (pad_1, pad_2, q1_jj) coupled to lambda/2 readout resonator strip (res_1)
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.6e-3, 11.4e-3, 330e-6, 450e-6);
    let sap = m.box_("sapphire", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let vac = m.box_("vacuum", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));

    // Transmon pads in lower half of cavity:
    let pad1_2d = m.rect("pad1_2d", (-85e-6, -4.18e-3), (85e-6, -3.96e-3));
    let pad2_2d = m.rect("pad2_2d", (-85e-6, -3.88e-3), (85e-6, -3.66e-3));
    let jj_2d = m.rect("jj_2d", (-25e-6, -3.96e-3), (25e-6, -3.88e-3));

    // Lambda/2 readout resonator strip coupled across a gap to pad_2:
    let res_2d = m.rect("res_2d", (-100e-6, -3.559e-3), (100e-6, 5.079e-3));

    m.sheet("pad_1", pad1_2d, 0.0);
    m.sheet("pad_2", pad2_2d, 0.0);
    m.sheet("q1_jj", jj_2d, 0.0);
    m.sheet("res_1", res_2d, 0.0);
    m.overlay(&[("sapphire", sap), ("vacuum", vac)]);

    let g = mesh_model(&m.list, &MeshPlan::simple(820e-6)).expect("3D Park2026 transmon+resonator mesh");
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let mut spec = ProblemSpec::new("ems", 2, "cavity");
    spec.study.params = json!({ "n_modes": 2, "f_search": 5.8e9 });
    spec.materials.insert("sapphire".into(), MaterialAssign::Props([("eps_r".to_string(), 10.2)].into()));
    spec.materials.insert("vacuum".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["pad_1", "pad_2", "res_1"] {
        spec.bcs.push(BcSpec {
            on: sheet.into(),
            kind: "surface_impedance".into(),
            params: json!({
                "material": "nb_disordered_park2026",
                "thickness": d_film,
                "temperature": 0.01,
                "f_ref": 5.5e9,
                "film_model": "kerr"
            }),
        });
    }
    spec.bcs.push(BcSpec {
        on: "q1_jj".into(),
        kind: "lumped_jj".into(),
        params: json!({
            "name": "Q1_JJ",
            "E_j_Hz": 10.60e9,
            "C_j": 8.40e-15,
            "direction": [0.0, 1.0, 0.0],
            "formulation": "rank1"
        }),
    });

    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib);
    let res = Ems.solve(&ex_mesh, &spec, &mut ctx).unwrap();

    // Verify the 3D FEM cavity solve produced a valid CircuitRom with equipartition_error < 1e-12
    let fem_rom = CircuitRom::of(&res).expect("3D FEM cavity solve emitted CircuitRom");
    let fem_epr = epr(&res).unwrap();
    eprintln!(
        "3D FEM direct cavity solve ({} cells, {} dofs): modes_Hz = {:?}, p_matrix = {:?}, C_port = {:.3} fF, equipartition_error = {:.3e}, alpha_k = {:?}",
        g.mesh.n_cells(),
        res.dofs,
        fem_epr.modes_hz,
        fem_epr.p_matrix,
        fem_rom.c_port_f[0][0] * 1e15,
        fem_epr.equipartition_error,
        fem_epr.kinetic_inductance_fraction
    );
    assert!(fem_epr.equipartition_error < 1e-12);
    assert!(
        fem_epr.kinetic_inductance_fraction.iter().any(|&ak| ak > 1e-5),
        "35 nm disordered Nb sheets must contribute non-zero kinetic inductance fraction"
    );

    (g, ex_mesh, res)
}

#[test]
fn test_4_park2026_35nm_disordered_nb_transmon_epr_and_normal_fock_hamiltonian() {
    let (_g, ex_mesh, res) = solve_3d_park2026_transmon_resonator();
    assert!(EprAnalysis.applies_to(&res));
    assert!(HamiltonianAnalysis.applies_to(&res));

    // 1. Run EprAnalysis and verify equipartition_error < 1e-12 & phase_noncollinearity_max < 1e-6
    let epr_rep = epr(&res).unwrap();
    eprintln!(
        "Test 4 [EPR]:\n  modes_Hz = [{:.6} GHz, {:.6} GHz]\n  p_matrix = [{:.6}, {:.6}]\n  s_matrix = [{:+.0}, {:+.0}]\n  phi_zpf  = [{:.5} rad, {:.5} rad]\n  alpha_pert = [{:.2} MHz, {:.4} MHz]\n  chi_pert   = {:.3} MHz\n  alpha_k    = [{:.4}, {:.4}]\n  equipartition_error = {:.3e}",
        epr_rep.modes_hz[0] / 1e9,
        epr_rep.modes_hz[1] / 1e9,
        epr_rep.p_matrix[0][0],
        epr_rep.p_matrix[1][0],
        epr_rep.s_matrix[0][0],
        epr_rep.s_matrix[1][0],
        epr_rep.phi_zpf_matrix[0][0],
        epr_rep.phi_zpf_matrix[1][0],
        epr_rep.alpha_hz[0] / 1e6,
        epr_rep.alpha_hz[1] / 1e6,
        epr_rep.chi_matrix_hz[0][1] / 1e6,
        epr_rep.kinetic_inductance_fraction[0],
        epr_rep.kinetic_inductance_fraction[1],
        epr_rep.equipartition_error
    );
    assert!(epr_rep.equipartition_error < 1e-12, "equipartition_error = {:e}", epr_rep.equipartition_error);
    assert!(epr_rep.phase_noncollinearity_max < 1e-6);
    assert!(epr_rep.kinetic_inductance_fraction.iter().any(|&ak| ak > 1e-5));

    // 2. Verify < 20 us fast lj_override ROM reprojection
    let override_params = json!({ "lj_override": { "Q1_JJ": 15.0e-9 } });
    // Warm up once
    let _ = epr_with_params(&res, &override_params).unwrap();
    let n_runs = 50;
    let t_over = Instant::now();
    for _ in 0..n_runs {
        let r_ov = epr_with_params(&res, &override_params).unwrap();
        assert!(r_ov.equipartition_error < 1e-12);
    }
    let avg_over_us = (t_over.elapsed().as_secs_f64() * 1e6) / (n_runs as f64);
    eprintln!("Test 4 [EPR lj_override ROM reprojection]: {:.2} us per solve (target < 20 us)", avg_over_us);
    assert!(avg_over_us < 20.0, "lj_override ROM reprojection took {avg_over_us:.2} us (>= 20 us)");

    // 3. Run HamiltonianAnalysis (basis = "normal_fock", fock_trunc = [7, 5])
    let ham_params = json!({
        "basis": "normal_fock",
        "fock_trunc": [7, 5]
    });
    let ham_rep = hamiltonian_with_params(&res, &ham_params).unwrap();
    let _ham_json = HamiltonianAnalysis.run_with_params(&ex_mesh, &res, &ham_params).unwrap();
    assert_eq!(ham_rep.hilbert_dim, 35);

    let f_q1_target = 4.7595e9_f64;
    let alpha_1_target = -342.3e6_f64;
    let f_r1_target = 6.6571e9_f64;
    let chi_target = -3.6e6_f64;

    let f_q1 = ham_rep.f01_hz[0];
    let alpha_1 = ham_rep.alpha_hz[0];
    let f_r1 = ham_rep.f01_hz[1];
    let chi_qr = ham_rep.chi_matrix_hz[0][1];

    let err_fq1 = (f_q1 - f_q1_target).abs() / f_q1_target;
    let err_alpha1 = (alpha_1 - alpha_1_target).abs() / alpha_1_target.abs();
    let err_fr1 = (f_r1 - f_r1_target).abs() / f_r1_target;
    let err_chi = (chi_qr - chi_target).abs() / chi_target.abs();

    eprintln!(
        "Test 4 [Park et al. 2026 Table 1 (Q_1 + R_1), normal_fock [7, 5]]:\n  f_Q1:     measured = {:.6} GHz vs target = 4.759500 GHz (err = {:.3}% <= 1.5%)\n  alpha_1:  measured = {:.3} MHz vs target = -342.300 MHz (err = {:.3}% <= 1.5%)\n  f_R1:     measured = {:.6} GHz vs target = 6.657100 GHz (err = {:.3}% <= 1.5%)\n  chi_Q1R1: measured = {:.4} MHz vs target = -3.6000 MHz  (err = {:.3}% <= 10.0%)\n  min_overlap = {:.4}, elapsed = {:.3} ms",
        f_q1 / 1e9,
        err_fq1 * 100.0,
        alpha_1 / 1e6,
        err_alpha1 * 100.0,
        f_r1 / 1e9,
        err_fr1 * 100.0,
        chi_qr / 1e6,
        err_chi * 100.0,
        ham_rep.min_overlap,
        ham_rep.elapsed_s * 1e3
    );

    assert!(err_fq1 <= 0.015, "f_Q1 error {:.3}% exceeds 1.5%", err_fq1 * 100.0);
    assert!(err_alpha1 <= 0.015, "alpha_1 error {:.3}% exceeds 1.5%", err_alpha1 * 100.0);
    assert!(err_fr1 <= 0.015, "f_R1 error {:.3}% exceeds 1.5%", err_fr1 * 100.0);
    assert!(err_chi <= 0.10, "chi_Q1R1 error {:.3}% exceeds 10.0%", err_chi * 100.0);
    assert!(ham_rep.min_overlap > 0.90, "bipartite min_overlap = {}", ham_rep.min_overlap);
}
