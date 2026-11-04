//! The kernels are three stages: assemble (matrices, constraints), solve (eigenpairs or solution), fields (plot fields and the
//! same quantities per cell). The one-call `solve_*` is their composition; these tests drive the stages on their own.
mod common;
use femcore::material::*;
use femcore::problem::*;

#[test]
fn optical_mode_stages_compose_to_the_one_call_solver() {
    use femcore::ems::*;
    let mesh = common::rect(1.0e-6, 0.5e-6, 12, 6);
    let nc = mesh.n_cells();
    let eps = vec![EpsCell { xx: 4.0, xy: 0.0, yy: 4.0, zz: 4.0 }; nc];
    let mu = vec![1.0; nc];
    let whole = solve_modes_tensor(&mesh, 2.0e14, &eps, &mu, 1.5, false, 3, 2).unwrap();
    let asm = assemble_modes(&mesh, 2.0e14, &eps, &mu, 2).unwrap();
    let staged = solve_mode_pairs(&asm, 1.5, false, 3).unwrap();
    assert_eq!(whole.neff, staged.neff);
    let (fields, cells) = mode_fields_cells(&mesh, &staged, 0, 2.0e14, 2);
    assert_eq!(fields.len(), 6);
    // the components carry the wafer names: x lateral, z normal, y along the guide
    for name in ["Ex", "Ez", "Ey"] {
        let avg = cells.node_average(&mesh, name).unwrap();
        let plot = fields.iter().find(|f| f.name == name).unwrap();
        let scale = plot.values.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        assert!(avg.iter().zip(&plot.values).all(|(a, b)| (a - b).abs() < 1e-12 * scale), "{name}");
    }
}

