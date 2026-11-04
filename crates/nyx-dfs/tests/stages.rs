//! The kernels are three stages (ADR 9): assemble (matrices, constraints), solve (the solution), fields (plot fields and the same quantities
//! per cell). The one-call `solve_*` is their composition; these tests drive the stages on their own.
mod common;
use common::*;
use nyx_dfs::{Bc, BcKind, Problem, Study};

#[test]
fn electrostatics_stages_compose_to_the_one_call_solver() {
    let mesh = rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![4], kind: BcKind::Ground }];
    let whole = nyx_dfs::solve_electrostatics(&mesh, &spec).unwrap();
    let asm = nyx_dfs::assemble_electrostatics(&mesh, &spec).unwrap();
    assert_eq!(asm.k.nrows(), mesh.n_dofs(spec.order));
    let (v, _) = nyx_dfs::solve_potential(&mesh, &spec, &asm).unwrap();
    let (fields, cells) = nyx_dfs::electrostatics_fields(&mesh, spec.order, &v);
    let want = &whole.solutions[0];
    assert_eq!(want.field("V").unwrap().values, v);
    assert_eq!(fields.len(), want.fields.len());
    // the flux is discontinuous across the interface between the two dielectrics only in the normal direction: the per-cell
    // field keeps the jump that the node average hides
    let ex = cells.get("Ex").unwrap();
    let avg = cells.node_average(&mesh, "Ex").unwrap();
    assert_eq!(avg, want.flat("Ex").unwrap());
    assert_eq!(ex.re.len(), mesh.n_cells() * mesh.nloc(spec.order));
    let spread =
        ex.re.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - ex.re.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(spread > 0.0);
    // no component named after a mesh axis: the cross-section has x and z
    assert!(cells.get("Ey").is_none() && cells.get("Ez").is_some());
}

#[test]
fn the_per_cell_flux_jumps_where_the_node_average_does_not() {
    // two dielectrics side by side across x: D_x is continuous, so E_x jumps by the ratio of the permittivities at the interface
    let mesh = rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![4], kind: BcKind::Ground }];
    let r = nyx_dfs::solve_electrostatics(&mesh, &spec).unwrap();
    let ex = r.solutions[0].cell_fields.iter().find(|f| f.name == "Ex").unwrap();
    let nloc = ex.nloc;
    // the largest and the smallest cell mean of |Ex| differ by about the permittivity ratio (4)
    let mean: Vec<f64> =
        (0..mesh.n_cells()).map(|c| ex.re[c * nloc..(c + 1) * nloc].iter().sum::<f64>() / nloc as f64).collect();
    let hi = mean.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let lo = mean.iter().filter(|x| x.abs() > 1e-9 * hi).fold(f64::INFINITY, |m, x| m.min(x.abs()));
    assert!((hi / lo - 4.0).abs() < 0.2, "ratio {}", hi / lo);
}

#[test]
fn heat_stages_give_the_temperature_and_the_flux_per_cell() {
    let mesh = rect(1.0, 0.4, 10, 4);
    let mut spec = Problem::new(Study::HeatStationary, 2);
    spec.materials.default.insert("k".into(), 3.0);
    spec.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![4], kind: BcKind::Temperature(280.0) },
    ];
    let sys = nyx_dfs::assemble_heat(&mesh, &spec).unwrap();
    let st = nyx_dfs::solve_heat_states(&mesh, &spec, &sys).unwrap();
    assert!(st.times.is_empty() && st.temps.len() == 1);
    let (fields, cells) = nyx_dfs::heat_fields(&mesh, &st.temps[0], &sys.ktensor, spec.order);
    assert!(fields.iter().any(|f| f.name == "q" && f.components == ["x", "z"]));
    // a linear temperature drop of 20 K over 1 m in k = 3: q_x = 60 W/m^2 in every cell
    let q = cells.get("qx").unwrap();
    assert!(q.re.iter().all(|v| (v - 60.0).abs() < 1e-6), "{:?}", &q.re[..4]);
}
