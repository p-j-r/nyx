//! The axisymmetric eigenfrequency study (`sms`, ADR 66): the envelope, the typed problem and the library entry give the same modes; the
//! study refuses what it cannot say; with an absorbing layer the modes are complex and carry Q, an eigenvalue and a layer fraction.
//! The meshes are rings off the axis (the axis conditions are the kernel's K1, tested elsewhere).
use femcore::mesh::SimplexMesh;
use nyx_physics::*;
use nyx_sms::{Bc, BcKind, Problem, Sms, Study};
use serde_json::json;

/// A ring section `w` wide and `h` high, its inner edge at radius `r0`, as the (r, z) half plane: domain 1 everywhere; `left` is boundary 1.
fn ring(r0: f64, w: f64, h: f64, nx: usize, ny: usize) -> SimplexMesh {
    let flat = femcore::testkit::rect(w, h, nx, ny);
    let mut pts = flat.points.clone();
    for i in 0..flat.n_points() {
        pts[2 * i] += r0;
    }
    SimplexMesh::new(
        2,
        pts,
        flat.cells.clone(),
        flat.cell_labels.clone(),
        flat.facets.clone(),
        flat.facet_labels.clone(),
    )
}

fn exchange(mesh: &SimplexMesh) -> Mesh {
    Mesh::from_solver_mesh(mesh).unwrap()
}

struct Table;
impl LabelResolver for Table {
    fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
        match r {
            "bar" => Ok(vec![1]),
            "layer" => Ok(vec![2]),
            _ => Err(Error::invalid(format!("no region {r}"))),
        }
    }
    fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
        match s {
            "left" => Ok(vec![1]),
            "bottom" => Ok(vec![2]),
            "top" => Ok(vec![3]),
            "right" => Ok(vec![4]),
            _ => Err(Error::invalid(format!("no selector {s}"))),
        }
    }
}

fn spec(m: i64, f_search: f64) -> ProblemSpec {
    let mut s = ProblemSpec::new("sms", 2, "axisymmetric_eigenfrequency");
    s.study.params = json!({ "n_modes": 4, "f_search": f_search, "azimuthal_order": m });
    s.materials.insert(
        "bar".into(),
        MaterialAssign::Props([("E".to_string(), 70e9), ("nu".to_string(), 0.33), ("rho".to_string(), 2700.0)].into()),
    );
    s.bcs.push(BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }) });
    s
}

fn typed(m: i64, f_search: f64) -> Problem {
    let mut p = Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 4, f_search, azimuthal_order: m }, 2);
    for (k, v) in [("E", 70e9), ("nu", 0.33), ("rho", 2700.0)] {
        p.materials.set(1, k, v);
    }
    p.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    p
}

fn frequencies(r: &SolveResult) -> Vec<f64> {
    r.solutions.iter().map(|s| s.value.unwrap()).collect()
}

fn scalar(r: &SolveResult, name: &str) -> Option<f64> {
    r.scalars.iter().find(|(n, _)| n == name).and_then(|(_, v)| match v {
        Scalar::Num(x) => Some(*x),
        _ => None,
    })
}

#[test]
fn the_envelope_the_typed_problem_and_the_library_entry_give_the_same_modes() {
    let mesh = ring(5.0e-3, 1.0e-3, 0.4e-3, 8, 3);
    let ex = exchange(&mesh);
    let table = Table;
    for m in [0i64, 2] {
        let env = Sms.solve(&ex, &spec(m, 0.0), &mut Ctx::new(&table)).unwrap();
        assert_eq!(
            (env.physics.as_str(), env.study.as_str(), env.dim, env.order),
            ("sms", "axisymmetric_eigenfrequency", 2, 2)
        );
        assert_eq!(env.dofs, 3 * mesh.n_dofs(2));
        let typ = nyx_sms::run(&ex, &typed(m, 0.0), &femcore::progress::Progress::new()).unwrap();
        // the library entry reads the plain eigenfrequency study and takes m as an argument
        let mut plain = typed(m, 0.0);
        plain.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
        let (lib, _) = nyx_sms::solve_solid_axisym(&mesh, &plain, m as f64).unwrap();
        assert_eq!(frequencies(&env), frequencies(&typ), "m = {m}: envelope against typed problem");
        assert_eq!(frequencies(&typ), frequencies(&lib), "m = {m}: typed problem against library entry");
        assert!(frequencies(&env).iter().all(|f| *f > 0.0));
        assert!(env.log.iter().any(|l| l.contains(&format!("azimuthal order m = {m}"))), "{:?}", env.log);
        // every mode has an eigenvalue [Re f, Im f] and a layer fraction of 0 (no layer)
        for (j, s) in env.solutions.iter().enumerate() {
            let [re, im] = s.eigenvalue.expect("eigenvalue");
            assert_eq!(Some(re), s.value);
            assert!(im.abs() < 1e-8 * re, "a lossless ring has a real frequency, got Im f = {im}");
            assert_eq!(scalar(&env, &format!("layer fraction {}", j + 1)), Some(0.0));
        }
    }
}

#[test]
fn an_absorbing_layer_makes_the_modes_complex_with_q_and_a_layer_fraction() {
    // the outer part of the ring section is a layer (domain 2, absorbing towards +r)
    let mut mesh = ring(5.0e-3, 2.0e-3, 0.4e-3, 16, 3);
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        let xc = (0..3).map(|k| mesh.point(v[k] as usize)[0]).sum::<f64>() / 3.0;
        mesh.cell_labels[ci] = if xc > 6.0e-3 { 2 } else { 1 };
    }
    let ex = exchange(&mesh);
    let table = Table;
    let mut s = spec(0, 0.0);
    s.materials.insert(
        "layer".into(),
        MaterialAssign::Props([("E".to_string(), 70e9), ("nu".to_string(), 0.33), ("rho".to_string(), 2700.0)].into()),
    );
    s.bcs.push(BcSpec {
        on: "layer".into(),
        kind: "pml".into(),
        params: json!({ "x": "plus", "strength": 3.0, "order": 2.0 }),
    });
    let r = Sms.solve(&ex, &s, &mut Ctx::new(&table)).unwrap();
    let mut seen_lossy = false;
    for (j, sol) in r.solutions.iter().enumerate() {
        let [re, im] = sol.eigenvalue.expect("eigenvalue");
        assert_eq!(Some(re), sol.value);
        let lf = scalar(&r, &format!("layer fraction {}", j + 1)).expect("layer fraction");
        assert!((0.0..=1.0).contains(&lf), "layer fraction {lf}");
        println!("mode {}: f = {re:.6e} {im:+.3e}i Hz, layer fraction {lf:.3}", j + 1);
        if im != 0.0 {
            seen_lossy = true;
            let q = scalar(&r, &format!("Q{}", j + 1)).expect("Q");
            assert!((q - re / (2.0 * im.abs())).abs() < 1e-9 * q.abs());
        }
    }
    assert!(seen_lossy, "no mode got an imaginary part from the layer");
    // some mode puts weight in the layer
    let fr: Vec<f64> = (1..=r.solutions.len()).map(|j| scalar(&r, &format!("layer fraction {j}")).unwrap()).collect();
    assert!(fr.iter().cloned().fold(0.0, f64::max) > 0.0, "{fr:?}");
}

#[test]
fn the_study_refuses_what_it_cannot_say() {
    let ring_mesh = ring(5.0e-3, 1.0e-3, 0.4e-3, 6, 2);
    let table = Table;
    let solve = |mesh: &SimplexMesh, s: &ProblemSpec| Sms.solve(&exchange(mesh), s, &mut Ctx::new(&table));
    let with = |kind: &str, on: &str, params: serde_json::Value| {
        let mut s = spec(0, 0.0);
        s.bcs.push(BcSpec { on: on.into(), kind: kind.into(), params });
        s
    };
    // 3D mesh
    let block = femcore::testkit::block(1.0, 0.5, 0.4, 3, 2, 2);
    let e = solve(&block, &spec(0, 0.0)).unwrap_err();
    assert!(e.to_string().contains("2D"), "{e}");
    // a mesh in x < 0
    let flat = femcore::testkit::rect(1.0e-3, 0.4e-3, 4, 2);
    let mut pts = flat.points.clone();
    for i in 0..flat.n_points() {
        pts[2 * i] -= 0.5e-3;
    }
    let neg = SimplexMesh::new(
        2,
        pts,
        flat.cells.clone(),
        flat.cell_labels.clone(),
        flat.facets.clone(),
        flat.facet_labels.clone(),
    );
    let e = solve(&neg, &spec(0, 0.0)).unwrap_err();
    assert!(e.to_string().contains("x >= 0"), "{e}");
    // a load, a Floquet pair, a low-reflecting boundary
    let e = solve(&ring_mesh, &with("load", "right", json!({ "traction": [0.0, 0.0, -1e5] }))).unwrap_err();
    assert!(e.to_string().contains("static"), "{e}");
    let e = solve(&ring_mesh, &with("periodic", "left", json!({ "dest": "right", "phase": 0.5 }))).unwrap_err();
    assert!(e.to_string().contains("Floquet"), "{e}");
    let e = solve(&ring_mesh, &with("low_reflecting", "right", json!({}))).unwrap_err();
    assert_eq!(e.code(), "unsupported", "{e}");
    assert!(e.to_string().contains("D1"), "{e}");
    // an absorbing layer along y (not in the plane) or towards the axis
    let e = solve(&ring_mesh, &with("pml", "bar", json!({ "y": "plus" }))).unwrap_err();
    assert!(e.to_string().contains("y"), "{e}");
    let e = solve(&ring_mesh, &with("pml", "bar", json!({ "x": "minus" }))).unwrap_err();
    assert!(e.to_string().contains("axis"), "{e}");
    // a constraint with a value other than zero
    let e = solve(&ring_mesh, &with("fixed", "right", json!({ "ux": 1e-6 }))).unwrap_err();
    assert!(e.to_string().contains("zero"), "{e}");
    // the typed problem is checked the same way
    let mut p = typed(0, 0.0);
    p.bcs.push(Bc { labels: vec![4], kind: BcKind::LowReflecting });
    assert!(nyx_sms::run(&exchange(&ring_mesh), &p, &femcore::progress::Progress::new()).is_err());
    // and a good one passes, so the refusals above are about their own cause
    solve(&ring_mesh, &spec(0, 0.0)).unwrap();
}

#[test]
fn the_study_is_listed_with_its_three_parameters_and_sizes_like_an_eigenfrequency_study() {
    let sd = Sms.studies().into_iter().find(|s| s.kind == "axisymmetric_eigenfrequency").expect("study");
    assert_eq!(sd.dims, Dims::Two);
    let names: Vec<&str> = sd.params.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(names, ["n_modes", "f_search", "azimuthal_order"]);
    assert!(!Sms.fields("axisymmetric_eigenfrequency").is_empty());
    // physics-driven sizing: a search frequency gives a wavelength per active region; none gives nothing
    let sz = Sms.sizing(&spec(0, 715e6), None);
    assert!(sz.wavelengths.contains_key("bar"), "{:?}", sz.wavelengths);
    assert!(Sms.sizing(&spec(0, 0.0), None).wavelengths.is_empty());
}

// ---- bodies that touch the axis (K1 of ADR 69), through the study -------------------------------------------------------------------------

/// `J_n(x)` by the periodic trapezoid rule (spectrally accurate), as in `axisym_axis.rs`.
fn bessel(n: i32, x: f64) -> f64 {
    let steps = 400;
    (0..steps)
        .map(|k| {
            let t = 2.0 * std::f64::consts::PI * k as f64 / steps as f64;
            (n as f64 * t - x * t.sin()).cos()
        })
        .sum::<f64>()
        / steps as f64
}

/// A thin solid disk (mesh from the axis, z from `-t/2` to `t/2`), domain 1, nothing fixed.
fn disk(rr: f64, t: f64) -> SimplexMesh {
    let h = t / 2.0;
    let m = femcore::testkit::parts_mesh(&[femcore::testkit::part(0.0, -t / 2.0, rr, t / 2.0, 1, Some(h))], h);
    SimplexMesh::new(2, m.points, m.cells, m.cell_labels, m.facets, m.facet_labels)
}

#[test]
fn a_solid_disk_touching_the_axis_is_solved_by_the_study_and_matches_the_contour_closed_form() {
    let (e, nu, rho, rr, t) = (70e9_f64, 0.33_f64, 2700.0_f64, 1.0e-3_f64, 2.0e-5_f64);
    // first root of lambda J0(lambda) / J1(lambda) = 1 - nu by bisection on (1, 3), then f = lambda / (2 pi R) sqrt(E / (rho (1 - nu^2)))
    let g = |x: f64| x * bessel(0, x) / bessel(1, x) - (1.0 - nu);
    let (mut a, mut b) = (1.0_f64, 3.0_f64);
    for _ in 0..60 {
        let mid = 0.5 * (a + b);
        if g(a) * g(mid) <= 0.0 {
            b = mid;
        } else {
            a = mid;
        }
    }
    let want = 0.5 * (a + b) / (2.0 * std::f64::consts::PI * rr) * (e / (rho * (1.0 - nu * nu))).sqrt();
    let mesh = disk(rr, t);
    let ex = exchange(&mesh);
    let table = Table;
    let mut s = spec(0, want);
    s.bcs.clear();
    let env = Sms.solve(&ex, &s, &mut Ctx::new(&table)).unwrap();
    let mut p = typed(0, want);
    p.bcs.clear();
    let typ = nyx_sms::run(&ex, &p, &femcore::progress::Progress::new()).unwrap();
    let mut plain = p.clone();
    plain.study = Study::Eigenfrequency { n_modes: 4, f_search: want };
    let (lib, _) = nyx_sms::solve_solid_axisym(&mesh, &plain, 0.0).unwrap();
    assert_eq!(frequencies(&env), frequencies(&typ));
    assert_eq!(frequencies(&typ), frequencies(&lib));
    let got = frequencies(&env).into_iter().min_by(|a, b| (a - want).abs().total_cmp(&(b - want).abs())).unwrap();
    println!("disk on the axis through the study: {got:.6e} Hz against the contour closed form {want:.6e} Hz");
    assert!((got / want - 1.0).abs() < 1e-3, "{got} against {want}");
    assert!(env.log.iter().any(|l| l.to_lowercase().contains("axis")), "the log names the axis condition: {:?}", env.log);
}
