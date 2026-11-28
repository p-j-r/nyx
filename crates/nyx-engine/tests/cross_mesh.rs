// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Per-physics meshing, active-domain subsets and domain-aware cross-mesh transfer through the engine (ADR 135). Two stand-in physics on the
//! SOI strip of the reduction tests: `src` meshes everything at the project element size and exports a quadratic temperature (`src.T`) and a
//! field that jumps 12:1 across the silicon / oxide interface; `dst` imports the temperature and is meshed on its own, coarser, with only the
//! oxide and silicon regions active. The quadratic is held exactly by `P2`, so every transferred value and every integral has a closed form.
use nyx_engine::*;
use nyx_physics::*;
use serde_json::{json, Value};
use std::sync::Arc;

const SILICON: ([f64; 2], [f64; 2]) = ([-2.5e-7, 0.0], [2.5e-7, 2.2e-7]);
const OXIDE: ([f64; 2], [f64; 2]) = ([-2e-6, -1.5e-6], [2e-6, 0.0]);

/// The quadratic both stand-ins agree on, in micrometres.
fn poly(x: &[f64]) -> f64 {
    let (a, c) = (x[0] * 1e6, x[1] * 1e6);
    1.0 + 2.0 * a - 3.0 * c + 0.5 * a * a + a * c - c * c
}

/// The coordinates of local node `a` of cell `c`: the vertices, then the midpoints of the vertex pairs.
fn node(sm: &Discretisation, c: usize, a: usize) -> Vec<f64> {
    let (v, dim) = (sm.cell(c), sm.dim);
    if a <= dim {
        return sm.point(v[a] as usize).to_vec();
    }
    let mut k = dim + 1;
    for i in 0..=dim {
        for j in i + 1..=dim {
            if k == a {
                return (0..dim).map(|d| 0.5 * (sm.point(v[i] as usize)[d] + sm.point(v[j] as usize)[d])).collect();
            }
            k += 1;
        }
    }
    sm.point(v[0] as usize).to_vec()
}

struct Stand {
    id: &'static str,
}

impl Physics for Stand {
    fn id(&self) -> &'static str {
        self.id
    }
    fn name(&self) -> &'static str {
        "Cross-mesh stand-in"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc {
            kind: "poly".into(),
            name: "Polynomial".into(),
            doc: "".into(),
            dims: Dims::Both,
            params: vec![],
        }]
    }
    fn bc_kinds(&self) -> Vec<BcDesc> {
        vec![]
    }
    fn properties(&self) -> Vec<PropDesc> {
        vec![]
    }
    fn params(&self) -> Vec<ParamDesc> {
        vec![]
    }
    fn fields(&self, _study: &str) -> Vec<FieldDesc> {
        vec![FieldDesc::new("T", "K", &[], "")]
    }
    fn exports(&self) -> Vec<FieldLink> {
        if self.id == "src" {
            vec![FieldLink::new("src.T", "temperature")]
        } else {
            vec![]
        }
    }
    fn imports(&self) -> Vec<FieldLink> {
        if self.id == "dst" {
            vec![FieldLink::new("src.T", "temperature")]
        } else {
            vec![]
        }
    }
    fn validate(&self, _s: &ProblemSpec, _m: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, _s: &ProblemSpec, _library: Option<&Arc<femcore::material::Library>>) -> Sizing {
        Sizing::default()
    }
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        let sm = discretise(mesh)?;
        let order = spec.order as usize;
        let nloc = sm.nloc(order);
        let silicon = ctx.resolver.domains("silicon")?;
        let air = ctx.resolver.domains("air")?;
        let t: Vec<f64> = sm.dof_coordinates(order).chunks(sm.dim).map(poly).collect();
        let (mut p, mut one, mut jump) = (vec![], vec![], vec![]);
        for c in 0..sm.n_cells() {
            let k = if silicon.contains(&sm.cell_labels[c]) { 12.0 } else { 1.0 };
            for a in 0..nloc {
                p.push(poly(&node(&sm, c, a)));
                one.push(1.0);
                jump.push(k);
            }
        }
        let cf = |name: &str, re: Vec<f64>| CellField { name: name.into(), unit: "a.u.".into(), nloc, re, im: vec![] };
        let mut scalars = vec![
            ("cells".to_string(), Scalar::Num(sm.n_cells() as f64)),
            ("air cells".to_string(), Scalar::Num(sm.cell_labels.iter().filter(|l| air.contains(l)).count() as f64)),
        ];
        if self.id == "dst" {
            let t_in = ctx.inputs.get("temperature").ok_or_else(|| Error::invalid("no imported temperature"))?;
            if t_in.len() != sm.n_dofs(order) {
                return Err(Error::invalid(format!("{} imported values for {} dofs", t_in.len(), sm.n_dofs(order))));
            }
            let err = t_in.iter().zip(&t).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            scalars.push(("input error".into(), Scalar::Num(err)));
            let tf = ctx.transferred_inputs.get("temperature");
            scalars.push(("transferred".into(), Scalar::Num(tf.is_some() as u8 as f64)));
            scalars.push(("same mesh".into(), Scalar::Num(tf.is_some_and(|f| f.is_same_mesh(mesh)) as u8 as f64)));
        }
        Ok(SolveResult {
            physics: self.id.into(),
            study: "poly".into(),
            order: spec.order,
            dim: sm.dim as u8,
            dofs: sm.n_dofs(order),
            seconds: 0.0,
            solutions: vec![Solution {
                label: "only".into(),
                fields: vec![Field::scalar("T", "K", t)],
                value: Some(0.0),
                cell_fields: vec![cf("p", p), cf("one", one), cf("jump", jump)],
                ..Default::default()
            }],
            scalars,
            log: vec![],
            ..Default::default()
        })
    }
}

fn call(e: &mut Engine, name: &str, args: Value) -> Value {
    e.call(name, &args).unwrap_or_else(|err| panic!("{name} failed: {err}"))
}

fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

/// The exact integral of a polynomial over a box: three-point Gauss-Legendre per axis.
fn exact(lo: &[f64; 2], hi: &[f64; 2], f: &dyn Fn(&[f64]) -> f64) -> f64 {
    let r = (3.0f64 / 5.0).sqrt();
    let nodes = [(-r, 5.0 / 9.0), (0.0, 8.0 / 9.0), (r, 5.0 / 9.0)];
    let mut total = 0.0;
    for (tx, wx) in nodes {
        for (tz, wz) in nodes {
            let x = [
                0.5 * (lo[0] + hi[0]) + 0.5 * (hi[0] - lo[0]) * tx,
                0.5 * (lo[1] + hi[1]) + 0.5 * (hi[1] - lo[1]) * tz,
            ];
            total += 0.25 * (hi[0] - lo[0]) * (hi[1] - lo[1]) * wx * wz * f(&x);
        }
    }
    total
}

fn close(a: f64, b: f64, rel: f64) {
    assert!((a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-300), "{a} vs {b}");
}

/// The strip with `src` on the project mesh and `dst` on its own coarser mesh of the oxide and the silicon only.
fn strip(own_mesh: bool) -> Engine {
    strip_orders(own_mesh, 2, 2)
}

fn strip_orders(own_mesh: bool, src_order: u8, dst_order: u8) -> Engine {
    let mut r = Registry::new();
    r.register(Arc::new(Stand { id: "src" })).unwrap();
    r.register(Arc::new(Stand { id: "dst" })).unwrap();
    let mut e = Engine::with(r, Arc::new(FemcoreMesher));
    call(&mut e, "project_new", json!({ "title": "Strip", "h": 5e-8 }));
    call(&mut e, "add_rect", json!({ "name": "window", "min": [-2e-6, -1.5e-6], "max": [2e-6, 1.5e-6] }));
    call(&mut e, "add_rect", json!({ "name": "box", "min": [-2e-6, -1.5e-6], "max": [2e-6, 0.0] }));
    call(&mut e, "add_rect", json!({ "name": "core", "min": [-2.5e-7, 0.0], "max": [2.5e-7, 2.2e-7] }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "air", "input": "window" }, { "name": "oxide", "input": "box" }, { "name": "silicon", "input": "core" }] }),
    );
    call(&mut e, "add_physics", json!({ "physics": "src", "study": "poly", "order": src_order }));
    call(&mut e, "add_physics", json!({ "physics": "dst", "study": "poly", "order": dst_order }));
    if own_mesh {
        call(&mut e, "set_mesh", json!({ "physics": "dst", "h": 8e-8, "active_domains": ["oxide", "silicon"] }));
        let p = e.project().clone();
        assert!(p.has_own_mesh("dst") && !p.has_own_mesh("src"));
        assert_eq!(p.mesh_for("dst").active_domains, Some(vec!["oxide".to_string(), "silicon".to_string()]));
        assert_eq!(p.mesh_for("src"), &p.mesh);
        let back = Project::from_json(&p.to_json()).unwrap();
        assert_eq!(back, p, "mesh_by_physics and active_domains survive a save");
        assert!(e.call("set_mesh", &json!({ "physics": "nope", "h": 1e-7 })).is_err());
    }
    call(&mut e, "solve_run", json!({ "physics": ["src", "dst"] }));
    e
}

fn scalar(e: &Engine, physics: &str, name: &str) -> f64 {
    e.result(physics)
        .unwrap_or_else(|| panic!("no result for {physics}"))
        .0
        .scalar(name)
        .unwrap_or_else(|| panic!("no {name}"))
}

#[test]
fn a_physics_on_its_own_active_domain_mesh_imports_the_field_by_region_name() {
    let e = strip(true);
    // dst was meshed on its own: coarser, and without the air
    assert_eq!(scalar(&e, "dst", "air cells"), 0.0);
    assert!(scalar(&e, "src", "air cells") > 0.0);
    assert!(scalar(&e, "dst", "cells") < scalar(&e, "src", "cells"));
    // the temperature came across meshes, sampled at dst's own P2 dofs: the quadratic is exact there
    assert_eq!(scalar(&e, "dst", "transferred"), 1.0);
    assert_eq!(scalar(&e, "dst", "same mesh"), 0.0);
    assert!(scalar(&e, "dst", "input error") < 1e-9, "transfer error {}", scalar(&e, "dst", "input error"));
}

#[test]
fn one_mesh_for_every_physics_keeps_the_zero_cost_path() {
    let e = strip(false);
    assert_eq!(scalar(&e, "dst", "transferred"), 0.0);
    assert_eq!(scalar(&e, "dst", "input error"), 0.0);
    assert_eq!(scalar(&e, "dst", "cells"), scalar(&e, "src", "cells"));
}

fn mul(a: Value, b: Value) -> Value {
    json!({ "op": "mul", "args": [a, b] })
}

fn leaf(name: &str, physics: &str) -> Value {
    json!({ "field": name, "physics": physics, "solution": 0 })
}

#[test]
fn reductions_read_fields_of_two_meshes_exactly_and_keep_the_jump() {
    let mut e = strip(true);
    // dst's mesh is primary (first leaf); src's quadratic is sampled at its P2 nodes: the degree-4 product is integrated exactly
    let args = json!({ "expr": mul(leaf("p", "dst"), leaf("T", "src")), "region": "silicon" });
    let got = call(&mut e, "integrate", args.clone());
    assert_eq!(got["cross_mesh_transfer"], json!(true));
    assert_eq!(got["cross_mesh"]["missing_nodes"], json!(0));
    close(num(&got["value"]), exact(&SILICON.0, &SILICON.1, &|x| poly(x) * poly(x)), 1e-9);
    // the 12:1 jump of src is read only in source cells of the same region: no smearing into the oxide
    for (region, k, b) in [("silicon", 12.0, SILICON), ("oxide", 1.0, OXIDE)] {
        let got = call(
            &mut e,
            "integrate",
            json!({ "expr": mul(leaf("one", "dst"), leaf("jump", "src")), "region": region }),
        );
        let area = (b.1[0] - b.0[0]) * (b.1[1] - b.0[1]);
        close(num(&got["value"]), k * area, 1e-12);
    }
    // src primary: dst has no air, so its field reads 0 there and the transfer says how many nodes had no source region
    let got = call(&mut e, "integrate", json!({ "expr": mul(leaf("one", "src"), leaf("one", "dst")) }));
    assert!(num(&got["cross_mesh"]["missing_nodes"]) > 0.0);
    let area = |b: ([f64; 2], [f64; 2])| (b.1[0] - b.0[0]) * (b.1[1] - b.0[1]);
    close(num(&got["value"]), area(SILICON) + area(OXIDE), 1e-9);
    // field_op and overlap take the same path
    let h = call(&mut e, "field_op", json!({ "op": "mul", "args": [leaf("p", "dst"), leaf("T", "src")] }));
    assert_eq!(h["cross_mesh_transfer"], json!(true));
    let ov = call(&mut e, "overlap", json!({ "a": leaf("p", "dst"), "b": leaf("T", "src"), "region": "silicon" }));
    assert_eq!(ov["cross_mesh_transfer"], json!(true));
    // direct quadrature-point sampling: even when dst is P1 (linear) on a coarse mesh, src's P2 quadratic is integrated exactly
    let mut e_p1 = strip_orders(true, 2, 1);
    let q_p1 =
        call(&mut e_p1, "integrate", json!({ "expr": mul(leaf("one", "dst"), leaf("T", "src")), "region": "silicon" }));
    assert_eq!(q_p1["cross_mesh_transfer"], json!(true));
    close(num(&q_p1["value"]), exact(&SILICON.0, &SILICON.1, &poly), 1e-9);
    // the same numbers on one thread and on the pool
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| call(&mut e, "integrate", args.clone()));
    assert_eq!(num(&one["value"]).to_bits(), num(&got_value(&mut e, &args)).to_bits());
    println!("cross-mesh checksum: {:016x}", num(&one["value"]).to_bits());
}

fn got_value(e: &mut Engine, args: &Value) -> Value {
    call(e, "integrate", args.clone())["value"].clone()
}

#[test]
fn cross_mesh_chain_and_reductions_are_bit_exact_between_nyx_threads_1_and_8() {
    if let Ok(out_path) = std::env::var("NYX_FR44_DET_OUT") {
        let mut e = strip(true);
        let mut words: Vec<u64> = Vec::new();
        for code in ["src", "dst"] {
            let (res, _) = e.result(code).unwrap();
            words.push(res.dofs as u64);
            for (_, s) in &res.scalars {
                if let Scalar::Num(v) = s {
                    words.push(v.to_bits());
                }
            }
            for f in &res.solutions[0].fields {
                words.extend(f.values.iter().map(|v| v.to_bits()));
            }
        }
        let ig =
            call(&mut e, "integrate", json!({ "expr": mul(leaf("p", "dst"), leaf("T", "src")), "region": "silicon" }));
        words.push(num(&ig["value"]).to_bits());
        let ov = call(&mut e, "overlap", json!({ "a": leaf("p", "dst"), "b": leaf("T", "src"), "region": "silicon" }));
        words.push(num(&ov["coefficient"]).to_bits());
        let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
        std::fs::write(&out_path, bytes).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_child = |nt: &str| -> Vec<u8> {
        let out_path = std::env::temp_dir().join(format!("nyx-fr44-det-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args([
                "--exact",
                "cross_mesh_chain_and_reductions_are_bit_exact_between_nyx_threads_1_and_8",
                "--nocapture",
            ])
            .env("NYX_THREADS", nt)
            .env("NYX_FR44_DET_OUT", &out_path)
            .output()
            .expect("spawn cross_mesh determinism subprocess");
        assert!(
            out.status.success(),
            "child NYX_THREADS={nt} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        bytes
    };

    let b1 = run_child("1");
    let b8 = run_child("8");
    assert!(!b1.is_empty());
    assert_eq!(b1, b8, "cross-mesh solve and reductions must be 0-ULP bit-identical for NYX_THREADS=1 and 8");
}
