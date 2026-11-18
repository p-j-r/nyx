//! The whole-array reductions (productivity plan N1, ADR 90) against two oracles. A stand-in physics `fld` writes polynomial fields of degree 2
//! at the nodes of the mesh, which a P2 field holds exactly, so every integral has a closed form (Gauss-Legendre over a box, which is exact for
//! these degrees). And each reduction is written out again here as a plain loop on the raw arrays of the result, in the order the module
//! documents, and must agree to the bit. The point location of `probe` is checked against the closed form and a loop to 1e-12 only, because the
//! barycentric inversion is not specified bit for bit.
#![allow(clippy::needless_range_loop)] // the reference loops are plain index loops on purpose: they are the oracle
use femcore::fem::{geometry, quad_rule, shape};
use nyx_engine::*;
use nyx_physics::*;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;

// ------------------------------------------------------------------------------------------------------ the stand-in physics

/// The polynomials the stand-in writes, in units of `1/scale` metres. `u` are the scaled coordinates (2D: x, z; 3D: x, y, z).
fn polys(u: &[f64]) -> (f64, (f64, f64), (f64, f64)) {
    let (a, b, c) = if u.len() == 2 { (u[0], 0.0, u[1]) } else { (u[0], u[1], u[2]) };
    let p = 1.0 + 2.0 * a - 3.0 * c + 0.5 * a * a + a * c - c * c + 0.25 * b;
    let ex = (a + 2.0 * c + 0.5 * b, a * a - c);
    let ez = (1.0 + 0.5 * a * c, 0.25 * a - c + 0.5 * b * c);
    (p, ex, ez)
}

struct Fld {
    scale: f64,
}

impl Physics for Fld {
    fn id(&self) -> &'static str {
        "fld"
    }
    fn name(&self) -> &'static str {
        "Field stand-in"
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
        vec![]
    }
    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn validate(&self, _s: &ProblemSpec, _m: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, _s: &ProblemSpec, _library: Option<&Arc<femcore::material::Library>>) -> Sizing {
        Sizing::default()
    }
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, _ctx: &mut Ctx) -> Result<SolveResult, Error> {
        let sm = mesh.discretise().map_err(Error::Invalid)?;
        let (dim, order) = (sm.dim, spec.order as usize);
        let nloc = sm.nloc(order);
        let scaled = |p: &[f64]| -> Vec<f64> { p.iter().map(|x| x * self.scale).collect() };
        // the coordinates of local node a of cell c: the vertices, then the midpoints of the vertex pairs
        let node = |c: usize, a: usize| -> Vec<f64> {
            let v = sm.cell(c);
            if a <= dim {
                return sm.point(v[a] as usize).to_vec();
            }
            let mut k = dim + 1;
            for i in 0..=dim {
                for j in i + 1..=dim {
                    if k == a {
                        return (0..dim)
                            .map(|d| 0.5 * (sm.point(v[i] as usize)[d] + sm.point(v[j] as usize)[d]))
                            .collect();
                    }
                    k += 1;
                }
            }
            unreachable!()
        };
        let (mut p, mut one, mut exr, mut exi, mut ezr, mut ezi) = (vec![], vec![], vec![], vec![], vec![], vec![]);
        for c in 0..sm.n_cells() {
            for a in 0..nloc {
                let (pv, ex, ez) = polys(&scaled(&node(c, a)));
                p.push(pv);
                one.push(1.0);
                exr.push(ex.0);
                exi.push(ex.1);
                ezr.push(ez.0);
                ezi.push(ez.1);
            }
        }
        let cf =
            |name: &str, re: Vec<f64>, im: Vec<f64>| CellField { name: name.into(), unit: "a.u.".into(), nloc, re, im };
        let mut t: Vec<f64> = (0..sm.n_points()).map(|i| polys(&scaled(sm.point(i))).0).collect();
        for &(a, b) in sm.edges() {
            let mid: Vec<f64> = (0..dim).map(|d| 0.5 * (sm.point(a as usize)[d] + sm.point(b as usize)[d])).collect();
            t.push(polys(&scaled(&mid)).0);
        }
        t.truncate(sm.n_dofs(order));
        Ok(SolveResult {
            physics: "fld".into(),
            study: "poly".into(),
            order: spec.order,
            dim: dim as u8,
            dofs: sm.n_dofs(order),
            seconds: 0.0,
            solutions: vec![Solution {
                label: "only".into(),
                fields: vec![Field::scalar("T", "K", t)],
                value: Some(0.0),
                cell_fields: vec![cf("p", p, vec![]), cf("one", one, vec![]), cf("Ex", exr, exi), cf("Ez", ezr, ezi)],
                ..Default::default()
            }],
            scalars: vec![],
            log: vec![],
        })
    }
}

fn engine(scale: f64) -> Engine {
    let mut r = Registry::new();
    r.register(Arc::new(Fld { scale })).unwrap();
    Engine::with(r, Arc::new(FemcoreMesher))
}

fn call(e: &mut Engine, name: &str, args: Value) -> Value {
    e.call(name, &args).unwrap_or_else(|err| panic!("{name} failed: {err}"))
}

fn err(e: &mut Engine, name: &str, args: Value) -> String {
    match e.call(name, &args) {
        Ok(v) => panic!("{name} should have failed, gave {v}"),
        Err(x) => x.to_string(),
    }
}

/// The SOI strip of the engine tests (window 4 um x 3 um; oxide below z = 0; the silicon core 0.5 um x 0.22 um), meshed finely enough to
/// have several blocks of cells, with the stand-in solved at order 2.
fn strip() -> Engine {
    let mut e = engine(1e6);
    call(&mut e, "project_new", json!({ "title": "Strip", "h": 5e-8 }));
    call(&mut e, "add_rect", json!({ "name": "window", "min": [-2e-6, -1.5e-6], "max": [2e-6, 1.5e-6] }));
    call(&mut e, "add_rect", json!({ "name": "box", "min": [-2e-6, -1.5e-6], "max": [2e-6, 0.0] }));
    call(&mut e, "add_rect", json!({ "name": "core", "min": [-2.5e-7, 0.0], "max": [2.5e-7, 2.2e-7] }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "air", "input": "window" }, { "name": "oxide", "input": "box" }, { "name": "silicon", "input": "core" }] }),
    );
    call(&mut e, "add_physics", json!({ "physics": "fld", "study": "poly", "order": 2 }));
    call(&mut e, "solve_run", json!({}));
    e
}

// ------------------------------------------------------------------------------------------------------ oracles

/// The exact integral of a polynomial integrand over a box: three-point Gauss-Legendre per axis (exact to degree 5 per axis).
fn exact(lo: &[f64], hi: &[f64], f: &dyn Fn(&[f64]) -> f64) -> f64 {
    let r = (3.0f64 / 5.0).sqrt();
    let nodes = [(-r, 5.0 / 9.0), (0.0, 8.0 / 9.0), (r, 5.0 / 9.0)];
    let dim = lo.len();
    let mut total = 0.0;
    let count = 3usize.pow(dim as u32);
    for k in 0..count {
        let (mut x, mut w, mut kk) = (vec![0.0; dim], 1.0, k);
        for d in 0..dim {
            let (t, wt) = nodes[kk % 3];
            kk /= 3;
            x[d] = 0.5 * (lo[d] + hi[d]) + 0.5 * (hi[d] - lo[d]) * t;
            w *= 0.5 * (hi[d] - lo[d]) * wt;
        }
        total += w * f(&x);
    }
    total
}

fn close(a: f64, b: f64, rel: f64) {
    assert!((a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-300), "{a} vs {b}");
}

fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

fn pair(v: &Value) -> (f64, f64) {
    (num(&v[0]), num(&v[1]))
}

/// The integrand of a reference loop: the interpolated (re, im) of each named field at a point in, one (re, im) value out.
type Integrand<'a> = &'a dyn Fn(&[(f64, f64)]) -> (f64, f64);

/// The raw arrays a reference loop reads.
struct Raw {
    mesh: Arc<femcore::mesh::SimplexMesh>,
    sol: Solution,
    order: usize,
    nloc: usize,
    mask: Vec<bool>,
}

fn raw(e: &Engine, region: Option<&str>) -> Raw {
    let (res, _) = e.result("fld").expect("a result");
    let model = e.model().expect("a mesh").clone();
    let mesh = model.solver_mesh();
    let order = res.order as usize;
    let mask = match region {
        None => vec![true; mesh.n_cells()],
        Some(r) => {
            let labels = model.domains(r).unwrap();
            mesh.cell_labels.iter().map(|l| labels.contains(l)).collect()
        }
    };
    Raw { nloc: mesh.nloc(order), mesh, sol: res.solutions[0].clone(), order, mask }
}

impl Raw {
    fn cf(&self, name: &str) -> &CellField {
        self.sol.cell_fields.iter().find(|f| f.name == name).unwrap()
    }

    /// `sum_c vol * sum_q w_q f_q` over blocks of 2048 cells, `f` given the interpolated (re, im) of each named field at the point.
    fn integral(&self, names: &[&str], f: Integrand) -> (f64, f64) {
        let dim = self.mesh.dim;
        let quad = quad_rule(dim, 4);
        let ns: Vec<Vec<f64>> = quad.lam.iter().map(|l| shape(l, self.order).0).collect();
        let vol = geometry(&self.mesh).vol;
        let fields: Vec<&CellField> = names.iter().map(|n| self.cf(n)).collect();
        let (mut tre, mut tim) = (0.0, 0.0);
        for block in 0..self.mesh.n_cells().div_ceil(2048) {
            let (mut bre, mut bim) = (0.0, 0.0);
            for c in block * 2048..((block + 1) * 2048).min(self.mesh.n_cells()) {
                if !self.mask[c] {
                    continue;
                }
                let (mut cre, mut cim) = (0.0, 0.0);
                for q in 0..quad.w.len() {
                    let at: Vec<(f64, f64)> = fields
                        .iter()
                        .map(|fl| {
                            let (mut r, mut i) = (0.0, 0.0);
                            for a in 0..self.nloc {
                                r += ns[q][a] * fl.re[c * self.nloc + a];
                                if !fl.im.is_empty() {
                                    i += ns[q][a] * fl.im[c * self.nloc + a];
                                }
                            }
                            (r, i)
                        })
                        .collect();
                    let v = f(&at);
                    cre += quad.w[q] * v.0;
                    cim += quad.w[q] * v.1;
                }
                bre += vol[c] * cre;
                bim += vol[c] * cim;
            }
            tre += bre;
            tim += bim;
        }
        (tre, tim)
    }
}

const SILICON: ([f64; 2], [f64; 2]) = ([-2.5e-7, 0.0], [2.5e-7, 2.2e-7]);
const OXIDE: ([f64; 2], [f64; 2]) = ([-2e-6, -1.5e-6], [2e-6, 0.0]);
const WINDOW: ([f64; 2], [f64; 2]) = ([-2e-6, -1.5e-6], [2e-6, 1.5e-6]);

fn u(x: &[f64]) -> Vec<f64> {
    x.iter().map(|v| v * 1e6).collect()
}

// ------------------------------------------------------------------------------------------------------ tests

#[test]
fn integrals_equal_their_closed_forms_and_their_loops_to_the_bit() {
    let mut e = strip();
    let cells = e.model().unwrap().mesh().n_cells();
    assert!(cells > 2 * 2048, "the mesh must span several blocks of cells: {cells}");
    // a real field over the whole window and over two regions
    for (region, (lo, hi)) in [(None, WINDOW), (Some("oxide"), OXIDE), (Some("silicon"), SILICON)] {
        let mut args = json!({ "expr": "p", "physics": "fld" });
        if let Some(r) = region {
            args["region"] = json!(r);
        }
        let got = call(&mut e, "integrate", args);
        assert_eq!(got["complex"], json!(false));
        let closed = exact(&lo, &hi, &|x| polys(&u(x)).0);
        close(num(&got["value"]), closed, 1e-9);
        let area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
        close(num(&got["measure"]), area, 1e-9);
        let reference = raw(&e, region).integral(&["p"], &|a| (a[0].0, 0.0));
        assert_eq!(num(&got["value"]).to_bits(), reference.0.to_bits(), "{region:?}");
        assert_eq!(got["stale"], json!(false));
    }
    // |Ex|^2 + |Ez|^2 (degree 4: exact with the rule), as a handle built from commands
    let ex = call(&mut e, "field", json!({ "physics": "fld", "name": "Ex" }))["handle"].clone();
    let ez = call(&mut e, "field", json!({ "physics": "fld", "name": "Ez" }))["handle"].clone();
    let a = call(&mut e, "field_op", json!({ "op": "abs2", "args": [ex.clone()] }))["handle"].clone();
    let b = call(&mut e, "field_op", json!({ "op": "abs2", "args": [ez.clone()] }))["handle"].clone();
    let sum = call(&mut e, "field_op", json!({ "op": "add", "args": [a, b] }));
    assert_eq!(sum["complex"], json!(false));
    let got = call(&mut e, "integrate", json!({ "expr": sum["handle"], "region": "silicon" }));
    let closed = exact(&SILICON.0, &SILICON.1, &|x| {
        let (_, ex, ez) = polys(&u(x));
        ex.0 * ex.0 + ex.1 * ex.1 + ez.0 * ez.0 + ez.1 * ez.1
    });
    close(num(&got["value"]), closed, 1e-9);
    let reference = raw(&e, Some("silicon"))
        .integral(&["Ex", "Ez"], &|a| (a[0].0 * a[0].0 + a[0].1 * a[0].1 + (a[1].0 * a[1].0 + a[1].1 * a[1].1), 0.0));
    assert_eq!(num(&got["value"]).to_bits(), reference.0.to_bits());
    // a complex product: the result is a pair, and both parts agree to the bit
    let got = call(&mut e, "integrate", json!({ "expr": { "op": "mul", "args": ["Ex", "Ez"] }, "physics": "fld" }));
    assert_eq!(got["complex"], json!(true));
    let (re, im) = pair(&got["value"]);
    let closed_re = exact(&WINDOW.0, &WINDOW.1, &|x| {
        let (_, ex, ez) = polys(&u(x));
        ex.0 * ez.0 - ex.1 * ez.1
    });
    let closed_im = exact(&WINDOW.0, &WINDOW.1, &|x| {
        let (_, ex, ez) = polys(&u(x));
        ex.0 * ez.1 + ex.1 * ez.0
    });
    close(re, closed_re, 1e-8);
    close(im, closed_im, 1e-8);
    let reference = raw(&e, None)
        .integral(&["Ex", "Ez"], &|a| (a[0].0 * a[1].0 - a[0].1 * a[1].1, a[0].0 * a[1].1 + a[0].1 * a[1].0));
    assert_eq!((re.to_bits(), im.to_bits()), (reference.0.to_bits(), reference.1.to_bits()));
    // the same through the node-averaged flat field `T` (one value per dof): the dof path
    let got = call(&mut e, "integrate", json!({ "expr": "T", "physics": "fld", "region": "oxide" }));
    close(num(&got["value"]), exact(&OXIDE.0, &OXIDE.1, &|x| polys(&u(x)).0), 1e-9);
}

#[test]
fn the_sums_do_not_depend_on_the_number_of_threads() {
    let mut e = strip();
    let args = json!({ "expr": { "op": "mul", "args": ["Ex", { "op": "conj", "args": ["Ez"] }] }, "physics": "fld" });
    let many = call(&mut e, "integrate", args.clone());
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| call(&mut e, "integrate", args.clone()));
    let three =
        rayon::ThreadPoolBuilder::new().num_threads(3).build().unwrap().install(|| call(&mut e, "integrate", args));
    assert_eq!(many["value"], one["value"]);
    assert_eq!(many["value"], three["value"]);
}

#[test]
fn extrema_probes_and_line_cuts_read_the_polynomial() {
    let mut e = strip();
    // the largest and smallest p at the nodes of the silicon: the loop on the raw arrays, and the point it reports
    for (op, max) in [("max", true), ("min", false)] {
        let got = call(&mut e, "extremum", json!({ "expr": "p", "physics": "fld", "region": "silicon", "op": op }));
        let r = raw(&e, Some("silicon"));
        let f = r.cf("p");
        let mut best: Option<(f64, usize, usize)> = None;
        for c in 0..r.mesh.n_cells() {
            if !r.mask[c] {
                continue;
            }
            for a in 0..r.nloc {
                let v = f.re[c * r.nloc + a];
                if best.is_none_or(|(b, _, _)| if max { v > b } else { v < b }) {
                    best = Some((v, c, a));
                }
            }
        }
        let (v, c, a) = best.unwrap();
        assert_eq!(num(&got["value"]).to_bits(), v.to_bits());
        assert_eq!((got["cell"].as_u64().unwrap() as usize, got["node"].as_u64().unwrap() as usize), (c, a));
        let pt: Vec<f64> = got["point"].as_array().unwrap().iter().map(num).collect();
        close(polys(&u(&pt)).0, v, 1e-9);
    }
    // probes: inside each region and outside the window
    let pts = [[1e-7, 5e-8], [-1e-6, -7e-7], [1.2e-6, 9e-7], [5e-6, 0.0]];
    let got = call(&mut e, "probe", json!({ "expr": "p", "physics": "fld", "points": pts }));
    assert_eq!(got["inside"], json!([true, true, true, false]));
    assert_eq!(got["values"][3], Value::Null);
    for (i, p) in pts.iter().take(3).enumerate() {
        close(num(&got["values"][i]), polys(&u(p)).0, 1e-9);
    }
    // a complex field answers `values` (re) and `im`
    let got = call(&mut e, "probe", json!({ "expr": "Ex", "physics": "fld", "points": [pts[0]] }));
    let (_, ex, _) = polys(&u(&pts[0]));
    close(num(&got["values"][0]), ex.0, 1e-9);
    close(num(&got["im"][0]), ex.1, 1e-9);
    // a line cut crossing the oxide, the silicon and the air, with a few points outside
    let cut = call(
        &mut e,
        "line_cut",
        json!({ "expr": "p", "physics": "fld", "from": [0.0, -1.4e-6], "to": [0.0, 1.4e-6], "n": 57 }),
    );
    assert_eq!(cut["inside"].as_array().unwrap().len(), 57);
    assert!(cut["inside"].as_array().unwrap().iter().all(|b| *b == json!(true)));
    close(num(&cut["s"][56]), 2.8e-6, 1e-12);
    for i in 0..57 {
        let z = -1.4e-6 + 2.8e-6 * (i as f64 / 56.0);
        close(num(&cut["values"][i]), polys(&u(&[0.0, z])).0, 1e-9);
    }
    let out = call(
        &mut e,
        "line_cut",
        json!({ "expr": "p", "physics": "fld", "from": [1.5e-6, 0.0], "to": [3.0e-6, 0.0], "n": 4 }),
    );
    assert_eq!(out["inside"], json!([true, true, false, false]));
}

#[test]
fn overlap_and_confinement_follow_their_definitions() {
    let mut e = strip();
    // the fraction of |Ex|^2 + |Ez|^2 in the silicon is the ratio of two `integrate` answers, to the bit (the plan's example)
    let power = json!({ "op": "add", "args": [{ "op": "abs2", "args": ["Ex"] }, { "op": "abs2", "args": ["Ez"] }] });
    let inside = call(&mut e, "integrate", json!({ "expr": power, "physics": "fld", "region": "silicon" }));
    let total = call(&mut e, "integrate", json!({ "expr": power, "physics": "fld" }));
    let conf = call(&mut e, "confinement", json!({ "fields": ["Ex", "Ez"], "physics": "fld", "region": "silicon" }));
    assert_eq!(num(&conf["fraction"]).to_bits(), (num(&inside["value"]) / num(&total["value"])).to_bits());
    let closed_in = exact(&SILICON.0, &SILICON.1, &|x| {
        let (_, ex, ez) = polys(&u(x));
        ex.0 * ex.0 + ex.1 * ex.1 + ez.0 * ez.0 + ez.1 * ez.1
    });
    let closed_all = exact(&WINDOW.0, &WINDOW.1, &|x| {
        let (_, ex, ez) = polys(&u(x));
        ex.0 * ex.0 + ex.1 * ex.1 + ez.0 * ez.0 + ez.1 * ez.1
    });
    close(num(&conf["fraction"]), closed_in / closed_all, 1e-9);
    // a uniform field: the fraction is the ratio of the areas
    let one = call(&mut e, "confinement", json!({ "fields": "one", "physics": "fld", "region": "oxide" }));
    close(num(&one["fraction"]), (4e-6 * 1.5e-6) / (4e-6 * 3e-6), 1e-9);
    // the overlap of a field with itself has coefficient 1; with 2i times itself, the overlap is 2i times the norm
    let same =
        call(&mut e, "overlap", json!({ "a": ["Ex", "Ez"], "b": ["Ex", "Ez"], "physics": "fld", "region": "silicon" }));
    close(num(&same["coefficient"]), 1.0, 1e-12);
    close(num(&same["norm_a"]), closed_in, 1e-9);
    let twice = json!([{ "op": "mul", "args": [{ "num": [0.0, 2.0] }, "Ex"] }, { "op": "mul", "args": [{ "num": [0.0, 2.0] }, "Ez"] }]);
    let scaled =
        call(&mut e, "overlap", json!({ "a": ["Ex", "Ez"], "b": twice, "physics": "fld", "region": "silicon" }));
    close(num(&scaled["coefficient"]), 1.0, 1e-12);
    let (re, im) = pair(&scaled["overlap"]);
    assert!(re.abs() < 1e-9 * closed_in);
    close(im, 2.0 * closed_in, 1e-9);
    // two components that share nothing: Ex against Ez of the same position is not orthogonal, but the coefficient lies in [0, 1]
    let mixed = call(&mut e, "overlap", json!({ "a": "Ex", "b": "Ez", "physics": "fld" }));
    let c = num(&mixed["coefficient"]);
    assert!((0.0..=1.0 + 1e-12).contains(&c), "{c}");
}

#[test]
fn a_body_of_revolution_is_weighted_by_two_pi_r() {
    let mut e = engine(1e6);
    call(&mut e, "project_new", json!({ "title": "Ring", "h": 1e-7 }));
    call(&mut e, "add_rect", json!({ "name": "ring", "min": [1e-6, 0.0], "max": [3e-6, 1e-6] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "ring", "input": "ring" }] }));
    call(&mut e, "add_physics", json!({ "physics": "fld", "study": "poly", "order": 2 }));
    call(&mut e, "solve_run", json!({}));
    let pi = std::f64::consts::PI;
    let got = call(&mut e, "integrate", json!({ "expr": "one", "physics": "fld", "revolution": true }));
    close(num(&got["value"]), pi * (9e-12 - 1e-12) * 1e-6, 1e-9);
    close(num(&got["measure"]), pi * (9e-12 - 1e-12) * 1e-6, 1e-9);
    assert_eq!(got["revolution"], json!(true));
    let got = call(&mut e, "integrate", json!({ "expr": "p", "physics": "fld", "revolution": true }));
    let closed = exact(&[1e-6, 0.0], &[3e-6, 1e-6], &|x| 2.0 * pi * x[0] * polys(&u(x)).0);
    close(num(&got["value"]), closed, 1e-9);
    // the plane measure is the default and is told
    let plain = call(&mut e, "integrate", json!({ "expr": "one", "physics": "fld" }));
    close(num(&plain["value"]), 2e-6 * 1e-6, 1e-9);
    assert_eq!(plain["revolution"], json!(false));
}

#[test]
fn a_3d_mesh_integrates_probes_and_finds_extrema() {
    let mut e = engine(1.0);
    call(&mut e, "project_new", json!({ "dim": 3, "h": 0.12 }));
    call(&mut e, "add_box", json!({ "name": "b", "min": [0.0, 0.0, 0.0], "max": [1.0, 2.0, 0.5] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "b", "input": "b" }] }));
    call(&mut e, "add_physics", json!({ "physics": "fld", "study": "poly", "order": 2 }));
    call(&mut e, "solve_run", json!({}));
    let got = call(&mut e, "integrate", json!({ "expr": { "op": "abs2", "args": ["Ex"] }, "physics": "fld" }));
    let closed = exact(&[0.0, 0.0, 0.0], &[1.0, 2.0, 0.5], &|x| {
        let (_, ex, _) = polys(x);
        ex.0 * ex.0 + ex.1 * ex.1
    });
    close(num(&got["value"]), closed, 1e-9);
    let reference = raw(&e, None).integral(&["Ex"], &|a| (a[0].0 * a[0].0 + a[0].1 * a[0].1, 0.0));
    assert_eq!(num(&got["value"]).to_bits(), reference.0.to_bits());
    close(num(&got["measure"]), 1.0, 1e-9);
    let pt = [0.3, 1.1, 0.2];
    let probe = call(&mut e, "probe", json!({ "expr": "p", "physics": "fld", "points": [pt, [2.0, 0.0, 0.0]] }));
    assert_eq!(probe["inside"], json!([true, false]));
    close(num(&probe["values"][0]), polys(&pt).0, 1e-9);
    let top = call(&mut e, "extremum", json!({ "expr": "p", "physics": "fld" }));
    let pt: Vec<f64> = top["point"].as_array().unwrap().iter().map(num).collect();
    assert_eq!(pt.len(), 3);
    close(polys(&pt).0, num(&top["value"]), 1e-9);
}

#[test]
fn what_cannot_be_computed_is_said() {
    let mut e = engine(1e6);
    assert!(err(&mut e, "integrate", json!({ "expr": "p", "physics": "fld" })).contains("no result for 'fld'"));
    let mut e = strip();
    let msg = err(&mut e, "integrate", json!({ "expr": "nope", "physics": "fld" }));
    assert!(msg.contains("no field 'nope'") && msg.contains("Ex") && msg.contains("T"), "{msg}");
    assert!(err(&mut e, "extremum", json!({ "expr": "Ex", "physics": "fld" })).contains("complex"));
    assert!(
        err(&mut e, "integrate", json!({ "expr": "p", "physics": "fld", "region": "nowhere" })).contains("no region")
    );
    call(&mut e, "integrate", json!({ "expr": "p" })); // a lone result is the default physics
    assert!(err(&mut e, "integrate", json!({ "expr": { "op": "frobnicate", "args": [1] }, "physics": "fld" }))
        .contains("no operation"));
    assert!(err(&mut e, "integrate", json!({ "expr": { "op": "pow", "args": ["p"], "n": 99 }, "physics": "fld" }))
        .contains("1 to 16"));
    assert!(err(&mut e, "integrate", json!({ "expr": { "op": "div", "args": ["p", 0.0] }, "physics": "fld" }))
        .contains("not finite"));
    assert!(
        err(&mut e, "probe", json!({ "expr": "p", "physics": "fld", "points": [[0.0, 0.0, 0.0]] })).contains("2D mesh")
    );
    assert!(err(
        &mut e,
        "line_cut",
        json!({ "expr": "p", "physics": "fld", "from": [0.0, 0.0], "to": [1.0, 1.0], "n": 1 })
    )
    .contains("2 to"));
    assert!(err(&mut e, "overlap", json!({ "a": ["Ex", "Ez"], "b": ["Ex"], "physics": "fld" })).contains("must match"));
    // the vector name gives one handle per component
    let v = call(&mut e, "field", json!({ "physics": "fld", "name": "E" }));
    assert_eq!(v["components"]["x"]["handle"]["field"], json!("Ex"));
    assert_eq!(v["components"]["z"]["complex"], json!(true));
    assert!(v["components"].get("y").is_none());
    let f = call(&mut e, "field", json!({ "physics": "fld", "name": "T" }));
    assert_eq!(f["source"], json!("dof"));
    assert_eq!(f["unit"], json!("K"));
}

#[test]
fn a_result_is_read_on_its_own_mesh_and_says_when_it_is_stale() {
    let mut e = strip();
    let before = call(&mut e, "integrate", json!({ "expr": "p", "physics": "fld", "region": "oxide" }));
    assert_eq!(before["stale"], json!(false));
    // a different element size is a different mesh; the result still lives on the cells it was solved on
    call(&mut e, "set_mesh", json!({ "h": 2e-7 }));
    let after = call(&mut e, "integrate", json!({ "expr": "p", "physics": "fld", "region": "oxide" }));
    assert_eq!(after["stale"], json!(true));
    assert_eq!(after["value"], before["value"]);
    assert_eq!(after["cells"], before["cells"]);
    assert!(after["age_s"].as_f64().unwrap() >= 0.0);
}

// ------------------------------------------------------------------------------------------------------ real physics

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_str()
        .unwrap()
        .to_string()
}

#[test]
fn the_confinement_of_a_real_mode_is_the_ratio_the_plan_writes() {
    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
    call(&mut e, "solve_run", json!({}));
    let info = call(&mut e, "result_info", json!({ "physics": "ems" }));
    assert!(!info[0]["solutions"].as_array().unwrap().is_empty());
    // gamma = integrate(|Ex|^2 + |Ez|^2, silicon) / integrate(|Ex|^2 + |Ez|^2): one command or two, the same bits
    let power = json!({ "op": "add", "args": [{ "op": "abs2", "args": ["Ex"] }, { "op": "abs2", "args": ["Ez"] }] });
    let inside = call(&mut e, "integrate", json!({ "expr": power, "physics": "ems", "region": "Si" }));
    let total = call(&mut e, "integrate", json!({ "expr": power, "physics": "ems" }));
    let conf = call(&mut e, "confinement", json!({ "fields": ["Ex", "Ez"], "physics": "ems", "region": "Si" }));
    let gamma = num(&conf["fraction"]);
    assert_eq!(gamma.to_bits(), (num(&inside["value"]) / num(&total["value"])).to_bits());
    assert!(gamma > 0.05 && gamma < 1.0, "the mode is partly in the core: {gamma}");
    // the mode overlaps itself completely
    let own = call(&mut e, "overlap", json!({ "a": ["Ex", "Ez"], "b": ["Ex", "Ez"], "physics": "ems" }));
    close(num(&own["coefficient"]), 1.0, 1e-12);
    // the largest |E| of the mode is at a node of the mesh and is the maximum of the loop over the raw field
    let top = call(&mut e, "extremum", json!({ "expr": { "op": "abs", "args": ["Ex"] }, "physics": "ems" }));
    let (res, _) = e.result("ems").unwrap();
    let f = res.solutions[0].cell_fields.iter().find(|f| f.name == "Ex").unwrap();
    let want = (0..f.re.len())
        .map(|i| (f.re[i] * f.re[i] + if f.im.is_empty() { 0.0 } else { f.im[i] * f.im[i] }).sqrt())
        .fold(f64::MIN, f64::max);
    assert_eq!(num(&top["value"]).to_bits(), want.to_bits());
}

#[test]
fn a_flat_field_of_a_real_physics_reads_through_the_dofs() {
    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("04_thermal_bar") }));
    call(&mut e, "solve_run", json!({}));
    let info = call(&mut e, "result_info", json!({}));
    let physics = info[0]["physics"].as_str().unwrap().to_string();
    let name = info[0]["solutions"][0]["fields"][0]["name"].as_str().unwrap().to_string();
    let got = call(&mut e, "integrate", json!({ "expr": name, "physics": physics }));
    assert!(num(&got["value"]).is_finite() && num(&got["measure"]) > 0.0);
    // against the loop: the field's dof values gathered cell by cell
    let (res, _) = e.result(&physics).unwrap();
    let res = res.clone();
    let model = e.model().unwrap().clone();
    let mesh = model.solver_mesh();
    let order = res.order as usize;
    let vals = res.solutions[0].flat(&name).unwrap();
    let quad = quad_rule(mesh.dim, 4);
    let ns: Vec<Vec<f64>> = quad.lam.iter().map(|l| shape(l, order).0).collect();
    let vol = geometry(&mesh).vol;
    let (mut total, mut block) = (0.0, 0.0);
    let mut dofs = vec![];
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        let mut cell = 0.0;
        for q in 0..quad.w.len() {
            let mut r = 0.0;
            for a in 0..dofs.len() {
                r += ns[q][a] * vals[dofs[a]];
            }
            cell += quad.w[q] * r;
        }
        block += vol[c] * cell;
        if (c + 1) % 2048 == 0 {
            total += block;
            block = 0.0;
        }
    }
    total += block;
    assert_eq!(num(&got["value"]).to_bits(), total.to_bits());
}

#[test]
fn a_result_that_came_back_from_a_worker_reduces_to_the_same_bits() {
    let solve = |isolate: bool| {
        let mut e = Engine::new();
        call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
        if isolate {
            e.set_worker(PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker")), false);
        }
        let reply = call(&mut e, "solve_run", json!({ "isolate": isolate }));
        assert_eq!(reply["worker"].is_object(), isolate);
        e
    };
    let (mut here, mut there) = (solve(false), solve(true));
    let power = json!({ "op": "add", "args": [{ "op": "abs2", "args": ["Ex"] }, { "op": "abs2", "args": ["Ez"] }] });
    for region in [None, Some("Si")] {
        let mut args = json!({ "expr": power, "physics": "ems" });
        if let Some(r) = region {
            args["region"] = json!(r);
        }
        let (a, b) = (call(&mut here, "integrate", args.clone()), call(&mut there, "integrate", args));
        assert_eq!(a["value"], b["value"], "{region:?}");
        assert_eq!(a["measure"], b["measure"]);
    }
    let args = json!({ "fields": ["Ex", "Ez"], "physics": "ems", "region": "Si" });
    assert_eq!(
        call(&mut here, "confinement", args.clone())["fraction"],
        call(&mut there, "confinement", args)["fraction"]
    );
}
