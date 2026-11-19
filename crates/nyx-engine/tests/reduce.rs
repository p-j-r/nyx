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
        vec![
            PropDesc::new("rho", "Density", "kg/m3", PropKind::Scalar, &[], Some(1.0)),
            PropDesc::new("eps_r", "Relative permittivity", "", PropKind::Scalar, &[], Some(1.0)),
        ]
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
/// The same with the interpolated coordinates (x, z) of the point.
type IntegrandAt<'a> = &'a dyn Fn(&[(f64, f64)], [f64; 2]) -> (f64, f64);

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

    fn integral(&self, names: &[&str], f: Integrand) -> (f64, f64) {
        self.integral_at(names, &|a, _| f(a))
    }

    /// The wafer coordinates (x, z) of local node `a` of cell `c`: the vertices, then the midpoints of the vertex pairs.
    fn node(&self, c: usize, a: usize) -> [f64; 2] {
        let (sm, dim) = (&self.mesh, self.mesh.dim);
        let v = sm.cell(c);
        if a <= dim {
            let p = sm.point(v[a] as usize);
            return [p[0], p[1]];
        }
        let mut k = dim + 1;
        for i in 0..=dim {
            for j in i + 1..=dim {
                if k == a {
                    let (p, q) = (sm.point(v[i] as usize), sm.point(v[j] as usize));
                    return [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])];
                }
                k += 1;
            }
        }
        unreachable!()
    }

    /// `sum_c vol * sum_q w_q f_q` over blocks of 2048 cells, `f` given the interpolated (re, im) of each named field and the interpolated
    /// coordinates (x, z) at the point.
    fn integral_at(&self, names: &[&str], f: IntegrandAt) -> (f64, f64) {
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
                    let (mut x, mut z) = (0.0, 0.0);
                    for a in 0..self.nloc {
                        let p = self.node(c, a);
                        x += ns[q][a] * p[0];
                        z += ns[q][a] * p[1];
                    }
                    let v = f(&at, [x, z]);
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

// ------------------------------------------------------------------------------------------------------ N2: formulas

fn value(e: &mut Engine, expr: Value, region: Option<&str>) -> Value {
    let mut args = json!({ "expr": expr, "physics": "fld" });
    if let Some(r) = region {
        args["region"] = json!(r);
    }
    call(e, "integrate", args)
}

fn bits(v: &Value) -> (u64, u64) {
    if v.is_array() {
        let (re, im) = pair(v);
        (re.to_bits(), im.to_bits())
    } else {
        (num(v).to_bits(), 0)
    }
}

#[test]
fn a_formula_equals_the_handle_form_to_the_bit() {
    let mut e = strip();
    let h = |op: &str, args: Vec<Value>| json!({ "op": op, "args": args });
    let pow2 = |x: Value| json!({ "op": "pow", "args": [x], "n": 2 });
    // |Ex|^2 + |Ez|^2 written three ways
    let handle = h("add", vec![pow2(h("abs", vec![json!("Ex")])), pow2(h("abs", vec![json!("Ez")]))]);
    let by_handle = value(&mut e, handle, Some("silicon"));
    let by_text = value(&mut e, json!("abs(Ex)^2 + abs(Ez)^2"), Some("silicon"));
    assert_eq!(bits(&by_handle["value"]), bits(&by_text["value"]));
    let handle2 = h("add", vec![h("abs2", vec![json!("Ex")]), h("abs2", vec![json!("Ez")])]);
    let by_handle2 = value(&mut e, handle2, Some("silicon"));
    let by_text2 = value(&mut e, json!("abs2(Ex) + abs2(Ez)"), Some("silicon"));
    assert_eq!(bits(&by_handle2["value"]), bits(&by_text2["value"]));
    // a formula and a handle mix inside one operation
    let mixed = value(&mut e, h("add", vec![json!("abs2(Ex)"), h("abs2", vec![json!("Ez")])]), Some("silicon"));
    assert_eq!(bits(&mixed["value"]), bits(&by_handle2["value"]));
    // and the loop on the raw arrays, and the closed form (degree 4: exact)
    let reference = raw(&e, Some("silicon"))
        .integral(&["Ex", "Ez"], &|a| (a[0].0 * a[0].0 + a[0].1 * a[0].1 + (a[1].0 * a[1].0 + a[1].1 * a[1].1), 0.0));
    assert_eq!(num(&by_text2["value"]).to_bits(), reference.0.to_bits());
    let closed = exact(&SILICON.0, &SILICON.1, &|x| {
        let (_, ex, ez) = polys(&u(x));
        ex.0 * ex.0 + ex.1 * ex.1 + ez.0 * ez.0 + ez.1 * ez.1
    });
    close(num(&by_text["value"]), closed, 1e-9);
    // a complex product: a pair, equal to the handle `mul` to the bit
    let prod = value(&mut e, json!("Ex * Ez"), None);
    assert_eq!(prod["complex"], json!(true));
    let by_handle = value(&mut e, h("mul", vec![json!("Ex"), json!("Ez")]), None);
    assert_eq!(bits(&prod["value"]), bits(&by_handle["value"]));
    // a plain field name still reads the field, and a parenthesised one is the same formula
    let plain = value(&mut e, json!("p"), None);
    let as_formula = value(&mut e, json!("(p)"), None);
    assert_eq!(bits(&plain["value"]), bits(&as_formula["value"]));
    let q1 = value(&mut e, json!("p / 4"), None);
    let q2 = value(&mut e, json!("p * 0.25"), None);
    close(num(&q1["value"]), num(&q2["value"]), 1e-13);
    // the prefixes name the physics and the solution
    let a = value(&mut e, json!("`fld:Ex` * `fld:0:Ez`"), None);
    assert_eq!(bits(&a["value"]), bits(&prod["value"]));
    assert!(err(&mut e, "integrate", json!({ "expr": "`fld:5:Ex`", "physics": "fld" })).contains("there is no 5"));
}

#[test]
fn coordinates_where_and_regions_read_each_point() {
    let mut e = strip();
    // x * p has degree 3, exact with the rule: the closed form and the loop (with interpolated coordinates) agree
    let xp = value(&mut e, json!("x * p"), None);
    close(num(&xp["value"]), exact(&WINDOW.0, &WINDOW.1, &|x| x[0] * polys(&u(x)).0), 1e-9);
    let reference = raw(&e, None).integral_at(&["p"], &|a, xz| (xz[0] * a[0].0, 0.0));
    assert_eq!(num(&xp["value"]).to_bits(), reference.0.to_bits());
    // z = 0 is the edge between oxide and air, so each cell is on one side of it and where(z > 0, ..) is the integral over the upper half
    let upper = value(&mut e, json!("where(z > 0, p, 0)"), None);
    let half = ([-2e-6, 0.0], [2e-6, 1.5e-6]);
    close(num(&upper["value"]), exact(&half.0, &half.1, &|x| polys(&u(x)).0), 1e-9);
    let reference = raw(&e, None).integral_at(&["p"], &|a, xz| (if xz[1] > 0.0 { a[0].0 } else { 0.0 }, 0.0));
    assert_eq!(num(&upper["value"]).to_bits(), reference.0.to_bits());
    // the region indicator is the area of the region; `y` is 0 on a cross-section
    let area = value(&mut e, json!("region(\"silicon\")"), None);
    close(num(&area["value"]), 0.5e-6 * 0.22e-6, 1e-9);
    let ox = value(&mut e, json!("where(region(\"oxide\"), 1, 0)"), None);
    close(num(&ox["value"]), 4e-6 * 1.5e-6, 1e-9);
    let y = value(&mut e, json!("y + 1"), None);
    close(num(&y["value"]), 4e-6 * 3e-6, 1e-9);
    // probes, extrema and line cuts read the same formulas
    let pr = call(&mut e, "probe", json!({ "expr": "x * z", "physics": "fld", "points": [[1e-7, 1e-7], [3e-6, 0.0]] }));
    close(num(&pr["values"][0]), 1e-14, 1e-9);
    assert_eq!(pr["inside"], json!([true, false]));
    let m = call(&mut e, "extremum", json!({ "expr": "where(z > 0, p, -1e9)", "physics": "fld" }));
    assert!(num(&m["value"]) > -1e9);
    let lc = call(
        &mut e,
        "line_cut",
        json!({ "expr": "x", "physics": "fld", "from": [-1e-6, 0.5e-6], "to": [1e-6, 0.5e-6], "n": 3 }),
    );
    close(num(&lc["values"][2]), 1e-6, 1e-9);
}

#[test]
fn vectors_are_the_components_of_a_name() {
    let mut e = strip();
    // norm2(E) is abs2(Ex) + abs2(Ez), the same tree, so the same bits
    let n2 = value(&mut e, json!("norm2(E)"), Some("silicon"));
    let by_hand = value(&mut e, json!("abs2(Ex) + abs2(Ez)"), Some("silicon"));
    assert_eq!(bits(&n2["value"]), bits(&by_hand["value"]));
    assert_eq!(n2["complex"], json!(false));
    // inner(E, E) is the same sum written with conjugates: complex-typed, the same real part to rounding and a zero imaginary part
    let inner = value(&mut e, json!("inner(E, E)"), Some("silicon"));
    let (re, im) = pair(&inner["value"]);
    close(re, num(&n2["value"]), 1e-12);
    assert_eq!(im, 0.0);
    // dot(E, E) has no conjugate: the sum of the squares, equal to the loop to the bit
    let dot = value(&mut e, json!("dot(E, E)"), Some("silicon"));
    let reference = raw(&e, Some("silicon")).integral(&["Ex", "Ez"], &|a| {
        (
            (a[0].0 * a[0].0 - a[0].1 * a[0].1) + (a[1].0 * a[1].0 - a[1].1 * a[1].1),
            (a[0].0 * a[0].1 + a[0].1 * a[0].0) + (a[1].0 * a[1].1 + a[1].1 * a[1].0),
        )
    });
    assert_eq!(bits(&dot["value"]), (reference.0.to_bits(), reference.1.to_bits()));
    let comp = value(&mut e, json!("E.x * conj(E.x)"), Some("silicon"));
    let ex2 = value(&mut e, json!("abs2(Ex)"), Some("silicon"));
    close(pair(&comp["value"]).0, num(&ex2["value"]), 1e-12);
    // a vector can be scaled before it is reduced
    let scaled = value(&mut e, json!("norm2(2 * E)"), Some("silicon"));
    close(num(&scaled["value"]), 4.0 * num(&n2["value"]), 1e-12);
    let norm = value(&mut e, json!("where(norm(E) > 0, 1, 0)"), Some("silicon"));
    assert!(num(&norm["value"]) > 0.0);
}

#[test]
fn material_values_are_those_the_solver_read() {
    let mut e = strip();
    call(&mut e, "set_material", json!({ "physics": "fld", "region": "oxide", "material": { "rho": 2200.0 } }));
    call(&mut e, "set_material", json!({ "physics": "fld", "region": "silicon", "material": { "rho": 2330.0 } }));
    call(&mut e, "solve_run", json!({}));
    let (a_ox, a_si) = (4e-6 * 1.5e-6, 0.5e-6 * 0.22e-6);
    let a_air = 4e-6 * 3e-6 - a_ox - a_si;
    let rho = value(&mut e, json!("rho"), None);
    close(num(&rho["value"]), 2200.0 * a_ox + 2330.0 * a_si + 1.0 * a_air, 1e-12);
    let si = value(&mut e, json!("material(\"rho\") * region(\"silicon\")"), None);
    close(num(&si["value"]), 2330.0 * a_si, 1e-12);
    // a material value with a field: rho * p over the oxide is 2200 times the integral of p
    let rp = value(&mut e, json!("rho * p"), Some("oxide"));
    let p = value(&mut e, json!("p"), Some("oxide"));
    close(num(&rp["value"]), 2200.0 * num(&p["value"]), 1e-12);
    // `eps` is `eps_r`; every region has the default 1
    let eps = value(&mut e, json!("eps"), None);
    close(num(&eps["value"]), 4e-6 * 3e-6, 1e-12);
    // not a property of this physics
    let msg = err(&mut e, "integrate", json!({ "expr": "material(\"E\")", "physics": "fld" }));
    assert!(msg.contains("no material property 'E'") && msg.contains("rho"), "{msg}");
    // a changed material makes the result stale, and the answer says so
    call(&mut e, "set_material", json!({ "physics": "fld", "region": "oxide", "material": { "rho": 1000.0 } }));
    let stale = value(&mut e, json!("rho"), Some("oxide"));
    assert_eq!(stale["stale"], json!(true));
}

#[test]
fn a_formula_that_cannot_be_read_says_where_and_why() {
    let mut e = strip();
    let msg = err(&mut e, "integrate", json!({ "expr": "abs(Ex) * foo", "physics": "fld" }));
    assert!(
        msg.contains("unknown name 'foo'") && msg.contains("Ex") && msg.contains("rho") && msg.contains("character 14"),
        "{msg}"
    );
    // an order of a complex value is refused before any point is read, naming the formula
    let msg = err(&mut e, "integrate", json!({ "expr": "where(Ex > 0, 1, 2)", "physics": "fld" }));
    assert!(msg.contains("where(Ex > 0, 1, 2)") && msg.contains("complex") && msg.contains("abs"), "{msg}");
    assert!(err(&mut e, "integrate", json!({ "expr": "log(Ex)", "physics": "fld" })).contains("complex"));
    // a vector is not a number
    assert!(err(&mut e, "integrate", json!({ "expr": "E", "physics": "fld" })).contains("vector"));
    // syntax
    assert!(err(&mut e, "integrate", json!({ "expr": "abs(Ex", "physics": "fld" })).contains("expected"));
    assert!(err(&mut e, "integrate", json!({ "expr": "__import__('os')", "physics": "fld" })).contains("does not read"));
    // a division by zero that the condition guards is fine, and an unguarded one is the same error as before
    let ok = value(&mut e, json!("where(z != 0, 1 / z, 0)"), None);
    assert!(num(&ok["value"]).is_finite());
    assert!(err(&mut e, "integrate", json!({ "expr": "1 / (z - z)", "physics": "fld" })).contains("not finite"));
    // a name that is no field and no formula keeps the message that lists the fields
    let msg = err(&mut e, "integrate", json!({ "expr": "nope", "physics": "fld" }));
    assert!(msg.contains("no field 'nope'") && msg.contains("Ex"), "{msg}");
}

#[test]
fn the_permittivity_of_a_real_mode_comes_from_its_materials() {
    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
    call(&mut e, "solve_run", json!({}));
    // eps_r in a region is the index squared of the material at the study wavelength: a constant, so its integral is it times the area
    let area = call(&mut e, "integrate", json!({ "expr": 1.0, "physics": "ems", "region": "Si" }));
    let eps = call(&mut e, "integrate", json!({ "expr": "eps_r", "physics": "ems", "region": "Si" }));
    let n_si = (num(&eps["value"]) / num(&area["value"])).sqrt();
    assert!((3.4..3.6).contains(&n_si), "silicon index at the wavelength of the example: {n_si}");
    // the energy fraction of the plan: integrate(eps |E|^2, core) / integrate(eps |E|^2), as one formula and as the handle form
    let by_text = call(&mut e, "integrate", json!({ "expr": "eps_r * norm2(E)", "physics": "ems", "region": "Si" }));
    let eps_handle = json!("eps_r");
    let power = json!({ "op": "add", "args": [{ "op": "abs2", "args": ["Ex"] }, { "op": "abs2", "args": ["Ez"] }] });
    let by_hand = call(
        &mut e,
        "integrate",
        json!({ "expr": { "op": "mul", "args": [eps_handle, power] }, "physics": "ems", "region": "Si" }),
    );
    // the vector E has Ey too, so norm2(E) adds |Ey|^2 and is larger
    assert!(num(&by_text["value"]) >= num(&by_hand["value"]));
    let ey = call(&mut e, "integrate", json!({ "expr": "eps_r * abs2(Ey)", "physics": "ems", "region": "Si" }));
    close(num(&by_text["value"]), num(&by_hand["value"]) + num(&ey["value"]), 1e-12);
    let gamma = num(&by_text["value"])
        / num(&call(&mut e, "integrate", json!({ "expr": "eps_r * norm2(E)", "physics": "ems" }))["value"]);
    assert!(gamma > 0.05 && gamma < 1.0, "{gamma}");
    // a material value that is a tensor is refused, not read as its default
    let msg = err(&mut e, "integrate", json!({ "expr": "material(\"rho\")", "physics": "ems" }));
    assert!(msg.contains("no material property 'rho'") && msg.contains("eps_r"), "{msg}");
}

#[test]
fn the_entries_of_a_rotated_crystal_are_read_in_wafer_axes() {
    // r-plane sapphire under a silicon strip (example 13): the c axis is tilted 57.6 degrees in the y-z plane, so the permittivity has a
    // non-zero yz entry. A swapped index or a frame change would show in that entry, which a diagonal crystal could not.
    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("13_silicon_on_sapphire") }));
    call(&mut e, "solve_run", json!({ "physics": "ems" }));
    let area =
        num(&call(&mut e, "integrate", json!({ "expr": 1.0, "physics": "ems", "region": "sapphire*" }))["value"]);
    let t = 57.6f64.to_radians();
    let euler =
        femcore::material::euler_from_rotation(&[[1.0, 0.0, 0.0], [0.0, t.cos(), -t.sin()], [0.0, t.sin(), t.cos()]]);
    let lib = femcore::material::Library::load_dir(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../materials"))
        .unwrap();
    let lambda = 299792458.0 / (299792458.0 / 1550e-9);
    let want = lib.get("sapphire").unwrap().optical_permittivity(lambda, None, euler).unwrap();
    let axes = ["x", "y", "z"];
    for i in 0..3 {
        for j in 0..3 {
            let got = call(
                &mut e,
                "integrate",
                json!({ "expr": format!("eps_r_{}{}", axes[i], axes[j]), "physics": "ems", "region": "sapphire*" }),
            );
            close(num(&got["value"]), want[i][j].re * area, 1e-9);
        }
    }
    assert!(want[1][2].re.abs() > 1e-3, "the crystal is tilted, so yz is not zero: {}", want[1][2].re);
    // eps_r_zy is the same entry as eps_r_yz, and `eps_xy` is `eps_r_xy`
    let yz = call(&mut e, "integrate", json!({ "expr": "eps_r_yz", "physics": "ems", "region": "sapphire*" }));
    let zy = call(&mut e, "integrate", json!({ "expr": "eps_r_zy", "physics": "ems", "region": "sapphire*" }));
    assert_eq!(yz["value"], zy["value"]);
    let a = call(&mut e, "integrate", json!({ "expr": "eps_xy", "physics": "ems", "region": "sapphire*" }));
    let b = call(&mut e, "integrate", json!({ "expr": "eps_r_xy", "physics": "ems", "region": "sapphire*" }));
    assert_eq!(a["value"], b["value"]);
    // an isotropic region of the same mesh gives v I: silicon has the same value on the diagonal and zero off it
    let sxx = call(&mut e, "integrate", json!({ "expr": "eps_r_xx", "physics": "ems", "region": "silicon" }));
    let szz = call(&mut e, "integrate", json!({ "expr": "eps_r_zz", "physics": "ems", "region": "silicon" }));
    let sxy = call(&mut e, "integrate", json!({ "expr": "eps_r_xy", "physics": "ems", "region": "silicon" }));
    assert_eq!(sxx["value"], szz["value"]);
    assert_eq!(num(&sxy["value"]), 0.0);
    assert!(num(&sxx["value"]) > 0.0);
    // the plain name is refused on this mesh, and the message names the entries; so does a formula that combines the entries with a field
    let msg = err(&mut e, "integrate", json!({ "expr": "eps_r", "physics": "ems" }));
    assert!(msg.contains("eps_r_xx") && msg.contains("tensor"), "{msg}");
    let msg = err(&mut e, "integrate", json!({ "expr": "eps_r_xq * 1", "physics": "ems" }));
    assert!(msg.contains("eps_r_xx"), "the names of the tensor are listed: {msg}");
    // a field with a tensor entry: the permittivity-weighted intensity along the propagation, as a formula
    let w =
        call(&mut e, "integrate", json!({ "expr": "eps_r_zz * abs2(Ey)", "physics": "ems", "region": "sapphire*" }));
    assert!(num(&w["value"]).is_finite());
}

#[test]
fn the_wavelength_argument_moves_the_material_values_and_nothing_else() {
    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
    call(&mut e, "solve_run", json!({}));
    let lib = femcore::material::Library::load_dir(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../materials"))
        .unwrap();
    let area = num(&call(&mut e, "integrate", json!({ "expr": 1.0, "physics": "ems", "region": "Si" }))["value"]);
    let eps_at = |e: &mut Engine, w: Option<f64>| {
        let mut args = json!({ "expr": "eps_r", "physics": "ems", "region": "Si" });
        if let Some(w) = w {
            args["wavelength"] = json!(w);
        }
        num(&call(e, "integrate", args)["value"]) / area
    };
    // the study's own wavelength gives the same as no argument, bit for bit
    let study = eps_at(&mut e, None);
    let own = eps_at(&mut e, Some(299792458.0 / 193414489032258.06));
    assert_eq!(study.to_bits(), own.to_bits());
    // another wavelength reads the model of the material there: the tabulated silicon of the library
    for lambda in [1.31e-6, 1.7e-6] {
        let want = lib.get("silicon").unwrap().optical_permittivity(lambda, None, [0.0; 3]).unwrap()[0][0].re;
        close(eps_at(&mut e, Some(lambda)), want, 1e-9);
    }
    assert!((eps_at(&mut e, Some(1.31e-6)) - study).abs() > 1e-4, "silicon is dispersive here");
    // the fields are those of the solve whatever the wavelength: a formula without a material value does not move
    let a = call(&mut e, "integrate", json!({ "expr": "abs2(Ex)", "physics": "ems", "region": "Si" }));
    let b = call(
        &mut e,
        "integrate",
        json!({ "expr": "abs2(Ex)", "physics": "ems", "region": "Si", "wavelength": 1.31e-6 }),
    );
    assert_eq!(a["value"], b["value"]);
    // outside the range of the table the material says so; a bad value is refused
    assert!(
        err(&mut e, "integrate", json!({ "expr": "eps_r", "physics": "ems", "wavelength": 5e-7 })).contains("range")
    );
    assert!(err(&mut e, "integrate", json!({ "expr": "eps_r", "physics": "ems", "wavelength": -1.0 }))
        .contains("wavelength"));
}

#[test]
fn a_weight_makes_the_confinement_the_energy_fraction() {
    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
    call(&mut e, "solve_run", json!({}));
    let power = json!({ "op": "add", "args": [{ "op": "add", "args": [{ "op": "abs2", "args": ["Ex"] }, { "op": "abs2", "args": ["Ey"] }] }, { "op": "abs2", "args": ["Ez"] }] });
    let ratio = |e: &mut Engine, expr: Value| -> f64 {
        let inside = num(&call(e, "integrate", json!({ "expr": expr, "physics": "ems", "region": "Si" }))["value"]);
        let total = num(&call(e, "integrate", json!({ "expr": expr, "physics": "ems" }))["value"]);
        inside / total
    };
    let conf = |e: &mut Engine, weight: Option<&str>| -> Value {
        let mut a = json!({ "fields": ["Ex", "Ey", "Ez"], "physics": "ems", "region": "Si" });
        if let Some(w) = weight {
            a["weight"] = json!(w);
        }
        call(e, "confinement", a)
    };
    // no weight: the ratio of the integrals of the sum of squared moduli, to the bit, as before
    let plain = conf(&mut e, None);
    assert_eq!(num(&plain["fraction"]).to_bits(), ratio(&mut e, power.clone()).to_bits());
    assert_eq!(plain["weighted"], json!(false));
    // a weight: the ratio of the integrals of weight * power, the energy fraction of the plan, to the bit
    let weighted = conf(&mut e, Some("eps_r"));
    let by_hand = ratio(&mut e, json!({ "op": "mul", "args": ["eps_r", power] }));
    assert_eq!(num(&weighted["fraction"]).to_bits(), by_hand.to_bits());
    assert_eq!(weighted["weighted"], json!(true));
    // silicon is the dense material, so weighting by the permittivity raises the fraction of the core
    assert!(num(&weighted["fraction"]) > num(&plain["fraction"]));
    assert!(num(&weighted["fraction"]) < 1.0);
    // a weight with a region: 1 inside the region is the same as no weight
    let one = conf(&mut e, Some("1"));
    assert_eq!(num(&one["fraction"]).to_bits(), num(&plain["fraction"]).to_bits());
    // a complex weight is refused (the stand-in physics has a complex Ex; the real mode above has none)
    let mut s = strip();
    let msg =
        err(&mut s, "confinement", json!({ "fields": ["Ex"], "physics": "fld", "region": "silicon", "weight": "Ex" }));
    assert!(msg.contains("weight") && msg.contains("real"), "{msg}");
}

const OUTER_2D: [&str; 6] = ["window.top", "window.left", "window.right", "box.left", "box.right", "box.bottom"];

#[test]
fn a_line_integral_follows_the_boundary_with_its_normal() {
    let mut e = strip();
    let line = |e: &mut Engine, expr: Value, boundary: &str, side: Option<&str>| -> Value {
        let mut a = json!({ "expr": expr, "physics": "fld", "boundary": boundary });
        if let Some(s) = side {
            a["side"] = json!(s);
        }
        call(e, "integrate", a)
    };
    // the length of known boundaries (expr 1): an outer edge, an interface read from one side, the edges of the core
    let top = line(&mut e, json!(1.0), "window.top", None);
    close(num(&top["value"]), 4e-6, 1e-12);
    close(num(&top["measure"]), 4e-6, 1e-12);
    assert!(top["facets"].as_u64().unwrap() >= 60);
    close(num(&line(&mut e, json!(1.0), "oxide|silicon", Some("silicon"))["value"]), 0.5e-6, 1e-12);
    close(num(&line(&mut e, json!(1.0), "air|silicon", Some("silicon"))["value"]), 0.5e-6 + 2.0 * 0.22e-6, 1e-12);
    close(num(&line(&mut e, json!(1.0), "air|oxide", Some("oxide"))["value"]), 4e-6 - 0.5e-6, 1e-12);
    // the normal points out of the cell each facet is read in: up on top, and either way at an interface by the side
    close(num(&line(&mut e, json!("nz"), "window.top", None)["value"]), 4e-6, 1e-12);
    close(num(&line(&mut e, json!("nx"), "window.left", None)["value"]), -1.5e-6, 1e-12);
    close(num(&line(&mut e, json!("nz"), "air|oxide", Some("oxide"))["value"]), 3.5e-6, 1e-12);
    close(num(&line(&mut e, json!("nz"), "air|oxide", Some("air"))["value"]), -3.5e-6, 1e-12);
    assert!(
        num(&line(&mut e, json!("ny"), "window.top", None)["value"]).abs() < 1e-18,
        "a cross-section has no y normal"
    );
    // a closed form: the integral of the quadratic p along the top edge is exact in the rule
    let want = exact(&[-2e-6], &[2e-6], &|x| polys(&[x[0] * 1e6, 1.5]).0);
    close(num(&line(&mut e, json!("p"), "window.top", None)["value"]), want, 1e-9);
    // the divergence theorem: the flux of F = (Im Ex, p) through the closed outer boundary equals the area integral of div F
    let flux: f64 = OUTER_2D.iter().map(|b| num(&line(&mut e, json!("im(Ex) * nx + p * nz"), b, None)["value"])).sum();
    let (lo, hi) = WINDOW;
    let div = exact(&lo, &hi, &|x| {
        let a = x[0] * 1e6;
        let c = x[1] * 1e6;
        1e6 * (2.0 * a + (-3.0 + a - 2.0 * c))
    });
    close(flux, div, 1e-9);
    // the same by the engine's own area integral of the divergence, written as a formula of the coordinates
    let area = call(
        &mut e,
        "integrate",
        json!({ "expr": "1e6 * (2 * x * 1e6 + (-3 + x * 1e6 - 2 * z * 1e6))", "physics": "fld" }),
    );
    close(num(&area["value"]), div, 1e-9);
    // a field and a material differ across an interface: each side reads its own
    call(&mut e, "set_material", json!({ "physics": "fld", "region": "oxide", "material": { "rho": 2200.0 } }));
    call(&mut e, "set_material", json!({ "physics": "fld", "region": "silicon", "material": { "rho": 2330.0 } }));
    call(&mut e, "solve_run", json!({}));
    close(num(&line(&mut e, json!("rho"), "oxide|silicon", Some("silicon"))["value"]), 2330.0 * 0.5e-6, 1e-12);
    close(num(&line(&mut e, json!("rho"), "oxide|silicon", Some("oxide"))["value"]), 2200.0 * 0.5e-6, 1e-12);
    // an interface without a side is an error that names both regions; so is the wrong side, and a region next to a boundary
    let msg = err(&mut e, "integrate", json!({ "expr": 1.0, "physics": "fld", "boundary": "oxide|silicon" }));
    assert!(
        msg.contains("interface") && msg.contains("oxide") && msg.contains("silicon") && msg.contains("side"),
        "{msg}"
    );
    let msg =
        err(&mut e, "integrate", json!({ "expr": 1.0, "physics": "fld", "boundary": "oxide|silicon", "side": "air" }));
    assert!(msg.contains("side"), "{msg}");
    let msg =
        err(&mut e, "integrate", json!({ "expr": 1.0, "physics": "fld", "boundary": "window.top", "region": "air" }));
    assert!(msg.contains("not both"), "{msg}");
    assert!(err(&mut e, "integrate", json!({ "expr": 1.0, "physics": "fld", "side": "air" })).contains("boundary"));
    assert!(
        err(&mut e, "integrate", json!({ "expr": 1.0, "physics": "fld", "boundary": "nowhere" })).contains("nowhere")
    );
    // the normal is a boundary's: elsewhere it is refused, in a plain integral and in a probe
    let msg = err(&mut e, "integrate", json!({ "expr": "nx", "physics": "fld" }));
    assert!(msg.contains("normal") && msg.contains("boundary"), "{msg}");
    assert!(
        err(&mut e, "probe", json!({ "expr": "nz", "physics": "fld", "points": [[0.0, 0.0]] })).contains("boundary")
    );
    // the boundary of the world has one cell per facet: the answer does not depend on a side being given for it
    let some = line(&mut e, json!(1.0), "window.top", Some("air"));
    assert_eq!(some["value"], top["value"]);
}

#[test]
fn a_surface_integral_in_3d_and_a_revolution_weight() {
    let mut e = engine(1.0);
    call(&mut e, "project_new", json!({ "dim": 3, "h": 0.12 }));
    call(&mut e, "add_box", json!({ "name": "b", "min": [0.0, 0.0, 0.0], "max": [1.0, 2.0, 0.5] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "b", "input": "b" }] }));
    call(&mut e, "add_physics", json!({ "physics": "fld", "study": "poly", "order": 2 }));
    call(&mut e, "solve_run", json!({}));
    let face = |e: &mut Engine, expr: &str, b: &str| {
        num(&call(e, "integrate", json!({ "expr": expr, "physics": "fld", "boundary": b }))["value"])
    };
    // areas and signed normals of the six faces
    close(face(&mut e, "1", "b.zmax"), 2.0, 1e-12);
    close(face(&mut e, "nz", "b.zmax"), 2.0, 1e-12);
    close(face(&mut e, "nz", "b.zmin"), -2.0, 1e-12);
    close(face(&mut e, "nx", "b.xmin"), -1.0, 1e-12);
    close(face(&mut e, "ny", "b.ymax"), 0.5, 1e-12);
    assert!(face(&mut e, "nx", "b.zmax").abs() < 1e-12, "the top face has no x normal");
    // a closed form: p over the top face (z = 0.5)
    close(face(&mut e, "p", "b.zmax"), exact(&[0.0, 0.0], &[1.0, 2.0], &|x| polys(&[x[0], x[1], 0.5]).0), 1e-9);
    // the divergence theorem for F = (p, 2p, p): the flux through the six faces equals the integral of div F = -0.5 + 2a - c over the box
    let flux: f64 = ["b.xmax", "b.xmin", "b.ymax", "b.ymin", "b.zmax", "b.zmin"]
        .iter()
        .map(|b| face(&mut e, "p * (nx + 2 * ny + nz)", b))
        .sum();
    let div = exact(&[0.0, 0.0, 0.0], &[1.0, 2.0, 0.5], &|x| -0.5 + 2.0 * x[0] - x[2]);
    close(flux, div, 1e-9);
    // `revolution` is for a (r, z) mesh
    assert!(err(
        &mut e,
        "integrate",
        json!({ "expr": 1.0, "physics": "fld", "boundary": "b.zmax", "revolution": true })
    )
    .contains("2D"));
}

#[test]
fn a_boundary_integral_of_revolution_weights_by_two_pi_r() {
    let mut e = strip();
    // read as an (r, z) mesh: the edge at x = 2e-6 above the oxide, weighted by 2 pi r (z from 0 to 1.5e-6)
    let got = call(
        &mut e,
        "integrate",
        json!({ "expr": 1.0, "physics": "fld", "boundary": "window.right", "revolution": true }),
    );
    let want = 2.0 * std::f64::consts::PI * 2e-6 * 1.5e-6;
    close(num(&got["value"]), want, 1e-12);
    close(num(&got["measure"]), want, 1e-12);
    assert_eq!(got["revolution"], json!(true));
}

#[test]
fn an_interface_in_3d_needs_a_side_and_has_the_area_of_the_contact() {
    let mut e = engine(1.0);
    call(&mut e, "project_new", json!({ "dim": 3, "h": 0.12 }));
    call(&mut e, "add_box", json!({ "name": "lo", "min": [0.0, 0.0, 0.0], "max": [1.0, 2.0, 0.25] }));
    call(&mut e, "add_box", json!({ "name": "hi", "min": [0.0, 0.0, 0.25], "max": [1.0, 2.0, 0.5] }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "lo", "input": "lo" }, { "name": "hi", "input": "hi" }] }),
    );
    call(&mut e, "add_physics", json!({ "physics": "fld", "study": "poly", "order": 2 }));
    call(&mut e, "solve_run", json!({}));
    let names = call(&mut e, "geometry_names", json!({}));
    let name = ["lo|hi", "hi|lo"]
        .into_iter()
        .find(|n| names["boundaries"].as_array().unwrap().iter().any(|b| b == n))
        .unwrap_or_else(|| panic!("the contact has a name: {}", names["boundaries"]));
    let at = |e: &mut Engine, expr: &str, side: Option<&str>| -> Result<Value, String> {
        let mut a = json!({ "expr": expr, "physics": "fld", "boundary": name });
        if let Some(s) = side {
            a["side"] = json!(s);
        }
        e.call("integrate", &a).map_err(|x| x.to_string())
    };
    // two cells share each face: without a side it is an error that says so
    let msg = at(&mut e, "1", None).unwrap_err();
    assert!(msg.contains("interface") && msg.contains("side"), "{msg}");
    // the contact is the square 1 x 2 at z = 0.25; the normal points out of the side read: up out of `lo`, down out of `hi`
    close(num(&at(&mut e, "1", Some("lo")).unwrap()["value"]), 2.0, 1e-12);
    close(num(&at(&mut e, "nz", Some("lo")).unwrap()["value"]), 2.0, 1e-12);
    close(num(&at(&mut e, "nz", Some("hi")).unwrap()["value"]), -2.0, 1e-12);
    close(
        num(&at(&mut e, "p", Some("lo")).unwrap()["value"]),
        exact(&[0.0, 0.0], &[1.0, 2.0], &|x| polys(&[x[0], x[1], 0.25]).0),
        1e-9,
    );
}
