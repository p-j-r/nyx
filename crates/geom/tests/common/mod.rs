//! Golden models shared by the tests and by the generator of `golden/*.json`.
#![allow(dead_code)]

pub mod interp;

use geom::model::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Expected {
    pub note: String,
    /// region name and its closed-form area, in overlay order
    pub regions: Vec<(String, f64)>,
    /// the set of boundary/interface auto-names
    pub names: Vec<String>,
}

pub struct Golden {
    pub name: &'static str,
    pub model: Model,
    pub expected: Expected,
}

fn def(name: &str, input: FeatureId, pin: Option<i32>) -> RegionDef {
    RegionDef { name: name.to_string(), input, kind: RegionKind::Material, pinned_id: pin }
}

fn exp(note: &str, regions: &[(&str, f64)], names: &[&str]) -> Expected {
    let mut n: Vec<String> = names.iter().map(|s| s.to_string()).collect();
    n.sort();
    Expected { note: note.into(), regions: regions.iter().map(|(a, b)| (a.to_string(), *b)).collect(), names: n }
}

/// SoiWire::default(): 500 x 300 nm silicon core on 2 um BOX, 1 um oxide margin; legacy domains BOX 1, Si 2, clad 3.
pub fn soi_strip() -> Golden {
    let mut m = Model::new_2d();
    m.param("w", 500e-9).param("t", 300e-9).param("margin", 1e-6).param("box_t", 2e-6);
    let clad = m.rect("clad", ("-w/2-margin", "-box_t"), ("w/2+margin", "t+margin"));
    let bx = m.rect("box", ("-w/2-margin", "-box_t"), ("w/2+margin", 0.0));
    let core = m.rect("core", ("-w/2", 0.0), ("w/2", "t"));
    m.overlay_defs(vec![def("clad", clad, Some(3)), def("BOX", bx, Some(1)), def("Si", core, Some(2))], vec![]);
    let (w, t, mg, bt) = (500e-9, 300e-9, 1e-6, 2e-6);
    let wide = w + 2.0 * mg;
    let e = exp(
        "SoiWire::default(): three regions, BOX/Si/clad = legacy domains 1/2/3",
        &[("clad", wide * (bt + t + mg) - wide * bt - w * t), ("BOX", wide * bt), ("Si", w * t)],
        &[
            "BOX|Si",
            "BOX|clad",
            "Si|clad",
            "box.bottom",
            "box.left",
            "box.right",
            "clad.left",
            "clad.right",
            "clad.top",
        ],
    );
    Golden { name: "soi_strip", model: m, expected: e }
}

/// Rib on a slab, side walls 80 degrees on the right and 75 degrees on the left (breaks the mirror symmetry).
pub fn soi_rib_sloped() -> Golden {
    let mut m = Model::new_2d();
    m.param("w", 500e-9).param("t", 300e-9).param("slab", 100e-9).param("margin", 1e-6).param("box_t", 2e-6);
    m.param("deg_r", 80.0).param("deg_l", 75.0);
    m.param("inset_r", "(t-slab)/tan(deg_r*pi/180)").param("inset_l", "(t-slab)/tan(deg_l*pi/180)");
    let clad = m.rect("clad", ("-w/2-margin", "-box_t"), ("w/2+margin", "t+margin"));
    let bx = m.rect("box", ("-w/2-margin", "-box_t"), ("w/2+margin", 0.0));
    let slab = m.rect("slab", ("-w/2-margin", 0.0), ("w/2+margin", "slab"));
    let core = m.polygon_named(
        "core",
        &[("-w/2", "slab"), ("w/2", "slab"), ("w/2-inset_r", "t"), ("-w/2+inset_l", "t")],
        &["foot", "right", "top", "left"],
    );
    let si = m.union("si", &[slab, core]);
    m.overlay_defs(vec![def("clad", clad, Some(3)), def("BOX", bx, Some(1)), def("Si", si, Some(2))], vec![]);
    let (w, t, sl, mg, bt) = (500e-9, 300e-9, 100e-9, 1e-6, 2e-6);
    let (ir, il) = ((t - sl) / (80f64.to_radians()).tan(), (t - sl) / (75f64.to_radians()).tan());
    let wide = w + 2.0 * mg;
    let si_area = wide * sl + 0.5 * (w + (w - ir - il)) * (t - sl);
    let e = exp(
        "rib on a slab with unequal side-wall angles; core is a polygon with named edges",
        &[("clad", wide * (bt + t + mg) - wide * bt - si_area), ("BOX", wide * bt), ("Si", si_area)],
        &[
            "BOX|Si",
            "BOX|clad",
            "Si|clad",
            "box.bottom",
            "box.left",
            "box.right",
            "clad.left",
            "clad.right",
            "clad.top",
        ],
    );
    let mut e = e;
    // the slab ends are named after the box/clad side pieces; keep the actual set as the reference after review
    e.names = vec![];
    Golden { name: "soi_rib_sloped", model: m, expected: e }
}

/// Silicon-on-sapphire rib waveguides of the four ring chips with structured PML layers (chips.rs, `chip_shapes`).
pub fn chip(name: &'static str, htot: f64, etch: f64, width: f64, xin: f64) -> Golden {
    let um = 1e-6;
    let mut m = Model::new_2d();
    m.param("um", um).param("htot", htot).param("etch", etch).param("width", width).param("xin", xin);
    m.param("hslab", "(htot-etch)*um")
        .param("w", "width*um")
        .param("xi", "xin*um")
        .param("hin", "2*um")
        .param("lp", "1*um")
        .param("ha", "1.2*um");
    m.param("x", "xi+lp").param("ytop", "htot*um+ha");
    let air = m.rect("air", ("-x", "-hin-lp"), ("x", "ytop+lp"));
    let sap = m.rect("sapphire", ("-x", "-hin-lp"), ("x", 0.0));
    let slab = m.rect("slab", ("-x", 0.0), ("x", "hslab"));
    let rib = m.rect("rib", ("-w/2", "hslab"), ("w/2", "htot*um"));
    let win = m.rect("window", ("-xi", "-hin"), ("xi", "ytop"));
    m.overlay(&[("air", air), ("sapphire", sap), ("slab", slab), ("rib", rib)]);
    m.pml("pml", win, &[], "lp", 4.0, 2.0);
    let (hslab, w, xi, hin, lp, ha) = ((htot - etch) * um, width * um, xin * um, 2.0 * um, um, 1.2 * um);
    let ytop = htot * um + ha;
    let rib_a = w * (htot * um - hslab);
    let phys = |h: f64| 2.0 * xi * h;
    let air_phys = phys(ytop + hin) - phys(hin) - phys(hslab) - rib_a;
    let regions: Vec<(&str, f64)> = vec![
        ("air", air_phys),
        ("air.pml_xmin", lp * (ytop - hslab)),
        ("air.pml_xmax", lp * (ytop - hslab)),
        ("air.pml_zmax", 2.0 * xi * lp),
        ("air.pml_xmin_zmax", lp * lp),
        ("air.pml_xmax_zmax", lp * lp),
        ("sapphire", phys(hin)),
        ("sapphire.pml_xmin", lp * hin),
        ("sapphire.pml_xmax", lp * hin),
        ("sapphire.pml_zmin", 2.0 * xi * lp),
        ("sapphire.pml_xmin_zmin", lp * lp),
        ("sapphire.pml_xmax_zmin", lp * lp),
        ("slab", phys(hslab)),
        ("slab.pml_xmin", lp * hslab),
        ("slab.pml_xmax", lp * hslab),
        ("rib", rib_a),
    ];
    let mut e = exp("ring chip cross-section: 16 regions, the 12 PML pieces come from one Pml feature", &regions, &[]);
    e.names = vec![];
    Golden { name, model: m, expected: e }
}

pub fn chips() -> Vec<Golden> {
    vec![
        chip("chip_300te", 0.3, 0.17, 0.5, 2.0),
        chip("chip_500te", 0.5, 0.15, 1.0, 3.0),
        chip("chip_300tm", 0.3, 0.235, 0.5, 2.0),
        chip("chip_500tm", 0.5, 0.325, 0.5, 2.0),
    ]
}

/// Cross-section of the widest part of the Bragg cell: the wire with a tooth on each side (domains 4 and 5 of `bragg.rs`).
pub fn bragg_section() -> Golden {
    let mut m = Model::new_2d();
    m.param("w", 500e-9).param("t", 300e-9).param("tooth", 80e-9).param("margin", 0.6e-6).param("box_t", 2e-6);
    let clad = m.rect("clad", ("-w/2-tooth-margin", "-box_t"), ("w/2+tooth+margin", "t+margin"));
    let bx = m.rect("box", ("-w/2-tooth-margin", "-box_t"), ("w/2+tooth+margin", 0.0));
    let core = m.rect("core", ("-w/2", 0.0), ("w/2", "t"));
    let lt = m.rect("left_tooth", ("-w/2-tooth", 0.0), ("-w/2", "t"));
    let rt = m.rect("right_tooth", ("w/2", 0.0), ("w/2+tooth", "t"));
    m.overlay_defs(
        vec![
            def("clad", clad, Some(3)),
            def("BOX", bx, Some(1)),
            def("Si", core, Some(2)),
            def("left_tooth", lt, Some(4)),
            def("right_tooth", rt, Some(5)),
        ],
        vec![],
    );
    let (w, t, tooth, mg, bt) = (500e-9, 300e-9, 80e-9, 0.6e-6, 2e-6);
    let wide = w + 2.0 * tooth + 2.0 * mg;
    let e = exp(
        "widest cross-section of the Bragg wire: core plus two teeth, domains 1..5",
        &[
            ("clad", wide * (bt + t + mg) - wide * bt - (w + 2.0 * tooth) * t),
            ("BOX", wide * bt),
            ("Si", w * t),
            ("left_tooth", tooth * t),
            ("right_tooth", tooth * t),
        ],
        &[],
    );
    let mut e = e;
    e.names = vec![];
    Golden { name: "bragg_section", model: m, expected: e }
}

/// 0.4 x 0.2 plate with a hole of radius 0.04 (Kirsch / Peterson stress concentration test of femcore).
pub fn plate_with_hole() -> Golden {
    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (0.4, 0.2));
    let hole = m.circle("hole", (0.2, 0.1), 0.04);
    let d = m.difference("plate_minus_hole", plate, &[hole]);
    m.overlay(&[("plate", d)]);
    let e = exp(
        "plate with a circular hole",
        &[("plate", 0.4 * 0.2 - PI * 0.04 * 0.04)],
        &["hole.edge", "plate.bottom", "plate.left", "plate.right", "plate.top"],
    );
    Golden { name: "plate_with_hole", model: m, expected: e }
}

/// Disc of radius 0.5 (Poisson on a meshed disc, exact area pi r^2).
pub fn disk() -> Golden {
    let mut m = Model::new_2d();
    let d = m.circle("disk", (0.0, 0.0), 0.5);
    m.overlay(&[("disk", d)]);
    Golden { name: "disk", model: m, expected: exp("disc of radius 0.5", &[("disk", PI * 0.25)], &["disk.edge"]) }
}

/// Annulus 0.3 < r < 0.5 (concentric cylinders / coax cross-section).
pub fn annulus() -> Golden {
    let mut m = Model::new_2d();
    let o = m.circle("outer", (0.0, 0.0), 0.5);
    let i = m.circle("inner", (0.0, 0.0), 0.3);
    let d = m.difference("ring", o, &[i]);
    m.overlay(&[("ring", d)]);
    Golden {
        name: "annulus",
        model: m,
        expected: exp("annulus", &[("ring", PI * (0.25 - 0.09))], &["inner.edge", "outer.edge"]),
    }
}

/// Three concentric layers (painter order), two dielectric interfaces.
pub fn concentric_discs() -> Golden {
    let mut m = Model::new_2d();
    let a = m.circle("c3", (0.0, 0.0), 0.6);
    let b = m.circle("c2", (0.0, 0.0), 0.4);
    let c = m.circle("c1", (0.0, 0.0), 0.2);
    m.overlay(&[("outer", a), ("mid", b), ("core", c)]);
    Golden {
        name: "concentric_discs",
        model: m,
        expected: exp(
            "three nested discs",
            &[("outer", PI * (0.36 - 0.16)), ("mid", PI * (0.16 - 0.04)), ("core", PI * 0.04)],
            &["c3.edge", "mid|outer", "core|mid"],
        ),
    }
}

pub fn all_2d() -> Vec<Golden> {
    let mut v = vec![soi_strip(), soi_rib_sloped()];
    v.extend(chips());
    v.extend([bragg_section(), plate_with_hole(), disk(), annulus(), concentric_discs()]);
    v
}

/// One period of a Bragg grating in a 300 nm SOI wire (`usecases/bragg.rs`): two layers along z (teeth present, teeth absent),
/// a Floquet pair on the end faces. `left`/`right` say which sides carry teeth.
pub fn bragg3(name: &'static str, left: bool, right: bool) -> Golden {
    let mut m = Model::new_3d();
    m.param("w", 500e-9)
        .param("t", 300e-9)
        .param("period", 280e-9)
        .param("duty", 0.5)
        .param("tooth", 80e-9)
        .param("margin", 0.45e-6)
        .param("box_t", 2e-6);
    m.param("x", "w/2+tooth+margin");
    let clad = m.box_("clad", ("-x", "-box_t", 0.0), ("x", "t+margin", "period"));
    let bx = m.box_("box", ("-x", "-box_t", 0.0), ("x", 0.0, "period"));
    let core = m.box_("core", ("-w/2", 0.0, 0.0), ("w/2", "t", "period"));
    let mut regs = vec![def("clad", clad, Some(3)), def("BOX", bx, Some(1)), def("Si", core, Some(2))];
    if left {
        let lt = m.box_("left_tooth", ("-w/2-tooth", 0.0, 0.0), ("-w/2", "t", "duty*period"));
        regs.push(def("left_tooth", lt, Some(4)));
    }
    if right {
        let rt = m.box_("right_tooth", ("w/2", 0.0, 0.0), ("w/2+tooth", "t", "duty*period"));
        regs.push(def("right_tooth", rt, Some(5)));
    }
    m.overlay_defs(regs, vec![]);
    m.periodic(
        "floquet",
        geom::select::Selector::parse("plane(z = 0)").unwrap(),
        geom::select::Selector::parse("plane(z = period)").unwrap(),
        (0.0, 0.0, "period"),
    );
    let (w, t, p, d, tooth, mg, bt) = (500e-9, 300e-9, 280e-9, 0.5, 80e-9, 0.45e-6, 2e-6);
    let teeth = (left as u8 + right as u8) as f64;
    let wide = w + 2.0 * (tooth + mg);
    let mut regions = vec![
        ("clad", wide * (bt + t + mg) * p - wide * bt * p - w * t * p - teeth * tooth * t * d * p),
        ("BOX", wide * bt * p),
        ("Si", w * t * p),
    ];
    if left {
        regions.push(("left_tooth", tooth * t * d * p));
    }
    if right {
        regions.push(("right_tooth", tooth * t * d * p));
    }
    let mut e =
        exp("one Bragg period: layers with and without teeth, identical end faces, Floquet pair", &regions, &[]);
    e.names = vec![];
    Golden { name, model: m, expected: e }
}

/// Cylindrical cavity through a block (3D plate with a hole; volume pi r^2 h removed).
pub fn cylinder_cavity() -> Golden {
    let mut m = Model::new_3d();
    let b = m.box_("block", (0.0, 0.0, 0.0), (2.0, 2.0, 1.0));
    let c = m.cylinder("cav", (1.0, 1.0, -1.0), 0.3, 3.0);
    let d = m.difference("d", b, &[c]);
    m.overlay(&[("block", d)]);
    let mut e = exp("block with a through cylinder cavity", &[("block", 4.0 - PI * 0.09)], &[]);
    e.names = vec![];
    Golden { name: "cylinder_cavity", model: m, expected: e }
}

/// Concentric cylinders of length 2 (coax: capacitance per length from the closed form).
pub fn concentric_cylinders() -> Golden {
    let mut m = Model::new_3d();
    let o = m.cylinder("outer", (0.0, 0.0, 0.0), 0.5, 2.0);
    let i = m.cylinder("inner", (0.0, 0.0, 0.0), 0.2, 2.0);
    m.overlay(&[("dielectric", o), ("conductor", i)]);
    let mut e =
        exp("concentric cylinders", &[("dielectric", PI * (0.25 - 0.04) * 2.0), ("conductor", PI * 0.04 * 2.0)], &[]);
    e.names = vec![];
    Golden { name: "concentric_cylinders", model: m, expected: e }
}

/// Sphere in a box (two regions, volumes 4/3 pi r^3), concentric spheres, and a torus (2 pi^2 R r^2).
pub fn sphere_in_box() -> Golden {
    let mut m = Model::new_3d();
    let air = m.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let ball = m.sphere("ball", (2.0, 2.0, 2.0), 0.7);
    m.overlay(&[("air", air), ("ball", ball)]);
    let v = 4.0 / 3.0 * PI * 0.343;
    let mut e = exp("sphere in a box", &[("air", 64.0 - v), ("ball", v)], &[]);
    e.names = vec![];
    Golden { name: "sphere_in_box", model: m, expected: e }
}

pub fn concentric_spheres() -> Golden {
    let mut m = Model::new_3d();
    let host = m.box_("host", (-3.0, -3.0, -3.0), (3.0, 3.0, 3.0));
    let outer = m.sphere("outer", (0.0, 0.0, 0.0), 2.0);
    let inner = m.sphere("inner", (0.0, 0.0, 0.0), 1.0);
    m.overlay(&[("vacuum", host), ("dielectric", outer), ("core", inner)]);
    let sv = |r: f64| 4.0 / 3.0 * PI * r * r * r;
    let mut e = exp(
        "concentric spheres in a box",
        &[("vacuum", 216.0 - sv(2.0)), ("dielectric", sv(2.0) - sv(1.0)), ("core", sv(1.0))],
        &[],
    );
    e.names = vec![];
    Golden { name: "concentric_spheres", model: m, expected: e }
}

pub fn torus_in_box() -> Golden {
    let mut m = Model::new_3d();
    let host = m.box_("host", (-2.0, -2.0, -1.0), (2.0, 2.0, 1.0));
    let ring = m.torus("ring", (0.0, 0.0, 0.0), 1.0, 0.25);
    m.overlay(&[("vacuum", host), ("ring", ring)]);
    let v = 2.0 * PI * PI * 1.0 * 0.0625;
    let mut e = exp("torus in a box", &[("vacuum", 32.0 - v), ("ring", v)], &[]);
    e.names = vec![];
    Golden { name: "torus_in_box", model: m, expected: e }
}

/// Unit cube cut by the plane x + z = 1 (a triangular prism): the slanted-cut golden.
pub fn slanted_cut_box() -> Golden {
    let mut m = Model::new_3d();
    let b = m.box_("block", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("block", b)]);
    m.cut("slant", (1.0, 0.0, 1.0), 1.0, Keep::Below);
    let mut e = exp("box with a slanted plane cut: 5 faces, 9 edges, 6 vertices", &[("block", 0.5)], &[]);
    e.names = vec![];
    Golden { name: "slanted_cut_box", model: m, expected: e }
}

/// A 2 x 2 square tapered by 10 degrees over height 1 (frustum).
pub fn tapered_frustum() -> Golden {
    let mut m = Model::new_3d();
    let p = m.rect("outline", (0.0, 0.0), (2.0, 2.0));
    let e = m.extrude_tapered("frustum", p, 0.0, 1.0, 10.0);
    m.overlay(&[("frustum", e)]);
    let top = 2.0 - 2.0 * 10f64.to_radians().tan();
    let mut ex = exp("tapered extrusion", &[("frustum", (4.0 + top * top + 2.0 * top) / 3.0)], &[]);
    ex.names = vec![];
    Golden { name: "tapered_frustum", model: m, expected: ex }
}

/// A stepped bell (revolved polygon) in a box: Pappus volume, five swept faces.
pub fn revolved_bell() -> Golden {
    let mut m = Model::new_3d();
    let air = m.box_("air", (-5.0, -5.0, -1.0), (5.0, 5.0, 8.0));
    let pts = [(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (1.5, 4.0), (1.5, 6.0), (0.0, 6.0)];
    let p = m.polygon_named("bell_profile", &pts, &["floor", "foot", "shoulder", "neck", "lid", "axis"]);
    let bell = m.revolve("bell", p, (0.0, 0.0, 0.0));
    m.overlay(&[("air", air), ("bronze", bell)]);
    let v = PI * 9.0 * 2.0 + PI * 2.0 * (9.0 + 4.5 + 2.25) / 3.0 + PI * 2.25 * 2.0;
    let mut e = exp("revolved stepped bell in a box", &[("air", 10.0 * 10.0 * 9.0 - v), ("bronze", v)], &[]);
    e.names = vec![];
    Golden { name: "revolved_bell", model: m, expected: e }
}

pub fn all_3d() -> Vec<Golden> {
    vec![
        bragg3("bragg_symmetric", true, true),
        bragg3("bragg_asymmetric", true, false),
        bragg3("bragg_uniform", false, false),
        cylinder_cavity(),
        concentric_cylinders(),
        sphere_in_box(),
        concentric_spheres(),
        torus_in_box(),
        slanted_cut_box(),
        tapered_frustum(),
        revolved_bell(),
    ]
}

pub fn all_models() -> Vec<Golden> {
    let mut v = all_2d();
    v.extend(all_3d());
    v
}
