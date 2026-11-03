//! Elastic stiffness of a material: the full 6x6 matrix in the textbook Voigt order (1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy, engineering
//! shear strain) in the crystal frame, in Pa. It is read from a crystal class with its independent constants, from the reduced matrix
//! (upper triangle or full), from the compliance (constants or matrix), or from isotropic constants (E and nu, lambda and mu, K and G).
//! The checks use plain messages: symmetry, positive definiteness (stability), consistency with the declared class, derived entries.

use super::reader::{MaterialError, Section, Tree};
use super::units::Dim;
use super::M6;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Unique {
    Y,
    Z,
}

/// Crystal class of the stiffness (the elastic symmetry group, not the full point group).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Class {
    Isotropic,
    Cubic,
    /// transversely isotropic about the crystal z axis
    Hexagonal,
    /// classes 32, 3m and -3m about z (sapphire, LiNbO3): c14 = -c24 = c56
    Trigonal,
    /// classes 4mm, 422, -42m, 4/mmm about z
    Tetragonal,
    Orthotropic,
    /// two-fold axis along y (default setting) or z
    Monoclinic(Unique),
    Triclinic,
}

/// How one nonzero entry follows from the independent constants.
#[derive(Clone, Copy)]
enum E {
    C(&'static str),
    Neg(&'static str),
    HalfDiff(&'static str, &'static str),
}

impl E {
    fn describe(self) -> String {
        match self {
            E::C(n) => n.to_string(),
            E::Neg(n) => format!("-{n}"),
            E::HalfDiff(a, b) => format!("({a} - {b})/2"),
        }
    }
}

pub(super) const NAMES: [[&str; 6]; 6] = [
    ["c11", "c12", "c13", "c14", "c15", "c16"],
    ["c21", "c22", "c23", "c24", "c25", "c26"],
    ["c31", "c32", "c33", "c34", "c35", "c36"],
    ["c41", "c42", "c43", "c44", "c45", "c46"],
    ["c51", "c52", "c53", "c54", "c55", "c56"],
    ["c61", "c62", "c63", "c64", "c65", "c66"],
];

fn pos(name: &str) -> (usize, usize) {
    let b = name.as_bytes();
    ((b[1] - b'1') as usize, (b[2] - b'1') as usize)
}

fn cname(i: usize, j: usize) -> &'static str {
    NAMES[i.min(j)][i.max(j)]
}

impl Class {
    pub const ALL: [Class; 9] = [
        Class::Isotropic,
        Class::Cubic,
        Class::Hexagonal,
        Class::Trigonal,
        Class::Tetragonal,
        Class::Orthotropic,
        Class::Monoclinic(Unique::Y),
        Class::Monoclinic(Unique::Z),
        Class::Triclinic,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Class::Isotropic => "isotropic",
            Class::Cubic => "cubic",
            Class::Hexagonal => "hexagonal",
            Class::Trigonal => "trigonal",
            Class::Tetragonal => "tetragonal",
            Class::Orthotropic => "orthotropic",
            Class::Monoclinic(_) => "monoclinic",
            Class::Triclinic => "triclinic",
        }
    }

    pub fn parse(s: &str, unique: Unique) -> Option<Class> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "isotropic" => Class::Isotropic,
            "cubic" => Class::Cubic,
            "hexagonal" => Class::Hexagonal,
            "trigonal" => Class::Trigonal,
            "tetragonal" => Class::Tetragonal,
            "orthotropic" | "orthorhombic" => Class::Orthotropic,
            "monoclinic" => Class::Monoclinic(unique),
            "triclinic" => Class::Triclinic,
            _ => return None,
        })
    }

    /// Every nonzero entry of the upper triangle and where it comes from.
    fn layout(self) -> Vec<((usize, usize), E)> {
        use E::*;
        let diag3 = |a: &'static str| vec![((0, 0), C(a)), ((1, 1), C(a)), ((2, 2), C(a))];
        match self {
            Class::Isotropic => {
                let mut v = diag3("c11");
                v.extend([((0, 1), C("c12")), ((0, 2), C("c12")), ((1, 2), C("c12"))]);
                v.extend([((3, 3), HalfDiff("c11", "c12")), ((4, 4), HalfDiff("c11", "c12")), ((5, 5), HalfDiff("c11", "c12"))]);
                v
            }
            Class::Cubic => {
                let mut v = diag3("c11");
                v.extend([((0, 1), C("c12")), ((0, 2), C("c12")), ((1, 2), C("c12"))]);
                v.extend([((3, 3), C("c44")), ((4, 4), C("c44")), ((5, 5), C("c44"))]);
                v
            }
            Class::Hexagonal | Class::Trigonal => {
                let mut v = vec![((0, 0), C("c11")), ((1, 1), C("c11")), ((2, 2), C("c33"))];
                v.extend([((0, 1), C("c12")), ((0, 2), C("c13")), ((1, 2), C("c13"))]);
                v.extend([((3, 3), C("c44")), ((4, 4), C("c44")), ((5, 5), HalfDiff("c11", "c12"))]);
                if self == Class::Trigonal {
                    v.extend([((0, 3), C("c14")), ((1, 3), Neg("c14")), ((4, 5), C("c14"))]);
                }
                v
            }
            Class::Tetragonal => vec![
                ((0, 0), C("c11")),
                ((1, 1), C("c11")),
                ((2, 2), C("c33")),
                ((0, 1), C("c12")),
                ((0, 2), C("c13")),
                ((1, 2), C("c13")),
                ((3, 3), C("c44")),
                ((4, 4), C("c44")),
                ((5, 5), C("c66")),
            ],
            Class::Orthotropic => {
                let mut v = Vec::new();
                for n in ["c11", "c22", "c33", "c12", "c13", "c23", "c44", "c55", "c66"] {
                    v.push((pos(n), C(n)));
                }
                v
            }
            Class::Monoclinic(u) => {
                let mut v = Vec::new();
                let extra: [&'static str; 4] = match u {
                    Unique::Y => ["c15", "c25", "c35", "c46"],
                    Unique::Z => ["c16", "c26", "c36", "c45"],
                };
                for n in ["c11", "c22", "c33", "c12", "c13", "c23", "c44", "c55", "c66"].into_iter().chain(extra) {
                    v.push((pos(n), C(n)));
                }
                v
            }
            Class::Triclinic => {
                let mut v = Vec::new();
                for i in 0..6 {
                    for j in i..6 {
                        v.push(((i, j), C(NAMES[i][j])));
                    }
                }
                v
            }
        }
    }

    /// Names of the independent constants, in the order of the class definition.
    pub fn constants(self) -> Vec<&'static str> {
        let mut v: Vec<&'static str> = Vec::new();
        for (_, e) in self.layout() {
            if let E::C(n) = e {
                if !v.contains(&n) {
                    v.push(n);
                }
            }
        }
        v.sort_by_key(|n| pos(n));
        v
    }

    /// The matrix (as `c_ijkl` arranged in Voigt without shear factors) from the independent constants.
    fn expand(self, get: &dyn Fn(&str) -> f64) -> M6 {
        let mut c = [[0.0; 6]; 6];
        for ((i, j), e) in self.layout() {
            let v = match e {
                E::C(n) => get(n),
                E::Neg(n) => -get(n),
                E::HalfDiff(a, b) => 0.5 * (get(a) - get(b)),
            };
            c[i][j] = v;
            c[j][i] = v;
        }
        c
    }

    fn rule_of(self, i: usize, j: usize) -> Option<E> {
        let (i, j) = (i.min(j), i.max(j));
        self.layout().into_iter().find(|(p, _)| *p == (i, j)).map(|(_, e)| e)
    }
}

/// Engineering-shear factor of a Voigt index (compliance entries carry one factor 2 per shear index).
fn shear(i: usize) -> f64 {
    if i >= 3 {
        2.0
    } else {
        1.0
    }
}

fn gpa(x: f64) -> String {
    let s = format!("{:.4}", x / 1e9);
    let s = s.trim_end_matches('0').trim_end_matches('.');
    format!("{} GPa", if s.is_empty() || s == "-0" { "0" } else { s })
}

/// Eigenvalues of a symmetric matrix by cyclic Jacobi rotations, ascending.
pub(super) fn sym_eigenvalues<const N: usize>(mut a: [[f64; N]; N]) -> [f64; N] {
    for _ in 0..100 {
        let mut off = 0.0;
        for (i, row) in a.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                if i != j {
                    off += v * v;
                }
            }
        }
        let diag: f64 = (0..N).map(|i| a[i][i] * a[i][i]).sum();
        if off <= 1e-30 * diag.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..N {
            for q in p + 1..N {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..N {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..N {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut e = [0.0; N];
    for (i, v) in e.iter_mut().enumerate() {
        *v = a[i][i];
    }
    e.sort_by(|x, y| x.partial_cmp(y).unwrap());
    e
}

/// Inverse of a 6x6 matrix by Gauss-Jordan elimination with partial pivoting; `None` if singular.
pub fn invert6(m: &M6) -> Option<M6> {
    let mut a = *m;
    let mut inv = [[0.0; 6]; 6];
    for (i, row) in inv.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for col in 0..6 {
        let piv = (col..6).max_by(|&x, &y| a[x][col].abs().partial_cmp(&a[y][col].abs()).unwrap())?;
        if a[piv][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, piv);
        inv.swap(col, piv);
        let d = a[col][col];
        for k in 0..6 {
            a[col][k] /= d;
            inv[col][k] /= d;
        }
        for r in 0..6 {
            if r != col {
                let f = a[r][col];
                if f != 0.0 {
                    for k in 0..6 {
                        a[r][k] -= f * a[col][k];
                        inv[r][k] -= f * inv[col][k];
                    }
                }
            }
        }
    }
    Some(inv)
}

fn scale(m: &M6) -> f64 {
    m.iter().flatten().fold(0.0_f64, |a, v| a.max(v.abs()))
}

/// First entry where `a` and `b` differ by more than `tol` times the largest entry.
fn first_difference(a: &M6, b: &M6, tol: f64) -> Option<(usize, usize)> {
    let s = scale(a).max(scale(b)).max(f64::MIN_POSITIVE);
    for i in 0..6 {
        for j in i..6 {
            if (a[i][j] - b[i][j]).abs() > tol * s {
                return Some((i, j));
            }
        }
    }
    None
}

const CLASS_TOL: f64 = 1e-6;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Elastic {
    /// stiffness in Pa, textbook Voigt order, crystal frame
    pub c: M6,
    /// the declared class, if the file named one
    pub class: Option<Class>,
    /// true if the file gave the independent constants of `class` (so the file can be written back as constants)
    pub exact_constants: bool,
}

impl Elastic {
    pub fn from_matrix(c: M6) -> Result<Elastic, String> {
        check_matrix(&c, None)?;
        Ok(Elastic { c, class: None, exact_constants: false })
    }

    /// Class and independent constants in Pa, for example `(Class::Cubic, [("c11", 165.7e9), ...])`.
    pub fn from_constants(class: Class, constants: &[(&str, f64)]) -> Result<Elastic, String> {
        let get = |n: &str| constants.iter().find(|(k, _)| *k == n).map(|(_, v)| *v);
        for n in class.constants() {
            if get(n).is_none() {
                return Err(format!("class '{}' needs {}; missing: {n}", class.name(), class.constants().join(", ")));
            }
        }
        let c = class.expand(&|n| get(n).unwrap_or(0.0));
        check_matrix(&c, Some(class))?;
        Ok(Elastic { c, class: Some(class), exact_constants: true })
    }

    /// Isotropic from Young's modulus and Poisson's ratio.
    pub fn isotropic(e: f64, nu: f64) -> Result<Elastic, String> {
        if e <= 0.0 {
            return Err("Young's modulus must be positive".into());
        }
        if nu <= -1.0 || nu >= 0.5 {
            return Err(format!("Poisson's ratio {nu} is outside the stable range -1 to 0.5"));
        }
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        Elastic::from_constants(Class::Isotropic, &[("c11", lam + 2.0 * mu), ("c12", lam)])
    }

    pub fn compliance(&self) -> M6 {
        invert6(&self.c).unwrap_or([[f64::NAN; 6]; 6])
    }

    /// The highest-symmetry class whose constants reproduce the matrix (to 1e-6 of its largest entry).
    pub fn detect_class(&self) -> Class {
        for class in Class::ALL {
            let exp = class.expand(&|n| {
                let (i, j) = pos(n);
                self.c[i][j]
            });
            if first_difference(&self.c, &exp, CLASS_TOL).is_none() {
                return class;
            }
        }
        Class::Triclinic
    }

    /// Stiffness in the solver's Voigt order (xx, yy, zz, xy, yz, xz) after rotating the crystal with Bunge angles (degrees).
    pub fn global(&self, euler_deg: [f64; 3]) -> M6 {
        super::stiffness_global(&super::Stiffness::General(self.c), euler_deg)
    }

    /// Smallest and largest eigenvalue of the stiffness in the normalised (Mandel) form, Pa.
    pub fn eigenvalue_range(&self) -> (f64, f64) {
        let e = mandel_eigenvalues(&self.c);
        (e[0], e[5])
    }

    /// Reads the `[elastic]` table.
    pub fn read(s: &mut Section) -> Result<Elastic, MaterialError> {
        let unit = s.string("unit")?;
        let cunit = s.string("compliance_unit")?;
        let unique = match s.string("unique_axis")?.as_deref() {
            None | Some("y") | Some("Y") => Unique::Y,
            Some("z") | Some("Z") => Unique::Z,
            Some(o) => return Err(s.err("unique_axis", format!("'{o}' is not an axis; write \"y\" or \"z\" (the two-fold axis of a monoclinic crystal)"))),
        };
        let class = match s.string("class")? {
            None => None,
            Some(name) => Some(Class::parse(&name, unique).ok_or_else(|| {
                s.err("class", format!("unknown class '{name}'; use isotropic, cubic, hexagonal, trigonal, tetragonal, orthotropic, monoclinic or triclinic"))
            })?),
        };
        // which kind of data is present
        let keys = s.keys();
        let ckeys: Vec<&String> = keys.iter().filter(|k| is_const_key(k, 'c')).collect();
        let skeys: Vec<&String> = keys.iter().filter(|k| is_const_key(k, 's')).collect();
        let iso_keys: Vec<&str> = ["E", "nu", "lambda", "mu", "K", "G"].into_iter().filter(|k| s.has(k)).collect();
        let kinds = [
            (!ckeys.is_empty(), "constants c11, c12, ..."),
            (!skeys.is_empty(), "compliance constants s11, s12, ..."),
            (!iso_keys.is_empty(), "isotropic constants (E, nu, ...)"),
            (s.has("matrix"), "a matrix"),
            (s.has("compliance"), "a compliance matrix"),
        ];
        let present: Vec<&str> = kinds.iter().filter(|(p, _)| *p).map(|(_, n)| *n).collect();
        if present.is_empty() {
            return Err(s.err(
                "",
                "no stiffness given; write a class with its constants (class = \"cubic\", c11 = ..., c12 = ..., c44 = ...), or matrix = [...], or E and nu",
            ));
        }
        if present.len() > 1 {
            return Err(s.err("", format!("give the stiffness in one way only, but found {}", present.join(" and "))));
        }

        let mut el = if !iso_keys.is_empty() {
            if let Some(c) = class {
                if c != Class::Isotropic {
                    return Err(s.err("class", format!("E and nu describe an isotropic material but the class is '{}'", c.name())));
                }
            }
            let q = |s: &mut Section, k: &str| s.quantity(k, Dim::Pressure, unit.as_deref());
            let (lam, mu) = if s.has("E") || s.has("nu") {
                let e = q(s, "E")?.ok_or_else(|| s.err("E", "E is given without nu"))?;
                let nu = s.number("nu")?.ok_or_else(|| s.err("nu", "nu is missing; give E together with nu"))?;
                if e <= 0.0 {
                    return Err(s.err("E", "Young's modulus must be positive"));
                }
                if nu <= -1.0 || nu >= 0.5 {
                    return Err(s.err("nu", format!("Poisson's ratio {nu} is outside the stable range -1 to 0.5")));
                }
                (e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), e / (2.0 * (1.0 + nu)))
            } else if s.has("lambda") || s.has("mu") {
                let l = q(s, "lambda")?.ok_or_else(|| s.err("lambda", "mu is given without lambda"))?;
                let m = q(s, "mu")?.ok_or_else(|| s.err("mu", "lambda is given without mu"))?;
                (l, m)
            } else {
                let k = q(s, "K")?.ok_or_else(|| s.err("K", "G is given without K"))?;
                let g = q(s, "G")?.ok_or_else(|| s.err("G", "K is given without G"))?;
                (k - 2.0 * g / 3.0, g)
            };
            Elastic::from_constants(Class::Isotropic, &[("c11", lam + 2.0 * mu), ("c12", lam)]).map_err(|m| s.err("", m))?
        } else if !ckeys.is_empty() || !skeys.is_empty() {
            let letter = if ckeys.is_empty() { 's' } else { 'c' };
            let given: Vec<String> = if ckeys.is_empty() { skeys.iter().map(|k| k.to_string()).collect() } else { ckeys.iter().map(|k| k.to_string()).collect() };
            let class = class.ok_or_else(|| {
                s.err(
                    "class",
                    format!(
                        "constants {} need a class so that the missing entries are known; add class = \"cubic\" (or hexagonal, trigonal, ...) or write the whole matrix",
                        given.join(", ")
                    ),
                )
            })?;
            let need: Vec<String> = class.constants().iter().map(|n| n.replacen('c', &letter.to_string(), 1)).collect();
            let missing: Vec<&String> = need.iter().filter(|n| !given.contains(n)).collect();
            if !missing.is_empty() {
                return Err(s.err(
                    "",
                    format!(
                        "class '{}' needs the constants {}; missing: {}",
                        class.name(),
                        need.join(", "),
                        missing.iter().map(|m| m.as_str()).collect::<Vec<_>>().join(", ")
                    ),
                ));
            }
            // read the values in SI (stiffness in Pa, compliance in 1/Pa)
            let mut vals: Vec<(String, f64)> = Vec::new();
            for k in &given {
                let (dim, u) = if letter == 'c' { (Dim::Pressure, unit.as_deref()) } else { (Dim::Compliance, cunit.as_deref()) };
                let v = s.quantity(k, dim, u)?.unwrap_or(0.0);
                vals.push((k.clone(), v));
            }
            let get_cn = |cn: &str| -> f64 {
                let key = cn.replacen('c', &letter.to_string(), 1);
                let (i, j) = pos(cn);
                let v = vals.iter().find(|(k, _)| *k == key).map(|(_, v)| *v).unwrap_or(0.0);
                if letter == 's' {
                    v / (shear(i) * shear(j))
                } else {
                    v
                }
            };
            let m = class.expand(&get_cn);
            // derived entries the user wrote anyway must agree
            for (k, _) in &vals {
                let cn = k.replacen(letter, "c", 1);
                if class.constants().contains(&cn.as_str()) {
                    continue;
                }
                let (i, j) = pos(&cn);
                let want = m[i][j];
                let got = get_cn(&cn);
                let tol = CLASS_TOL * scale(&m).max(f64::MIN_POSITIVE);
                if (got - want).abs() > tol {
                    let why = match class.rule_of(i, j) {
                        Some(e) => format!("must equal {}", e.describe().replacen('c', &letter.to_string(), 3)),
                        None => "must be zero".to_string(),
                    };
                    return Err(s.err(
                        k,
                        format!("{k} is not independent in a {} crystal: it {why}, which gives {}, but the file says {}", class.name(), fmt_entry(want, letter), fmt_entry(got, letter)),
                    ));
                }
            }
            let cm = if letter == 'c' {
                m
            } else {
                // m holds s_ijkl; scale back to the engineering compliance and invert
                let mut sm = m;
                for (i, row) in sm.iter_mut().enumerate() {
                    for (j, v) in row.iter_mut().enumerate() {
                        *v *= shear(i) * shear(j);
                    }
                }
                invert6(&sm).ok_or_else(|| s.err("", "the compliance matrix is singular, so no stiffness exists"))?
            };
            check_matrix(&cm, Some(class)).map_err(|m| s.err("", m))?;
            Elastic { c: cm, class: Some(class), exact_constants: letter == 'c' }
        } else {
            let (key, dim, u) = if s.has("matrix") { ("matrix", Dim::Pressure, unit.as_deref()) } else { ("compliance", Dim::Compliance, cunit.as_deref()) };
            let tree = s.quantities(key, dim, u)?.unwrap();
            let m = matrix_from_tree(s, key, &tree, key == "matrix")?;
            let c = if key == "matrix" {
                m
            } else {
                invert6(&m).ok_or_else(|| s.err(key, "the compliance matrix is singular, so no stiffness exists"))?
            };
            check_matrix(&c, class).map_err(|m| s.err(key, m))?;
            Elastic { c, class, exact_constants: false }
        };
        if let Some(c) = class {
            el.class = Some(c);
        }
        Ok(el)
    }

    /// The `[elastic]` table written from this value: constants if the file gave them, otherwise the matrix (upper triangle), in GPa.
    pub fn write(&self) -> toml::Table {
        let mut t = toml::Table::new();
        let f = |x: f64| toml::Value::Float(x / 1e9);
        t.insert("unit".into(), "GPa".into());
        match (self.class, self.exact_constants) {
            (Some(class), true) => {
                t.insert("class".into(), class.name().into());
                if let Class::Monoclinic(u) = class {
                    t.insert("unique_axis".into(), if u == Unique::Y { "y" } else { "z" }.into());
                }
                for n in class.constants() {
                    let (i, j) = pos(n);
                    t.insert(n.into(), f(self.c[i][j]));
                }
            }
            _ => {
                if let Some(class) = self.class {
                    t.insert("class".into(), class.name().into());
                    if let Class::Monoclinic(u) = class {
                        t.insert("unique_axis".into(), if u == Unique::Y { "y" } else { "z" }.into());
                    }
                }
                let rows: Vec<toml::Value> = (0..6).map(|i| toml::Value::Array((i..6).map(|j| f(self.c[i][j])).collect())).collect();
                t.insert("matrix".into(), toml::Value::Array(rows));
            }
        }
        t
    }
}

fn fmt_entry(x: f64, letter: char) -> String {
    if letter == 'c' {
        gpa(x)
    } else {
        format!("{:.4e} 1/Pa", x)
    }
}

fn is_const_key(k: &str, letter: char) -> bool {
    let b = k.as_bytes();
    b.len() == 3 && b[0] == letter as u8 && (b'1'..=b'6').contains(&b[1]) && (b'1'..=b'6').contains(&b[2])
}

fn mandel_eigenvalues(c: &M6) -> [f64; 6] {
    let m = |i: usize| if i >= 3 { std::f64::consts::SQRT_2 } else { 1.0 };
    let mut a = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            a[i][j] = c[i][j] * m(i) * m(j);
        }
    }
    sym_eigenvalues(a)
}

/// Symmetry (already ensured by the callers), stability and, if a class is given, consistency with it.
fn check_matrix(c: &M6, class: Option<Class>) -> Result<(), String> {
    if c.iter().flatten().any(|v| !v.is_finite()) {
        return Err("the stiffness contains a value that is not a finite number".into());
    }
    if let Some(class) = class {
        let exp = class.expand(&|n| {
            let (i, j) = pos(n);
            c[i][j]
        });
        if let Some((i, j)) = first_difference(c, &exp, CLASS_TOL) {
            let why = match class.rule_of(i, j) {
                Some(e) => format!("in a {} crystal {} must equal {} = {}", class.name(), cname(i, j), e.describe(), gpa(exp[i][j])),
                None => format!("in a {} crystal {} must be zero", class.name(), cname(i, j)),
            };
            return Err(format!(
                "the matrix does not have the symmetry of the declared class '{}': {} is {}, but {why}. Correct the entry, or declare a lower-symmetry class",
                class.name(),
                cname(i, j),
                gpa(c[i][j])
            ));
        }
    }
    let e = mandel_eigenvalues(c);
    let tiny = 1e-9 * e[5].abs().max(f64::MIN_POSITIVE);
    if e[0] <= tiny {
        let mut msg = format!(
            "the stiffness is not positive definite, so the crystal would be mechanically unstable (smallest eigenvalue {}). Typical causes: a wrong sign or a mistyped unit",
            gpa(e[0])
        );
        if class == Some(Class::Cubic) || class == Some(Class::Isotropic) {
            msg.push_str("; a cubic crystal is stable only if c11 > |c12|, c11 + 2*c12 > 0 and c44 > 0");
        } else if class.is_some() {
            msg.push_str("; check c11 > |c12|, c33 > 0 and the shear constants c44, c66 > 0");
        }
        return Err(msg);
    }
    Ok(())
}

fn matrix_from_tree(s: &Section, key: &str, tree: &Tree, symmetric_required: bool) -> Result<M6, MaterialError> {
    let rows = match tree {
        Tree::List(r) => r,
        _ => return Err(s.err(key, "must be a list of 6 rows")),
    };
    if rows.len() != 6 {
        return Err(s.err(key, format!("must have 6 rows (found {}); give all 6 x 6 entries or the upper triangle with rows of 6, 5, 4, 3, 2, 1 entries", rows.len())));
    }
    let mut m = [[f64::NAN; 6]; 6];
    let mut full = true;
    for (i, row) in rows.iter().enumerate() {
        let Tree::List(r) = row else { return Err(s.err(key, format!("row {} must be a list of numbers", i + 1))) };
        let nums: Vec<f64> = r
            .iter()
            .map(|t| match t {
                Tree::Num(x) => Ok(*x),
                _ => Err(s.err(key, format!("row {} holds something that is not a number", i + 1))),
            })
            .collect::<Result<_, _>>()?;
        if nums.len() == 6 {
            for (j, v) in nums.iter().enumerate() {
                m[i][j] = *v;
            }
        } else if nums.len() == 6 - i {
            full = false;
            for (k, v) in nums.iter().enumerate() {
                m[i][i + k] = *v;
            }
        } else {
            return Err(s.err(key, format!("row {} has {} entries; expected 6 (full matrix) or {} (upper triangle starting at the diagonal)", i + 1, nums.len(), 6 - i)));
        }
    }
    if !full {
        for (i, row) in rows.iter().enumerate() {
            if let Tree::List(r) = row {
                if r.len() == 6 {
                    return Err(s.err(key, format!("row {} has 6 entries but the rows before it are triangular; use either all full rows or the triangle", i + 1)));
                }
            }
        }
        for i in 0..6 {
            for j in 0..i {
                m[i][j] = m[j][i];
            }
        }
    } else if symmetric_required {
        let sc = scale(&m).max(f64::MIN_POSITIVE);
        for i in 0..6 {
            for j in i + 1..6 {
                if (m[i][j] - m[j][i]).abs() > 1e-6 * sc {
                    return Err(s.err(
                        key,
                        format!(
                            "the matrix is not symmetric: {} = {} but {} = {}. A stiffness matrix must be symmetric; check the entries or give only the upper triangle",
                            NAMES[i][j],
                            gpa(m[i][j]),
                            NAMES[j][i],
                            gpa(m[j][i])
                        ),
                    ));
                }
            }
        }
        for i in 0..6 {
            for j in 0..i {
                let v = 0.5 * (m[i][j] + m[j][i]);
                m[i][j] = v;
                m[j][i] = v;
            }
        }
    } else {
        let sc = scale(&m).max(f64::MIN_POSITIVE);
        for i in 0..6 {
            for j in i + 1..6 {
                if (m[i][j] - m[j][i]).abs() > 1e-6 * sc {
                    return Err(s.err(key, format!("the compliance matrix is not symmetric: s{}{} differs from s{}{}", i + 1, j + 1, j + 1, i + 1)));
                }
                let v = 0.5 * (m[i][j] + m[j][i]);
                m[i][j] = v;
                m[j][i] = v;
            }
        }
    }
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigenvalues_of_a_known_matrix() {
        let e = sym_eigenvalues([[2.0, 1.0], [1.0, 2.0]]);
        assert!((e[0] - 1.0).abs() < 1e-12 && (e[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_times_matrix_is_identity() {
        let el = Elastic::from_constants(Class::Cubic, &[("c11", 165.7e9), ("c12", 63.9e9), ("c44", 79.6e9)]).unwrap();
        let s = el.compliance();
        for i in 0..6 {
            for j in 0..6 {
                let p: f64 = (0..6).map(|k| el.c[i][k] * s[k][j]).sum();
                assert!((p - (i == j) as u8 as f64).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn classes_have_the_right_counts() {
        let n: Vec<usize> = Class::ALL.iter().map(|c| c.constants().len()).collect();
        assert_eq!(n, vec![2, 3, 5, 6, 6, 9, 13, 13, 21]);
    }

    #[test]
    fn detect_class_finds_the_symmetry() {
        let hex = Elastic::from_constants(Class::Hexagonal, &[("c11", 160e9), ("c12", 60e9), ("c13", 50e9), ("c33", 180e9), ("c44", 40e9)]).unwrap();
        assert_eq!(hex.detect_class(), Class::Hexagonal);
        let cubic = Elastic::from_constants(Class::Cubic, &[("c11", 165.7e9), ("c12", 63.9e9), ("c44", 79.6e9)]).unwrap();
        assert_eq!(cubic.detect_class(), Class::Cubic);
    }

    #[test]
    fn unstable_cubic_is_refused_with_the_criteria() {
        let e = Elastic::from_constants(Class::Cubic, &[("c11", 10e9), ("c12", 60e9), ("c44", 5e9)]).unwrap_err();
        assert!(e.contains("not positive definite") && e.contains("c11 > |c12|"), "{e}");
    }
}
