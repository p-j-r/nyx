//! Analysis of optical modes: the entry point of nyx's "analysis section".
//!
//! Takes the modes of a mode study (complex Ex, Ey, Ez at sample points, with area weights) and produces, per mode, the
//! polarization report of `polarization` (integrated Stokes vector, Poincaré point, degree of polarization, TE and
//! longitudinal fractions, class), a table, a Poincaré-sphere picture, and an analysis of degenerate modes.
//!
//! Degenerate modes (effective indices equal within a tolerance) are the one place where a mode's polarization is not
//! well defined. If two modes share an effective index, any unitary mixture of the two is an equally good pair of
//! modes, and the eigensolver returns whichever mixture it happens to find. The polarization of each returned mode then
//! depends on that arbitrary choice, but two things do not:
//!
//! * the *group Stokes vector*: the sum of the Stokes vectors of the (equally normalized) members, which is unchanged
//!   by any unitary mixing within the group, and
//! * the *resolved pair*: for a pair, the mixture in which one chosen Stokes operator (x/y, +-45 degrees, or circular)
//!   is diagonal, i.e. the two modes of definite polarization of that kind. That pair is the same whatever the solver
//!   returned.
//!
//! The mixing assumes that the members are orthogonal in the plain sum of |E|^2 over the samples (true for lossless
//! modes of a symmetric structure); the overlap is reported so that a violation is visible. Only pairs are resolved; a
//! larger group gets its group Stokes vector only.

use crate::polarization::*;
use crate::optical::{mode_samples, OpticalSet};
use crate::mesh::SimplexMesh;
use crate::problem::{Scalar, Solution, SolveResult};
use std::fmt::Write;

fn cmul(a: Cx, b: Cx) -> Cx {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn cconj_mul(a: Cx, b: Cx) -> Cx {
    cmul([a[0], -a[1]], b)
}

fn abs2(a: Cx) -> f64 {
    a[0] * a[0] + a[1] * a[1]
}

/// Names of the field components in a mode `Solution` (real and imaginary parts as separate fields).
#[derive(Clone, Debug)]
pub struct FieldNames {
    pub ex: (String, String),
    pub ey: (String, String),
    /// `None`: the longitudinal component is not available and counts as zero
    pub ez: Option<(String, String)>,
}

impl Default for FieldNames {
    /// A guess at the naming; `modes_from_solutions` lists the fields that exist when a name is not found.
    fn default() -> Self {
        let p = |a: &str, b: &str| (a.to_string(), b.to_string());
        FieldNames { ex: p("Ex_re", "Ex_im"), ey: p("Ey_re", "Ey_im"), ez: Some(p("Ez_re", "Ez_im")) }
    }
}

#[derive(Clone, Debug)]
pub struct ModeInput {
    pub label: String,
    /// effective index (re, im)
    pub neff: (f64, f64),
    pub samples: Vec<Sample>,
}

/// Mode inputs from the solutions of a mode study. `weights` has one entry per sample (for P1 fields one per vertex, see
/// `p1_vertex_weights`); the first `weights.len()` values of each field are used, which for P2 fields are the vertices.
pub fn modes_from_solutions(solutions: &[Solution], names: &FieldNames, weights: &[f64]) -> Result<Vec<ModeInput>, String> {
    let n = weights.len();
    let mut out = vec![];
    for (k, sol) in solutions.iter().enumerate() {
        let get = |name: &str| -> Result<&[f64], String> {
            let f = sol.field(name).ok_or_else(|| {
                let have: Vec<&str> = sol.fields.iter().map(|f| f.name.as_str()).collect();
                format!("mode {} ('{}') has no field '{}'; the fields are: {}", k + 1, sol.label, name, have.join(", "))
            })?;
            if f.values.len() < n {
                return Err(format!("field '{}' of mode {} has {} values but there are {} weights", name, k + 1, f.values.len(), n));
            }
            Ok(&f.values[..n])
        };
        let (exr, exi) = (get(&names.ex.0)?, get(&names.ex.1)?);
        let (eyr, eyi) = (get(&names.ey.0)?, get(&names.ey.1)?);
        let ez = match &names.ez {
            Some((r, i)) => Some((get(r)?, get(i)?)),
            None => None,
        };
        let samples = (0..n)
            .map(|i| Sample { ex: [exr[i], exi[i]], ey: [eyr[i], eyi[i]], ez: ez.map_or([0.0; 2], |(r, im)| [r[i], im[i]]), w: weights[i] })
            .collect();
        out.push(ModeInput { label: sol.label.clone(), neff: sol.neff.unwrap_or((f64::NAN, 0.0)), samples });
    }
    Ok(out)
}

/// Which Stokes operator is diagonalized to resolve a degenerate pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ResolveAxis {
    /// x versus y: the pair becomes quasi-TE and quasi-TM
    S1,
    /// +45 versus -45 degrees
    S2,
    /// counter-clockwise versus clockwise circular
    S3,
}

#[derive(Clone, Debug)]
pub struct AnalysisOptions {
    pub class: ClassOptions,
    /// what the number beside each mode is called in the table: "neff" for waveguide modes, "f (THz)" for cavity modes
    pub value_label: String,
    /// modes whose effective indices differ by at most this are treated as degenerate; mesh asymmetry splits an exactly degenerate pair by
    /// about 1e-4 (a square core in uniform oxide), so the default is larger than the eigensolver tolerance
    pub degeneracy_tol: f64,
    pub resolve_axis: ResolveAxis,
}

impl Default for AnalysisOptions {
    fn default() -> Self {
        AnalysisOptions { class: ClassOptions::default(), value_label: "neff".into(), degeneracy_tol: 5e-4, resolve_axis: ResolveAxis::S1 }
    }
}

#[derive(Clone, Debug)]
pub struct ModeRow {
    pub label: String,
    pub neff: (f64, f64),
    /// `None` when the mode has no intensity
    pub report: Option<ModeReport>,
    /// index into `Analysis::groups` when the mode is degenerate with others
    pub group: Option<usize>,
}

/// The pair of modes of definite polarization inside a degenerate pair.
#[derive(Clone, Debug)]
pub struct Resolved {
    pub axis: ResolveAxis,
    /// new mode r = `mixing[r][0]` * (member 0) + `mixing[r][1]` * (member 1), members normalized to unit power
    pub mixing: [[Cx; 2]; 2],
    /// eigenvalue of the Stokes operator for each new mode (+1 and -1 for fully polarized members)
    pub eigenvalues: [f64; 2],
    pub reports: [Option<ModeReport>; 2],
}

#[derive(Clone, Debug)]
pub struct GroupReport {
    /// indices into `Analysis::rows`
    pub members: Vec<usize>,
    /// Stokes vector of the group, independent of how the solver mixed the modes
    pub invariant: Option<ModeReport>,
    /// largest `|<m|n>|` between two members for unit-power modes (0 for orthogonal modes)
    pub overlap: f64,
    pub resolved: Option<Resolved>,
    pub note: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Analysis {
    pub value_label: String,
    pub rows: Vec<ModeRow>,
    pub groups: Vec<GroupReport>,
}

fn power(m: &ModeInput) -> f64 {
    m.samples.iter().map(|s| s.w * (abs2(s.ex) + abs2(s.ey) + abs2(s.ez))).sum()
}

fn scaled(m: &ModeInput, k: f64) -> Vec<Sample> {
    m.samples.iter().map(|s| Sample { w: s.w * k, ..*s }).collect()
}

/// `conj(a) O b` for the chosen Stokes operator O, summed with the weights, where `a` and `b` are samples of two modes.
fn operator_element(axis: ResolveAxis, conv: TimeConvention, a: &[Sample], b: &[Sample]) -> Cx {
    let sgn = match conv {
        TimeConvention::MinusI => 1.0,
        TimeConvention::PlusJ => -1.0,
    };
    let mut acc = [0.0, 0.0];
    for (p, q) in a.iter().zip(b) {
        let (xx, yy) = (cconj_mul(p.ex, q.ex), cconj_mul(p.ey, q.ey));
        let (xy, yx) = (cconj_mul(p.ex, q.ey), cconj_mul(p.ey, q.ex));
        let v = match axis {
            ResolveAxis::S1 => [xx[0] - yy[0], xx[1] - yy[1]],
            ResolveAxis::S2 => [xy[0] + yx[0], xy[1] + yx[1]],
            // sgn * (-i) * (xy - yx)
            ResolveAxis::S3 => {
                let d = [xy[0] - yx[0], xy[1] - yx[1]];
                [sgn * d[1], -sgn * d[0]]
            }
        };
        acc[0] += p.w * v[0];
        acc[1] += p.w * v[1];
    }
    acc
}

fn resolve_pair(a: &ModeInput, b: &ModeInput, axis: ResolveAxis, opts: &AnalysisOptions) -> Result<Resolved, String> {
    if a.samples.len() != b.samples.len() {
        return Err("the two modes are not sampled at the same points".into());
    }
    let (pa, pb) = (power(a), power(b));
    if !(pa > 0.0 && pb > 0.0) {
        return Err("a member has no intensity".into());
    }
    let ua = scaled(a, 1.0 / pa);
    let ub = scaled(b, 1.0 / pb);
    let conv = opts.class.time;
    let h00 = operator_element(axis, conv, &ua, &ua)[0];
    let h11 = operator_element(axis, conv, &ub, &ub)[0];
    let h01 = operator_element(axis, conv, &ua, &ub);
    let mean = 0.5 * (h00 + h11);
    let rad = (0.25 * (h00 - h11).powi(2) + abs2(h01)).sqrt();
    let lam = [mean + rad, mean - rad];
    // eigenvector of the larger eigenvalue; the other one is orthogonal to it
    let (mut x, mut y) = if abs2(h01) > 1e-30 * (1.0 + rad * rad) { (h01, [lam[0] - h00, 0.0]) } else if h00 >= h11 { ([1.0, 0.0], [0.0, 0.0]) } else { ([0.0, 0.0], [1.0, 0.0]) };
    let n = (abs2(x) + abs2(y)).sqrt();
    x = [x[0] / n, x[1] / n];
    y = [y[0] / n, y[1] / n];
    let v0 = [x, y];
    let v1 = [[-y[0], y[1]], [x[0], -x[1]]];
    let mixing = [v0, v1];
    let mut reports: [Option<ModeReport>; 2] = [None, None];
    for (r, c) in mixing.iter().enumerate() {
        let mixed: Vec<Sample> = ua
            .iter()
            .zip(&ub)
            .map(|(p, q)| {
                let add = |u: Cx, v: Cx| {
                    let (s, t) = (cmul(c[0], u), cmul(c[1], v));
                    [s[0] + t[0], s[1] + t[1]]
                };
                Sample { ex: add(p.ex, q.ex), ey: add(p.ey, q.ey), ez: add(p.ez, q.ez), w: p.w }
            })
            .collect();
        reports[r] = analyze(&mixed, &opts.class);
    }
    Ok(Resolved { axis, mixing, eigenvalues: lam, reports })
}

/// Analyze all modes. Rows are in the order of `modes`.
pub fn analyze_modes(modes: &[ModeInput], opts: &AnalysisOptions) -> Analysis {
    let mut rows: Vec<ModeRow> = modes.iter().map(|m| ModeRow { label: m.label.clone(), neff: m.neff, report: analyze(&m.samples, &opts.class), group: None }).collect();
    let neff: Vec<(f64, f64)> = modes.iter().map(|m| m.neff).collect();
    let mut groups = vec![];
    for members in degenerate_groups(&neff, opts.degeneracy_tol) {
        let gi = groups.len();
        for &i in &members {
            rows[i].group = Some(gi);
        }
        let powers: Vec<f64> = members.iter().map(|&i| power(&modes[i])).collect();
        let mut note = None;
        let mut all: Vec<Sample> = vec![];
        for (&i, &p) in members.iter().zip(&powers) {
            if p > 0.0 {
                all.extend(scaled(&modes[i], 1.0 / p));
            } else {
                note = Some(format!("mode '{}' has no intensity and is left out", modes[i].label));
            }
        }
        let invariant = analyze(&all, &opts.class);
        let mut overlap = 0.0f64;
        for (a, &i) in members.iter().enumerate() {
            for (b, &j) in members.iter().enumerate().skip(a + 1) {
                if modes[i].samples.len() == modes[j].samples.len() && powers[a] > 0.0 && powers[b] > 0.0 {
                    let mut g = [0.0, 0.0];
                    for (p, q) in modes[i].samples.iter().zip(&modes[j].samples) {
                        for (u, v) in [(p.ex, q.ex), (p.ey, q.ey), (p.ez, q.ez)] {
                            let t = cconj_mul(u, v);
                            g[0] += p.w * t[0];
                            g[1] += p.w * t[1];
                        }
                    }
                    overlap = overlap.max(abs2(g).sqrt() / (powers[a] * powers[b]).sqrt());
                }
            }
        }
        let resolved = if members.len() == 2 {
            match resolve_pair(&modes[members[0]], &modes[members[1]], opts.resolve_axis, opts) {
                Ok(r) => Some(r),
                Err(e) => {
                    note = Some(e);
                    None
                }
            }
        } else {
            note = Some(format!("a group of {} modes: only the group Stokes vector is given (pairs are resolved)", members.len()));
            None
        };
        groups.push(GroupReport { members, invariant, overlap, resolved, note });
    }
    Analysis { value_label: opts.value_label.clone(), rows, groups }
}

impl Analysis {
    /// Plain-text table, one line per mode, then one block per degenerate group.
    #[allow(clippy::write_literal)] // the header cells are literals next to the label of the first column
    pub fn table(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "{:>3}  {:<14} {:>10} {:>9} {:>6} {:>6} {:>6} {:>7} {:>7} {:>7} {:>6}  {}", "#", "mode", self.value_label, format!("Im {}", self.value_label), "TE", "Ey", "P", "s1", "s2", "s3", "group", "class");
        for (i, r) in self.rows.iter().enumerate() {
            let g = r.group.map_or("-".to_string(), |g| format!("G{}", g + 1));
            match &r.report {
                Some(m) => {
                    let _ = writeln!(s, "{:>3}  {:<14} {:>10.6} {:>9.2e} {:>6.3} {:>6.3} {:>6.3} {:>+7.3} {:>+7.3} {:>+7.3} {:>6}  {}", i + 1, r.label, r.neff.0, r.neff.1, m.te_fraction, m.long_fraction, m.dop, m.sphere[0], m.sphere[1], m.sphere[2], g, m.class.label());
                }
                None => {
                    let _ = writeln!(s, "{:>3}  {:<14} {:>10.6} {:>9.2e}  no intensity", i + 1, r.label, r.neff.0, r.neff.1);
                }
            }
        }
        for (gi, g) in self.groups.iter().enumerate() {
            let names: Vec<String> = g.members.iter().map(|&i| format!("{}", i + 1)).collect();
            let _ = writeln!(s, "\nG{}: modes {} are degenerate; their individual polarization depends on the solver's choice of mixture.", gi + 1, names.join(", "));
            if let Some(m) = &g.invariant {
                let _ = writeln!(s, "  group (independent of the mixture): {}", m.summary_line());
            }
            if g.overlap > 0.05 {
                let _ = writeln!(s, "  warning: the members overlap by {:.3}; the resolved pair assumes orthogonal modes", g.overlap);
            }
            if let Some(r) = &g.resolved {
                let axis = match r.axis {
                    ResolveAxis::S1 => "x/y",
                    ResolveAxis::S2 => "+-45 deg",
                    ResolveAxis::S3 => "circular",
                };
                for (k, rep) in r.reports.iter().enumerate() {
                    if let Some(m) = rep {
                        let _ = writeln!(s, "  resolved ({axis}) {}: {}", k + 1, m.summary_line());
                    }
                }
            }
            if let Some(n) = &g.note {
                let _ = writeln!(s, "  note: {n}");
            }
        }
        s
    }

    pub fn csv(&self) -> String {
        let name = self.value_label.split_whitespace().next().unwrap_or("neff");
        let mut s = format!("mode,label,{name}_re,{name}_im,s1,s2,s3,dop,chi_deg,psi_deg,te_fraction,long_fraction,class,group\n");
        for (i, r) in self.rows.iter().enumerate() {
            let g = r.group.map_or(String::new(), |g| format!("G{}", g + 1));
            match &r.report {
                Some(m) => {
                    let _ = writeln!(s, "{},\"{}\",{:e},{:e},{:.6},{:.6},{:.6},{:.6},{:.4},{:.4},{:.6},{:.6},{},{}", i + 1, r.label.replace('"', "'"), r.neff.0, r.neff.1, m.sphere[0], m.sphere[1], m.sphere[2], m.dop, m.chi_deg, m.psi_deg, m.te_fraction, m.long_fraction, m.class.label(), g);
                }
                None => {
                    let _ = writeln!(s, "{},\"{}\",{:e},{:e},,,,,,,,,no intensity,{}", i + 1, r.label.replace('"', "'"), r.neff.0, r.neff.1, g);
                }
            }
        }
        s
    }

    /// Poincaré sphere with one numbered marker per mode on a line from the centre to its integrated state (a hybrid mode lies inside the
    /// sphere); members of a degenerate group are dashed.
    pub fn sphere_svg(&self, size: f64) -> String {
        let states: Vec<SphereState> = self.rows.iter().enumerate().filter_map(|(i, r)| r.report.as_ref().map(|m| SphereState { label: (i + 1).to_string(), s: m.sphere, dashed: r.group.is_some() })).collect();
        poincare_modes_svg(&states, size)
    }
}

// ------------------------------------------------------------------ from the mode solvers

/// Polarization analysis of the modes of a 2D mode study, from the complex fields the solver found (not from the plotted
/// fields, which are real and lose the phase between `Ex` and `Ey`). Cells inside a PML are left out. The axes are those of
/// the cross-section in the wafer frame: `x` lateral (quasi-TE), `z` the surface normal (quasi-TM), `y` the propagation axis.
pub fn analyze_optical(mesh: &SimplexMesh, set: &OpticalSet, opts: &AnalysisOptions) -> Analysis {
    let modes: Vec<ModeInput> = (0..set.modes.len())
        .map(|j| ModeInput { label: format!("mode {}", j + 1), neff: (set.modes[j].neff.re, set.modes[j].neff.im), samples: mode_samples(mesh, set, j) })
        .collect();
    let mut opts = opts.clone();
    if set.left_handed {
        opts.class.time = crate::polarization::TimeConvention::MinusI;
    }
    analyze_modes(&modes, &opts)
}

impl Analysis {
    /// Adds the result to a solver's output: one summary entry per mode in the scalars, the table in the log and the
    /// whole analysis for the analysis tab and the exports.
    pub fn attach(self, result: &mut SolveResult) {
        let (scalars, log) = self.summary();
        result.scalars.extend(scalars.into_iter().map(|(n, t)| (n, Scalar::Text(t))));
        result.log.extend(log);
        result.analysis = Some(self);
    }

    /// What `attach` adds to the scalars (one text per mode) and to the log (the table), for the physics crates whose result is not the old
    /// `SolveResult`.
    pub fn summary(&self) -> (Vec<(String, String)>, Vec<String>) {
        let mut scalars = vec![];
        for (i, r) in self.rows.iter().enumerate() {
            if let Some(m) = &r.report {
                let text = format!("{}  TE {:.3}  Ey {:.3}  P {:.3}  s = ({:+.3}, {:+.3}, {:+.3})", m.class.label(), m.te_fraction, m.long_fraction, m.dop, m.sphere[0], m.sphere[1], m.sphere[2]);
                scalars.push((format!("polarization {}", i + 1), text));
            }
        }
        let mut log = vec!["polarization of the modes (Stokes parameters; s1 > 0 is x, the lateral direction, s1 < 0 is z, the surface normal):".to_string()];
        log.extend(self.table().lines().filter(|l| !l.is_empty()).map(|l| l.to_string()));
        (scalars, log)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::problem::Field;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    /// Two orthogonal modes on 64 points: A is x-polarized (quasi-TE), B is y-polarized (quasi-TM), as for a square
    /// waveguide where both have the same effective index. Different bump shapes make the example less trivial.
    fn te_tm() -> (Vec<Sample>, Vec<Sample>) {
        let n = 64;
        let g = |k: usize, c: f64| (-((k as f64 / n as f64 - c).powi(2)) * 30.0).exp();
        let a = (0..n).map(|k| Sample { ex: [g(k, 0.5), 0.0], ey: [0.0; 2], ez: [0.0; 2], w: 1.0 }).collect();
        let b = (0..n).map(|k| Sample { ex: [0.0; 2], ey: [g(k, 0.55), 0.0], ez: [0.2 * g(k, 0.5), 0.0], w: 1.0 }).collect();
        // unit power each, so that the mixtures below are unitary mixtures
        let unit = |v: Vec<Sample>| {
            let p: f64 = v.iter().map(|s| s.w * (abs2(s.ex) + abs2(s.ey) + abs2(s.ez))).sum();
            let k = 1.0 / p.sqrt();
            let f = |a: Cx| [k * a[0], k * a[1]];
            v.into_iter().map(|s| Sample { ex: f(s.ex), ey: f(s.ey), ez: f(s.ez), ..s }).collect::<Vec<_>>()
        };
        (unit(a), unit(b))
    }

    /// c0 * u + c1 * v for two sample sets
    fn mix(u: &[Sample], v: &[Sample], c0: Cx, c1: Cx) -> Vec<Sample> {
        let add = |a: Cx, b: Cx| {
            let (s, t) = (cmul(c0, a), cmul(c1, b));
            [s[0] + t[0], s[1] + t[1]]
        };
        u.iter().zip(v).map(|(p, q)| Sample { ex: add(p.ex, q.ex), ey: add(p.ey, q.ey), ez: add(p.ez, q.ez), w: p.w }).collect()
    }

    fn input(label: &str, neff: f64, samples: Vec<Sample>) -> ModeInput {
        ModeInput { label: label.into(), neff: (neff, 0.0), samples }
    }

    #[test]
    fn ordinary_modes_get_rows_and_no_groups() {
        let (a, b) = te_tm();
        let an = analyze_modes(&[input("TE", 2.4, a), input("TM", 1.9, b)], &AnalysisOptions::default());
        assert!(an.groups.is_empty());
        assert_eq!(an.rows[0].report.as_ref().unwrap().class, PolClass::QuasiTe);
        assert_eq!(an.rows[1].report.as_ref().unwrap().class, PolClass::QuasiTm);
        assert!(an.rows[1].report.as_ref().unwrap().long_fraction > 0.0);
        assert!(an.table().contains("quasi-TE") && an.csv().lines().count() == 3);
    }

    /// The same two physical modes, returned by the solver as +-45 degree mixtures or as circular mixtures
    #[test]
    fn degenerate_pair_depends_on_the_mixture_but_the_group_and_resolved_pair_do_not() {
        let (a, b) = te_tm();
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let opts = AnalysisOptions { class: ClassOptions { time: TimeConvention::MinusI, ..Default::default() }, ..Default::default() };
        let neff = 2.2;
        let mixes: [(Cx, Cx, Cx, Cx); 3] = [
            // as they are
            ([1.0, 0.0], [0.0; 2], [0.0; 2], [1.0, 0.0]),
            // (A + B)/sqrt2 and (A - B)/sqrt2
            ([r, 0.0], [r, 0.0], [r, 0.0], [-r, 0.0]),
            // (A + iB)/sqrt2 and (A - iB)/sqrt2
            ([r, 0.0], [0.0, r], [r, 0.0], [0.0, -r]),
        ];
        let mut classes = vec![];
        let mut groups = vec![];
        for (c00, c01, c10, c11) in mixes {
            let m0 = mix(&a, &b, c00, c01);
            let m1 = mix(&a, &b, c10, c11);
            let an = analyze_modes(&[input("m1", neff, m0), input("m2", neff, m1)], &opts);
            assert_eq!(an.groups.len(), 1);
            classes.push((an.rows[0].report.as_ref().unwrap().class, an.rows[1].report.as_ref().unwrap().class));
            groups.push(an.groups[0].clone());
        }
        // what each returned mode looks like changes with the mixture
        assert_eq!(classes[0], (PolClass::QuasiTe, PolClass::QuasiTm));
        assert_eq!(classes[1], (PolClass::LinearPlus45, PolClass::LinearMinus45));
        assert_eq!(classes[2], (PolClass::Circular { ccw: true }, PolClass::Circular { ccw: false }));
        // the group Stokes vector and the resolved pair do not
        let g0 = groups[0].invariant.as_ref().unwrap();
        for g in &groups {
            let inv = g.invariant.as_ref().unwrap();
            // not exactly 0 and 0.5 only because only the second sample mode has a longitudinal part
            assert!(inv.dop < 0.02 && (inv.te_fraction - 0.5).abs() < 0.02);
            for k in 0..3 {
                assert!(close(inv.sphere[k], g0.sphere[k]), "k {k} {:?} vs {:?}", inv.sphere, g0.sphere);
            }
            assert!(g.overlap < 1e-9);
            let res = g.resolved.as_ref().unwrap();
            let c0 = res.reports[0].as_ref().unwrap();
            let c1 = res.reports[1].as_ref().unwrap();
            assert_eq!(c0.class, PolClass::QuasiTe);
            assert_eq!(c1.class, PolClass::QuasiTm);
            assert!(close(c0.te_fraction, 1.0) && close(c1.te_fraction, 0.0));
            // +1 and -1 up to the longitudinal share of the power, which is not part of the Stokes operator
            assert!((res.eigenvalues[0] - 1.0).abs() < 0.05 && (res.eigenvalues[1] + 1.0).abs() < 0.05);
        }
    }

    #[test]
    fn resolving_along_other_axes() {
        let (a, b) = te_tm();
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let m0 = mix(&a, &b, [r, 0.0], [r, 0.0]);
        let m1 = mix(&a, &b, [r, 0.0], [-r, 0.0]);
        for (axis, c0, c1) in [
            (ResolveAxis::S2, PolClass::LinearPlus45, PolClass::LinearMinus45),
            (ResolveAxis::S3, PolClass::Circular { ccw: true }, PolClass::Circular { ccw: false }),
        ] {
            let opts = AnalysisOptions { class: ClassOptions { time: TimeConvention::MinusI, ..Default::default() }, resolve_axis: axis, ..Default::default() };
            let an = analyze_modes(&[input("m1", 2.0, m0.clone()), input("m2", 2.0, m1.clone())], &opts);
            let res = an.groups[0].resolved.as_ref().unwrap();
            assert_eq!(res.reports[0].as_ref().unwrap().class, c0, "{axis:?}");
            assert_eq!(res.reports[1].as_ref().unwrap().class, c1, "{axis:?}");
        }
    }

    /// The CCW mode is always the first (eigenvalue +1); what the time convention changes is its mixture of TE and TM:
    /// (TE + i TM)/sqrt2 for exp(-iwt) and (TE - i TM)/sqrt2 for exp(+iwt)
    #[test]
    fn the_time_convention_changes_the_mixture_of_the_ccw_mode() {
        let (a, b) = te_tm();
        let (m0, m1) = (mix(&a, &b, [1.0, 0.0], [0.0; 2]), mix(&a, &b, [0.0; 2], [1.0, 0.0]));
        for (conv, sign) in [(TimeConvention::MinusI, 1.0), (TimeConvention::PlusJ, -1.0)] {
            let opts = AnalysisOptions { class: ClassOptions { time: conv, ..Default::default() }, resolve_axis: ResolveAxis::S3, ..Default::default() };
            let an = analyze_modes(&[input("m1", 2.0, m0.clone()), input("m2", 2.0, m1.clone())], &opts);
            let res = an.groups[0].resolved.as_ref().unwrap();
            assert_eq!(res.reports[0].as_ref().unwrap().class, PolClass::Circular { ccw: true });
            let ratio = cconj_mul(res.mixing[0][0], res.mixing[0][1]);
            let norm = abs2(ratio).sqrt();
            assert!(close(ratio[0] / norm, 0.0) && close(ratio[1] / norm, sign));
        }
    }

    #[test]
    fn groups_of_three_get_only_the_group_vector_and_missing_intensity_is_noted() {
        let (a, b) = te_tm();
        let zero: Vec<Sample> = a.iter().map(|s| Sample { ex: [0.0; 2], ey: [0.0; 2], ez: [0.0; 2], ..*s }).collect();
        let an = analyze_modes(&[input("a", 2.0, a.clone()), input("b", 2.0, b), input("c", 2.0, a)], &AnalysisOptions::default());
        assert_eq!(an.groups[0].members, vec![0, 1, 2]);
        assert!(an.groups[0].resolved.is_none() && an.groups[0].note.is_some());
        let an = analyze_modes(&[input("z", 2.0, zero)], &AnalysisOptions::default());
        assert!(an.rows[0].report.is_none());
        assert!(an.table().contains("no intensity") && an.csv().contains("no intensity"));
    }

    #[test]
    fn overlapping_members_are_flagged() {
        let (a, _) = te_tm();
        let an = analyze_modes(&[input("a", 2.0, a.clone()), input("a2", 2.0, a)], &AnalysisOptions::default());
        assert!(an.groups[0].overlap > 0.99);
        assert!(an.table().contains("overlap"));
    }

    #[test]
    fn solutions_are_read_by_field_name() {
        let f = |name: &str, v: Vec<f64>| Field { name: name.into(), unit: String::new(), values: v };
        let sol = Solution {
            label: "mode 1".into(),
            fields: vec![f("Ex_re", vec![1.0, 2.0, 9.0]), f("Ex_im", vec![0.0; 3]), f("Ey_re", vec![0.0; 3]), f("Ey_im", vec![0.0, 1.0, 9.0]), f("Ez_re", vec![0.0; 3]), f("Ez_im", vec![0.0; 3])],
            neff: Some((2.5, 1e-4)),
            ..Default::default()
        };
        let w = [1.0, 1.0];
        let modes = modes_from_solutions(std::slice::from_ref(&sol), &FieldNames::default(), &w).unwrap();
        assert_eq!(modes[0].samples.len(), 2);
        assert_eq!(modes[0].samples[1].ex, [2.0, 0.0]);
        assert_eq!(modes[0].samples[1].ey, [0.0, 1.0]);
        assert_eq!(modes[0].neff, (2.5, 1e-4));
        let err = modes_from_solutions(std::slice::from_ref(&sol), &FieldNames { ex: ("Er".into(), "Ei".into()), ..Default::default() }, &w).unwrap_err();
        assert!(err.contains("no field 'Er'") && err.contains("Ex_re"));
        let none_ez = FieldNames { ez: None, ..Default::default() };
        assert_eq!(modes_from_solutions(std::slice::from_ref(&sol), &none_ez, &w).unwrap()[0].samples[0].ez, [0.0; 2]);
        assert!(modes_from_solutions(&[sol], &FieldNames::default(), &[1.0; 5]).is_err());
    }

    #[test]
    fn sphere_picture_has_one_point_per_mode() {
        let (a, b) = te_tm();
        let an = analyze_modes(&[input("TE", 2.4, a), input("TM", 1.9, b)], &AnalysisOptions::default());
        let svg = an.sphere_svg(300.0);
        assert!(svg.starts_with("<svg") && svg.ends_with("</svg>") && !svg.contains("NaN"));
        assert_eq!(svg.matches("font-weight=\"700\"").count(), 2);
    }
}
