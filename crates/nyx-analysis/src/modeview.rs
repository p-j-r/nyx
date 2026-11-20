//! A mode shape on a section in a fixed phase, so that it can be drawn and compared (ADR 99).
//!
//! An eigensolver returns a mode up to a complex factor `exp(i phi)`: every `phi` gives an equally good eigenvector. A picture shows one instant,
//! `Re(u)`, so a mode returned with `phi` near a quarter turn would look empty. The mode is therefore put in a phase that does not depend on what
//! the solver returned:
//!
//! 1. **Phase.** The in-plane components (all but the propagation axis) are rotated by `theta = -1/2 arg( sum_nodes sum_c u_c^2 )`, with the
//!    plain square, not the modulus square: `(exp(i phi) w)^2 = exp(2 i phi) w^2` for a real `w`, so this is the phase in which the in-plane motion
//!    is largest in `Re u`.
//! 2. **Sign.** `theta` leaves a factor of -1 open (`theta` and `theta + pi`); the mode is turned so that its largest in-plane component, at the
//!    node of largest in-plane amplitude, is positive. Modes then keep their orientation from run to run.
//! 3. **Scale.** The peak of `|u|` over the core is 1.
//!
//! Why the propagation component is then in `Im`: a guided mode goes as `u(x, z) exp(-i q y)`, and every `d/dy` of the elasticity (or Maxwell)
//! equations brings a factor `-i q`, so the in-plane components and the propagation component are a quarter period apart. After step 1 the
//! propagation component is in `Im`, and its `Re` is zero for a lossless guide whose section has a mirror symmetry along the propagation axis
//! (time reversal turns `q` into `-q`, the mirror turns it back). The mode reports how much of that component is still in `Re`
//! (`quadrature`): a number near 0 says the picture is in the expected phase; a mode with loss, an anisotropic cut without that symmetry or a
//! degenerate pair can leave more, and the figure says so instead of hiding it.
//!
//! The in-plane components of a mode with a complex frequency (a resonance with an absorbing layer) are put in phase the same way; the
//! propagation component of such a mode is 0 or not in quadrature, which `quadrature` also tells.

use crate::section::Section;
use nyx_physics::{Error, Solution};

/// A vector field kept as cell fields `<prefix>x`, `<prefix>y`, `<prefix>z` (`ux uy uz`, `Ex Ey Ez`): the cell fields of one prefix whose last
/// letter is an axis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VectorFamily {
    pub prefix: String,
    /// the axes present, in the order the fields appear
    pub axes: Vec<char>,
    /// the index of each component in the cell fields of the solution
    pub fields: Vec<usize>,
}

/// The vector families of a solution, in the order of their first component.
pub fn vector_families(sol: &Solution) -> Vec<VectorFamily> {
    let mut out: Vec<VectorFamily> = vec![];
    for (i, f) in sol.cell_fields.iter().enumerate() {
        let mut ch = f.name.chars();
        let (Some(axis), prefix) = (ch.next_back(), ch.as_str()) else { continue };
        if !matches!(axis, 'x' | 'y' | 'z') || prefix.is_empty() {
            continue;
        }
        match out.iter_mut().find(|v| v.prefix == prefix) {
            Some(v) => {
                v.axes.push(axis);
                v.fields.push(i);
            }
            None => out.push(VectorFamily { prefix: prefix.to_string(), axes: vec![axis], fields: vec![i] }),
        }
    }
    out.retain(|v| v.axes.len() >= 2);
    out
}

#[derive(Clone, Debug)]
pub struct AlignOptions {
    /// the vector field (its prefix); the first one of the solution when `None`
    pub family: Option<String>,
    /// the domain numbers that carry the mode: the phase, the sign and the scale are taken over their cells (the guide, not the substrate); all
    /// cells when `None`
    pub core: Option<Vec<i32>>,
    /// the propagation axis; the wafer frame of nyx propagates along y
    pub propagation: char,
}

impl Default for AlignOptions {
    fn default() -> Self {
        AlignOptions { family: None, core: None, propagation: 'y' }
    }
}

#[derive(Clone, Debug)]
pub struct AlignedMode {
    pub prefix: String,
    pub axes: Vec<char>,
    /// the components at the nodes of the section: `re[component][node]`, `im` the same (empty for a real field)
    pub re: Vec<Vec<f64>>,
    pub im: Vec<Vec<f64>>,
    /// the angle the mode was rotated by, radians in (-pi, pi]
    pub theta: f64,
    /// whether the sign step turned it by pi
    pub flipped: bool,
    /// the peak of `|u|` over the core before scaling
    pub peak: f64,
    /// the components that are not the propagation axis, and the propagation component
    pub in_plane: Vec<usize>,
    pub along: Option<usize>,
    /// share of the propagation component's energy that is in `Re` after the phase step: near 0 when it is in quadrature; `None` for a real field
    pub quadrature: Option<f64>,
    pub complex: bool,
}

impl AlignedMode {
    /// `|u|` at every node.
    pub fn magnitude(&self) -> Vec<f64> {
        let n = self.re.first().map_or(0, |c| c.len());
        (0..n)
            .map(|i| {
                let s: f64 = (0..self.re.len())
                    .map(|c| self.re[c][i].powi(2) + self.im.get(c).and_then(|v| v.get(i)).map_or(0.0, |x| x * x))
                    .sum();
                s.sqrt()
            })
            .collect()
    }

    /// The part of component `c` to draw as its second view: `Im` when it holds the larger share of the energy (the propagation component of a
    /// guided mode), `Re` otherwise; with the name of the part.
    pub fn quadrature_part(&self, c: usize) -> (&[f64], &'static str) {
        let en = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>();
        match self.im.get(c) {
            Some(im) if en(im) >= en(&self.re[c]) => (im, "Im"),
            _ => (&self.re[c], "Re"),
        }
    }
}

/// Puts solution `sol` on `section` in the phase, sign and scale described above.
pub fn align(section: &Section, sol: &Solution, opts: &AlignOptions) -> Result<AlignedMode, Error> {
    let fams = vector_families(sol);
    let fam = match &opts.family {
        Some(p) => fams.iter().find(|f| &f.prefix == p),
        None => fams.first(),
    }
    .ok_or_else(|| {
        let have: Vec<&str> = sol.cell_fields.iter().map(|f| f.name.as_str()).collect();
        Error::invalid(format!(
            "the solution ('{}') has no vector field in its cell fields (names ending in x, y, z, like `ux uy uz`); it has: {}",
            sol.label,
            have.join(", ")
        ))
    })?;
    let n = section.n_nodes();
    let nc = fam.fields.len();
    let complex = fam.fields.iter().any(|&i| !sol.cell_fields[i].im.is_empty());
    let mut re = vec![vec![0.0; n]; nc];
    let mut im = if complex { vec![vec![0.0; n]; nc] } else { vec![] };
    for (c, &fi) in fam.fields.iter().enumerate() {
        let cf = &sol.cell_fields[fi];
        for (k, &node) in section.tris.iter().enumerate() {
            let s = section.src[k] as usize;
            let v = cf.re.get(s).ok_or_else(|| {
                Error::invalid(format!("cell field '{}' is shorter than the mesh it is drawn on", cf.name))
            })?;
            re[c][node as usize] = *v;
            if complex && !cf.im.is_empty() {
                im[c][node as usize] = cf.im[s];
            }
        }
    }
    // the nodes of the core
    let mut core = vec![false; n];
    for t in 0..section.n_tris() {
        if opts.core.as_ref().is_none_or(|l| l.contains(&section.labels[t])) {
            for &node in section.tri(t) {
                core[node as usize] = true;
            }
        }
    }
    if !core.iter().any(|&b| b) {
        return Err(Error::invalid("no cell of the section carries the domains named as the core of the mode"));
    }
    let along = fam.axes.iter().position(|&a| a == opts.propagation);
    let in_plane: Vec<usize> = (0..nc).filter(|&c| Some(c) != along).collect();

    // 1. phase
    let mut theta = 0.0;
    if complex {
        let (mut sr, mut si) = (0.0f64, 0.0f64);
        for i in (0..n).filter(|&i| core[i]) {
            for &c in &in_plane {
                let (a, b) = (re[c][i], im[c][i]);
                sr += a * a - b * b;
                si += 2.0 * a * b;
            }
        }
        if sr.hypot(si) > 1e-300 {
            theta = -0.5 * si.atan2(sr);
        }
        let (ct, st) = (theta.cos(), theta.sin());
        for c in 0..nc {
            for i in 0..n {
                let (a, b) = (re[c][i], im[c][i]);
                re[c][i] = a * ct - b * st;
                im[c][i] = a * st + b * ct;
            }
        }
    }
    // 2. sign: the largest in-plane component at the node of largest in-plane amplitude is positive
    let amp = |i: usize| -> f64 {
        in_plane.iter().map(|&c| re[c][i].powi(2) + if complex { im[c][i].powi(2) } else { 0.0 }).sum()
    };
    let top = (0..n).filter(|&i| core[i]).max_by(|&a, &b| amp(a).total_cmp(&amp(b))).unwrap_or(0);
    let lead = in_plane.iter().copied().max_by(|&a, &b| re[a][top].abs().total_cmp(&re[b][top].abs()));
    let mut flipped = false;
    if let Some(c) = lead {
        if re[c][top] < 0.0 {
            flipped = true;
            re.iter_mut().flatten().for_each(|x| *x = -*x);
            im.iter_mut().flatten().for_each(|x| *x = -*x);
            theta += if theta > 0.0 { -std::f64::consts::PI } else { std::f64::consts::PI };
        }
    }
    // 3. scale
    let mag2 =
        |i: usize| -> f64 { (0..nc).map(|c| re[c][i].powi(2) + if complex { im[c][i].powi(2) } else { 0.0 }).sum() };
    let peak = (0..n).filter(|&i| core[i]).map(mag2).fold(0.0f64, f64::max).sqrt();
    if peak > 0.0 {
        re.iter_mut().flatten().for_each(|x| *x /= peak);
        im.iter_mut().flatten().for_each(|x| *x /= peak);
    }
    // a real field (the real solver of a mode study keeps one real amplitude per component, the propagation one without its factor of i) has no
    // phase to look at: there is nothing to say about quadrature
    let quadrature = along.filter(|_| complex).and_then(|c| {
        let er: f64 = (0..n).filter(|&i| core[i]).map(|i| re[c][i].powi(2)).sum();
        let ei: f64 = (0..n).filter(|&i| core[i]).map(|i| im[c][i].powi(2)).sum();
        (er + ei > 0.0).then_some(er / (er + ei))
    });
    Ok(AlignedMode {
        prefix: fam.prefix.clone(),
        axes: fam.axes.clone(),
        re,
        im,
        theta,
        flipped,
        peak,
        in_plane,
        along,
        quadrature,
        complex,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use femcore::mesh::SimplexMesh;
    use nyx_physics::CellField;

    /// A P1 square with a "mode": in-plane (x, z) real shape times `exp(i phi)`, propagation component in quadrature.
    fn mode(phi: f64, sign: f64) -> (Section, Solution) {
        let m = SimplexMesh::new(
            2,
            vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
            vec![0, 1, 2, 0, 2, 3],
            vec![1, 1],
            vec![],
            vec![],
        );
        let s = Section::new(&m, 1, None).unwrap();
        // cell 0 holds the vertices 0, 1, 2 and cell 1 the vertices 0, 2, 3 (a field that is continuous across their shared edge)
        let xs = [0.0, 1.0, 1.0, 0.0, 1.0, 0.0];
        let zs = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0];
        let (c, sn) = (phi.cos(), phi.sin());
        let mk = |name: &str, w: &[f64], quad: bool| {
            // real shape w, in quadrature when `quad`: i w
            let (re, im): (Vec<f64>, Vec<f64>) = w
                .iter()
                .map(|v| {
                    let z = if quad { (0.0, sign * v) } else { (sign * v, 0.0) };
                    (z.0 * c - z.1 * sn, z.0 * sn + z.1 * c)
                })
                .unzip();
            CellField { name: name.into(), unit: "a.u.".into(), nloc: 3, re, im }
        };
        let ux: Vec<f64> = xs.iter().zip(&zs).map(|(x, z)| 0.3 + 0.5 * x - 0.2 * z).collect();
        let uz: Vec<f64> = xs.iter().zip(&zs).map(|(x, z)| 1.0 - 0.6 * x + 0.1 * z).collect();
        let uy: Vec<f64> = xs.iter().zip(&zs).map(|(x, z)| 0.4 * x + 0.2 * z).collect();
        let sol = Solution {
            label: "mode".into(),
            fields: vec![],
            value: None,
            eigenvalue: None,
            cell_fields: vec![mk("ux", &ux, false), mk("uy", &uy, true), mk("uz", &uz, false)],
        };
        (s, sol)
    }

    #[test]
    fn whatever_phase_the_solver_returned_the_mode_comes_out_the_same() {
        let (s, sol) = mode(0.0, 1.0);
        let base = align(&s, &sol, &AlignOptions::default()).unwrap();
        // in-plane motion is real, the propagation component is in Im, and the peak is 1
        for &c in &base.in_plane {
            assert!(base.im[c].iter().all(|x| x.abs() < 1e-12), "in-plane Im is not 0");
        }
        let y = base.along.unwrap();
        assert!(base.re[y].iter().all(|x| x.abs() < 1e-12));
        assert!(base.quadrature.unwrap() < 1e-20);
        assert_eq!(base.quadrature_part(y).1, "Im");
        let top = base.magnitude().into_iter().fold(0.0f64, f64::max);
        assert!((top - 1.0).abs() < 1e-12);
        for (phi, sign) in [(0.7, 1.0), (std::f64::consts::FRAC_PI_2, 1.0), (2.9, 1.0), (-1.1, -1.0), (0.3, -1.0)] {
            let (s2, sol2) = mode(phi, sign);
            let a = align(&s2, &sol2, &AlignOptions::default()).unwrap();
            for c in 0..3 {
                for i in 0..base.re[c].len() {
                    assert!(
                        (a.re[c][i] - base.re[c][i]).abs() < 1e-9 && (a.im[c][i] - base.im[c][i]).abs() < 1e-9,
                        "phi {phi}, sign {sign}: component {c}, node {i}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_vector_field_is_found_by_the_names_of_its_cell_fields() {
        let (_, sol) = mode(0.0, 1.0);
        let f = vector_families(&sol);
        assert_eq!(f.len(), 1);
        assert_eq!((f[0].prefix.as_str(), f[0].axes.clone()), ("u", vec!['x', 'y', 'z']));
    }
}
