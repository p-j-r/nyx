//! The only place where the wafer frame meets the 2D kernels.
//!
//! Every public description (problems, geometry, results) is in the axes of the wafer: x lateral, y the direction of propagation, z the
//! surface normal (right-handed). The 2D kernels were written for a mesh in the x-y plane with z out of it; for the cross-section of a
//! waveguide that is the mesh `(x, z)` with the out-of-plane axis `-y`. This module turns a wafer-frame problem into that model frame
//! (`to_model_frame`), gives the result fields their wafer names again (`rename_wafer_fields`) and holds the one remaining frame
//! decision of the 3D mode analysis (`axes_of_3d_problem`). Nothing else in the crate knows both vocabularies, so the day the kernels are
//! ported to wafer axes this module goes away.
//!
//! Two kernel layouts exist, and both are tables of "model axis -> (wafer axis, sign)":
//!
//! * `Kernel::Section`, a straight waveguide cross-section: model axes (x, y, z) = wafer (x, z, -y);
//! * `Kernel::Ring`, the (r, z, theta) half-plane of a ring (`brillouin.ring_radius`): x = r, y = z, and the third component is
//!   `+theta`, which is wafer `+y` (r, theta, z is right-handed). Crystal data are averaged about the ring axis, which makes the sign of
//!   the third axis irrelevant for orientations (test `ring_frame_adapter`); only the sign of the third result component differs.
//!
//! Tensor orientations (Euler angles) and boundary-condition components follow the table; result fields are renamed component by component
//! from the same table (`wafer_field_name`), so a field the table cannot place is reported rather than guessed.

use crate::material::{self, M3};
use crate::mesh::SimplexMesh;
use crate::problem::*;

/// Which family of 2D kernels a problem goes to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kernel {
    Section,
    Ring,
}

impl Kernel {
    pub fn of(spec: &ProblemSpec) -> Kernel {
        if spec.physics == Physics::Brillouin && spec.brillouin.as_ref().is_some_and(|b| b.ring_radius.is_some()) {
            Kernel::Ring
        } else {
            Kernel::Section
        }
    }

    /// For each model axis (x, y, z of the kernel) the wafer axis (0 = x, 1 = y, 2 = z) it lies along and the sign of that axis.
    pub fn axes(self) -> [(usize, f64); 3] {
        match self {
            Kernel::Section => [(0, 1.0), (2, 1.0), (1, -1.0)],
            Kernel::Ring => [(0, 1.0), (2, 1.0), (1, 1.0)],
        }
    }
}

/// The rotation from the wafer axes (x, y, z) to the axes of a straight cross-section (mesh x, mesh y = z, out of plane = -y):
/// the rows of `Kernel::Section.axes()`.
pub const WAFER_TO_MODEL: M3 = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]];

fn mul3(a: &M3, b: &M3) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

/// Euler angles of a crystal whose orientation is `euler_deg` in the wafer frame, expressed in the 2D model frame (for a ring this
/// is the same after the averaging about the ring axis).
pub fn euler_to_model_frame(euler_deg: [f64; 3]) -> [f64; 3] {
    material::euler_from_rotation(&mul3(&WAFER_TO_MODEL, &material::rotation(euler_deg)))
}

/// Euler angles in the wafer frame of a crystal given in the 2D model frame (the inverse of `euler_to_model_frame`).
pub fn euler_from_model_frame(euler_deg: [f64; 3]) -> [f64; 3] {
    let st: M3 = [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]; // transpose of WAFER_TO_MODEL
    material::euler_from_rotation(&mul3(&st, &material::rotation(euler_deg)))
}

/// A layer of a 2D cross-section written along the wafer axes x and z, along the mesh axes x and y of the solvers.
fn pml_to_model_frame(p: &Pml) -> Pml {
    Pml { y: p.z, z: PmlDir::Off, ..p.clone() }
}

/// The problem in the axes the solvers use. On a 2D mesh the crystal orientations, the components of constraints and
/// loads, the layers and the sign of the out-of-plane wave number are turned into the model frame; everything else is unchanged.
pub fn to_model_frame(spec: &ProblemSpec, mesh_dim: usize) -> ProblemSpec {
    let mut s = spec.clone();
    if mesh_dim != 2 {
        return s;
    }
    let axes = Kernel::of(spec).axes();
    // exp(-i k y) with the model's third axis z = sign * y is exp(-i (sign k) z): a positive k is a wave along +y
    s.out_of_plane = s.out_of_plane.map(|k| k * axes[2].1);
    for p in s.materials.pml.values_mut() {
        *p = pml_to_model_frame(p);
    }
    for a in s.materials.aniso.values_mut() {
        a.euler_deg = euler_to_model_frame(a.euler_deg);
    }
    for bc in &mut s.bcs {
        match &mut bc.kind {
            BcKind::Fixed { fixed, values } => {
                let (f0, v0) = (*fixed, *values);
                for i in 0..3 {
                    fixed[i] = f0[axes[i].0];
                    values[i] = axes[i].1 * v0[axes[i].0];
                }
            }
            BcKind::Load(f) => {
                let f0 = *f;
                for i in 0..3 {
                    f[i] = axes[i].1 * f0[axes[i].0];
                }
            }
            _ => {}
        }
    }
    if let Some(b) = s.brillouin.as_mut() {
        for p in b.optical_pml.values_mut() {
            *p = pml_to_model_frame(p);
        }
        for p in b.photoelastic.values_mut() {
            p.euler_deg = euler_to_model_frame(p.euler_deg);
        }
        for a in b.viscosity.values_mut() {
            a.euler_deg = euler_to_model_frame(a.euler_deg);
        }
    }
    s
}

/// Refuses what a wafer-frame 2D problem cannot say: there is no y axis in the plane of a cross-section.
pub fn check_wafer_input(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<(), String> {
    if mesh.dim == 2 {
        let bad = spec
            .materials
            .pml
            .iter()
            .map(|(d, p)| ("domain", *d, p))
            .chain(spec.brillouin.iter().flat_map(|b| b.optical_pml.iter().map(|(d, p)| ("optical domain", *d, p))))
            .find(|(_, _, p)| p.y != PmlDir::Off);
        if let Some((what, d, _)) = bad {
            return Err(format!("the perfectly matched layer of {what} {d} absorbs along y, which is out of the plane of a 2D cross-section (x-z): use z for the direction of the surface normal"));
        }
    }
    Ok(())
}

// ------------------------------------------------------------------ result fields

const LETTERS: [char; 3] = ['x', 'y', 'z'];

fn model_axis(c: char) -> Option<usize> {
    LETTERS.iter().position(|&l| l == c)
}

/// Wafer letter and sign of a model component letter.
fn component(axes: &[(usize, f64); 3], c: char) -> Option<(char, f64)> {
    model_axis(c).map(|i| (LETTERS[axes[i].0], axes[i].1))
}

/// The wafer name and sign of a field of the 2D kernels, or `None` when the table cannot place it. The vocabulary of the kernels:
///
/// * frame-free scalars and norms: `V`, `T`, `|E|`, `|q|`, `|u|`, `von Mises`;
/// * a vector component `u?`, `E?`, `q?` (with an optional `Re ` in front, or `|...|` around it for the magnitude of that component);
/// * the transverse electric components `Etx`, `Ety` (in the plane of the cross-section) and their magnitude `|Et|`;
/// * a stress component `s??`.
pub fn wafer_field_name(model: &str, kernel: Kernel) -> Option<(String, f64)> {
    let axes = kernel.axes();
    if matches!(model, "V" | "T" | "|E|" | "|q|" | "|u|" | "von Mises") {
        return Some((model.to_string(), 1.0));
    }
    if model == "|Et|" {
        return Some(("|Exz|".to_string(), 1.0));
    }
    // |uz|: magnitude of one component
    if let Some(inner) = model.strip_prefix('|').and_then(|r| r.strip_suffix('|')) {
        let (v, c) = split_vector(inner)?;
        let (w, _) = component(&axes, c)?;
        return Some((format!("|{v}{w}|"), 1.0));
    }
    let (re, rest) = match model.strip_prefix("Re ") {
        Some(r) => ("Re ", r),
        None => ("", model),
    };
    // transverse components Etx, Ety
    if let Some(c) = rest.strip_prefix("Et") {
        let mut it = c.chars();
        let (c, None) = (it.next()?, it.next()) else { return None };
        if model_axis(c)? > 1 {
            return None;
        }
        let (w, s) = component(&axes, c)?;
        return Some((format!("{re}E{w}"), s));
    }
    if let Some((v, c)) = split_vector(rest) {
        let (w, s) = component(&axes, c)?;
        return Some((format!("{re}{v}{w}"), s));
    }
    // stress sxy
    if re.is_empty() {
        let cs: Vec<char> = rest.strip_prefix('s')?.chars().collect();
        if cs.len() == 2 {
            let ((a, sa), (b, sb)) = (component(&axes, cs[0])?, component(&axes, cs[1])?);
            let (a, b) = canonical_pair(a, b);
            return Some((format!("s{a}{b}"), sa * sb));
        }
    }
    None
}

/// `u`, `E` or `q` followed by one axis letter
fn split_vector(s: &str) -> Option<(char, char)> {
    let mut it = s.chars();
    let (v, c) = (it.next()?, it.next()?);
    (it.next().is_none() && matches!(v, 'u' | 'E' | 'q') && model_axis(c).is_some()).then_some((v, c))
}

/// the order the solvers name the off-diagonal stresses in: xx yy zz, xy yz xz
fn canonical_pair(a: char, b: char) -> (char, char) {
    match (a, b) {
        ('y', 'x') => ('x', 'y'),
        ('z', 'y') => ('y', 'z'),
        ('z', 'x') => ('x', 'z'),
        p => p,
    }
}

/// Gives the result fields of a 2D kernel their wafer names (and signs). Returns the names the table could not place; those fields are
/// left as they are.
pub fn rename_wafer_fields(r: &mut SolveResult, kernel: Kernel) -> Vec<String> {
    let mut unmapped: Vec<String> = vec![];
    for sol in &mut r.solutions {
        for f in &mut sol.fields {
            match wafer_field_name(&f.name, kernel) {
                Some((name, sign)) => {
                    f.name = name;
                    if sign < 0.0 {
                        f.values.iter_mut().for_each(|v| *v = -*v);
                    }
                }
                None => {
                    if !unmapped.contains(&f.name) {
                        unmapped.push(f.name.clone());
                    }
                }
            }
        }
    }
    unmapped
}

// ------------------------------------------------------------------ 3D mode analysis

/// The three axes of a 3D problem that play the roles of the cross-section: `(lateral, normal, propagation)` as indices of x, y, z,
/// and whether that triple is left-handed (then the sense of rotation of the field must be read the other way round).
///
/// The propagation axis is the axis a Floquet pair of the problem is periodic along (the first one found); without periodic
/// conditions it is `y`. The surface normal is `z` unless the problem is periodic along z, then `y`.
pub fn axes_of_3d_problem(spec: &ProblemSpec, mesh: &SimplexMesh) -> ([usize; 3], bool) {
    let mut long = 1;
    // a Floquet pair maps the facets of `source` onto those of `dest`: the periodic axis is the one the two groups differ in
    if let Some(p) = spec.periodic.first() {
        let centre = |labels: &[i32]| -> [f64; 3] {
            let (mut c, mut n) = ([0.0f64; 3], 0.0f64);
            for f in 0..mesh.n_facets() {
                if labels.contains(&mesh.facet_labels[f]) {
                    for &v in mesh.facet(f) {
                        for k in 0..3 {
                            c[k] += mesh.point(v as usize)[k];
                        }
                        n += 1.0;
                    }
                }
            }
            c.map(|x| x / n.max(1.0))
        };
        let (a, b) = (centre(&p.source), centre(&p.dest));
        let d = [(a[0] - b[0]).abs(), (a[1] - b[1]).abs(), (a[2] - b[2]).abs()];
        if d.iter().cloned().fold(0.0, f64::max) > 0.0 {
            long = (0..3).max_by(|&i, &j| d[i].total_cmp(&d[j])).unwrap_or(long);
        }
    }
    let normal = if long != 2 { 2 } else { 1 };
    let lateral = (0..3).find(|&k| k != long && k != normal).unwrap_or(0);
    // (lateral, normal, long) is right-handed when it is an even permutation of (x, y, z)
    let even = matches!((lateral, normal, long), (0, 1, 2) | (1, 2, 0) | (2, 0, 1));
    ([lateral, normal, long], !even)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn model_frame_orientations_come_back() {
        let e = [30.0, 40.0, 10.0];
        let back = euler_from_model_frame(euler_to_model_frame(e));
        let (r0, r1) = (material::rotation(e), material::rotation(back));
        assert!((0..3).all(|i| (0..3).all(|j| (r0[i][j] - r1[i][j]).abs() < 1e-12)));
    }

    #[test]
    fn the_rotation_is_the_table() {
        let axes = Kernel::Section.axes();
        for i in 0..3 {
            for j in 0..3 {
                let want = if axes[i].0 == j { axes[i].1 } else { 0.0 };
                assert_eq!(WAFER_TO_MODEL[i][j], want);
            }
        }
    }

    #[test]
    fn the_vocabulary_of_the_kernels_has_wafer_names() {
        let s = Kernel::Section;
        let name = |m: &str| wafer_field_name(m, s).map(|(n, g)| (n, g as i32));
        assert_eq!(name("uy"), Some(("uz".into(), 1)));
        assert_eq!(name("uz"), Some(("uy".into(), -1)));
        assert_eq!(name("|uz|"), Some(("|uy|".into(), 1)));
        assert_eq!(name("Re uz"), Some(("Re uy".into(), -1)));
        assert_eq!(name("Etx"), Some(("Ex".into(), 1)));
        assert_eq!(name("Ety"), Some(("Ez".into(), 1)));
        assert_eq!(name("Re Ety"), Some(("Re Ez".into(), 1)));
        assert_eq!(name("Ez"), Some(("Ey".into(), -1)));
        assert_eq!(name("|Et|"), Some(("|Exz|".into(), 1)));
        assert_eq!(name("sxy"), Some(("sxz".into(), 1)));
        assert_eq!(name("sxz"), Some(("sxy".into(), -1)));
        assert_eq!(name("syz"), Some(("syz".into(), -1)));
        assert_eq!(name("szz"), Some(("syy".into(), 1)));
        assert_eq!(name("von Mises"), Some(("von Mises".into(), 1)));
        assert_eq!(name("qy"), Some(("qz".into(), 1)));
        assert_eq!(name("Etz"), None);
        assert_eq!(name("something else"), None);
        // the third component of a ring is +theta = +y
        assert_eq!(wafer_field_name("uz", Kernel::Ring), Some(("uy".into(), 1.0)));
    }
}
