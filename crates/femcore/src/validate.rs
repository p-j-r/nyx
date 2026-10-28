//! Checks that run before and after a solve so that bad numbers produce a clear error instead of a silent
//! "solution" full of NaN.
//!
//! Before: every number the user can type (mesh coordinates, material properties, boundary-condition values, study
//! parameters, thickness, ...) must be finite. NaN or infinity in a stiffness or a coordinate does not stop the
//! factorisation or the Krylov iteration; it just flows into every field. After: every result field and displacement
//! must be finite as well; the message names the solution, the field and how many values are bad. Scalar results
//! (for example Q = infinity of a lossless mode) may legitimately be infinite and are not checked.

use crate::mesh::SimplexMesh;
use crate::problem::*;

/// Looks at `spec` and `mesh` for NaN or infinite numbers; the error names the offending item.
pub fn check_inputs(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<(), String> {
    if let Some(i) = mesh.points.iter().position(|v| !v.is_finite()) {
        return Err(format!("the mesh has a non-finite coordinate (point {}, axis {}): NaN or infinity in the geometry", i / mesh.dim.max(1), i % mesh.dim.max(1)));
    }
    let m = &spec.materials;
    for (name, v) in &m.default {
        if !v.is_finite() {
            return Err(format!("material property '{name}' (default for all domains) is {v}: it must be a finite number"));
        }
    }
    for (dom, props) in &m.domains {
        for (name, v) in props {
            if !v.is_finite() {
                return Err(format!("material property '{name}' of domain {dom} is {v}: it must be a finite number"));
            }
        }
    }
    for (i, bc) in spec.bcs.iter().enumerate() {
        let vals: Vec<f64> = match &bc.kind {
            BcKind::Potential(v) | BcKind::SurfaceCharge(v) | BcKind::Temperature(v) | BcKind::Flux(v) => vec![*v],
            BcKind::Convection { h, text } => vec![*h, *text],
            BcKind::Fixed { values, .. } => values.to_vec(),
            BcKind::Load(f) => f.to_vec(),
            BcKind::Ground | BcKind::LowReflecting => vec![],
        };
        if vals.iter().any(|v| !v.is_finite()) {
            return Err(format!("boundary condition {} ({} on boundaries {:?}) has a NaN or infinite value", i + 1, bc.kind.name(), bc.labels));
        }
    }
    let study_bad = match &spec.study {
        Study::Stationary | Study::Static => None,
        Study::Transient { t_end, .. } => (!t_end.is_finite()).then_some("end time"),
        Study::Eigenfrequency { f_search, .. } => (!f_search.is_finite()).then_some("search frequency"),
        Study::Modes { freq, n_target, .. } => (!freq.is_finite()).then_some("frequency").or((!n_target.is_finite()).then_some("target effective index")),
    };
    if let Some(what) = study_bad {
        return Err(format!("the study's {what} is NaN or infinite"));
    }
    if !spec.thickness.is_finite() {
        return Err("the 2D thickness is NaN or infinite".into());
    }
    if let Some(kz) = spec.out_of_plane {
        if !kz.is_finite() {
            return Err("the out-of-plane wave number is NaN or infinite".into());
        }
    }
    // everything else (tensors, PML data, thermal expansion, Brillouin settings): scan the printed form of the spec,
    // outside of quoted strings, for the words NaN and inf
    let text = format!("{spec:?}");
    if let Some(at) = non_finite_token(&text) {
        let from = text[..at].char_indices().rev().nth(70).map_or(0, |(i, _)| i);
        return Err(format!("the problem description contains a NaN or infinite number, near: ...{}", &text[from..(at + 5).min(text.len())]));
    }
    Ok(())
}

/// Byte position of the first bare `NaN` / `inf` token outside quotes in a `Debug` string.
fn non_finite_token(text: &str) -> Option<usize> {
    let b = text.as_bytes();
    let mut in_str = false;
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'\\' if in_str => i += 1,
            b'"' => in_str = !in_str,
            b'N' | b'i' if !in_str => {
                let word = |w: &str| text[i..].starts_with(w) && !b.get(i + w.len()).is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_') && (i == 0 || !(b[i - 1].is_ascii_alphanumeric() || b[i - 1] == b'_'));
                if word("NaN") || word("inf") {
                    return Some(i);
                }
            }
            _ => {}
        }
        i += 1;
    }
    None
}

/// Every field and displacement array of the result must be finite.
pub fn check_result(result: &SolveResult) -> Result<(), String> {
    for (k, sol) in result.solutions.iter().enumerate() {
        for f in &sol.fields {
            let bad = f.values.iter().filter(|v| !v.is_finite()).count();
            if bad > 0 {
                return Err(format!("the solution is not finite: field '{}' of solution {} ('{}') has {bad} NaN/infinite values out of {}; check materials, boundary conditions and the mesh", f.name, k + 1, sol.label, f.values.len()));
            }
        }
        if let Some(d) = &sol.disp {
            let bad = d.iter().filter(|v| !v.is_finite()).count();
            if bad > 0 {
                return Err(format!("the solution is not finite: the displacement of solution {} ('{}') has {bad} NaN/infinite values out of {}", k + 1, sol.label, d.len()));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_bare_nan_and_inf_but_not_inside_strings() {
        assert_eq!(non_finite_token("a: 1.0, b: NaN"), Some(11));
        assert!(non_finite_token("t0: \"inf NaN\", x: 2.0").is_none());
        assert!(non_finite_token("name: \"a\\\"NaN\"").is_none());
        assert!(non_finite_token("Info { x: 1.0 }").is_none());
        assert!(non_finite_token("v: -inf").is_some());
    }
}
