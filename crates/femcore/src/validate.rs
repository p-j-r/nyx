// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Checks that run before and after a solve so that bad numbers produce a clear error instead of a silent
//! "solution" full of NaN.
//!
//! Before: every number the user can type (mesh coordinates, material properties, boundary-condition values, study
//! parameters, thickness, ...) must be finite. NaN or infinity in a stiffness or a coordinate does not stop the
//! factorisation or the Krylov iteration; it just flows into every field. After: every result field and displacement
//! must be finite as well; the message names the solution, the field and how many values are bad. Scalar results
//! (for example Q = infinity of a lossless mode) may legitimately be infinite and are not checked.

/// Byte position of the first bare `NaN` / `inf` token outside quotes in a `Debug` string.
pub fn non_finite_token(text: &str) -> Option<usize> {
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
