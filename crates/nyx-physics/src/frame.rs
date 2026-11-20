// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The frame a result's vector fields are to be read in by an analysis (restructure phase 5, ADR 56).
//!
//! An analysis is a function of a `SolveResult`, so everything it needs that is not in the numbers has to be in the result. For the
//! polarization of optical modes that is: which wafer axis is the lateral, the normal and the propagation direction of the modes (a cavity
//! with a Floquet pair is periodic along its own propagation axis), whether that triple is left-handed (then the sense of rotation is read the
//! other way), and which domains of the mesh are absorbing layers and are left out. A physics writes this as one text scalar,
//! `FRAME_KEY`; `FieldFrame::of` reads it back. The text is meant to be read by a person too:
//! `lateral x, normal z, propagation y, left-handed, outside domains 3 4`.

use crate::result::{Scalar, SolveResult};

/// Name of the scalar that carries the frame.
pub const FRAME_KEY: &str = "field frame";

const AXES: [&str; 3] = ["x", "y", "z"];

/// Wafer axes (0 = x, 1 = y, 2 = z) of the lateral, normal and propagation direction of the modes of a result, and the domains the
/// analysis leaves out.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FieldFrame {
    pub lateral: usize,
    pub normal: usize,
    pub propagation: usize,
    /// whether (lateral, normal, propagation) is a left-handed triple
    pub left_handed: bool,
    /// domain numbers (1-based) of absorbing layers: their cells are not part of the analysed mode
    pub outside: Vec<i32>,
}

impl FieldFrame {
    /// The text written to the result.
    pub fn to_text(&self) -> String {
        let mut s = format!(
            "lateral {}, normal {}, propagation {}, {}",
            AXES[self.lateral],
            AXES[self.normal],
            AXES[self.propagation],
            if self.left_handed { "left-handed" } else { "right-handed" }
        );
        if !self.outside.is_empty() {
            let d: Vec<String> = self.outside.iter().map(|d| d.to_string()).collect();
            s.push_str(&format!(", outside domains {}", d.join(" ")));
        }
        s
    }

    /// Reads `to_text` back; `None` for anything else.
    pub fn from_text(text: &str) -> Option<FieldFrame> {
        let axis = |part: &str, word: &str| -> Option<usize> {
            let rest = part.trim().strip_prefix(word)?.trim();
            AXES.iter().position(|a| *a == rest)
        };
        let parts: Vec<&str> = text.split(',').collect();
        if parts.len() < 4 {
            return None;
        }
        let left_handed = match parts[3].trim() {
            "left-handed" => true,
            "right-handed" => false,
            _ => return None,
        };
        let mut outside = vec![];
        if let Some(p) = parts.get(4) {
            for d in p.trim().strip_prefix("outside domains")?.split_whitespace() {
                outside.push(d.parse().ok()?);
            }
        }
        Some(FieldFrame {
            lateral: axis(parts[0], "lateral")?,
            normal: axis(parts[1], "normal")?,
            propagation: axis(parts[2], "propagation")?,
            left_handed,
            outside,
        })
    }

    /// The frame a result states, if it states one.
    pub fn of(result: &SolveResult) -> Option<FieldFrame> {
        result.scalars.iter().find(|(n, _)| n == FRAME_KEY).and_then(|(_, v)| match v {
            Scalar::Text(t) => FieldFrame::from_text(t),
            Scalar::Num(_) => None,
        })
    }

    /// The scalar to push on a result.
    pub fn scalar(&self) -> (String, Scalar) {
        (FRAME_KEY.to_string(), Scalar::Text(self.to_text()))
    }

    /// Wafer axis name (`x`, `y` or `z`) of the lateral, normal and propagation direction.
    pub fn names(&self) -> [&'static str; 3] {
        [AXES[self.lateral], AXES[self.normal], AXES[self.propagation]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_text_round_trips() {
        let a = FieldFrame { lateral: 0, normal: 2, propagation: 1, left_handed: true, outside: vec![] };
        assert_eq!(a.to_text(), "lateral x, normal z, propagation y, left-handed");
        assert_eq!(FieldFrame::from_text(&a.to_text()), Some(a));
        let b = FieldFrame { lateral: 2, normal: 1, propagation: 0, left_handed: false, outside: vec![3, 4] };
        assert_eq!(FieldFrame::from_text(&b.to_text()), Some(b));
        assert_eq!(FieldFrame::from_text("something else"), None);
        assert_eq!(FieldFrame::from_text("lateral x, normal q, propagation y, left-handed"), None);
    }
}
