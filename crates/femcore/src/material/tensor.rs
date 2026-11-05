//! Second-rank tensor properties of a material (relative permittivity, thermal and electrical conductivity, thermal expansion) in the
//! crystal frame. A file gives one value (isotropic), three (diagonal), six (xx, yy, zz, yz, xz, xy as for the stiffness), a 3x3
//! array, or `{ perp = .., par = .. }` for a uniaxial crystal with the unique axis along crystal z. Symmetry and, where physics
//! demands it, positive definiteness are checked with plain messages.

use super::elastic::sym_eigenvalues;
use super::reader::{MaterialError, Section, Tree};
use super::units::Dim;
use super::M3;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum TensorProp {
    Scalar(f64),
    Diagonal([f64; 3]),
    /// `perp` in the crystal xy plane, `par` along the crystal z axis
    Uniaxial { perp: f64, par: f64 },
    /// xx, yy, zz, yz, xz, xy
    Symmetric([f64; 6]),
}

/// What a physical tensor must satisfy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Requirement {
    /// every principal value above zero (permittivity, thermal conductivity)
    PositiveDefinite,
    /// principal values not below zero (electrical conductivity)
    NonNegative,
    /// only symmetry (thermal expansion may be negative)
    SymmetricOnly,
}

impl TensorProp {
    pub fn matrix(&self) -> M3 {
        match *self {
            TensorProp::Scalar(v) => [[v, 0.0, 0.0], [0.0, v, 0.0], [0.0, 0.0, v]],
            TensorProp::Diagonal([x, y, z]) => [[x, 0.0, 0.0], [0.0, y, 0.0], [0.0, 0.0, z]],
            TensorProp::Uniaxial { perp, par } => [[perp, 0.0, 0.0], [0.0, perp, 0.0], [0.0, 0.0, par]],
            TensorProp::Symmetric([xx, yy, zz, yz, xz, xy]) => [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]],
        }
    }

    pub fn principal_values(&self) -> [f64; 3] {
        sym_eigenvalues(self.matrix())
    }

    /// The scalar if the tensor is isotropic.
    pub fn as_scalar(&self) -> Option<f64> {
        match *self {
            TensorProp::Scalar(v) => Some(v),
            _ => None,
        }
    }

    /// The tensor in the older `Anisotropy` form.
    pub fn to_legacy(&self) -> super::Tensor {
        match *self {
            TensorProp::Scalar(v) => super::Tensor::Diagonal { x: v, y: v, z: v },
            TensorProp::Diagonal([x, y, z]) => super::Tensor::Diagonal { x, y, z },
            TensorProp::Uniaxial { perp, par } => super::Tensor::Uniaxial { perp, par },
            TensorProp::Symmetric([xx, yy, zz, yz, xz, xy]) => super::Tensor::Symmetric { xx, yy, zz, xy, yz, xz },
        }
    }

    /// Checks the requirement; `label` names the property in the message.
    pub fn check(&self, req: Requirement, label: &str) -> Result<(), String> {
        let m = self.matrix();
        if m.iter().flatten().any(|v| !v.is_finite()) {
            return Err(format!("the {label} contains a value that is not a finite number"));
        }
        let e = self.principal_values();
        let tiny = 1e-12 * e[2].abs().max(e[0].abs()).max(f64::MIN_POSITIVE);
        match req {
            Requirement::PositiveDefinite if e[0] <= tiny => Err(format!(
                "the {label} must be positive in every direction, but its smallest principal value is {}. Check the signs and the number of entries (1 = isotropic, 3 = x y z, 6 = xx yy zz yz xz xy)",
                e[0]
            )),
            Requirement::NonNegative if e[0] < -tiny => Err(format!("the {label} must not be negative in any direction, but its smallest principal value is {}", e[0])),
            _ => Ok(()),
        }
    }

    pub fn from_tree(s: &Section, key: &str, tree: &Tree) -> Result<TensorProp, MaterialError> {
        let num = |t: &Tree| match t {
            Tree::Num(x) => Some(*x),
            _ => None,
        };
        match tree {
            Tree::Num(x) => Ok(TensorProp::Scalar(*x)),
            Tree::Uniaxial(perp, par) => Ok(TensorProp::Uniaxial { perp: *perp, par: *par }),
            Tree::List(l) => {
                let flat: Option<Vec<f64>> = l.iter().map(num).collect();
                match (flat, l.len()) {
                    (Some(v), 3) => Ok(TensorProp::Diagonal([v[0], v[1], v[2]])),
                    (Some(v), 6) => Ok(TensorProp::Symmetric([v[0], v[1], v[2], v[3], v[4], v[5]])),
                    (Some(v), n) => Err(s.err(key, format!("has {n} entries; give 1 value (isotropic), 3 (x y z), 6 (xx yy zz yz xz xy), a 3x3 array, or {{ perp = .., par = .. }} (got {:?})", v))),
                    (None, 3) => {
                        let mut m = [[0.0; 3]; 3];
                        for (i, row) in l.iter().enumerate() {
                            let Tree::List(r) = row else { return Err(s.err(key, "a 3x3 array needs three rows of three numbers")) };
                            let r: Option<Vec<f64>> = r.iter().map(num).collect();
                            match r {
                                Some(r) if r.len() == 3 => m[i] = [r[0], r[1], r[2]],
                                _ => return Err(s.err(key, "a 3x3 array needs three rows of three numbers")),
                            }
                        }
                        let sc = m.iter().flatten().fold(0.0_f64, |a, v| a.max(v.abs())).max(f64::MIN_POSITIVE);
                        for i in 0..3 {
                            for j in i + 1..3 {
                                if (m[i][j] - m[j][i]).abs() > 1e-9 * sc {
                                    return Err(s.err(key, format!("the 3x3 array is not symmetric: entry ({},{}) = {} but ({},{}) = {}; the tensor must be symmetric", i + 1, j + 1, m[i][j], j + 1, i + 1, m[j][i])));
                                }
                            }
                        }
                        Ok(TensorProp::Symmetric([m[0][0], m[1][1], m[2][2], m[1][2], m[0][2], m[0][1]]))
                    }
                    _ => Err(s.err(key, "must be numbers, or three rows of three numbers")),
                }
            }
        }
    }

    /// Reads `key` of the section, converts to SI, checks the requirement.
    pub fn read(s: &mut Section, key: &str, dim: Dim, req: Requirement, label: &str) -> Result<Option<TensorProp>, MaterialError> {
        let Some(tree) = s.quantities(key, dim, None)? else { return Ok(None) };
        let t = Self::from_tree(s, key, &tree)?;
        t.check(req, label).map_err(|m| s.err(key, m))?;
        Ok(Some(t))
    }

    /// The value as TOML, numbers divided by `unit` (SI factor of the unit written next to it).
    pub fn write(&self, unit: f64) -> toml::Value {
        let f = |x: f64| toml::Value::Float(x / unit);
        match *self {
            TensorProp::Scalar(v) => f(v),
            TensorProp::Diagonal(v) => toml::Value::Array(v.iter().map(|x| f(*x)).collect()),
            TensorProp::Symmetric(v) => toml::Value::Array(v.iter().map(|x| f(*x)).collect()),
            TensorProp::Uniaxial { perp, par } => {
                let mut t = toml::Table::new();
                t.insert("perp".into(), f(perp));
                t.insert("par".into(), f(par));
                toml::Value::Table(t)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn negative_permittivity_is_refused_in_words() {
        let t = TensorProp::Diagonal([2.0, -1.0, 3.0]);
        let e = t.check(Requirement::PositiveDefinite, "permittivity").unwrap_err();
        assert!(e.contains("permittivity") && e.contains("positive"), "{e}");
        assert!(t.check(Requirement::SymmetricOnly, "expansion").is_ok());
    }

    #[test]
    fn symmetric_six_entry_order_is_voigt() {
        let t = TensorProp::Symmetric([1.0, 2.0, 3.0, 0.4, 0.5, 0.6]);
        let m = t.matrix();
        assert_eq!(m[1][2], 0.4);
        assert_eq!(m[0][2], 0.5);
        assert_eq!(m[0][1], 0.6);
    }
}
