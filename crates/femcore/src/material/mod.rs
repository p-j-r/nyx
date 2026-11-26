// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Anisotropic material data: crystal stiffness (cubic, hexagonal, trigonal, orthotropic),
//! symmetric tensors for conductivity, permittivity and thermal expansion, and rotation of the
//! crystal axes into the global frame.
//!
//! Orientation: `euler_deg = [phi1, Phi, phi2]` (Bunge, Z-X-Z). The crystal frame is rotated
//! into the global frame by `R = Rz(phi1) Rx(Phi) Rz(phi2)`, i.e. the crystal x axis points along
//! the first column of `R`. Stiffness matrices are returned in the solver's Voigt order
//! (xx, yy, zz, xy, yz, xz) with engineering shear strains.

use serde::{Deserialize, Serialize};

pub mod bcs;
pub mod coupled;
pub mod elastic;
pub mod import;
pub mod optics;
pub mod reader;
pub mod schema;
pub mod source;
pub mod spec;
pub mod state_t;
pub mod tensor;
pub mod units;

pub use bcs::{
    bcs_gap, bcs_surface_impedance, surface_impedance_from_params, BcsEval, FilmModel, Superconducting,
    BCS_GAP_RATIO, EV_TO_J, HBAR, K_B, MU0,
};
pub use coupled::{
    Dielectric, ElectroOptic, InterfaceLayer, Interfaces, Photoelastic, PiezoForm, Piezoelectric, Viscoelastic, EPSILON_0,
};
pub use elastic::{Class, Elastic, Unique};
pub use optics::{AxisModel, Dispersion, Optical, OpticsError};
pub use reader::MaterialError;
pub use source::Source;
pub use spec::{merge_over, Library, MaterialSpec};
pub use state_t::{AxisClamp, CheckedTempEval, ComponentClamp, CrystalSetting, MaterialKind, PointGroup, TempClampInfo, TempCurve, TempTensor, TEMP_CURVE_KEYS};
pub use tensor::TensorProp;

pub type M3 = [[f64; 3]; 3];
pub type M6 = [[f64; 6]; 6];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[allow(clippy::large_enum_variant)] // a stiffness is made once per material, not per cell
pub enum Stiffness {
    Cubic { c11: f64, c12: f64, c44: f64 },
    /// transversely isotropic about the crystal z axis
    Hexagonal { c11: f64, c12: f64, c13: f64, c33: f64, c44: f64 },
    /// classes 3m and -3m (sapphire, LiNbO3): c14 = -c24 = c56, c66 = (c11 - c12) / 2
    Trigonal { c11: f64, c12: f64, c13: f64, c14: f64, c33: f64, c44: f64 },
    Orthotropic { c11: f64, c22: f64, c33: f64, c12: f64, c13: f64, c23: f64, c44: f64, c55: f64, c66: f64 },
    /// any symmetric 6x6 matrix in the textbook Voigt order (xx, yy, zz, yz, xz, xy) and the crystal frame, in Pa (upper triangle is read)
    General(M6),
}

impl Stiffness {
    pub fn name(&self) -> &'static str {
        match self {
            Stiffness::Cubic { .. } => "Cubic",
            Stiffness::Hexagonal { .. } => "Hexagonal",
            Stiffness::Trigonal { .. } => "Trigonal",
            Stiffness::Orthotropic { .. } => "Orthotropic",
            Stiffness::General(_) => "General",
        }
    }
    /// Standard Voigt order (xx, yy, zz, yz, xz, xy) in the crystal frame.
    fn voigt(&self) -> M6 {
        let mut c = [[0.0; 6]; 6];
        let mut set = |i: usize, j: usize, v: f64| {
            c[i][j] = v;
            c[j][i] = v;
        };
        match *self {
            Stiffness::Cubic { c11, c12, c44 } => {
                for i in 0..3 {
                    set(i, i, c11);
                    set(i + 3, i + 3, c44);
                    for j in i + 1..3 {
                        set(i, j, c12);
                    }
                }
            }
            Stiffness::Hexagonal { c11, c12, c13, c33, c44 } => {
                set(0, 0, c11);
                set(1, 1, c11);
                set(2, 2, c33);
                set(0, 1, c12);
                set(0, 2, c13);
                set(1, 2, c13);
                set(3, 3, c44);
                set(4, 4, c44);
                set(5, 5, 0.5 * (c11 - c12));
            }
            Stiffness::Trigonal { c11, c12, c13, c14, c33, c44 } => {
                set(0, 0, c11);
                set(1, 1, c11);
                set(2, 2, c33);
                set(0, 1, c12);
                set(0, 2, c13);
                set(1, 2, c13);
                set(0, 3, c14);
                set(1, 3, -c14);
                set(4, 5, c14);
                set(3, 3, c44);
                set(4, 4, c44);
                set(5, 5, 0.5 * (c11 - c12));
            }
            Stiffness::General(m) => {
                for i in 0..6 {
                    for j in i..6 {
                        set(i, j, m[i][j]);
                    }
                }
            }
            Stiffness::Orthotropic { c11, c22, c33, c12, c13, c23, c44, c55, c66 } => {
                set(0, 0, c11);
                set(1, 1, c22);
                set(2, 2, c33);
                set(0, 1, c12);
                set(0, 2, c13);
                set(1, 2, c23);
                set(3, 3, c44);
                set(4, 4, c55);
                set(5, 5, c66);
            }
        }
        c
    }

    /// The constants as `(label, value)` for editing.
    pub fn constants_mut(&mut self) -> Vec<(&'static str, &mut f64)> {
        match self {
            Stiffness::Cubic { c11, c12, c44 } => vec![("c11", c11), ("c12", c12), ("c44", c44)],
            Stiffness::Hexagonal { c11, c12, c13, c33, c44 } => vec![("c11", c11), ("c12", c12), ("c13", c13), ("c33", c33), ("c44", c44)],
            Stiffness::Trigonal { c11, c12, c13, c14, c33, c44 } => vec![("c11", c11), ("c12", c12), ("c13", c13), ("c14", c14), ("c33", c33), ("c44", c44)],
            Stiffness::General(m) => {
                let mut v = Vec::new();
                for (i, row) in m.iter_mut().enumerate() {
                    for (j, x) in row.iter_mut().enumerate() {
                        if j >= i {
                            v.push((elastic::NAMES[i][j], x));
                        }
                    }
                }
                v
            }
            Stiffness::Orthotropic { c11, c22, c33, c12, c13, c23, c44, c55, c66 } => {
                vec![("c11", c11), ("c22", c22), ("c33", c33), ("c12", c12), ("c13", c13), ("c23", c23), ("c44", c44), ("c55", c55), ("c66", c66)]
            }
        }
    }
}

/// A symmetric second-rank tensor property in the crystal frame.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Tensor {
    /// `perp` in the crystal xy plane, `par` along the crystal z axis
    Uniaxial { perp: f64, par: f64 },
    Diagonal { x: f64, y: f64, z: f64 },
    Symmetric { xx: f64, yy: f64, zz: f64, xy: f64, yz: f64, xz: f64 },
}

impl Tensor {
    pub fn name(&self) -> &'static str {
        match self {
            Tensor::Uniaxial { .. } => "Uniaxial",
            Tensor::Diagonal { .. } => "Diagonal",
            Tensor::Symmetric { .. } => "Symmetric",
        }
    }
    pub fn matrix(&self) -> M3 {
        match *self {
            Tensor::Uniaxial { perp, par } => [[perp, 0.0, 0.0], [0.0, perp, 0.0], [0.0, 0.0, par]],
            Tensor::Diagonal { x, y, z } => [[x, 0.0, 0.0], [0.0, y, 0.0], [0.0, 0.0, z]],
            Tensor::Symmetric { xx, yy, zz, xy, yz, xz } => [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]],
        }
    }
    pub fn constants_mut(&mut self) -> Vec<(&'static str, &mut f64)> {
        match self {
            Tensor::Uniaxial { perp, par } => vec![("⊥", perp), ("∥ (z)", par)],
            Tensor::Diagonal { x, y, z } => vec![("x", x), ("y", y), ("z", z)],
            Tensor::Symmetric { xx, yy, zz, xy, yz, xz } => vec![("xx", xx), ("yy", yy), ("zz", zz), ("xy", xy), ("yz", yz), ("xz", xz)],
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Anisotropy {
    pub stiffness: Option<Stiffness>,
    pub conductivity: Option<Tensor>,
    pub permittivity: Option<Tensor>,
    pub expansion: Option<Tensor>,
    /// Bunge Z-X-Z angles of the crystal axes in degrees
    pub euler_deg: [f64; 3],
    /// Imaginary part magnitude of relative permittivity (microwave/optical dielectric loss in crystal frame).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub permittivity_im: Option<Tensor>,
    /// Phonon/structural viscosity tensor eta in Pa*s (textbook Voigt order in crystal frame).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<M6>,
    /// Piezoelectric stress-charge coupling e_iJ in C/m^2 (textbook Voigt in crystal frame).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub piezo_e: Option<[[f64; 6]; 3]>,
    /// Photoelastic Pockels tensor p_IJ (dimensionless, textbook Voigt in crystal frame).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub photoelastic: Option<M6>,
}

pub fn rotation(euler_deg: [f64; 3]) -> M3 {
    let [a, b, c] = euler_deg.map(|d| d.to_radians());
    let rz = |t: f64| {
        let (s, c) = t.sin_cos();
        [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
    };
    let rx = |t: f64| {
        let (s, c) = t.sin_cos();
        [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]
    };
    mul3(&mul3(&rz(a), &rx(b)), &rz(c))
}

/// Bunge angles (degrees) of a rotation matrix whose columns are the crystal axes in the global
/// frame: the inverse of `rotation`.
pub fn euler_from_rotation(m: &M3) -> [f64; 3] {
    let b = m[2][2].clamp(-1.0, 1.0).acos();
    if b.sin().abs() < 1e-9 {
        // the first and last rotation are about the same axis: put it all in the first
        let a = m[1][0].atan2(m[0][0]);
        return [a.to_degrees(), b.to_degrees(), 0.0];
    }
    let a = m[0][2].atan2(-m[1][2]);
    let c = m[2][0].atan2(m[2][1]);
    [a.to_degrees(), b.to_degrees(), c.to_degrees()]
}

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

pub fn rotate_tensor(t: &M3, r: &M3) -> M3 {
    let rt = [[r[0][0], r[1][0], r[2][0]], [r[0][1], r[1][1], r[2][1]], [r[0][2], r[1][2], r[2][2]]];
    mul3(&mul3(r, t), &rt)
}

const VP: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];
/// position in the solver's Voigt order (xx, yy, zz, xy, yz, xz) of each standard Voigt component (xx, yy, zz, yz, xz, xy)
const TO_SOLVER: [usize; 6] = [0, 1, 2, 4, 5, 3];

fn rotate_voigt_raw(c: &M6, r: &M3, to_solver: bool) -> M6 {
    let mut t = [[[[0.0; 3]; 3]; 3]; 3];
    for (i, &(a, b)) in VP.iter().enumerate() {
        for (j, &(cc, d)) in VP.iter().enumerate() {
            let v = c[i][j];
            for (p, q) in [(a, b), (b, a)] {
                for (s, u) in [(cc, d), (d, cc)] {
                    t[p][q][s][u] = v;
                }
            }
        }
    }
    let mut rot = [[[[0.0; 3]; 3]; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let mut s = 0.0;
                    for a in 0..3 {
                        for b in 0..3 {
                            for cc in 0..3 {
                                for d in 0..3 {
                                    s += r[i][a] * r[j][b] * r[k][cc] * r[l][d] * t[a][b][cc][d];
                                }
                            }
                        }
                    }
                    rot[i][j][k][l] = s;
                }
            }
        }
    }
    let mut out = [[0.0; 6]; 6];
    for (i, &(a, b)) in VP.iter().enumerate() {
        let oi = if to_solver { TO_SOLVER[i] } else { i };
        for (j, &(cc, d)) in VP.iter().enumerate() {
            let oj = if to_solver { TO_SOLVER[j] } else { j };
            out[oi][oj] = rot[a][b][cc][d];
        }
    }
    out
}

/// Rotates a stiffness given in standard Voigt order and returns it in the solver's Voigt order.
fn rotate_stiffness(c: &M6, r: &M3) -> M6 {
    rotate_voigt_raw(c, r, true)
}

/// Stiffness in the global frame (solver Voigt order).
pub fn stiffness_global(s: &Stiffness, euler_deg: [f64; 3]) -> M6 {
    rotate_stiffness(&s.voigt(), &rotation(euler_deg))
}

/// Rotates a symmetric `6x6` fourth-order tensor given in textbook Voigt order `(xx, yy, zz, yz, xz, xy)`
/// from crystal frame to global frame by Bunge Z-X-Z Euler angles (degrees), returning solver Voigt order `(xx, yy, zz, xy, yz, xz)`.
pub fn rotate_m6_global(m: &M6, euler_deg: [f64; 3]) -> M6 {
    rotate_stiffness(m, &rotation(euler_deg))
}

/// Rotates a symmetric `6x6` fourth-order tensor given in textbook Voigt order `(xx, yy, zz, yz, xz, xy)`
/// from crystal frame to global frame by Bunge Z-X-Z Euler angles (degrees), preserving textbook Voigt order.
pub fn rotate_m6_textbook(m: &M6, euler_deg: [f64; 3]) -> M6 {
    rotate_voigt_raw(m, &rotation(euler_deg), false)
}

/// Rotates a `3x6` piezoelectric stress-charge coupling matrix `e_iJ` (`C/m^2`) from the crystal frame
/// (textbook Voigt columns `xx, yy, zz, yz, xz, xy`) into the global frame in solver Voigt column order
/// `(xx, yy, zz, xy, yz, xz)`.
pub fn piezo_e_global(e_textbook: &[[f64; 6]; 3], euler_deg: [f64; 3]) -> [[f64; 6]; 3] {
    let pz = Piezoelectric { source: Vec::new(), form: PiezoForm::StressCharge, matrix: *e_textbook };
    let rot_textbook = pz.rotate(&rotation(euler_deg));
    let mut out = [[0.0; 6]; 3];
    for row in 0..3 {
        for (j, &s) in TO_SOLVER.iter().enumerate() {
            out[row][s] = rot_textbook[row][j];
        }
    }
    out
}

/// Isotropic stiffness (solver Voigt order).
pub fn isotropic_stiffness(e: f64, nu: f64) -> M6 {
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = lam;
        }
        c[i][i] += 2.0 * mu;
        c[3 + i][3 + i] = mu;
    }
    c
}

pub fn tensor_global(t: &Tensor, euler_deg: [f64; 3]) -> M3 {
    rotate_tensor(&t.matrix(), &rotation(euler_deg))
}

pub fn iso_tensor(v: f64) -> M3 {
    [[v, 0.0, 0.0], [0.0, v, 0.0], [0.0, 0.0, v]]
}

/// Positive definiteness check by Cholesky on the 6x6 matrix.
pub fn is_positive_definite(c: &M6) -> bool {
    let mut l = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let mut s = c[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if s <= 0.0 {
                    return false;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    true
}

/// Typical room-temperature crystal data for the material library (GPa converted to Pa).
pub fn sapphire() -> (Stiffness, Tensor, Tensor, Tensor, f64) {
    (
        Stiffness::Trigonal { c11: 496.8e9, c12: 163.6e9, c13: 110.9e9, c14: -23.5e9, c33: 498.1e9, c44: 147.4e9 },
        Tensor::Uniaxial { perp: 9.4, par: 11.6 },
        Tensor::Uniaxial { perp: 23.0, par: 25.2 },
        Tensor::Uniaxial { perp: 5.0e-6, par: 5.3e-6 },
        3980.0,
    )
}

pub fn silicon() -> Stiffness {
    Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_with_isotropic_constants_is_isotropic() {
        let (e, nu) = (200e9, 0.3);
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let c = stiffness_global(&Stiffness::Cubic { c11: lam + 2.0 * mu, c12: lam, c44: mu }, [17.0, 33.0, 71.0]);
        let iso = isotropic_stiffness(e, nu);
        for i in 0..6 {
            for j in 0..6 {
                assert!((c[i][j] - iso[i][j]).abs() < 1e-12 * iso[0][0], "{i}{j}: {} vs {}", c[i][j], iso[i][j]);
            }
        }
    }

    #[test]
    fn trigonal_is_invariant_under_threefold_rotation_about_c() {
        let (s, ..) = sapphire();
        let a = stiffness_global(&s, [0.0, 0.0, 0.0]);
        let b = stiffness_global(&s, [120.0, 0.0, 0.0]);
        let scale = a[0][0];
        for i in 0..6 {
            for j in 0..6 {
                assert!((a[i][j] - b[i][j]).abs() < 1e-9 * scale, "{i}{j}");
            }
        }
        // but not under a 60 degree rotation
        let c = stiffness_global(&s, [60.0, 0.0, 0.0]);
        let diff = (0..6).flat_map(|i| (0..6).map(move |j| (i, j))).map(|(i, j)| (a[i][j] - c[i][j]).abs()).fold(0.0, f64::max);
        assert!(diff > 1e9, "a 60 degree rotation must change the trigonal stiffness, diff {diff}");
        assert!(is_positive_definite(&a));
    }

    #[test]
    fn rotation_matrix_is_orthogonal() {
        let r = rotation([30.0, 50.0, 70.0]);
        let rt = [[r[0][0], r[1][0], r[2][0]], [r[0][1], r[1][1], r[2][1]], [r[0][2], r[1][2], r[2][2]]];
        let p = mul3(&r, &rt);
        for i in 0..3 {
            for j in 0..3 {
                assert!((p[i][j] - (i == j) as u8 as f64).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn tensor_rotation_swaps_axes() {
        let t = tensor_global(&Tensor::Diagonal { x: 1.0, y: 2.0, z: 3.0 }, [90.0, 0.0, 0.0]);
        assert!((t[0][0] - 2.0).abs() < 1e-14 && (t[1][1] - 1.0).abs() < 1e-14 && (t[2][2] - 3.0).abs() < 1e-14);
    }

    #[test]
    fn piezo_e_global_maps_voigt_and_preserves_symmetry() {
        // Wurtzite 6mm (AlN): textbook e31=-0.58, e32=-0.58, e33=1.55, e24=-0.48 (yz), e15=-0.48 (xz)
        let e_tb = [
            [0.0, 0.0, 0.0, 0.0, -0.48, 0.0],
            [0.0, 0.0, 0.0, -0.48, 0.0, 0.0],
            [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0],
        ];
        let e0 = piezo_e_global(&e_tb, [0.0, 0.0, 0.0]);
        // Solver Voigt cols: (xx, yy, zz, xy, yz, xz) => e15 at [0][5], e24 at [1][4]
        assert!((e0[0][5] - (-0.48)).abs() < 1e-15);
        assert!((e0[1][4] - (-0.48)).abs() < 1e-15);
        assert!((e0[2][0] - (-0.58)).abs() < 1e-15);
        assert!((e0[2][2] - 1.55).abs() < 1e-15);
        // Hexagonal 6mm is invariant under arbitrary rotation about z
        let e_rot = piezo_e_global(&e_tb, [37.5, 0.0, 0.0]);
        for i in 0..3 {
            for j in 0..6 {
                assert!((e0[i][j] - e_rot[i][j]).abs() < 1e-14, "{i}{j}: {} vs {}", e0[i][j], e_rot[i][j]);
            }
        }
    }
}
