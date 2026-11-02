//! The cross-section of a waveguide in wafer axes: the one place where the axis numbers of the 2D kernels are written down.
//!
//! The wafer frame is x lateral, y the direction of propagation, z the surface normal (right-handed). A 2D mesh is the x-z plane of the
//! wafer, so mesh axis 0 is wafer x and mesh axis 1 is wafer z; y is the axis out of the plane. A kernel never writes a wafer axis number
//! itself: it asks this module for the wafer axis behind a mesh axis, for the in-plane block of a tensor, and for the wafer name of a
//! component.
//!
//! Two views of a wafer tensor are offered:
//!
//! * `plane_tensor`: axes (x, z, y), no signs. The leading 2x2 block is the x-z block that diffusion and in-plane elasticity use, entry
//!   `[2][2]` is the yy (out-of-plane) value and `[0][2]`, `[1][2]` are the couplings of x and z to y. For a 3D mesh it is the tensor itself.
//! * `local_tensor`: the axes (x, z, `PROP_SIGN` y) the waveguide mode solvers work in (lateral, normal, propagation; right-handed because
//!   `PROP_SIGN` is -1). A mode solver's own equations are written in these local axes, so the optical kernels read their tensors and
//!   name their fields through this module and nothing else.
//!
//! `PROP_SIGN` is the one open frame decision of the optical modes: the mode solvers carry the field as `exp(-j beta z_local)`, which with
//! `z_local = -y` is a wave along `-y` in wafer terms. Everything that depends on it is written in terms of `PROP_SIGN`, so flipping it is
//! a one-line change that the golden tests (`tests/golden_frames.rs`) will show.

use crate::material::M3;

/// Wafer axis (0 = x, 1 = y, 2 = z) of each axis of a mesh of dimension `dim`.
pub fn wafer_axes(dim: usize) -> &'static [usize] {
    if dim == 2 {
        &[0, 2]
    } else {
        &[0, 1, 2]
    }
}

/// Wafer name of each component of a vector field on a mesh of dimension `dim`: `["x", "z"]` on a cross-section.
pub fn component_names(dim: usize) -> &'static [&'static str] {
    if dim == 2 {
        &["x", "z"]
    } else {
        &["x", "y", "z"]
    }
}

/// The wafer coordinates `(x, y, z)` of a mesh point (`y = 0` on a cross-section).
pub fn wafer_point(p: &[f64], dim: usize) -> [f64; 3] {
    if dim == 2 {
        [p[0], 0.0, p[1]]
    } else {
        [p[0], p[1], p[2]]
    }
}

/// A wafer tensor in the axes (x, z, y) of a cross-section (the identity on a 3D mesh), see the module documentation.
pub fn plane_tensor(t: &M3, dim: usize) -> M3 {
    if dim != 2 {
        return *t;
    }
    permute(t, &[0, 2, 1], &[1.0, 1.0, 1.0])
}

/// Sign of the propagation axis of the mode solvers against wafer y: everything propagates along `+y`.
pub const PROP_SIGN: f64 = 1.0;

/// The propagation axis of the Brillouin kernels, which still work in the model frame (x, z, -y) until phase 5 of the kernel port.
pub const LEGACY_PROP_SIGN: f64 = -1.0;

/// The wafer axis behind each axis of the local frame (lateral, normal, propagation) and its sign.
pub const LOCAL_AXES: [(usize, f64); 3] = [(0, 1.0), (2, 1.0), (1, PROP_SIGN)];

/// A wafer tensor in the local axes (lateral, normal, propagation) of the mode solvers.
pub fn local_tensor(t: &M3) -> M3 {
    local_tensor_signed(t, PROP_SIGN)
}

/// As `local_tensor` for a propagation axis along `prop_sign * y`.
pub fn local_tensor_signed(t: &M3, prop_sign: f64) -> M3 {
    permute(t, &[0, 2, 1], &[1.0, 1.0, prop_sign])
}

/// Whether the local triple (lateral, normal, propagation) of a mode solver with propagation `prop_sign * y` is left-handed.
pub fn is_left_handed(prop_sign: f64) -> bool {
    prop_sign > 0.0
}

fn permute(t: &M3, axis: &[usize; 3], sign: &[f64; 3]) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = sign[i] * sign[j] * t[axis[i]][axis[j]];
        }
    }
    r
}

/// Wafer name and sign of a component of a field written in the local axes: `0` lateral -> `x`, `1` normal -> `z`, `2` propagation -> `y`.
pub fn local_component(i: usize) -> (&'static str, f64) {
    (["x", "z", "y"][i], LOCAL_AXES[i].1)
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: M3 = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]];

    #[test]
    fn plane_tensor_puts_x_z_first() {
        let p = plane_tensor(&T, 2);
        assert_eq!([p[0][0], p[0][1], p[1][1]], [T[0][0], T[0][2], T[2][2]]);
        assert_eq!([p[2][2], p[0][2], p[1][2]], [T[1][1], T[0][1], T[2][1]]);
        assert_eq!(plane_tensor(&T, 3), T);
    }

    #[test]
    fn local_tensor_follows_the_table() {
        let l = local_tensor(&T);
        // lateral, normal, propagation = x, z, PROP_SIGN y
        assert_eq!([l[0][0], l[1][1], l[2][2], l[0][1]], [T[0][0], T[2][2], T[1][1], T[0][2]]);
        assert_eq!([l[0][2], l[1][2]], [PROP_SIGN * T[0][1], PROP_SIGN * T[2][1]]);
        let b = local_tensor_signed(&T, LEGACY_PROP_SIGN);
        assert_eq!([b[0][2], b[1][2]], [-T[0][1], -T[2][1]]);
    }

    #[test]
    fn handedness_of_the_local_axes_is_known() {
        for sign in [PROP_SIGN, LEGACY_PROP_SIGN] {
            let a = [(0usize, 1.0), (2, 1.0), (1, sign)];
            let e = |i: usize| {
                let mut v = [0.0; 3];
                v[a[i].0] = a[i].1;
                v
            };
            let (x, z) = (e(0), e(1));
            let cross = [x[1] * z[2] - x[2] * z[1], x[2] * z[0] - x[0] * z[2], x[0] * z[1] - x[1] * z[0]];
            let right = cross == e(2);
            assert_eq!(is_left_handed(sign), !right, "sign {sign}");
        }
        assert_eq!(LOCAL_AXES[2].1, PROP_SIGN);
    }
}
