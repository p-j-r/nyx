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
//! * `local_tensor`: the axes (x, z, `PROP_SIGN` y) the waveguide mode solvers work in (lateral, normal, propagation). Everything
//!   propagates along `+y` (`PROP_SIGN = 1`), so the local triple is left-handed: the curl-curl equations of the optical kernels do not
//!   care, a Poynting vector written with the same component formulas keeps its sign, and only the sense of rotation of a polarization
//!   is read the other way (`OpticalSet::left_handed`). The Brillouin kernels run in the same axes: the optical fields, the elastic modes
//!   (displacement blocks x, z, y, stiffness in the Voigt order of `plane_stiffness`) and the photoelastic tensor all use (x, z, y).
//!
//! The mode solvers carry the field as `exp(-j beta y)`, the elastic kz kernel as `exp(-j k y)`; a positive wave number is a wave along `+y`.

use crate::material::{M3, M6};
use crate::mesh::SimplexMesh;
use crate::props::*;

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

/// For each entry of the Voigt order of the axes (x, z, y) of a cross-section, (xx, zz, yy, xz, zy, xy), the entry of the wafer Voigt order
/// (xx, yy, zz, xy, yz, xz) it is. The solid kernels work in the first order: their in-plane block is `[0, 1, 3]` (xx, zz, xz) and the
/// out-of-plane normal strain (plane strain: zero, plane stress: condensed out) is entry 2 (yy).
pub const PLANE_VOIGT: [usize; 6] = [0, 2, 1, 5, 4, 3];

/// A wafer stiffness (Voigt order xx, yy, zz, xy, yz, xz, engineering shear) in the Voigt order of the axes (x, z, y) of a cross-section;
/// the identity on a 3D mesh.
pub fn plane_stiffness(c: &M6, dim: usize) -> M6 {
    if dim != 2 {
        return *c;
    }
    let mut r = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            r[i][j] = c[PLANE_VOIGT[i]][PLANE_VOIGT[j]];
        }
    }
    r
}

/// Wafer name of stress component `idx` of the solid kernels (their Voigt order: wafer order in 3D, `PLANE_VOIGT` on a cross-section).
pub fn stress_name(dim: usize, idx: usize) -> &'static str {
    const WAFER: [&str; 6] = ["sxx", "syy", "szz", "sxy", "syz", "sxz"];
    WAFER[if dim == 2 { PLANE_VOIGT[idx] } else { idx }]
}

/// The directions of a perfectly matched layer along the axes of a mesh: (x, z) on a cross-section, (x, y, z) in 3D.
pub fn mesh_pml_dirs(p: &Pml, dim: usize) -> [PmlDir; 3] {
    if dim == 2 {
        [p.x, p.z, PmlDir::Off]
    } else {
        [p.x, p.y, p.z]
    }
}

/// Sign of the propagation axis of the mode solvers against wafer y: everything propagates along `+y`.
pub const PROP_SIGN: f64 = 1.0;

/// The wafer axis behind each axis of the local frame (lateral, normal, propagation) and its sign.
pub const LOCAL_AXES: [(usize, f64); 3] = [(0, 1.0), (2, 1.0), (1, PROP_SIGN)];

/// A wafer tensor in the local axes (lateral, normal, propagation) of the mode solvers.
pub fn local_tensor(t: &M3) -> M3 {
    local_tensor_signed(t, PROP_SIGN)
}

/// As `local_tensor` for a propagation axis along `prop_sign * y` (the reciprocity test runs a guide both ways).
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

/// The three axes of a 3D problem that play the roles of the cross-section: `(lateral, normal, propagation)` as indices of x, y, z,
/// and whether that triple is left-handed (then the sense of rotation of the field must be read the other way round).
///
/// The propagation axis is the axis a Floquet pair of the problem is periodic along (the first one found); without periodic
/// conditions it is `y`. The surface normal is `z` unless the problem is periodic along z, then `y`.
pub fn axes_of_3d_periodic(periodic: &[Periodic], mesh: &SimplexMesh) -> ([usize; 3], bool) {
    let mut long = 1;
    // a Floquet pair maps the facets of `source` onto those of `dest`: the periodic axis is the one the two groups differ in
    if let Some(p) = periodic.first() {
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

    const T: M3 = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]];

    #[test]
    fn plane_tensor_puts_x_z_first() {
        let p = plane_tensor(&T, 2);
        assert_eq!([p[0][0], p[0][1], p[1][1]], [T[0][0], T[0][2], T[2][2]]);
        assert_eq!([p[2][2], p[0][2], p[1][2]], [T[1][1], T[0][1], T[2][1]]);
        assert_eq!(plane_tensor(&T, 3), T);
    }

    #[test]
    fn plane_voigt_is_the_voigt_order_of_x_z_y() {
        // an isotropic-like pattern: c[i][j] labelled by the tensor indices of Voigt entries i and j
        let vp = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]; // wafer order
        let pl = [(0, 0), (2, 2), (1, 1), (0, 2), (2, 1), (0, 1)]; // plane order, as wafer axes
        for k in 0..6 {
            let (a, b) = pl[k];
            let w = PLANE_VOIGT[k];
            assert!(vp[w] == (a, b) || vp[w] == (b, a), "entry {k}");
        }
        assert_eq!(stress_name(2, 1), "szz");
        assert_eq!(stress_name(2, 2), "syy");
        assert_eq!(stress_name(2, 3), "sxz");
        assert_eq!(stress_name(3, 5), "sxz");
    }

    #[test]
    fn local_tensor_follows_the_table() {
        let l = local_tensor(&T);
        // lateral, normal, propagation = x, z, PROP_SIGN y
        assert_eq!([l[0][0], l[1][1], l[2][2], l[0][1]], [T[0][0], T[2][2], T[1][1], T[0][2]]);
        assert_eq!([l[0][2], l[1][2]], [PROP_SIGN * T[0][1], PROP_SIGN * T[2][1]]);
    }

    #[test]
    fn handedness_of_the_local_axes_is_known() {
        for sign in [PROP_SIGN, -PROP_SIGN] {
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
