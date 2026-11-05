//! Anisotropic materials: tensor permittivity in the waveguide modes and the rotation of a crystal stiffness. The solid mechanics tests (stiffness,
//! expansion) moved to `nyx-sms` with the kernel.
mod common;
use common::*;
use femcore::material::*;
use femcore::props::*;

#[test]
fn shear_components_of_a_rotated_stiffness_are_in_the_solver_order() {
    // solver Voigt order: xx, yy, zz, xy, yz, xz. An unrotated orthotropic material has c66 on xy, c44 on yz, c55 on xz.
    let st = Stiffness::Orthotropic { c11: 200e9, c22: 180e9, c33: 160e9, c12: 60e9, c13: 55e9, c23: 50e9, c44: 40e9, c55: 45e9, c66: 50e9 };
    let c = stiffness_global(&st, [0.0, 0.0, 0.0]);
    assert_eq!((c[3][3], c[4][4], c[5][5]), (50e9, 40e9, 45e9));
    // a rotated crystal against a direct fourth-rank transform C'_ijkl = R_ia R_jb R_kc R_ld C_abcd (columns of R = crystal axes)
    let st = Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 };
    let c0 = stiffness_global(&st, [0.0, 0.0, 0.0]);
    let vp = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    let mut t = [[[[0.0; 3]; 3]; 3]; 3];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(c, d)) in vp.iter().enumerate() {
            for (p, q) in [(a, b), (b, a)] {
                for (u, v) in [(c, d), (d, c)] {
                    t[p][q][u][v] = c0[i][j];
                }
            }
        }
    }
    for psi in [0.0, 30.0, 45.0, 90.0] {
        // a proper rotation about the surface normal followed by a tilt of 57.6 degrees (any rotation tests the tensor transformation)
        let (cp, sp, ct, st_) = ((psi as f64).to_radians().cos(), (psi as f64).to_radians().sin(), 57.6f64.to_radians().cos(), 57.6f64.to_radians().sin());
        let rz = [[cp, -sp, 0.0], [sp, cp, 0.0], [0.0, 0.0, 1.0]];
        let rx = [[1.0, 0.0, 0.0], [0.0, ct, -st_], [0.0, st_, ct]];
        let mut r = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = (0..3).map(|k| rz[i][k] * rx[k][j]).sum();
            }
        }
        let m = stiffness_global(&st, euler_from_rotation(&r));
        for (i, &(a, b)) in vp.iter().enumerate() {
            for (j, &(c, d)) in vp.iter().enumerate() {
                let mut want = 0.0;
                for p in 0..3 {
                    for q in 0..3 {
                        for u in 0..3 {
                            for v in 0..3 {
                                want += r[a][p] * r[b][q] * r[c][u] * r[d][v] * t[p][q][u][v];
                            }
                        }
                    }
                }
                assert!((m[i][j] - want).abs() < 1e-6 * 500e9, "psi {psi} ({i},{j}): {} against {want}", m[i][j]);
            }
        }
    }
}
