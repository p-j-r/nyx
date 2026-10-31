//! Floquet eigenmodes against an exact band edge. A stack of two dielectrics seen by a plane wave at normal incidence is the textbook
//! Bragg mirror: with the cell periodic in x and y (phase 0) and in z (phase `pi`, the band edge of the stack) the eigenmodes of a 3D cell
//! are the plane waves of the stack, and the Bloch relation
//! `cos(K L) = cos p1 cos p2 - (n1/n2 + n2/n1)/2 sin p1 sin p2`, `p_i = n_i k0 d_i`, gives the frequencies at which `K L = pi`. The edges
//! of the first gap are what the Bragg grating use case measures, here on a structure whose answer is known. Modes come in pairs
//! (two polarizations).
use femcore::mesh::SimplexMesh;
use femcore::problem::*;

const C0: f64 = 299_792_458.0;

fn bloch_cos(f: f64, n: [f64; 2], d: [f64; 2]) -> f64 {
    let k0 = 2.0 * std::f64::consts::PI * f / C0;
    let (p1, p2) = (n[0] * k0 * d[0], n[1] * k0 * d[1]);
    p1.cos() * p2.cos() - 0.5 * (n[0] / n[1] + n[1] / n[0]) * p1.sin() * p2.sin()
}

/// Edges of the first gap at `K L = pi`: the first frequencies at which `cos(K L)` reaches -1, found by a scan and bisection.
fn first_gap(n: [f64; 2], d: [f64; 2]) -> (f64, f64) {
    let g = |f: f64| bloch_cos(f, n, d) + 1.0;
    let l = d[0] + d[1];
    let df = C0 / (2.0 * l * n[0].max(n[1])) * 1e-4;
    let mut f = df;
    let mut edges = vec![];
    let mut prev = g(f);
    while edges.len() < 2 {
        let next = g(f + df);
        if (prev > 0.0 && next <= 0.0 && edges.is_empty()) || (prev <= 0.0 && next > 0.0 && edges.len() == 1) {
            let (mut lo, mut hi) = (f, f + df);
            let s = g(lo) > 0.0;
            for _ in 0..80 {
                let m = 0.5 * (lo + hi);
                if (g(m) > 0.0) == s {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            edges.push(0.5 * (lo + hi));
        }
        prev = next;
        f += df;
    }
    (edges[0], edges[1])
}

#[test]
fn the_floquet_cell_reproduces_the_band_edges_of_a_dielectric_stack() {
    let (n, d) = ([12.0f64.sqrt(), 2.0f64.sqrt()], [0.0075, 0.0125]);
    let (a, l) = (0.01, d[0] + d[1]);
    let (f_lo, f_hi) = first_gap(n, d);
    assert!(f_hi > 1.1 * f_lo, "there is a gap: {f_lo} {f_hi}");
    let mut errs = vec![];
    for nz in [20, 40] {
        let mut mesh = SimplexMesh::block(a, a, l, 2, 2, nz);
        mesh.relabel(&[(vec![0.0, 0.0, 0.0, a, a, d[0]], 1), (vec![0.0, 0.0, d[0], a, a, l], 2)], &[]).unwrap();
        let mut s = ProblemSpec::new(Physics::Emw);
        s.order = 2;
        for (dom, nn) in [(1, n[0]), (2, n[1])] {
            s.materials.set(dom, "eps_r", nn * nn);
            s.materials.set(dom, "mu_r", 1.0);
        }
        s.periodic = vec![
            Periodic { source: vec![1], dest: vec![6], phase: 0.0 },
            Periodic { source: vec![2], dest: vec![5], phase: 0.0 },
            Periodic { source: vec![3], dest: vec![4], phase: std::f64::consts::PI },
        ];
        s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.97 * f_lo };
        let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
        let mut f: Vec<f64> = r.solutions.iter().map(|x| x.value.unwrap()).collect();
        f.sort_by(|x, y| x.total_cmp(y));
        // two polarizations: the lower edge twice, then the upper edge twice
        let mut worst = 0.0f64;
        for (got, exact) in [(f[0], f_lo), (f[1], f_lo), (f[2], f_hi), (f[3], f_hi)] {
            worst = worst.max((got / exact - 1.0).abs());
            assert!((got / exact - 1.0).abs() < 0.015, "nz = {nz}: {got:e} against {exact:e} (modes {f:?}, exact edges {f_lo:e} and {f_hi:e})");
        }
        errs.push(worst);
        // the gap itself (a small difference of two frequencies) is right to better than 2 %
        assert!(((f[2] - f[0]) / (f_hi - f_lo) - 1.0).abs() < 0.02, "gap {} against {}", f[2] - f[0], f_hi - f_lo);
    }
    assert!(errs[1] < 0.6 * errs[0], "the edges converge with refinement: {errs:?}");
}
