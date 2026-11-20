// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Warm starts of the eigensolvers (productivity plan P5b, ADR 96): the previous solve's modes as the start block give the same modes in
//! fewer iterations when the problem moved a little; a solve with no context, or a context with nothing to start from, is bit for bit the
//! solve it always was; a start of the wrong size is a cold start.

use femcore::linalg::csr::Csr;
use femcore::linalg::driver::{eigs_sym_indefinite, EigOpts};
use femcore::warm::{self, Warm};

/// `-u'' + v(x) u` on a line of `n` points (a Schrodinger operator): the potential moves with `shift`, the modes with it.
fn operator(n: usize, shift: f64) -> Csr<f64> {
    let mut t = Vec::new();
    for i in 0..n {
        let x = i as f64 / n as f64;
        let v = 400.0 * (1.0 + shift) * (x - 0.5).powi(2) + 30.0 * shift * (7.0 * x).sin();
        t.push((i, i, 2.0 + v / (n * n) as f64 * 1.0e4));
        if i + 1 < n {
            t.push((i, i + 1, -1.0));
            t.push((i + 1, i, -1.0));
        }
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

fn identity(n: usize) -> Csr<f64> {
    Csr::from_triplets(n, n, &(0..n).map(|i| (i, i, 1.0)).collect::<Vec<_>>()).unwrap()
}

fn solve(shift: f64) -> femcore::linalg::driver::EigPairs {
    let n = 1500;
    eigs_sym_indefinite(&operator(n, shift), &identity(n), 0.002, 5, &EigOpts::default(), None).unwrap()
}

#[test]
fn a_solve_with_no_context_and_one_with_nothing_to_start_from_are_the_same_bits() {
    let plain = solve(0.0);
    let w = Warm::new(vec![]);
    let recorded = {
        let _g = warm::install(&w);
        solve(0.0)
    };
    assert_eq!(plain.values, recorded.values);
    assert_eq!(plain.vectors, recorded.vectors);
    assert_eq!(w.calls().len(), 1);
    assert!(!w.was_used(), "the first solve of a sweep is cold");
}

#[test]
fn the_previous_modes_as_a_start_give_the_same_modes_in_fewer_iterations() {
    let first = Warm::new(vec![]);
    {
        let _g = warm::install(&first);
        solve(0.0);
    }
    let cold_iterations = first.calls()[0].iterations;
    let second = Warm::new(first.output());
    let warm_values = {
        let _g = warm::install(&second);
        solve(0.02)
    };
    assert!(second.was_used(), "the second solve started from the first's modes");
    let warm_iterations = second.calls()[0].iterations;
    let cold_values = solve(0.02);
    assert_eq!(warm_values.values.len(), cold_values.values.len());
    for (a, b) in warm_values.values.iter().zip(&cold_values.values) {
        assert!((a - b).abs() <= 1e-8 * b.abs().max(1e-12), "{a} against {b}");
    }
    assert!(warm_iterations < cold_iterations, "{warm_iterations} restarts warm, {cold_iterations} cold");
}

#[test]
fn a_start_of_another_size_or_count_is_a_cold_start_and_not_an_error() {
    let small = Warm::new(vec![]);
    {
        let _g = warm::install(&small);
        eigs_sym_indefinite(&operator(900, 0.0), &identity(900), 0.002, 5, &EigOpts::default(), None).unwrap();
    }
    let w = Warm::new(small.output());
    let _g = warm::install(&w);
    let p = solve(0.0);
    assert!(!w.was_used(), "another number of unknowns");
    assert_eq!(p.values.len(), 5);
}
