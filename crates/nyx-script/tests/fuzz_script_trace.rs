// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified `bolero` fuzz and regression tests for `nyx_script::trace::guard`, `femcore::formula::parse`, and `nyx_script::kernel` limits.

use bolero::TypeGenerator;
use femcore::fieldexpr::Tree;
use femcore::formula::{parse as parse_formula, Named, Names};
use nyx_script::kernel::{parse_timeout_secs, push_bounded_output, MAX_TIMEOUT};
use nyx_script::trace;

struct FuzzNames;

impl Names for FuzzNames {
    fn name(&mut self, name: &str) -> Result<Option<Named>, String> {
        Ok(match name {
            "x" => Some(Named::Scalar(Tree::Leaf(0))),
            "y" => Some(Named::Scalar(Tree::Leaf(1))),
            "z" => Some(Named::Scalar(Tree::Leaf(2))),
            "V" => Some(Named::Scalar(Tree::Leaf(3))),
            "Ex" => Some(Named::Scalar(Tree::Leaf(4))),
            "Ey" => Some(Named::Scalar(Tree::Leaf(5))),
            "Ez" => Some(Named::Scalar(Tree::Leaf(6))),
            "E" => Some(Named::Vector(vec![
                (0, Tree::Leaf(4)),
                (1, Tree::Leaf(5)),
                (2, Tree::Leaf(6)),
            ])),
            _ => None,
        })
    }

    fn material(&mut self, name: &str) -> Result<Tree, String> {
        if name.trim().is_empty() {
            Err("empty material".into())
        } else {
            Ok(Tree::num(1.0))
        }
    }

    fn region(&mut self, name: &str) -> Result<Tree, String> {
        if name.trim().is_empty() {
            Err("empty region".into())
        } else {
            Ok(Tree::Leaf(7))
        }
    }
}

#[derive(Debug, TypeGenerator)]
struct KernelFuzzCase {
    timeout_bits: u64,
    cap: u16,
    chunk_splits: [u8; 4],
    payload: Vec<u8>,
}

#[test]
fn fuzz_trace_guard_and_formula_parse_no_panic() {
    let mut wrapped = String::new();
    bolero::check!().with_iterations(512).for_each(|bytes: &[u8]| {
        let text = String::from_utf8_lossy(bytes);
        let snippet: &str = match text.char_indices().nth(2048) {
            Some((idx, _)) => &text[..idx],
            None => &text,
        };

        let _ = trace::guard("fuzz.py", snippet);

        wrapped.clear();
        wrapped.push_str("def helper(x):\n    return x\n");
        wrapped.push_str("def target(x, y, z):\n    return (");
        wrapped.push_str(snippet);
        wrapped.push_str(")\nnyx.ems(eps_r=target)\n");
        let _ = trace::guard("fuzz_wrapped.py", &wrapped);

        let mut names = FuzzNames;
        if let Ok(tree) = parse_formula(snippet, &mut names) {
            assert!(tree.leaves() <= 8, "unexpected leaf index: {}", tree.leaves());
        }
    });
}

#[test]
fn fuzz_kernel_timeout_and_bounded_printer_no_panic() {
    bolero::check!()
        .with_type::<KernelFuzzCase>()
        .with_iterations(512)
        .for_each(|case: &KernelFuzzCase| {
            let secs = f64::from_bits(case.timeout_bits);
            let res = parse_timeout_secs(secs, "timeout_s");
            if !secs.is_finite() || secs < 0.0 {
                assert!(res.is_err(), "non-finite or negative timeout {secs} must fail");
            } else if let Ok(d) = res {
                assert!(d <= MAX_TIMEOUT, "clamped duration {d:?} exceeded MAX_TIMEOUT");
            }

            let text = String::from_utf8_lossy(&case.payload);
            let cap = usize::from(case.cap);
            let mut out = String::new();
            let total = text.len();
            let mut prev = 0usize;
            for &raw_pos in &case.chunk_splits {
                if total == 0 {
                    push_bounded_output(&mut out, cap, "");
                    continue;
                }
                let raw = (usize::from(raw_pos) * total / 256).clamp(prev, total);
                let mut end = raw;
                while end > prev && !text.is_char_boundary(end) {
                    end -= 1;
                }
                push_bounded_output(&mut out, cap, &text[prev..end]);
                assert!(out.len() <= cap, "output len {} exceeded cap {cap}", out.len());
                prev = end;
            }
            push_bounded_output(&mut out, cap, &text[prev..]);
            assert!(out.len() <= cap, "output len {} exceeded cap {cap}", out.len());
            assert!(std::str::from_utf8(out.as_bytes()).is_ok());
        });
}

#[test]
fn regression_trace_guard_rejects_over_500_rounds() {
    use std::fmt::Write as _;

    let mut src = String::from(
        "def h0(p):\n    return 1.0 if p > 0.0 else 0.0\n",
    );
    for i in 1..505usize {
        let _ = writeln!(src, "def h{i}(p):\n    return p + 1.0");
    }
    src.push_str("def f(x, y, z):\n    return ");
    for i in 0..505usize {
        if i > 0 {
            src.push_str(" + ");
        }
        let _ = write!(src, "h{i}(x)");
    }
    src.push_str("\nnyx.ems(eps_r=f)\n");

    let err = trace::guard("deep_helpers.py", &src)
        .expect_err("trace::guard must reject >500 helper call rounds with explicit Err");
    assert!(err.contains("500"), "unexpected error message: {err}");
}
