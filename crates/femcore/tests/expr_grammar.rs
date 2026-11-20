// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! `geom::expr` and `femcore::expr` are two copies of one grammar (ROUGH_SPOTS R06): geom must not depend on femcore, so the evaluator is
//! copied rather than shared. This test keeps the copies equal: on one corpus both accept the same inputs, both refuse the same inputs,
//! and the values agree to the last bit. (`femcore` depends on `geom`, so the test lives here.)

const X: f64 = 0.37;
const Y: f64 = -1.25;
const Z: f64 = 2.5;

fn femcore_value(text: &str) -> Result<f64, String> {
    femcore::expr::Expr::parse(text).map(|e| e.eval(X, Y, Z))
}

fn geom_value(text: &str) -> Result<f64, String> {
    let env = |n: &str| match n {
        "x" => Some(X),
        "y" => Some(Y),
        "z" => Some(Z),
        _ => None,
    };
    geom::expr::eval_text(text, &env)
}

/// Inputs both grammars accept: numbers, the four operations and powers, unary minus, parentheses, every function, the constants.
const ACCEPTED: &[&str] = &[
    "1",
    "1.5e-7",
    ".5",
    "2 + 3 * 4",
    "(2 + 3) * 4",
    "2 ^ 3 ^ 2",
    "-x",
    "--x",
    "-(x + y)",
    "x - y - z",
    "x / y / z",
    "x * y + z",
    "2 ** 3",
    "pi",
    "e",
    "x + pi * y - e * z",
    "sin(x) + cos(y) + tan(z)",
    "exp(-((x - 0.5) ^ 2) / 0.01)",
    "log(2.5) + ln(2.5) + log10(2.5)",
    "sqrt(z) + abs(y)",
    "sinh(x) + cosh(y) + tanh(z)",
    "atan(x * y)",
    "min(x, y) + max(y, z)",
    "step(x) + step(-x)",
    "293.15 + 50*exp(-((x-0.5)^2)/0.01)",
    "  sin( x )  ",
    "1e3 * x",
    "2^-1",
];

/// Inputs both grammars refuse: unbalanced or empty input, unknown names and functions, wrong argument counts, stray characters.
const REFUSED: &[&str] = &[
    "",
    "   ",
    "1 +",
    "* 2",
    "(1 + 2",
    "1 + 2)",
    "foo",
    "foo(1)",
    "sin()",
    "sin(1, 2)",
    "min(1)",
    "min(1, 2, 3)",
    "1 2",
    "x y",
    "1 $ 2",
    "system(\"ls\")",
    "x = 1",
    "1..2",
    "sin x",
];

#[test]
fn both_grammars_accept_the_same_inputs_with_the_same_values() {
    for text in ACCEPTED {
        let a = femcore_value(text).unwrap_or_else(|e| panic!("femcore refuses '{text}': {e}"));
        let b = geom_value(text).unwrap_or_else(|e| panic!("geom refuses '{text}': {e}"));
        assert!(a == b || (a.is_nan() && b.is_nan()), "'{text}': femcore {a:?}, geom {b:?}");
    }
}

#[test]
fn both_grammars_refuse_the_same_inputs() {
    for text in REFUSED {
        assert!(femcore_value(text).is_err(), "femcore accepts '{text}'");
        assert!(geom_value(text).is_err(), "geom accepts '{text}'");
    }
}
