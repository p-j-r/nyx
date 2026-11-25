// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Property-based (`proptest`) checks for `ProblemSpec` canonical JSON round-trips and exhaustive `check_result_finite` fault injection.

use nyx_physics::result::{CellField, Field, Scalar, Solution, SolveResult};
use nyx_physics::{
    canonical_json, check_result_finite, BcSpec, MaterialAssign, ProblemSpec, PropValue, StudySpec, PROBLEM_FORMAT,
};
use proptest::prelude::*;
use serde_json::{Map, Number, Value};

fn arb_ident() -> impl Strategy<Value = String> {
    "[a-zA-Z][a-zA-Z0-9_-]{0,15}"
}

fn arb_finite_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        Just(0.0_f64),
        (-1e6_f64..1e6_f64).prop_filter("non-zero finite", |x| x.is_finite() && *x != 0.0),
        prop_oneof![Just(1e-12), Just(-1e-12), Just(std::f64::consts::PI), Just(1e9)],
    ]
    .prop_map(|x| if x == 0.0 { 0.0 } else { x })
}

fn arb_json_value() -> impl Strategy<Value = Value> {
    let leaf = prop_oneof![
        Just(Value::Null),
        any::<bool>().prop_map(Value::Bool),
        (-10_000i64..10_000).prop_map(|n| Value::Number(Number::from(n))),
        arb_finite_f64().prop_map(|f| Value::Number(Number::from_f64(f).expect("finite f64"))),
        "[a-zA-Z0-9 _.:/-]{0,24}".prop_map(Value::String),
    ];
    leaf.prop_recursive(3, 24, 6, |inner| {
        prop_oneof![
            prop::collection::vec(inner.clone(), 0..5).prop_map(Value::Array),
            prop::collection::vec((arb_ident(), inner), 0..5).prop_map(|entries| {
                let mut m = Map::new();
                for (k, v) in entries {
                    m.insert(k, v);
                }
                Value::Object(m)
            }),
        ]
    })
}

fn arb_json_object_or_null() -> impl Strategy<Value = Value> {
    prop_oneof![
        Just(Value::Null),
        prop::collection::vec((arb_ident(), arb_json_value()), 0..6).prop_map(|entries| {
            let mut m = Map::new();
            for (k, v) in entries {
                m.insert(k, v);
            }
            Value::Object(m)
        }),
    ]
}

fn arb_prop_value_tensor() -> impl Strategy<Value = PropValue> {
    let f = arb_finite_f64;
    prop_oneof![
        proptest::array::uniform3(f()).prop_map(PropValue::Diag),
        proptest::array::uniform6(f()).prop_map(PropValue::Sym6),
        proptest::array::uniform3(proptest::array::uniform3(f())).prop_map(PropValue::Matrix),
        proptest::array::uniform3(proptest::array::uniform6(f())).prop_map(PropValue::Matrix3x6),
    ]
}

fn arb_prop_value_any() -> impl Strategy<Value = PropValue> {
    prop_oneof![arb_finite_f64().prop_map(PropValue::Scalar), arb_prop_value_tensor(),]
}

fn arb_material_assign() -> impl Strategy<Value = MaterialAssign> {
    prop_oneof![
        arb_ident().prop_map(MaterialAssign::Named),
        (arb_ident(), proptest::array::uniform3(arb_finite_f64()))
            .prop_map(|(material, euler_deg)| MaterialAssign::Oriented { material, euler_deg }),
        prop::collection::btree_map(
            arb_ident().prop_filter("avoid oriented key", |k| k != "material"),
            arb_finite_f64(),
            1..5
        )
        .prop_map(MaterialAssign::Props),
        (
            arb_ident().prop_filter("avoid oriented key", |k| k != "material"),
            arb_prop_value_tensor(),
            prop::collection::btree_map(
                arb_ident().prop_filter("avoid oriented key", |k| k != "material"),
                arb_prop_value_any(),
                0..4,
            ),
        )
            .prop_map(|(k0, v0, mut rest)| {
                rest.insert(k0, v0);
                MaterialAssign::PropsTensor(rest)
            }),
    ]
}

fn arb_problem_spec() -> impl Strategy<Value = ProblemSpec> {
    (
        prop_oneof![Just("dfs".to_string()), Just("sms".to_string()), Just("ems".to_string()), arb_ident()],
        1u8..=4,
        (arb_ident(), arb_json_object_or_null()).prop_map(|(kind, params)| StudySpec { kind, params }),
        prop::collection::btree_map(arb_ident(), arb_material_assign(), 0..5),
        prop::collection::vec(
            (arb_ident(), arb_ident(), arb_json_object_or_null()).prop_map(|(on, kind, params)| BcSpec {
                on,
                kind,
                params,
            }),
            0..5,
        ),
        arb_json_object_or_null(),
    )
        .prop_map(|(physics, order, study, materials, bcs, params)| ProblemSpec {
            format: PROBLEM_FORMAT,
            physics,
            order,
            study,
            materials,
            bcs,
            params,
        })
}

fn shuffle_json_keys(v: &Value, mut state: u64) -> Value {
    match v {
        Value::Object(map) => {
            let mut pairs: Vec<(String, Value)> = map
                .iter()
                .map(|(k, val)| {
                    (k.clone(), shuffle_json_keys(val, state.wrapping_mul(6364136223846793005).wrapping_add(1)))
                })
                .collect();
            let n = pairs.len();
            if n > 1 {
                for i in (1..n).rev() {
                    state = state.wrapping_mul(2862933555777941757).wrapping_add(3037000493);
                    let j = (state as usize) % (i + 1);
                    pairs.swap(i, j);
                }
            }
            let mut out = Map::new();
            for (k, val) in pairs {
                out.insert(k, val);
            }
            Value::Object(out)
        }
        Value::Array(arr) => {
            Value::Array(arr.iter().map(|item| shuffle_json_keys(item, state ^ 0x9e37_79b9_7f4a_7c15)).collect())
        }
        other => other.clone(),
    }
}

fn arb_solution() -> impl Strategy<Value = Solution> {
    (
        arb_ident(),
        prop::collection::vec(
            (arb_ident(), arb_ident(), prop::collection::vec(arb_finite_f64(), 1..6))
                .prop_map(|(name, unit, values)| Field::scalar(name, unit, values)),
            1..4,
        ),
        prop::option::of(arb_finite_f64()),
        prop::option::of(proptest::array::uniform2(arb_finite_f64())),
        prop::collection::vec(
            (arb_ident(), arb_ident(), prop::collection::vec(arb_finite_f64(), 1..6), any::<bool>()).prop_flat_map(
                |(name, unit, re, with_im)| {
                    let n = re.len();
                    let im_strat = if with_im {
                        prop::collection::vec(arb_finite_f64(), n..=n).boxed()
                    } else {
                        Just(Vec::new()).boxed()
                    };
                    (Just(name), Just(unit), Just(re), im_strat).prop_map(|(name, unit, re, im)| CellField {
                        name,
                        unit,
                        nloc: 1,
                        re,
                        im,
                    })
                },
            ),
            1..4,
        ),
    )
        .prop_map(|(label, fields, value, eigenvalue, cell_fields)| Solution {
            label,
            fields,
            value,
            eigenvalue,
            cell_fields,
        })
}

fn arb_clean_solve_result() -> impl Strategy<Value = SolveResult> {
    (
        arb_ident(),
        arb_ident(),
        1u8..=3,
        2u8..=3,
        1usize..100,
        arb_finite_f64(),
        prop::collection::vec(arb_solution(), 1..4),
        prop::collection::vec(
            (
                arb_ident(),
                prop_oneof![
                    arb_finite_f64().prop_map(Scalar::Num),
                    Just(Scalar::Num(f64::INFINITY)),
                    Just(Scalar::Num(f64::NEG_INFINITY)),
                    arb_ident().prop_map(Scalar::Text),
                ],
            ),
            1..6,
        ),
        prop::collection::vec(arb_ident(), 0..4),
    )
        .prop_map(|(physics, study, order, dim, dofs, seconds, solutions, scalars, log)| SolveResult {
            physics,
            study,
            order,
            dim,
            dofs,
            seconds,
            solutions,
            scalars,
            log,
        })
}

proptest! {
    #[test]
    fn problem_spec_canonical_json_roundtrip_and_key_order_invariance(
        spec in arb_problem_spec(),
        shuffle_seed in any::<u64>(),
    ) {
        let val = serde_json::to_value(&spec).expect("serialize ProblemSpec");
        let canon = canonical_json(&val);
        prop_assert_eq!(canonical_json(&canon), canon.clone(), "canonical_json must be idempotent");

        let shuffled = shuffle_json_keys(&val, shuffle_seed);
        let canon_shuffled = canonical_json(&shuffled);
        prop_assert_eq!(&canon_shuffled, &canon, "canonical_json must ignore key insertion order");
        prop_assert_eq!(
            serde_json::to_string(&canon_shuffled).unwrap(),
            serde_json::to_string(&canon).unwrap()
        );

        let wire = serde_json::to_string(&canon).unwrap();
        let decoded: ProblemSpec = serde_json::from_str(&wire).expect("deserialize ProblemSpec");
        prop_assert_eq!(&decoded, &spec);

        let re_val = serde_json::to_value(&decoded).unwrap();
        prop_assert_eq!(canonical_json(&re_val), canon);
    }

    #[test]
    fn problem_spec_scalar_props_tensor_preserves_canonical_json(
        entries in prop::collection::btree_map(
            arb_ident().prop_filter("avoid oriented key", |k| k != "material"),
            arb_prop_value_any(),
            0..5,
        ),
        seed in any::<u64>(),
    ) {
        let mut spec = ProblemSpec::new("ems", 2, "eigen");
        spec.materials.insert("core".into(), MaterialAssign::PropsTensor(entries));
        let val = serde_json::to_value(&spec).unwrap();
        let canon = canonical_json(&val);
        let shuffled = canonical_json(&shuffle_json_keys(&val, seed));
        prop_assert_eq!(&shuffled, &canon);

        let roundtripped: ProblemSpec = serde_json::from_value(canon.clone()).unwrap();
        let rt_canon = canonical_json(&serde_json::to_value(&roundtripped).unwrap());
        prop_assert_eq!(rt_canon, canon);
    }

    #[test]
    fn check_result_finite_exhaustive_fault_injection(
        clean in arb_clean_solve_result(),
        sol_sel in any::<prop::sample::Index>(),
        field_sel in any::<prop::sample::Index>(),
        val_sel in any::<prop::sample::Index>(),
        cf_sel in any::<prop::sample::Index>(),
        cf_idx_sel in any::<prop::sample::Index>(),
        scalar_sel in any::<prop::sample::Index>(),
    ) {
        prop_assert!(check_result_finite(&clean).is_ok(), "clean finite result must pass");

        let non_finite = [
            f64::NAN,
            f64::from_bits(0x7ff8_0000_0000_0042),
            f64::INFINITY,
            f64::NEG_INFINITY,
        ];
        let si = sol_sel.index(clean.solutions.len());
        let fi = field_sel.index(clean.solutions[si].fields.len());
        let vi = val_sel.index(clean.solutions[si].fields[fi].values.len());
        let cfi = cf_sel.index(clean.solutions[si].cell_fields.len());
        let ri = cf_idx_sel.index(clean.solutions[si].cell_fields[cfi].re.len());

        for bad in non_finite {
            let mut r = clean.clone();
            r.solutions[si].value = Some(bad);
            prop_assert!(check_result_finite(&r).is_err(), "solution.value {bad} must fail");

            let mut r = clean.clone();
            r.solutions[si].eigenvalue = Some([bad, 1.0]);
            prop_assert!(check_result_finite(&r).is_err(), "eigenvalue[0] {bad} must fail");

            let mut r = clean.clone();
            r.solutions[si].eigenvalue = Some([1.0, bad]);
            prop_assert!(check_result_finite(&r).is_err(), "eigenvalue[1] {bad} must fail");

            let mut r = clean.clone();
            r.solutions[si].fields[fi].values[vi] = bad;
            prop_assert!(check_result_finite(&r).is_err(), "field value {bad} must fail");

            let mut r = clean.clone();
            r.solutions[si].cell_fields[cfi].re[ri] = bad;
            prop_assert!(check_result_finite(&r).is_err(), "cell_field re {bad} must fail");

            if !clean.solutions[si].cell_fields[cfi].im.is_empty() {
                let ii = cf_idx_sel.index(clean.solutions[si].cell_fields[cfi].im.len());
                let mut r = clean.clone();
                r.solutions[si].cell_fields[cfi].im[ii] = bad;
                prop_assert!(check_result_finite(&r).is_err(), "cell_field im {bad} must fail");
            }
        }

        let sc_idx = scalar_sel.index(clean.scalars.len());
        let mut r = clean.clone();
        r.scalars[sc_idx].1 = Scalar::Num(f64::NAN);
        prop_assert!(check_result_finite(&r).is_err(), "scalar NaN must fail");

        let mut r = clean.clone();
        r.scalars[sc_idx].1 = Scalar::Num(f64::INFINITY);
        prop_assert!(check_result_finite(&r).is_ok(), "scalar +inf represents lossless Q and must pass");
    }
}
