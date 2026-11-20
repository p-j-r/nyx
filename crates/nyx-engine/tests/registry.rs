// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The engine builds its registry from cargo features; every physics it holds follows the interface rules.
use nyx_engine::{physics_for, registry};
use nyx_physics::*;

#[test]
fn the_default_build_has_dfs_sms_ems_in_alphabetical_order() {
    let r = registry();
    if cfg!(all(feature = "dfs", feature = "sms", feature = "ems")) {
        assert_eq!(r.ids(), vec!["dfs", "ems", "sms"]);
    }
    for id in r.ids() {
        let p = physics_for(&r, id).unwrap();
        assert_eq!(p.id(), id);
        assert!(registry::is_valid_code(id));
    }
    let e = physics_for(&r, "xyz").err().unwrap();
    assert_eq!(e.code(), "unsupported");
    assert!(e.to_string().contains("xyz"));
}

#[test]
fn every_registered_physics_describes_itself_completely() {
    let r = registry();
    for id in r.ids() {
        let p = r.get(id).unwrap();
        assert!(!p.name().is_empty());
        let studies = p.studies();
        assert!(!studies.is_empty(), "{id} has no study");
        let mut kinds: Vec<_> = studies.iter().map(|s| s.kind.clone()).collect();
        kinds.sort();
        kinds.dedup();
        assert_eq!(kinds.len(), studies.len(), "{id}: study kinds must be unique");
        for s in &studies {
            assert!(!p.fields(&s.kind).is_empty(), "{id}/{}: no fields", s.kind);
            // a study kind is lowercase words joined by underscores (it is a key in files)
            assert!(s.kind.bytes().all(|b| b.is_ascii_lowercase() || b == b'_'), "{id}: study kind '{}'", s.kind);
        }
        // every property a study needs is one the physics lists
        let props = p.properties();
        for pr in &props {
            for k in &pr.needed_by {
                assert!(kinds.contains(k), "{id}: property {} is needed by the unknown study '{k}'", pr.name);
            }
        }
        // exports are qualified with the physics code
        for l in p.exports().iter().chain(p.imports().iter()) {
            assert!(l.field.starts_with(&format!("{id}.")), "{id}: link field '{}' is not qualified", l.field);
        }
    }
}

#[test]
fn heat_then_solid_is_a_one_way_chain_the_registry_can_order() {
    let r = registry();
    if r.get("dfs").is_some() && r.get("sms").is_some() {
        // sms imports a temperature that dfs exports; the order must put dfs first even though "sms" > "dfs" is also alphabetical,
        // so test the opposite request order
        assert_eq!(r.solve_order(&["sms", "dfs"]).unwrap(), vec!["dfs", "sms"]);
    }
    if r.get("sms").is_some() && r.get("ems").is_some() {
        // no coupling between them: both are solved, in any order
        assert_eq!(r.solve_order(&["sms", "ems"]).unwrap().len(), 2);
    }
}

#[test]
fn an_engine_without_features_is_empty_not_broken() {
    // compiled with `--no-default-features` the registry is empty and lookups fail with the list of what the build has (nothing)
    if r_is_empty() {
        let e = physics_for(&registry(), "dfs").err().unwrap();
        assert!(e.to_string().contains("this build has no physics"));
    }
}

fn r_is_empty() -> bool {
    registry().is_empty()
}
