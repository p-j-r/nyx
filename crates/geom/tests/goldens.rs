//! Golden models: the JSON feature lists in `golden/` evaluate to the regions and closed-form areas in `*.expected.json`.
//! `GEOM_WRITE_GOLDEN=1 cargo test --test goldens` rewrites the files from `tests/common`.

mod common;
use common::*;
use geom::eval::Tier;
use geom::model::FeatureList;
use std::path::PathBuf;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden")
}

#[test]
fn write_goldens_when_asked() {
    if std::env::var("GEOM_WRITE_GOLDEN").is_err() {
        return;
    }
    for g in all_models() {
        let ev = g.model.eval(Tier::Complex).unwrap();
        let mut e = g.expected.clone();
        if e.names.is_empty() {
            let mut n: Vec<String> = if g.model.list.dim == 3 {
                ev.faces().iter().map(|x| x.name.clone()).collect()
            } else {
                ev.edges().iter().map(|x| x.name.clone()).collect()
            };
            n.sort();
            n.dedup();
            e.names = n;
        }
        std::fs::write(
            dir().join(format!("{}.json", g.name)),
            serde_json::to_string_pretty(&g.model.list).unwrap() + "\n",
        )
        .unwrap();
        std::fs::write(
            dir().join(format!("{}.expected.json", g.name)),
            serde_json::to_string_pretty(&e).unwrap() + "\n",
        )
        .unwrap();
    }
}

#[test]
fn golden_models_match_their_invariants() {
    for g in all_models() {
        let list: FeatureList =
            serde_json::from_str(&std::fs::read_to_string(dir().join(format!("{}.json", g.name))).unwrap()).unwrap();
        assert_eq!(list, g.model.list, "{}: file differs from the builder", g.name);
        let exp: Expected =
            serde_json::from_str(&std::fs::read_to_string(dir().join(format!("{}.expected.json", g.name))).unwrap())
                .unwrap();
        let ev = geom::eval::Evaluator::new().eval(&list, Tier::Complex).unwrap();
        let cx = ev.complex();
        cx.check().unwrap();
        let names: Vec<&str> = cx.regions.iter().map(|r| r.name.as_str()).collect();
        let want: Vec<&str> = exp.regions.iter().map(|r| r.0.as_str()).collect();
        assert_eq!(names, want, "{}: regions", g.name);
        for (r, (n, a)) in cx.regions.iter().zip(&exp.regions) {
            assert!(
                (r.measure - a).abs() <= 1e-10 * a,
                "{}: region {n}: area {} vs closed form {a}",
                g.name,
                r.measure
            );
        }
        let mut en: Vec<String> = (0..cx.facet_count() as u32).map(|i| cx.facet_name(i).to_string()).collect();
        en.sort();
        en.dedup();
        assert_eq!(en, exp.names, "{}: boundary names", g.name);
        assert!(ev.warnings.is_empty(), "{}: {:?}", g.name, ev.warnings);
    }
}
