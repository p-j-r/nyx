//! No new colour literal in the window code (ADR 88, rule "no Color32 literal outside theme.rs, tokens.rs and colormap.rs"). The literals that were
//! there before the colour tokens are counted per file and may only go DOWN: when you remove one, lower the number here; a file that is not in the
//! table has none and may not get one. The table is the to-do list of plan steps U4 (data and canvas colours: `viewport.rs`, `app.rs`, `tabs.rs`,
//! `dashboard.rs`, `wfview.rs`) and later; `geomui.rs` and `scene.rs` hold feature colours of the geometry preview and test fixtures.

use std::collections::BTreeMap;

const PATTERNS: [&str; 15] = [
    "Color32::from_rgb", "Color32::from_rgba_unmultiplied", "Color32::from_rgba_premultiplied", "Color32::from_gray", "Color32::from_black_alpha", "Color32::from_white_alpha",
    "Color32::WHITE", "Color32::BLACK", "Color32::RED", "Color32::GREEN", "Color32::BLUE", "Color32::GRAY", "Color32::YELLOW", "Color32::LIGHT_GRAY", "Color32::DARK_GRAY",
];

/// Where colour literals are still allowed, and how many (counted 2026-11-17, after plan step U1). `TRANSPARENT` is not a colour.
const ALLOWED: [(&str, usize); 7] = [("app.rs", 13), ("dashboard.rs", 3), ("geomui.rs", 3), ("scene.rs", 3), ("tabs.rs", 20), ("viewport.rs", 18), ("wfview.rs", 3)];
const EXEMPT: [&str; 3] = ["theme.rs", "tokens.rs", "colormap.rs"];

/// Occurrences of `pat` that are not the start of a longer name (`from_rgb` must not count `from_rgba_unmultiplied`).
fn count(src: &str, pat: &str) -> usize {
    src.match_indices(pat).filter(|(i, _)| !src[i + pat.len()..].starts_with(|c: char| c.is_alphanumeric() || c == '_')).count()
}

#[test]
fn colour_literals_only_go_down() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let allowed: BTreeMap<&str, usize> = ALLOWED.into_iter().collect();
    let mut bad = vec![];
    for e in std::fs::read_dir(&dir).expect("src folder") {
        let p = e.expect("entry").path();
        let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
        if !name.ends_with(".rs") || EXEMPT.contains(&name.as_str()) {
            continue;
        }
        let src = std::fs::read_to_string(&p).expect("read source");
        let n: usize = PATTERNS.iter().map(|pat| count(&src, pat)).sum();
        let limit = allowed.get(name.as_str()).copied().unwrap_or(0);
        if n > limit {
            bad.push(format!("{name}: {n} colour literals, {limit} allowed: take the colour from theme.rs (a role of tokens.rs) instead"));
        }
        if n < limit {
            bad.push(format!("{name}: {n} colour literals, the table allows {limit}: lower the number in tests/colour_literals.rs"));
        }
    }
    assert!(bad.is_empty(), "{}", bad.join("\n"));
}
