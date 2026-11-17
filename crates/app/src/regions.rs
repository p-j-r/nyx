//! Which colour a region gets on the canvas (ADR 88, plan U4, section 5 of docs/ux/UX_MERGED_PLAN.md). The colour follows the entity, not the rank:
//! a region is keyed by its material when it has one (regions of the same material share a colour), otherwise by its own number, and slots are
//! handed out in order of first appearance over the region numbers in ascending order, so hiding, filtering or reordering never repaints. Slots 1
//! to 8 are the eight validated solid colours (`tokens.rs`), slots 9 to 16 the same eight hues again with a hatch (a texture of thin lines: the
//! one cue that does not depend on colour), slot 17 and up a neutral grey with only the number chip. A ninth hue is never generated.

use crate::tokens::{self, Mode};
use eframe::egui::Color32;
use std::collections::{BTreeMap, HashMap};

/// How a slot is drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Style {
    /// one of the eight validated colours (index 0 to 7)
    Solid(usize),
    /// the same eight hues again, with 1 px diagonal lines over the fill
    Hatched(usize),
    /// slot 17 and up: neutral, told apart by the number chip only
    Grey,
}

/// The style of slot `n` (0-based). Without a hatch pass (3D surfaces) slots 9 and up are grey, so a colour is never used twice without a cue.
pub fn style_of(n: usize, hatch: bool) -> Style {
    match n {
        0..=7 => Style::Solid(n),
        8..=15 if hatch => Style::Hatched(n - 8),
        _ => Style::Grey,
    }
}

/// The fill of a style on the canvas of `canvas` (opaque).
pub fn fill(style: Style, canvas: Mode) -> Color32 {
    match style {
        Style::Solid(i) | Style::Hatched(i) => crate::theme::rgb(tokens::regions(canvas)[i]),
        Style::Grey => crate::theme::rgb(tokens::REGION_GREY.get(canvas)),
    }
}

/// One legend row: the regions that share a slot.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    /// 0-based slot
    pub slot: usize,
    /// the material (or the region's name when it has none)
    pub name: String,
    pub labels: Vec<i32>,
}

/// The colour plan of one model: a slot per region number.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Paint {
    pub slot_of: HashMap<i32, usize>,
    /// slot order, the first appearance first
    pub rows: Vec<Row>,
    /// whether slots 9 to 16 may be hatched (2D)
    pub hatch: bool,
    /// changes whenever a colour of any region would
    pub hash: u64,
}

impl Paint {
    pub fn style(&self, label: i32) -> Style {
        style_of(self.slot_of.get(&label).copied().unwrap_or(usize::MAX), self.hatch)
    }
    pub fn slot(&self, label: i32) -> Option<usize> {
        self.slot_of.get(&label).copied()
    }
}

/// Hands out the slots: `labels` in any order, `key_of` the colour key of a label (material name or the region's own name).
pub fn assign(labels: &[i32], key_of: impl Fn(i32) -> String, hatch: bool) -> Paint {
    let mut sorted: Vec<i32> = labels.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut slot_of = HashMap::new();
    let mut rows: Vec<Row> = vec![];
    let mut slot_by_key: HashMap<String, usize> = HashMap::new();
    for l in sorted {
        let key = key_of(l);
        let slot = *slot_by_key.entry(key.clone()).or_insert_with(|| {
            rows.push(Row { slot: rows.len(), name: key, labels: vec![] });
            rows.len() - 1
        });
        rows[slot].labels.push(l);
        slot_of.insert(l, slot);
    }
    let mut hash = 0xcbf2_9ce4_8422_2325u64 ^ u64::from(hatch);
    for r in &rows {
        for &l in &r.labels {
            hash = (hash ^ (l as u32 as u64 | (r.slot as u64) << 32)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    Paint { slot_of, rows, hatch, hash }
}

/// The key of every region number by its material: `regions` are (number, name), `materials` the region patterns of the active problem with their
/// material (`None` for properties given as numbers: the region keeps a colour of its own). A later pattern overrides an earlier one, as in the
/// solver (`materials_by_label`); `matches` says which numbers a pattern names.
pub fn keys_by_material(regions: &[(i32, String)], materials: &BTreeMap<String, Option<String>>, matches: impl Fn(&str) -> Vec<i32>) -> HashMap<i32, String> {
    let mut key: HashMap<i32, String> = regions.iter().map(|(l, n)| (*l, n.clone())).collect();
    for (pattern, material) in materials {
        if let Some(m) = material {
            for l in matches(pattern) {
                key.insert(l, m.clone());
            }
        }
    }
    key
}

#[cfg(test)]
mod tests {
    use super::*;

    fn by_id(l: i32) -> String {
        format!("region {l}")
    }

    #[test]
    fn slots_follow_first_appearance_in_ascending_label_order() {
        let p = assign(&[7, 3, 5], by_id, true);
        assert_eq!((p.slot(3), p.slot(5), p.slot(7)), (Some(0), Some(1), Some(2)));
        // the same regions in another order get the same slots: nothing repaints when a list is reordered
        assert_eq!(assign(&[5, 7, 3], by_id, true), p);
    }

    #[test]
    fn regions_of_one_material_share_a_slot() {
        let key = |l: i32| if l == 2 || l == 4 { "silicon".to_string() } else { by_id(l) };
        let p = assign(&[1, 2, 3, 4], key, true);
        assert_eq!(p.slot(2), p.slot(4));
        assert_ne!(p.slot(1), p.slot(2));
        assert_eq!(p.rows.len(), 3);
        assert_eq!(p.rows[1].labels, vec![2, 4]);
    }

    #[test]
    fn eight_solid_then_hatched_then_grey() {
        let labels: Vec<i32> = (1..=20).collect();
        let p = assign(&labels, by_id, true);
        for l in 1..=8 {
            assert_eq!(p.style(l), Style::Solid(l as usize - 1));
        }
        for l in 9..=16 {
            assert_eq!(p.style(l), Style::Hatched(l as usize - 9));
        }
        for l in 17..=20 {
            assert_eq!(p.style(l), Style::Grey);
        }
        // a surface without a hatch pass never repeats a hue: slot 9 and up are grey
        let p3 = assign(&labels, by_id, false);
        assert_eq!(p3.style(9), Style::Grey);
        assert_eq!(p3.style(8), Style::Solid(7));
    }

    #[test]
    fn a_hatched_slot_is_the_same_hue_as_its_solid_twin_and_grey_is_neither() {
        for mode in [Mode::Dark, Mode::Light] {
            for i in 0..8 {
                assert_eq!(fill(Style::Solid(i), mode), fill(Style::Hatched(i), mode));
                assert_ne!(fill(Style::Grey, mode), fill(Style::Solid(i), mode));
            }
        }
    }

    #[test]
    fn material_keys_override_region_names_and_later_patterns_win() {
        let regions = vec![(1, "core".to_string()), (2, "clad".to_string()), (3, "box".to_string())];
        let mut m: BTreeMap<String, Option<String>> = BTreeMap::new();
        m.insert("clad".into(), Some("silica".into()));
        m.insert("box".into(), Some("silica".into()));
        m.insert("core".into(), None); // properties as numbers: keeps its own name
        let matches = |p: &str| regions.iter().filter(|(_, n)| n == p).map(|(l, _)| *l).collect::<Vec<_>>();
        let key = keys_by_material(&regions, &m, matches);
        assert_eq!(key[&1], "core");
        assert_eq!(key[&2], "silica");
        assert_eq!(key[&3], "silica");
        // the paint of this model: core, then silica twice
        let p = assign(&[1, 2, 3], |l| key[&l].clone(), true);
        assert_eq!(p.rows.len(), 2);
        assert_eq!(p.slot(2), p.slot(3));
    }

    #[test]
    fn the_hash_changes_when_a_colour_would() {
        let a = assign(&[1, 2, 3], by_id, true);
        let b = assign(&[1, 2, 3], |l| if l == 3 { "x".into() } else { by_id(l) }, true);
        let c = assign(&[1, 2, 3], |_| "all".into(), true);
        assert_ne!(a.hash, c.hash);
        assert_eq!(a.hash, b.hash, "same slots, same colours");
        assert_ne!(a.hash, assign(&[1, 2, 3], by_id, false).hash);
    }
}
