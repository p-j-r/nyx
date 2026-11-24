// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Which colour a region gets on the canvas (ADR 88, plan U4, section 5 of docs/ux/UX_MERGED_PLAN.md). The colour follows the entity, not the rank:
//! a region is keyed by its material when it has one (regions of the same material share a colour), otherwise by its own number, and slots are
//! handed out in order of first appearance over the region numbers in ascending order, so hiding, filtering or reordering never repaints. Slots 1
//! to 8 are the eight validated solid colours (`tokens.rs`), slots 9 to 16 the same eight hues again with a hatch (a texture of thin lines: the
//! one cue that does not depend on colour), slot 17 and up a neutral grey with only the number chip. A ninth hue is never generated.
//! Absorbing layers (PML) are a role, not a material (ADR 123): every region that any physics of the project marks as a PML shares one row, drawn
//! in one grey graded by the strength of its layer (the lighter, the weaker) and hatched, whichever physics is on show.

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
    /// an absorbing layer: the neutral grey graded by strength (0 weakest to 255 strongest), hatched
    Pml(u8),
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
        Style::Pml(level) => {
            // from 55 percent grey (weakest) to full grey (strongest) over the canvas: opaque, and the order of the grades is the order of lightness
            let (g, bg) = (tokens::REGION_GREY.get(canvas), tokens::CANVAS.get(canvas));
            let t = 0.55 + 0.45 * f32::from(level) / 255.0;
            let m = |a: u8, b: u8| (f32::from(a) * t + f32::from(b) * (1.0 - t)).round() as u8;
            crate::theme::rgb([m(g[0], bg[0]), m(g[1], bg[1]), m(g[2], bg[2])])
        }
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
    /// how the row is drawn (its swatch in the legend)
    pub style: Style,
}

/// The colour plan of one model: a slot per region number.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Paint {
    pub slot_of: HashMap<i32, usize>,
    /// the absorbing layers: the grade (0 to 255) of every region that is one; they share the last row
    pub pml: HashMap<i32, u8>,
    /// slot order, the first appearance first
    pub rows: Vec<Row>,
    /// whether slots 9 to 16 may be hatched (2D)
    pub hatch: bool,
    /// changes whenever a colour of any region would
    pub hash: u64,
}

impl Paint {
    pub fn style(&self, label: i32) -> Style {
        match self.pml.get(&label) {
            Some(&level) => Style::Pml(level),
            None => style_of(self.slot_of.get(&label).copied().unwrap_or(usize::MAX), self.hatch),
        }
    }
    pub fn slot(&self, label: i32) -> Option<usize> {
        self.slot_of.get(&label).copied()
    }
}

/// Hands out the slots: `labels` in any order, `key_of` the colour key of a label (material name or the region's own name).
pub fn assign(labels: &[i32], key_of: impl Fn(i32) -> String, hatch: bool) -> Paint {
    assign_with_pml(labels, key_of, hatch, &BTreeMap::new())
}

/// The same, with the absorbing layers: `pml` holds the strength of the layer of every region number that is one (the largest, when several
/// physics mark it). Those regions take no slot of a material: they share the last row, and the grade of each is its strength over the largest.
pub fn assign_with_pml(labels: &[i32], key_of: impl Fn(i32) -> String, hatch: bool, pml: &BTreeMap<i32, f64>) -> Paint {
    let mut sorted: Vec<i32> = labels.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut slot_of = HashMap::new();
    let mut rows: Vec<Row> = vec![];
    let mut slot_by_key: HashMap<String, usize> = HashMap::new();
    let mut layers: Vec<i32> = vec![];
    for l in sorted {
        if pml.contains_key(&l) {
            layers.push(l);
            continue;
        }
        let key = key_of(l);
        let slot = *slot_by_key.entry(key.clone()).or_insert_with(|| {
            let slot = rows.len();
            rows.push(Row { slot, name: key, labels: vec![], style: style_of(slot, hatch) });
            slot
        });
        rows[slot].labels.push(l);
        slot_of.insert(l, slot);
    }
    let mut grades = HashMap::new();
    if !layers.is_empty() {
        let top = layers.iter().map(|l| pml[l]).fold(0.0f64, f64::max);
        let slot = rows.len();
        for &l in &layers {
            // the grade is the strength over the strongest layer of the model; without a positive strength every layer is the full grey
            let grade = if top > 0.0 { (pml[&l] / top).clamp(0.0, 1.0) } else { 1.0 };
            grades.insert(l, (grade * 255.0).round() as u8);
            slot_of.insert(l, slot);
        }
        rows.push(Row { slot, name: "PML".into(), labels: layers, style: Style::Pml(255) });
    }
    let mut hash = 0xcbf2_9ce4_8422_2325u64 ^ u64::from(hatch);
    for r in &rows {
        for &l in &r.labels {
            hash = (hash ^ (l as u32 as u64 | (r.slot as u64) << 32 | u64::from(grades.get(&l).copied().unwrap_or(0)) << 48)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    Paint { slot_of, pml: grades, rows, hatch, hash }
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

    #[test]
    fn absorbing_layers_share_one_row_graded_by_strength_and_take_no_material_slot() {
        // air, slab, two layers of strength 4 and 2: the layers are one row "PML" after the materials, and the materials keep slots 0 and 1
        let pml: BTreeMap<i32, f64> = [(3, 4.0), (4, 2.0)].into();
        let p = assign_with_pml(&[1, 2, 3, 4], by_id, true, &pml);
        assert_eq!((p.slot(1), p.slot(2)), (Some(0), Some(1)));
        assert_eq!(p.rows.len(), 3);
        assert_eq!((p.rows[2].name.as_str(), p.rows[2].labels.clone()), ("PML", vec![3, 4]));
        assert_eq!(p.slot(3), p.slot(4));
        assert_eq!((p.style(3), p.style(4)), (Style::Pml(255), Style::Pml(128)));
        // one hue: the same grey, graded; the strongest layer is the full grey, a weaker one is nearer to the canvas
        for mode in [Mode::Dark, Mode::Light] {
            let (strong, weak) = (fill(p.style(3), mode), fill(p.style(4), mode));
            assert_eq!(strong, crate::theme::rgb(tokens::REGION_GREY.get(mode)));
            assert_ne!(strong, weak);
            assert_eq!(strong.a(), 255);
            assert_eq!(weak.a(), 255);
        }
    }

    #[test]
    fn the_plan_of_the_materials_does_not_move_when_a_layer_appears_or_the_physics_changes() {
        // the same materials with and without layers: the slots of the materials are the same, so switching the physics shown never repaints them
        let key = |l: i32| format!("m{}", l.min(2));
        let without = assign_with_pml(&[1, 2], key, true, &BTreeMap::new());
        let with = assign_with_pml(&[1, 2, 7], key, true, &[(7, 3.0)].into());
        assert_eq!((without.slot(1), without.slot(2)), (with.slot(1), with.slot(2)));
        assert_eq!(with.rows.len(), 3);
        // a layer of no strength is still a layer, at the full grey
        let flat = assign_with_pml(&[1, 7], by_id, true, &[(7, 0.0)].into());
        assert_eq!(flat.style(7), Style::Pml(255));
        assert_ne!(with.hash, without.hash);
    }
}
