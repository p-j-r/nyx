//! `LabelMap`: the compiled, deterministic mapping between names and the integer labels the solvers consume. Domain label =
//! region (1-based in overlay order, pins honoured and checked); boundary label = one per facet (ordered by provenance, then
//! position), except pinned groups, whose fragments share the pinned integer. Ints are output; users address names and rules.

use crate::error::{GeomError, GeomResult};
use crate::eval::Evaluated;
use crate::model::{Op, RegionKind};
use crate::select::{EntityKind, Selector};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DomainLabel {
    pub id: i32,
    pub name: String,
    pub region: u32,
    pub kind: RegionKind,
    pub pinned: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BoundaryLabel {
    pub id: i32,
    pub name: String,
    /// facets (2D: edges) carrying this label
    pub entities: Vec<u32>,
    pub pinned: bool,
}

/// A named group of boundary labels and the rule that produced it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Group {
    pub name: String,
    pub rule: String,
    pub labels: Vec<i32>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LabelMap {
    pub domains: Vec<DomainLabel>,
    pub boundaries: Vec<BoundaryLabel>,
    pub groups: Vec<Group>,
}

impl LabelMap {
    /// All boundary labels whose name is `name` (a face split into fragments has several).
    pub fn boundary(&self, name: &str) -> Vec<i32> {
        let mut v: Vec<i32> = self.boundaries.iter().filter(|b| b.name == name).map(|b| b.id).collect();
        for g in self.groups.iter().filter(|g| g.name == name) {
            v.extend(&g.labels);
        }
        v.sort_unstable();
        v.dedup();
        v
    }
    pub fn name_of_boundary(&self, id: i32) -> Option<&str> {
        self.boundaries.iter().find(|b| b.id == id).map(|b| b.name.as_str())
    }
    pub fn domain(&self, name: &str) -> Option<i32> {
        self.domains.iter().find(|d| d.name == name).map(|d| d.id)
    }
    pub fn name_of_domain(&self, id: i32) -> Option<&str> {
        self.domains.iter().find(|d| d.id == id).map(|d| d.name.as_str())
    }
    pub fn edges_of_label(&self, id: i32) -> Vec<u32> {
        self.boundaries.iter().filter(|b| b.id == id).flat_map(|b| b.entities.clone()).collect()
    }
    pub fn label_of_edge(&self, e: u32) -> Option<i32> {
        self.boundaries.iter().find(|b| b.entities.contains(&e)).map(|b| b.id)
    }
}

/// Smallest positive integer not in `used`, starting the search at `from`.
fn free(used: &[i32], from: &mut i32) -> i32 {
    while used.contains(from) {
        *from += 1;
    }
    let v = *from;
    *from += 1;
    v
}

pub fn compile(ev: &Evaluated) -> GeomResult<LabelMap> {
    let cx = ev.complex();
    // domains
    let mut used: Vec<i32> = vec![];
    for r in &cx.regions {
        if let Some(p) = r.pinned_id {
            if p <= 0 {
                return Err(GeomError(format!(
                    "region '{}' is pinned to the label {p}; domain labels must be positive",
                    r.name
                )));
            }
            if let Some(o) = cx.regions.iter().find(|o| o.pinned_id == Some(p) && o.name != r.name) {
                return Err(GeomError(format!(
                    "regions '{}' and '{}' are both pinned to the domain label {p}",
                    r.name, o.name
                )));
            }
            used.push(p);
        }
    }
    let mut next = 1;
    let domains: Vec<DomainLabel> = cx
        .regions
        .iter()
        .enumerate()
        .map(|(i, r)| DomainLabel {
            id: r.pinned_id.unwrap_or_else(|| free(&used, &mut next)),
            name: r.name.clone(),
            region: i as u32,
            kind: r.kind,
            pinned: r.pinned_id.is_some(),
        })
        .collect();

    // boundaries: pins first
    let n = cx.facet_count();
    let mut label_of: Vec<Option<(i32, String, bool)>> = vec![None; n];
    let mut bused: Vec<i32> = vec![];
    for pin in &ev.scene().pins {
        if pin.id <= 0 {
            return Err(GeomError(format!(
                "the pin of '{}' uses the boundary label {}; labels must be positive",
                pin.sel, pin.id
            )));
        }
        let r = ev.resolve_inner(&pin.sel, true)?;
        if r.set.kind == EntityKind::Regions || r.set.kind == EntityKind::Vertices {
            return Err(GeomError(format!(
                "the pin of '{}' selects {:?}; only boundary facets can be pinned",
                pin.sel, r.set.kind
            )));
        }
        if r.set.ids.is_empty() {
            return Err(GeomError(format!("the pin of '{}' to label {} matches nothing", pin.sel, pin.id)));
        }
        for &e in &r.set.ids {
            if let Some((o, _, _)) = &label_of[e as usize] {
                if *o != pin.id {
                    return Err(GeomError(format!(
                        "facet '{}' is selected by pins of both {} and {}",
                        cx.facet_name(e),
                        o,
                        pin.id
                    )));
                }
            }
            label_of[e as usize] = Some((pin.id, pin.sel.to_string(), true));
        }
        bused.push(pin.id);
    }
    let mut nb = 1;
    for (e, slot) in label_of.iter_mut().enumerate() {
        if slot.is_none() {
            *slot = Some((free(&bused, &mut nb), cx.facet_name(e as u32).to_string(), false));
        }
    }
    let mut boundaries: Vec<BoundaryLabel> = vec![];
    for (e, l) in label_of.into_iter().enumerate() {
        let (id, name, pinned) = l.unwrap();
        match boundaries.iter_mut().find(|b| b.id == id) {
            Some(b) => b.entities.push(e as u32),
            None => boundaries.push(BoundaryLabel { id, name, entities: vec![e as u32], pinned }),
        }
    }

    // groups: one per auto-name, one per tag
    let mut groups: Vec<Group> = vec![];
    let mut names: Vec<&str> = (0..n as u32).map(|e| cx.facet_name(e)).collect();
    names.sort_unstable();
    names.dedup();
    let labels_of = |ids: &[u32]| -> Vec<i32> {
        let mut v: Vec<i32> =
            boundaries.iter().filter(|b| b.entities.iter().any(|e| ids.contains(e))).map(|b| b.id).collect();
        v.sort_unstable();
        v
    };
    for nm in names {
        let ids: Vec<u32> = (0..n as u32).filter(|&e| cx.facet_name(e) == nm).collect();
        groups.push(Group {
            name: nm.to_string(),
            rule: Selector::Name(nm.to_string()).to_string(),
            labels: labels_of(&ids),
        });
    }
    let mut tags: Vec<String> = ev.list.features.iter().flat_map(|f| f.tags.iter().cloned()).collect();
    tags.extend(ev.list.features.iter().filter_map(|f| {
        if let Op::Tag { name, .. } = &f.op {
            Some(name.clone())
        } else {
            None
        }
    }));
    tags.extend(ev.scene().auto_tags.iter().map(|t| t.2.clone()));
    tags.sort();
    tags.dedup();
    for t in tags {
        let sel = Selector::Tag(t.clone());
        if let Ok(r) = ev.resolve_inner(&sel, true) {
            if r.set.kind != EntityKind::Regions && r.set.kind != EntityKind::Vertices {
                groups.push(Group { name: t, rule: sel.to_string(), labels: labels_of(&r.set.ids) });
            }
        }
    }
    Ok(LabelMap { domains, boundaries, groups })
}

impl Evaluated {
    pub fn label_map(&self) -> GeomResult<LabelMap> {
        compile(self)
    }

    /// Integer labels (boundary labels for facets, domain labels for regions) a selector resolves to.
    pub fn select_labels(&self, sel: &Selector) -> GeomResult<Vec<i32>> {
        let lm = compile(self)?;
        let r = self.resolve(sel)?;
        let mut v: Vec<i32> = match r.set.kind {
            EntityKind::Regions => r.set.ids.iter().map(|&i| lm.domains[i as usize].id).collect(),
            EntityKind::Vertices => {
                return Err(GeomError(format!("'{sel}' selects vertices, which have no integer label")))
            }
            _ => r.set.ids.iter().filter_map(|&e| lm.label_of_edge(e)).collect(),
        };
        v.sort_unstable();
        v.dedup();
        Ok(v)
    }
}
