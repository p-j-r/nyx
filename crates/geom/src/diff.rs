//! Change detection between two evaluations of related feature lists, for parameter sweeps: per-feature content hashes (the operation
//! with its expressions evaluated), propagated through the input references, plus the topology signature.

use crate::eval::Evaluated;
use crate::math::Hasher64;
use crate::model::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChangeKind {
    /// nothing that matters changed
    None,
    /// the 2D overlay of a layer stack is unchanged and only z extents differ
    ZOnly,
    /// equal topology signature, moved geometry: a mesh can be morphed
    Morph,
    /// the topology changed: remesh
    Remesh,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChangeReport {
    /// features whose own content or whose inputs changed (features present in only one list included), in id order
    pub changed_features: Vec<FeatureId>,
    pub signature_equal: bool,
    pub kind: ChangeKind,
}

/// Content hash of every feature: the operation, its name, its suppression flag and the numeric value of each expression.
pub(crate) fn feature_hashes(list: &FeatureList, env: &crate::expr::Env) -> Vec<(FeatureId, u64)> {
    list.features
        .iter()
        .map(|f| {
            let mut h = Hasher64::new();
            h.str(&format!("{:?}", f.op)).str(&f.name).u64(f.suppressed as u64);
            for e in f.op.exprs() {
                // an expression that cannot be evaluated hashes as NaN; the evaluation itself reports the error
                h.f64(env.eval(e, "").unwrap_or(f64::NAN));
            }
            (f.id, h.finish())
        })
        .collect()
}

pub fn diff(old: &Evaluated, new: &Evaluated) -> ChangeReport {
    let (a, b) = (&old.feature_hashes, &new.feature_hashes);
    let mut changed: Vec<FeatureId> = vec![];
    let find = |v: &Vec<(FeatureId, u64)>, id: FeatureId| v.iter().find(|x| x.0 == id).map(|x| x.1);
    for f in &new.list.features {
        let own = find(a, f.id) != find(b, f.id);
        let inputs = f.op.inputs().iter().any(|i| changed.contains(i));
        if own || inputs {
            changed.push(f.id);
        }
    }
    for f in &old.list.features {
        if find(b, f.id).is_none() {
            changed.push(f.id);
        }
    }
    changed.sort_unstable();
    changed.dedup();
    let signature_equal = old.topology_signature() == new.topology_signature();
    let kind = if !signature_equal {
        ChangeKind::Remesh
    } else if changed.is_empty() {
        ChangeKind::None
    } else if old.overlay_geometry_hash().is_some() && old.overlay_geometry_hash() == new.overlay_geometry_hash() {
        ChangeKind::ZOnly
    } else {
        ChangeKind::Morph
    };
    ChangeReport { changed_features: changed, signature_equal, kind }
}
