//! Editing the geometry of a project: every edit of a `geom` feature list the window or a host can make, as a function that either does it
//! or leaves the list untouched and says why. The operations and their arguments come from `geomdesc`, so a new operation in geom
//! is available here (and as a tool) without a change in this file. Features are named by name, never by id, in this interface.

use crate::geomdesc::{args_of, build_op, find, kind_of};
use geom::expr::Expr;
use geom::model::{check_name, Feature, FeatureList, Origin, Param};
use nyx_physics::Error;
use serde_json::{json, Value};

fn bad(e: geom::GeomError) -> Error {
    Error::Invalid(e.to_string())
}

/// Applies `edit` to a copy; the copy replaces the list only if it is still a valid feature list.
fn transact<R>(list: &mut FeatureList, edit: impl FnOnce(&mut FeatureList) -> Result<R, Error>) -> Result<R, Error> {
    let mut copy = list.clone();
    let r = edit(&mut copy)?;
    copy.validate().map_err(bad)?;
    *list = copy;
    Ok(r)
}

fn index_of(list: &FeatureList, name: &str) -> Result<usize, Error> {
    list.features.iter().position(|f| f.name == name).ok_or_else(|| {
        Error::invalid(format!(
            "there is no feature named '{name}' (the features are: {})",
            list.features.iter().map(|f| f.name.as_str()).collect::<Vec<_>>().join(", ")
        ))
    })
}

fn next_id(list: &FeatureList) -> u32 {
    list.features.iter().map(|f| f.id).max().map_or(1, |m| m + 1)
}

/// A name for a feature the caller did not name: the kind in lower case and the first number not yet used (`rect1`, `rect2`).
fn auto_name(list: &FeatureList, kind: &str) -> String {
    let base = kind.to_lowercase();
    (1..).map(|n| format!("{base}{n}")).find(|c| list.feature_by_name(c).is_none()).expect("a free name exists")
}

/// Adds the feature `kind` built from `args` (see `geomdesc`); returns the name it was given.
pub fn add_feature(list: &mut FeatureList, kind: &str, name: Option<&str>, args: &Value) -> Result<String, Error> {
    if find(kind).is_none() {
        return Err(Error::invalid(format!("there is no operation '{kind}'")));
    }
    transact(list, |l| {
        let name = match name {
            Some(n) => n.to_string(),
            None => auto_name(l, kind),
        };
        check_name("feature", &name).map_err(bad)?;
        if l.feature_by_name(&name).is_some() {
            return Err(Error::invalid(format!("a feature named '{name}' exists; names must be unique")));
        }
        let op = build_op(l, kind, args)?;
        let id = next_id(l);
        l.features.push(Feature { id, name: name.clone(), op, suppressed: false, origin: Origin::Code, tags: vec![] });
        Ok(name)
    })
}

/// Changes some fields of a feature; `fields` holds only the fields to change, in the same form `add_feature` takes.
pub fn update_feature(list: &mut FeatureList, name: &str, fields: &Value) -> Result<(), Error> {
    let changes =
        fields.as_object().ok_or_else(|| Error::invalid("`fields` must be an object of the fields to change"))?;
    transact(list, |l| {
        let i = index_of(l, name)?;
        let kind = kind_of(&l.features[i].op);
        let mut args = args_of(l, &l.features[i].op);
        for (k, v) in changes {
            args[k] = v.clone();
        }
        // a feature may refer only to earlier ones: build against the part of the list before it
        let mut before = l.clone();
        before.features.truncate(i);
        l.features[i].op = build_op(&before, kind, &args)?;
        Ok(())
    })
}

/// Removes a feature that no later feature uses.
pub fn remove_feature(list: &mut FeatureList, name: &str) -> Result<(), Error> {
    transact(list, |l| {
        let i = index_of(l, name)?;
        let id = l.features[i].id;
        let users: Vec<&str> =
            l.features.iter().filter(|f| f.op.inputs().contains(&id)).map(|f| f.name.as_str()).collect();
        if !users.is_empty() {
            return Err(Error::invalid(format!(
                "'{name}' is used by {}; remove or change those first",
                users.join(", ")
            )));
        }
        l.features.remove(i);
        Ok(())
    })
}

pub fn rename_feature(list: &mut FeatureList, name: &str, new_name: &str) -> Result<(), Error> {
    transact(list, |l| {
        let i = index_of(l, name)?;
        check_name("feature", new_name).map_err(bad)?;
        if name != new_name && l.feature_by_name(new_name).is_some() {
            return Err(Error::invalid(format!("a feature named '{new_name}' exists")));
        }
        l.features[i].name = new_name.into();
        Ok(())
    })
}

/// Suppressed features stay in the list and are skipped when the model is evaluated.
pub fn suppress_feature(list: &mut FeatureList, name: &str, suppressed: bool) -> Result<(), Error> {
    transact(list, |l| {
        let i = index_of(l, name)?;
        l.features[i].suppressed = suppressed;
        Ok(())
    })
}

/// Adds or removes a tag of a feature (`remove` true removes it).
pub fn tag_feature(list: &mut FeatureList, name: &str, tag: &str, remove: bool) -> Result<(), Error> {
    transact(list, |l| {
        let i = index_of(l, name)?;
        let tags = &mut l.features[i].tags;
        if remove {
            tags.retain(|t| t != tag);
        } else if !tags.iter().any(|t| t == tag) {
            tags.push(tag.into());
        }
        Ok(())
    })
}

/// Sets a model parameter (a number, or an expression over the other parameters); creates it when new.
pub fn set_param(list: &mut FeatureList, name: &str, value: &Value) -> Result<(), Error> {
    let expr: Expr = match value {
        Value::Number(n) => Expr::Num(n.as_f64().ok_or_else(|| Error::invalid("the value is not a finite number"))?),
        Value::String(s) => Expr::Text(s.clone()),
        other => return Err(Error::invalid(format!("a parameter is a number or an expression text, not {other}"))),
    };
    transact(list, |l| {
        check_name("parameter", name).map_err(bad)?;
        match l.params.iter_mut().find(|p| p.name == name) {
            Some(p) => p.value = expr,
            None => l.params.push(Param { name: name.into(), value: expr }),
        }
        Ok(())
    })
}

pub fn remove_param(list: &mut FeatureList, name: &str) -> Result<(), Error> {
    transact(list, |l| {
        let n = l.params.len();
        l.params.retain(|p| p.name != name);
        if l.params.len() == n {
            return Err(Error::invalid(format!("there is no parameter '{name}'")));
        }
        Ok(())
    })
}

/// The feature list as a host reads it: each feature with its kind and its arguments in the form the add tools take.
pub fn describe(list: &FeatureList) -> Value {
    json!({
        "dim": list.dim,
        "params": list.params.iter().map(|p| json!({ "name": p.name, "value": p.value })).collect::<Vec<_>>(),
        "features": list.features.iter().map(|f| json!({
            "name": f.name,
            "kind": kind_of(&f.op),
            "suppressed": f.suppressed,
            "tags": f.tags,
            "args": args_of(list, &f.op),
        })).collect::<Vec<_>>(),
    })
}

/// Evaluates the model and returns the names a problem can refer to: regions and boundaries, before any mesh exists.
pub fn names(list: &FeatureList) -> Result<Value, Error> {
    let ev = geom::eval::Evaluator::new()
        .eval(list, geom::eval::Tier::Complex)
        .map_err(|e| Error::invalid(format!("the geometry does not evaluate: {e}")))?;
    let labels = ev.label_map().map_err(|e| Error::invalid(format!("the geometry has no consistent names: {e}")))?;
    let mut boundaries: Vec<&str> = labels.boundaries.iter().map(|b| b.name.as_str()).collect();
    boundaries.sort_unstable();
    boundaries.dedup();
    let b = ev.bounds();
    Ok(json!({
        "regions": labels.domains.iter().map(|d| json!({ "label": d.id, "name": d.name, "kind": format!("{:?}", d.kind) })).collect::<Vec<_>>(),
        "boundaries": boundaries,
        "groups": labels.groups.iter().map(|g| json!({ "name": g.name, "rule": g.rule })).collect::<Vec<_>>(),
        "bounds": { "min": b.lo, "max": b.hi },
    }))
}

/// Replaces the whole feature list (the window edits it as a whole); it must be a valid list of the same dimension class (2 or 3).
pub fn set_all(list: &mut FeatureList, new: FeatureList) -> Result<(), Error> {
    new.validate().map_err(bad)?;
    *list = new;
    Ok(())
}
