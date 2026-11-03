//! The registry: the physics a build contains, by code. The engine fills it from cargo features (physics are compiled in); a project
//! loads a few of them. Nothing outside the physics crates names a physics: everything goes through `Registry::get`.

use crate::error::Error;
use crate::physics::Physics;
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Default, Clone)]
pub struct Registry {
    by_id: BTreeMap<&'static str, Arc<dyn Physics>>,
}

/// A code is exactly three lowercase ASCII letters (ADR 12).
pub fn is_valid_code(id: &str) -> bool {
    id.len() == 3 && id.bytes().all(|b| b.is_ascii_lowercase())
}

impl Registry {
    pub fn new() -> Registry {
        Registry::default()
    }

    /// Adds a physics. Refuses a code that is not three lowercase letters and a code that is already taken.
    pub fn register(&mut self, p: Arc<dyn Physics>) -> Result<(), Error> {
        let id = p.id();
        if !is_valid_code(id) {
            return Err(Error::invalid(format!("physics code '{id}' is not three lowercase letters (ADR 12)")));
        }
        if self.by_id.contains_key(id) {
            return Err(Error::invalid(format!("physics '{id}' is registered twice")));
        }
        self.by_id.insert(id, p);
        Ok(())
    }

    pub fn get(&self, id: &str) -> Option<&dyn Physics> {
        self.by_id.get(id).map(|p| p.as_ref())
    }

    /// The codes in alphabetical order.
    pub fn ids(&self) -> Vec<&'static str> {
        self.by_id.keys().copied().collect()
    }

    pub fn len(&self) -> usize {
        self.by_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_id.is_empty()
    }

    /// The order in which `ids` must be solved so that every import finds its export first. An import with no matching export among `ids` is
    /// left to the user (the field simply is not available); a cycle is an error because version 1 supports one-way chains only.
    pub fn solve_order(&self, ids: &[&str]) -> Result<Vec<String>, Error> {
        let mut deps: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for &id in ids {
            let p = self.get(id).ok_or_else(|| Error::invalid(format!("no physics '{id}' in this build")))?;
            let mut d = vec![];
            for imp in p.imports() {
                for &other in ids {
                    if other == id {
                        continue;
                    }
                    if self.get(other).is_some_and(|q| q.exports().iter().any(|e| e.role == imp.role))
                        && !d.contains(&other.to_string())
                    {
                        d.push(other.to_string());
                    }
                }
            }
            deps.insert(id.to_string(), d);
        }
        let mut order: Vec<String> = vec![];
        while order.len() < deps.len() {
            let next = deps
                .iter()
                .find(|(k, d)| !order.contains(k) && d.iter().all(|x| order.contains(x)))
                .map(|(k, _)| k.clone());
            match next {
                Some(k) => order.push(k),
                None => return Err(Error::invalid(
                    "the physics of this project depend on each other in a cycle; only one-way chains are supported",
                )),
            }
        }
        Ok(order)
    }
}
