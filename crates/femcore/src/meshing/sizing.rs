//! The neutral request a physics makes of the mesher: what to resolve, with no knowledge of any physics or of `ProblemSpec`.
//!
//! A physics (`nyx_physics::Physics::sizing`) states the element order it will use, the shortest wavelength it wants resolved in each
//! region, and a `headroom` factor (an eigenfrequency study finds modes above the search frequency, so the wavelengths are divided by it).
//! The mesher turns that into element sizes (`physics_size`), with its own points-per-wavelength calibration. Regions are named; the
//! caller resolves a name to a domain label (`Sizing::resolve`) with the label map of the meshed model, or, for a mesh that only has
//! numbers, with [`label_names`] (the key is the decimal label). The serde shape is the one of the `physics_descriptors` contract.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// What the mesher needs from a physics to size the elements before any solution exists.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Sizing {
    /// element order the problem asks for
    pub order: u8,
    /// factor on the search frequency of an eigenfrequency study (the answer may lie above the guess)
    pub headroom: f64,
    /// the shortest wavelength (m) the physics wants resolved, by region name; regions without a wave are absent
    pub wavelengths: BTreeMap<String, f64>,
}

impl Sizing {
    /// A request for a mesh that only has domain numbers: the keys are the decimal labels.
    pub fn by_label(order: u8, headroom: f64, wavelengths: &BTreeMap<i32, f64>) -> Sizing {
        Sizing { order, headroom, wavelengths: wavelengths.iter().map(|(d, l)| (d.to_string(), *l)).collect() }
    }

    /// The same request with a name that ends in `*` replaced by the regions it stands for: `Air*` is `Air` and its absorbing pieces
    /// `Air.pml_xmin`, `Air.pml_xmax` (the rule `set_material` and the solvers use). `names` are the region names of the model. A name that is
    /// given in full keeps its own wavelength; a `*` that matches nothing adds nothing (the solver refuses it when it reads the materials).
    pub fn expand_wildcards(&self, names: &[&str]) -> Sizing {
        let mut out: BTreeMap<String, f64> = self.wavelengths.iter().filter(|(k, _)| !k.ends_with('*')).map(|(k, l)| (k.clone(), *l)).collect();
        for (k, lam) in self.wavelengths.iter().filter(|(k, _)| k.ends_with('*')) {
            let base = &k[..k.len() - 1];
            for n in names.iter().filter(|n| **n == base || n.strip_prefix(base).is_some_and(|rest| rest.starts_with('.'))) {
                out.entry((*n).to_string()).or_insert(*lam);
            }
        }
        Sizing { order: self.order, headroom: self.headroom, wavelengths: out }
    }

    /// The wavelengths by domain label, the names resolved by `domain`. A region the model does not have is an error.
    pub fn resolve(&self, domain: &dyn Fn(&str) -> Option<i32>) -> Result<BTreeMap<i32, f64>, String> {
        let mut out = BTreeMap::new();
        for (name, &lam) in &self.wavelengths {
            let d = domain(name).ok_or_else(|| format!("the sizing names region '{name}', which the model does not have"))?;
            out.insert(d, lam);
        }
        Ok(out)
    }
}

/// Resolver for meshes without names: the key of a region is its decimal domain label.
pub fn label_names(name: &str) -> Option<i32> {
    name.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_star_stands_for_the_region_and_its_absorbing_pieces_and_a_full_name_keeps_its_own_wavelength() {
        let s = Sizing {
            order: 2,
            headroom: 1.0,
            wavelengths: [("Air*", 1.5e-6), ("Si", 4e-7), ("Sapphire*", 8e-7), ("Sapphire.pml_xmin", 7e-7), ("Gone*", 1e-6)].iter().map(|(k, l)| (k.to_string(), *l)).collect(),
        };
        let names = ["Air", "Air.pml_xmin", "Air.pml_xmax", "Airy", "Si", "Sapphire", "Sapphire.pml_xmin", "Sapphire.pml_xmax"];
        let e = s.expand_wildcards(&names);
        let keys: Vec<&str> = e.wavelengths.keys().map(String::as_str).collect();
        assert_eq!(keys, vec!["Air", "Air.pml_xmax", "Air.pml_xmin", "Sapphire", "Sapphire.pml_xmax", "Sapphire.pml_xmin", "Si"], "`Airy` is not `Air*`, `Gone*` adds nothing");
        assert_eq!(e.wavelengths["Sapphire.pml_xmin"], 7e-7, "a name given in full keeps its own value");
        assert_eq!(e.wavelengths["Sapphire.pml_xmax"], 8e-7);
        // every name resolves now
        let id = |n: &str| names.iter().position(|m| *m == n).map(|i| i as i32 + 1);
        assert!(e.resolve(&id).is_ok());
        assert!(s.resolve(&id).is_err(), "before, `Air*` named a region the model does not have");
    }

    #[test]
    fn labels_round_trip_and_unknown_regions_are_refused() {
        let s = Sizing::by_label(2, 1.5, &[(1, 5e-7), (3, 2e-7)].into_iter().collect());
        assert_eq!(s.resolve(&label_names).unwrap(), [(1, 5e-7), (3, 2e-7)].into_iter().collect());
        let named = Sizing { order: 1, headroom: 1.0, wavelengths: [("core".to_string(), 1e-6)].into_iter().collect() };
        assert!(named.resolve(&label_names).is_err());
        assert_eq!(named.resolve(&|n| (n == "core").then_some(7)).unwrap()[&7], 1e-6);
        let json = serde_json::to_value(&s).unwrap();
        assert_eq!(serde_json::from_value::<Sizing>(json).unwrap(), s);
    }
}
