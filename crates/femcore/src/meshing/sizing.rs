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
