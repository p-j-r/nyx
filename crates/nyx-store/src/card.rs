//! The concept card: a summary of an experiment generated from its spec, so nothing has to be typed to make it findable.
//!
//! The spec is read as a `serde_json::Value` of the saved project (`.nyx.xml` through `femcore::xmlio`), tolerantly: a field that is not
//! there is skipped, so a spec written by a newer nyx still gets a card. Materials carry numbers only (no names), so a name such as
//! "sapphire" is added only when the refractive index of an optical domain matches a known material within 1 %, and the card says
//! it is inferred ("Al2O3 (inferred from n 1.746)").

use crate::model::{Card, Dim, GeoCard};
use serde_json::Value;

const C0: f64 = 299_792_458.0;

/// (name, synonyms for search, refractive index near 1550 nm)
const KNOWN: &[(&str, &str, f64)] = &[
    ("Si", "silicon", 3.4757),
    ("SiO2", "silica oxide glass", 1.444),
    ("Al2O3", "sapphire alumina", 1.746),
    ("Si3N4", "silicon nitride nitride", 1.996),
    ("air", "vacuum", 1.0),
];

fn num(v: &Value) -> Option<f64> {
    v.as_f64()
}

fn get<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    path.iter().try_fold(v, |v, k| v.get(k))
}

/// `{ "Variant": {...} }` as (variant, payload)
fn variant(v: &Value) -> Option<(&str, &Value)> {
    let o = v.as_object()?;
    if o.len() == 1 {
        o.iter().next().map(|(k, v)| (k.as_str(), v))
    } else {
        None
    }
}

fn arr3(v: &Value) -> Option<[f64; 3]> {
    let a = v.as_array()?;
    Some([num(a.first()?)?, num(a.get(1)?)?, num(a.get(2)?)?])
}

fn dim(out: &mut Vec<Dim>, name: &str, of: &str, m: f64) {
    if m.is_finite() && m > 0.0 {
        out.push(Dim { name: name.into(), of: of.into(), m });
    }
}

fn shape_dims(out: &mut Vec<Dim>, shape: &Value, of: &str) {
    let Some((kind, p)) = variant(shape) else { return };
    match kind {
        "Rect" => {
            let g = |k: &str| p.get(k).and_then(num);
            if let (Some(x0), Some(x1), Some(y0), Some(y1)) = (g("x0"), g("x1"), g("y0"), g("y1")) {
                dim(out, "width", of, (x1 - x0).abs());
                dim(out, "height", of, (y1 - y0).abs());
            }
        }
        "Circle" => {
            if let Some(r) = p.get("r").and_then(num) {
                dim(out, "diameter", of, 2.0 * r);
            }
        }
        _ => {}
    }
}

fn geometry(proj: &Value) -> (GeoCard, u8) {
    let Some((family, p)) = get(proj, &["geometry", "geometry"]).and_then(variant) else {
        return (GeoCard { family: "unknown".into(), summary: "geometry not recognised".into(), dims: vec![] }, 2);
    };
    let mut dims = vec![];
    let mut dim_n = 2u8;
    let summary;
    match family {
        "Rectangle" => {
            let g = |k: &str| p.get(k).and_then(num).unwrap_or(0.0);
            dim(&mut dims, "lx", "frame", g("lx"));
            dim(&mut dims, "ly", "frame", g("ly"));
            summary = "rectangle".to_string();
        }
        "Block" => {
            dim_n = 3;
            let g = |k: &str| p.get(k).and_then(num).unwrap_or(0.0);
            for k in ["lx", "ly", "lz"] {
                dim(&mut dims, k, "frame", g(k));
            }
            summary = "block".into();
        }
        "Shapes" => {
            let parts = get(p, &["spec", "parts"]).and_then(|v| v.as_array()).cloned().unwrap_or_default();
            for part in &parts {
                let of = format!("domain {}", part.get("domain").and_then(|d| d.as_i64()).unwrap_or(0));
                if let Some(s) = part.get("shape") {
                    shape_dims(&mut dims, s, &of);
                }
            }
            let layers = p.get("extrude").and_then(|e| e.as_array()).cloned().unwrap_or_default();
            if !layers.is_empty() {
                dim_n = 3;
                for (i, l) in layers.iter().enumerate() {
                    if let Some(t) = l.get("thickness").and_then(num) {
                        dim(&mut dims, &format!("layer {} thickness", i + 1), "extrusion", t);
                    }
                }
            }
            summary = format!("{} shape{}{}", parts.len(), if parts.len() == 1 { "" } else { "s" }, if dim_n == 3 { ", extruded" } else { "" });
        }
        "Solids" => {
            dim_n = 3;
            let bodies = get(p, &["spec", "bodies"]).and_then(|v| v.as_array()).cloned().unwrap_or_default();
            for b in &bodies {
                let of = format!("domain {}", b.get("domain").and_then(|d| d.as_i64()).unwrap_or(0));
                let Some((kind, s)) = b.get("solid").and_then(variant) else { continue };
                match kind {
                    "Box" => {
                        if let (Some(a), Some(z)) = (s.get("min").and_then(arr3), s.get("max").and_then(arr3)) {
                            for (k, n) in ["x", "y", "z"].iter().enumerate() {
                                dim(&mut dims, &format!("box {n}"), &of, (z[k] - a[k]).abs());
                            }
                        }
                    }
                    "Sphere" => {
                        if let Some(r) = s.get("r").and_then(num) {
                            dim(&mut dims, "diameter", &of, 2.0 * r);
                        }
                    }
                    "Cylinder" => {
                        if let Some(r) = s.get("r").and_then(num) {
                            dim(&mut dims, "diameter", &of, 2.0 * r);
                        }
                        if let Some(l) = s.get("len").and_then(num) {
                            dim(&mut dims, "length", &of, l);
                        }
                    }
                    "Extrusion" => {
                        if let Some(sh) = s.get("profile") {
                            shape_dims(&mut dims, sh, &of);
                        }
                        if let (Some(a), Some(b)) = (s.get("z0").and_then(num), s.get("z1").and_then(num)) {
                            dim(&mut dims, "height", &of, (b - a).abs());
                        }
                    }
                    _ => {}
                }
            }
            summary = format!("{} solid{}", bodies.len(), if bodies.len() == 1 { "" } else { "s" });
        }
        other => summary = other.to_lowercase(),
    }
    let family = match family {
        "Rectangle" => "rectangle with a structured mesh",
        "Block" => "block with a structured mesh",
        "Shapes" => "2D shapes",
        "Solids" => "3D solids",
        o => o,
    };
    (GeoCard { family: family.into(), summary, dims }, dim_n)
}

fn physics_words(name: &str) -> &'static str {
    match name {
        "Electrostatics" => "electrostatics potential field capacitance electric",
        "Heat" => "heat thermal temperature conduction",
        "Solid" => "solid mechanics structure stress vibration eigenfrequency elastic",
        "Emw" => "em waves electromagnetic optical waveguide modes neff photonics",
        _ => "",
    }
}

fn physics_label(name: &str) -> &'static str {
    match name {
        "Electrostatics" => "electrostatics",
        "Heat" => "heat transfer",
        "Solid" => "solid mechanics",
        "Emw" => "EM waves",
        _ => "unknown physics",
    }
}

fn study_label(study: &Value) -> (String, Option<f64>) {
    let Some((k, p)) = variant(study).map(|(k, p)| (k, p.clone())).or_else(|| study.as_str().map(|s| (s, Value::Null))) else {
        return ("unknown".into(), None);
    };
    let wavelength = if k == "Modes" { p.get("freq").and_then(num).filter(|f| *f > 0.0).map(|f| C0 / f) } else { None };
    let label = match k {
        "Modes" => "modes",
        "Eigenfrequency" => "eigenfrequencies",
        "Stationary" => "stationary",
        "Transient" => "transient",
        "Static" => "static",
        o => o,
    };
    (label.into(), wavelength)
}

/// The card of an experiment from its project value (`Project` serialised to JSON, or read from `.nyx.xml`).
pub fn from_value(proj: &Value) -> Card {
    let active = proj.get("active").and_then(|a| a.as_str()).unwrap_or("Emw").to_string();
    let spec = proj.get("specs").and_then(|s| s.as_array()).and_then(|a| a.iter().find(|s| s.get("physics").and_then(|p| p.as_str()) == Some(&active)));
    let (geo, dimn) = geometry(proj);
    let (study, wavelength) = spec.and_then(|s| s.get("study")).map(study_label).unwrap_or(("unknown".into(), None));
    // materials: refractive index per domain for the optical physics, name when it matches a known material
    let mut materials = vec![];
    let mut keywords: Vec<String> = physics_words(&active).split_whitespace().map(String::from).collect();
    if active == "Emw" {
        if let Some(domains) = spec.and_then(|s| get(s, &["materials", "domains"])).and_then(|d| d.as_object()) {
            for (d, props) in domains {
                let Some(eps) = props.get("eps_r").and_then(num).filter(|e| *e > 0.0) else { continue };
                let n = eps.sqrt();
                match KNOWN.iter().find(|(_, _, kn)| ((n - kn) / kn).abs() < 0.01) {
                    Some((name, syn, _)) => {
                        materials.push(format!("domain {d}: {name} (inferred from n {n:.3})"));
                        keywords.push((*name).to_string());
                        keywords.extend(syn.split_whitespace().map(String::from));
                    }
                    None => materials.push(format!("domain {d}: n {n:.3}")),
                }
            }
        }
    }
    keywords.extend(geo.family.split_whitespace().map(|s| s.to_lowercase()));
    keywords.sort();
    keywords.dedup();
    let wl = wavelength.map(|w| format!(", {:.0} nm", w * 1e9)).unwrap_or_default();
    let summary = format!("{} {}D, {} {}{}", geo.summary, dimn, physics_label(&active), study, wl);
    Card { physics: active, study, dim: dimn, wavelength_m: wavelength, geometry: geo, materials, keywords, summary }
}

/// The card of a `.nyx.xml` text.
pub fn from_xml(xml: &str) -> Result<Card, String> {
    let (root, v) = femcore::xmlio::xml_to_value(xml)?;
    if root != "nyx_project" {
        return Err(format!("not a nyx project: root element is <{root}>"));
    }
    Ok(from_value(&v))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn project() -> Value {
        json!({
            "geometry": {"geometry": {"Shapes": {"spec": {"parts": [
                {"shape": {"Rect": {"x0": 0.0, "y0": 0.0, "x1": 5e-7, "y1": 3e-7}}, "domain": 1},
                {"shape": {"Circle": {"cx": 0.0, "cy": 0.0, "r": 1.5e-7}}, "domain": 2}]}, "extrude": null}}},
            "active": "Emw",
            "specs": [{"physics": "Emw", "study": {"Modes": {"freq": 1.9341e14, "n_target": 2.4, "evanescent": false, "n_modes": 4}},
                       "materials": {"domains": {"1": {"eps_r": 12.081}, "2": {"eps_r": 3.05}, "3": {"eps_r": 7.7}}}}]
        })
    }

    #[test]
    fn the_card_reads_dims_study_wavelength_and_infers_materials() {
        let c = from_value(&project());
        assert_eq!(c.physics, "Emw");
        assert_eq!(c.study, "modes");
        assert_eq!(c.dim, 2);
        assert!((c.wavelength_m.unwrap() - 1.55e-6).abs() < 1e-9);
        assert!(c.geometry.dims.iter().any(|d| d.name == "width" && (d.m - 5e-7).abs() < 1e-15));
        assert!(c.geometry.dims.iter().any(|d| d.name == "diameter" && (d.m - 3e-7).abs() < 1e-15));
        assert!(c.materials.iter().any(|m| m.contains("Si") && m.contains("inferred")));
        assert!(c.materials.iter().any(|m| m.contains("Al2O3")));
        assert!(c.materials.iter().any(|m| m.starts_with("domain 3: n 2.775")));
        assert!(c.keywords.contains(&"sapphire".to_string()));
        assert!(c.summary.contains("1550 nm"), "{}", c.summary);
    }

    #[test]
    fn a_strange_spec_still_gets_a_card() {
        let c = from_value(&json!({"geometry": 3}));
        assert_eq!(c.geometry.family, "unknown");
        assert!(!c.summary.is_empty());
    }
}
