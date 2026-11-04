//! The versioned JSON Schema of a material file (a TOML file read as JSON data). `SCHEMA_VERSION` is the schema's own version and
//! changes with the contract exemplar `material_schema.v1.json`; the file format number `format = 1` is `spec::FORMAT`.
//! The schema is documentation and an editor aid; the loader in `spec.rs` is the authority and reports errors in plain words.

use serde_json::{json, Value};

pub const SCHEMA_VERSION: &str = "1.0";

fn quantity(desc: &str) -> Value {
    json!({
        "description": desc,
        "oneOf": [
            { "type": "string", "description": "a number with its unit, for example \"496.8 GPa\"" },
            { "type": "number", "description": "a bare number: needs the unit declared in a sibling key `<name>_unit` (or `unit` in [elastic]), unless dimensionless" },
            { "type": "object", "required": ["value"], "properties": { "value": {}, "unit": { "type": "string" } }, "additionalProperties": false }
        ]
    })
}

fn tensor(desc: &str) -> Value {
    json!({
        "description": format!("{desc}. One value (isotropic), 3 (x y z), 6 (xx yy zz yz xz xy), a symmetric 3x3 array, or {{ perp, par }} (unique axis along crystal z)"),
        "oneOf": [
            { "type": ["string", "number"] },
            { "type": "array", "items": {} },
            { "type": "object", "properties": { "perp": {}, "par": {}, "value": {}, "unit": { "type": "string" } }, "additionalProperties": false }
        ]
    })
}

fn model() -> Value {
    json!({
        "type": "object",
        "required": ["model"],
        "properties": {
            "model": { "enum": ["constant", "sellmeier", "cauchy", "tabulated"] },
            "n": { "type": "number", "description": "constant model: refractive index" },
            "eps": { "type": "number", "description": "constant model: relative permittivity (alternative to n)" },
            "k": { "type": "number", "minimum": 0, "description": "constant model: extinction coefficient; the complex index is n - j k" },
            "A": { "type": "number", "description": "sellmeier: n^2 = A + sum B x^2/(x^2 - C), default 1" },
            "B": { "type": "array", "items": { "type": "number" } },
            "C": { "type": "array", "items": { "type": "number", "exclusiveMinimum": 0 }, "description": "sellmeier: squared resonance wavelengths in wavelength_unit^2" },
            "lambda0": { "type": "array", "items": { "type": "number", "exclusiveMinimum": 0 }, "description": "sellmeier: resonance wavelengths in wavelength_unit (instead of C)" },
            "coefficients": { "type": "array", "items": { "type": "number" }, "description": "cauchy: n = c0 + c1/x^2 + c2/x^4 + ..." },
            "data": { "type": "array", "items": { "type": "array", "minItems": 2, "maxItems": 3, "items": { "type": "number" } }, "minItems": 2, "description": "tabulated: rows [wavelength, n] or [wavelength, n, k], wavelength increasing" },
            "interpolation": { "enum": ["pchip", "linear"], "default": "pchip" },
            "extrapolate": { "enum": ["error", "linear"], "default": "error" },
            "dn_dT": quantity("thermo-optic coefficient, 1/K"),
            "dn_dT_unit": { "type": "string" }
        }
    })
}

/// The JSON Schema (draft 2020-12) of format 1.
pub fn material_schema() -> Value {
    let mut optical_props = model()["properties"].as_object().cloned().unwrap();
    optical_props.insert("model".into(), json!({ "enum": ["constant", "sellmeier", "cauchy", "tabulated"] }));
    optical_props.insert("wavelength_unit".into(), json!({ "type": "string", "description": "unit of every wavelength, range and coefficient of this section, for example \"um\"; required unless all models are constant" }));
    optical_props.insert("range".into(), json!({ "type": "array", "minItems": 2, "maxItems": 2, "description": "validity range [min, max] in wavelength_unit; required for formula models; a wavelength outside it is an error" }));
    optical_props.insert("reference_temperature".into(), json!({ "type": ["string", "number"], "description": "temperature the models hold at, default 293.15 K; with dn_dT the index shifts by dn_dT (T - reference)" }));
    for axis in ["ordinary", "extraordinary", "x", "y", "z"] {
        optical_props.insert(axis.into(), json!({ "$ref": "#/$defs/model", "description": "model along this crystal axis (ordinary/extraordinary: uniaxial about z; x/y/z: biaxial)" }));
    }
    let mut elastic_props = serde_json::Map::new();
    elastic_props.insert("unit".into(), json!({ "type": "string", "description": "unit of the stiffness values written as bare numbers, for example \"GPa\"" }));
    elastic_props.insert("compliance_unit".into(), json!({ "type": "string", "description": "unit of compliance values, for example \"1/GPa\"" }));
    elastic_props.insert("class".into(), json!({ "enum": ["isotropic", "cubic", "hexagonal", "trigonal", "tetragonal", "orthotropic", "orthorhombic", "monoclinic", "triclinic"] }));
    elastic_props.insert("unique_axis".into(), json!({ "enum": ["y", "z"], "default": "y", "description": "two-fold axis of a monoclinic crystal" }));
    elastic_props.insert("matrix".into(), json!({ "type": "array", "description": "stiffness in Voigt order 1=xx 2=yy 3=zz 4=yz 5=xz 6=xy (engineering shear): 6 rows of 6, or the upper triangle as rows of 6, 5, 4, 3, 2, 1 entries" }));
    elastic_props.insert("compliance".into(), json!({ "type": "array", "description": "compliance matrix, same shapes as matrix (Voigt, engineering shear: the matrix inverse of the stiffness)" }));
    for k in ["E", "lambda", "mu", "K", "G"] {
        elastic_props.insert(k.into(), quantity("isotropic constant"));
    }
    elastic_props.insert("nu".into(), json!({ "type": "number", "exclusiveMinimum": -1, "exclusiveMaximum": 0.5, "description": "Poisson's ratio" }));
    for i in 1..=6 {
        for j in i..=6 {
            elastic_props.insert(format!("c{i}{j}"), quantity("independent stiffness constant of the class"));
            elastic_props.insert(format!("s{i}{j}"), quantity("independent compliance constant of the class"));
        }
    }
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://nyx.local/schemas/material.v1.json",
        "title": "nyx material file",
        "x-schema-version": SCHEMA_VERSION,
        "x-file-format": super::spec::FORMAT,
        "description": "A TOML material file read as data. Units are never assumed: write a quantity with its unit.",
        "type": "object",
        "required": ["format"],
        "additionalProperties": false,
        "properties": {
            "format": { "const": super::spec::FORMAT },
            "name": { "type": "string" },
            "description": { "type": "string" },
            "source": { "oneOf": [{ "type": "string" }, { "type": "array", "items": { "type": "string" } }], "description": "references the numbers come from" },
            "notes": { "type": "string" },
            "density": quantity("mass density, kg/m^3"),
            "density_unit": { "type": "string" },
            "elastic": { "type": "object", "additionalProperties": false, "properties": elastic_props, "description": "stiffness in the crystal frame; give it in exactly one way: class + constants, matrix, compliance, or isotropic constants" },
            "thermal": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "conductivity": tensor("thermal conductivity, W/(m K)"),
                    "conductivity_unit": { "type": "string" },
                    "specific_heat": quantity("specific heat, J/(kg K)"),
                    "specific_heat_unit": { "type": "string" },
                    "expansion": tensor("thermal expansion, 1/K"),
                    "expansion_unit": { "type": "string" }
                }
            },
            "electrical": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "permittivity": tensor("relative permittivity at low frequency"),
                    "permeability": { "type": "number", "exclusiveMinimum": 0 },
                    "conductivity": tensor("electrical conductivity, S/m"),
                    "conductivity_unit": { "type": "string" }
                }
            },
            "optical": {
                "type": "object",
                "additionalProperties": false,
                "description": "wavelength-dependent index; time convention e^{j omega t}, complex index n - j k. Put the model keys here (isotropic) or in per-axis tables",
                "properties": optical_props
            }
        },
        "$defs": { "model": model() }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_declares_its_versions() {
        let s = material_schema();
        assert_eq!(s["x-schema-version"], SCHEMA_VERSION);
        assert_eq!(s["properties"]["format"]["const"], 1);
    }
}
