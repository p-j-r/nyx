// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The versioned JSON Schemas of a material file (`format = 2` current schema `2.0`, plus `format = 1`
//! read-only compatibility schema `1.1`). The loader in `spec.rs` is the authority and reports errors
//! in plain words.

use super::state_t::PointGroup;
use serde_json::{json, Value};

pub const SCHEMA_VERSION: &str = "2.0";
pub const SCHEMA_VERSION_V1: &str = "1.1";

fn source_ref() -> Value {
    json!({ "$ref": "#/$defs/source" })
}

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

fn temp_scalar_or_curve(desc: &str) -> Value {
    json!({
        "description": format!("{desc}. Either a constant quantity or a temperature-curve table {{ form = \"tabulated_T\" | \"log10_poly\" | \"einstein_shift\", ... }}"),
        "oneOf": [
            quantity(desc),
            { "$ref": "#/$defs/temp_curve" }
        ]
    })
}

fn temp_tensor_or_curve(desc: &str) -> Value {
    json!({
        "description": format!("{desc}. A constant tensor, a scalar temperature-curve table, or {{ perp = <temp_curve>, par = <temp_curve> }}"),
        "oneOf": [
            tensor(desc),
            { "$ref": "#/$defs/temp_curve" },
            {
                "type": "object",
                "required": ["perp", "par"],
                "additionalProperties": false,
                "properties": {
                    "perp": temp_scalar_or_curve("in-plane component"),
                    "par": temp_scalar_or_curve("c-axis component")
                }
            }
        ]
    })
}

fn source_def() -> Value {
    json!({
        "description": "where the numbers come from: a citation text, a table { cite, url, doi, verified } (verified = YYYY-MM-DD), or a list of those; data that were not checked say so in the text",
        "oneOf": [
            { "type": "string" },
            { "type": "object", "required": ["cite"], "additionalProperties": false, "properties": { "cite": { "type": "string" }, "url": { "type": "string", "pattern": "^https?://" }, "doi": { "type": "string", "pattern": "^10[.]" }, "verified": { "type": "string", "pattern": "^[0-9]{4}-[0-9]{2}-[0-9]{2}$" } } },
            { "type": "array", "items": { "$ref": "#/$defs/source" } }
        ]
    })
}

fn model_v1() -> Value {
    json!({
        "type": "object",
        "required": ["model"],
        "properties": {
            "model": { "enum": ["constant", "sellmeier", "cauchy", "schott", "herzberger", "drude", "lorentz", "drude_lorentz", "tabulated"] },
            "source": source_ref(),
            "frequency_unit": { "type": "string", "description": "drude/lorentz: unit of the frequencies, for example \"eV\", \"THz\", \"rad/s\", \"1/cm\"" },
            "eps_inf": { "type": "number", "default": 1, "description": "drude/lorentz: permittivity at high frequency" },
            "plasma_frequency": quantity("drude/lorentz: the scale Wp of the oscillator strengths, in frequency_unit"),
            "damping": quantity("drude: damping (collision frequency) of the free carriers, in frequency_unit"),
            "drude_weight": { "type": "number", "default": 1, "description": "drude: weight f0 of the Drude term" },
            "oscillators": { "type": "array", "items": { "type": "array", "minItems": 3, "maxItems": 3, "items": { "type": "number" } }, "description": "lorentz: rows [f, resonance, damping], eps += f Wp^2/(resonance^2 - w^2 + j damping w)" },
            "n": { "type": "number", "description": "constant model: refractive index" },
            "eps": { "type": "number", "description": "constant model: relative permittivity (alternative to n)" },
            "k": { "type": "number", "minimum": 0, "description": "constant model: extinction coefficient; the complex index is n - j k" },
            "A": { "type": "number", "description": "sellmeier: n^2 = A + sum B x^2/(x^2 - C), default 1" },
            "B": { "type": "array", "items": { "type": "number" } },
            "C": { "type": "array", "items": { "type": "number", "exclusiveMinimum": 0 }, "description": "sellmeier: squared resonance wavelengths in wavelength_unit^2" },
            "lambda0": { "type": "array", "items": { "type": "number", "exclusiveMinimum": 0 }, "description": "sellmeier: resonance wavelengths in wavelength_unit (instead of C)" },
            "coefficients": { "type": "array", "items": { "type": "number" }, "description": "cauchy: n = c0 + c1/x^2 + c2/x^4 + ...; schott: a0..a5 of n^2 = a0 + a1 x^2 + a2 x^-2 + a3 x^-4 + a4 x^-6 + a5 x^-8; herzberger: A..F of n = A + B L + C L^2 + D x^2 + E x^4 + F x^6, L = 1/(x^2 - 0.028 um^2)" },
            "data": { "type": "array", "items": { "type": "array", "minItems": 2, "maxItems": 3, "items": { "type": "number" } }, "minItems": 2, "description": "tabulated: rows [wavelength, n] or [wavelength, n, k], wavelength increasing" },
            "interpolation": { "enum": ["pchip", "linear"], "default": "pchip" },
            "extrapolate": { "enum": ["error", "linear"], "default": "error" },
            "dn_dT": quantity("thermo-optic coefficient, 1/K"),
            "dn_dT_unit": { "type": "string" }
        }
    })
}

fn model_v2() -> Value {
    let mut props = model_v1()["properties"].as_object().cloned().unwrap();
    props.insert(
        "model".into(),
        json!({ "enum": ["constant", "sellmeier", "cauchy", "schott", "herzberger", "rii_formula_4", "rii_formula_5", "rii_formula_6", "drude", "lorentz", "drude_lorentz", "tabulated"] }),
    );
    props.insert(
        "k_data".into(),
        json!({ "type": "array", "items": { "type": "array", "minItems": 2, "maxItems": 2, "items": { "type": "number" } }, "minItems": 2, "description": "optional separate extinction table [wavelength, k] paired with an analytical n dispersion model" }),
    );
    props.insert("dn_dT_einstein_theta".into(), quantity("Einstein temperature Theta_E [K] saturating dn/dT -> 0 as T -> 0 K"));
    props.insert("dn_dT_einstein_theta_unit".into(), json!({ "type": "string" }));
    props.insert("dn_dT_range".into(), json!({ "type": "array", "minItems": 2, "maxItems": 2, "items": { "type": "number" }, "description": "valid temperature interval [T_min, T_max] for thermo-optic shifting" }));
    props.insert("dn_dT_range_unit".into(), json!({ "type": "string" }));
    json!({
        "type": "object",
        "required": ["model"],
        "properties": props
    })
}

fn elastic_props() -> serde_json::Map<String, Value> {
    let mut elastic_props = serde_json::Map::new();
    elastic_props.insert("source".into(), source_ref());
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
    elastic_props
}

/// The JSON Schema (draft 2020-12) of legacy format 1 (retained verbatim for v1 contract tests).
pub fn material_schema_v1() -> Value {
    let mut optical_props = model_v1()["properties"].as_object().cloned().unwrap();
    optical_props.insert("wavelength_unit".into(), json!({ "type": "string", "description": "unit of every wavelength, range and coefficient of this section, for example \"um\"; required unless all models are constant" }));
    optical_props.insert("range".into(), json!({ "type": "array", "minItems": 2, "maxItems": 2, "description": "validity range [min, max] in wavelength_unit; required for formula models; a wavelength outside it is an error" }));
    optical_props.insert("reference_temperature".into(), json!({ "type": ["string", "number"], "description": "temperature the models hold at, default 293.15 K; with dn_dT the index shifts by dn_dT (T - reference)" }));
    for axis in ["ordinary", "extraordinary", "x", "y", "z"] {
        optical_props.insert(axis.into(), json!({ "$ref": "#/$defs/model", "description": "model along this crystal axis (ordinary/extraordinary: uniaxial about z; x/y/z: biaxial)" }));
    }
    let elastic = elastic_props();
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://nyx.local/schemas/material.v1.json",
        "title": "nyx material file",
        "x-schema-version": SCHEMA_VERSION_V1,
        "x-file-format": super::spec::FORMAT_V1,
        "description": "A TOML material file read as data. Units are never assumed: write a quantity with its unit.",
        "type": "object",
        "required": ["format"],
        "additionalProperties": false,
        "properties": {
            "format": { "const": super::spec::FORMAT_V1 },
            "name": { "type": "string" },
            "description": { "type": "string" },
            "source": source_ref(),
            "extends": { "type": "string", "description": "name of a library material to start from; the tables of this file replace the base's tables as a whole" },
            "notes": { "type": "string" },
            "density": quantity("mass density, kg/m^3"),
            "density_unit": { "type": "string" },
            "elastic": { "type": "object", "additionalProperties": false, "properties": elastic, "description": "stiffness in the crystal frame; give it in exactly one way: class + constants, matrix, compliance, or isotropic constants" },
            "thermal": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "source": source_ref(),
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
                    "source": source_ref(),
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
        "$defs": {
            "model": model_v1(),
            "source": source_def()
        }
    })
}

fn thermal_v2_schema() -> Value {
    let mut p = serde_json::Map::new();
    p.insert("source".into(), source_ref());
    p.insert("conductivity".into(), temp_tensor_or_curve("thermal conductivity, W/(m K)"));
    p.insert("conductivity_unit".into(), json!({ "type": "string" }));
    p.insert("conductivity_T".into(), temp_tensor_or_curve("explicit k(T) curve, W/(m K)"));
    p.insert("conductivity_T_unit".into(), json!({ "type": "string" }));
    p.insert("specific_heat".into(), temp_scalar_or_curve("specific heat, J/(kg K)"));
    p.insert("specific_heat_unit".into(), json!({ "type": "string" }));
    p.insert("specific_heat_T".into(), temp_scalar_or_curve("explicit c_p(T) curve, J/(kg K)"));
    p.insert("specific_heat_T_unit".into(), json!({ "type": "string" }));
    p.insert("expansion".into(), temp_tensor_or_curve("thermal expansion alpha, 1/K"));
    p.insert("expansion_unit".into(), json!({ "type": "string" }));
    p.insert("expansion_T".into(), temp_tensor_or_curve("explicit alpha(T) curve, 1/K"));
    p.insert("expansion_T_unit".into(), json!({ "type": "string" }));
    p.insert("contraction".into(), temp_tensor_or_curve("integrated thermal contraction (L(T) - L_ref)/L_ref"));
    p.insert("contraction_unit".into(), json!({ "type": "string" }));
    json!({ "type": "object", "additionalProperties": false, "properties": p })
}

fn dielectric_schema() -> Value {
    let mut p = serde_json::Map::new();
    p.insert("source".into(), source_ref());
    p.insert("eps_r".into(), tensor("static/microwave relative permittivity alias"));
    p.insert("eps_clamped".into(), tensor("constant-strain microwave relative permittivity eps^S / eps_0"));
    p.insert("eps_free".into(), tensor("constant-stress static relative permittivity eps^T / eps_0"));
    p.insert("eps_inf".into(), tensor("high-frequency clamped electronic permittivity eps_inf"));
    p.insert("eps_T".into(), temp_tensor_or_curve("temperature-dependent relative permittivity"));
    p.insert("loss_tangent".into(), tensor("bulk microwave loss tangent tan(delta) >= 0"));
    p.insert("loss_ref_freq".into(), quantity("reference frequency for bulk loss tangent, Hz"));
    p.insert("loss_ref_freq_unit".into(), json!({ "type": "string" }));
    p.insert("tls_tan_delta_0".into(), json!({ "type": "number", "minimum": 0 }));
    p.insert("tls_e_crit".into(), quantity("critical saturation field E_c for TLS dielectric loss, V/m"));
    p.insert("tls_e_crit_unit".into(), json!({ "type": "string" }));
    json!({ "type": "object", "additionalProperties": false, "properties": p })
}

fn superconducting_schema() -> Value {
    let mut p = serde_json::Map::new();
    p.insert("source".into(), source_ref());
    p.insert("tc".into(), quantity("superconducting critical temperature T_c, K"));
    p.insert("tc_unit".into(), json!({ "type": "string" }));
    p.insert("gap_0".into(), quantity("zero-temperature BCS energy gap Delta(0), J or ueV"));
    p.insert("gap_0_unit".into(), json!({ "type": "string" }));
    p.insert("lambda_L".into(), quantity("clean-limit London penetration depth, m or nm"));
    p.insert("lambda_L_unit".into(), json!({ "type": "string" }));
    p.insert("xi_0".into(), quantity("Pippard/BCS coherence length, m or nm"));
    p.insert("xi_0_unit".into(), json!({ "type": "string" }));
    p.insert("mean_free_path".into(), quantity("electron mean free path ell, m or nm"));
    p.insert("mean_free_path_unit".into(), json!({ "type": "string" }));
    p.insert("sigma_n".into(), quantity("normal-state conductivity sigma_n, S/m"));
    p.insert("sigma_n_unit".into(), json!({ "type": "string" }));
    p.insert("dynes_gamma".into(), quantity("Dynes sub-gap pair-breaking broadening Gamma, J or ueV"));
    p.insert("dynes_gamma_unit".into(), json!({ "type": "string" }));
    p.insert("r_res".into(), quantity("non-BCS residual sheet resistance R_res, Ohm/sq"));
    p.insert("r_res_unit".into(), json!({ "type": "string" }));
    p.insert("x_qp_floor".into(), json!({ "type": "number", "minimum": 0 }));
    p.insert("nominal_thickness".into(), quantity("default film thickness, m or nm"));
    p.insert("nominal_thickness_unit".into(), json!({ "type": "string" }));
    p.insert("film_model".into(), json!({ "enum": ["kerr", "two_sided", "thin_limit"] }));
    json!({ "type": "object", "additionalProperties": false, "required": ["tc"], "properties": p })
}

fn defs_v2_schema() -> Value {
    let mut iface_p = serde_json::Map::new();
    iface_p.insert("source".into(), source_ref());
    iface_p.insert("oxide".into(), json!({ "type": "string" }));
    iface_p.insert("thickness".into(), quantity("native oxide layer thickness, nm"));
    iface_p.insert("thickness_unit".into(), json!({ "type": "string" }));
    iface_p.insert("eps_r".into(), json!({ "type": "number", "exclusiveMinimum": 0 }));
    iface_p.insert("tan_delta_tls_0".into(), json!({ "type": "number", "minimum": 0 }));
    iface_p.insert("e_crit".into(), quantity("TLS critical saturation electric field, V/m"));
    iface_p.insert("e_crit_unit".into(), json!({ "type": "string" }));

    let mut tc_p = serde_json::Map::new();
    tc_p.insert("form".into(), json!({ "enum": ["tabulated_T", "log10_poly", "einstein_shift"] }));
    tc_p.insert("source".into(), source_ref());
    tc_p.insert("T".into(), json!({ "type": "array", "items": { "type": "number" } }));
    tc_p.insert("T_unit".into(), json!({ "type": "string" }));
    tc_p.insert("T_range".into(), json!({ "type": "array", "minItems": 2, "maxItems": 2, "items": { "type": "number" } }));
    tc_p.insert("values".into(), json!({ "type": "array" }));
    tc_p.insert("values_unit".into(), json!({ "type": "string" }));
    tc_p.insert("coeffs".into(), json!({ "type": "array", "items": { "type": "number" } }));
    tc_p.insert("ref_value".into(), quantity("value at t_ref"));
    tc_p.insert("t_ref".into(), quantity("reference temperature"));
    tc_p.insert("theta_e".into(), quantity("Einstein characteristic temperature [K]"));
    tc_p.insert("alpha".into(), json!({ "type": "number" }));
    tc_p.insert("unit".into(), json!({ "type": "string" }));

    let mut defs = serde_json::Map::new();
    defs.insert("model".into(), model_v2());
    defs.insert("source".into(), source_def());
    defs.insert(
        "interface_layer".into(),
        json!({
            "type": "object",
            "additionalProperties": false,
            "required": ["thickness", "eps_r", "tan_delta_tls_0"],
            "properties": iface_p
        }),
    );
    defs.insert(
        "temp_curve".into(),
        json!({ "type": "object", "required": ["form"], "properties": tc_p }),
    );
    Value::Object(defs)
}

/// The JSON Schema (draft 2020-12) of format 2 (`SCHEMA_VERSION = "2.0"`).
pub fn material_schema() -> Value {
    let mut optical_props = model_v2()["properties"].as_object().cloned().unwrap();
    optical_props.insert("wavelength_unit".into(), json!({ "type": "string", "description": "unit of every wavelength, range and coefficient of this section, for example \"um\"; required unless all models are constant" }));
    optical_props.insert("range".into(), json!({ "type": "array", "minItems": 2, "maxItems": 2, "description": "validity range [min, max] in wavelength_unit; required for formula models" }));
    optical_props.insert("reference_temperature".into(), json!({ "type": ["string", "number"], "description": "temperature the models hold at, default 293.15 K" }));
    for axis in ["ordinary", "extraordinary", "x", "y", "z"] {
        optical_props.insert(axis.into(), json!({ "$ref": "#/$defs/model", "description": "model along this crystal axis" }));
    }

    let mut piezo_props = serde_json::Map::new();
    piezo_props.insert("source".into(), source_ref());
    piezo_props.insert("form".into(), json!({ "enum": ["stress_charge", "strain_charge"] }));
    piezo_props.insert("unit".into(), json!({ "type": "string", "description": "C/m^2 for stress_charge (e), pC/N or pm/V for strain_charge (d)" }));
    piezo_props.insert("matrix".into(), json!({ "type": "array", "description": "3x6 coupling matrix in textbook Voigt (1=xx..6=xy)" }));
    for i in 1..=3 {
        for j in 1..=6 {
            piezo_props.insert(format!("e{i}{j}"), quantity("stress-charge piezoelectric entry, C/m^2"));
            piezo_props.insert(format!("d{i}{j}"), quantity("strain-charge piezoelectric entry, pC/N or m/V"));
        }
    }

    let mut pe_props = serde_json::Map::new();
    pe_props.insert("source".into(), source_ref());
    pe_props.insert("wavelength".into(), quantity("reference optical wavelength"));
    pe_props.insert("wavelength_unit".into(), json!({ "type": "string" }));
    pe_props.insert("matrix".into(), json!({ "type": "array", "description": "6x6 dimensionless Pockels photoelastic matrix (may be asymmetric and indefinite)" }));
    for i in 1..=6 {
        for j in 1..=6 {
            pe_props.insert(format!("p{i}{j}"), json!({ "type": "number" }));
        }
    }

    let mut eo_sub_props = serde_json::Map::new();
    eo_sub_props.insert("unit".into(), json!({ "type": "string" }));
    eo_sub_props.insert("matrix".into(), json!({ "type": "array", "description": "6x3 Pockels electro-optic matrix in m/V (default pm/V)" }));
    for i in 1..=6 {
        for j in 1..=3 {
            eo_sub_props.insert(format!("r{i}{j}"), quantity("linear electro-optic entry, pm/V"));
        }
    }
    let eo_sub_schema = json!({ "type": "object", "additionalProperties": false, "properties": eo_sub_props });

    let mut eo_props = serde_json::Map::new();
    eo_props.insert("source".into(), source_ref());
    eo_props.insert("wavelength".into(), quantity("reference optical wavelength"));
    eo_props.insert("wavelength_unit".into(), json!({ "type": "string" }));
    eo_props.insert("unit".into(), json!({ "type": "string" }));
    eo_props.insert("clamped".into(), eo_sub_schema.clone());
    eo_props.insert("free".into(), eo_sub_schema);
    eo_props.insert("matrix".into(), json!({ "type": "array" }));
    eo_props.insert("state".into(), json!({ "enum": ["clamped", "free"] }));
    for i in 1..=6 {
        for j in 1..=3 {
            eo_props.insert(format!("r{i}{j}"), quantity("linear electro-optic entry, pm/V"));
        }
    }

    let mut ve_props = serde_json::Map::new();
    ve_props.insert("source".into(), source_ref());
    ve_props.insert("unit".into(), json!({ "type": "string", "description": "dynamic viscosity unit, e.g. mPa*s or Pa*s" }));
    ve_props.insert("matrix".into(), json!({ "type": "array", "description": "symmetric 6x6 phonon viscosity matrix eta_IJ" }));
    ve_props.insert("cryo_qm_floor".into(), json!({ "type": "number", "exclusiveMinimum": 0 }));
    for i in 1..=6 {
        for j in i..=6 {
            ve_props.insert(format!("eta{i}{j}"), quantity("viscosity entry, Pa*s"));
        }
    }

    let pg_enum: Vec<&str> = PointGroup::ALL.iter().map(|p| p.as_str()).collect();
    let elastic = elastic_props();

    let mut top = serde_json::Map::new();
    top.insert("format".into(), json!({ "const": super::spec::FORMAT }));
    top.insert("name".into(), json!({ "type": "string" }));
    top.insert("description".into(), json!({ "type": "string" }));
    top.insert("source".into(), source_ref());
    top.insert("extends".into(), json!({ "type": "string", "description": "base material name; format 2 merges nested tables recursively with representation reset guards" }));
    top.insert("notes".into(), json!({ "type": "string" }));
    top.insert("kind".into(), json!({ "enum": ["crystalline", "polycrystalline", "amorphous", "superconductor", "semiconductor", "dielectric", "metal", "polymer", "fluid"] }));
    top.insert("point_group".into(), json!({ "enum": pg_enum, "description": "Hermann-Mauguin point group; centrosymmetric groups forbid piezoelectric and electro_optic" }));
    top.insert("crystal_setting".into(), json!({ "enum": ["IEEE-1987", "IRE-1949"] }));
    top.insert("reference_temperature".into(), quantity("reference temperature T_ref, default 293.15 K"));
    top.insert("reference_temperature_unit".into(), json!({ "type": "string" }));
    top.insert("density".into(), quantity("mass density, kg/m^3"));
    top.insert("density_unit".into(), json!({ "type": "string" }));
    top.insert("residual_stress".into(), tensor("built-in thin-film residual stress tensor, Pa"));
    top.insert("residual_stress_unit".into(), json!({ "type": "string" }));
    top.insert("elastic".into(), json!({ "type": "object", "additionalProperties": false, "properties": elastic }));
    top.insert("thermal".into(), thermal_v2_schema());
    top.insert(
        "electrical".into(),
        json!({
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "source": source_ref(),
                "permittivity": tensor("low-frequency relative permittivity (fallback when [dielectric] is absent)"),
                "permeability": { "type": "number", "exclusiveMinimum": 0 },
                "conductivity": tensor("electrical conductivity, S/m"),
                "conductivity_unit": { "type": "string" }
            }
        }),
    );
    top.insert("dielectric".into(), dielectric_schema());
    top.insert("superconducting".into(), superconducting_schema());
    top.insert("piezoelectric".into(), json!({ "type": "object", "additionalProperties": false, "properties": piezo_props }));
    top.insert("photoelastic".into(), json!({ "type": "object", "additionalProperties": false, "properties": pe_props }));
    top.insert("electro_optic".into(), json!({ "type": "object", "additionalProperties": false, "properties": eo_props }));
    top.insert("viscoelastic".into(), json!({ "type": "object", "additionalProperties": false, "properties": ve_props }));
    top.insert(
        "interfaces".into(),
        json!({
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "ma": { "$ref": "#/$defs/interface_layer" },
                "ms": { "$ref": "#/$defs/interface_layer" },
                "sa": { "$ref": "#/$defs/interface_layer" }
            }
        }),
    );
    top.insert(
        "optical".into(),
        json!({
            "type": "object",
            "additionalProperties": false,
            "description": "wavelength-dependent complex refractive index n - j k (e^{j omega t} convention)",
            "properties": optical_props
        }),
    );

    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://nyx.local/schemas/material.v2.json",
        "title": "nyx material file",
        "x-schema-version": SCHEMA_VERSION,
        "x-file-format": super::spec::FORMAT,
        "description": "Multi-physics material definition (format 2). Supports cryogenic T-curves, superconductors + surface oxide TLS interfaces, microwave dielectric loss, piezoelectric Legendre duals, photoelastic, electro-optic, and viscoelastic tensors.",
        "type": "object",
        "required": ["format"],
        "additionalProperties": false,
        "properties": Value::Object(top),
        "$defs": defs_v2_schema()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_declares_its_versions() {
        let v1 = material_schema_v1();
        assert_eq!(v1["x-schema-version"], SCHEMA_VERSION_V1);
        assert_eq!(v1["properties"]["format"]["const"], 1);
        let v2 = material_schema();
        assert_eq!(v2["x-schema-version"], SCHEMA_VERSION);
        assert_eq!(v2["properties"]["format"]["const"], 2);
    }
}
