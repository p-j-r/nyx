// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-Rust importers from open material databases (`refractiveindex.info`, `AFLOW`, `Materials Project`, `NumBAT`) into validated
//! `format = 1` material TOML text (`MaterialSpec`), plus a curated catalog of remote presets for on-demand download.
//!
//! This module contains no networking code: callers (such as `nyx-engine` or CLI front-ends) fetch HTTPS payloads or read local files and
//! hand the raw text to [`import_to_toml`]. Every converted file marks automatic imports with `"UNVERIFIED"` in its source citation so the
//! user knows to check bulk 0 K DFT or unverified handbook values against their process before quantitative use.

use super::elastic::Elastic;
use super::optics::{Axes, AxisModel, Dispersion, Extrapolate, Interp, Optical};
use super::reader::MaterialError;
use super::source::Source;
use super::spec::{Library, MaterialSpec, FORMAT};
use super::tensor::TensorProp;
use super::M6;
use serde_json::Value;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImportFormat {
    RefractiveIndexInfo,
    Aflow,
    MaterialsProject,
    Numbat,
}

impl ImportFormat {
    pub fn parse(s: &str) -> Result<Self, MaterialError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "rii" | "refractiveindex" | "refractiveindex.info" => Ok(Self::RefractiveIndexInfo),
            "aflow" => Ok(Self::Aflow),
            "mp" | "materials_project" | "materials-project" | "materialsproject" => Ok(Self::MaterialsProject),
            "numbat" => Ok(Self::Numbat),
            other => Err(err("format", format!("unknown import format '{other}'; use rii, aflow, mp or numbat"))),
        }
    }

    pub fn id(self) -> &'static str {
        match self {
            Self::RefractiveIndexInfo => "rii",
            Self::Aflow => "aflow",
            Self::MaterialsProject => "mp",
            Self::Numbat => "numbat",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::RefractiveIndexInfo => "refractiveindex.info (YAML)",
            Self::Aflow => "AFLOW AEL/AGL (JSON)",
            Self::MaterialsProject => "Materials Project (JSON)",
            Self::Numbat => "NumBAT material_data (JSONC)",
        }
    }
}

pub struct ImportOptions<'a> {
    pub name: &'a str,
    pub description: Option<&'a str>,
    pub extends: Option<&'a str>,
    pub source_url: Option<&'a str>,
}

#[derive(Clone, Copy, Debug)]
pub struct CatalogEntry {
    pub id: &'static str,
    pub domain: &'static str,
    pub format: ImportFormat,
    pub target: &'static str,
    pub extends: Option<&'static str>,
    pub summary: &'static str,
}

/// Curated remote presets across photonics, phononics and superconducting/cryogenic metals that are not shipped directly in `materials/*.toml`.
pub const REMOTE_CATALOG: &[CatalogEntry] = &[
    CatalogEntry {
        id: "si_green_2008",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/Si/nk/Green-2008",
        extends: Some("silicon"),
        summary: "Intrinsic c-Si 300 K KK-consistent n,k table (Green 2008, 0.25-1.45 um)",
    },
    CatalogEntry {
        id: "si_franta_10k",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/Si/nk/Franta-10K",
        extends: Some("silicon"),
        summary: "Cryogenic FZ c-Si at 10 K incorporating Frey 2006 (Franta 2018, 0.03-20 um)",
    },
    CatalogEntry {
        id: "ge_burnett",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/Ge/nk/Burnett",
        extends: None,
        summary: "Single-crystal germanium NIR/LWIR Sellmeier (Burnett 2016, 1.8-14 um, n~4.0)",
    },
    CatalogEntry {
        id: "inp_pettit",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/InP/nk/Pettit",
        extends: None,
        summary: "Indium phosphide InP NIR/MIR Sellmeier (Pettit & Turner 1965, 0.95-10 um)",
    },
    CatalogEntry {
        id: "sic_4h_o",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/SiC/nk/Wang-4H-o",
        extends: None,
        summary: "4H-SiC ordinary Sellmeier for integrated SiC-on-insulator (Wang 2013, 0.40-2.32 um)",
    },
    CatalogEntry {
        id: "tio2_rutile_o",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/TiO2/nk/Devore-o",
        extends: None,
        summary: "Rutile TiO2 ordinary Sellmeier (Devore 1951, 0.43-1.53 um)",
    },
    CatalogEntry {
        id: "bbo_o",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/BaB2O4/nk/Eimerl-o",
        extends: None,
        summary: "Beta-barium borate BBO ordinary Sellmeier (Eimerl 1987, 0.22-1.06 um)",
    },
    CatalogEntry {
        id: "schott_nbk7",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "specs/schott/optical/N-BK7",
        extends: None,
        summary: "Schott N-BK7 borosilicate crown optical glass (Schott formula)",
    },
    CatalogEntry {
        id: "pmma_sultanova",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "organic/(C5H8O2)n - poly(methyl methacrylate)/nk/Sultanova",
        extends: None,
        summary: "PMMA polymer Sellmeier for e-beam resist / cladding (Sultanova 2009, 0.44-1.05 um)",
    },
    CatalogEntry {
        id: "hfo2_alkuhaili",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/HfO2/nk/Al-Kuhaili",
        extends: None,
        summary: "Hafnium dioxide HfO2 thin-film optical constants (Al-Kuhaili 2004)",
    },
    CatalogEntry {
        id: "sio2_franta",
        domain: "photonics",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/SiO2/nk/Franta",
        extends: Some("silica"),
        summary: "Broadband fused silica EUV-MIR optical table (Franta 2016, 0.025-50 um)",
    },
    CatalogEntry {
        id: "aflow_si",
        domain: "phononics",
        format: ImportFormat::Aflow,
        target: "species(Si),nspecies(1)",
        extends: None,
        summary: "AFLOW 0 K DFT single-crystal Si stiffness + AGL thermal values",
    },
    CatalogEntry {
        id: "aflow_sic",
        domain: "phononics",
        format: ImportFormat::Aflow,
        target: "species(Si,C),nspecies(2)",
        extends: None,
        summary: "AFLOW 0 K DFT hexagonal SiC stiffness tensor + AGL thermal values",
    },
    CatalogEntry {
        id: "aflow_gan",
        domain: "phononics",
        format: ImportFormat::Aflow,
        target: "species(Ga,N),nspecies(2)",
        extends: None,
        summary: "AFLOW 0 K DFT GaN stiffness tensor + AGL thermal values",
    },
    CatalogEntry {
        id: "aflow_aln",
        domain: "phononics",
        format: ImportFormat::Aflow,
        target: "species(Al,N),nspecies(2)",
        extends: None,
        summary: "AFLOW 0 K DFT wurtzite AlN stiffness tensor + AGL thermal values",
    },
    CatalogEntry {
        id: "aflow_zno",
        domain: "phononics",
        format: ImportFormat::Aflow,
        target: "species(Zn,O),nspecies(2)",
        extends: None,
        summary: "AFLOW 0 K DFT ZnO stiffness tensor + AGL thermal values",
    },
    CatalogEntry {
        id: "aflow_ge",
        domain: "phononics",
        format: ImportFormat::Aflow,
        target: "species(Ge),nspecies(1)",
        extends: None,
        summary: "AFLOW 0 K DFT diamond-cubic Ge stiffness tensor + AGL thermal values",
    },
    CatalogEntry {
        id: "numbat_si_poulton",
        domain: "phononics",
        format: ImportFormat::Numbat,
        target: "Si_2021_Poulton",
        extends: None,
        summary: "NumBAT calibrated Si [100] elastic, photoelastic and acoustic viscosity notes",
    },
    CatalogEntry {
        id: "numbat_gaas_poulton",
        domain: "phononics",
        format: ImportFormat::Numbat,
        target: "GaAs_2021_Poulton",
        extends: None,
        summary: "NumBAT calibrated GaAs [100] elastic, photoelastic and acoustic viscosity notes",
    },
    CatalogEntry {
        id: "numbat_sio2_poulton",
        domain: "phononics",
        format: ImportFormat::Numbat,
        target: "SiO2_2021_Poulton",
        extends: None,
        summary: "NumBAT calibrated fused silica elastic, photoelastic and acoustic viscosity notes",
    },
    CatalogEntry {
        id: "numbat_as2s3_poulton",
        domain: "phononics",
        format: ImportFormat::Numbat,
        target: "As2S3_2021_Poulton",
        extends: None,
        summary: "NumBAT calibrated As2S3 chalcogenide glass for backward SBS waveguides",
    },
    CatalogEntry {
        id: "nb_weaver",
        domain: "superconducting",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/Nb/nk/Weaver",
        extends: Some("niobium_sc"),
        summary: "Niobium optical n,k 0.026-12.4 um (Weaver, Lynch & Olson 1973)",
    },
    CatalogEntry {
        id: "al_rakic_tabulated",
        domain: "superconducting",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/Al/nk/Rakic",
        extends: Some("aluminium"),
        summary: "Evaporated Al KK self-consistent table 0.12 pm - 200 um (Rakic 1995)",
    },
    CatalogEntry {
        id: "tin_pfluger",
        domain: "superconducting",
        format: ImportFormat::RefractiveIndexInfo,
        target: "main/TiN/nk/Pfluger",
        extends: Some("titanium_nitride_sc"),
        summary: "TiN1.0 thin-film optical constants from EELS/KK (Pfluger & Fink 1984/1991)",
    },
];

pub fn catalog_entry(id: &str) -> Option<&'static CatalogEntry> {
    let id = id.trim();
    REMOTE_CATALOG.iter().find(|c| c.id.eq_ignore_ascii_case(id))
}

fn err(path: &str, msg: impl Into<String>) -> MaterialError {
    MaterialError { material: "import".into(), path: path.into(), message: msg.into() }
}

fn encode_url_path(path: &str) -> String {
    let mut out = String::with_capacity(path.len() + 8);
    for ch in path.chars() {
        match ch {
            ' ' => out.push_str("%20"),
            '(' => out.push_str("%28"),
            ')' => out.push_str("%29"),
            c => out.push(c),
        }
    }
    out
}

/// Resolves a curated identifier, database coordinate or raw URL into an HTTPS download URL.
pub fn resolve_url(format: ImportFormat, target: &str) -> Result<String, MaterialError> {
    let t = target.trim();
    if t.is_empty() {
        return Err(err("target", "target must not be empty"));
    }
    if t.starts_with("https://") || t.starts_with("http://") {
        return Ok(t.to_string());
    }
    match format {
        ImportFormat::RefractiveIndexInfo => {
            let trimmed = t.trim_start_matches('/');
            let no_ext = trimmed.strip_suffix(".yml").or_else(|| trimmed.strip_suffix(".yaml")).unwrap_or(trimmed);
            let mut parts: Vec<&str> = no_ext.split('/').filter(|s| !s.is_empty()).collect();
            if parts.len() >= 3 {
                let is_specs = parts[0].eq_ignore_ascii_case("specs");
                let has_kind = parts.iter().any(|p| *p == "nk" || *p == "optical");
                if !is_specs && !has_kind {
                    let last = parts.pop().unwrap();
                    parts.push("nk");
                    parts.push(last);
                }
            }
            let rel = encode_url_path(&parts.join("/"));
            Ok(format!("https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/database/data/{rel}.yml"))
        }
        ImportFormat::Aflow => {
            let clean: String = t.chars().filter(|c| !c.is_whitespace()).collect();
            let mut query = clean;
            for clause in [
                "ael_stiffness_tensor(*)",
                "density(*)",
                "agl_thermal_conductivity_300K(*)",
                "agl_thermal_expansion_300K(*)",
            ] {
                let stem = clause.trim_end_matches("(*)");
                if !query.contains(stem) {
                    query.push(',');
                    query.push_str(clause);
                }
            }
            if !query.contains("paging(") {
                query.push_str(",paging(1,1)");
            }
            if !query.contains("format(") {
                query.push_str(",format(json)");
            }
            Ok(format!("https://aflow.org/API/aflux/?{query}"))
        }
        ImportFormat::Numbat => {
            let clean = t.trim_start_matches('/').strip_suffix(".json").unwrap_or(t.trim_start_matches('/'));
            Ok(format!("https://raw.githubusercontent.com/bjornsturmberg/NumBAT/master/backend/material_data/{clean}.json"))
        }
        ImportFormat::MaterialsProject => Err(err(
            "format",
            "Materials Project requires an API key header (X-API-KEY); download the summary/elasticity JSON from api.materialsproject.org locally or paste it into `import_material` with format = \"mp\"",
        )),
    }
}

/// Converts raw text from an external database format into validated `format = 1` material TOML text.
pub fn import_to_toml(format: ImportFormat, raw: &str, opts: &ImportOptions<'_>) -> Result<String, MaterialError> {
    if opts.name.trim().is_empty() || opts.name.contains(['/', '\\', '*', ',']) {
        return Err(err("name", "a material name is plain text without slashes, commas or stars"));
    }
    let toml_text = match format {
        ImportFormat::RefractiveIndexInfo => import_rii_yaml(raw, opts)?,
        ImportFormat::Aflow => import_aflow_json(raw, opts)?,
        ImportFormat::MaterialsProject => import_mp_json(raw, opts)?,
        ImportFormat::Numbat => import_numbat_json(raw, opts)?,
    };
    validate_converted(opts, &toml_text)?;
    Ok(toml_text)
}

fn validate_converted(opts: &ImportOptions<'_>, toml_text: &str) -> Result<(), MaterialError> {
    if opts.extends.is_some() {
        let mut table: toml::Table = toml_text.parse().map_err(|e: toml::de::Error| err("toml", e.to_string()))?;
        table.remove("extends");
        MaterialSpec::from_table(&table, opts.name)?;
    } else {
        MaterialSpec::from_toml_named(toml_text, opts.name)?;
    }
    Ok(())
}

/// Validates converted TOML together with an optional base `Library` (for `extends` checks in `nyx-engine`).
pub fn validate_against_library(lib: &Library, name: &str, toml_text: &str) -> Result<MaterialSpec, MaterialError> {
    lib.read(toml_text, name)
}

// ------------------------------------------------------------------ shared helpers

fn strip_html(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut in_tag = false;
    for ch in s.chars() {
        match ch {
            '<' => {
                if !out.ends_with(' ') {
                    out.push(' ');
                }
                in_tag = true;
            }
            '>' => in_tag = false,
            '\r' | '\n' if !in_tag => {
                if !out.ends_with(' ') {
                    out.push(' ');
                }
            }
            c if !in_tag => out.push(c),
            _ => {}
        }
    }
    out.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn extract_doi(s: &str) -> Option<String> {
    let mut rest = s;
    while let Some(idx) = rest.find("10.") {
        let cand = &rest[idx..];
        let end = cand
            .find(|c: char| c.is_whitespace() || matches!(c, '"' | '\'' | '<' | '>' | ')' | ']' | ',' | ';'))
            .unwrap_or(cand.len());
        let token = cand[..end].trim_end_matches('.');
        if let Some((_, rhs)) = token.split_once('/') {
            if !rhs.is_empty() && !token.contains("xx.xx") && !token.contains("XXXX") {
                return Some(token.to_string());
            }
        }
        rest = &cand[3..];
    }
    None
}

fn clean_http_url(s: &str) -> Option<String> {
    let s = s.trim();
    (s.starts_with("https://") || s.starts_with("http://")).then(|| s.to_string())
}

fn push_common_header(t: &mut toml::Table, opts: &ImportOptions<'_>, description: &str, notes: Option<&str>, source: Source) {
    t.insert("format".into(), toml::Value::Integer(FORMAT));
    t.insert("name".into(), opts.name.to_string().into());
    if let Some(ext) = opts.extends {
        t.insert("extends".into(), ext.to_string().into());
    }
    let desc = opts.description.map(|s| s.trim()).filter(|s| !s.is_empty()).unwrap_or(description);
    if !desc.is_empty() {
        t.insert("description".into(), desc.to_string().into());
    }
    if let Some(n) = notes.filter(|s| !s.trim().is_empty()) {
        t.insert("notes".into(), n.to_string().into());
    }
    t.insert("source".into(), source_to_toml(&source));
}

fn source_to_toml(s: &Source) -> toml::Value {
    let mut t = toml::Table::new();
    t.insert("cite".into(), s.cite.clone().into());
    if let Some(u) = &s.url {
        t.insert("url".into(), u.clone().into());
    }
    if let Some(d) = &s.doi {
        t.insert("doi".into(), d.clone().into());
    }
    if let Some(v) = &s.verified {
        t.insert("verified".into(), v.clone().into());
    }
    toml::Value::Table(t)
}

fn render_toml_table(t: &toml::Table) -> Result<String, MaterialError> {
    toml::to_string_pretty(t).map_err(|e| err("toml", e.to_string()))
}

fn build_spec_toml(spec: &MaterialSpec, extends: Option<&str>, notes: Option<&str>) -> Result<String, MaterialError> {
    let mut raw = spec.to_toml();
    if let Some(ext) = extends {
        raw = raw.replacen(&format!("name = \"{}\"\n", spec.name), &format!("name = \"{}\"\nextends = \"{ext}\"\n", spec.name), 1);
    }
    if let Some(n) = notes.filter(|s| !s.trim().is_empty()) {
        if spec.notes.is_none() {
            let mut t: toml::Table = raw.parse().map_err(|e: toml::de::Error| err("toml", e.to_string()))?;
            t.insert("notes".into(), n.to_string().into());
            if let Some(ext) = extends {
                t.insert("extends".into(), ext.to_string().into());
            }
            return render_toml_table(&t);
        }
    }
    Ok(raw)
}

// ------------------------------------------------------------------ 1. refractiveindex.info YAML

#[derive(Default)]
struct RiiBlock {
    kind: String,
    range: Option<(f64, f64)>,
    coeffs: Vec<f64>,
    rows: Vec<[f64; 3]>,
    has_n: bool,
    has_k: bool,
}

fn eval_rii_n(kind: &str, coeffs: &[f64], x: f64) -> Option<f64> {
    let x2 = x * x;
    match kind {
        "formula 1" => {
            if coeffs.is_empty() {
                return None;
            }
            let mut n2 = 1.0 + coeffs[0];
            for pair in coeffs[1..].as_chunks::<2>().0 {
                let d = x2 - pair[1] * pair[1];
                if d.abs() < 1e-18 {
                    return None;
                }
                n2 += pair[0] * x2 / d;
            }
            (n2 > 0.0).then(|| n2.sqrt())
        }
        "formula 2" => {
            if coeffs.is_empty() {
                return None;
            }
            let mut n2 = 1.0 + coeffs[0];
            for pair in coeffs[1..].as_chunks::<2>().0 {
                let d = x2 - pair[1];
                if d.abs() < 1e-18 {
                    return None;
                }
                n2 += pair[0] * x2 / d;
            }
            (n2 > 0.0).then(|| n2.sqrt())
        }
        "formula 3" => {
            if coeffs.is_empty() {
                return None;
            }
            let mut n2 = coeffs[0];
            for pair in coeffs[1..].as_chunks::<2>().0 {
                n2 += pair[0] * x.powf(pair[1]);
            }
            (n2 > 0.0).then(|| n2.sqrt())
        }
        "formula 4" => {
            // n^2 = c0 + c1*x^c2/(x^2 - c3^c4) + c5*x^c6/(x^2 - c7^c8) + c9*x^c10 + ...
            let mut n2 = *coeffs.first()?;
            let mut i = 1;
            for _ in 0..2 {
                if i + 3 < coeffs.len() {
                    if coeffs[i] != 0.0 {
                        let denom = x2 - coeffs[i + 2].powf(coeffs[i + 3]);
                        if denom.abs() < 1e-18 {
                            return None;
                        }
                        n2 += coeffs[i] * x.powf(coeffs[i + 1]) / denom;
                    }
                    i += 4;
                }
            }
            while i + 1 < coeffs.len() {
                n2 += coeffs[i] * x.powf(coeffs[i + 1]);
                i += 2;
            }
            (n2 > 0.0).then(|| n2.sqrt())
        }
        "formula 5" => {
            if coeffs.is_empty() {
                return None;
            }
            let mut n = coeffs[0];
            for pair in coeffs[1..].as_chunks::<2>().0 {
                n += pair[0] * x.powf(pair[1]);
            }
            (n > 0.0).then_some(n)
        }
        "formula 6" => {
            if coeffs.is_empty() {
                return None;
            }
            let mut n = 1.0 + coeffs[0];
            let inv = 1.0 / x2;
            for pair in coeffs[1..].as_chunks::<2>().0 {
                let d = pair[1] - inv;
                if d.abs() < 1e-18 {
                    return None;
                }
                n += pair[0] / d;
            }
            (n > 0.0).then_some(n)
        }
        "formula 7" => {
            let mut c = [0.0; 6];
            for (i, v) in coeffs.iter().take(6).enumerate() {
                c[i] = *v;
            }
            let l = 1.0 / (x2 - 0.028);
            let n = c[0] + c[1] * l + c[2] * l * l + c[3] * x2 + c[4] * x2 * x2 + c[5] * x2 * x2 * x2;
            (n > 0.0).then_some(n)
        }
        _ => None,
    }
}

fn rii_block_to_dispersion(block: &RiiBlock, range_out: &mut Option<(f64, f64)>) -> Result<Dispersion, MaterialError> {
    let k = block.kind.as_str();
    if k.starts_with("tabulated") {
        if block.rows.len() < 2 {
            return Err(err("optical.data", "refractiveindex.info tabulated block has fewer than 2 valid rows"));
        }
        let mut lambda = Vec::with_capacity(block.rows.len());
        let mut n = Vec::with_capacity(block.rows.len());
        let mut kv = Vec::with_capacity(block.rows.len());
        for r in &block.rows {
            if r[0] <= 0.0 || r[1] <= 0.0 || r[2] < 0.0 {
                continue;
            }
            if lambda.last().is_some_and(|prev: &f64| r[0] <= *prev) {
                continue;
            }
            lambda.push(r[0]);
            n.push(r[1]);
            kv.push(r[2]);
        }
        if lambda.len() < 2 {
            return Err(err("optical.data", "tabulated rows must have increasing positive wavelengths and positive n"));
        }
        *range_out = Some((lambda[0], *lambda.last().unwrap()));
        return Ok(Dispersion::Tabulated { lambda, n, k: kv, interpolation: Interp::Pchip, extrapolate: Extrapolate::Error });
    }

    let range = block.range.ok_or_else(|| err("optical.range", format!("RII `{k}` lacks wavelength_range")))?;
    *range_out = Some(range);
    let c = &block.coeffs;

    match k {
        "formula 1" => {
            let mut a = 1.0 + *c.first().ok_or_else(|| err("optical", "missing formula 1 coefficients"))?;
            let (mut b, mut l0) = (Vec::new(), Vec::new());
            for p in c[1..].as_chunks::<2>().0 {
                if p[1] == 0.0 {
                    a += p[0];
                } else {
                    b.push(p[0]);
                    l0.push(p[1] * p[1]);
                }
            }
            if !b.is_empty() && l0.iter().all(|ci| *ci > 0.0 && !(ci.sqrt() >= range.0 && ci.sqrt() <= range.1)) {
                return Ok(Dispersion::Sellmeier { a, b, c: l0 });
            }
            sample_formula_to_table(k, c, range)
        }
        "formula 2" => {
            let mut a = 1.0 + *c.first().ok_or_else(|| err("optical", "missing formula 2 coefficients"))?;
            let (mut b, mut cv) = (Vec::new(), Vec::new());
            for p in c[1..].as_chunks::<2>().0 {
                if p[1] == 0.0 {
                    a += p[0];
                } else {
                    b.push(p[0]);
                    cv.push(p[1]);
                }
            }
            if !b.is_empty() && cv.iter().all(|ci| *ci > 0.0 && !(ci.sqrt() >= range.0 && ci.sqrt() <= range.1)) {
                return Ok(Dispersion::Sellmeier { a, b, c: cv });
            }
            sample_formula_to_table(k, c, range)
        }
        "formula 3" => {
            // Check exact Schott exponents: c0, c1, +2, c3, -2, c5, -4, c7, -6, c9, -8
            if c.len() == 11 && [2.0, -2.0, -4.0, -6.0, -8.0].iter().enumerate().all(|(idx, exp)| (c[2 + 2 * idx] - exp).abs() < 1e-12) {
                return Ok(Dispersion::Schott { coefficients: vec![c[0], c[1], c[3], c[5], c[7], c[9]] });
            }
            sample_formula_to_table(k, c, range)
        }
        "formula 7" => {
            if c.is_empty() || c.len() > 6 {
                return Err(err("optical", "RII formula 7 (Herzberger) expects 1..6 coefficients"));
            }
            let mut out = vec![0.0; 6];
            for (i, v) in c.iter().enumerate() {
                out[i] = *v;
            }
            Ok(Dispersion::Herzberger { coefficients: out })
        }
        "formula 4" => {
            if c.is_empty() {
                return Err(err("optical", "missing RII formula 4 coefficients"));
            }
            let _ = eval_rii_n(k, c, 0.5 * (range.0 + range.1))
                .ok_or_else(|| err("optical", "RII formula 4 fails to evaluate cleanly in wavelength_range"))?;
            Ok(Dispersion::RiiFormula4 { coefficients: c.to_vec() })
        }
        "formula 5" => {
            if c.is_empty() || c.len() % 2 != 1 {
                return Err(err("optical", "RII formula 5 expects 1 + 2*m coefficients"));
            }
            Ok(Dispersion::RiiFormula5 { coefficients: c.to_vec() })
        }
        "formula 6" => {
            if c.is_empty() || c.len() % 2 != 1 {
                return Err(err("optical", "RII formula 6 expects 1 + 2*m coefficients"));
            }
            Ok(Dispersion::RiiFormula6 { coefficients: c.to_vec() })
        }
        other => Err(err("optical", format!("unsupported refractiveindex.info DATA type `{other}`"))),
    }
}

fn extract_k_data(k_block: &RiiBlock, n_range: Option<(f64, f64)>) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut kl = Vec::new();
    let mut kv = Vec::new();
    for r in &k_block.rows {
        let wl = r[0];
        let k = r[2];
        if wl <= 0.0 || k < 0.0 || !k.is_finite() {
            continue;
        }
        if let Some((lo, hi)) = n_range {
            if wl < lo * (1.0 - 1e-12) || wl > hi * (1.0 + 1e-12) {
                continue;
            }
        }
        if kl.last().is_some_and(|prev: &f64| wl <= *prev) {
            continue;
        }
        kl.push(wl);
        kv.push(k);
    }
    (kl.len() >= 2).then_some((kl, kv))
}

fn sample_formula_to_table(kind: &str, coeffs: &[f64], range: (f64, f64)) -> Result<Dispersion, MaterialError> {
    let n_pts = 81_usize;
    let (lo, hi) = range;
    if lo <= 0.0 || hi <= lo {
        return Err(err("optical.range", format!("invalid wavelength_range [{lo}, {hi}]")));
    }
    let mut lambda = Vec::with_capacity(n_pts);
    let mut n = Vec::with_capacity(n_pts);
    for i in 0..n_pts {
        let t = i as f64 / (n_pts - 1) as f64;
        let x = lo * (hi / lo).powf(t);
        let nv = eval_rii_n(kind, coeffs, x).ok_or_else(|| err("optical", format!("cannot evaluate `{kind}` at {x} um")))?;
        lambda.push(x);
        n.push(nv);
    }
    let k = vec![0.0; n_pts];
    Ok(Dispersion::Tabulated { lambda, n, k, interpolation: Interp::Pchip, extrapolate: Extrapolate::Error })
}

fn import_rii_yaml(raw: &str, opts: &ImportOptions<'_>) -> Result<String, MaterialError> {
    let mut section = "";
    let mut refs = String::new();
    let mut comments = String::new();
    let mut temp_k = 293.15_f64;
    let mut blocks: Vec<RiiBlock> = Vec::new();
    let mut in_block_data = false;

    for raw_line in raw.lines() {
        let stripped = raw_line.trim_end();
        let trimmed = stripped.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Top-level YAML keys have 0 leading indent (and do not start with '-')
        if !stripped.starts_with(char::is_whitespace) && !trimmed.starts_with('-') {
            in_block_data = false;
            if let Some((k, v)) = trimmed.split_once(':') {
                section = k.trim();
                let v = v.trim();
                if v != "|" && !v.is_empty() {
                    if section == "REFERENCES" {
                        refs.push_str(v);
                        refs.push(' ');
                    } else if section == "COMMENTS" {
                        comments.push_str(v);
                        comments.push(' ');
                    }
                }
            }
            continue;
        }

        match section {
            "REFERENCES" => {
                refs.push_str(trimmed);
                refs.push(' ');
            }
            "COMMENTS" => {
                comments.push_str(trimmed);
                comments.push(' ');
            }
            "CONDITIONS" => {
                if let Some((k, v)) = trimmed.split_once(':') {
                    if k.trim() == "temperature" {
                        let num_str = v.split_whitespace().next().unwrap_or("");
                        if let Ok(t) = num_str.parse::<f64>() {
                            if t > 0.0 {
                                temp_k = t;
                            }
                        }
                    }
                }
            }
            "DATA" => {
                if let Some(after) = trimmed.strip_prefix('-') {
                    blocks.push(RiiBlock::default());
                    in_block_data = false;
                    let after = after.trim();
                    if let Some((k, v)) = after.split_once(':') {
                        if k.trim() == "type" {
                            let cur = blocks.last_mut().unwrap();
                            cur.kind = v.trim().to_ascii_lowercase();
                            cur.has_n = cur.kind.contains("formula") || cur.kind.ends_with("nk") || cur.kind.ends_with(" n");
                            cur.has_k = cur.kind.ends_with("nk") || cur.kind.ends_with(" k");
                        }
                    }
                    continue;
                }
                let Some(cur) = blocks.last_mut() else { continue };
                if in_block_data && !trimmed.starts_with("type:") && !trimmed.starts_with("wavelength_range:") && !trimmed.starts_with("coefficients:") {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    let parse = |idx: usize| parts.get(idx).and_then(|s| s.parse::<f64>().ok());
                    if cur.kind == "tabulated nk" {
                        if let (Some(w), Some(n), Some(k)) = (parse(0), parse(1), parse(2)) {
                            cur.rows.push([w, n, k]);
                        }
                    } else if cur.kind == "tabulated n" {
                        if let (Some(w), Some(n)) = (parse(0), parse(1)) {
                            cur.rows.push([w, n, 0.0]);
                        }
                    } else if cur.kind == "tabulated k" {
                        if let (Some(w), Some(k)) = (parse(0), parse(1)) {
                            cur.rows.push([w, 1.0, k]);
                        }
                    }
                    continue;
                }
                if let Some((k, v)) = trimmed.split_once(':') {
                    in_block_data = false;
                    let key = k.trim();
                    let val = v.trim();
                    match key {
                        "type" => {
                            cur.kind = val.to_ascii_lowercase();
                            cur.has_n = cur.kind.contains("formula") || cur.kind.ends_with("nk") || cur.kind.ends_with(" n");
                            cur.has_k = cur.kind.ends_with("nk") || cur.kind.ends_with(" k");
                        }
                        "wavelength_range" => {
                            let nums: Vec<f64> = val.split_whitespace().filter_map(|s| s.parse().ok()).collect();
                            if nums.len() >= 2 {
                                cur.range = Some((nums[0], nums[1]));
                            }
                        }
                        "coefficients" => {
                            cur.coeffs = val.split_whitespace().filter_map(|s| s.parse().ok()).collect();
                        }
                        "data" => {
                            in_block_data = true;
                        }
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }

    if blocks.is_empty() {
        return Err(err("optical", "no DATA blocks found in refractiveindex.info YAML"));
    }

    let n_block = blocks.iter().find(|b| b.has_n).ok_or_else(|| err("optical", "refractiveindex.info file has no refractive index (`n` or `formula`) DATA block"))?;
    let k_only_block = blocks.iter().find(|b| !b.has_n && b.has_k);

    let mut range = None;
    let dispersion = rii_block_to_dispersion(n_block, &mut range)?;
    let k_data = k_only_block.and_then(|kb| extract_k_data(kb, range));

    let clean_ref = strip_html(&refs);
    let clean_cmt = strip_html(&comments);
    let doi = extract_doi(&refs);
    let url = opts.source_url.and_then(clean_http_url);
    let cite = if clean_ref.is_empty() {
        "UNVERIFIED automatic import from refractiveindex.info (CC0)".to_string()
    } else {
        format!("{clean_ref} — UNVERIFIED automatic import from refractiveindex.info")
    };
    let src = Source { cite, url, doi, verified: None };

    let mut desc = format!("Optical model {} imported from refractiveindex.info", n_block.kind);
    if !clean_cmt.is_empty() {
        desc.push_str(&format!(" ({})", clean_cmt.chars().take(110).collect::<String>()));
    }

    let spec = MaterialSpec {
        name: opts.name.to_string(),
        description: Some(opts.description.map(|s| s.to_string()).unwrap_or(desc)),
        source: vec![src.clone()],
        notes: (!clean_cmt.is_empty()).then_some(clean_cmt),
        optical: Some(Optical {
            source: vec![src],
            wavelength_unit: "um".into(),
            unit_m: 1e-6,
            range,
            reference_temperature: temp_k,
            axes: Axes::Isotropic(AxisModel {
                dispersion,
                source: Vec::new(),
                dn_dt: 0.0,
                k_data,
                dn_dt_einstein_theta: None,
                dn_dt_range: None,
            }),
            gyrotropy: None,
        }),
        ..Default::default()
    };
    build_spec_toml(&spec, opts.extends, None)
}

// ------------------------------------------------------------------ 2. shared elastic matrix importer

fn symmetrize_6x6_gpa(raw: &[[f64; 6]; 6], label: &str) -> Result<Elastic, MaterialError> {
    let max_val = raw.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs())).max(1e-12);
    let mut max_asym = 0.0_f64;
    let mut sym = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let d = (raw[i][j] - raw[j][i]).abs() / max_val;
            max_asym = max_asym.max(d);
            sym[i][j] = 0.5 * (raw[i][j] + raw[j][i]) * 1e9;
        }
    }
    if max_asym > 0.05 {
        return Err(err("elastic.matrix", format!("{label} stiffness tensor asymmetry exceeds 5% ({max_asym:.3})")));
    }
    Elastic::from_matrix(sym).map_err(|e| err("elastic.matrix", e))
}

fn parse_6x6_json(v: &Value, path: &str) -> Result<[[f64; 6]; 6], MaterialError> {
    let rows = v.as_array().ok_or_else(|| err(path, "expected a 6x6 array"))?;
    if rows.len() != 6 {
        return Err(err(path, format!("expected 6 rows, found {}", rows.len())));
    }
    let mut out = [[0.0; 6]; 6];
    for (i, row) in rows.iter().enumerate() {
        let cols = row.as_array().ok_or_else(|| err(path, format!("row {i} is not an array")))?;
        if cols.len() != 6 {
            return Err(err(path, format!("row {i} has {} entries, expected 6", cols.len())));
        }
        for (j, val) in cols.iter().enumerate() {
            out[i][j] = val.as_f64().ok_or_else(|| err(path, format!("non-numeric entry [{i}][{j}]")))?;
        }
    }
    Ok(out)
}

fn matrix_gpa_toml(c_pa: &M6) -> toml::Value {
    toml::Value::Array(
        c_pa.iter()
            .map(|r| {
                toml::Value::Array(
                    r.iter()
                        .map(|v| {
                            let gpa = ((v / 1e9) * 10000.0).round() / 10000.0;
                            toml::Value::Float(gpa)
                        })
                        .collect(),
                )
            })
            .collect(),
    )
}

// ------------------------------------------------------------------ 3. AFLOW JSON

fn pick_aflow_entry(v: &Value) -> Result<&Value, MaterialError> {
    if let Some(arr) = v.as_array() {
        return arr.first().ok_or_else(|| err("aflow", "empty AFLOW response array"));
    }
    if let Some(obj) = v.as_object() {
        if obj.contains_key("ael_stiffness_tensor") || obj.contains_key("compound") || obj.contains_key("auid") {
            return Ok(v);
        }
        if let Some((_, first)) = obj.iter().find(|(_, item)| item.is_object()) {
            return Ok(first);
        }
    }
    Err(err("aflow", "unrecognized AFLOW JSON structure (expected object, array or {\"1 of N\": {...}})"))
}

fn import_aflow_json(raw: &str, opts: &ImportOptions<'_>) -> Result<String, MaterialError> {
    let root: Value = serde_json::from_str(raw).map_err(|e| err("json", format!("invalid AFLOW JSON: {e}")))?;
    let entry = pick_aflow_entry(&root)?;

    let compound = entry.get("compound").and_then(Value::as_str).unwrap_or("unknown");
    let auid = entry.get("auid").and_then(Value::as_str).unwrap_or("aflow:unknown");
    let sg = entry.get("spacegroup_relax").and_then(Value::as_u64);

    let mut t = toml::Table::new();
    let cite = format!(
        "UNVERIFIED automatic 0 K DFT import from AFLOW AEL/AGL ({compound}, {auid}{}); single-crystal bulk calculation, do not use uncalibrated for room-temperature or amorphous/polycrystalline thin films",
        sg.map(|s| format!(", spacegroup {s}")).unwrap_or_default()
    );
    let src = Source {
        cite: cite.clone(),
        url: opts.source_url.and_then(clean_http_url).or_else(|| Some("https://aflow.org/".into())),
        doi: Some("10.1016/j.commatsci.2012.02.005".into()),
        verified: None,
    };
    let desc = format!("AFLOW 0 K DFT bulk {compound} ({auid}) — uncalibrated single-crystal reference; check thin-film/temperature offsets.");
    let notes = "0 K DFT AEL/AGL calculation (Volumetric thermal expansion divided by 3 for linear isotropic coefficient; thin-film stress, grain boundaries and finite temperatures are excluded).";
    push_common_header(&mut t, opts, &desc, Some(notes), src.clone());

    if let Some(rho_gcc) = entry.get("density").and_then(Value::as_f64).filter(|d| *d > 0.0) {
        t.insert("density".into(), format!("{:.2} kg/m^3", rho_gcc * 1000.0).into());
    }

    if let Some(tensor) = entry.get("ael_stiffness_tensor").filter(|v| !v.is_null()) {
        let raw6 = parse_6x6_json(tensor, "elastic.matrix")?;
        let el = symmetrize_6x6_gpa(&raw6, "AFLOW AEL")?;
        let mut e_tab = toml::Table::new();
        e_tab.insert("unit".into(), "GPa".into());
        e_tab.insert("matrix".into(), matrix_gpa_toml(&el.c));
        e_tab.insert("source".into(), source_to_toml(&src));
        t.insert("elastic".into(), toml::Value::Table(e_tab));
    } else if opts.extends.is_none() {
        return Err(err("ael_stiffness_tensor", "AFLOW entry has no `ael_stiffness_tensor`"));
    }

    let k_300 = entry.get("agl_thermal_conductivity_300K").and_then(Value::as_f64).filter(|x| *x > 0.0);
    let a_300 = entry.get("agl_thermal_expansion_300K").and_then(Value::as_f64);
    if k_300.is_some() || a_300.is_some() {
        let mut th = toml::Table::new();
        if let Some(k) = k_300 {
            th.insert("conductivity".into(), format!("{k:.4} W/(m*K)").into());
        }
        if let Some(a_vol) = a_300 {
            th.insert("expansion".into(), format!("{:.6e} 1/K", a_vol / 3.0).into());
        }
        th.insert("source".into(), source_to_toml(&src));
        t.insert("thermal".into(), toml::Value::Table(th));
    }

    render_toml_table(&t)
}

// ------------------------------------------------------------------ 4. Materials Project JSON

fn mp_3x3_toml(arr: &Value) -> Option<toml::Value> {
    let diel = arr.as_array()?;
    if diel.len() != 3 || !diel.iter().all(|r| r.as_array().is_some_and(|c| c.len() == 3)) {
        return None;
    }
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = diel[i][j].as_f64()?;
        }
    }
    for i in 0..3 {
        for j in (i + 1)..3 {
            let avg = 0.5 * (m[i][j] + m[j][i]);
            m[i][j] = avg;
            m[j][i] = avg;
        }
    }
    let pd = TensorProp::Symmetric([m[0][0], m[1][1], m[2][2], m[1][2], m[0][2], m[0][1]])
        .principal_values()
        .into_iter()
        .all(|v| v > 0.0);
    if !pd {
        return None;
    }
    let off_zero = m[0][1].abs() <= 1e-8 && m[0][2].abs() <= 1e-8 && m[1][2].abs() <= 1e-8;
    let iso_xx_yy = (m[0][0] - m[1][1]).abs() <= 1e-6 * m[0][0].abs();
    let iso_xx_zz = (m[0][0] - m[2][2]).abs() <= 1e-6 * m[0][0].abs();
    Some(if off_zero && iso_xx_yy && iso_xx_zz {
        let s = if m[0][0] == m[1][1] && m[0][0] == m[2][2] {
            m[0][0]
        } else {
            (m[0][0] + m[1][1] + m[2][2]) / 3.0
        };
        toml::Value::Float(s)
    } else if off_zero && iso_xx_yy {
        let mut uni = toml::Table::new();
        let perp = if m[0][0] == m[1][1] { m[0][0] } else { 0.5 * (m[0][0] + m[1][1]) };
        uni.insert("perp".into(), toml::Value::Float(perp));
        uni.insert("par".into(), toml::Value::Float(m[2][2]));
        toml::Value::Table(uni)
    } else if off_zero {
        toml::Value::Array(vec![
            toml::Value::Float(m[0][0]),
            toml::Value::Float(m[1][1]),
            toml::Value::Float(m[2][2]),
        ])
    } else {
        toml::Value::Array(m.iter().map(|r| toml::Value::Array(r.iter().map(|v| toml::Value::Float(*v)).collect())).collect())
    })
}

fn import_mp_json(raw: &str, opts: &ImportOptions<'_>) -> Result<String, MaterialError> {
    let root: Value = serde_json::from_str(raw).map_err(|e| err("json", format!("invalid Materials Project JSON: {e}")))?;
    let doc = root
        .get("data")
        .and_then(Value::as_array)
        .and_then(|a| a.first())
        .or_else(|| root.as_array().and_then(|a| a.first()))
        .unwrap_or(&root);

    let mp_id = doc.get("material_id").and_then(Value::as_str).unwrap_or("mp-unknown");
    let formula = doc.get("formula_pretty").and_then(Value::as_str).unwrap_or(opts.name);

    let cite = format!("UNVERIFIED automatic 0 K DFT import from Materials Project ({formula}, {mp_id}); single-crystal PBE IEEE Voigt tensor, not calibrated for thin films or room-temperature anharmonicity");
    let src = Source {
        cite,
        url: opts.source_url.and_then(clean_http_url).or_else(|| Some(format!("https://next-gen.materialsproject.org/materials/{mp_id}"))),
        doi: Some("10.1038/sdata.2015.9".into()),
        verified: None,
    };

    let mut t = toml::Table::new();
    let desc = format!("Materials Project 0 K DFT {formula} ({mp_id}) — IEEE Voigt single-crystal tensor; calibrate for film/T.");
    let notes = "Exported in IEEE Std 176 Voigt order (xx, yy, zz, yz, xz, xy), matching nyx format = 2 directly.";
    push_common_header(&mut t, opts, &desc, Some(notes), src.clone());

    let pg = doc
        .get("symmetry")
        .and_then(|s| s.get("point_group"))
        .and_then(Value::as_str)
        .and_then(super::state_t::PointGroup::parse);
    if let Some(p) = pg {
        t.insert("point_group".into(), p.as_str().into());
    }

    if let Some(rho) = doc.get("density").and_then(Value::as_f64).filter(|d| *d > 0.0) {
        t.insert("density".into(), format!("{:.2} kg/m^3", rho * 1000.0).into());
    }

    let el_node = doc.get("elasticity").unwrap_or(doc);
    let tensor_val = el_node
        .get("elastic_tensor")
        .and_then(|et| et.get("ieee_format").or_else(|| et.get("raw")).or(Some(et)))
        .filter(|v| v.is_array())
        .or_else(|| doc.get("homogeneous_elastic_tensor").filter(|v| v.is_array()));

    if let Some(tv) = tensor_val {
        let raw6 = parse_6x6_json(tv, "elastic.matrix")?;
        let el = symmetrize_6x6_gpa(&raw6, "Materials Project")?;
        let mut e_tab = toml::Table::new();
        e_tab.insert("unit".into(), "GPa".into());
        e_tab.insert("matrix".into(), matrix_gpa_toml(&el.c));
        e_tab.insert("source".into(), source_to_toml(&src));
        t.insert("elastic".into(), toml::Value::Table(e_tab));
    }

    let diel_node = doc.get("dielectric");
    let e_total = diel_node.and_then(|d| d.get("e_total")).and_then(mp_3x3_toml);
    let e_elec = diel_node.and_then(|d| d.get("e_electronic")).and_then(mp_3x3_toml);
    let dfpt_src = Source {
        cite: format!("UNVERIFIED DFPT dielectric tensor from Materials Project ({mp_id})"),
        url: src.url.clone(),
        doi: Some("10.1038/sdata.2016.134".into()),
        verified: None,
    };
    if let Some(perm_val) = e_total.clone() {
        let mut el_tab = toml::Table::new();
        el_tab.insert("permittivity".into(), perm_val);
        el_tab.insert("source".into(), source_to_toml(&dfpt_src));
        t.insert("electrical".into(), toml::Value::Table(el_tab));
    }
    if e_total.is_some() || e_elec.is_some() {
        let mut dl_tab = toml::Table::new();
        if let Some(et) = e_total {
            dl_tab.insert("eps_clamped".into(), et);
        }
        if let Some(ee) = e_elec {
            dl_tab.insert("eps_inf".into(), ee);
        }
        dl_tab.insert("source".into(), source_to_toml(&dfpt_src));
        t.insert("dielectric".into(), toml::Value::Table(dl_tab));
    }

    if !pg.is_some_and(|p| p.forbids_piezo_and_pockels()) {
        if let Some(pz_arr) = doc
            .get("piezoelectric")
            .and_then(|p| p.get("total").or_else(|| p.get("e_ij")).or(Some(p)))
            .and_then(Value::as_array)
        {
            if pz_arr.len() == 3 && pz_arr.iter().all(|r| r.as_array().is_some_and(|c| c.len() == 6)) {
                let rows: Vec<toml::Value> = pz_arr
                    .iter()
                    .map(|r| {
                        toml::Value::Array(
                            r.as_array()
                                .unwrap()
                                .iter()
                                .map(|v| toml::Value::Float(v.as_f64().unwrap_or(0.0)))
                                .collect(),
                        )
                    })
                    .collect();
                let mut pz_tab = toml::Table::new();
                pz_tab.insert("form".into(), "stress_charge".into());
                pz_tab.insert("unit".into(), "C/m^2".into());
                pz_tab.insert("matrix".into(), toml::Value::Array(rows));
                pz_tab.insert("source".into(), source_to_toml(&src));
                t.insert("piezoelectric".into(), toml::Value::Table(pz_tab));
            }
        }
    }

    if !t.contains_key("elastic") && !t.contains_key("electrical") && !t.contains_key("dielectric") {
        return Err(err("mp", "Materials Project JSON contained neither `elastic_tensor` nor `dielectric.e_total`"));
    }
    render_toml_table(&t)
}

// ------------------------------------------------------------------ 5. NumBAT JSONC

/// Strips `//` line comments outside string literals so NumBAT `.json` files parse with standard `serde_json`.
pub fn strip_jsonc_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut in_str = false;
    let mut esc = false;
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if in_str {
            out.push(c);
            if esc {
                esc = false;
            } else if c == '\\' {
                esc = true;
            } else if c == '"' {
                in_str = false;
            }
            continue;
        }
        if c == '"' {
            in_str = true;
            out.push(c);
        } else if c == '/' && chars.peek() == Some(&'/') {
            for nc in chars.by_ref() {
                if nc == '\n' {
                    out.push('\n');
                    break;
                }
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn import_numbat_json(raw: &str, opts: &ImportOptions<'_>) -> Result<String, MaterialError> {
    let clean = strip_jsonc_comments(raw);
    let doc: Value = serde_json::from_str(&clean).map_err(|e| err("json", format!("invalid NumBAT JSON: {e}")))?;

    let chem = doc.get("chemical").and_then(Value::as_str).unwrap_or(opts.name);
    let author = doc.get("author").and_then(Value::as_str).unwrap_or("NumBAT");
    let date = doc.get("date").and_then(Value::as_str).unwrap_or("");
    let raw_doi = doc.get("doi").and_then(Value::as_str).unwrap_or("");
    let doi = extract_doi(raw_doi);

    let src = Source {
        cite: format!("UNVERIFIED automatic import from NumBAT material_data ({chem}, {author} {date})"),
        url: opts.source_url.and_then(clean_http_url).or_else(|| Some("https://github.com/bjornsturmberg/NumBAT".into())),
        doi,
        verified: None,
    };

    let get_opt = |k: &str| doc.get(k).and_then(Value::as_f64);
    let mut notes_parts: Vec<String> = vec![format!("Imported from NumBAT ({chem}).")];
    let mut pe: Vec<String> = vec![];
    for k in ["p_11", "p_12", "p_13", "p_14", "p_31", "p_33", "p_44", "p_66"] {
        if let Some(v) = get_opt(k) {
            pe.push(format!("{k}={v}"));
        }
    }
    if !pe.is_empty() {
        notes_parts.push(format!("Photoelastic: {}.", pe.join(", ")));
    }
    let mut visc: Vec<String> = vec![];
    for k in ["eta_11", "eta_12", "eta_44"] {
        if let Some(v) = get_opt(k) {
            visc.push(format!("{k}={v} Pa*s"));
        }
    }
    if !visc.is_empty() {
        notes_parts.push(format!("Acoustic viscosity: {}.", visc.join(", ")));
    }

    let mut t = toml::Table::new();
    let desc = format!("NumBAT calibrated Brillouin dataset ({chem}, {author} {date}) — single-wavelength n & bulk stiffness.");
    push_common_header(&mut t, opts, &desc, Some(&notes_parts.join(" ")), src.clone());

    if let Some(rho) = get_opt("s").filter(|d| *d > 0.0) {
        t.insert("density".into(), format!("{rho} kg/m^3").into());
    }

    if let (Some(c11), Some(c12), Some(c44)) = (get_opt("c_11"), get_opt("c_12"), get_opt("c_44")) {
        let mut m = [[0.0; 6]; 6];
        let c22 = get_opt("c_22").unwrap_or(c11);
        let c33 = get_opt("c_33").unwrap_or(c11);
        let c13 = get_opt("c_13").unwrap_or(c12);
        let c23 = get_opt("c_23").unwrap_or(c13);
        let c55 = get_opt("c_55").unwrap_or(c44);
        let c66 = get_opt("c_66").unwrap_or(if c13 == c12 && c33 == c11 { c44 } else { 0.5 * (c11 - c12) });
        m[0][0] = c11;
        m[1][1] = c22;
        m[2][2] = c33;
        m[0][1] = c12;
        m[1][0] = c12;
        m[0][2] = c13;
        m[2][0] = c13;
        m[1][2] = c23;
        m[2][1] = c23;
        m[3][3] = c44;
        m[4][4] = c55;
        m[5][5] = c66;
        if let Some(c14) = get_opt("c_14") {
            m[0][3] = c14;
            m[3][0] = c14;
            m[1][3] = -c14;
            m[3][1] = -c14;
            m[4][5] = c14;
            m[5][4] = c14;
        }
        // Convert Pa -> GPa view check and construct positive-definite Elastic
        Elastic::from_matrix(m).map_err(|e| err("elastic.matrix", e))?;
        let mut e_tab = toml::Table::new();
        e_tab.insert("unit".into(), "GPa".into());
        e_tab.insert("matrix".into(), matrix_gpa_toml(&m));
        e_tab.insert("source".into(), source_to_toml(&src));
        t.insert("elastic".into(), toml::Value::Table(e_tab));
    }

    if let Some(re_n) = get_opt("Re_n").filter(|x| *x > 0.0) {
        let im_k = get_opt("Im_n").unwrap_or(0.0).max(0.0);
        let mut o_tab = toml::Table::new();
        o_tab.insert("model".into(), "constant".into());
        o_tab.insert("n".into(), toml::Value::Float(re_n));
        if im_k > 0.0 {
            o_tab.insert("k".into(), toml::Value::Float(im_k));
        }
        o_tab.insert("source".into(), source_to_toml(&src));
        t.insert("optical".into(), toml::Value::Table(o_tab));
    }

    if let Some(p11) = get_opt("p_11") {
        let p12 = get_opt("p_12").unwrap_or(0.0);
        let p13 = get_opt("p_13").unwrap_or(p12);
        let p31 = get_opt("p_31").unwrap_or(p12);
        let p33 = get_opt("p_33").unwrap_or(p11);
        let p44 = get_opt("p_44").unwrap_or(0.5 * (p11 - p12));
        let p66 = get_opt("p_66").unwrap_or(0.5 * (p11 - p12));
        let p14 = get_opt("p_14").unwrap_or(0.0);
        let p41 = get_opt("p_41").unwrap_or(p14);
        let mut pm = [
            [p11, p12, p13, p14, 0.0, 0.0],
            [p12, p11, p13, -p14, 0.0, 0.0],
            [p31, p31, p33, 0.0, 0.0, 0.0],
            [p41, -p41, 0.0, p44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, p44, p41],
            [0.0, 0.0, 0.0, 0.0, p14, p66],
        ];
        if p14 == 0.0 && p41 == 0.0 && (p31 - p12).abs() < 1e-12 && (p33 - p11).abs() < 1e-12 && p66 == p44 {
            // Cubic or isotropic exact zeros
            pm[3][0] = 0.0;
        }
        let rows = pm
            .iter()
            .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
            .collect();
        let mut pe_tab = toml::Table::new();
        pe_tab.insert("matrix".into(), toml::Value::Array(rows));
        pe_tab.insert("source".into(), source_to_toml(&src));
        t.insert("photoelastic".into(), toml::Value::Table(pe_tab));
    }

    if let Some(eta11) = get_opt("eta_11").filter(|&v| v > 0.0) {
        let eta12 = get_opt("eta_12").unwrap_or(0.0);
        let eta13 = get_opt("eta_13").unwrap_or(eta12);
        let eta33 = get_opt("eta_33").unwrap_or(eta11);
        let eta44 = get_opt("eta_44").unwrap_or((0.5 * (eta11 - eta12)).max(0.0));
        let eta66 = get_opt("eta_66").unwrap_or((0.5 * (eta11 - eta12)).max(0.0));
        let em = [
            [eta11, eta12, eta13, 0.0, 0.0, 0.0],
            [eta12, eta11, eta13, 0.0, 0.0, 0.0],
            [eta13, eta13, eta33, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, eta44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, eta44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, eta66],
        ];
        let evals = super::elastic::sym_eigenvalues(em);
        if evals[0] >= -1e-12 * eta11 {
            let rows = em
                .iter()
                .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
                .collect();
            let mut ve_tab = toml::Table::new();
            ve_tab.insert("unit".into(), "Pa*s".into());
            ve_tab.insert("matrix".into(), toml::Value::Array(rows));
            ve_tab.insert("source".into(), source_to_toml(&src));
            t.insert("viscoelastic".into(), toml::Value::Table(ve_tab));
        }
    }

    render_toml_table(&t)
}
