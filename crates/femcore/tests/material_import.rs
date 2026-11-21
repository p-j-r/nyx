// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Offline verification of the pure-Rust importers (`femcore::material::import`) for refractiveindex.info, AFLOW, Materials Project and NumBAT.

use femcore::material::import::{catalog_entry, import_to_toml, resolve_url, strip_jsonc_comments, ImportFormat, ImportOptions, REMOTE_CATALOG};
use femcore::material::{Dispersion, Library, MaterialSpec, PointGroup};
use std::collections::BTreeSet;
use std::path::Path;

fn lib() -> Library {
    Library::load_dir(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials")).unwrap()
}

#[test]
fn url_resolution_handles_shorthands_specs_and_queries() {
    assert_eq!(
        resolve_url(ImportFormat::RefractiveIndexInfo, "main/Si/Salzberg").unwrap(),
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/database/data/main/Si/nk/Salzberg.yml"
    );
    assert_eq!(
        resolve_url(ImportFormat::RefractiveIndexInfo, "other/doped crystals/MgO-LiNbO3/Zelmon-o").unwrap(),
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/database/data/other/doped%20crystals/MgO-LiNbO3/nk/Zelmon-o.yml"
    );
    assert_eq!(
        resolve_url(ImportFormat::RefractiveIndexInfo, "specs/schott/optical/N-BK7").unwrap(),
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/database/data/specs/schott/optical/N-BK7.yml"
    );
    assert_eq!(
        resolve_url(ImportFormat::Aflow, "species(Si,C),nspecies(2)").unwrap(),
        "https://aflow.org/API/aflux/?species(Si,C),nspecies(2),ael_stiffness_tensor(*),density(*),agl_thermal_conductivity_300K(*),agl_thermal_expansion_300K(*),paging(1,1),format(json)"
    );
    assert_eq!(
        resolve_url(ImportFormat::Numbat, "Si_2021_Poulton").unwrap(),
        "https://raw.githubusercontent.com/bjornsturmberg/NumBAT/master/backend/material_data/Si_2021_Poulton.json"
    );
    assert!(resolve_url(ImportFormat::MaterialsProject, "mp-149").is_err());
}

#[test]
fn remote_catalog_entries_have_unique_ids_valid_urls_and_known_bases() {
    let base_lib = lib();
    let mut ids = BTreeSet::new();
    for entry in REMOTE_CATALOG {
        assert!(ids.insert(entry.id), "duplicate catalog id {}", entry.id);
        assert!(catalog_entry(entry.id).is_some());
        assert!(matches!(entry.domain, "photonics" | "phononics" | "superconducting"));
        let u = resolve_url(entry.format, entry.target).unwrap();
        assert!(u.starts_with("https://"), "{u}");
        if let Some(ext) = entry.extends {
            assert!(base_lib.get(ext).is_some(), "{} extends missing base {ext}", entry.id);
        }
    }
    assert!(REMOTE_CATALOG.len() >= 20);
}

#[test]
fn rii_yaml_imports_formula1_with_c0_herzberger_padding_and_out_of_order_conditions() {
    // 1) formula 1 with nonzero c0 (GaAs Skauli)
    let skauli = r#"
REFERENCES: |
    T. Skauli et al. <a href="https://doi.org/10.1063/1.1621740"><i>J. Appl. Phys.</i></a>
DATA:
  - type: formula 1
    wavelength_range: 0.97 17
    coefficients: 4.372514 5.466742 0.4431307 0.02429960 0.8746453 1.957522 36.9166
"#;
    let out = import_to_toml(ImportFormat::RefractiveIndexInfo, skauli, &ImportOptions { name: "gaas_test", description: None, extends: None, source_url: None }).unwrap();
    let spec = MaterialSpec::from_toml(&out).unwrap();
    assert!((spec.optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0 - 3.3702).abs() < 1e-3);
    assert!(spec.source[0].flags_uncertainty() && spec.source[0].doi.as_deref() == Some("10.1063/1.1621740"));

    // 2) formula 7 (Edwards 5 coefficients -> padded to 6) + CONDITIONS after DATA + extends = "silicon"
    let edwards = r#"
REFERENCES: |
    D. F. Edwards, <a href="https://doi.org/10.1364/AO.19.004130">Appl. Opt.</a>
DATA:
  - type: formula 7
    wavelength_range: 2.4373 25
    coefficients: 3.41983 0.159906 -0.123109 1.26878E-6 -1.95104E-9
CONDITIONS:
    temperature: 299
"#;
    let out = import_to_toml(
        ImportFormat::RefractiveIndexInfo,
        edwards,
        &ImportOptions { name: "si_edwards", description: None, extends: Some("silicon"), source_url: None },
    )
    .unwrap();
    let merged = lib().read(&out, "si_edwards").unwrap();
    let n5 = merged.optical.as_ref().unwrap().nk(0, 5.0e-6, None).unwrap().0;
    assert!((n5 - 3.4223).abs() < 5e-3, "{n5}");
    assert_eq!(merged.optical.as_ref().unwrap().reference_temperature, 299.0);
}

#[test]
fn aflow_afux_keyed_json_symmetrizes_and_warns_about_0k_dft() {
    let aflow = r#"{
        "1 of 13": {
            "compound": "C6Si6",
            "auid": "aflow:059e25f57cc58053",
            "spacegroup_relax": 186,
            "density": 3.17182,
            "agl_thermal_conductivity_300K": 34.5239,
            "agl_thermal_expansion_300K": 1.93848e-5,
            "ael_stiffness_tensor": [
                [488.417, 101.2, 50.5333, 0, 0, 0],
                [100.867, 488.4, 50.5333, 0, 0, 0],
                [50.4, 50.5, 533.133, 0, 0, 0],
                [0, 0, 0, 160.4, 0, 0],
                [0, 0, 0, 0, 160.4, 0],
                [0, 0, 0, 0, 0, 192.9]
            ]
        }
    }"#;
    let out = import_to_toml(ImportFormat::Aflow, aflow, &ImportOptions { name: "sic_aflow", description: None, extends: None, source_url: None }).unwrap();
    let spec = MaterialSpec::from_toml(&out).unwrap();
    assert!(spec.source[0].flags_uncertainty());
    assert!(spec.notes.as_deref().unwrap().contains("0 K DFT"));
    assert!((spec.density.unwrap() - 3171.82).abs() < 0.1);
    assert!((spec.elastic.as_ref().unwrap().c[2][2] - 533.133e9).abs() < 1e7);
    assert!((spec.legacy_scalars()["alpha"] - 1.93848e-5 / 3.0).abs() < 1e-8);
}

#[test]
fn materials_project_json_imports_ieee_elasticity_and_dielectric() {
    let mp = r#"{
        "data": [{
            "material_id": "mp-149",
            "formula_pretty": "Si",
            "density": 2.329,
            "elasticity": {
                "elastic_tensor": {
                    "ieee_format": [
                        [159.0, 62.0, 62.0, 0, 0, 0],
                        [62.0, 159.0, 62.0, 0, 0, 0],
                        [62.0, 62.0, 159.0, 0, 0, 0],
                        [0, 0, 0, 76.0, 0, 0],
                        [0, 0, 0, 0, 76.0, 0],
                        [0, 0, 0, 0, 0, 76.0]
                    ]
                }
            },
            "dielectric": {
                "e_total": [[12.1, 0, 0], [0, 12.1, 0], [0, 0, 12.1]]
            }
        }]
    }"#;
    let out = import_to_toml(ImportFormat::MaterialsProject, mp, &ImportOptions { name: "si_mp", description: None, extends: None, source_url: None }).unwrap();
    let spec = MaterialSpec::from_toml(&out).unwrap();
    assert!((spec.elastic.unwrap().c[0][0] - 159e9).abs() < 1e6);
    assert_eq!(spec.electrical.permittivity.unwrap().as_scalar(), Some(12.1));
}

#[test]
fn numbat_jsonc_strips_comments_and_ignores_placeholder_dois() {
    let numbat = r#"{
        "file_name": "Si_2021_Poulton", // Name
        "chemical": "Si [100]", // Composition
        "author": "Poulton",
        "date": "2021",
        "doi": "dx.doi.org/xx.xxxx/xx.xx.xxxxxx", // placeholder!
        "Re_n": 3.4757,
        "Im_n": 0.0,
        "s": 2331,
        "c_11": 165.8e9,
        "c_12": 63.9e9,
        "c_44": 79.6e9,
        "p_11": -0.094,
        "eta_11": 5.9e-3
    }"#;
    assert!(!strip_jsonc_comments(numbat).contains("// Name"));
    let out = import_to_toml(ImportFormat::Numbat, numbat, &ImportOptions { name: "si_numbat", description: None, extends: None, source_url: None }).unwrap();
    let spec = MaterialSpec::from_toml(&out).unwrap();
    assert!(spec.source[0].doi.is_none(), "placeholder xx.xxxx must not become a DOI");
    assert_eq!(spec.density, Some(2331.0));
    assert!(spec.notes.as_deref().unwrap().contains("p_11=-0.094"));
    assert!(spec.photoelastic.is_some());
    assert!((spec.viscoelastic.as_ref().unwrap().matrix[0][0] - 5.9e-3).abs() < 1e-12);
    assert!((spec.optical.unwrap().nk(0, 1.55e-6, None).unwrap().0 - 3.4757).abs() < 1e-6);
}

#[test]
fn format2_importers_preserve_rii_formula4_kdata_and_mp_piezo() {
    let rii = r#"
REFERENCES: "Test Zelmon + k"
DATA:
  - type: formula 4
    wavelength_range: 0.4 4.0
    coefficients: 4.9048 0.11768 0 0.0475 1 0 0 1 1 -0.027169 2
  - type: tabulated k
    data: |
      0.40 0.002
      1.00 0.0001
      2.00 0.0
"#;
    let out = import_to_toml(ImportFormat::RefractiveIndexInfo, rii, &ImportOptions { name: "ln_rii", description: None, extends: None, source_url: None }).unwrap();
    let spec = MaterialSpec::from_toml(&out).unwrap();
    assert!(out.starts_with("format = 2") || out.contains("\nformat = 2"));
    let opt = spec.optical.as_ref().unwrap();
    let (n, k, _, _) = opt.nk(0, 1.0e-6, None).unwrap();
    assert!(n > 2.2 && (k - 1e-4).abs() < 1e-12);
    assert!(matches!(opt.axes.model(0).dispersion, Dispersion::RiiFormula4 { .. }));

    let mp = r#"{
        "data": [{
            "material_id": "mp-aln",
            "formula_pretty": "AlN",
            "symmetry": { "point_group": "6mm" },
            "density": 3.26,
            "elasticity": { "elastic_tensor": { "ieee_format": [[410,149,99,0,0,0],[149,410,99,0,0,0],[99,99,389,0,0,0],[0,0,0,125,0,0],[0,0,0,0,125,0],[0,0,0,0,0,130.5]] } },
            "dielectric": { "e_total": [[8.5,0,0],[0,8.5,0],[0,0,9.5]], "e_electronic": [[4.6,0,0],[0,4.6,0],[0,0,4.8]] },
            "piezoelectric": { "total": [[0,0,0,0,-0.48,0],[0,0,0,-0.48,0,0],[-0.58,-0.58,1.55,0,0,0]] }
        }]
    }"#;
    let mp_toml = import_to_toml(ImportFormat::MaterialsProject, mp, &ImportOptions { name: "aln_mp", description: None, extends: None, source_url: None }).unwrap();
    let mp_spec = MaterialSpec::from_toml(&mp_toml).unwrap();
    assert_eq!(mp_spec.point_group, Some(PointGroup::C6v));
    assert_eq!(mp_spec.dielectric.as_ref().unwrap().eps_inf.as_ref().unwrap().matrix()[2][2], 4.8);
    assert!((mp_spec.piezoelectric.as_ref().unwrap().matrix[2][2] - 1.55).abs() < 1e-12);
}
