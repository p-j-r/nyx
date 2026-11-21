// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! `extends`: a material that starts from another of the library and replaces whole tables. The combined result is validated.

use femcore::material::*;
use std::path::Path;

fn root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials")
}

fn library() -> Library {
    Library::load_dir(&root()).unwrap_or_else(|e| panic!("{}", e.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\n")))
}

fn temp_dir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-extends-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn the_shipped_example_replaces_only_the_optical_table_and_the_source() {
    let lib = library();
    let text = std::fs::read_to_string(root().join("examples/silicon_measured.toml")).unwrap();
    let m = lib.read(&text, "silicon_measured").unwrap();
    let base = lib.get("silicon").unwrap();
    assert_eq!(m.name, "silicon_measured", "the name is not inherited");
    assert!(m.source[0].cite.contains("Own ellipsometry"));
    assert_eq!(m.density, base.density);
    assert_eq!(m.elastic, base.elastic);
    assert_eq!(m.thermal, base.thermal);
    assert_eq!(m.electrical, base.electrical);
    let o = m.optical.as_ref().unwrap();
    assert_eq!(o.nk(0, 1.55e-6, None).unwrap().0, 3.4761);
    assert!(o.nk(0, 1.3e-6, None).is_err(), "the measured range is 1.5 to 1.6 um");
    assert_eq!(base.optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0, 3.4757, "the base is untouched");
}

#[test]
fn library_files_may_extend_each_other_in_a_chain() {
    let d = temp_dir("chain");
    std::fs::write(d.join("a.toml"), "format = 1\ndensity = \"1000 kg/m^3\"\n[electrical]\npermittivity = 2.0\n").unwrap();
    std::fs::write(d.join("b.toml"), "format = 1\nextends = \"a\"\ndescription = \"b\"\n").unwrap();
    std::fs::write(d.join("c.toml"), "format = 1\nextends = \"B\"\ndensity = \"2000 kg/m^3\"\n").unwrap();
    let lib = Library::load_dir(&d).unwrap();
    std::fs::remove_dir_all(&d).ok();
    let c = lib.get("c").unwrap();
    assert_eq!(c.density, Some(2000.0));
    assert_eq!(c.electrical.permittivity.as_ref().unwrap().as_scalar(), Some(2.0));
    assert_eq!(c.name, "c");
    assert_eq!(lib.get("b").unwrap().density, Some(1000.0));
}

#[test]
fn a_table_is_replaced_as_a_whole() {
    let lib = library();
    // silicon has a cubic stiffness; replacing [elastic] by isotropic constants leaves nothing of the class behind
    let m = lib.read("format = 1\nextends = \"silicon\"\n[elastic]\nunit = \"GPa\"\nE = 70\nnu = 0.2\n", "soft").unwrap();
    assert_eq!(m.elastic.as_ref().unwrap().detect_class(), Class::Isotropic);
    assert_eq!(m.thermal, lib.get("silicon").unwrap().thermal);
    // a unit of the base does not follow a replaced value
    let m = lib.read("format = 1\nextends = \"silicon\"\ndensity = \"2.3 g/cm^3\"\n", "x").unwrap();
    assert_eq!(m.density, Some(2300.0));
    let e = lib.read("format = 1\nextends = \"silicon\"\n[thermal]\nconductivity = 150\n", "x").unwrap_err().to_string();
    assert!(e.contains("no unit"), "{e}");
}

#[test]
fn the_combined_result_is_validated() {
    let lib = library();
    let e = lib.read("format = 1\nextends = \"silicon\"\n[optical]\nwavelength_unit = \"um\"\nrange = [1.0, 2.0]\nmodel = \"tabulated\"\ndata = [[1.2, 3.5], [1.8, 3.4]]\n", "x").unwrap_err().to_string();
    assert!(e.contains("wider than the table"), "{e}");
    let e = lib.read("format = 1\nextends = \"silicon\"\n[elastic]\nclass = \"cubic\"\nunit = \"GPa\"\nc11 = 10\nc12 = 60\nc44 = 5\n", "x").unwrap_err().to_string();
    assert!(e.contains("not positive definite"), "{e}");
    let e = lib.read("format = 1\nextends = \"silicon\"\n[electrical]\npermittivity = -3\n", "x").unwrap_err().to_string();
    assert!(e.contains("relative permittivity must be positive"), "{e}");
}

#[test]
fn a_missing_base_is_named_with_a_suggestion() {
    let lib = library();
    let e = lib.read("format = 1\nextends = \"silicn\"\n", "mine").unwrap_err().to_string();
    assert!(e.contains("material 'mine'") && e.contains("extends = \"silicn\"") && e.contains("\"silicon\"") && e.contains("available:"), "{e}");
    let d = temp_dir("missing");
    std::fs::write(d.join("x.toml"), "format = 1\nextends = \"nothing\"\n").unwrap();
    let errs = Library::load_dir(&d).unwrap_err();
    std::fs::remove_dir_all(&d).ok();
    assert!(errs[0].to_string().contains("no material of that name"), "{}", errs[0]);
}

#[test]
fn a_loop_is_reported_as_a_chain() {
    let d = temp_dir("loop");
    std::fs::write(d.join("a.toml"), "format = 1\nextends = \"b\"\n").unwrap();
    std::fs::write(d.join("b.toml"), "format = 1\nextends = \"c\"\n").unwrap();
    std::fs::write(d.join("c.toml"), "format = 1\nextends = \"a\"\n").unwrap();
    std::fs::write(d.join("self.toml"), "format = 1\nextends = \"self\"\n").unwrap();
    let errs = Library::load_dir(&d).unwrap_err();
    std::fs::remove_dir_all(&d).ok();
    let all: Vec<String> = errs.iter().map(|e| e.to_string()).collect();
    assert!(all.iter().any(|e| e.contains("loop: a -> b -> c -> a")), "{all:?}");
    assert!(all.iter().any(|e| e.contains("loop: self -> self")), "{all:?}");
}

#[test]
fn extends_without_a_library_is_explained() {
    let e = MaterialSpec::from_toml("format = 1\nextends = \"silicon\"\n").unwrap_err().to_string();
    assert!(e.contains("together with its base material") && e.contains("Library::load_dir"), "{e}");
}

#[test]
fn format_2_deep_merges_subtables_and_resets_conflicting_models() {
    let d = temp_dir("v2deep");
    std::fs::write(
        d.join("base.toml"),
        "format = 2\nkind = \"crystalline\"\npoint_group = \"m-3m\"\ndensity = \"2329 kg/m^3\"\n[dielectric]\neps_clamped = 11.45\nloss_tangent = 1e-6\ntls_tan_delta_0 = 2e-6\n[optical]\nwavelength_unit = \"um\"\nrange = [1.2, 6.0]\nmodel = \"sellmeier\"\nB = [10.6684293]\nC = [0.09091319]\n",
    )
    .unwrap();
    std::fs::write(
        d.join("child.toml"),
        "format = 2\nextends = \"base\"\n[dielectric]\nloss_tangent = 5e-8\n[optical]\nmodel = \"constant\"\nn = 3.45\n",
    )
    .unwrap();
    let lib = Library::load_dir(&d).unwrap();
    std::fs::remove_dir_all(&d).ok();

    let child = lib.get("child").unwrap();
    let dl = child.dielectric.as_ref().unwrap();
    assert_eq!(dl.eps_clamped.as_ref().unwrap().as_scalar(), Some(11.45));
    assert_eq!(dl.loss_tangent.as_ref().unwrap().as_scalar(), Some(5e-8));
    assert_eq!(dl.tls_tan_delta_0, Some(2e-6));
    assert!((child.optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0 - 3.45).abs() < 1e-12);
}

#[test]
fn centrosymmetric_point_group_forbids_piezo_and_pockels() {
    let e_pz = MaterialSpec::from_toml(
        "format = 2\npoint_group = \"m-3m\"\n[piezoelectric]\nform = \"stress_charge\"\nunit = \"C/m^2\"\ne14 = 0.5\n",
    )
    .unwrap_err()
    .to_string();
    assert!(e_pz.contains("forbids linear piezoelectricity"), "{e_pz}");

    let e_eo = MaterialSpec::from_toml(
        "format = 2\npoint_group = \"-3m\"\n[electro_optic]\nunit = \"pm/V\"\nr33 = 10.0\n",
    )
    .unwrap_err()
    .to_string();
    assert!(e_eo.contains("forbids the linear electro-optic"), "{e_eo}");
}
