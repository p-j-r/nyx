// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The materials of a session: the library of material files (`materials/*.toml`, ADR 29 to 35) and the materials a project carries as TOML
//! text, which may extend a library material. A problem names a material of either (`MaterialAssign::Named` or `Oriented`); the engine
//! hands the merged library to the physics through `Ctx` and works out the wavelength an optical model is evaluated at from the study.

use crate::engine::Engine;
use femcore::material::{Library, MaterialSpec};
use nyx_physics::{Error, Physics, ProblemSpec};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// The folders searched for bundled/configured material files, first valid match wins when loading the base library:
/// `NYX_MATERIALS` (a path list), split/bundle layouts near `current_exe()` (and its canonical target for symlinks),
/// and an ancestor `materials` folder of the program or current working directory (a build tree).
pub fn default_dirs() -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = vec![];
    let mut push_unique = |p: PathBuf| {
        if !v.contains(&p) {
            v.push(p);
        }
    };
    if let Some(p) = std::env::var_os("NYX_MATERIALS") {
        for dir in std::env::split_paths(&p) {
            push_unique(dir);
        }
    }
    let has_files = |d: &Path| {
        std::fs::read_dir(d)
            .map(|r| r.filter_map(|e| e.ok()).any(|e| e.path().extension().and_then(|x| x.to_str()) == Some("toml")))
            .unwrap_or(false)
    };
    let mut exe_dirs: Vec<PathBuf> = vec![];
    if let Ok(exe) = std::env::current_exe() {
        for cand in [Some(exe.clone()), exe.canonicalize().ok()] {
            if let Some(d) = cand.as_deref().and_then(Path::parent) {
                let dp = d.to_path_buf();
                if !exe_dirs.contains(&dp) {
                    exe_dirs.push(dp);
                }
            }
        }
    }
    for d in &exe_dirs {
        for rel in [
            d.join("materials"),
            d.join("..").join("materials"),
            d.join("..").join("Resources").join("materials"),
            d.join("..").join("share").join("nyx").join("materials"),
        ] {
            if has_files(&rel) {
                push_unique(rel);
            }
        }
    }
    let mut starts = exe_dirs;
    if let Ok(cwd) = std::env::current_dir() {
        starts.push(cwd);
    }
    for s in starts {
        let mut d: Option<&Path> = Some(&s);
        while let Some(dir) = d {
            let m = dir.join("materials");
            if has_files(&m) {
                push_unique(m);
                break;
            }
            d = dir.parent();
        }
    }
    v
}

/// The user's writable material overlay folder (`<workspace>/materials`, i.e. `$NYX_WORKSPACE/materials` or `Documents/nyx/materials`).
/// Installed builds write fetched or custom material `.toml` files here instead of mutating read-only or signed install directories.
pub fn user_materials_dir() -> Option<PathBuf> {
    Some(crate::shared_queue::default_workspace().join("materials"))
}

fn read_toml_map(dir: &Path, notes: &mut Vec<String>) -> std::collections::BTreeMap<String, String> {
    let mut files: Vec<PathBuf> = match std::fs::read_dir(dir) {
        Ok(rd) => rd
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| p.is_file() && p.extension().and_then(|x| x.to_str()) == Some("toml"))
            .collect(),
        Err(_) => return Default::default(),
    };
    files.sort();
    let mut texts = std::collections::BTreeMap::new();
    for path in files {
        let Some(stem) = path.file_stem().and_then(|s| s.to_str()) else { continue };
        match std::fs::read_to_string(&path) {
            Ok(txt) => {
                texts.insert(stem.to_string(), txt);
            }
            Err(e) => notes.push(format!("{}: {e}", path.display())),
        }
    }
    texts
}

/// Loads the first folder that has material files as the base library, then overlays any user `*.toml` files
/// from [`user_materials_dir`] via `Library::with_texts` so user materials can `extends` bundled ones without hiding them.
pub fn load_library(dirs: &[PathBuf]) -> (Library, Vec<String>) {
    let mut notes = vec![];
    let mut loaded_dir: Option<&Path> = None;
    let mut base = None;
    for d in dirs {
        match Library::load_dir(d) {
            Ok(l) => {
                base = Some(l);
                loaded_dir = Some(d.as_path());
                break;
            }
            Err(errs) => notes.extend(errs.iter().map(|e| format!("{}: {e}", d.display()))),
        }
    }
    let mut lib = match base {
        Some(l) => l,
        None => {
            if dirs.is_empty() {
                notes.push(
                    "no material folder found (set NYX_MATERIALS, or keep a `materials` folder next to the program)"
                        .into(),
                );
            }
            Library::default()
        }
    };
    if let Some(user_dir) = user_materials_dir().filter(|p| p.is_dir()) {
        let same = loaded_dir
            .and_then(|b| b.canonicalize().ok())
            .zip(user_dir.canonicalize().ok())
            .is_some_and(|(a, b)| a == b);
        if !same {
            let texts = read_toml_map(&user_dir, &mut notes);
            if !texts.is_empty() {
                match lib.with_texts(&texts) {
                    Ok(merged) => lib = merged,
                    Err(_) => {
                        for (name, txt) in texts {
                            let mut one = std::collections::BTreeMap::new();
                            one.insert(name, txt);
                            match lib.with_texts(&one) {
                                Ok(merged) => lib = merged,
                                Err(errs) => notes.extend(errs.iter().map(|e| format!("{}: {e}", user_dir.display()))),
                            }
                        }
                    }
                }
            }
        }
    }
    (lib, notes)
}

/// The frequency (Hz) of a study, from the first parameter of the study descriptor with the unit Hz (given, else its default),
/// or an explicit `"freq"` / `"frequency"` study parameter.
pub fn study_frequency(p: &dyn Physics, spec: &ProblemSpec) -> Option<f64> {
    use nyx_physics::ParamKind;
    let studies = p.studies();
    let sd = studies.iter().find(|s| s.kind == spec.study.kind)?;
    if spec.study.kind == "driven" {
        if let Some(arr) = spec.study.params.get("freqs_hz").and_then(Value::as_array) {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            let mut count = 0usize;
            for v in arr {
                if let Some(f) = v.as_f64().filter(|&x| x > 0.0 && x.is_finite()) {
                    lo = lo.min(f);
                    hi = hi.max(f);
                    count += 1;
                }
            }
            if count > 0 {
                return Some(0.5 * (lo + hi));
            }
        }
        let param_num = |name: &str| -> Option<f64> {
            spec.study
                .params
                .get(name)
                .and_then(Value::as_f64)
                .or_else(|| {
                    sd.params.iter().find(|d| d.name == name).and_then(|d| match d.kind {
                        ParamKind::Number { default, .. } => Some(default),
                        _ => None,
                    })
                })
                .filter(|&f| f > 0.0 && f.is_finite())
        };
        if let (Some(f_min), Some(f_max)) = (param_num("f_min"), param_num("f_max")) {
            return Some(0.5 * (f_min + f_max));
        }
    }
    if let Some(d) = sd.params.iter().find(|d| d.unit == "Hz") {
        let f = spec.study.params.get(&d.name).and_then(Value::as_f64).or(match d.kind {
            ParamKind::Number { default, .. } => Some(default),
            _ => None,
        })?;
        return (f >= 0.0 && f.is_finite()).then_some(f);
    }
    spec.study
        .params
        .get("freq")
        .or_else(|| spec.study.params.get("frequency"))
        .and_then(Value::as_f64)
        .filter(|&f| f >= 0.0 && f.is_finite())
}

/// The frequency (Hz) of a study, from the first parameter of the study descriptor with the unit Hz (given, else its default); the free-space
/// wavelength follows. `None` for a study without a frequency, or for quasi-static `dfs` `electrostatic` studies where `freq` evaluates
/// `[dielectric]` / `[electrical]` AC loss rather than `[optical]` dispersion.
pub fn study_wavelength(p: &dyn Physics, spec: &ProblemSpec) -> Option<f64> {
    if (p.id() == "dfs" && spec.study.kind == "electrostatic") || p.id() == "sms" {
        return None;
    }
    let f = study_frequency(p, spec)?;
    (f > 0.0 && f.is_finite()).then(|| 299_792_458.0 / f)
}

pub(crate) fn describe(m: &MaterialSpec) -> Value {
    json!({
        "format": m.format,
        "name": m.name,
        "description": m.description,
        "notes": m.notes,
        "kind": m.kind.map(|k| k.as_str()),
        "point_group": m.point_group.map(|p| p.as_str()),
        "crystal_setting": m.crystal_setting.map(|c| c.as_str()),
        "reference_temperature": m.reference_temperature,
        "density": m.density,
        "residual_stress": m.residual_stress,
        "elastic": m.elastic.as_ref().map(|e| json!({ "class": format!("{:?}", e.detect_class()) })),
        "optical": m.optical.is_some(),
        "thermal": m.thermal.conductivity.is_some()
            || m.thermal.specific_heat.is_some()
            || m.thermal.expansion.is_some()
            || m.thermal.conductivity_t.is_some()
            || m.thermal.specific_heat_t.is_some()
            || m.thermal.expansion_t.is_some(),
        "electrical": m.electrical.permittivity.is_some()
            || m.electrical.permeability.is_some()
            || m.electrical.conductivity.is_some(),
        "dielectric": m.dielectric.is_some(),
        "superconducting": m.superconducting.is_some(),
        "piezoelectric": m.piezoelectric.is_some(),
        "photoelastic": m.photoelastic.is_some(),
        "electro_optic": m.electro_optic.is_some(),
        "viscoelastic": m.viscoelastic.is_some(),
        "interfaces": !m.interfaces.is_empty(),
        "toml": m.to_toml(),
    })
}

impl Engine {
    /// The library a solve sees: the folder's materials and the project's own, which may extend them.
    pub fn library(&self) -> Result<Arc<Library>, Error> {
        if self.project.materials.is_empty() {
            return Ok(self.base_library().clone());
        }
        let lib = self.base_library().with_texts(&self.project.materials).map_err(|errs| {
            Error::invalid(format!(
                "the materials of the project do not read: {}",
                errs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join("; ")
            ))
        })?;
        Ok(Arc::new(lib))
    }

    /// Replaces the folder library (tests, or a front-end that lets the person pick another folder).
    pub fn set_base_library(&mut self, library: Library) {
        self.set_base_library_arc(Arc::new(library));
    }

    pub fn list_materials(&self) -> Result<Value, Error> {
        let lib = self.library()?;
        Ok(json!(lib
            .names()
            .iter()
            .filter_map(|n| lib.get(n))
            .map(|m| json!({
                "format": m.format,
                "name": m.name,
                "description": m.description,
                "kind": m.kind.map(|k| k.as_str()),
                "point_group": m.point_group.map(|p| p.as_str()),
                "in_project": self.project.materials.keys().any(|k| k.eq_ignore_ascii_case(&m.name)),
                "optical": m.optical.is_some(),
                "elastic": m.elastic.is_some(),
                "dielectric": m.dielectric.is_some(),
                "superconducting": m.superconducting.is_some(),
                "piezoelectric": m.piezoelectric.is_some(),
                "photoelastic": m.photoelastic.is_some(),
                "electro_optic": m.electro_optic.is_some(),
                "viscoelastic": m.viscoelastic.is_some(),
            }))
            .collect::<Vec<_>>()))
    }

    pub fn describe_material(&self, name: &str) -> Result<Value, Error> {
        let lib = self.library()?;
        let m = lib.require(name).map_err(|e| Error::invalid(e.to_string()))?;
        Ok(describe(m))
    }

    /// Evaluates a material's temperature-, frequency-, wavelength-, and field-dependent multi-physics properties.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_material(
        &self,
        name: &str,
        temperature_k: Option<f64>,
        frequency_hz: Option<f64>,
        wavelength_m: Option<f64>,
        e_field_v_per_m: Option<f64>,
        thickness_m: Option<f64>,
        film_model: Option<&str>,
        euler_deg: Option<[f64; 3]>,
    ) -> Result<Value, Error> {
        use femcore::material::{bcs_surface_impedance, FilmModel};
        let lib = self.library()?;
        let m = lib.require(name).map_err(|e| Error::invalid(e.to_string()))?;
        let t_k = temperature_k.or(m.reference_temperature);
        let euler = euler_deg.unwrap_or([0.0; 3]);

        let mut out = serde_json::Map::new();
        out.insert("name".into(), json!(m.name));
        out.insert("format".into(), json!(m.format));
        out.insert("kind".into(), json!(m.kind.map(|k| k.as_str())));
        out.insert("point_group".into(), json!(m.point_group.map(|p| p.as_str())));
        out.insert("temperature_k".into(), json!(t_k));
        out.insert("euler_deg".into(), json!(euler));
        if let Some(rho) = m.density {
            out.insert("density_kg_m3".into(), json!(rho));
        }

        // Thermal evaluation at T
        let k_t = m.thermal_conductivity_at(t_k).map(|tp| tp.matrix());
        let cp_t = m.specific_heat_at(t_k);
        let alpha_t = m.thermal_expansion_at(t_k).map(|tp| tp.matrix());
        if k_t.is_some() || cp_t.is_some() || alpha_t.is_some() {
            out.insert(
                "thermal".into(),
                json!({
                    "conductivity_w_mk": k_t,
                    "specific_heat_j_kgk": cp_t,
                    "expansion_1_k": alpha_t,
                }),
            );
        }

        // Optical evaluation at (wavelength, T)
        let opt_lambda =
            wavelength_m.or_else(|| frequency_hz.filter(|&f| f >= 1.0e12 && f.is_finite()).map(|f| 299_792_458.0 / f));
        if let (Some(opt), Some(lam)) = (&m.optical, opt_lambda) {
            if let Ok(eps_c) = m.optical_permittivity(lam, t_k, euler) {
                let nk0 = opt.nk(0, lam, t_k).ok();
                let re = [
                    [eps_c[0][0].re, eps_c[0][1].re, eps_c[0][2].re],
                    [eps_c[1][0].re, eps_c[1][1].re, eps_c[1][2].re],
                    [eps_c[2][0].re, eps_c[2][1].re, eps_c[2][2].re],
                ];
                let im = [
                    [eps_c[0][0].im, eps_c[0][1].im, eps_c[0][2].im],
                    [eps_c[1][0].im, eps_c[1][1].im, eps_c[1][2].im],
                    [eps_c[2][0].im, eps_c[2][1].im, eps_c[2][2].im],
                ];
                out.insert(
                    "optical".into(),
                    json!({
                        "wavelength_m": lam,
                        "n": nk0.map(|v| v.0),
                        "k": nk0.map(|v| v.1),
                        "eps_re": re,
                        "eps_im": im,
                    }),
                );
            }
        }

        // Microwave / RF complex permittivity at (f, T, E)
        if let Some((re, im_loss)) = m.complex_microwave_permittivity(frequency_hz, t_k, e_field_v_per_m, false, euler)
        {
            let tan_delta_xx = if re[0][0].abs() > 0.0 { im_loss[0][0] / re[0][0] } else { 0.0 };
            let tan_delta_zz = if re[2][2].abs() > 0.0 { im_loss[2][2] / re[2][2] } else { 0.0 };
            out.insert(
                "dielectric".into(),
                json!({
                    "frequency_hz": frequency_hz,
                    "e_field_v_per_m": e_field_v_per_m,
                    "eps_re": re,
                    "eps_im_loss": im_loss,
                    "tan_delta_xx": tan_delta_xx,
                    "tan_delta_zz": tan_delta_zz,
                }),
            );
        }

        // Piezoelectric Legendre duals (eps^S, eps^T, C^E, C^D) rotated into the wafer frame
        if let Some(pz) = &m.piezoelectric {
            use femcore::material::{rotate_m6_textbook, rotate_tensor, rotation, PiezoForm, Piezoelectric};
            let ce_voigt = m.elastic.as_ref().map(|e| e.c);
            let rot = rotation(euler);
            let rotated = euler != [0.0; 3];
            let rot3 = |t: femcore::material::M3| if rotated { rotate_tensor(&t, &rot) } else { t };
            let rot6 = |c: femcore::material::M6| if rotated { rotate_m6_textbook(&c, euler) } else { c };
            let e_rot = pz.stress_charge(ce_voigt.as_ref()).ok().map(|e| {
                if rotated {
                    Piezoelectric { source: Vec::new(), form: PiezoForm::StressCharge, matrix: e }.rotate(&rot)
                } else {
                    e
                }
            });
            let d_rot = pz.strain_charge(ce_voigt.as_ref()).ok().map(|d| {
                if rotated {
                    Piezoelectric { source: Vec::new(), form: PiezoForm::StrainCharge, matrix: d }.rotate(&rot)
                } else {
                    d
                }
            });
            out.insert(
                "piezoelectric".into(),
                json!({
                    "eps_clamped": m.permittivity_clamped().map(rot3),
                    "eps_free": m.permittivity_free().map(rot3),
                    "c_e_pa": m.stiffness_constant_e().map(rot6),
                    "c_d_pa": m.stiffness_constant_d().map(rot6),
                    "e_c_m2": e_rot,
                    "d_c_n": d_rot,
                }),
            );
        }

        if let Some(pe) = &m.photoelastic {
            let p_mat =
                if euler != [0.0; 3] { pe.rotate_voigt(&femcore::material::rotation(euler)) } else { pe.matrix };
            out.insert("photoelastic".into(), json!({ "wavelength_m": pe.wavelength, "p_voigt": p_mat }));
        }
        if let Some(eo) = &m.electro_optic {
            let ce_voigt = m.elastic.as_ref().map(|e| e.c);
            let d_mat = m.piezoelectric.as_ref().and_then(|p| p.strain_charge(ce_voigt.as_ref()).ok());
            out.insert(
                "electro_optic".into(),
                json!({
                    "wavelength_m": eo.wavelength,
                    "r_clamped_m_per_v": eo.resolved_clamped(m.photoelastic.as_ref(), d_mat.as_ref()),
                    "r_free_m_per_v": eo.resolved_free(m.photoelastic.as_ref(), d_mat.as_ref()),
                }),
            );
        }
        if let Some(ve) = &m.viscoelastic {
            let eta_mat =
                if euler != [0.0; 3] { femcore::material::rotate_m6_textbook(&ve.matrix, euler) } else { ve.matrix };
            out.insert("viscoelastic".into(), json!({ "eta_pa_s": eta_mat, "cryo_qm_floor": ve.cryo_qm_floor }));
        }

        // Superconducting BCS evaluation if [superconducting] is present
        if let Some(sc) = &m.superconducting {
            let f_eval = frequency_hz.unwrap_or(5.0e9);
            let t_eval = temperature_k.unwrap_or(0.02);
            let d_eval = thickness_m.or(sc.nominal_thickness).unwrap_or(100.0e-9);
            let fm = match film_model {
                Some(s) => FilmModel::parse(s).ok_or_else(|| {
                    Error::invalid(format!("unknown film_model '{s}' (expected 'kerr', 'two_sided' or 'thin_limit')"))
                })?,
                None => sc.film_model.unwrap_or(FilmModel::Kerr),
            };
            if f_eval > 0.0 && t_eval >= 0.0 && d_eval > 0.0 {
                if let Ok(ev) = bcs_surface_impedance(sc, f_eval, t_eval, d_eval, fm) {
                    out.insert(
                        "superconducting".into(),
                        json!({
                            "tc_k": sc.tc,
                            "frequency_hz": f_eval,
                            "temperature_k": t_eval,
                            "thickness_m": d_eval,
                            "film_model": fm.as_str(),
                            "r_s_ohm_sq": ev.r_s,
                            "x_s_ohm_sq": ev.z_s.im,
                            "l_k_sq_h_sq": ev.l_k_sq,
                            "lambda_eff_m": ev.lambda_eff,
                            "delta_t_j": ev.delta_j,
                            "delta_t_mev": ev.delta_j / (1.602_176_634e-19 * 1.0e-3),
                            "sigma_1_s_per_m": ev.sigma.re,
                            "sigma_2_s_per_m": -ev.sigma.im,
                        }),
                    );
                }
            }
        }

        Ok(Value::Object(out))
    }

    /// Evaluates Zimmermann arbitrary-purity BCS + sub-gap mK surface impedance for a named superconducting
    /// material or explicit parameters (`tc`, `gap_0`, `gap_0_eV`, `lambda_L`, `xi_0`, `mean_free_path`/`ell`, `sigma_n`, `dynes_gamma`, `r_res`, `x_qp_floor`).
    pub fn bcs_surface_impedance_cmd(&self, args: &Value) -> Result<Value, Error> {
        use femcore::material::{bcs_surface_impedance, FilmModel, Superconducting, EV_TO_J};
        let freq = args
            .get("freq")
            .or_else(|| args.get("frequency"))
            .and_then(Value::as_f64)
            .ok_or_else(|| Error::invalid("`freq` (Hz) is required and must be a number"))?;
        let temp = args.get("temperature").and_then(Value::as_f64).unwrap_or(0.01);

        let mat_name = args.get("material").and_then(Value::as_str).map(str::trim).filter(|s| !s.is_empty());
        let mut sc = if let Some(name) = mat_name {
            let lib = self.library()?;
            let m = lib.require(name).map_err(|e| Error::invalid(e.to_string()))?;
            m.superconducting
                .clone()
                .ok_or_else(|| Error::invalid(format!("material '{name}' has no [superconducting] block")))?
        } else {
            let tc = args
                .get("tc")
                .and_then(Value::as_f64)
                .ok_or_else(|| Error::invalid("provide either `material` or `tc` (K) for bcs_surface_impedance"))?;
            Superconducting { tc, ..Default::default() }
        };

        if let Some(v) = args.get("tc").and_then(Value::as_f64) {
            sc.tc = v;
        }
        if let Some(v) = args
            .get("gap_0")
            .and_then(Value::as_f64)
            .or_else(|| args.get("gap_0_eV").and_then(Value::as_f64).map(|e| e * EV_TO_J))
        {
            sc.gap_0 = Some(v);
        }
        if let Some(v) = args.get("lambda_L").and_then(Value::as_f64) {
            sc.lambda_l = Some(v);
        }
        if let Some(v) = args.get("xi_0").and_then(Value::as_f64) {
            sc.xi_0 = Some(v);
        }
        if let Some(v) = args.get("mean_free_path").or_else(|| args.get("ell")).and_then(Value::as_f64) {
            sc.mean_free_path = Some(v);
        }
        if let Some(v) = args.get("sigma_n").and_then(Value::as_f64) {
            sc.sigma_n = Some(v);
        }
        if let Some(v) = args.get("dynes_gamma").and_then(Value::as_f64) {
            sc.dynes_gamma = Some(v);
        }
        if let Some(v) = args.get("r_res").and_then(Value::as_f64) {
            sc.r_res = Some(v);
        }
        if let Some(v) = args.get("x_qp_floor").and_then(Value::as_f64) {
            sc.x_qp_floor = Some(v);
        }

        let thickness = args.get("thickness").and_then(Value::as_f64).or(sc.nominal_thickness).unwrap_or(100.0e-9);
        let fm = match args.get("film_model").and_then(Value::as_str) {
            Some(s) => FilmModel::parse(s).ok_or_else(|| {
                Error::invalid(format!("unknown film_model '{s}'; use two_sided, kerr or thin_limit"))
            })?,
            None => sc.film_model.unwrap_or(FilmModel::Kerr),
        };

        let ev = bcs_surface_impedance(&sc, freq, temp, thickness, fm).map_err(Error::invalid)?;
        Ok(json!({
            "material": mat_name,
            "tc_k": sc.tc,
            "frequency_hz": freq,
            "temperature_k": temp,
            "thickness_m": thickness,
            "film_model": fm.as_str(),
            "R_s": ev.r_s,
            "L_k_sq": ev.l_k_sq,
            "L_eff": ev.l_eff,
            "G_eff": ev.g_eff,
            "sigma_re": ev.sigma.re,
            "sigma_im": ev.sigma.im,
            "sigma_n": ev.sigma_n,
            "delta_J": ev.delta_j,
            "delta_meV": ev.delta_j / (EV_TO_J * 1e-3),
            "lambda_eff": ev.lambda_eff,
            "r_s_ohm_sq": ev.r_s,
            "x_s_ohm_sq": ev.z_s.im,
            "l_k_sq_h_sq": ev.l_k_sq,
            "lambda_eff_m": ev.lambda_eff,
            "delta_t_j": ev.delta_j,
            "delta_t_mev": ev.delta_j / (EV_TO_J * 1e-3),
            "sigma_1_s_per_m": ev.sigma.re,
            "sigma_2_s_per_m": -ev.sigma.im,
        }))
    }

    /// Adds (or replaces) a material the project carries, given as TOML text. It is read now, so a mistake is reported here.
    pub fn add_material(&mut self, name: &str, toml_text: &str) -> Result<Value, Error> {
        if name.trim().is_empty() || name.contains(['/', '\\', '*', ',']) {
            return Err(Error::invalid("a material name is plain text without slashes, commas or stars"));
        }
        let mut texts = self.project.materials.clone();
        texts.insert(name.to_string(), toml_text.to_string());
        let lib = self
            .base_library()
            .with_texts(&texts)
            .map_err(|errs| Error::invalid(errs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join("; ")))?;
        let m = lib.get(name).expect("just added");
        let v = describe(m);
        self.project.materials = texts;
        Ok(v)
    }

    pub fn remove_material(&mut self, name: &str) -> Result<(), Error> {
        if self.project.materials.remove(name).is_none() {
            return Err(Error::invalid(format!("the project carries no material '{name}'")));
        }
        Ok(())
    }

    /// Evaluates Zimmermann arbitrary-purity BCS complex conductivity and thin-film surface impedance.
    pub fn bcs_surface_impedance(&self, args: &Value) -> Result<Value, Error> {
        self.bcs_surface_impedance_cmd(args)
    }
}
