//! Inspector tabs (geometry, materials, physics, study, results) and the results panel.

use crate::app::App;
use crate::colormap::Colormap;
use crate::state::*;
use eframe::egui::{self, Color32, DragValue, RichText, Ui};
use femcore::material::{self, Anisotropy, Stiffness, Tensor};
use femcore::problem::*;

// ------------------------------------------------------------------ helpers

pub fn fmt_num(n: f64) -> String {
    if n == 0.0 {
        "0".into()
    } else if n.abs() >= 1e5 || n.abs() < 1e-3 {
        format!("{n:.5e}")
    } else {
        let s = format!("{n:.6}");
        let s = s.trim_end_matches('0').trim_end_matches('.');
        s.to_string()
    }
}

/// A number field: drag to change, or double-click and type (scientific notation accepted).
pub fn num(ui: &mut Ui, v: &mut f64) -> egui::Response {
    let speed = if *v == 0.0 { 0.01 } else { v.abs() * 0.005 };
    ui.add(DragValue::new(v).speed(speed).custom_formatter(|n, _| fmt_num(n)).custom_parser(|s| s.trim().replace(',', ".").parse::<f64>().ok()))
}

fn labels_text(l: &[i32]) -> String {
    l.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
}

fn parse_labels(s: &str) -> Option<Vec<i32>> {
    let mut out = vec![];
    for part in s.split(|c: char| c == ',' || c == ';' || c.is_whitespace()).filter(|p| !p.is_empty()) {
        if let Some((a, b)) = part.split_once('-') {
            let (a, b) = (a.parse::<i32>().ok()?, b.parse::<i32>().ok()?);
            out.extend(a..=b);
        } else {
            out.push(part.parse().ok()?);
        }
    }
    out.sort();
    out.dedup();
    Some(out)
}

struct Material {
    name: &'static str,
    props: &'static [(&'static str, f64)],
}

const MATERIALS: &[Material] = &[
    Material { name: "Vacuum / air", props: &[("eps_r", 1.0), ("mu_r", 1.0), ("k", 0.026), ("rho", 1.2), ("Cp", 1005.0)] },
    Material { name: "Silicon", props: &[("eps_r", 11.7), ("mu_r", 1.0), ("E", 170e9), ("nu", 0.28), ("rho", 2329.0), ("k", 130.0), ("Cp", 700.0), ("alpha", 2.6e-6)] },
    Material { name: "Silicon dioxide (fused silica)", props: &[("eps_r", 3.75), ("mu_r", 1.0), ("E", 73.1e9), ("nu", 0.17), ("rho", 2203.0), ("k", 1.38), ("Cp", 703.0), ("alpha", 0.55e-6)] },
    Material { name: "Silicon nitride", props: &[("eps_r", 7.5), ("mu_r", 1.0), ("E", 250e9), ("nu", 0.23), ("rho", 3100.0), ("k", 20.0), ("Cp", 700.0), ("alpha", 2.3e-6)] },
    Material { name: "Sapphire (isotropic approximation)", props: &[("eps_r", 9.4), ("mu_r", 1.0), ("E", 345e9), ("nu", 0.28), ("rho", 3980.0), ("k", 35.0), ("Cp", 760.0), ("alpha", 5.3e-6)] },
    Material { name: "Aluminium", props: &[("E", 70e9), ("nu", 0.33), ("rho", 2700.0), ("k", 238.0), ("Cp", 900.0), ("alpha", 23e-6)] },
    Material { name: "Copper", props: &[("E", 110e9), ("nu", 0.35), ("rho", 8960.0), ("k", 400.0), ("Cp", 385.0), ("alpha", 17e-6)] },
    Material { name: "Structural steel", props: &[("E", 200e9), ("nu", 0.3), ("rho", 7850.0), ("k", 44.5), ("Cp", 475.0), ("alpha", 12.3e-6)] },
    Material { name: "Water", props: &[("eps_r", 80.0), ("mu_r", 1.0), ("rho", 997.0), ("k", 0.6), ("Cp", 4181.0)] },
    Material { name: "FR-4", props: &[("eps_r", 4.4), ("mu_r", 1.0), ("E", 22e9), ("nu", 0.15), ("rho", 1850.0), ("k", 0.3), ("Cp", 1100.0), ("alpha", 14e-6)] },
];

/// (key, label, unit, hint)
fn props_for(p: Physics) -> &'static [(&'static str, &'static str, &'static str)] {
    match p {
        Physics::Electrostatics => &[("eps_r", "relative permittivity ε_r", ""), ("rho_v", "space charge density ρ_v", "C/m³")],
        Physics::Heat => &[("k", "thermal conductivity k", "W/(m·K)"), ("rho", "density ρ", "kg/m³"), ("Cp", "heat capacity C_p", "J/(kg·K)"), ("Q", "heat source Q", "W/m³")],
        Physics::Solid | Physics::Brillouin => &[("E", "Young's modulus E", "Pa"), ("nu", "Poisson's ratio ν", ""), ("rho", "density ρ", "kg/m³"), ("alpha", "thermal expansion α", "1/K")],
        Physics::Emw => &[("eps_r", "relative permittivity ε_r", ""), ("mu_r", "relative permeability μ_r", "")],
    }
}

fn default_prop(key: &str) -> f64 {
    match key {
        "eps_r" | "mu_r" => 1.0,
        "rho_v" | "Q" | "alpha" => 0.0,
        "nu" => 0.3,
        "E" => 200e9,
        "k" => 1.0,
        "rho" => 1000.0,
        "Cp" => 500.0,
        _ => 0.0,
    }
}

struct Crystal {
    name: &'static str,
    make: fn() -> (Anisotropy, f64),
}

const CRYSTALS: &[Crystal] = &[
    Crystal {
        name: "Sapphire (anisotropic, c axis = z)",
        make: || {
            let (st, eps, k, alpha, rho) = material::sapphire();
            (Anisotropy { stiffness: Some(st), conductivity: Some(k), permittivity: Some(eps), expansion: Some(alpha), euler_deg: [0.0; 3] }, rho)
        },
    },
    Crystal {
        name: "Silicon (cubic, [100] along x, y, z)",
        make: || {
            (
                Anisotropy {
                    stiffness: Some(material::silicon()),
                    conductivity: Some(Tensor::Diagonal { x: 130.0, y: 130.0, z: 130.0 }),
                    permittivity: Some(Tensor::Diagonal { x: 11.7, y: 11.7, z: 11.7 }),
                    expansion: Some(Tensor::Diagonal { x: 2.6e-6, y: 2.6e-6, z: 2.6e-6 }),
                    euler_deg: [0.0; 3],
                },
                2329.0,
            )
        },
    },
    Crystal {
        name: "GaAs (cubic)",
        make: || {
            (
                Anisotropy {
                    stiffness: Some(Stiffness::Cubic { c11: 118.8e9, c12: 53.8e9, c44: 59.4e9 }),
                    conductivity: Some(Tensor::Diagonal { x: 55.0, y: 55.0, z: 55.0 }),
                    permittivity: Some(Tensor::Diagonal { x: 12.9, y: 12.9, z: 12.9 }),
                    expansion: Some(Tensor::Diagonal { x: 5.7e-6, y: 5.7e-6, z: 5.7e-6 }),
                    euler_deg: [0.0; 3],
                },
                5316.0,
            )
        },
    },
];

fn tensor_combo(ui: &mut Ui, id: &str, label: &str, unit: &str, t: &mut Option<Tensor>, scalar: f64) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label(label);
        let cur = t.as_ref().map(|x| x.name()).unwrap_or("Isotropic");
        egui::ComboBox::from_id_salt(id).selected_text(cur).show_ui(ui, |ui| {
            if ui.selectable_label(t.is_none(), "Isotropic (use the value above)").clicked() && t.is_some() {
                *t = None;
                changed = true;
            }
            for (name, make) in [
                ("Uniaxial", Tensor::Uniaxial { perp: scalar, par: scalar }),
                ("Diagonal", Tensor::Diagonal { x: scalar, y: scalar, z: scalar }),
                ("Symmetric", Tensor::Symmetric { xx: scalar, yy: scalar, zz: scalar, xy: 0.0, yz: 0.0, xz: 0.0 }),
            ] {
                if ui.selectable_label(t.as_ref().map(|x| x.name()) == Some(name), name).clicked() && t.as_ref().map(|x| x.name()) != Some(name) {
                    *t = Some(make);
                    changed = true;
                }
            }
        });
        ui.weak(unit);
    });
    if let Some(tt) = t {
        ui.horizontal_wrapped(|ui| {
            for (n, v) in tt.constants_mut() {
                ui.label(n);
                changed |= num(ui, v).changed();
            }
        });
    }
    changed
}

/// Editor of the anisotropic data of one domain for the active physics. Returns true on change.
fn aniso_editor(ui: &mut Ui, phys: Physics, d: i32, an: &mut Anisotropy, props: &Materials) -> bool {
    let mut changed = false;
    if phys == Physics::Solid || phys == Physics::Brillouin {
        ui.horizontal(|ui| {
            ui.label("stiffness");
            let cur = an.stiffness.as_ref().map(|x| x.name()).unwrap_or("Isotropic (E, ν)");
            egui::ComboBox::from_id_salt(format!("stiff_{d}")).selected_text(cur).show_ui(ui, |ui| {
                if ui.selectable_label(an.stiffness.is_none(), "Isotropic (E, ν)").clicked() && an.stiffness.is_some() {
                    an.stiffness = None;
                    changed = true;
                }
                for (name, make) in [
                    ("Cubic", Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }),
                    ("Hexagonal", Stiffness::Hexagonal { c11: 160e9, c12: 60e9, c13: 50e9, c33: 180e9, c44: 40e9 }),
                    ("Trigonal", Stiffness::Trigonal { c11: 496.8e9, c12: 163.6e9, c13: 110.9e9, c14: -23.5e9, c33: 498.1e9, c44: 147.4e9 }),
                    ("Orthotropic", Stiffness::Orthotropic { c11: 200e9, c22: 180e9, c33: 160e9, c12: 60e9, c13: 55e9, c23: 50e9, c44: 40e9, c55: 45e9, c66: 50e9 }),
                ] {
                    if ui.selectable_label(an.stiffness.as_ref().map(|x| x.name()) == Some(name), name).clicked() && an.stiffness.as_ref().map(|x| x.name()) != Some(name) {
                        an.stiffness = Some(make);
                        changed = true;
                    }
                }
            });
            ui.weak("GPa");
        });
        if let Some(st) = &mut an.stiffness {
            ui.horizontal_wrapped(|ui| {
                for (n, v) in st.constants_mut() {
                    ui.label(n);
                    let mut g = *v / 1e9;
                    if num(ui, &mut g).changed() {
                        *v = g * 1e9;
                        changed = true;
                    }
                }
            });
        }
        let a0 = props.get(d, "alpha").unwrap_or(0.0);
        changed |= tensor_combo(ui, &format!("exp_{d}"), "thermal expansion", "1/K", &mut an.expansion, a0);
    }
    if phys == Physics::Heat {
        let k0 = props.get(d, "k").unwrap_or(1.0);
        changed |= tensor_combo(ui, &format!("cond_{d}"), "conductivity", "W/(m·K)", &mut an.conductivity, k0);
    }
    if phys == Physics::Electrostatics || phys == Physics::Emw || phys == Physics::Brillouin {
        let e0 = props.get(d, "eps_r").unwrap_or(1.0);
        changed |= tensor_combo(ui, &format!("perm_{d}"), if phys == Physics::Brillouin { "optical permittivity (replaces the refractive index)" } else { "permittivity" }, "", &mut an.permittivity, e0);
    }
    let any = an.stiffness.is_some() || an.conductivity.is_some() || an.permittivity.is_some() || an.expansion.is_some();
    if any {
        ui.horizontal(|ui| {
            ui.label("crystal orientation (Euler Z-X-Z, °)");
            for a in an.euler_deg.iter_mut() {
                changed |= ui.add(DragValue::new(a).speed(1.0).range(-360.0..=360.0)).changed();
            }
        });
    }
    changed
}

/// Photoelastic presets (constants of the crystal frame; orientation follows the elastic crystal data).
type PhotoPreset = (&'static str, fn() -> PhotoKind);
const PHOTO_PRESETS: &[PhotoPreset] = &[
    ("Silicon (-0.094, 0.017, -0.051)", || PhotoKind::Cubic { p11: -0.094, p12: 0.017, p44: -0.051 }),
    ("Sapphire (-0.23, 0.03, 0.02, 0, -0.2, -0.1)", || PhotoKind::Trigonal { p11: -0.23, p12: 0.03, p13: 0.02, p14: 0.0, p33: -0.2, p44: -0.1 }),
    ("Silica (0.121, 0.270)", || PhotoKind::Isotropic { p11: 0.121, p12: 0.270 }),
    ("Silicon nitride (0.2, 0.2) isotropic", || PhotoKind::Isotropic { p11: 0.2, p12: 0.2 }),
];

fn photo_kind_name(k: &PhotoKind) -> &'static str {
    match k {
        PhotoKind::Isotropic { .. } => "isotropic",
        PhotoKind::Cubic { .. } => "cubic",
        PhotoKind::Trigonal { .. } => "trigonal",
    }
}

fn photo_constants(k: &mut PhotoKind) -> Vec<(&'static str, &mut f64)> {
    match k {
        PhotoKind::Isotropic { p11, p12 } => vec![("p11", p11), ("p12", p12)],
        PhotoKind::Cubic { p11, p12, p44 } => vec![("p11", p11), ("p12", p12), ("p44", p44)],
        PhotoKind::Trigonal { p11, p12, p13, p14, p33, p44 } => vec![("p11", p11), ("p12", p12), ("p13", p13), ("p14", p14), ("p33", p33), ("p44", p44)],
    }
}

/// The optical side of a Brillouin study: wavelength, modes, backward or forward scattering, and per
/// domain the refractive index, the photoelastic data, the moving boundary and the optical PML.
fn brillouin_editor(ui: &mut Ui, spec: &mut ProblemSpec, domains: &[i32], bufs: &mut std::collections::HashMap<(Physics, usize), String>) {
    let Some(b) = spec.brillouin.as_mut() else {
        spec.brillouin = Some(BrillouinSpec::default());
        return;
    };
    ui.label(RichText::new("Optical modes").strong());
    egui::Grid::new("brill_opt").num_columns(2).show(ui, |ui| {
        let mut nm = b.wavelength * 1e9;
        ui.label("pump wavelength (nm)");
        if num(ui, &mut nm).changed() {
            b.wavelength = nm * 1e-9;
        }
        ui.end_row();
        ui.label("search near effective index");
        num(ui, &mut b.n_target);
        ui.end_row();
        ui.label("number of optical modes");
        ui.add(DragValue::new(&mut b.n_optical_modes).range(1..=20));
        ui.end_row();
    });
    ui.horizontal(|ui| {
        ui.label("scattering");
        ui.selectable_value(&mut b.backward, true, "backward (q = β_pump + β_Stokes)");
        ui.selectable_value(&mut b.backward, false, "forward (q = β_pump − β_Stokes)");
    });
    let mut ring = b.ring_radius.is_some();
    ui.horizontal(|ui| {
        if ui.checkbox(&mut ring, "ring resonator (axisymmetric)").changed() {
            b.ring_radius = if ring { Some(100e-6) } else { None };
        }
        if let Some(r) = b.ring_radius.as_mut() {
            let mut um = *r * 1e6;
            ui.label("reference radius (µm)");
            if num(ui, &mut um).changed() {
                *r = um * 1e-6;
            }
        }
    });
    if b.ring_radius.is_some() {
        ui.weak("The mesh is the (r, z) half plane: x = r, y = z, away from the axis. Optical modes are whispering-gallery modes of azimuthal order m = round(n_target k0 R) (perfectly conducting walls; optical PML layers absorb the radiation, their stretch is applied to r and z), the elastic modes have order m_pump + m_Stokes (backward) or the difference (forward). Crystals are averaged over rotations about the ring axis. Viscosity is not used for rings.");
    }
    let list_field = |ui: &mut Ui, bufs: &mut std::collections::HashMap<(Physics, usize), String>, key: usize, value: &mut Vec<usize>, label: &str| {
        let as_i32: Vec<i32> = value.iter().map(|&x| x as i32).collect();
        let shown = labels_text(&as_i32);
        let entry = bufs.entry((Physics::Brillouin, key)).or_insert_with(|| shown.clone());
        let id = ui.make_persistent_id(("brill_list", key));
        let focused = ui.memory(|m| m.has_focus(id));
        if parse_labels(entry).as_ref() != Some(&as_i32) && !focused {
            *entry = shown;
        }
        ui.horizontal(|ui| {
            ui.label(label);
            let r = ui.add(egui::TextEdit::singleline(entry).id(id).desired_width(110.0).hint_text("empty = all"));
            if r.changed() {
                if let Some(l) = parse_labels(entry) {
                    *value = l.into_iter().filter(|x| *x >= 1).map(|x| x as usize).collect();
                }
            }
        });
    };
    list_field(ui, bufs, 9200, &mut b.pump_modes, "pump modes (numbers from 1)");
    list_field(ui, bufs, 9201, &mut b.stokes_modes, "Stokes modes (numbers from 1)");
    ui.weak("Every pump/Stokes pair gets its own elastic problem at the phase-matched wave number q. Intra-mode scattering is pump = Stokes. The elastic side (materials, crystal data, PML, fixed boundaries, domains without solid, number of modes, search frequency) is set in the Materials, Boundary and this Study tab like a Solid mechanics study.");
    ui.add_space(4.0);
    ui.label(RichText::new("Domains").strong());
    ui.weak("n: refractive index (air 1). Photoelastic: the domains that couple light and strain. Moving boundary: the domains whose outline moves (the body). Optical PML: absorbing layers for the light.");
    for &d in domains {
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("domain {d}")).strong());
                let mut n = b.refractive_index.get(&d).copied().unwrap_or(1.0);
                ui.label("n");
                if num(ui, &mut n).changed() {
                    b.refractive_index.insert(d, n);
                }
                let mut moving = b.moving.contains(&d);
                if ui.checkbox(&mut moving, "moving boundary").changed() {
                    if moving {
                        b.moving.push(d);
                        b.moving.sort();
                    } else {
                        b.moving.retain(|x| *x != d);
                    }
                }
            });
            let cur = b.photoelastic.get(&d).map(|p| photo_kind_name(&p.kind)).unwrap_or("none");
            ui.horizontal(|ui| {
                ui.label("photoelastic");
                egui::ComboBox::from_id_salt(format!("pe_{d}")).selected_text(cur).show_ui(ui, |ui| {
                    if ui.selectable_label(!b.photoelastic.contains_key(&d), "none").clicked() {
                        b.photoelastic.remove(&d);
                    }
                    for (name, make) in PHOTO_PRESETS {
                        if ui.selectable_label(false, *name).clicked() {
                            let follow = b.photoelastic.get(&d).is_none_or(|p| p.follow_stiffness);
                            let euler = b.photoelastic.get(&d).map_or([0.0; 3], |p| p.euler_deg);
                            b.photoelastic.insert(d, PhotoSpec { kind: make(), euler_deg: euler, follow_stiffness: follow });
                        }
                    }
                });
            });
            if let Some(p) = b.photoelastic.get_mut(&d) {
                ui.horizontal_wrapped(|ui| {
                    for (n, v) in photo_constants(&mut p.kind) {
                        ui.label(n);
                        num(ui, v);
                    }
                });
                ui.horizontal(|ui| {
                    ui.checkbox(&mut p.follow_stiffness, "orientation of the elastic crystal data");
                    if !p.follow_stiffness {
                        ui.label("Euler Z-X-Z (°)");
                        for a in p.euler_deg.iter_mut() {
                            ui.add(DragValue::new(a).speed(1.0).range(-360.0..=360.0));
                        }
                    }
                });
            }
            let mut pml = b.optical_pml.get(&d).cloned().unwrap_or_default();
            let active = pml.x != PmlDir::Off || pml.y != PmlDir::Off;
            let mut edited = false;
            egui::CollapsingHeader::new(if active { "Optical PML (active)" } else { "Optical PML" }).id_salt(format!("brill_pml_{d}")).default_open(active).show(ui, |ui| {
                for (axis, dir) in [("x", &mut pml.x), ("y", &mut pml.y)] {
                    ui.horizontal(|ui| {
                        ui.label(format!("absorb along {axis}"));
                        for (v, name) in [(PmlDir::Off, "off"), (PmlDir::Plus, "towards +"), (PmlDir::Minus, "towards −")] {
                            if ui.selectable_label(*dir == v, name).clicked() && *dir != v {
                                *dir = v;
                                edited = true;
                            }
                        }
                    });
                }
                ui.horizontal(|ui| {
                    ui.label("strength σ");
                    edited |= num(ui, &mut pml.strength).changed();
                    ui.label("order p");
                    edited |= num(ui, &mut pml.order).changed();
                });
            });
            if edited {
                if pml.x == PmlDir::Off && pml.y == PmlDir::Off {
                    b.optical_pml.remove(&d);
                } else {
                    b.optical_pml.insert(d, pml);
                }
            }
        });
    }
    ui.add_space(4.0);
    ui.label(RichText::new("Gain").strong());
    egui::Grid::new("brill_gain").num_columns(2).show(ui, |ui| {
        ui.label("mechanical Q cap");
        num(ui, &mut b.q_cap);
        ui.end_row();
        ui.label("waveguide length for g₀ (m)");
        num(ui, &mut b.length);
        ui.end_row();
    });
    ui.horizontal(|ui| {
        ui.label("mode parity under the mirror");
        for (v, name) in [(None, "off"), (Some(0usize), "x → −x"), (Some(1usize), "y → −y")] {
            ui.selectable_value(&mut b.mirror_axis, v, name);
        }
    });
    ui.weak("For cross-sections whose mesh vertices are symmetric about the middle (structured or mirrored meshes): +1 even, −1 odd, in between for mixed (degenerate) modes.");
    ui.checkbox(&mut b.conjugate_strain, "phase-matched photoelastic overlap (conjugate strain)");
    ui.weak("Off uses the strain as it is (the same expression for a phonon of the opposite wave-number sign); the two differ when the strain components have different phases. Gain = 2 ω_p (Q/k_eff) |f|², power P = 2∫Re(E*×H)·z, Q = min(Q_eig, cap).");
}

fn solids_editor(ui: &mut Ui, spec: &mut femcore::mesher3d::Spec3d) {
    use femcore::mesher::Shape;
    use femcore::mesher3d::{Body3d, Solid};
    ui.horizontal(|ui| {
        ui.label("element size h");
        num(ui, &mut spec.h);
        ui.label("alpha");
        ui.add(DragValue::new(&mut spec.alpha).speed(0.005).range(0.05..=0.45));
        ui.checkbox(&mut spec.snap_features, "fit sharp edges");
        ui.label("min layers");
        ui.add(DragValue::new(&mut spec.min_layers).speed(0.1).range(0.0..=10.0));
    });
    ui.weak("Bodies are painted in order: a later body overrides the earlier ones where they overlap; domain 0 makes a hole. The mesh is a tetrahedral lattice cut along the surfaces (very fast); interfaces between domains and sharp edges are fitted exactly. A body with its own element size is refined locally; a body thinner than min layers elements makes the lattice finer across it.");
    ui.weak("Boundary numbers: box 1 x=min, 2 x=max, 3 y=min, 4 y=max, 5 z=min, 6 z=max; cylinder 1 side, 2 base, 3 top; sphere 1; each body continues after the previous one.");
    let mut remove = None;
    let mut swap = None;
    let nb = spec.bodies.len();
    let spec_h = spec.h;
    for (i, body) in spec.bodies.iter_mut().enumerate() {
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("{}. {}", i + 1, body.solid.name())).strong());
                ui.label("domain");
                ui.add(DragValue::new(&mut body.domain).range(0..=999));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("✕").clicked() {
                        remove = Some(i);
                    }
                    if i + 1 < nb && ui.small_button("↓").clicked() {
                        swap = Some((i, i + 1));
                    }
                    if i > 0 && ui.small_button("↑").clicked() {
                        swap = Some((i - 1, i));
                    }
                });
            });
            if body.domain > 0 {
                ui.horizontal(|ui| {
                    let mut local = body.h.is_some();
                    if ui.checkbox(&mut local, "own element size").changed() {
                        body.h = if local { Some(spec_h * 0.5) } else { None };
                    }
                    if let Some(h) = body.h.as_mut() {
                        num(ui, h);
                    }
                });
            }
            match &mut body.solid {
                Solid::Box { min, max } => {
                    for k in 0..3 {
                        ui.horizontal(|ui| {
                            ui.label(["x", "y", "z"][k]);
                            num(ui, &mut min[k]);
                            num(ui, &mut max[k]);
                        });
                    }
                }
                Solid::Sphere { c, r } => {
                    ui.horizontal(|ui| {
                        ui.label("centre");
                        for v in c.iter_mut() {
                            num(ui, v);
                        }
                        ui.label("radius");
                        num(ui, r);
                    });
                }
                Solid::Cylinder { base, axis, r, len } => {
                    ui.horizontal(|ui| {
                        ui.label("base centre");
                        for v in base.iter_mut() {
                            num(ui, v);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("axis");
                        for (k, n) in ["x", "y", "z"].iter().enumerate() {
                            ui.selectable_value(axis, k, *n);
                        }
                        ui.label("radius");
                        num(ui, r);
                        ui.label("length");
                        num(ui, len);
                    });
                }
                Solid::Extrusion { profile, z0, z1 } => {
                    ui.horizontal(|ui| {
                        ui.label("z");
                        num(ui, z0);
                        num(ui, z1);
                    });
                    match profile {
                        Shape::Rect { x0, y0, x1, y1 } => {
                            ui.horizontal(|ui| {
                                ui.label("x");
                                num(ui, x0);
                                num(ui, x1);
                                ui.label("y");
                                num(ui, y0);
                                num(ui, y1);
                            });
                        }
                        Shape::Circle { cx, cy, r } => {
                            ui.horizontal(|ui| {
                                ui.label("centre");
                                num(ui, cx);
                                num(ui, cy);
                                ui.label("radius");
                                num(ui, r);
                            });
                        }
                        Shape::Polygon(v) => {
                            let mut rm = None;
                            for (k, p) in v.iter_mut().enumerate() {
                                ui.horizontal(|ui| {
                                    ui.label(format!("vertex {}", k + 1));
                                    num(ui, &mut p[0]);
                                    num(ui, &mut p[1]);
                                    if ui.small_button("✕").clicked() {
                                        rm = Some(k);
                                    }
                                });
                            }
                            if let Some(k) = rm {
                                if v.len() > 3 {
                                    v.remove(k);
                                }
                            }
                            if ui.small_button("+ vertex").clicked() {
                                let last = *v.last().unwrap_or(&[0.0, 0.0]);
                                v.push([last[0] + 0.1, last[1] + 0.1]);
                            }
                        }
                    }
                }
            }
        });
    }
    if let Some(i) = remove {
        spec.bodies.remove(i);
    }
    if let Some((a, b)) = swap {
        spec.bodies.swap(a, b);
    }
    ui.horizontal(|ui| {
        let next = spec.bodies.len() as i32 + 1;
        if ui.button("+ Box").clicked() {
            spec.bodies.push(Body3d { solid: Solid::Box { min: [0.0, 0.0, 0.0], max: [0.5, 0.5, 0.5] }, domain: next, h: None });
        }
        if ui.button("+ Sphere").clicked() {
            spec.bodies.push(Body3d { solid: Solid::Sphere { c: [0.5, 0.25, 0.125], r: 0.1 }, domain: next, h: None });
        }
        if ui.button("+ Cylinder").clicked() {
            spec.bodies.push(Body3d { solid: Solid::Cylinder { base: [0.5, 0.25, -0.1], axis: 2, r: 0.1, len: 0.5 }, domain: next, h: None });
        }
        if ui.button("+ Extrusion").clicked() {
            spec.bodies.push(Body3d { solid: Solid::Extrusion { profile: Shape::Polygon(vec![[0.0, 0.0], [0.5, 0.0], [0.25, 0.4]]), z0: 0.0, z1: 0.2 }, domain: next, h: None });
        }
    });
}

fn shapes_editor(ui: &mut Ui, spec: &mut femcore::mesher::Spec2d, extrude: &mut Option<Vec<femcore::mesher::Layer>>, wafer: bool) {
    let second = if wafer { "z" } else { "y" };
    use femcore::mesher::{Part, Shape};
    ui.horizontal(|ui| {
        ui.label("element size h");
        num(ui, &mut spec.h);
        ui.label("min angle °");
        ui.add(DragValue::new(&mut spec.min_angle_deg).range(0.0..=33.0));
    });
    ui.weak("Shapes are painted in order: a later shape overrides the earlier ones where they overlap; domain 0 makes a hole.");
    let mut remove = None;
    let mut swap = None;
    let np = spec.parts.len();
    for (i, part) in spec.parts.iter_mut().enumerate() {
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("{}. {}", i + 1, part.shape.name())).strong());
                ui.label("domain");
                ui.add(DragValue::new(&mut part.domain).range(0..=999));
                let mut own = part.h.is_some();
                if ui.checkbox(&mut own, "own h").changed() {
                    part.h = if own { Some(spec.h / 2.0) } else { None };
                }
                if let Some(h) = &mut part.h {
                    num(ui, h);
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("✕").clicked() {
                        remove = Some(i);
                    }
                    if i + 1 < np && ui.small_button("↓").clicked() {
                        swap = Some((i, i + 1));
                    }
                    if i > 0 && ui.small_button("↑").clicked() {
                        swap = Some((i - 1, i));
                    }
                });
            });
            match &mut part.shape {
                Shape::Rect { x0, y0, x1, y1 } => {
                    ui.horizontal(|ui| {
                        ui.label("x");
                        num(ui, x0);
                        num(ui, x1);
                        ui.label(second);
                        num(ui, y0);
                        num(ui, y1);
                    });
                }
                Shape::Circle { cx, cy, r } => {
                    ui.horizontal(|ui| {
                        ui.label("centre");
                        num(ui, cx);
                        num(ui, cy);
                        ui.label("radius");
                        num(ui, r);
                    });
                }
                Shape::Polygon(v) => {
                    let mut rm = None;
                    for (k, p) in v.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("vertex {}", k + 1));
                            num(ui, &mut p[0]);
                            num(ui, &mut p[1]);
                            if v_len_ok(k) && ui.small_button("✕").clicked() {
                                rm = Some(k);
                            }
                        });
                    }
                    if let Some(k) = rm {
                        if v.len() > 3 {
                            v.remove(k);
                        }
                    }
                    if ui.small_button("+ vertex").clicked() {
                        let last = *v.last().unwrap_or(&[0.0, 0.0]);
                        v.push([last[0] + 0.1, last[1] + 0.1]);
                    }
                }
            }
        });
    }
    if let Some(i) = remove {
        spec.parts.remove(i);
    }
    if let Some((a, b)) = swap {
        spec.parts.swap(a, b);
    }
    ui.horizontal(|ui| {
        if ui.button("+ Rectangle").clicked() {
            spec.parts.push(Part { shape: Shape::Rect { x0: 0.0, y0: 0.0, x1: 0.5, y1: 0.5 }, domain: spec.parts.len() as i32 + 1, h: None });
        }
        if ui.button("+ Circle").clicked() {
            spec.parts.push(Part { shape: Shape::Circle { cx: 0.5, cy: 0.25, r: 0.1 }, domain: spec.parts.len() as i32 + 1, h: None });
        }
        if ui.button("+ Polygon").clicked() {
            spec.parts.push(Part { shape: Shape::Polygon(vec![[0.0, 0.0], [0.5, 0.0], [0.25, 0.4]]), domain: spec.parts.len() as i32 + 1, h: None });
        }
    });
    extrude_editor(ui, extrude);
}

/// Layers along y (the propagation axis) of a swept cross-section.
fn extrude_editor(ui: &mut Ui, extrude: &mut Option<Vec<femcore::mesher::Layer>>) {
    use femcore::mesher::Layer;
    let mut on = extrude.is_some();
    if ui.checkbox(&mut on, "Sweep along y into 3D layers").changed() {
        *extrude = if on { Some(vec![Layer { thickness: 0.1, n: 2, domain_offset: 0 }]) } else { None };
    }
    if let Some(layers) = extrude {
        let mut rm = None;
        for (i, l) in layers.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("layer {}", i + 1));
                ui.label("thickness");
                num(ui, &mut l.thickness);
                ui.label("elements");
                ui.add(DragValue::new(&mut l.n).range(1..=200));
                ui.label("domain +");
                ui.add(DragValue::new(&mut l.domain_offset).range(0..=900));
                if ui.small_button("✕").clicked() {
                    rm = Some(i);
                }
            });
        }
        if let Some(i) = rm {
            if layers.len() > 1 {
                layers.remove(i);
            }
        }
        if ui.small_button("+ layer").clicked() {
            layers.push(Layer { thickness: 0.1, n: 2, domain_offset: 0 });
        }
        ui.weak("Side boundaries keep the 2D numbers; the end at y = 0 is the next number, the far end the one after.");
    }
}

/// A `geom` model: element size, sweep, the names of its regions and edges with the integers the solvers use, and the feature list as JSON.
fn model_editor(ui: &mut Ui, list: &mut geom::model::FeatureList, h: &mut f64, extrude: &mut Option<Vec<femcore::mesher::Layer>>) {
    ui.horizontal(|ui| {
        ui.label("element size h (m)");
        num(ui, h);
        ui.weak(format!("{}D model, x-z cross-section; the size rules of the model refine it", list.dim));
    });
    extrude_editor(ui, extrude);
    ui.collapsing("Names and numbers", |ui| match femcore::geomdrive::evaluate(list).and_then(|ev| ev.label_map().map_err(|e| e.to_string())) {
        Ok(lm) => {
            ui.weak("Problems refer to these names (materials, boundary conditions, selectors); the integers are what the solvers store. Read-only.");
            egui::Grid::new("model_domains").num_columns(2).show(ui, |ui| {
                ui.label(RichText::new("domain").strong());
                ui.label(RichText::new("number").strong());
                ui.end_row();
                for d in &lm.domains {
                    ui.label(&d.name);
                    ui.label(d.id.to_string());
                    ui.end_row();
                }
            });
            ui.add_space(4.0);
            egui::Grid::new("model_boundaries").num_columns(2).show(ui, |ui| {
                ui.label(RichText::new("boundary").strong());
                ui.label(RichText::new("number").strong());
                ui.end_row();
                for b in &lm.boundaries {
                    ui.label(&b.name);
                    ui.label(b.id.to_string());
                    ui.end_row();
                }
            });
        }
        Err(e) => {
            ui.colored_label(Color32::from_rgb(220, 100, 100), e);
        }
    });
    ui.collapsing("Feature list (JSON)", |ui| {
        let id = ui.next_auto_id().with("model_json");
        let mut text: String = ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| serde_json::to_string_pretty(&*list).unwrap_or_default());
        ui.add(egui::TextEdit::multiline(&mut text).code_editor().desired_rows(12).desired_width(f32::INFINITY));
        let parsed = serde_json::from_str::<geom::model::FeatureList>(&text);
        ui.horizontal(|ui| {
            if ui.add_enabled(parsed.is_ok(), egui::Button::new("Apply")).clicked() {
                if let Ok(l) = &parsed {
                    *list = l.clone();
                }
                ui.data_mut(|d| d.remove::<String>(id));
                return;
            }
            if ui.button("Reset").clicked() {
                ui.data_mut(|d| d.remove::<String>(id));
                return;
            }
            if let Err(e) = &parsed {
                ui.colored_label(Color32::from_rgb(220, 100, 100), e.to_string());
            }
        });
        ui.data_mut(|d| d.insert_temp(id, text));
    });
}

fn v_len_ok(_k: usize) -> bool {
    true
}

impl App {
    pub fn domain_labels(&self) -> Vec<i32> {
        if let Some(m) = &self.doc.mesh {
            m.domain_labels()
        } else if let Some(v) = &self.doc.view {
            v.domain_labels.clone()
        } else {
            vec![1]
        }
    }

    fn set_material(&mut self, domain: i32, key: &str, v: f64) {
        for s in &mut self.proj.specs {
            s.materials.set(domain, key, v);
        }
    }

    // ------------------------------------------------------------------ inspector

    pub fn inspector(&mut self, ui: &mut Ui) {
        ui.horizontal_wrapped(|ui| {
            for t in Tab::ALL {
                if t == Tab::Results && self.results.is_empty() {
                    continue;
                }
                if t == Tab::Analysis && !self.has_analysis() {
                    continue;
                }
                if ui.selectable_label(self.tab == t, RichText::new(t.name()).strong()).clicked() {
                    self.tab = t;
                    if t != Tab::Physics {
                        self.picking = None;
                    }
                }
            }
        });
        ui.separator();
        egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| match self.tab {
            Tab::Geometry => self.tab_geometry(ui),
            Tab::Materials => self.tab_materials(ui),
            Tab::Physics => self.tab_physics(ui),
            Tab::Study => self.tab_study(ui),
            Tab::Results => self.tab_results(ui),
            Tab::Analysis => self.tab_analysis(ui),
        });
    }

    // ------------------------------------------------------------------ geometry

    fn tab_geometry(&mut self, ui: &mut Ui) {
        ui.heading("Geometry");
        {
            let wafer = self.proj.active_spec().materials.frame == Frame::Wafer;
            let mut w = wafer;
            ui.horizontal(|ui| {
                ui.label("Axes");
                ui.selectable_value(&mut w, true, "wafer: x lateral, y propagation, z normal (2D = x-z)");
                ui.selectable_value(&mut w, false, "legacy model axes: 2D = x-y, z out of plane");
            });
            if w != wafer {
                for sp in &mut self.proj.specs {
                    sp.materials.frame = if w { Frame::Wafer } else { Frame::Model };
                }
            }
            ui.weak(if w { "Crystal orientations (Euler angles), photoelastic tensors, permittivities, constraint and load components and the names of result fields use the right-handed axes x, y, z of the wafer; a 2D cross-section is drawn in x-z and y (the direction of propagation) points into the page. A wave travelling towards the viewer travels along -y." } else { "Model axes: a 2D mesh is the x-y plane and z is the out-of-plane direction." });
        }
        let mut g = self.geom_edit.clone();
        let dim = g.geometry.dim();
        ui.horizontal(|ui| {
            ui.label("Shape");
            if ui.selectable_label(dim == 2, "Rectangle (2D)").clicked() && dim != 2 {
                g = GeometrySpec::rectangle(1.0, 0.5, 40, 20);
            }
            if ui.selectable_label(matches!(g.geometry, Geometry::Block { .. }), "Block (3D)").clicked() && !matches!(g.geometry, Geometry::Block { .. }) {
                g = GeometrySpec::block(1.0, 0.1, 0.1, 16, 3, 3);
            }
            if ui.selectable_label(matches!(g.geometry, Geometry::Solids { .. }), "Solids (3D mesher)").clicked() && !matches!(g.geometry, Geometry::Solids { .. }) {
                g = GeometrySpec::solids();
            }
            if ui.selectable_label(matches!(g.geometry, Geometry::Shapes { .. }), "Shapes (mesher)").clicked() && !matches!(g.geometry, Geometry::Shapes { .. }) {
                g = GeometrySpec::shapes();
            }
            if ui.selectable_label(matches!(g.geometry, Geometry::Model { .. }), "Model (geom)").clicked() && !matches!(g.geometry, Geometry::Model { .. }) {
                let mut m = geom::model::Model::new_2d();
                let r = m.rect("plate", (0.0, 0.0), (1.0, 0.5));
                m.overlay(&[("plate", r)]);
                g = GeometrySpec::model(m.list, 0.05);
            }
        });
        egui::Grid::new("geom_grid").num_columns(4).show(ui, |ui| match &mut g.geometry {
            Geometry::Rectangle { lx, lz, nx, nz } => {
                ui.label("width x (m)");
                num(ui, lx);
                ui.label("elements");
                ui.add(DragValue::new(nx).range(1..=2000));
                ui.end_row();
                ui.label("height z (m)");
                num(ui, lz);
                ui.label("elements");
                ui.add(DragValue::new(nz).range(1..=2000));
                ui.end_row();
            }
            Geometry::Shapes { .. } | Geometry::Solids { .. } | Geometry::Model { .. } => {}
            Geometry::Block { lx, ly, lz, nx, ny, nz } => {
                for (name, l, n) in [("length x (m)", lx, nx), ("width y (m)", ly, ny), ("height z (m)", lz, nz)] {
                    ui.label(name);
                    num(ui, l);
                    ui.label("elements");
                    ui.add(DragValue::new(n).range(1..=200));
                    ui.end_row();
                }
            }
        });
        if let Geometry::Shapes { spec, extrude } = &mut g.geometry {
            shapes_editor(ui, spec, extrude, self.proj.active_spec().materials.frame == Frame::Wafer);
        }
        if let Geometry::Solids { spec } = &mut g.geometry {
            solids_editor(ui, spec);
        }
        if let Geometry::Model { list, h, extrude } = &mut g.geometry {
            model_editor(ui, list, h, extrude);
        }
        ui.label(RichText::new("Boundary numbers: rectangle 1 x=0, 2 z=0, 3 z=lz, 4 x=lx; block 1 x=0, 2 y=0, 3 z=0, 4 z=lz, 5 y=ly, 6 x=lx.").weak().small());
        ui.add_space(6.0);

        // regions
        for (title, which) in [("Domain regions (relabel cells whose centre lies in a box)", 0), ("Boundary regions (relabel boundary facets whose midpoint lies in a box)", 1)] {
            ui.collapsing(title, |ui| {
                let list = if which == 0 { &mut g.domains } else { &mut g.boundaries };
                let mut remove = None;
                for (i, r) in list.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        ui.label(if which == 0 { "domain" } else { "boundary" });
                        ui.add(DragValue::new(&mut r.label).range(1..=999));
                        if ui.small_button("✕").on_hover_text("remove").clicked() {
                            remove = Some(i);
                        }
                    });
                    egui::Grid::new(format!("region_{which}_{i}")).num_columns(5).show(ui, |ui| {
                        for k in 0..dim {
                            ui.label(["x", "y", "z"][k]);
                            ui.label("from");
                            num(ui, &mut r.bbox[k]);
                            ui.label("to");
                            num(ui, &mut r.bbox[dim + k]);
                            ui.end_row();
                        }
                    });
                    ui.separator();
                }
                if let Some(i) = remove {
                    list.remove(i);
                }
                if ui.button("+ add region").clicked() {
                    let (lo, hi) = match &g.geometry {
                        Geometry::Rectangle { lx, lz, .. } => (vec![0.0, 0.0], vec![*lx, *lz]),
                        Geometry::Block { lx, ly, lz, .. } => (vec![0.0, 0.0, 0.0], vec![*lx, *ly, *lz]),
                        Geometry::Shapes { .. } | Geometry::Solids { .. } | Geometry::Model { .. } => (vec![0.0; dim], vec![1.0; dim]),
                    };
                    let mut b = lo.clone();
                    b.extend(hi.iter().enumerate().map(|(k, h)| if k == 0 { 0.5 * h } else { *h }));
                    let list = if which == 0 { &mut g.domains } else { &mut g.boundaries };
                    let next = if which == 0 { 2 + list.len() as i32 } else { 5 + list.len() as i32 };
                    list.push(Region { bbox: b, label: next });
                }
            });
        }
        let changed = g != self.geom_edit;
        self.geom_edit = g;
        if changed {
            self.geometry_dirty = true;
        }
        ui.add_space(6.0);
        let dirty = self.geometry_dirty || self.geom_edit != self.proj.geometry;
        ui.horizontal(|ui| {
            let b = ui.add_enabled(self.busy.is_none(), egui::Button::new(RichText::new("Build mesh").strong()).fill(if dirty { Color32::from_rgb(40, 90, 150) } else { Color32::from_rgb(52, 56, 66) }));
            if b.clicked() {
                self.rebuild_mesh();
            }
            if dirty {
                ui.weak("(geometry changed)");
            }
        });
        let elements = match &self.geom_edit.geometry {
            Geometry::Rectangle { nx, nz, .. } => 2 * nx * nz,
            Geometry::Block { nx, ny, nz, .. } => 6 * nx * ny * nz,
            Geometry::Shapes { .. } | Geometry::Solids { .. } | Geometry::Model { .. } => 0,
        };
        if elements > 0 {
            ui.weak(format!("{elements} elements when built"));
        } else {
            ui.weak("the element count follows from the sizes when the mesh is built");
        }

        ui.add_space(8.0);
        ui.separator();
        ui.heading("Mesh");
        if let Some(v) = self.doc.view.clone() {
            ui.label(format!("{}D, {} points, {} cells, {} boundary facets/segments", v.dim, v.points.len(), v.n_cells, if v.dim == 2 { v.segs.len() } else { v.tris.len() }));
            ui.label(format!("domains {:?}", v.domain_labels));
            ui.label(format!("boundaries {:?}", v.boundary_labels));
            ui.label(format!("element types: {}", v.types.join(", ")));
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.show_edges, "edges");
                ui.checkbox(&mut self.show_numbers, "boundary numbers");
            });
        } else {
            ui.weak("no mesh yet");
        }
        ui.add_space(8.0);
        ui.separator();
        ui.heading("Examples");
        for (i, n) in crate::examples::NAMES.iter().enumerate() {
            if ui.link(*n).clicked() {
                self.load_example(i);
            }
        }
    }

    // ------------------------------------------------------------------ materials

    fn tab_materials(&mut self, ui: &mut Ui) {
        let phys = self.proj.active;
        ui.heading("Materials");
        ui.label(format!("Properties used by {} for each domain. Values are shared between physics that use the same property.", phys.name()));
        let props = props_for(phys);
        let mesh_dim = self.doc.mesh.as_ref().map(|m| m.dim).or_else(|| self.doc.view.as_ref().map(|v| v.dim)).unwrap_or(2);
        let domains = self.domain_labels();
        let mut apply_all: Option<usize> = None;
        ui.horizontal(|ui| {
            ui.label("All domains:");
            egui::ComboBox::from_id_salt("mat_all").selected_text("apply material...").show_ui(ui, |ui| {
                for (i, m) in MATERIALS.iter().enumerate() {
                    if ui.selectable_label(false, m.name).clicked() {
                        apply_all = Some(i);
                    }
                }
            });
        });
        let mut edits: Vec<(i32, &'static str, f64)> = vec![];
        let mut preset: Vec<(i32, usize)> = vec![];
        let mut aniso_edits: Vec<(i32, Option<Anisotropy>)> = vec![];
        let mut pml_edits: Vec<(i32, Option<Pml>)> = vec![];
        let mut crystal: Vec<(i32, usize)> = vec![];
        for &d in &domains {
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("Domain {d}")).strong());
                    egui::ComboBox::from_id_salt(format!("mat_{d}")).selected_text("material...").show_ui(ui, |ui| {
                        for (i, m) in MATERIALS.iter().enumerate() {
                            if ui.selectable_label(false, m.name).clicked() {
                                preset.push((d, i));
                            }
                        }
                        ui.separator();
                        for (i, c) in CRYSTALS.iter().enumerate() {
                            if ui.selectable_label(false, c.name).clicked() {
                                crystal.push((d, i));
                            }
                        }
                    });
                });
                egui::Grid::new(format!("props_{d}")).num_columns(3).show(ui, |ui| {
                    for (key, label, unit) in props {
                        let mut v = self.proj.active_spec().materials.get(d, key).unwrap_or_else(|| default_prop(key));
                        ui.label(*label);
                        if num(ui, &mut v).changed() {
                            edits.push((d, key, v));
                        }
                        ui.weak(*unit);
                        ui.end_row();
                    }
                });
                let mats = self.proj.active_spec().materials.clone();
                if phys == Physics::Emw || phys == Physics::Solid || phys == Physics::Brillouin {
                    let mut pml = mats.pml.get(&d).cloned().unwrap_or_default();
                    let active = pml.x != PmlDir::Off || pml.y != PmlDir::Off || pml.z != PmlDir::Off;
                    let mut edited = false;
                    egui::CollapsingHeader::new(if active { "Perfectly matched layer (active)" } else { "Perfectly matched layer" }).id_salt(format!("pml_{d}")).default_open(active).show(ui, |ui| {
                        let nax = if mesh_dim == 3 { 3 } else { 2 };
                        let second = if mesh_dim == 2 && mats.frame == Frame::Wafer { "z" } else { "y" };
                        for (axis, dir) in [("x", &mut pml.x), (second, &mut pml.y), ("z", &mut pml.z)].into_iter().take(nax) {
                            ui.horizontal(|ui| {
                                ui.label(format!("absorb along {axis}"));
                                for (v, name) in [(PmlDir::Off, "off"), (PmlDir::Plus, "towards +"), (PmlDir::Minus, "towards −")] {
                                    if ui.selectable_label(*dir == v, name).clicked() && *dir != v {
                                        *dir = v;
                                        edited = true;
                                    }
                                }
                            });
                        }
                        ui.horizontal(|ui| {
                            ui.label("strength σ");
                            edited |= num(ui, &mut pml.strength).changed();
                            ui.label("order p");
                            edited |= num(ui, &mut pml.order).changed();
                            ui.label("real stretch");
                            edited |= num(ui, &mut pml.re).changed();
                        });
                        ui.weak(if phys != Physics::Emw { "Stretch s = re − jσ·t^p, t from 0 at the interface to 1 at the outer edge of this domain. Used by the Eigenfrequency study (complex frequencies, Q); fix the outer edge of the layer." } else { "Stretch s = re − jσ·t^p, t from 0 at the interface to 1 at the outer edge of this domain. The layer covers the whole domain: make it a separate domain region and put the PEC wall at its outer edge." });
                    });
                    if edited {
                        let off = pml.x == PmlDir::Off && pml.y == PmlDir::Off && pml.z == PmlDir::Off;
                        pml_edits.push((d, if off { None } else { Some(pml) }));
                    }
                }
                let mut an = mats.aniso.get(&d).cloned().unwrap_or_default();
                let has = mats.aniso.contains_key(&d);
                egui::CollapsingHeader::new(if has { "Anisotropy / crystal data (active)" } else { "Anisotropy / crystal data" }).id_salt(format!("an_{d}")).default_open(has).show(ui, |ui| {
                    if aniso_editor(ui, phys, d, &mut an, &mats) {
                        let empty = an.stiffness.is_none() && an.conductivity.is_none() && an.permittivity.is_none() && an.expansion.is_none();
                        aniso_edits.push((d, if empty { None } else { Some(an.clone()) }));
                    }
                });
            });
        }
        for (d, p) in pml_edits {
            for sp in &mut self.proj.specs {
                match &p {
                    Some(p) => {
                        sp.materials.pml.insert(d, p.clone());
                    }
                    None => {
                        sp.materials.pml.remove(&d);
                    }
                }
            }
        }
        for (d, a) in aniso_edits {
            for sp in &mut self.proj.specs {
                match &a {
                    Some(a) => {
                        sp.materials.aniso.insert(d, a.clone());
                    }
                    None => {
                        sp.materials.aniso.remove(&d);
                    }
                }
            }
        }
        for (d, i) in crystal {
            let (a, rho) = (CRYSTALS[i].make)();
            for sp in &mut self.proj.specs {
                sp.materials.aniso.insert(d, a.clone());
                sp.materials.set(d, "rho", rho);
            }
        }
        for (d, k, v) in edits {
            self.set_material(d, k, v);
        }
        for (d, i) in preset {
            for (k, v) in MATERIALS[i].props {
                self.set_material(d, k, *v);
            }
        }
        if let Some(i) = apply_all {
            for &d in &domains {
                for (k, v) in MATERIALS[i].props {
                    self.set_material(d, k, *v);
                }
            }
        }
        ui.add_space(6.0);
        ui.weak("Library values are typical room-temperature numbers; check them for your application. Anisotropic data replaces the scalar value it covers in that domain.");
    }

    // ------------------------------------------------------------------ physics

    fn tab_physics(&mut self, ui: &mut Ui) {
        let phys = self.proj.active;
        let mut spec = self.proj.active_spec().clone();
        ui.heading(phys.name());
        ui.horizontal(|ui| {
            ui.label("Element order");
            ui.selectable_value(&mut spec.order, 1, "P1 (linear)");
            ui.selectable_value(&mut spec.order, 2, "P2 (quadratic, default)");
        });
        if phys != Physics::Emw {
            ui.horizontal(|ui| {
                ui.label("Linear solver");
                egui::ComboBox::from_id_salt("linsolver").selected_text(spec.solver.name()).show_ui(ui, |ui| {
                    for sv in [LinearSolver::Auto, LinearSolver::Direct, LinearSolver::Iterative] {
                        ui.selectable_value(&mut spec.solver, sv, sv.name());
                    }
                });
            });
            ui.weak("Automatic uses the sparse direct solver for small problems and algebraic multigrid (AMG-CG) for large ones: from about 30 000 free dofs in 3D, 400 000 in 2D. It applies to static solid, stationary and transient heat, electrostatics, and the lowest eigenfrequencies (LOBPCG with a multigrid preconditioner, chosen when multigrid converges fast; otherwise a factorisation). If multigrid fails to converge, the direct solver is used.");
        }
        let dim = self.doc.mesh.as_ref().map(|m| m.dim).or_else(|| self.doc.view.as_ref().map(|v| v.dim)).unwrap_or(2);
        if phys == Physics::Solid {
            if dim == 2 {
                ui.horizontal(|ui| {
                    ui.label("2D approximation");
                    ui.selectable_value(&mut spec.plane, Plane::Strain, "plane strain (default)");
                    ui.selectable_value(&mut spec.plane, Plane::Stress, "plane stress");
                });
                ui.horizontal(|ui| {
                    ui.label("thickness (m)");
                    num(ui, &mut spec.thickness);
                });
            }
            let mut on = spec.thermal.is_some();
            if ui.checkbox(&mut on, "Thermal expansion").changed() {
                spec.thermal = if on { Some(ThermalExpansion { source: ThermalSource::Uniform(293.15), t_ref: 293.15 }) } else { None };
            }
            if let Some(th) = &mut spec.thermal {
                let mut from_heat = th.source == ThermalSource::FromHeat;
                ui.horizontal(|ui| {
                    ui.label("temperature from");
                    if ui.selectable_label(!from_heat, "uniform").clicked() {
                        from_heat = false;
                        if matches!(th.source, ThermalSource::FromHeat) {
                            th.source = ThermalSource::Uniform(th.t_ref + 100.0);
                        }
                    }
                    if ui.selectable_label(from_heat, "last heat-transfer result").clicked() {
                        from_heat = true;
                        th.source = ThermalSource::FromHeat;
                    }
                });
                if let ThermalSource::Uniform(t) = &mut th.source {
                    ui.horizontal(|ui| {
                        ui.label("T (K)");
                        num(ui, t);
                    });
                } else if self.solve_state.heat_t.is_none() {
                    ui.colored_label(Color32::from_rgb(235, 170, 90), "No heat result yet: switch Physics to Heat transfer and Compute first.");
                }
                ui.horizontal(|ui| {
                    ui.label("strain reference T_ref (K)");
                    num(ui, &mut th.t_ref);
                });
            }
        }
        ui.add_space(4.0);
        ui.separator();
        ui.label(RichText::new("Boundary conditions").strong());
        let note = match phys {
            Physics::Electrostatics => "Boundaries without a condition have zero normal D (zero charge).",
            Physics::Heat => "Boundaries without a condition are thermally insulated.",
            Physics::Solid => "Boundaries without a condition are free.",
            Physics::Emw => "All outer boundaries are perfect electric conductors (except periodic ones).",
            Physics::Brillouin => "Elastic boundary conditions. Boundaries without a condition are free; fix the outer edges of the substrate and its PML. The optical problem has perfectly conducting walls.",
        };
        ui.weak(note);
        let mut remove = None;
        for i in 0..spec.bcs.len() {
            let buf_key = (phys, i);
            let labels = spec.bcs[i].labels.clone();
            let entry = self.label_buffers.entry(buf_key).or_insert_with(|| labels_text(&labels));
            if parse_labels(entry).as_ref() != Some(&labels) && !ui.memory(|m| m.has_focus(ui.make_persistent_id(("bclabels", phys as u8, i)))) {
                *entry = labels_text(&labels);
            }
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("{}. {}", i + 1, spec.bcs[i].kind.name())).strong().color(crate::colormap::palette(i + 3)));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.small_button("✕").on_hover_text("remove").clicked() {
                            remove = Some(i);
                        }
                        let picking = self.picking == Some(i);
                        if ui.selectable_label(picking, "Pick on mesh").on_hover_text("click boundaries in the viewport to add or remove them").clicked() {
                            self.picking = if picking { None } else { Some(i) };
                            self.show_results = false;
                        }
                    });
                });
                ui.horizontal(|ui| {
                    ui.label("boundaries");
                    let r = ui.add(egui::TextEdit::singleline(entry).id(ui.make_persistent_id(("bclabels", phys as u8, i))).desired_width(120.0).hint_text("e.g. 1, 4 or 2-4"));
                    if r.changed() {
                        if let Some(l) = parse_labels(entry) {
                            spec.bcs[i].labels = l;
                        }
                    }
                    if parse_labels(entry).is_none() {
                        ui.colored_label(Color32::from_rgb(235, 110, 100), "numbers separated by commas");
                    }
                });
                match &mut spec.bcs[i].kind {
                    BcKind::Potential(v) => {
                        ui.horizontal(|ui| {
                            ui.label("V₀ (V)");
                            num(ui, v);
                        });
                    }
                    BcKind::Ground => {}
                    BcKind::SurfaceCharge(s) => {
                        ui.horizontal(|ui| {
                            ui.label("ρ_s (C/m²)");
                            num(ui, s);
                        });
                    }
                    BcKind::Temperature(t) => {
                        ui.horizontal(|ui| {
                            ui.label("T₀ (K)");
                            num(ui, t);
                        });
                    }
                    BcKind::Flux(q) => {
                        ui.horizontal(|ui| {
                            ui.label("inward flux q₀ (W/m²)");
                            num(ui, q);
                        });
                    }
                    BcKind::Convection { h, text } => {
                        ui.horizontal(|ui| {
                            ui.label("h (W/m²·K)");
                            num(ui, h);
                            ui.label("T_ext (K)");
                            num(ui, text);
                        });
                    }
                    BcKind::Fixed { fixed, values } => {
                        for k in 0..dim.min(3) {
                            ui.horizontal(|ui| {
                                ui.checkbox(&mut fixed[k], format!("u{} =", ["x", "y", "z"][k]));
                                ui.add_enabled_ui(fixed[k], |ui| {
                                    num(ui, &mut values[k]);
                                });
                                ui.weak("m");
                            });
                        }
                    }
                    BcKind::Load(f) => {
                        for k in 0..dim.min(3) {
                            ui.horizontal(|ui| {
                                ui.label(format!("f{} (N/m²)", ["x", "y", "z"][k]));
                                num(ui, &mut f[k]);
                            });
                        }
                    }
                    BcKind::LowReflecting => {
                        ui.weak("Absorbs outgoing waves (traction −ρc_p(v·n)n − ρc_s v_t with the speeds of the adjacent material). Used by the Eigenfrequency study, which then gives complex frequencies and Q.");
                    }
                }
            });
        }
        if let Some(i) = remove {
            spec.bcs.remove(i);
            self.label_buffers.clear();
            if self.picking == Some(i) {
                self.picking = None;
            }
        }
        let avail = BcKind::available(phys);
        if !avail.is_empty() {
            ui.menu_button("+ Add boundary condition", |ui| {
                for k in avail {
                    if ui.button(k.name()).clicked() {
                        spec.bcs.push(Bc { labels: vec![], kind: k });
                        self.label_buffers.clear();
                        ui.close();
                    }
                }
            });
        }
        // domains that do not take part, and Floquet-periodic boundary pairs
        if phys == Physics::Solid || phys == Physics::Emw || phys == Physics::Brillouin {
            ui.add_space(6.0);
            ui.separator();
            let text_field = |ui: &mut Ui, bufs: &mut std::collections::HashMap<(Physics, usize), String>, key: (Physics, usize), value: &mut Vec<i32>, id: &str| {
                let shown = labels_text(value);
                let entry = bufs.entry(key).or_insert_with(|| shown.clone());
                let focused = ui.memory(|m| m.has_focus(ui.make_persistent_id((id, key.1))));
                if parse_labels(entry).as_ref() != Some(&*value) && !focused {
                    *entry = shown;
                }
                let r = ui.add(egui::TextEdit::singleline(entry).id(ui.make_persistent_id((id, key.1))).desired_width(110.0).hint_text("e.g. 1, 4"));
                if r.changed() {
                    if let Some(l) = parse_labels(entry) {
                        *value = l;
                    }
                }
            };
            if phys == Physics::Solid || phys == Physics::Brillouin {
                ui.horizontal(|ui| {
                    ui.label("Domains without solid (e.g. air)");
                    text_field(ui, &mut self.label_buffers, (phys, 9000), &mut spec.inactive, "inactive");
                });
            }
            ui.label(RichText::new("Periodic conditions (Floquet)").strong());
            ui.weak("The destination boundaries repeat the source boundaries, shifted by the distance between the two groups, times exp(−i·phase). Both must carry identical meshes. Used by the Eigenfrequency study.");
            let mut remove_p = None;
            for i in 0..spec.periodic.len() {
                ui.group(|ui| {
                    ui.horizontal(|ui| {
                        ui.label(format!("{}. source", i + 1));
                        text_field(ui, &mut self.label_buffers, (phys, 9100 + 2 * i), &mut spec.periodic[i].source, "persrc");
                        ui.label("→ destination");
                        text_field(ui, &mut self.label_buffers, (phys, 9101 + 2 * i), &mut spec.periodic[i].dest, "perdst");
                        if ui.small_button("✕").clicked() {
                            remove_p = Some(i);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("phase (degrees)");
                        let mut deg = spec.periodic[i].phase.to_degrees();
                        if num(ui, &mut deg).changed() {
                            spec.periodic[i].phase = deg.to_radians();
                        }
                    });
                });
            }
            if let Some(i) = remove_p {
                spec.periodic.remove(i);
                self.label_buffers.retain(|k, _| k.1 < 9100);
            }
            if ui.button("+ periodic pair").clicked() {
                spec.periodic.push(Periodic::default());
            }
        }
        if spec != *self.proj.active_spec() {
            *self.proj.spec_mut(phys) = spec;
        }
        if self.picking.is_some() {
            ui.add_space(4.0);
            ui.colored_label(Color32::from_rgb(120, 190, 255), "Picking: click boundaries in the viewport to add or remove them from the selected condition.");
        }
    }

    // ------------------------------------------------------------------ study

    fn tab_study(&mut self, ui: &mut Ui) {
        let phys = self.proj.active;
        let mut spec = self.proj.active_spec().clone();
        ui.heading(format!("Study: {}", phys.name()));
        match phys {
            Physics::Electrostatics => {
                ui.label("Stationary: solves −∇·(ε₀ε_r∇V) = ρ_v. Reports energy, maximum field and, with two different potentials, the capacitance.");
            }
            Physics::Heat => {
                let mut transient = matches!(spec.study, Study::Transient { .. });
                ui.horizontal(|ui| {
                    if ui.selectable_label(!transient, "Stationary").clicked() {
                        transient = false;
                    }
                    if ui.selectable_label(transient, "Transient").clicked() {
                        transient = true;
                    }
                });
                if transient && !matches!(spec.study, Study::Transient { .. }) {
                    spec.study = Study::Transient { t_end: 100.0, n_out: 20, t0: "293.15".into() };
                } else if !transient {
                    spec.study = Study::Stationary;
                }
                if let Study::Transient { t_end, n_out, t0 } = &mut spec.study {
                    egui::Grid::new("tr").num_columns(2).show(ui, |ui| {
                        ui.label("end time (s)");
                        num(ui, t_end);
                        ui.end_row();
                        ui.label("output intervals");
                        ui.add(DragValue::new(n_out).range(1..=500));
                        ui.end_row();
                        ui.label("initial T (K), expression in x, y, z");
                        ui.text_edit_singleline(t0);
                        ui.end_row();
                    });
                    ui.weak("Examples: 293.15   or   293.15 + 80*exp(-((x-0.5)^2)/0.01). Up to about 2500 free dofs the semi-discrete system is solved exactly (no time-stepping error); larger models use BDF2.");
                }
            }
            Physics::Solid => {
                let mut eigen = matches!(spec.study, Study::Eigenfrequency { .. });
                ui.horizontal(|ui| {
                    if ui.selectable_label(!eigen, "Static").clicked() {
                        eigen = false;
                    }
                    if ui.selectable_label(eigen, "Eigenfrequencies").clicked() {
                        eigen = true;
                    }
                });
                if eigen && !matches!(spec.study, Study::Eigenfrequency { .. }) {
                    spec.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
                } else if !eigen {
                    spec.study = Study::Static;
                }
                if let Study::Eigenfrequency { n_modes, f_search } = &mut spec.study {
                    egui::Grid::new("eig").num_columns(2).show(ui, |ui| {
                        ui.label("number of modes");
                        ui.add(DragValue::new(n_modes).range(1..=60));
                        ui.end_row();
                        ui.label("search around (Hz, 0 = lowest)");
                        num(ui, f_search);
                        ui.end_row();
                    });
                    ui.weak("Boundary loads are ignored; only fixed constraints matter. Without a fixed boundary the first modes are rigid-body motions at 0 Hz.");
                    if self.doc.mesh.as_ref().is_some_and(|m| m.dim == 2) {
                        let mut on = spec.out_of_plane.is_some();
                        if ui.checkbox(&mut on, "Out-of-plane wave number (waveguide cross-section)").changed() {
                            spec.out_of_plane = if on { Some(0.0) } else { None };
                        }
                        let wafer = spec.materials.frame == Frame::Wafer;
                        if let Some(kz) = &mut spec.out_of_plane {
                            ui.horizontal(|ui| {
                                ui.label(if wafer { "k (rad/m), propagation along the page normal: positive towards the viewer (-y)" } else { "kz (rad/m)" });
                                num(ui, kz);
                            });
                            ui.weak("Three displacement components on the cross-section, u(x, y) exp(-i kz z); anisotropic crystals use their full stiffness. For backward Brillouin scattering kz = 2 n_eff 2π/λ; for forward scattering kz = n_eff 2π f_acoustic / c.");
                        }
                    }
                }
            }
            Physics::Brillouin => {
                if !matches!(spec.study, Study::Eigenfrequency { .. }) {
                    spec.study = Study::Eigenfrequency { n_modes: 30, f_search: 10e9 };
                }
                if let Study::Eigenfrequency { n_modes, f_search } = &mut spec.study {
                    ui.label(RichText::new("Elastic modes").strong());
                    egui::Grid::new("brill_eig").num_columns(2).show(ui, |ui| {
                        ui.label("modes per pump/Stokes pair");
                        ui.add(DragValue::new(n_modes).range(1..=200));
                        ui.end_row();
                        ui.label("search around (Hz)");
                        num(ui, f_search);
                        ui.end_row();
                    });
                    ui.weak("About 11 GHz for forward-like and 15 to 25 GHz for backward scattering in silicon; the search is a shift-and-invert around this frequency, so several tens of modes are needed to see the guided ones among the substrate modes.");
                }
                if self.doc.mesh.as_ref().is_some_and(|m| m.dim == 3) {
                    ui.colored_label(Color32::from_rgb(235, 110, 100), "The Brillouin study works on a 2D cross-section.");
                }
                let domains = self.domain_labels();
                brillouin_editor(ui, &mut spec, &domains, &mut self.label_buffers);
            }
            Physics::Emw if self.doc.mesh.as_ref().is_some_and(|m| m.dim == 3) => {
                if !matches!(spec.study, Study::Eigenfrequency { .. }) {
                    spec.study = Study::Eigenfrequency { n_modes: 6, f_search: 10e9 };
                }
                if let Study::Eigenfrequency { n_modes, f_search } = &mut spec.study {
                    let mut ghz = *f_search / 1e9;
                    egui::Grid::new("emeig").num_columns(2).show(ui, |ui| {
                        ui.label("number of modes");
                        ui.add(DragValue::new(n_modes).range(1..=40));
                        ui.end_row();
                        ui.label("search around (GHz)");
                        if num(ui, &mut ghz).changed() {
                            *f_search = ghz * 1e9;
                        }
                        ui.end_row();
                    });
                    ui.weak("Eigenmodes of the 3D structure with perfect electric conductor on the outer boundary, found with quadratic Nedelec elements (1 for first order in the Physics tab); PML domains give complex frequencies and Q. Modes nearest the search frequency are returned; the shift must lie in the band you want (modes far above it may be missed).");
                }
            }
            Physics::Emw => {
                if let Study::Modes { freq, n_target, evanescent, n_modes } = &mut spec.study {
                    let mut ghz = *freq / 1e9;
                    egui::Grid::new("modes").num_columns(2).show(ui, |ui| {
                        ui.label("frequency (GHz; 1550 nm = 193414.5)");
                        if num(ui, &mut ghz).changed() {
                            *freq = ghz * 1e9;
                        }
                        ui.end_row();
                        ui.label("search near effective index");
                        num(ui, n_target);
                        ui.end_row();
                        ui.label("number of modes");
                        ui.add(DragValue::new(n_modes).range(1..=40));
                        ui.end_row();
                    });
                    ui.checkbox(evanescent, "target is imaginary (find evanescent modes)");
                    ui.weak("Finds the modes with n_eff closest to the target. For a guide filled with ε_r, n_eff lies between 0 and √ε_r.");
                }
            }
        }
        ui.add_space(6.0);
        ui.separator();
        ui.label(RichText::new("Run folder").strong());
        ui.checkbox(&mut self.autosave, "save every run (results .hdf5 + XDMF + project .nyx.xml) as temp_untitled_<date>_<time>");
        if self.autosave {
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.autosave_dir).desired_width(300.0));
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new().pick_folder() {
                        self.autosave_dir = p.display().to_string();
                    }
                }
            });
        }
        ui.add_space(4.0);
        ui.label(RichText::new("Results file").strong());
        ui.horizontal(|ui| {
            ui.label("Write after Compute");
            egui::ComboBox::from_id_salt("outmode").selected_text(spec.output.mode.name()).show_ui(ui, |ui| {
                for m in [OutputMode::None, OutputMode::Hdf5Serial, OutputMode::Hdf5Parallel] {
                    ui.selectable_value(&mut spec.output.mode, m, m.name());
                }
            });
        });
        if spec.output.mode != OutputMode::None {
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut spec.output.path).hint_text("results.h5").desired_width(240.0));
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("HDF5", &["h5", "hdf5"]).set_file_name("results.h5").save_file() {
                        spec.output.path = p.display().to_string();
                    }
                }
            });
            ui.weak("Mesh, problem description and every solution with all fields. Serial uses the sequential HDF5 driver. Parallel (MPI-IO) is not available in this build.");
        }
        if spec != *self.proj.active_spec() {
            *self.proj.spec_mut(phys) = spec;
        }
        ui.add_space(8.0);
        let ok = self.busy.is_none() && self.doc.mesh.is_some();
        if ui.add_enabled(ok, egui::Button::new(RichText::new("  Compute  ").strong()).fill(Color32::from_rgb(38, 110, 70))).clicked() {
            self.compute();
        }
    }

    // ------------------------------------------------------------------ results

    fn tab_results(&mut self, ui: &mut Ui) {
        let Some(mut cur) = self.current else {
            ui.weak("No results yet. Press Compute.");
            return;
        };
        ui.heading("Results");
        if self.results.len() > 1 {
            egui::ComboBox::from_id_salt("result_sel").selected_text(self.results[cur].name.clone()).show_ui(ui, |ui| {
                for (i, r) in self.results.iter().enumerate() {
                    ui.selectable_value(&mut cur, i, r.name.clone());
                }
            });
            if cur != self.current.unwrap() {
                self.current = Some(cur);
                self.show_results = true;
            }
        } else {
            ui.label(self.results[cur].name.clone());
        }
        let mut cmap = self.proj.colormap;
        let mut fit = false;
        {
            let rv = &mut self.results[cur];
            let nsol = rv.result.solutions.len();
            if nsol > 1 {
                ui.horizontal(|ui| {
                    if ui.small_button("<").clicked() && rv.sol > 0 {
                        rv.sol -= 1;
                    }
                    let mut s = rv.sol;
                    ui.add(egui::Slider::new(&mut s, 0..=nsol - 1).show_value(false));
                    rv.sol = s;
                    if ui.small_button(">").clicked() && rv.sol + 1 < nsol {
                        rv.sol += 1;
                    }
                    ui.label(rv.result.solutions[rv.sol].label.clone());
                });
            } else {
                ui.label(rv.result.solutions[0].label.clone());
            }
            let sol = rv.sol.min(nsol - 1);
            let names: Vec<String> = rv.result.solutions[sol].fields.iter().map(|f| format!("{} [{}]", f.name, f.unit)).collect();
            if rv.field >= names.len() {
                rv.field = 0;
            }
            ui.horizontal(|ui| {
                ui.label("Field");
                egui::ComboBox::from_id_salt("field_sel").selected_text(names[rv.field].clone()).show_ui(ui, |ui| {
                    for (i, n) in names.iter().enumerate() {
                        ui.selectable_value(&mut rv.field, i, n);
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.label("Colours");
                egui::ComboBox::from_id_salt("cmap").selected_text(cmap.name()).show_ui(ui, |ui| {
                    for c in Colormap::ALL {
                        ui.selectable_value(&mut cmap, c, c.name());
                    }
                });
            });
            // range
            let auto = rv.manual_range.is_none();
            let mut auto_sel = auto;
            ui.horizontal(|ui| {
                ui.checkbox(&mut auto_sel, "auto range");
                if !auto_sel {
                    let (lo, hi) = rv.manual_range.unwrap_or_else(|| {
                        rv.manual_range = None;
                        rv.range()
                    });
                    let (mut a, mut b) = (lo as f64, hi as f64);
                    let c1 = num(ui, &mut a).changed();
                    ui.label("to");
                    let c2 = num(ui, &mut b).changed();
                    if c1 || c2 || rv.manual_range.is_none() {
                        rv.manual_range = Some((a as f32, b as f32));
                    }
                }
            });
            if auto_sel {
                rv.manual_range = None;
            } else if rv.manual_range.is_none() {
                let r = rv.range();
                rv.manual_range = Some(r);
            }
            if rv.result.solutions.iter().any(|s| s.disp.is_some()) {
                ui.horizontal(|ui| {
                    ui.label("deformation scale");
                    let max = (rv.auto_warp * 10.0).max(1.0);
                    ui.add(egui::Slider::new(&mut rv.warp, 0.0..=max).logarithmic(false).custom_formatter(|v, _| format!("{:.3e}", v)));
                    if ui.small_button("auto").clicked() {
                        rv.warp = rv.auto_warp;
                    }
                    if ui.small_button("0").clicked() {
                        rv.warp = 0.0;
                    }
                });
            }
            ui.checkbox(&mut self.show_edges, "mesh edges");
            if ui.button("Fit view").clicked() {
                fit = true;
            }
        }
        if fit {
            self.vp.fit_token += 1;
        }
        self.proj.colormap = cmap;
        ui.separator();
        ui.label(RichText::new("Summary").strong());
        let rv = &self.results[cur];
        egui::Grid::new("scalars").num_columns(2).striped(true).show(ui, |ui| {
            for (n, v) in &rv.result.scalars {
                ui.label(n);
                match v {
                    Scalar::Num(x) => ui.monospace(format!("{x:.6e}")),
                    Scalar::Text(t) => ui.monospace(t.clone()),
                };
                ui.end_row();
            }
        });
        ui.weak(format!("{} unknowns, solved in {:.2} s", rv.result.dofs, rv.result.seconds));
        egui::CollapsingHeader::new(format!("What the solver did ({} lines)", rv.result.log.len())).default_open(true).show(ui, |ui| {
            egui::ScrollArea::vertical().max_height(260.0).id_salt("solver_log").show(ui, |ui| {
                for l in &rv.result.log {
                    ui.add(egui::Label::new(RichText::new(l).monospace().weak().size(11.0)).wrap());
                }
            });
        });
        ui.add_space(6.0);
        if ui.button("Export CSV...").clicked() {
            self.export_csv();
        }
    }

    // ------------------------------------------------------------------ analysis

    /// Whether the result on show carries a polarization analysis.
    pub fn has_analysis(&self) -> bool {
        self.current.is_some_and(|i| self.results[i].result.analysis.is_some())
    }

    /// The analysis section: Stokes parameters of every optical mode, the Poincare sphere and the table of `femcore::analysis`.
    fn tab_analysis(&mut self, ui: &mut Ui) {
        use femcore::polarization::{mode_markers, sphere_drawing, SphereState, MODE_COLOURS};
        let Some(cur) = self.current else {
            ui.weak("No results yet. Press Compute.");
            return;
        };
        ui.heading("Analysis");
        let mut export: Option<&'static str> = None;
        {
            let rv = &mut self.results[cur];
            let selectable = rv.physics == Some(Physics::Emw);
            let Some(a) = rv.result.analysis.as_ref() else {
                ui.weak("This result has no polarization analysis. Mode studies of the electromagnetic waves physics and the Brillouin study (optical modes) have one.");
                return;
            };
            ui.label(RichText::new(if rv.physics == Some(Physics::Brillouin) { "Polarization of the optical modes (pump and Stokes)" } else { "Polarization of the modes" }).strong());
            ui.weak("Stokes parameters from the complex field of each mode, summed over the cross-section. s1 > 0 is x (lateral, quasi-TE), s1 < 0 is z (normal, quasi-TM), s2 the diagonals, s3 the sense of rotation (+ counter-clockwise). P is the degree of polarization: below 1 where the polarization changes across the mode.");
            // the picture
            let side = ui.available_width().clamp(220.0, 360.0);
            let (resp, painter) = ui.allocate_painter(egui::vec2(side, side), egui::Sense::hover());
            let dark = ui.visuals().dark_mode;
            let ink = if dark { Color32::from_rgb(155, 183, 217) } else { Color32::from_rgb(59, 91, 130) };
            let with_alpha = |c: Color32, o: f64| Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), (o.clamp(0.0, 1.0) * 255.0) as u8);
            let hex = |s: &str| {
                let v = u32::from_str_radix(s.trim_start_matches('#'), 16).unwrap_or(0);
                Color32::from_rgb((v >> 16) as u8, (v >> 8) as u8, v as u8)
            };
            let d = sphere_drawing(side as f64);
            let origin = resp.rect.min;
            let pt = |x: f64, y: f64| origin + egui::vec2(x as f32, y as f32);
            let centre = pt(d.centre.0, d.centre.1);
            let states: Vec<(usize, SphereState)> = a.rows.iter().enumerate().filter_map(|(i, r)| r.report.as_ref().map(|m| (i, SphereState { label: (i + 1).to_string(), s: m.sphere, dashed: r.group.is_some() }))).collect();
            let markers = mode_markers(&states.iter().map(|(_, s)| s.clone()).collect::<Vec<_>>(), side as f64);
            let selected = if selectable { Some(rv.sol) } else { None };
            let draw_markers = |front: bool| {
                for (m, (row, _)) in markers.iter().zip(&states) {
                    if m.front != front {
                        continue;
                    }
                    let col = with_alpha(hex(MODE_COLOURS[m.colour]), if front { 1.0 } else { 0.6 });
                    let (tip, c0) = (pt(m.tip.0, m.tip.1), pt(m.centre.0, m.centre.1));
                    if (tip - c0).length() >= 1.0 {
                        let stroke = egui::Stroke::new(1.8, col);
                        if m.dashed {
                            painter.extend(egui::Shape::dashed_line(&[c0, tip], stroke, 4.0, 3.0));
                        } else {
                            painter.line_segment([c0, tip], stroke);
                        }
                    }
                    let big = selected == Some(*row);
                    painter.circle(tip, side * if big { 0.040 } else { 0.030 }, col, egui::Stroke::new(if big { 2.0 } else { 0.8 }, Color32::WHITE));
                    painter.text(tip, egui::Align2::CENTER_CENTER, &m.label, egui::FontId::proportional(side * 0.036), Color32::WHITE);
                }
            };
            let line = |l: &femcore::polarization::WireLine| -> Vec<egui::Shape> {
                let pts: Vec<egui::Pos2> = l.pts.iter().map(|p| pt(p.0, p.1)).collect();
                if l.front {
                    vec![egui::Shape::line(pts, egui::Stroke::new(l.width as f32, with_alpha(ink, l.opacity)))]
                } else {
                    egui::Shape::dashed_line(&pts, egui::Stroke::new(l.width as f32, with_alpha(ink, l.opacity * 0.6)), 2.0, 3.0)
                }
            };
            // far side, the translucent ball, near side, labels
            for l in d.lines.iter().filter(|l| !l.front) {
                painter.extend(line(l));
            }
            draw_markers(false);
            painter.circle(centre, d.radius as f32, with_alpha(Color32::from_rgb(141, 179, 221), if dark { 0.10 } else { 0.22 }), egui::Stroke::new(1.0, with_alpha(ink, 0.6)));
            for l in d.lines.iter().filter(|l| l.front) {
                painter.extend(line(l));
            }
            draw_markers(true);
            for lab in &d.labels {
                painter.text(pt(lab.x, lab.y), egui::Align2::CENTER_CENTER, &lab.text, egui::FontId::proportional(side * 0.05), with_alpha(ui.visuals().text_color(), if lab.front { 0.9 } else { 0.5 }));
            }
            // the modes
            ui.add_space(4.0);
            let mut pick = None;
            for (i, r) in a.rows.iter().enumerate() {
                let colour = states.iter().position(|(k, _)| *k == i).map(|p| hex(MODE_COLOURS[p % MODE_COLOURS.len()])).unwrap_or(Color32::GRAY);
                let text = match &r.report {
                    Some(m) => format!("{:>2}  {} {:.5}  {}  P {:.2}  TE {:.2}{}", i + 1, a.value_label, r.neff.0, m.class.label(), m.dop, m.te_fraction, r.group.map_or(String::new(), |g| format!("  G{}", g + 1))),
                    None => format!("{:>2}  {} {:.5}  no intensity", i + 1, a.value_label, r.neff.0),
                };
                ui.horizontal(|ui| {
                    ui.label(RichText::new("\u{25cf}").color(colour));
                    let on = selected == Some(i);
                    let resp = ui.selectable_label(on, RichText::new(text).monospace().size(11.5));
                    if selectable && resp.clicked() {
                        pick = Some(i);
                    }
                });
            }
            if let Some(i) = pick {
                if i < rv.result.solutions.len() {
                    rv.sol = i;
                    rv.values_cache = None;
                }
            }
            for (gi, g) in a.groups.iter().enumerate() {
                let names: Vec<String> = g.members.iter().map(|i| (i + 1).to_string()).collect();
                egui::CollapsingHeader::new(format!("G{}: modes {} are degenerate", gi + 1, names.join(", "))).default_open(true).show(ui, |ui| {
                    ui.weak("The eigensolver returns some mixture of degenerate modes, so their own points on the sphere are not unique. The group state and the resolved pair are.");
                    if let Some(m) = &g.invariant {
                        ui.label(RichText::new(format!("group: {}", m.summary_line())).monospace().size(11.0));
                    }
                    if let Some(r) = &g.resolved {
                        for (k, rep) in r.reports.iter().enumerate() {
                            if let Some(m) = rep {
                                ui.label(RichText::new(format!("resolved {}: {}", k + 1, m.summary_line())).monospace().size(11.0));
                            }
                        }
                    }
                    if let Some(n) = &g.note {
                        ui.weak(n.clone());
                    }
                });
            }
            egui::CollapsingHeader::new("Table").default_open(false).show(ui, |ui| {
                egui::ScrollArea::horizontal().id_salt("pol_table").show(ui, |ui| {
                    ui.add(egui::Label::new(RichText::new(a.table()).monospace().size(10.5)).extend());
                });
            });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button("Export CSV...").clicked() {
                    export = Some("csv");
                }
                if ui.button("Export sphere (SVG)...").clicked() {
                    export = Some("svg");
                }
                if ui.button("Copy table").clicked() {
                    ui.ctx().copy_text(a.table());
                }
            });
        }
        if let Some(kind) = export {
            self.export_analysis(kind);
        }
    }

    /// Bottom panel: how the shown field's extremes (or the mode frequencies) vary over the
    /// solutions of the current result.
    pub fn results_panel(&mut self, ui: &mut Ui) {
        let Some(cur) = self.current else { return };
        let rv = &self.results[cur];
        let nsol = rv.result.solutions.len();
        ui.horizontal(|ui| {
            ui.label(RichText::new(rv.name.clone()).strong());
            ui.weak(format!("{} solution(s)", nsol));
        });
        if nsol < 2 {
            ui.weak("Single solution: use the viewport probe (hover) for values.");
            return;
        }
        let fname = rv.result.solutions[rv.sol.min(nsol - 1)].fields.get(rv.field).map(|f| f.name.clone()).unwrap_or_default();
        let xs: Vec<f64> = rv.result.solutions.iter().enumerate().map(|(i, s)| s.value.unwrap_or(i as f64 + 1.0)).collect();
        let is_modes = rv.result.solutions[0].neff.is_some();
        // eigenfrequencies: solids, Brillouin, and the 3D electromagnetic eigenmodes (frequency, no effective index)
        let is_eigen_plot = (rv.physics == Some(Physics::Solid) || rv.physics == Some(Physics::Brillouin) || (rv.physics == Some(Physics::Emw) && rv.result.solutions[0].neff.is_none())) && rv.result.solutions[0].value.is_some();
        let xlabel = if is_eigen_plot || is_modes { "mode number" } else if rv.result.solutions[0].value.is_some() { "time (s)" } else { "solution" };
        let mut max_pts = vec![];
        let mut min_pts = vec![];
        for (i, s) in rv.result.solutions.iter().enumerate() {
            if let Some(f) = s.fields.iter().find(|f| f.name == fname) {
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for &v in &f.values {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
                max_pts.push([xs[i], hi]);
                min_pts.push([xs[i], lo]);
            }
        }
        let current = xs[rv.sol.min(nsol - 1)];
        let is_eigen = is_eigen_plot;
        egui_plot::Plot::new("summary_plot").height(ui.available_height().max(90.0)).x_axis_label(xlabel).y_axis_label(if is_eigen { "frequency (Hz)" } else if is_modes { "|neff|" } else { fname.as_str() }).legend(egui_plot::Legend::default()).show(ui, |p| {
            if is_eigen {
                p.points(egui_plot::Points::new("f", egui_plot::PlotPoints::from(xs.iter().enumerate().map(|(i, _)| [i as f64 + 1.0, xs[i]]).collect::<Vec<_>>())).radius(4.0));
            } else if is_modes {
                // neff per mode
                let pts: Vec<[f64; 2]> = rv.result.solutions.iter().enumerate().map(|(i, s)| [i as f64 + 1.0, s.neff.map(|n| n.0.hypot(n.1)).unwrap_or(0.0)]).collect();
                p.points(egui_plot::Points::new("|neff|", egui_plot::PlotPoints::from(pts)).radius(4.0));
            } else {
                p.line(egui_plot::Line::new(format!("max {fname}"), egui_plot::PlotPoints::from(max_pts.clone())));
                p.line(egui_plot::Line::new(format!("min {fname}"), egui_plot::PlotPoints::from(min_pts.clone())));
                p.vline(egui_plot::VLine::new("now", current).color(Color32::from_rgb(230, 190, 80)));
            }
        });
    }
}
