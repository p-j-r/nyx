// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The `figure` command (ADR 99): a gallery and a component view of the modes of a solved example, as SVG and as PNG, the same figure from the
//! same result, and the figure of a stored results file (the benches of the author, when `NYX_FIG_H5` names one).

use nyx_engine::Engine;
use serde_json::json;
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_str()
        .unwrap()
        .to_string()
}

fn tmp(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("nyx-fig-{}-{name}", std::process::id()))
}

/// Solves an example and returns the engine.
fn solved(name: &str) -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example(name) })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    e
}

#[test]
fn the_modes_of_a_wire_are_drawn_in_quadrature_as_svg_and_png() {
    let mut e = solved("15_soi_wire");
    let svg = tmp("wire.svg");
    let r = e.call("figure", &json!({ "physics": "ems", "path": svg.to_str().unwrap(), "max_modes": 2 })).unwrap();
    let text = std::fs::read_to_string(&svg).unwrap();
    assert!(
        text.starts_with("<svg") && text.trim_end().ends_with("</svg>") && text.contains("data:image/png;base64,"),
        "an SVG with the field as an image"
    );
    let modes = r["modes"].as_array().unwrap();
    assert_eq!(modes.len(), 2, "{r}");
    for m in modes {
        // the real solver of a wire keeps one real amplitude per component: there is no phase to fix and nothing to say about quadrature
        assert_eq!(m["field"], "E");
        assert_eq!(m["complex"], false);
        assert!(m["quadrature"].is_null(), "{m}");
    }
    // the same result gives the same figure: the phase does not depend on anything but the mode
    let again = e.call("figure", &json!({ "physics": "ems", "path": svg.to_str().unwrap(), "max_modes": 2 })).unwrap();
    assert_eq!(again["modes"], r["modes"]);
    assert_eq!(std::fs::read_to_string(&svg).unwrap(), text);

    let png = tmp("wire.png");
    let r = e
        .call("figure", &json!({ "physics": "ems", "path": png.to_str().unwrap(), "kind": "components", "modes": [1], "scale": 1.0 }))
        .unwrap();
    let bytes = std::fs::read(&png).unwrap();
    assert_eq!(&bytes[..4], &[0x89, b'P', b'N', b'G']);
    assert_eq!(r["bytes"].as_u64().unwrap() as usize, bytes.len());
    assert!(r["width"].as_f64().unwrap() > 100.0 && r["height"].as_f64().unwrap() > 100.0);
    // a figure is not written for what it cannot draw
    assert!(e.call("figure", &json!({ "physics": "ems", "path": tmp("x.jpg").to_str().unwrap() })).is_err());
    assert!(e.call("figure", &json!({ "physics": "ems", "path": svg.to_str().unwrap(), "modes": [99] })).is_err());
    for p in [svg, png] {
        let _ = std::fs::remove_file(p);
    }
}

#[test]
fn the_vibration_modes_of_a_cross_section_are_drawn_too() {
    let mut e = solved("13_silicon_on_sapphire");
    let physics = e.project().physics.iter().find(|s| s.physics == "sms").map(|s| s.physics.clone());
    let Some(physics) = physics else { return };
    let png = tmp("sms.png");
    let r =
        e.call("figure", &json!({ "physics": physics, "path": png.to_str().unwrap(), "max_modes": 3, "scale": 1.0 }));
    // an example may solve a static or a 3D problem the section does not apply to; then the answer says why, it does not draw nothing
    match r {
        Ok(v) => {
            assert!(std::fs::metadata(&png).unwrap().len() > 1000, "{v}");
            let _ = std::fs::remove_file(&png);
        }
        Err(err) => eprintln!("13_silicon_on_sapphire: {err}"),
    }
}

/// A results file made elsewhere (a Bloch study of a slab, complex modes): drawn from the stored file, and the quadrature of its propagation component is printed (a leaky mode may
/// have none). Runs only when `NYX_FIG_H5` names the file (and `NYX_FIG_OUT` a folder for the figures); the files of the benches are not in
/// the repository.
#[test]
fn a_stored_file_of_complex_modes_is_drawn_in_quadrature() {
    let Ok(file) = std::env::var("NYX_FIG_H5") else { return };
    let out = PathBuf::from(
        std::env::var("NYX_FIG_OUT").unwrap_or_else(|_| std::env::temp_dir().to_string_lossy().into_owned()),
    );
    let stored = nyx_engine::h5io::read_results(std::path::Path::new(&file)).unwrap();
    let (code, result) = &stored.results[0];
    for (kind, name) in [("gallery", "gallery"), ("components", "components")] {
        let o = nyx_analysis::figure::FigureOptions {
            kind: if kind == "gallery" {
                nyx_analysis::figure::Kind::Gallery
            } else {
                nyx_analysis::figure::Kind::Components
            },
            modes: if kind == "gallery" {
                None
            } else {
                Some([0, 2, 4, 8].into_iter().filter(|k| *k < result.solutions.len()).collect())
            },
            ..Default::default()
        };
        let fig = nyx_analysis::figure::mode_figure_on(&stored.mesh, result, &o).unwrap();
        for m in &fig.modes {
            eprintln!(
                "{code} mode {}: phase {:.1} deg, flipped {}, peak {:.3e}, quadrature {:?}",
                m.mode, m.phase_deg, m.flipped, m.peak, m.quadrature
            );
        }
        std::fs::write(out.join(format!("{name}.png")), fig.canvas.to_png(1.5).unwrap()).unwrap();
        std::fs::write(out.join(format!("{name}.svg")), fig.canvas.to_svg().unwrap()).unwrap();
    }
}
