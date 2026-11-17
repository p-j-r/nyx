//! Colour maps for field plots (Viridis, Cividis, Turbo, Inferno, cool-warm, gray) as functions of a value in [0, 1], and the
//! palette used to tell boundary conditions apart (the eight region colours of `tokens.rs`).

use crate::tokens::{self, Mode};
use eframe::egui::Color32;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Colormap {
    Viridis,
    /// designed to read the same for every kind of colour-vision difference (Nunez, Anderton and Renslow, 2018)
    Cividis,
    Turbo,
    Inferno,
    Coolwarm,
    Gray,
}

impl Colormap {
    pub const ALL: [Colormap; 6] = [Colormap::Viridis, Colormap::Cividis, Colormap::Turbo, Colormap::Inferno, Colormap::Coolwarm, Colormap::Gray];

    /// Turbo is a rainbow: its lightness goes up and down, so equal steps of the value do not look like equal steps. It says so in its name.
    pub fn name(self) -> &'static str {
        match self {
            Colormap::Viridis => "Viridis",
            Colormap::Cividis => "Cividis (colour-vision safe)",
            Colormap::Turbo => "Turbo (rainbow)",
            Colormap::Inferno => "Inferno",
            Colormap::Coolwarm => "Cool-warm",
            Colormap::Gray => "Gray",
        }
    }

    /// Colour at `t` in [0, 1]. (The Turbo polynomial coefficients are the published ones, with all their digits.)
    #[allow(clippy::excessive_precision)]
    pub fn sample(self, t: f32) -> Color32 {
        let t = if t.is_finite() { t.clamp(0.0, 1.0) } else { 0.0 };
        let rgb = match self {
            Colormap::Viridis => lerp_table(&VIRIDIS, t),
            Colormap::Inferno => lerp_table(&INFERNO, t),
            Colormap::Cividis => lerp_table(&CIVIDIS, t),
            Colormap::Coolwarm => lerp_table(&COOLWARM, t),
            Colormap::Gray => [t, t, t],
            Colormap::Turbo => {
                let r = 0.13572138 + t * (4.6153926 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))));
                let g = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))));
                let b = 0.1066733 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))));
                [r, g, b]
            }
        };
        Color32::from_rgb(
            (rgb[0].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
            (rgb[1].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
            (rgb[2].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
        )
    }
}

fn lerp_table(table: &[[f32; 3]], t: f32) -> [f32; 3] {
    let x = t * (table.len() - 1) as f32;
    let i = (x.floor() as usize).min(table.len() - 2);
    let f = x - i as f32;
    let (a, b) = (table[i], table[i + 1]);
    [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]
}

// sampled colour values (0.318 is a green component, not 1/pi)
#[allow(clippy::approx_constant)]
const VIRIDIS: [[f32; 3]; 11] = [
    [0.267, 0.005, 0.329],
    [0.283, 0.141, 0.458],
    [0.254, 0.265, 0.530],
    [0.207, 0.372, 0.553],
    [0.164, 0.471, 0.558],
    [0.128, 0.567, 0.551],
    [0.135, 0.659, 0.518],
    [0.267, 0.749, 0.441],
    [0.478, 0.821, 0.318],
    [0.741, 0.873, 0.150],
    [0.993, 0.906, 0.144],
];

const INFERNO: [[f32; 3]; 11] = [
    [0.001, 0.000, 0.014],
    [0.087, 0.044, 0.224],
    [0.258, 0.039, 0.406],
    [0.416, 0.090, 0.433],
    [0.578, 0.148, 0.404],
    [0.735, 0.216, 0.330],
    [0.865, 0.317, 0.226],
    [0.954, 0.469, 0.103],
    [0.988, 0.652, 0.040],
    [0.964, 0.859, 0.264],
    [0.988, 0.998, 0.645],
];

const COOLWARM: [[f32; 3]; 5] = [
    [0.230, 0.299, 0.754],
    [0.552, 0.691, 0.996],
    [0.866, 0.866, 0.866],
    [0.956, 0.604, 0.486],
    [0.706, 0.016, 0.150],
];

/// Eleven samples (every 25.5th of 256) of matplotlib's `cividis` table (`_cm_listed.py`, CC0): values read from that file, not typed from memory.
const CIVIDIS: [[f32; 3]; 11] = [
    [0.000000, 0.135112, 0.304751],
    [0.046205, 0.203903, 0.440196],
    [0.208926, 0.272546, 0.424809],
    [0.313287, 0.343120, 0.425120],
    [0.401418, 0.411790, 0.440708],
    [0.488697, 0.485318, 0.471008],
    [0.582087, 0.558670, 0.468118],
    [0.683950, 0.638793, 0.444251],
    [0.785965, 0.720438, 0.399613],
    [0.896818, 0.811030, 0.320982],
    [0.995737, 0.909344, 0.217772],
];

/// One of the eight validated categorical colours for the canvas of `mode` (any two are told apart under protan and deutan vision, see
/// `tokens.rs`). Index `i` wraps at eight: the policy for more than eight things (hatch, grey) is in `regions.rs`.
pub fn palette(i: usize, mode: Mode) -> Color32 {
    crate::theme::rgb(tokens::regions(mode)[i % 8])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn luma(c: Color32) -> f32 {
        0.2126 * f32::from(c.r()) + 0.7152 * f32::from(c.g()) + 0.0722 * f32::from(c.b())
    }

    #[test]
    fn cividis_gets_lighter_with_the_value_and_viridis_too() {
        for m in [Colormap::Cividis, Colormap::Viridis, Colormap::Gray] {
            let l: Vec<f32> = (0..=40).map(|i| luma(m.sample(i as f32 / 40.0))).collect();
            assert!(l.windows(2).all(|w| w[1] >= w[0] - 0.5), "{}: lightness is not monotone", m.name());
        }
        // the rainbow is not: its lightness goes up and down, and its name says so
        let t: Vec<f32> = (0..=40).map(|i| luma(Colormap::Turbo.sample(i as f32 / 40.0))).collect();
        assert!(!t.windows(2).all(|w| w[1] >= w[0] - 0.5));
        assert!(Colormap::Turbo.name().contains("rainbow"));
    }

    #[test]
    fn the_categorical_palette_is_the_validated_eight() {
        for mode in [Mode::Dark, Mode::Light] {
            for i in 0..8 {
                assert_eq!(palette(i, mode), crate::theme::rgb(tokens::regions(mode)[i]));
            }
            assert_eq!(palette(8, mode), palette(0, mode));
        }
    }
}
