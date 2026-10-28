//! Colour maps for field plots (Viridis, Turbo, Inferno, cool-warm, gray) as functions of a value in [0, 1], and the
//! palette used to tell domains and boundaries apart.

use eframe::egui::Color32;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Colormap {
    Viridis,
    Turbo,
    Inferno,
    Coolwarm,
    Gray,
}

impl Colormap {
    pub const ALL: [Colormap; 5] = [Colormap::Viridis, Colormap::Turbo, Colormap::Inferno, Colormap::Coolwarm, Colormap::Gray];

    pub fn name(self) -> &'static str {
        match self {
            Colormap::Viridis => "Viridis",
            Colormap::Turbo => "Turbo",
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

/// A categorical palette for domain and boundary numbers.
pub fn palette(i: usize) -> Color32 {
    const P: [[u8; 3]; 12] = [
        [86, 154, 211],
        [230, 145, 56],
        [106, 190, 108],
        [199, 99, 99],
        [164, 126, 214],
        [210, 180, 70],
        [80, 190, 190],
        [214, 120, 170],
        [150, 150, 150],
        [120, 170, 90],
        [190, 130, 90],
        [100, 120, 200],
    ];
    let c = P[i % P.len()];
    Color32::from_rgb(c[0], c[1], c[2])
}
