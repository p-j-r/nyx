// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Colour-vision simulation for the screenshots of the scripted runs (UX plan U6): `cvd protan|deutan|tritan` makes every later `shot` show the window
//! as a person with that kind of colour-vision difference sees it, so the canvas, the legend and the status colours can be looked at that way
//! without such a person at hand. It is the Machado, Oliveira and Fernandes (2009) model at full severity, applied in linear sRGB, the same
//! matrices the palette test of `tokens.rs` uses. A simulation is an approximation of one kind of difference, not a substitute for a person who has it.

/// Which cone type is missing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Protan,
    Deutan,
    Tritan,
}

impl Kind {
    pub fn parse(s: &str) -> Option<Kind> {
        match s.trim() {
            "protan" => Some(Kind::Protan),
            "deutan" => Some(Kind::Deutan),
            "tritan" => Some(Kind::Tritan),
            _ => None,
        }
    }

    fn matrix(self) -> [[f32; 3]; 3] {
        match self {
            Kind::Protan => [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
            Kind::Deutan => [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
            Kind::Tritan => [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
        }
    }
}

fn to_linear(c: u8) -> f32 {
    let s = c as f32 / 255.0;
    if s <= 0.04045 {
        s / 12.92
    } else {
        ((s + 0.055) / 1.055).powf(2.4)
    }
}

fn to_srgb(l: f32) -> u8 {
    let l = l.clamp(0.0, 1.0);
    let s = if l <= 0.0031308 { l * 12.92 } else { 1.055 * l.powf(1.0 / 2.4) - 0.055 };
    (s * 255.0).round().clamp(0.0, 255.0) as u8
}

/// One colour as the simulated viewer sees it.
pub fn simulate(rgb: [u8; 3], kind: Kind) -> [u8; 3] {
    let m = kind.matrix();
    let l = [to_linear(rgb[0]), to_linear(rgb[1]), to_linear(rgb[2])];
    let row = |r: usize| to_srgb(m[r][0] * l[0] + m[r][1] * l[1] + m[r][2] * l[2]);
    [row(0), row(1), row(2)]
}

/// A whole image, `rgba` bytes in place (alpha is kept).
pub fn simulate_rgba(raw: &mut [u8], kind: Kind) {
    for px in raw.as_chunks_mut::<4>().0.iter_mut() {
        let s = simulate([px[0], px[1], px[2]], kind);
        px[..3].copy_from_slice(&s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn greys_stay_grey_and_the_confusion_colours_come_together() {
        for kind in [Kind::Protan, Kind::Deutan, Kind::Tritan] {
            for g in [0u8, 90, 200, 255] {
                let s = simulate([g, g, g], kind);
                assert!(s.iter().all(|c| (*c as i32 - g as i32).abs() <= 2), "{kind:?} {g}: {s:?}");
            }
        }
        // a red and a green of one lightness are far apart for normal vision and close for a deuteranope
        let (red, green) = ([200u8, 50, 30], [30u8, 125, 60]);
        let dist = |a: [u8; 3], b: [u8; 3]| a.iter().zip(b).map(|(x, y)| (*x as f32 - y as f32).powi(2)).sum::<f32>().sqrt();
        assert!(dist(simulate(red, Kind::Deutan), simulate(green, Kind::Deutan)) < dist(red, green) * 0.6);
    }

    #[test]
    fn kinds_are_read_by_name() {
        assert_eq!(Kind::parse("deutan"), Some(Kind::Deutan));
        assert_eq!(Kind::parse("none"), None);
    }
}
