// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Small text helpers of the window (plan U3): plurals ("1 run", "2 runs"), ellipsis for long titles, SI prefixes ("201.7 µm"), and readable labels
//! for the metric and parameter keys of the runs table ("max_von_mises_pa" becomes "Max von Mises [Pa]"). Numbers are stored in SI everywhere; this is
//! only how they are shown.

/// "1 run", "2 runs", "0 runs": the noun gets an "s" unless the count is one.
pub fn plural(n: usize, noun: &str) -> String {
    format!("{n} {}", if n == 1 { noun.to_string() } else { format!("{noun}s") })
}

/// The text cut to `max` characters with an ellipsis; text that fits is returned as it is.
pub fn ellipsize(s: &str, max: usize) -> String {
    if s.chars().count() <= max {
        return s.to_string();
    }
    s.chars().take(max.saturating_sub(1)).collect::<String>().trim_end().to_string() + "…"
}

/// A time that has passed: "0.4 s", "12 s", "1:05", "1:02:09" (seconds below ten show a decimal; a minute and more show minutes and seconds).
pub fn seconds(s: f64) -> String {
    let s = if s.is_finite() { s.max(0.0) } else { 0.0 };
    if s < 10.0 {
        format!("{s:.1} s")
    } else if s < 60.0 {
        format!("{s:.0} s")
    } else {
        let t = s.round() as u64;
        if t < 3600 {
            format!("{}:{:02}", t / 60, t % 60)
        } else {
            format!("{}:{:02}:{:02}", t / 3600, (t / 60) % 60, t % 60)
        }
    }
}

/// A value with an SI prefix and four significant digits: `si(2.017e-4, "m")` is "201.7 µm", `si(6.175e6, "Pa")` is "6.175 MPa".
pub fn si(v: f64, unit: &str) -> String {
    if v == 0.0 || !v.is_finite() {
        return format!("{v} {unit}").replace("inf", "∞");
    }
    const PREFIXES: [&str; 17] = ["y", "z", "a", "f", "p", "n", "µ", "m", "", "k", "M", "G", "T", "P", "E", "Z", "Y"];
    let mut exp = ((v.abs().log10() / 3.0).floor() as i32).clamp(-8, 8);
    loop {
        let m = v / 10f64.powi(3 * exp);
        let digits = 3 - m.abs().log10().floor() as i32; // four significant digits
        let text = format!("{m:.*}", digits.max(0) as usize);
        let rounded: f64 = text.parse().unwrap_or(m);
        if rounded.abs() >= 1000.0 && exp < 8 {
            exp += 1; // 999.97 rounds to 1000: use the next prefix
            continue;
        }
        let text = if text.contains('.') { text.trim_end_matches('0').trim_end_matches('.').to_string() } else { text };
        return format!("{text} {}{unit}", PREFIXES[(exp + 8) as usize]);
    }
}

/// The unit a metric key ends with ("max_u_m" is in metres), when the last word is a physical unit.
/// Parameter keys ("parameter_w", "w", "r_w") never carry an implicit SI unit suffix.
pub fn unit_of(key: &str) -> Option<&'static str> {
    if key.starts_with("parameter_") {
        return None;
    }
    let (stem, last) = key.rsplit_once('_')?;
    match last {
        "pa" => return Some("Pa"),
        "hz" => return Some("Hz"),
        "kg" => return Some("kg"),
        "rad" => return Some("rad"),
        _ => {}
    }
    const SINGLE_UNITS: [(&str, &str); 8] = [("m", "m"), ("s", "s"), ("k", "K"), ("w", "W"), ("v", "V"), ("n", "N"), ("j", "J"), ("a", "A")];
    const METRIC_PREFIXES: [&str; 22] = [
        "max_", "min_", "mean_", "avg_", "total_", "rms_", "peak_", "net_", "freq_", "time_", "power_", "energy_", "force_", "temp_",
        "temperature_", "voltage_", "potential_", "current_", "disp_", "displacement_", "wavelength_", "charge_",
    ];
    const METRIC_STEMS: [&str; 20] = [
        "max", "min", "mean", "avg", "total", "rms", "peak", "net", "freq", "frequency", "time", "power", "energy", "force", "temp",
        "temperature", "voltage", "potential", "current", "displacement",
    ];
    let is_metric = METRIC_PREFIXES.iter().any(|p| stem.starts_with(p)) || METRIC_STEMS.contains(&stem);
    if !is_metric {
        return None;
    }
    SINGLE_UNITS.iter().find(|(w, _)| *w == last).map(|(_, u)| *u)
}

/// "max_von_mises_pa" as "Max von Mises [Pa]": first word capitalised, a few names spelled right, the unit in brackets.
pub fn key_label(key: &str) -> String {
    let clean = key.strip_prefix("parameter_").filter(|s| !s.is_empty()).unwrap_or(key);
    let unit = unit_of(key);
    let stem = if unit.is_some() { clean.rsplit_once('_').map(|(s, _)| s).unwrap_or(clean) } else { clean };
    let words: Vec<String> = stem
        .split('_')
        .filter(|w| !w.is_empty())
        .enumerate()
        .map(|(i, w)| {
            let special = match w {
                "mises" => Some("Mises"),
                "neff" => Some("n_eff"),
                "dofs" => Some("DOFs"),
                "fsr" => Some("FSR"),
                "pml" => Some("PML"),
                "te" => Some("TE"),
                "tm" => Some("TM"),
                "re" => Some("(real)"),
                "im" => Some("(imaginary)"),
                _ => None,
            };
            let word = special.unwrap_or(w);
            if i == 0 && !word.starts_with('(') {
                word[..1].to_uppercase() + &word[1..]
            } else {
                word.to_string()
            }
        })
        .collect();
    match unit {
        Some(u) => format!("{} [{u}]", words.join(" ")),
        None => words.join(" "),
    }
}

fn is_length_param(key: &str) -> bool {
    let k = key.strip_prefix("parameter_").filter(|s| !s.is_empty()).unwrap_or(key);
    ["width", "height", "thick", "diam", "wavelength", "layer", "box", "gap", "radius", "pitch", "period"].iter().any(|w| k.contains(w))
        || k.starts_with('l')
        || matches!(k, "w" | "t" | "h" | "r" | "d" | "g")
        || ["w_", "h_", "t_", "r_", "d_", "g_"].iter().any(|p| k.starts_with(p))
}

/// A number as the runs table shows it: SI prefixes when the key names a unit, whole numbers without decimals, lengths of a parameter in nm,
/// small and large values in scientific notation.
pub fn value(key: &str, v: f64) -> String {
    if let Some(u) = unit_of(key) {
        if v.is_finite() {
            return si(v, u);
        }
    }
    let lengthy = v != 0.0 && v.abs() < 1e-3 && is_length_param(key);
    if lengthy {
        let s = format!("{:.3}", v * 1e9);
        format!("{} nm", s.trim_end_matches('0').trim_end_matches('.'))
    } else if v.is_finite() && v.fract() == 0.0 && v.abs() < 1e9 {
        format!("{v:.0}")
    } else if v != 0.0 && (v.abs() < 1e-3 || v.abs() >= 1e6) {
        format!("{v:.3e}")
    } else {
        format!("{v:.4}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_time_that_has_passed_reads_naturally() {
        assert_eq!(seconds(0.43), "0.4 s");
        assert_eq!(seconds(12.4), "12 s");
        assert_eq!(seconds(65.0), "1:05");
        assert_eq!(seconds(3729.0), "1:02:09");
        assert_eq!(seconds(f64::NAN), "0.0 s");
    }

    #[test]
    fn plurals_follow_the_count() {
        assert_eq!(plural(0, "run"), "0 runs");
        assert_eq!(plural(1, "region"), "1 region");
        assert_eq!(plural(2, "experiment"), "2 experiments");
    }

    #[test]
    fn long_titles_get_an_ellipsis_and_short_ones_stay() {
        assert_eq!(ellipsize("short", 22), "short");
        let t = ellipsize("Parallel-plate capacitor with two dielectrics", 22);
        assert_eq!(t.chars().count(), 22);
        assert!(t.ends_with('…'));
    }

    #[test]
    fn si_prefixes_and_four_digits() {
        assert_eq!(si(2.017e-4, "m"), "201.7 µm");
        assert_eq!(si(6.175e6, "Pa"), "6.175 MPa");
        assert_eq!(si(193.4e12, "Hz"), "193.4 THz");
        assert_eq!(si(0.0, "m"), "0 m");
        assert_eq!(si(1.0, "s"), "1 s");
        assert_eq!(si(999.97, "m"), "1 km");
        assert_eq!(si(-1.5e-9, "m"), "-1.5 nm");
    }

    #[test]
    fn keys_read_as_labels_with_their_unit() {
        assert_eq!(key_label("max_u_m"), "Max u [m]");
        assert_eq!(key_label("max_von_mises_pa"), "Max von Mises [Pa]");
        assert_eq!(key_label("neff_re"), "N_eff (real)");
        assert_eq!(key_label("dofs"), "DOFs");
        assert_eq!(key_label("seconds"), "Seconds");
        assert_eq!(key_label("w"), "W");
        assert_eq!(key_label("parameter_w"), "W");
        assert_eq!(key_label("h_si"), "H si");
        assert_eq!(key_label("parameter_h_si"), "H si");
        assert_eq!(key_label("parameter_k"), "K");
        assert_eq!(unit_of("w"), None);
        assert_eq!(unit_of("parameter_w"), None);
        assert_eq!(unit_of("parameter_k"), None);
        assert_eq!(unit_of("r_w"), None);
    }

    #[test]
    fn values_show_units_whole_numbers_and_notation() {
        assert_eq!(value("max_u_m", 2.017e-4), "201.7 µm");
        assert_eq!(value("dofs", 3214.0), "3214");
        assert_eq!(value("width", 5e-7), "500 nm");
        assert_eq!(value("w", 5e-7), "500 nm");
        assert_eq!(value("parameter_w", 5e-7), "500 nm");
        assert_eq!(value("h_si", 3e-7), "300 nm");
        assert_eq!(value("parameter_h_si", 3e-7), "300 nm");
        assert_eq!(value("parameter_k", 3.5), "3.5000");
        assert_eq!(value("neff_re", 2.4), "2.4000");
    }
}
