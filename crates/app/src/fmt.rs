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

/// The unit a metric or parameter key ends with ("max_u_m" is in metres), when the last word is a unit.
pub fn unit_of(key: &str) -> Option<&'static str> {
    const UNITS: [(&str, &str); 12] = [("m", "m"), ("pa", "Pa"), ("hz", "Hz"), ("s", "s"), ("k", "K"), ("w", "W"), ("v", "V"), ("n", "N"), ("j", "J"), ("a", "A"), ("kg", "kg"), ("rad", "rad")];
    let (_, last) = key.rsplit_once('_')?;
    UNITS.iter().find(|(w, _)| *w == last).map(|(_, u)| *u)
}

/// "max_von_mises_pa" as "Max von Mises [Pa]": first word capitalised, a few names spelled right, the unit in brackets.
pub fn key_label(key: &str) -> String {
    let unit = unit_of(key);
    let stem = if unit.is_some() { key.rsplit_once('_').map(|(s, _)| s).unwrap_or(key) } else { key };
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

/// A number as the runs table shows it: SI prefixes when the key names a unit, whole numbers without decimals, lengths of a parameter in nm,
/// small and large values in scientific notation.
pub fn value(key: &str, v: f64) -> String {
    if let Some(u) = unit_of(key) {
        if v.is_finite() {
            return si(v, u);
        }
    }
    let lengthy = v != 0.0 && v.abs() < 1e-3 && ["width", "height", "thick", "diam", "wavelength", "layer", "box", "gap"].iter().any(|w| key.contains(w)) || (v != 0.0 && v.abs() < 1e-3 && key.starts_with('l'));
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
    }

    #[test]
    fn values_show_units_whole_numbers_and_notation() {
        assert_eq!(value("max_u_m", 2.017e-4), "201.7 µm");
        assert_eq!(value("dofs", 3214.0), "3214");
        assert_eq!(value("width", 5e-7), "500 nm");
        assert_eq!(value("neff_re", 2.4), "2.4000");
    }
}
