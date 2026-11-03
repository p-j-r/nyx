//! Units of material files: a quantity is written as `"496.8 GPa"`, as `{ value = 496.8, unit = "GPa" }`, or as a bare number
//! next to a declared unit. This module turns a unit name into a dimension and a factor to SI and refuses a unit of the wrong
//! kind with a message that says what was expected. All values leave this module in SI.

/// The kind of quantity an entry holds; a unit of another kind is refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dim {
    Dimensionless,
    Pressure,
    Density,
    Length,
    /// an absolute temperature (K or degC)
    Temperature,
    ThermalConductivity,
    SpecificHeat,
    /// thermal expansion and thermo-optic coefficients, 1/K
    PerKelvin,
    /// inverse pressure, m^2/N
    Compliance,
    ElectricalConductivity,
}

impl Dim {
    pub fn describe(self) -> &'static str {
        match self {
            Dim::Dimensionless => "a plain number (no unit)",
            Dim::Pressure => "a stiffness or pressure, for example GPa",
            Dim::Density => "a density, for example kg/m^3",
            Dim::Length => "a length, for example um or nm",
            Dim::Temperature => "a temperature, for example K or degC",
            Dim::ThermalConductivity => "a thermal conductivity, for example W/(m*K)",
            Dim::SpecificHeat => "a specific heat, for example J/(kg*K)",
            Dim::PerKelvin => "a coefficient per kelvin, for example 1/K or ppm/K",
            Dim::Compliance => "a compliance, for example 1/GPa",
            Dim::ElectricalConductivity => "an electrical conductivity, for example S/m",
        }
    }
}

fn normalise(u: &str) -> String {
    let mut s = String::new();
    for ch in u.trim().chars() {
        match ch {
            ' ' | '*' | '\u{b7}' => {}
            '\u{3bc}' | '\u{b5}' => s.push('u'),
            '\u{b0}' => s.push_str("deg"),
            '\u{c5}' => s.push_str("angstrom"),
            c => s.push(c),
        }
    }
    s
}

fn prefix(p: &str) -> Option<f64> {
    Some(match p {
        "p" => 1e-12,
        "n" => 1e-9,
        "u" => 1e-6,
        "m" => 1e-3,
        "c" => 1e-2,
        "" => 1.0,
        "k" => 1e3,
        "M" => 1e6,
        "G" => 1e9,
        "T" => 1e12,
        _ => return None,
    })
}

/// Splits `<prefix><base>` when `u` ends with `base`.
fn with_prefix(u: &str, base: &str) -> Option<f64> {
    u.strip_suffix(base).and_then(prefix)
}

/// Dimension and factor to SI of a unit name. `None` if the unit is unknown.
pub fn lookup(unit: &str) -> Option<(Dim, f64)> {
    let u = normalise(unit);
    match u.as_str() {
        "" | "1" => return Some((Dim::Dimensionless, 1.0)),
        "kg/m^3" | "kg/m3" => return Some((Dim::Density, 1.0)),
        "g/cm^3" | "g/cm3" | "g/cc" => return Some((Dim::Density, 1e3)),
        "W/(mK)" | "W/mK" | "W/m/K" => return Some((Dim::ThermalConductivity, 1.0)),
        "J/(kgK)" | "J/kgK" | "J/kg/K" => return Some((Dim::SpecificHeat, 1.0)),
        "kJ/(kgK)" | "kJ/kgK" | "kJ/kg/K" => return Some((Dim::SpecificHeat, 1e3)),
        "1/K" | "K^-1" | "/K" | "1/degC" => return Some((Dim::PerKelvin, 1.0)),
        "ppm/K" | "1e-6/K" | "10^-6/K" | "e-6/K" | "ppm/degC" | "1e-6/degC" => return Some((Dim::PerKelvin, 1e-6)),
        "1e-5/K" | "10^-5/K" => return Some((Dim::PerKelvin, 1e-5)),
        "m^2/N" | "m2/N" => return Some((Dim::Compliance, 1.0)),
        "K" => return Some((Dim::Temperature, 1.0)),
        "degC" | "C" => return Some((Dim::Temperature, 1.0)),
        "angstrom" => return Some((Dim::Length, 1e-10)),
        "S/m" => return Some((Dim::ElectricalConductivity, 1.0)),
        _ => {}
    }
    if let Some(f) = with_prefix(&u, "Pa") {
        return Some((Dim::Pressure, f));
    }
    if let Some(rest) = u.strip_prefix("1/") {
        if let Some(f) = with_prefix(rest, "Pa") {
            return Some((Dim::Compliance, 1.0 / f));
        }
    }
    if let Some(head) = u.strip_suffix("Pa^-1") {
        if let Some(f) = prefix(head) {
            return Some((Dim::Compliance, 1.0 / f));
        }
    }
    if let Some(head) = u.strip_suffix("S/m") {
        if let Some(f) = prefix(head) {
            return Some((Dim::ElectricalConductivity, f));
        }
    }
    if let Some(f) = with_prefix(&u, "m") {
        return Some((Dim::Length, f));
    }
    None
}

/// `value` in `unit` converted to SI, if `unit` is a unit of dimension `expect`.
pub fn to_si(value: f64, unit: &str, expect: Dim) -> Result<f64, String> {
    let (dim, f) = lookup(unit).ok_or_else(|| format!("the unit '{unit}' is not known; expected {}", expect.describe()))?;
    if dim != expect {
        return Err(format!("the unit '{unit}' does not fit here; expected {}", expect.describe()));
    }
    if expect == Dim::Temperature {
        let u = normalise(unit);
        return Ok(if u == "K" { value } else { value + 273.15 });
    }
    Ok(value * f)
}

/// Splits `"496.8 GPa"` (or `"496.8GPa"`) into the number and the unit text.
pub fn split_quantity(s: &str) -> Result<(f64, String), String> {
    let t = s.trim();
    if let Some((num, rest)) = t.split_once(char::is_whitespace) {
        if let Ok(v) = num.parse::<f64>() {
            return Ok((v, rest.trim().to_string()));
        }
    }
    for end in (1..=t.len()).rev().filter(|&i| t.is_char_boundary(i)) {
        if let Ok(v) = t[..end].parse::<f64>() {
            if end < t.len() && !t[end..].starts_with(['e', 'E']) {
                return Ok((v, t[end..].trim().to_string()));
            }
        }
    }
    t.parse::<f64>().map(|v| (v, String::new())).map_err(|_| format!("'{s}' is not a number with a unit; write it like \"496.8 GPa\""))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_units() {
        assert_eq!(to_si(496.8, "GPa", Dim::Pressure).unwrap(), 496.8e9);
        assert_eq!(to_si(2.2, "g/cm^3", Dim::Density).unwrap(), 2200.0);
        assert_eq!(to_si(1.55, "um", Dim::Length).unwrap(), 1.55e-6);
        assert_eq!(to_si(1.55, "\u{b5}m", Dim::Length).unwrap(), 1.55e-6);
        assert_eq!(to_si(5.0, "ppm/K", Dim::PerKelvin).unwrap(), 5.0e-6);
        assert!((to_si(20.0, "degC", Dim::Temperature).unwrap() - 293.15).abs() < 1e-12);
        assert_eq!(to_si(1.0, "1/GPa", Dim::Compliance).unwrap(), 1e-9);
        assert_eq!(to_si(3.0, "W/(m*K)", Dim::ThermalConductivity).unwrap(), 3.0);
    }

    #[test]
    fn wrong_kind_is_refused_with_a_hint() {
        let e = to_si(1.0, "GPa", Dim::Density).unwrap_err();
        assert!(e.contains("density"), "{e}");
        assert!(to_si(1.0, "furlong", Dim::Length).unwrap_err().contains("not known"));
    }

    #[test]
    fn quantity_strings() {
        assert_eq!(split_quantity("496.8 GPa").unwrap(), (496.8, "GPa".into()));
        assert_eq!(split_quantity("496.8GPa").unwrap(), (496.8, "GPa".into()));
        assert_eq!(split_quantity("1e-6 1/K").unwrap(), (1e-6, "1/K".into()));
        assert_eq!(split_quantity("12.5").unwrap(), (12.5, String::new()));
        assert!(split_quantity("GPa").is_err());
    }
}
