// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
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
    /// angular frequency, frequency or photon energy; the SI value is rad/s
    Frequency,
    /// linear frequency; the SI value is Hz
    FrequencyHz,
    /// energy; the SI value is J
    Energy,
    /// surface or sheet resistance; the SI value is Ohm (Ohm/sq)
    SheetResistance,
    /// sheet or kinetic inductance; the SI value is H (H/sq)
    SheetInductance,
    /// stress-charge piezoelectric coefficient; the SI value is C/m^2
    PiezoelectricStress,
    /// strain-charge piezoelectric coefficient; the SI value is C/N (== m/V)
    PiezoelectricStrain,
    /// linear electro-optic (Pockels) coefficient; the SI value is m/V
    ElectroOptic,
    /// dynamic viscosity; the SI value is Pa*s
    DynamicViscosity,
    /// electric field strength; the SI value is V/m
    ElectricField,
    /// 2D surface / sheet conductivity; the SI value is S (S*sq)
    SheetConductivity,
    /// time or relaxation time; the SI value is s
    Time,
    /// velocity; the SI value is m/s
    Velocity,
    /// magnetic flux density; the SI value is T
    MagneticFluxDensity,
    /// capacitance per unit area; the SI value is F/m^2
    CapacitancePerArea,
    /// nonlinear Kerr refractive index n_2; the SI value is m^2/W
    NonlinearIndex,
    /// third-order nonlinear susceptibility chi^(3); the SI value is m^2/V^2
    Chi3Susceptibility,
    /// two-photon absorption coefficient beta_TPA; the SI value is m/W
    TwoPhotonAbsorption,
    /// cross-section or area (e.g. free-carrier absorption sigma_FCA); the SI value is m^2
    Area,
    /// volume or free-carrier dispersion coefficient k_FCD; the SI value is m^3
    Volume,
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
            Dim::Frequency => "a frequency or energy, for example eV, THz or rad/s",
            Dim::FrequencyHz => "a frequency in Hz, for example GHz or MHz",
            Dim::Energy => "an energy, for example meV, ueV or J",
            Dim::SheetResistance => "a surface or sheet resistance, for example Ohm or nOhm/sq",
            Dim::SheetInductance => "a sheet inductance, for example pH/sq or nH/sq",
            Dim::PiezoelectricStress => "a stress-charge piezoelectric coefficient, for example C/m^2",
            Dim::PiezoelectricStrain => "a strain-charge piezoelectric coefficient, for example pC/N",
            Dim::ElectroOptic => "an electro-optic or chi^(2) coefficient, for example pm/V",
            Dim::DynamicViscosity => "a dynamic viscosity, for example Pa*s or mPa*s",
            Dim::ElectricField => "an electric field strength, for example V/m or kV/cm",
            Dim::SheetConductivity => "a 2D sheet conductivity, for example S, mS or uS",
            Dim::Time => "a time or relaxation time, for example ps, fs or ns",
            Dim::Velocity => "a velocity, for example m/s or cm/s",
            Dim::MagneticFluxDensity => "a magnetic flux density, for example T or mT",
            Dim::CapacitancePerArea => "a capacitance per unit area, for example F/m^2 or uF/cm^2",
            Dim::NonlinearIndex => "a nonlinear refractive index n_2, for example m^2/W or cm^2/W",
            Dim::Chi3Susceptibility => "a third-order susceptibility chi^(3), for example m^2/V^2 or pm^2/V^2",
            Dim::TwoPhotonAbsorption => "a two-photon absorption coefficient, for example m/W or cm/GW",
            Dim::Area => "an area or cross-section, for example m^2 or cm^2",
            Dim::Volume => "a volume or carrier-dispersion coefficient, for example m^3 or cm^3",
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
            '\u{3a9}' | '\u{2126}' => s.push_str("Ohm"),
            c => s.push(c),
        }
    }
    s
}

fn prefix(p: &str) -> Option<f64> {
    Some(match p {
        "f" => 1e-15,
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
        "P" => 1e15,
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
        "rad/s" => return Some((Dim::Frequency, 1.0)),
        "1/cm" | "cm^-1" | "cm-1" => return Some((Dim::Frequency, 2.0 * std::f64::consts::PI * 299_792_458.0 * 100.0)),
        _ => {}
    }
    if let Some(f) = with_prefix(&u, "eV") {
        return Some((Dim::Frequency, f * 1.602_176_634e-19 / 1.054_571_817e-34));
    }
    if let Some(f) = with_prefix(&u, "Hz") {
        return Some((Dim::Frequency, f * 2.0 * std::f64::consts::PI));
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
    if let Some(f) = with_prefix(&u, "Ohm/sq")
        .or_else(|| with_prefix(&u, "ohm/sq"))
        .or_else(|| with_prefix(&u, "Ohm"))
        .or_else(|| with_prefix(&u, "ohm"))
    {
        return Some((Dim::SheetResistance, f));
    }
    if let Some(f) = with_prefix(&u, "H/sq").or_else(|| with_prefix(&u, "H")) {
        return Some((Dim::SheetInductance, f));
    }
    if let Some(f) = with_prefix(&u, "C/m^2").or_else(|| with_prefix(&u, "C/m2")) {
        return Some((Dim::PiezoelectricStress, f));
    }
    if let Some(f) = with_prefix(&u, "C/N") {
        return Some((Dim::PiezoelectricStrain, f));
    }
    if let Some(f) = with_prefix(&u, "m/V") {
        return Some((Dim::ElectroOptic, f));
    }
    if let Some(f) = with_prefix(&u, "Pas") {
        return Some((Dim::DynamicViscosity, f));
    }
    if let Some(f) = with_prefix(&u, "V/m") {
        return Some((Dim::ElectricField, f));
    }
    if let Some(f) = with_prefix(&u, "V/cm") {
        return Some((Dim::ElectricField, f * 1e2));
    }
    if let Some(f) = with_prefix(&u, "Ssq")
        .or_else(|| with_prefix(&u, "S"))
        .or_else(|| with_prefix(&u, "Ohm^-1"))
        .or_else(|| with_prefix(&u, "ohm^-1"))
    {
        return Some((Dim::SheetConductivity, f));
    }
    if let Some(f) = with_prefix(&u, "m/s") {
        return Some((Dim::Velocity, f));
    }
    if let Some(f) = with_prefix(&u, "F/m^2").or_else(|| with_prefix(&u, "F/m2")) {
        return Some((Dim::CapacitancePerArea, f));
    }
    if let Some(f) = with_prefix(&u, "F/cm^2").or_else(|| with_prefix(&u, "F/cm2")) {
        return Some((Dim::CapacitancePerArea, f * 1e4));
    }
    if let Some(rest) = u.strip_suffix("/GW") {
        if rest == "cm^2" || rest == "cm2" {
            return Some((Dim::NonlinearIndex, 1e-13));
        }
        if rest == "cm" {
            return Some((Dim::TwoPhotonAbsorption, 1e-11));
        }
    }
    if let Some(rest) = u.strip_suffix("/W") {
        if let Some(head) = rest.strip_suffix("^2").or_else(|| rest.strip_suffix('2')) {
            if let Some(f) = with_prefix(head, "m") {
                return Some((Dim::NonlinearIndex, f * f));
            }
        }
        if let Some(f) = with_prefix(rest, "m") {
            return Some((Dim::TwoPhotonAbsorption, f));
        }
    }
    if let Some(rest) = u.strip_suffix("/V^2").or_else(|| u.strip_suffix("/V2")) {
        if let Some(head) = rest.strip_suffix("^2").or_else(|| rest.strip_suffix('2')) {
            if let Some(f) = with_prefix(head, "m") {
                return Some((Dim::Chi3Susceptibility, f * f));
            }
        }
    }
    if let Some(head) = u.strip_suffix("^2") {
        if let Some(f) = with_prefix(head, "m") {
            return Some((Dim::Area, f * f));
        }
    }
    if u == "m2" || u == "cm2" || u == "mm2" || u == "um2" || u == "nm2" || u == "pm2" {
        if let Some(f) = with_prefix(u.strip_suffix('2').unwrap_or(""), "m") {
            return Some((Dim::Area, f * f));
        }
    }
    if let Some(head) = u.strip_suffix("^3") {
        if let Some(f) = with_prefix(head, "m") {
            return Some((Dim::Volume, f * f * f));
        }
    }
    if u == "m3" || u == "cm3" || u == "mm3" || u == "um3" || u == "nm3" {
        if let Some(f) = with_prefix(u.strip_suffix('3').unwrap_or(""), "m") {
            return Some((Dim::Volume, f * f * f));
        }
    }
    if let Some(f) = with_prefix(&u, "T") {
        return Some((Dim::MagneticFluxDensity, f));
    }
    if let Some(f) = with_prefix(&u, "s") {
        return Some((Dim::Time, f));
    }
    if let Some(f) = with_prefix(&u, "J") {
        return Some((Dim::Energy, f));
    }
    if let Some(f) = with_prefix(&u, "K") {
        return Some((Dim::Temperature, f));
    }
    if let Some(f) = with_prefix(&u, "m") {
        return Some((Dim::Length, f));
    }
    None
}

/// `value` in `unit` converted to SI, if `unit` is a unit of dimension `expect`.
pub fn to_si(value: f64, unit: &str, expect: Dim) -> Result<f64, String> {
    let u = normalise(unit);
    match expect {
        Dim::Energy => {
            if let Some(f) = with_prefix(&u, "eV") {
                return Ok(value * f * 1.602_176_634e-19);
            }
            if let Some(f) = with_prefix(&u, "J") {
                return Ok(value * f);
            }
        }
        Dim::FrequencyHz => {
            if let Some(f) = with_prefix(&u, "Hz") {
                return Ok(value * f);
            }
        }
        Dim::PiezoelectricStrain => {
            // 1 C/N == 1 m/V in SI; accept both notations (e.g. pC/N and pm/V) for strain-charge d_iJ.
            if let Some(f) = with_prefix(&u, "C/N").or_else(|| with_prefix(&u, "m/V")) {
                return Ok(value * f);
            }
        }
        _ => {}
    }
    let (dim, f) = lookup(unit).ok_or_else(|| format!("the unit '{unit}' is not known; expected {}", expect.describe()))?;
    if dim != expect {
        return Err(format!("the unit '{unit}' does not fit here; expected {}", expect.describe()));
    }
    if expect == Dim::Temperature {
        return Ok(if u == "degC" || u == "C" { value + 273.15 } else { value * f });
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
    if let Ok(v) = t.parse::<f64>() {
        return Ok((v, String::new()));
    }
    for end in (1..t.len()).rev().filter(|&i| t.is_char_boundary(i)) {
        if let Ok(v) = t[..end].parse::<f64>() {
            let next = t[end..].chars().next().unwrap_or(' ');
            if !(next.is_ascii_digit() || next == '.' || next == 'e' || next == 'E') {
                return Ok((v, t[end..].trim().to_string()));
            }
        }
    }
    Err(format!("'{s}' is not a number with a unit; write it like \"496.8 GPa\""))
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
        assert!((to_si(5.0, "ppm/K", Dim::PerKelvin).unwrap() - 5.0e-6).abs() < 1e-18);
        assert!((to_si(20.0, "degC", Dim::Temperature).unwrap() - 293.15).abs() < 1e-12);
        assert_eq!(to_si(1.0, "1/GPa", Dim::Compliance).unwrap(), 1e-9);
        assert_eq!(to_si(3.0, "W/(m*K)", Dim::ThermalConductivity).unwrap(), 3.0);
        assert!((to_si(180.0, "ueV", Dim::Energy).unwrap() - 180e-6 * 1.602_176_634e-19).abs() < 1e-34);
        assert_eq!(to_si(10.0, "GHz", Dim::FrequencyHz).unwrap(), 10.0e9);
        assert!((to_si(12.0, "nOhm", Dim::SheetResistance).unwrap() - 12e-9).abs() < 1e-22);
        assert!((to_si(12.0, "nOhm/sq", Dim::SheetResistance).unwrap() - 12e-9).abs() < 1e-22);
        assert!((to_si(0.09, "pH/sq", Dim::SheetInductance).unwrap() - 0.09e-12).abs() < 1e-25);
        assert_eq!(to_si(2.45, "C/m^2", Dim::PiezoelectricStress).unwrap(), 2.45);
        assert!((to_si(69.2, "pC/N", Dim::PiezoelectricStrain).unwrap() - 69.2e-12).abs() < 1e-24);
        assert!((to_si(30.8, "pm/V", Dim::ElectroOptic).unwrap() - 30.8e-12).abs() < 1e-24);
        assert!((to_si(5.9, "mPa*s", Dim::DynamicViscosity).unwrap() - 5.9e-3).abs() < 1e-16);
        assert_eq!(to_si(1.0, "kV/cm", Dim::ElectricField).unwrap(), 1e5);
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
