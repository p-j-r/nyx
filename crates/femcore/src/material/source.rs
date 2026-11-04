//! Provenance of material data: a `source` key at the top of a file and in each table says where the numbers come from.
//! It is a citation text, or a table `{ cite, url, doi, verified }` (`verified` is the date `YYYY-MM-DD` on which someone compared the
//! values with the cited publication), or a list of those. Data that were not checked say so in the citation text.

use super::reader::{MaterialError, Section};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Source {
    pub cite: String,
    pub url: Option<String>,
    pub doi: Option<String>,
    /// date `YYYY-MM-DD` of the comparison with the publication
    pub verified: Option<String>,
}

impl Source {
    pub fn text(cite: &str) -> Source {
        Source { cite: cite.to_string(), url: None, doi: None, verified: None }
    }

    /// True if the entry has a verification date.
    pub fn is_verified(&self) -> bool {
        self.verified.is_some()
    }

    /// True if the citation text itself flags the data as unchecked.
    pub fn flags_uncertainty(&self) -> bool {
        let c = self.cite.to_ascii_lowercase();
        c.contains("unverified") || c.contains("approximate")
    }
}

fn valid_date(d: &str) -> bool {
    let b = d.as_bytes();
    if b.len() != 10 || b[4] != b'-' || b[7] != b'-' {
        return false;
    }
    let num = |r: std::ops::Range<usize>| d[r].parse::<u32>().ok();
    match (num(0..4), num(5..7), num(8..10)) {
        (Some(y), Some(m), Some(day)) => y >= 1900 && (1..=12).contains(&m) && (1..=31).contains(&day),
        _ => false,
    }
}

fn one(s: &Section, key: &str, v: &toml::Value) -> Result<Source, MaterialError> {
    match v {
        toml::Value::String(t) => {
            if t.trim().is_empty() {
                return Err(s.err(key, "must not be empty; cite where the numbers come from"));
            }
            Ok(Source::text(t))
        }
        toml::Value::Table(t) => {
            let get = |k: &str| -> Result<Option<String>, MaterialError> {
                match t.get(k) {
                    None => Ok(None),
                    Some(toml::Value::String(x)) => Ok(Some(x.clone())),
                    Some(_) => Err(s.err(key, format!("'{k}' must be text in quotes"))),
                }
            };
            for k in t.keys() {
                if !["cite", "url", "doi", "verified"].contains(&k.as_str()) {
                    return Err(s.err(key, format!("unknown entry '{k}' in a source; use cite, url, doi, verified")));
                }
            }
            let cite = get("cite")?.ok_or_else(|| s.err(key, "a source table needs cite = \"...\" (the citation text)"))?;
            let (url, doi, verified) = (get("url")?, get("doi")?, get("verified")?);
            if let Some(u) = &url {
                if !(u.starts_with("http://") || u.starts_with("https://")) {
                    return Err(s.err(key, format!("url \"{u}\" must start with http:// or https://")));
                }
            }
            if let Some(d) = &doi {
                if !d.starts_with("10.") {
                    return Err(s.err(key, format!("doi \"{d}\" must start with 10. (write the DOI itself, not a link)")));
                }
            }
            if let Some(v) = &verified {
                if !valid_date(v) {
                    return Err(s.err(key, format!("verified \"{v}\" must be a date written YYYY-MM-DD")));
                }
            }
            Ok(Source { cite, url, doi, verified })
        }
        _ => Err(s.err(key, "must be text, a table { cite, url, doi, verified }, or a list of those")),
    }
}

/// Reads `key` of the section as a list of sources (empty if the key is absent).
pub fn read_sources(s: &mut Section, key: &str) -> Result<Vec<Source>, MaterialError> {
    match s.value(key) {
        None => Ok(Vec::new()),
        Some(toml::Value::Array(a)) => a.iter().map(|v| one(s, key, v)).collect(),
        Some(v) => Ok(vec![one(s, key, v)?]),
    }
}

/// The TOML value for a list of sources, or `None` if it is empty.
pub fn write_sources(v: &[Source]) -> Option<toml::Value> {
    let one = |x: &Source| -> toml::Value {
        if x.url.is_none() && x.doi.is_none() && x.verified.is_none() {
            return x.cite.clone().into();
        }
        let mut t = toml::Table::new();
        t.insert("cite".into(), x.cite.clone().into());
        if let Some(u) = &x.url {
            t.insert("url".into(), u.clone().into());
        }
        if let Some(d) = &x.doi {
            t.insert("doi".into(), d.clone().into());
        }
        if let Some(d) = &x.verified {
            t.insert("verified".into(), d.clone().into());
        }
        toml::Value::Table(t)
    };
    match v.len() {
        0 => None,
        1 => Some(one(&v[0])),
        _ => Some(toml::Value::Array(v.iter().map(one).collect())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dates_are_checked() {
        assert!(valid_date("2026-11-04") && !valid_date("2026-13-04") && !valid_date("04-11-2026") && !valid_date("2026-11-4"));
    }
}
