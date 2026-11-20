// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Small shared helpers: UTC time without a time crate, short time-ordered ids, atomic file writes, hashes and slugs.

use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

/// Seconds since 1970-01-01 UTC.
pub fn now_secs() -> i64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0)
}

/// Days since 1970-01-01 to (year, month, day), proleptic Gregorian (Hinnant's algorithm).
fn civil(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn days_from_civil(y: i64, m: u32, d: u32) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y.rem_euclid(400);
    let mp = if m > 2 { m - 3 } else { m + 9 } as i64;
    let doy = (153 * mp + 2) / 5 + d as i64 - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// `2026-11-01T00:08:15Z`
pub fn iso(secs: i64) -> String {
    let (y, m, d) = civil(secs.div_euclid(86_400));
    let r = secs.rem_euclid(86_400);
    format!("{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z", r / 3600, r % 3600 / 60, r % 60)
}

/// Seconds of an `iso` string (`None` when it is not one).
pub fn parse_iso(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    if b.len() != 20 || b[4] != b'-' || b[7] != b'-' || b[10] != b'T' || b[13] != b':' || b[16] != b':' || b[19] != b'Z' {
        return None;
    }
    let n = |a: usize, z: usize| s.get(a..z)?.parse::<i64>().ok();
    let (y, mo, d, h, mi, se) = (n(0, 4)?, n(5, 7)?, n(8, 10)?, n(11, 13)?, n(14, 16)?, n(17, 19)?);
    if !(1..=12).contains(&mo) || !(1..=31).contains(&d) || h > 23 || mi > 59 || se > 60 {
        return None;
    }
    Some(days_from_civil(y, mo as u32, d as u32) * 86_400 + h * 3600 + mi * 60 + se)
}

/// "2 h ago", "yesterday", "3 weeks ago" for a time `then` seen from `now`.
pub fn ago(then: i64, now: i64) -> String {
    let s = (now - then).max(0);
    let plural = |n: i64, unit: &str| if n == 1 { format!("{n} {unit} ago") } else { format!("{n} {unit}s ago") };
    match s {
        0..=59 => "just now".into(),
        60..=3599 => format!("{} min ago", s / 60),
        3600..=86_399 => format!("{} h ago", s / 3600),
        86_400..=172_799 => "yesterday".into(),
        172_800..=1_209_599 => format!("{} days ago", s / 86_400),
        1_210_000..=2_591_999 => plural(s / 604_800, "week"),
        2_592_000..=31_535_999 => if s < 5_184_000 { "last month".into() } else { plural(s / 2_592_000, "month") },
        _ => plural(s / 31_536_000, "year"),
    }
}

const B32: &[u8; 32] = b"0123456789abcdefghjkmnpqrstvwxyz";

fn b32(mut v: u64, n: usize) -> String {
    let mut out = vec![b'0'; n];
    for i in (0..n).rev() {
        out[i] = B32[(v & 31) as usize];
        v >>= 5;
    }
    String::from_utf8(out).unwrap()
}

static COUNTER: AtomicU64 = AtomicU64::new(0);

/// A short id: 7 base32 characters of the time in seconds (so sorting ids sorts by creation) and 2 of entropy.
pub fn short_id() -> String {
    let t = now_secs().max(0) as u64;
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.subsec_nanos() as u64).unwrap_or(0);
    let c = COUNTER.fetch_add(1, Ordering::Relaxed);
    let h = blake3::hash(&[&nanos.to_le_bytes()[..], &c.to_le_bytes(), &(std::process::id() as u64).to_le_bytes()].concat());
    let e = u16::from_le_bytes([h.as_bytes()[0], h.as_bytes()[1]]) as u64 & 1023;
    format!("{}{}", b32(t, 7), b32(e, 2))
}

/// Three entropy characters for a run id (`r0007_7hb`).
pub fn run_suffix() -> String {
    let id = short_id();
    id[6..9].to_string()
}

/// Lower-case ascii slug for a folder name: letters and digits, single dashes, at most 40 characters.
pub fn slug(s: &str) -> String {
    let mut out = String::new();
    for ch in s.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
        } else if !out.ends_with('-') && !out.is_empty() {
            out.push('-');
        }
    }
    let out: String = out.trim_end_matches('-').chars().take(40).collect();
    let out = out.trim_end_matches('-').to_string();
    if out.is_empty() { "workflow".into() } else { out }
}

/// `b3:` and 32 hex digits of the BLAKE3 hash of `bytes`.
pub fn hash_bytes(bytes: &[u8]) -> String {
    format!("b3:{}", &blake3::hash(bytes).to_hex()[..32])
}

pub fn hash_file(path: &Path) -> std::io::Result<String> {
    Ok(hash_bytes(&std::fs::read(path)?))
}

/// Writes `data` to a temporary file next to `path` and renames it over `path`, so a crash never leaves half a file.
pub fn write_atomic(path: &Path, data: &[u8]) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    std::fs::write(&tmp, data)?;
    std::fs::rename(&tmp, path).inspect_err(|_| {
        let _ = std::fs::remove_file(&tmp);
    })
}

pub fn io_err(path: &Path, e: impl std::fmt::Display) -> String {
    format!("{}: {e}", path.display())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iso_round_trips_and_knows_leap_years() {
        for t in [0, 951_782_400, 1_782_864_495, 4_102_444_800] {
            assert_eq!(parse_iso(&iso(t)), Some(t), "{t}");
        }
        assert_eq!(iso(951_782_400), "2000-02-29T00:00:00Z");
        assert_eq!(parse_iso("2026-13-01T00:00:00Z"), None);
        assert_eq!(parse_iso("garbage"), None);
    }

    #[test]
    fn ids_sort_by_time_and_differ() {
        let a = short_id();
        let b = short_id();
        assert_ne!(a, b);
        assert_eq!(a.len(), 9);
        assert!(b32(1000, 7) < b32(1001, 7));
    }

    #[test]
    fn slugs_are_folder_safe() {
        assert_eq!(slug("SOI TE modes: sweep!"), "soi-te-modes-sweep");
        assert_eq!(slug("..//"), "workflow");
        assert_eq!(ago(0, 90_000), "yesterday");
    }
}
