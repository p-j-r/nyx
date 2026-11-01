//! The session: content-hashed handles for geometry, mesh, problem and result, a journal of every call, and the small header that goes
//! into every result so a host can tell whether the session changed under it.
//!
//! Handles are `g_`, `m_`, `p_`, `r_` plus hex digits of an FNV-1a hash of the inputs that make the artifact, so asking for the same thing
//! again returns the same handle without recomputing. Meshes and solver results stay in memory; what is written to the session folder is
//! the journal (`journal.jsonl`), the geometry model of every geometry (`g_*.geom.json`, the `geom` crate's own format) and a digest of every result.

use femcore::geomdrive::MeshedGeom;
use femcore::problem::ProblemSpec;
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::params::ToolError;

/// FNV-1a, 64 bit: stable across runs, platforms and compiler versions (std's hashers promise none of that).
pub fn fnv64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

#[derive(Clone, Debug, Serialize)]
pub struct GeomParams {
    pub thickness: f64,
    pub width: f64,
    pub box_thickness: f64,
    pub air_cladding: bool,
    pub margin: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct MeshParams {
    pub h_core: f64,
    pub h_near: f64,
    pub h_far: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ProblemParams {
    pub wavelength: f64,
    pub n_modes: usize,
    pub n_target: f64,
    pub auto_target: bool,
    pub pol: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Family {
    Te,
    Tm,
    Other,
}

impl Family {
    pub fn name(self) -> &'static str {
        match self {
            Family::Te => "TE",
            Family::Tm => "TM",
            Family::Other => "other",
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ModeRec {
    pub neff_re: f64,
    pub neff_im: f64,
    pub class: String,
    pub family: Family,
    pub te_fraction: f64,
    pub long_fraction: f64,
    pub dop: f64,
    pub degenerate: bool,
}

pub struct GeomItem {
    pub params: GeomParams,
}

pub struct MeshItem {
    pub geom: String,
    pub params: MeshParams,
    pub meshed: MeshedGeom,
}

pub struct ProblemItem {
    pub mesh: String,
    pub params: ProblemParams,
    pub spec: ProblemSpec,
}

pub struct ResultItem {
    pub problem: String,
    pub modes: Vec<ModeRec>,
    pub seconds: f64,
    pub dofs: usize,
}

pub enum Item {
    Geom(GeomItem),
    Mesh(MeshItem),
    Problem(ProblemItem),
    Result(ResultItem),
}

impl Item {
    fn kind(&self) -> &'static str {
        match self {
            Item::Geom(_) => "geometry",
            Item::Mesh(_) => "mesh",
            Item::Problem(_) => "problem",
            Item::Result(_) => "result",
        }
    }
}

pub struct Workspace {
    pub dir: PathBuf,
    pub session: String,
    revision: u64,
    items: BTreeMap<String, Item>,
    keys: BTreeMap<String, String>,
    journal_failed: bool,
}

impl Workspace {
    /// A new session in `dir/<session id>`; the folder is created on the first write.
    pub fn new(dir: PathBuf) -> Workspace {
        let secs = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
        let session = format!("s{secs:x}{:02x}", std::process::id() & 0xff);
        Workspace { dir: dir.join(&session), session, revision: 0, items: BTreeMap::new(), keys: BTreeMap::new(), journal_failed: false }
    }

    pub fn revision(&self) -> u64 {
        self.revision
    }

    /// Hash of the set of handles: changes when anything is added.
    pub fn project_hash(&self) -> String {
        let all: Vec<&str> = self.items.keys().map(String::as_str).collect();
        format!("{:08x}", fnv64(&all.join(",")) as u32)
    }

    pub fn header(&self) -> Value {
        json!({ "session": self.session, "revision": self.revision, "project_hash": self.project_hash() })
    }

    /// Adds an item under the handle derived from `key`, unless that handle already exists. Returns the handle and whether it was there before.
    pub fn insert(&mut self, prefix: char, key: &str, item: Item) -> (String, bool) {
        let h = fnv64(&format!("{prefix}|{key}"));
        for digits in [8usize, 12, 16] {
            let handle = format!("{prefix}_{:0width$x}", h >> (64 - 4 * digits as u32), width = digits);
            match self.keys.get(&handle) {
                Some(k) if k == key => return (handle, true),
                Some(_) => continue,
                None => {
                    self.keys.insert(handle.clone(), key.to_string());
                    self.items.insert(handle.clone(), item);
                    self.revision += 1;
                    return (handle, false);
                }
            }
        }
        unreachable!("64-bit hash collision on three widths")
    }

    fn missing(&self, handle: &str, want: &str) -> ToolError {
        let prefix = handle.chars().next().unwrap_or('?');
        let have: Vec<&String> = self.items.iter().filter(|(_, it)| it.kind() == want).map(|(k, _)| k).collect();
        let msg = if self.items.contains_key(handle) { format!("'{handle}' is not a {want} handle") } else { format!("no {want} '{handle}' in this session") };
        ToolError::new("unknown_handle", msg).with("got", handle).with("valid", json!(have)).with("kind", want).with("prefix", prefix.to_string())
    }

    pub fn geom(&self, h: &str) -> Result<&GeomItem, ToolError> {
        match self.items.get(h) {
            Some(Item::Geom(g)) => Ok(g),
            _ => Err(self.missing(h, "geometry")),
        }
    }
    pub fn mesh(&self, h: &str) -> Result<&MeshItem, ToolError> {
        match self.items.get(h) {
            Some(Item::Mesh(m)) => Ok(m),
            _ => Err(self.missing(h, "mesh")),
        }
    }
    pub fn problem(&self, h: &str) -> Result<&ProblemItem, ToolError> {
        match self.items.get(h) {
            Some(Item::Problem(p)) => Ok(p),
            _ => Err(self.missing(h, "problem")),
        }
    }
    pub fn result(&self, h: &str) -> Result<&ResultItem, ToolError> {
        match self.items.get(h) {
            Some(Item::Result(r)) => Ok(r),
            _ => Err(self.missing(h, "result")),
        }
    }

    pub fn handles(&self) -> Vec<(String, &'static str)> {
        self.items.iter().map(|(k, v)| (k.clone(), v.kind())).collect()
    }

    /// Writes a file into the session folder (best effort: a read-only disk must not stop a calculation).
    pub fn write_file(&mut self, name: &str, text: &str) {
        if std::fs::create_dir_all(&self.dir).and_then(|_| std::fs::write(self.dir.join(name), text)).is_err() {
            self.journal_failed = true;
        }
    }

    /// Appends one line to the journal: which tool, with which arguments, what it produced.
    pub fn journal(&mut self, tool: &str, args: &Value, outcome: &str, summary: Value) {
        let line = json!({ "rev": self.revision, "tool": tool, "args": args, "outcome": outcome, "summary": summary });
        let ok = std::fs::create_dir_all(&self.dir)
            .and_then(|_| std::fs::OpenOptions::new().create(true).append(true).open(self.dir.join("journal.jsonl")))
            .and_then(|mut f| writeln!(f, "{line}"))
            .is_ok();
        if !ok {
            self.journal_failed = true;
        }
    }

    pub fn journal_failed(&self) -> bool {
        self.journal_failed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv_is_stable() {
        // reference values of FNV-1a 64
        assert_eq!(fnv64(""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv64("a"), 0xaf63_dc4c_8601_ec8c);
    }

    #[test]
    fn the_same_inputs_give_the_same_handle_and_a_new_one_bumps_the_revision() {
        let mut ws = Workspace::new(std::env::temp_dir().join("nyx-session-test"));
        let item = || Item::Geom(GeomItem { params: GeomParams { thickness: 3e-7, width: 5e-7, box_thickness: 3e-6, air_cladding: true, margin: 1e-6 } });
        let (a, again) = ws.insert('g', "k1", item());
        assert!(!again && a.starts_with("g_") && a.len() == 10);
        let r = ws.revision();
        let (b, again) = ws.insert('g', "k1", item());
        assert!(again && a == b && ws.revision() == r);
        let (c, _) = ws.insert('g', "k2", item());
        assert_ne!(a, c);
        assert_eq!(ws.revision(), r + 1);
        assert!(ws.geom(&a).is_ok());
        assert_eq!(ws.mesh(&a).err().unwrap().code, "unknown_handle");
    }
}
