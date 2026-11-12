//! `.nyx` files: a workflow as one small, deterministic ZIP, and the safe way back in.
//!
//! **Export** writes the files in sorted order with a fixed timestamp, fixed permissions and fixed compression, plus `nyx-archive.json`
//! listing every file with its BLAKE3 hash, so the same workflow always gives the same bytes and a damaged file is detected. The share level
//! decides what goes in (`Recipe`: the set-up, run table, versions, notes and run logs; `Results`: also the heavy result files; `Meshes`: also
//! mesh files). Absolute paths in specs are reduced to file names when asked to. Scripts under `analysis/` travel but are never enabled.
//!
//! **Import** treats the file as hostile: only the layout nyx itself writes is accepted (anything else rejects the whole file), names are
//! checked against zip-slip (`..`, absolute paths, drive letters, backslashes, reserved Windows names, duplicates that differ only in case),
//! symbolic links are refused, sizes are counted while decompressing (the sizes in the headers are not trusted), entry, per-file and total
//! caps and a compression-ratio cap apply, and the content hashes must match. It extracts into a temporary folder and renames it into place
//! only when everything passed, under a fresh id, with scripts switched off.

use crate::flow::WorkflowDir;
use crate::model::*;
use crate::util::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::{Read, Seek, Write};
use std::path::{Path, PathBuf};

const MANIFEST: &str = "nyx-archive.json";

#[derive(Clone, Debug)]
pub struct Limits {
    pub max_entries: usize,
    pub max_file: u64,
    pub max_total: u64,
    /// uncompressed / compressed, checked for files above 1 MiB
    pub max_ratio: u64,
}

impl Default for Limits {
    fn default() -> Self {
        Limits { max_entries: 20_000, max_file: 2 << 30, max_total: 4 << 30, max_ratio: 1000 }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
struct FileEntry {
    path: String,
    size: u64,
    hash: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ArchiveManifest {
    format: u32,
    kind: String,
    level: ShareLevel,
    /// version of the nyx that wrote it (no names, paths or times: the file must not depend on who or when)
    nyx_version: String,
    files: Vec<FileEntry>,
}

// ---------------------------------------------------------------------------------------------------------- names

const RESERVED: &[&str] = &["con", "prn", "aux", "nul", "com1", "com2", "com3", "com4", "com5", "com6", "com7", "com8", "com9", "lpt1", "lpt2", "lpt3", "lpt4", "lpt5", "lpt6", "lpt7", "lpt8", "lpt9"];

/// `Ok` when `name` is a plain relative path that cannot leave the folder it is extracted into, on any platform.
pub fn check_name(name: &str) -> Result<(), String> {
    let bad = |why: &str| Err(format!("unsafe file name '{}': {why}", name.escape_debug()));
    if name.is_empty() || name.len() > 240 {
        return bad("empty or too long");
    }
    if name.contains('\\') || name.starts_with('/') || name.contains(':') || name.contains('\0') || name.chars().any(|c| c.is_control()) {
        return bad("not a plain relative path");
    }
    for part in name.split('/') {
        if part.is_empty() || part == "." || part == ".." {
            return bad("empty, '.' or '..' component");
        }
        if part.len() > 120 || part.ends_with('.') || part.ends_with(' ') {
            return bad("component too long or ends with a dot or space");
        }
        let stem = part.split('.').next().unwrap_or("").to_lowercase();
        if RESERVED.contains(&stem.as_str()) {
            return bad("reserved name");
        }
    }
    Ok(())
}

fn simple(s: &str, first: char) -> bool {
    s.starts_with(first) && s.len() <= 40 && s.chars().all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
}

fn simple_file(s: &str, exts: &[&str]) -> bool {
    let Some((stem, ext)) = s.rsplit_once('.') else { return false };
    !stem.is_empty() && stem.len() <= 80 && stem.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.')) && exts.contains(&ext)
}

/// The layout nyx writes; nothing else is accepted from a `.nyx` file.
pub fn allowed(name: &str) -> bool {
    let p: Vec<&str> = name.split('/').collect();
    match p.as_slice() {
        [MANIFEST] | ["workflow.json"] | ["history.jsonl"] | ["notes.md"] | [".gitignore"] => true,
        ["specs", f] => f.ends_with(".nyx.json") && f.len() < 100 && f.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '.' | '_')),
        ["analysis", f] => simple_file(f, &["rhai", "py", "md", "txt", "json", "csv"]),
        [e, "experiment.nyx.json"] | [e, "runs.jsonl"] => simple(e, 'e'),
        [e, "runs", r, f] => simple(e, 'e') && simple(r, 'r') && simple_file(f, &["txt", "h5", "hdf5", "json", "csv", "mesh", "png"]),
        _ => false,
    }
}

fn is_script(name: &str) -> bool {
    name.starts_with("analysis/")
}

fn is_mesh(name: &str) -> bool {
    name.ends_with(".mesh")
}

fn is_heavy(name: &str) -> bool {
    name.contains("/runs/") && !name.ends_with("/log.txt")
}

// ---------------------------------------------------------------------------------------------------------- export

/// Files an export of `level` would contain (relative name, path on disk), sorted by name. `only` restricts to some experiments.
fn plan(w: &WorkflowDir, level: ShareLevel, only: Option<&[String]>) -> Result<Vec<(String, PathBuf)>, String> {
    let keep = |eid: &str| only.is_none_or(|o| o.iter().any(|x| x == eid));
    let mut files: BTreeMap<String, PathBuf> = BTreeMap::new();
    for f in ["history.jsonl", "notes.md", ".gitignore"] {
        if w.root.join(f).is_file() {
            files.insert(f.into(), w.root.join(f));
        }
    }
    let mut hashes: Vec<String> = vec![];
    for e in w.wf.experiments.iter().filter(|e| keep(&e.id)) {
        files.insert(format!("{}/experiment.nyx.json", e.id), w.spec_path(&e.id));
        if w.runs_path(&e.id).is_file() {
            files.insert(format!("{}/runs.jsonl", e.id), w.runs_path(&e.id));
            hashes.extend(w.runs(&e.id).0.into_iter().map(|r| r.spec));
        }
        hashes.push(e.spec_hash.clone());
        let runs_dir = w.root.join(&e.id).join("runs");
        if let Ok(rd) = std::fs::read_dir(&runs_dir) {
            for run in rd.filter_map(|r| r.ok()).filter(|r| r.path().is_dir()) {
                let rid = run.file_name().to_string_lossy().to_string();
                for f in std::fs::read_dir(run.path()).map_err(|e| io_err(&run.path(), e))?.filter_map(|f| f.ok()).filter(|f| f.path().is_file()) {
                    let name = format!("{}/runs/{rid}/{}", e.id, f.file_name().to_string_lossy());
                    let wanted = match level {
                        ShareLevel::Recipe => !is_heavy(&name),
                        ShareLevel::Results => !is_mesh(&name),
                        ShareLevel::Meshes => true,
                    };
                    if wanted && allowed(&name) {
                        files.insert(name, f.path());
                    }
                }
            }
        }
    }
    let ids: Vec<&str> = w.wf.experiments.iter().filter(|e| keep(&e.id)).map(|e| e.id.as_str()).collect();
    hashes.extend(w.history().into_iter().filter(|v| ids.contains(&v.experiment.as_str())).flat_map(|v| [Some(v.spec_hash), v.prev_hash]).flatten());
    hashes.sort();
    hashes.dedup();
    for h in hashes {
        let p = w.frozen_path(&h);
        if p.is_file() {
            files.insert(format!("specs/{}", p.file_name().unwrap().to_string_lossy()), p);
        }
    }
    let an = w.root.join("analysis");
    if let Ok(rd) = std::fs::read_dir(&an) {
        for f in rd.filter_map(|f| f.ok()).filter(|f| f.path().is_file()) {
            let name = format!("analysis/{}", f.file_name().to_string_lossy());
            if allowed(&name) {
                files.insert(name, f.path());
            }
        }
    }
    Ok(files.into_iter().collect())
}

/// Sizes in bytes of the files each level would contain (uncompressed; the zip is smaller), for the Share screen.
pub fn sizes(w: &WorkflowDir, only: Option<&[String]>) -> Result<[(ShareLevel, u64); 3], String> {
    let mut out = [(ShareLevel::Recipe, 0), (ShareLevel::Results, 0), (ShareLevel::Meshes, 0)];
    for o in out.iter_mut() {
        o.1 = plan(w, o.0, only)?.iter().map(|(_, p)| std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)).sum();
    }
    Ok(out)
}

/// Absolute paths in string values of a spec reduced to their file names (`C:\Users\me\out\r.h5` becomes `r.h5`).
pub fn strip_paths(xml: &str) -> Result<(String, usize), String> {
    fn walk(v: &mut serde_json::Value, n: &mut usize) {
        match v {
            serde_json::Value::String(s) => {
                let b = s.as_bytes();
                let abs = s.starts_with('/') || s.starts_with("\\\\") || (b.len() > 2 && b[1] == b':' && (b[2] == b'\\' || b[2] == b'/') && b[0].is_ascii_alphabetic());
                if abs {
                    *s = s.rsplit(['/', '\\']).next().unwrap_or("").to_string();
                    *n += 1;
                }
            }
            serde_json::Value::Array(a) => a.iter_mut().for_each(|x| walk(x, n)),
            serde_json::Value::Object(o) => o.values_mut().for_each(|x| walk(x, n)),
            _ => {}
        }
    }
        let mut v: serde_json::Value = serde_json::from_str(xml).map_err(|e| format!("not a nyx project: {e}"))?;
        let mut n = 0;
        walk(&mut v, &mut n);
        if n == 0 {
            return Ok((xml.to_string(), 0));
        }
        Ok((serde_json::to_string_pretty(&v).map_err(|e| e.to_string())?, n))
}

/// How many personal paths an export would strip (for the "check before sending" panel).
pub fn count_personal_paths(w: &WorkflowDir) -> usize {
    plan(w, ShareLevel::Recipe, None).map(|p| p.iter().filter(|(n, _)| n.ends_with(".nyx.json")).filter_map(|(_, f)| std::fs::read_to_string(f).ok()).filter_map(|x| strip_paths(&x).ok()).map(|(_, n)| n).sum()).unwrap_or(0)
}

pub struct ExportOptions<'a> {
    pub level: ShareLevel,
    pub only: Option<&'a [String]>,
    pub strip_paths: bool,
    /// relative names left out of the file (scripts the sender chose not to share)
    pub exclude: &'a [String],
}

/// Writes the `.nyx` of `w` to `out`. The same workflow and options always give the same bytes.
pub fn export<W: Write + Seek>(w: &WorkflowDir, opt: &ExportOptions, out: W) -> Result<(), String> {
    let mut contents: Vec<(String, Vec<u8>)> = vec![];
    for (name, path) in plan(w, opt.level, opt.only)? {
        if opt.exclude.contains(&name) {
            continue;
        }
        let mut data = std::fs::read(&path).map_err(|e| io_err(&path, e))?;
        if opt.strip_paths && name.ends_with(".nyx.json") {
            let (x, n) = strip_paths(&String::from_utf8_lossy(&data))?;
            if n > 0 {
                data = x.into_bytes();
            }
        }
        contents.push((name, data));
    }
    // the manifest of the workflow: only the chosen experiments, scripts off, no pin
    let mut wf = w.wf.clone();
    if let Some(only) = opt.only {
        wf.experiments.retain(|e| only.contains(&e.id));
    }
    wf.pinned = false;
    wf.scripts_approved = false;
    wf.share = opt.level;
    contents.push(("workflow.json".into(), serde_json::to_string_pretty(&wf).map_err(|e| e.to_string())?.into_bytes()));
    contents.sort_by(|a, b| a.0.cmp(&b.0));
    let files = contents.iter().map(|(p, d)| FileEntry { path: p.clone(), size: d.len() as u64, hash: hash_bytes(d) }).collect();
    let am = ArchiveManifest { format: FORMAT, kind: "nyx-archive".into(), level: opt.level, nyx_version: env!("CARGO_PKG_VERSION").into(), files };
    contents.push((MANIFEST.into(), serde_json::to_string_pretty(&am).map_err(|e| e.to_string())?.into_bytes()));
    contents.sort_by(|a, b| a.0.cmp(&b.0));

    let mut z = zip::ZipWriter::new(out);
    let o = zip::write::SimpleFileOptions::default().compression_method(zip::CompressionMethod::Deflated).compression_level(Some(6)).last_modified_time(zip::DateTime::default()).unix_permissions(0o644);
    for (name, data) in contents {
        z.start_file(name, o).map_err(|e| e.to_string())?;
        z.write_all(&data).map_err(|e| e.to_string())?;
    }
    z.finish().map_err(|e| e.to_string())?;
    Ok(())
}

pub fn export_bytes(w: &WorkflowDir, opt: &ExportOptions) -> Result<Vec<u8>, String> {
    let mut buf = std::io::Cursor::new(Vec::new());
    export(w, opt, &mut buf)?;
    Ok(buf.into_inner())
}

pub fn export_file(w: &WorkflowDir, opt: &ExportOptions, path: &Path) -> Result<u64, String> {
    let bytes = export_bytes(w, opt)?;
    write_atomic(path, &bytes).map_err(|e| io_err(path, e))?;
    Ok(bytes.len() as u64)
}

// ---------------------------------------------------------------------------------------------------------- inspect and import

/// What a `.nyx` file contains, read without extracting anything.
#[derive(Clone, Debug, Default)]
pub struct Inspection {
    pub title: String,
    pub summary: String,
    pub tags: Vec<String>,
    pub experiments: usize,
    pub runs: usize,
    pub versions: usize,
    pub level: Option<ShareLevel>,
    pub nyx_version: String,
    pub file_bytes: u64,
    pub unpacked_bytes: u64,
    /// scripts inside (they stay off until approved)
    pub scripts: Vec<String>,
    /// the content hashes all match
    pub intact: bool,
    /// why the file cannot be opened; empty when it can
    pub problems: Vec<String>,
}

fn read_capped<R: Read>(r: &mut R, cap: u64) -> Result<Vec<u8>, String> {
    let mut buf = vec![];
    r.take(cap + 1).read_to_end(&mut buf).map_err(|e| e.to_string())?;
    if buf.len() as u64 > cap {
        return Err(format!("a file is larger than the {cap} byte limit"));
    }
    Ok(buf)
}

/// Checks every entry name, kind and size against `lim` before anything is read or written.
fn vet<R: Read + Seek>(z: &mut zip::ZipArchive<R>, lim: &Limits) -> Result<u64, String> {
    if z.len() > lim.max_entries {
        return Err(format!("{} entries; the limit is {}", z.len(), lim.max_entries));
    }
    let mut seen = std::collections::HashSet::new();
    let mut total = 0u64;
    for i in 0..z.len() {
        let f = z.by_index_raw(i).map_err(|e| e.to_string())?;
        let name = f.name().to_string();
        if f.is_dir() {
            continue;
        }
        check_name(&name)?;
        if !allowed(&name) {
            return Err(format!("'{name}' is not a file a nyx workflow contains"));
        }
        if !seen.insert(name.to_lowercase()) {
            return Err(format!("'{name}' appears twice"));
        }
        if f.unix_mode().is_some_and(|m| m & 0o170000 == 0o120000) || f.encrypted() {
            return Err(format!("'{name}' is a link or encrypted"));
        }
        if f.size() > lim.max_file {
            return Err(format!("'{name}' is {} bytes; the limit per file is {}", f.size(), lim.max_file));
        }
        if f.size() > (1 << 20) && f.size() / f.compressed_size().max(1) > lim.max_ratio {
            return Err(format!("'{name}' is compressed more than {}:1", lim.max_ratio));
        }
        total += f.size();
        if total > lim.max_total {
            return Err(format!("more than {} bytes in all", lim.max_total));
        }
    }
    Ok(total)
}

fn open_zip(path: &Path) -> Result<zip::ZipArchive<std::fs::File>, String> {
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    zip::ZipArchive::new(f).map_err(|e| io_err(path, format!("not a .nyx file: {e}")))
}

pub fn inspect(path: &Path, lim: &Limits) -> Inspection {
    let mut ins = Inspection { file_bytes: std::fs::metadata(path).map(|m| m.len()).unwrap_or(0), ..Default::default() };
    let fail = |ins: &mut Inspection, why: String| {
        ins.problems.push(why);
    };
    let mut z = match open_zip(path) {
        Ok(z) => z,
        Err(e) => {
            fail(&mut ins, e);
            return ins;
        }
    };
    match vet(&mut z, lim) {
        Ok(t) => ins.unpacked_bytes = t,
        Err(e) => {
            fail(&mut ins, e);
            return ins;
        }
    }
    let manifest: ArchiveManifest = match z.by_name(MANIFEST).map_err(|e| e.to_string()).and_then(|mut f| read_capped(&mut f, 16 << 20)).and_then(|b| serde_json::from_slice(&b).map_err(|e| e.to_string())) {
        Ok(m) => m,
        Err(e) => {
            fail(&mut ins, format!("no readable {MANIFEST}: {e}"));
            return ins;
        }
    };
    if manifest.kind != "nyx-archive" || manifest.format > FORMAT {
        fail(&mut ins, format!("made with a newer nyx (format {}); this one reads format {FORMAT}", manifest.format));
        return ins;
    }
    ins.level = Some(manifest.level);
    ins.nyx_version = manifest.nyx_version.clone();
    // every listed file must be present, hash-equal and of the listed size, and nothing unlisted may be there
    let mut intact = true;
    let listed: std::collections::HashSet<&str> = manifest.files.iter().map(|f| f.path.as_str()).collect();
    for i in 0..z.len() {
        let n = z.by_index_raw(i).map(|f| (f.name().to_string(), f.is_dir())).unwrap_or_default();
        if !n.1 && n.0 != MANIFEST && !listed.contains(n.0.as_str()) {
            intact = false;
            ins.problems.push(format!("'{}' is in the file but not in its list", n.0));
        }
    }
    for fe in &manifest.files {
        match z.by_name(&fe.path) {
            Ok(mut f) => match read_capped(&mut f, lim.max_file) {
                Ok(b) if b.len() as u64 == fe.size && hash_bytes(&b) == fe.hash => {
                    if fe.path == "workflow.json" {
                        if let Ok(wf) = serde_json::from_slice::<Workflow>(&b) {
                            ins.title = wf.title.clone();
                            ins.tags = wf.tags.clone();
                            ins.experiments = wf.experiments.len();
                            ins.runs = wf.run_count();
                            ins.summary = wf.summary();
                            if wf.nyx > FORMAT {
                                ins.problems.push(format!("made with a newer nyx (format {})", wf.nyx));
                            }
                        } else {
                            ins.problems.push("workflow.json is not a nyx workflow".into());
                        }
                    }
                    if fe.path == "history.jsonl" {
                        ins.versions = b.split(|c| *c == b'\n').filter(|l| !l.is_empty()).count();
                    }
                    if is_script(&fe.path) {
                        ins.scripts.push(fe.path.clone());
                    }
                }
                _ => {
                    intact = false;
                    ins.problems.push(format!("'{}' does not match its checksum", fe.path));
                }
            },
            Err(_) => {
                intact = false;
                ins.problems.push(format!("'{}' is listed but missing", fe.path));
            }
        }
    }
    if !listed.contains("workflow.json") {
        ins.problems.push("workflow.json is missing".into());
    }
    ins.intact = intact && ins.problems.is_empty();
    ins
}

/// Unpacks `path` as a new workflow in `workspace` (fresh id, scripts off, not pinned, origin "imported").
pub fn import(path: &Path, workspace: &Path, lim: &Limits) -> Result<WorkflowDir, String> {
    let ins = inspect(path, lim);
    if !ins.problems.is_empty() {
        return Err(format!("this .nyx file was not opened: {}", ins.problems.join("; ")));
    }
    let mut z = open_zip(path)?;
    let tmp = workspace.join(format!(".import-{}", short_id()));
    std::fs::create_dir_all(&tmp).map_err(|e| io_err(&tmp, e))?;
    let result = (|| -> Result<WorkflowDir, String> {
        let mut total = 0u64;
        for i in 0..z.len() {
            let mut f = z.by_index(i).map_err(|e| e.to_string())?;
            if f.is_dir() || f.name() == MANIFEST {
                continue;
            }
            let name = f.name().to_string();
            check_name(&name)?;
            let dest = tmp.join(&name);
            if let Some(d) = dest.parent() {
                std::fs::create_dir_all(d).map_err(|e| io_err(d, e))?;
            }
            // counted while decompressing: the header size is not trusted
            let mut out = std::fs::OpenOptions::new().write(true).create_new(true).open(&dest).map_err(|e| io_err(&dest, e))?;
            let mut buf = [0u8; 64 * 1024];
            let mut n_file = 0u64;
            loop {
                let n = f.read(&mut buf).map_err(|e| io_err(&dest, e))?;
                if n == 0 {
                    break;
                }
                n_file += n as u64;
                total += n as u64;
                if n_file > lim.max_file || total > lim.max_total {
                    return Err(format!("'{name}' is larger than its header said; the file was not opened"));
                }
                out.write_all(&buf[..n]).map_err(|e| io_err(&dest, e))?;
            }
        }
        let mut d = WorkflowDir::open(&tmp)?;
        d.wf.id = short_id();
        d.wf.origin = "imported".into();
        d.wf.pinned = false;
        d.wf.scripts_approved = false;
        d.wf.touched = iso(now_secs());
        d.save()?;
        let final_dir = workspace.join(format!("{}_{}", slug(&d.wf.title), d.wf.id));
        std::fs::rename(&tmp, &final_dir).map_err(|e| io_err(&final_dir, e))?;
        d.root = final_dir;
        Ok(d)
    })();
    if result.is_err() {
        let _ = std::fs::remove_dir_all(&tmp);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flow::tests::{spec_json, tmp};
    use crate::runs::RunRecord;

    fn workflow(ws: &Path) -> WorkflowDir {
        let mut d = WorkflowDir::create(ws, "SOI TE modes sweep", "gui", &spec_json(5e-7, 4)).unwrap();
        d.set_tags(vec!["soi".into(), "te modes".into()]).unwrap();
        d.set_notes("width sweep").unwrap();
        let spec = d.read_spec("e1").unwrap();
        for i in 0..3 {
            let r = d.record_run("e1", RunRecord { sec: 0.1, metrics: [("neff_re".to_string(), 2.4 + 0.01 * i as f64)].into(), ..Default::default() }, &spec).unwrap();
            let dir = d.run_dir("e1", &r.id);
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("log.txt"), "stage log").unwrap();
            std::fs::write(dir.join("results.h5"), vec![7u8; 5000]).unwrap();
        }
        std::fs::create_dir_all(d.root.join("analysis")).unwrap();
        std::fs::write(d.root.join("analysis/dispersion.rhai"), "let x = 1;").unwrap();
        d
    }

    fn names(bytes: &[u8]) -> Vec<String> {
        let mut z = zip::ZipArchive::new(std::io::Cursor::new(bytes)).unwrap();
        (0..z.len()).map(|i| z.by_index(i).unwrap().name().to_string()).collect()
    }

    fn opt(level: ShareLevel) -> ExportOptions<'static> {
        ExportOptions { level, only: None, strip_paths: true, exclude: &[] }
    }

    #[test]
    fn export_is_deterministic_and_levels_choose_the_contents() {
        let ws = tmp("exp");
        let d = workflow(&ws);
        let a = export_bytes(&d, &opt(ShareLevel::Recipe)).unwrap();
        // touching file times must not change the bytes
        std::thread::sleep(std::time::Duration::from_millis(30));
        std::fs::write(d.root.join("notes.md"), "width sweep").unwrap();
        let b = export_bytes(&d, &opt(ShareLevel::Recipe)).unwrap();
        assert_eq!(a, b, "same workflow, same bytes");
        let n = names(&a);
        let mut sorted = n.clone();
        sorted.sort();
        assert_eq!(n, sorted, "sorted order");
        assert!(n.iter().any(|x| x == "e1/runs.jsonl") && n.iter().any(|x| x == "workflow.json") && n.iter().any(|x| x == "analysis/dispersion.rhai") && n.iter().any(|x| x.ends_with("/log.txt")));
        assert!(!n.iter().any(|x| x.ends_with(".h5")), "recipe has no heavy files: {n:?}");
        let r = export_bytes(&d, &opt(ShareLevel::Results)).unwrap();
        assert_eq!(names(&r).iter().filter(|x| x.ends_with(".h5")).count(), 3);
        assert!(r.len() > a.len());
        let s = sizes(&d, None).unwrap();
        assert!(s[0].1 < s[1].1 && s[1].1 <= s[2].1);
        // subset: only e1 (the fork e2 is left out)
        let mut d2 = WorkflowDir::open(&d.root).unwrap();
        d2.fork("e1", "second").unwrap();
        let only = vec!["e1".to_string()];
        let sub = export_bytes(&d2, &ExportOptions { level: ShareLevel::Recipe, only: Some(&only), strip_paths: true, exclude: &[] }).unwrap();
        assert!(!names(&sub).iter().any(|x| x.starts_with("e2/")));
        let _ = std::fs::remove_dir_all(&ws);
    }

    #[test]
    fn a_workflow_goes_out_and_comes_back_whole_with_scripts_off() {
        let ws = tmp("rt");
        let d = workflow(&ws);
        let file = ws.join("out.nyx");
        let size = export_file(&d, &opt(ShareLevel::Recipe), &file).unwrap();
        assert_eq!(size, std::fs::metadata(&file).unwrap().len());
        let ins = inspect(&file, &Limits::default());
        assert!(ins.problems.is_empty(), "{:?}", ins.problems);
        assert!(ins.intact);
        assert_eq!((ins.title.as_str(), ins.experiments, ins.runs, ins.versions), ("SOI TE modes sweep", 1, 3, 1));
        assert_eq!(ins.scripts, ["analysis/dispersion.rhai"]);
        assert_eq!(ins.level, Some(ShareLevel::Recipe));

        let inbox = tmp("inbox");
        let mut got = import(&file, &inbox, &Limits::default()).unwrap();
        assert_ne!(got.wf.id, d.wf.id, "a new id, so the original is never overwritten");
        assert_eq!(got.wf.origin, "imported");
        assert!(!got.wf.scripts_approved);
        assert_eq!(got.runs("e1").0.len(), 3);
        assert_eq!(got.read_spec("e1").unwrap(), d.read_spec("e1").unwrap());
        assert_eq!(got.history().len(), d.history().len());
        assert!(got.root.join("analysis/dispersion.rhai").exists());
        assert!(!got.digest_is_stale());
        // and it can be shared on: the second export of the import has the same files
        got.touch().unwrap();
        assert_eq!(names(&export_bytes(&got, &opt(ShareLevel::Recipe)).unwrap()), names(&std::fs::read(&file).unwrap()));
        assert_eq!(std::fs::read_dir(&inbox).unwrap().count(), 1, "no temporary folder is left behind");
        let _ = std::fs::remove_dir_all(&ws);
        let _ = std::fs::remove_dir_all(&inbox);
    }

    #[test]
    fn personal_paths_are_reduced_to_file_names() {
        let mut v = serde_json::json!({"output": {"path": "C:\\Users\\dana\\out\\r.h5"}, "other": "relative/x.h5", "unix": "/home/dana/r.h5"});
        let (out, n) = strip_paths(&v.to_string()).unwrap();
        assert_eq!(n, 2);
        v = serde_json::from_str(&out).unwrap();
        assert_eq!(v["output"]["path"], "r.h5");
        assert_eq!(v["unix"], "r.h5");
        assert_eq!(v["other"], "relative/x.h5");
        assert!(!out.contains("dana"));
    }

    // ---- hostile files

    fn build(entries: &[(&str, Vec<u8>)]) -> Vec<u8> {
        let mut z = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        let o = zip::write::SimpleFileOptions::default().compression_method(zip::CompressionMethod::Deflated);
        for (n, d) in entries {
            z.start_file(*n, o).unwrap();
            z.write_all(d).unwrap();
        }
        z.finish().unwrap().into_inner()
    }

    fn try_open(bytes: &[u8], lim: &Limits) -> Result<WorkflowDir, String> {
        let ws = tmp("hostile");
        let f = ws.join("x.nyx");
        std::fs::write(&f, bytes).unwrap();
        let inbox = ws.join("inbox");
        std::fs::create_dir_all(&inbox).unwrap();
        let r = import(&f, &inbox, lim);
        if r.is_err() {
            assert_eq!(std::fs::read_dir(&inbox).unwrap().count(), 0, "a refused file leaves nothing behind");
        }
        let _ = std::fs::remove_dir_all(&ws);
        r
    }

    #[test]
    fn zip_slip_and_odd_names_are_refused() {
        for bad in ["../evil.txt", "e1/../../evil.txt", "/etc/passwd", "C:/evil.txt", "e1\\runs.jsonl", "e1/runs/r1_a/con.txt", "..", "a//b", "nul.txt", "e1/runs.jsonl "] {
            assert!(check_name(bad).is_err(), "{bad} must be refused");
            let bytes = build(&[(bad, b"x".to_vec())]);
            let err = try_open(&bytes, &Limits::default()).err().unwrap_or_else(|| panic!("{bad} was accepted"));
            assert!(err.contains("unsafe") || err.contains("not a file") || err.contains("not opened"), "{bad}: {err}");
        }
        assert!(check_name("e1/runs/r0001_abc/log.txt").is_ok());
        // names that are fine as paths but not part of a workflow are refused too
        for odd in ["evil.exe", "e1/evil.dll", "specs/../x.nyx.json", "e1/runs/r1_a/x.exe", "workflow.json.bak", "analysis/run.bat"] {
            assert!(!allowed(odd), "{odd}");
        }
        // case-insensitive duplicates would overwrite each other on Windows
        let bytes = build(&[("notes.md", b"a".to_vec()), ("NOTES.md", b"b".to_vec())]);
        assert!(try_open(&bytes, &Limits::default()).is_err());
    }

    #[test]
    fn size_ratio_and_count_caps_hold_and_corruption_is_detected() {
        let ws = tmp("caps");
        let d = workflow(&ws);
        let good = export_bytes(&d, &opt(ShareLevel::Recipe)).unwrap();
        assert!(try_open(&good, &Limits::default()).is_ok());
        assert!(try_open(&good, &Limits { max_entries: 3, ..Default::default() }).unwrap_err().contains("entries"));
        assert!(try_open(&good, &Limits { max_file: 10, ..Default::default() }).unwrap_err().contains("limit per file"));
        assert!(try_open(&good, &Limits { max_total: 100, ..Default::default() }).unwrap_err().contains("in all"));
        // a zip bomb: 8 MiB of zeros is a tiny file; the ratio cap refuses it
        let bomb = build(&[("notes.md", vec![0u8; 8 << 20])]);
        assert!(bomb.len() < 20_000);
        assert!(try_open(&bomb, &Limits { max_ratio: 100, ..Default::default() }).unwrap_err().contains("compressed more than"));
        // a changed byte inside a stored file is not an intact file
        let mut z = zip::ZipArchive::new(std::io::Cursor::new(good.clone())).unwrap();
        let mut entries: Vec<(String, Vec<u8>)> = (0..z.len()).map(|i| { let mut f = z.by_index(i).unwrap(); let mut b = vec![]; f.read_to_end(&mut b).unwrap(); (f.name().to_string(), b) }).collect();
        for e in entries.iter_mut() {
            if e.0 == "notes.md" {
                e.1 = b"tampered".to_vec();
            }
        }
        let refs: Vec<(&str, Vec<u8>)> = entries.iter().map(|(n, b)| (n.as_str(), b.clone())).collect();
        let err = try_open(&build(&refs), &Limits::default()).unwrap_err();
        assert!(err.contains("checksum"), "{err}");
        // not a zip at all
        assert!(try_open(b"hello", &Limits::default()).is_err());
        // a file from a newer format
        let mut m: serde_json::Value = serde_json::from_slice(&entries.iter().find(|e| e.0 == MANIFEST).unwrap().1).unwrap();
        m["format"] = 99.into();
        let refs2: Vec<(&str, Vec<u8>)> = entries.iter().map(|(n, b)| (n.as_str(), if n == MANIFEST { serde_json::to_vec(&m).unwrap() } else { b.clone() })).collect();
        assert!(try_open(&build(&refs2), &Limits::default()).unwrap_err().contains("newer nyx"));
        let _ = std::fs::remove_dir_all(&ws);
    }
}

/// The scripts a share would carry (relative names under `analysis/`).
pub fn scripts_of(w: &WorkflowDir) -> Vec<String> {
    plan(w, ShareLevel::Recipe, None).map(|p| p.into_iter().map(|(n, _)| n).filter(|n| is_script(n)).collect()).unwrap_or_default()
}

#[cfg(test)]
mod exclude_tests {
    use super::*;
    use crate::flow::tests::{spec_json, tmp};

    #[test]
    fn a_script_can_be_left_out() {
        let ws = tmp("excl");
        let d = WorkflowDir::create(&ws, "x", "gui", &spec_json(5e-7, 4)).unwrap();
        std::fs::create_dir_all(d.root.join("analysis")).unwrap();
        std::fs::write(d.root.join("analysis/a.py"), "print(1)").unwrap();
        assert_eq!(scripts_of(&d), ["analysis/a.py"]);
        let ex = vec!["analysis/a.py".to_string()];
        let b = export_bytes(&d, &ExportOptions { level: ShareLevel::Recipe, only: None, strip_paths: true, exclude: &ex }).unwrap();
        let mut z = zip::ZipArchive::new(std::io::Cursor::new(b)).unwrap();
        assert!((0..z.len()).all(|i| !z.by_index(i).unwrap().name().starts_with("analysis/")));
        let _ = std::fs::remove_dir_all(&ws);
    }
}
