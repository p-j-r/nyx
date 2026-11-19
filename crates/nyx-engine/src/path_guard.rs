//! The path guard (productivity plan S4, ADR [89]): a front-end that runs text it did not write (a script) asks the engine to keep every file
//! and folder argument of a command inside one folder. The guard is in the engine, not in the script host, so no language and no front-end has
//! to get it right: `Engine::set_path_root` turns it on and `Engine::call` then rewrites each path argument to an absolute path inside the
//! root, or refuses the command. A relative path is taken from the root, not from the program's working directory.
//!
//! What it stops: `..` out of the root, an absolute path elsewhere, and a link (symbolic link, junction) inside the root that points out of
//! it: the longest part of the path that exists is resolved on disk before it is compared with the root.

use nyx_physics::Error;
use serde_json::{Map, Value};
use std::path::{Component, Path, PathBuf};

/// The names of the arguments that hold a file or a folder, in every command of the catalog. `tests/path_guard.rs` fails when a command takes
/// a path-like argument that is not named here, so a new command cannot slip past the guard.
pub const PATH_ARGS: &[&str] = &["path", "folder", "workspace"];

/// A root made absolute, as the guard keeps it.
pub(crate) fn absolute_root(root: &Path) -> Result<PathBuf, String> {
    let abs = std::path::absolute(root).map_err(|e| format!("{}: {e}", root.display()))?;
    if !abs.is_dir() {
        return Err(format!("{} is not a folder", abs.display()));
    }
    Ok(abs)
}

/// `requested` as an absolute path inside `root`, or the reason it is not allowed.
pub fn confine(root: &Path, requested: &str) -> Result<PathBuf, String> {
    if requested.is_empty() || requested.contains('\0') {
        return Err("the path is empty or holds a NUL character".into());
    }
    let wanted = Path::new(requested);
    let joined = if wanted.is_absolute() { wanted.to_path_buf() } else { root.join(wanted) };
    // lexical normalisation: `.` and `..` resolved without touching the disk
    let mut norm = PathBuf::new();
    for c in joined.components() {
        match c {
            Component::CurDir => {}
            Component::ParentDir => {
                if !norm.pop() {
                    return Err(format!("{requested} goes above the root"));
                }
            }
            other => norm.push(other.as_os_str()),
        }
    }
    // links: resolve the longest part that exists, keep the rest as written
    let real_root = root.canonicalize().map_err(|e| format!("{}: {e}", root.display()))?;
    let mut existing = norm.clone();
    let mut rest: Vec<std::ffi::OsString> = Vec::new();
    while !existing.exists() {
        match existing.file_name() {
            Some(n) => rest.push(n.to_os_string()),
            None => return Err(format!("{requested} has no part that exists")),
        }
        if !existing.pop() {
            break;
        }
    }
    let mut real = existing.canonicalize().map_err(|e| format!("{}: {e}", existing.display()))?;
    for n in rest.iter().rev() {
        real.push(n);
    }
    let inside = real
        .strip_prefix(&real_root)
        .map_err(|_| format!("{requested} is outside the folder scripts may use ({})", root.display()))?;
    Ok(root.join(inside))
}

/// Rewrites the path arguments of one command call, or refuses it.
pub(crate) fn confine_args(root: &Path, command: &str, args: &mut Map<String, Value>) -> Result<(), Error> {
    for key in PATH_ARGS {
        if let Some(Value::String(s)) = args.get(*key) {
            let p = confine(root, s).map_err(|why| Error::invalid(format!("{command}: `{key}`: {why}")))?;
            args.insert((*key).to_string(), Value::String(p.to_string_lossy().into_owned()));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_root(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nyx-path-guard-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join("sub")).unwrap();
        d
    }

    #[test]
    fn a_relative_path_is_taken_from_the_root() {
        let root = temp_root("rel");
        let p = confine(&root, "model.nyx.json").unwrap();
        assert_eq!(p, root.join("model.nyx.json"));
        let p = confine(&root, "sub/../sub/out.h5").unwrap();
        assert_eq!(p, root.join("sub").join("out.h5"));
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn dot_dot_out_of_the_root_is_refused() {
        let root = temp_root("up");
        assert!(confine(&root, "../outside.json").is_err());
        assert!(confine(&root, "sub/../../outside.json").is_err());
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn an_absolute_path_is_allowed_only_inside_the_root() {
        let root = temp_root("abs");
        assert!(confine(&root, root.join("a.json").to_str().unwrap()).is_ok());
        let elsewhere = std::env::temp_dir().join("nyx-elsewhere.json");
        assert!(confine(&root, elsewhere.to_str().unwrap()).is_err());
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn a_sibling_folder_that_starts_with_the_same_name_is_outside() {
        let root = temp_root("sib");
        let sibling = PathBuf::from(format!("{}-evil", root.display()));
        std::fs::create_dir_all(&sibling).unwrap();
        assert!(confine(&root, sibling.join("x.json").to_str().unwrap()).is_err());
        let _ = std::fs::remove_dir_all(&root);
        let _ = std::fs::remove_dir_all(&sibling);
    }

    #[cfg(unix)]
    #[test]
    fn a_link_inside_the_root_that_points_out_is_refused() {
        let root = temp_root("link");
        let out = temp_root("link-target");
        std::os::unix::fs::symlink(&out, root.join("escape")).unwrap();
        assert!(confine(&root, "escape/secret.json").is_err());
        let _ = std::fs::remove_dir_all(&root);
        let _ = std::fs::remove_dir_all(&out);
    }

    #[cfg(windows)]
    #[test]
    fn a_junction_inside_the_root_that_points_out_is_refused() {
        // a directory junction needs no privilege on Windows: `mklink /J`
        let root = temp_root("junction");
        let out = temp_root("junction-target");
        let status = std::process::Command::new("cmd")
            .args(["/C", "mklink", "/J"])
            .arg(root.join("escape"))
            .arg(&out)
            .output()
            .unwrap();
        assert!(status.status.success(), "mklink failed: {}", String::from_utf8_lossy(&status.stdout));
        assert!(confine(&root, "escape/secret.json").is_err());
        let _ = std::fs::remove_dir_all(&root);
        let _ = std::fs::remove_dir_all(&out);
    }
}
