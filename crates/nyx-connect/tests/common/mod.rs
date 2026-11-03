//! Helpers shared by the test files of nyx-connect: a throw-away "machine" (home folder, program) and the fake MCP server.
#![allow(dead_code)]
use nyx_connect::{Env, Os};
use std::fs;
use std::path::PathBuf;

/// A temporary folder that is removed when it goes out of scope. Nothing in the tests reads the real home folder.
pub struct Sandbox(pub PathBuf);

impl Sandbox {
    pub fn new(tag: &str) -> Sandbox {
        let p = std::env::temp_dir().join(format!("nyx_connect_{tag}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Sandbox(p)
    }

    /// An `Env` whose home is inside the sandbox.
    pub fn env(&self, os: Os) -> Env {
        Env { home: self.0.join("home"), appdata: Some(self.0.join("home/AppData/Roaming")), local_appdata: Some(self.0.join("home/AppData/Local")), os }
    }

    /// A file that stands in for the `nyx-mcp` program (it only has to exist).
    pub fn exe(&self) -> PathBuf {
        let p = self.0.join(if cfg!(windows) { "nyx-mcp.exe" } else { "nyx-mcp" });
        fs::write(&p, "x").unwrap();
        p
    }
}

impl Drop for Sandbox {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The built `fake_mcp` example (`cargo test` builds the examples of the package before running its tests).
pub fn fake_mcp() -> PathBuf {
    let me = std::env::current_exe().unwrap();
    let dir = me.parent().unwrap().parent().unwrap().join("examples");
    let p = dir.join(if cfg!(windows) { "fake_mcp.exe" } else { "fake_mcp" });
    assert!(p.is_file(), "the fake_mcp example is not built at {}", p.display());
    p
}
