//! Behaviour of nyx-connect: what it does to a host's config file and how it finds and probes the server.
//! Everything runs in a temporary "machine" (see `common`); the real home folder is never read or written.
//! The key order of the user's config file is part of the behaviour (the crate enables `serde_json/preserve_order` itself, ADR 28),
//! which is why these tests may look at key order where the rule of ADR 26 forbids it for data formats.
mod common;
use common::*;
use nyx_connect::*;
use serde_json::{json, Value};
use std::fs;
use std::process::Command;
use std::time::Duration;

fn host<'a>(hs: &'a [Host], id: &str) -> &'a Host {
    hs.iter().find(|h| h.id == id).unwrap()
}

fn read(p: &std::path::Path) -> Value {
    serde_json::from_str(&fs::read_to_string(p).unwrap()).unwrap()
}

fn status_of(env: &Env, id: &str, exe: &std::path::Path) -> Status {
    status(host(&hosts(env), id), exe)
}

#[test]
fn a_host_that_is_not_there_is_not_installed_and_a_present_one_is_not_connected() {
    let sb = Sandbox::new("detect");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    assert!(list(&env, &exe).iter().all(|h| h.status == Status::NotInstalled));
    fs::create_dir_all(env.home.join(".cursor")).unwrap();
    let cursor = list(&env, &exe).into_iter().find(|h| h.id == "cursor").unwrap();
    assert_eq!(cursor.status, Status::NotConnected);
}

#[test]
fn connect_keeps_other_servers_and_settings_in_order_and_backs_up_the_original() {
    let sb = Sandbox::new("merge");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".claude.json");
    fs::create_dir_all(&env.home).unwrap();
    let original = r#"{"numStartups": 7, "mcpServers": {"paper": {"type": "http", "url": "http://x"}}, "theme": "dark"}"#;
    fs::write(&cfg, original).unwrap();
    let done = connect(&env, "claude-code", &exe).unwrap();
    assert!(done.changed);
    let backup = done.backup.expect("a backup of the old file");
    assert_eq!(backup, cfg.with_file_name(".claude.json.nyx-backup"));
    assert_eq!(fs::read_to_string(&backup).unwrap(), original);
    let v = read(&cfg);
    assert_eq!(v["numStartups"], 7);
    assert_eq!(v["theme"], "dark");
    assert_eq!(v["mcpServers"]["paper"]["url"], "http://x");
    assert_eq!(v["mcpServers"]["nyx"], json!({ "type": "stdio", "command": exe.to_string_lossy(), "args": [], "env": {} }));
    let keys: Vec<&String> = v.as_object().unwrap().keys().collect();
    assert_eq!(keys, ["numStartups", "mcpServers", "theme"]);
    let order: Vec<&String> = v["mcpServers"].as_object().unwrap().keys().collect();
    assert_eq!(order, ["paper", "nyx"]);
    assert_eq!(status_of(&env, "claude-code", &exe), Status::Connected);
}

#[test]
fn the_backup_is_always_the_file_as_it_was_before_nyx_first_touched_it() {
    let sb = Sandbox::new("backup");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".cursor/mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    let original = r#"{"mcpServers": {"other": {"command": "o"}}}"#;
    fs::write(&cfg, original).unwrap();
    let backup = cfg.with_file_name("mcp.json.nyx-backup");
    connect(&env, "cursor", &exe).unwrap();
    disconnect(&env, "cursor").unwrap(); // within the same second: a time-stamped name would collide
    connect(&env, "cursor", &exe).unwrap();
    fs::write(&cfg, r#"{"mcpServers": {"nyx": {"command": "/old/nyx-mcp"}}}"#).unwrap();
    connect(&env, "cursor", &exe).unwrap();
    assert_eq!(fs::read_to_string(&backup).unwrap(), original);
    // a connect that changes nothing still reports the backup that exists
    assert_eq!(connect(&env, "cursor", &exe).unwrap().backup, Some(backup));
}

#[test]
fn connecting_twice_changes_nothing_the_second_time() {
    let sb = Sandbox::new("twice");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let first = connect(&env, "cursor", &exe).unwrap();
    assert!(first.changed && first.backup.is_none(), "there was no file to back up");
    let again = connect(&env, "cursor", &exe).unwrap();
    assert!(!again.changed);
}

#[test]
fn an_entry_for_another_exe_is_stale_and_connect_repairs_it() {
    let sb = Sandbox::new("stale");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".cursor").join("mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    fs::write(&cfg, r#"{"mcpServers": {"nyx": {"command": "/old/nyx-mcp", "args": []}}}"#).unwrap();
    assert!(matches!(status_of(&env, "cursor", &exe), Status::Stale { ref found } if found == "/old/nyx-mcp"));
    assert!(connect(&env, "cursor", &exe).unwrap().changed);
    assert_eq!(status_of(&env, "cursor", &exe), Status::Connected);
}

#[test]
fn an_entry_without_a_command_is_stale_with_nothing_found() {
    let sb = Sandbox::new("nocommand");
    let env = sb.env(Os::Linux);
    let cfg = env.home.join(".cursor").join("mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    fs::write(&cfg, r#"{"mcpServers": {"nyx": {"url": "http://x"}}}"#).unwrap();
    assert_eq!(status_of(&env, "cursor", &sb.exe()), Status::Stale { found: String::new() });
}

#[test]
fn keys_the_user_added_to_the_nyx_entry_survive_a_reconnect() {
    let sb = Sandbox::new("userkeys");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".claude.json");
    fs::create_dir_all(&env.home).unwrap();
    fs::write(&cfg, r#"{"mcpServers": {"nyx": {"command": "/old", "env": {"NYX_THREADS": "4"}, "timeout": 30}}}"#).unwrap();
    connect(&env, "claude-code", &exe).unwrap();
    let e = &read(&cfg)["mcpServers"]["nyx"];
    assert_eq!(e["command"], json!(exe.to_string_lossy()));
    assert_eq!(e["env"], json!({ "NYX_THREADS": "4" }), "the user's environment variables are kept");
    assert_eq!(e["timeout"], 30);
    assert_eq!(e["type"], "stdio");
}

#[test]
fn disconnect_removes_only_nyx() {
    let sb = Sandbox::new("remove");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".cursor").join("mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    fs::write(&cfg, r#"{"mcpServers": {"first": {"command": "a"}, "other": {"command": "o"}}}"#).unwrap();
    connect(&env, "cursor", &exe).unwrap();
    assert!(disconnect(&env, "cursor").unwrap().changed);
    let v = read(&cfg);
    assert!(v["mcpServers"].get("nyx").is_none());
    assert_eq!(v["mcpServers"]["other"]["command"], "o");
    let order: Vec<&String> = v["mcpServers"].as_object().unwrap().keys().collect();
    assert_eq!(order, ["first", "other"], "removing an entry does not reorder the rest");
    assert!(!disconnect(&env, "cursor").unwrap().changed);
    assert_eq!(status_of(&env, "cursor", &exe), Status::NotConnected);
}

#[test]
fn disconnect_on_a_missing_file_writes_nothing() {
    let sb = Sandbox::new("noconfig");
    let env = sb.env(Os::Linux);
    let d = disconnect(&env, "cursor").unwrap();
    assert!(!d.changed);
    assert!(!env.home.join(".cursor").exists());
}

#[test]
fn a_file_with_comments_is_never_overwritten() {
    let sb = Sandbox::new("jsonc");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".config/Code/User/mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    let original = "{\n  // my servers\n  \"servers\": {}\n}\n";
    fs::write(&cfg, original).unwrap();
    assert!(matches!(status_of(&env, "vscode", &exe), Status::Unreadable { .. }));
    assert!(connect(&env, "vscode", &exe).is_err());
    assert!(disconnect(&env, "vscode").is_err());
    assert_eq!(fs::read_to_string(&cfg).unwrap(), original);
    let names: Vec<_> = fs::read_dir(cfg.parent().unwrap()).unwrap().map(|e| e.unwrap().file_name()).collect();
    assert_eq!(names.len(), 1, "no backup and no temporary file next to a file that was not changed: {names:?}");
}

#[test]
fn a_config_that_is_not_an_object_or_has_a_non_object_servers_key_is_left_alone() {
    let sb = Sandbox::new("shapes");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".cursor/mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    for text in ["[1, 2]", r#"{"mcpServers": []}"#, r#"{"mcpServers": "x"}"#] {
        fs::write(&cfg, text).unwrap();
        assert!(connect(&env, "cursor", &exe).is_err(), "{text}");
        assert_eq!(fs::read_to_string(&cfg).unwrap(), text);
    }
}

#[test]
fn an_empty_file_and_a_file_with_a_byte_order_mark_are_read() {
    let sb = Sandbox::new("bom");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    let cfg = env.home.join(".cursor/mcp.json");
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    fs::write(&cfg, "  \n").unwrap();
    assert_eq!(status_of(&env, "cursor", &exe), Status::NotConnected);
    assert!(connect(&env, "cursor", &exe).unwrap().changed);
    fs::write(&cfg, "\u{feff}{\"keep\": 1}").unwrap();
    connect(&env, "cursor", &exe).unwrap();
    assert_eq!(read(&cfg)["keep"], 1);
}

#[test]
fn vscode_uses_the_servers_key_and_no_env() {
    let sb = Sandbox::new("vscode");
    let env = sb.env(Os::Mac);
    let exe = sb.exe();
    connect(&env, "vscode", &exe).unwrap();
    let v = read(&env.home.join("Library/Application Support/Code/User/mcp.json"));
    assert_eq!(v["servers"]["nyx"], json!({ "type": "stdio", "command": exe.to_string_lossy(), "args": [] }));
    assert!(v.get("mcpServers").is_none());
}

#[test]
fn claude_desktop_is_written_untyped() {
    let sb = Sandbox::new("desktop");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    connect(&env, "claude-desktop", &exe).unwrap();
    let v = read(&env.home.join(".config/Claude/claude_desktop_config.json"));
    assert_eq!(v["mcpServers"]["nyx"], json!({ "command": exe.to_string_lossy(), "args": [] }));
}

#[test]
fn the_store_build_of_claude_desktop_is_found_on_windows() {
    let sb = Sandbox::new("store");
    let env = sb.env(Os::Windows);
    let pkg = env.local_appdata.clone().unwrap().join("Packages/Claude_pzs8sxrjxfjjc/LocalCache/Roaming/Claude");
    fs::create_dir_all(&pkg).unwrap();
    let hs = hosts(&env);
    assert_eq!(host(&hs, "claude-desktop").config, pkg.join("claude_desktop_config.json"));
}

#[test]
fn without_the_store_package_windows_uses_appdata() {
    let sb = Sandbox::new("appdata");
    let env = sb.env(Os::Windows);
    let hs = hosts(&env);
    assert_eq!(host(&hs, "claude-desktop").config, env.appdata.clone().unwrap().join("Claude/claude_desktop_config.json"));
    assert_eq!(host(&hs, "vscode").config, env.appdata.clone().unwrap().join("Code/User/mcp.json"));
}

#[test]
fn connect_refuses_when_the_server_program_is_missing_or_the_host_is_unknown() {
    let sb = Sandbox::new("noexe");
    let env = sb.env(Os::Linux);
    let e = connect(&env, "cursor", &sb.0.join("missing")).unwrap_err();
    assert!(e.contains("build it first") && !e.contains("cd ai"), "{e}");
    let e = connect(&env, "nope", &sb.0).unwrap_err();
    assert!(e.contains("valid: claude-code, claude-desktop, cursor, vscode, windsurf, antigravity"), "{e}");
    assert!(!env.home.exists(), "a refused connect creates nothing");
}

#[test]
fn every_host_has_its_own_id_and_config_file_and_a_hint() {
    let sb = Sandbox::new("table");
    for os in [Os::Windows, Os::Mac, Os::Linux] {
        let hs = hosts(&sb.env(os));
        let mut ids: Vec<_> = hs.iter().map(|h| h.id).collect();
        ids.dedup();
        assert_eq!(ids.len(), 6);
        let mut files: Vec<_> = hs.iter().map(|h| h.config.clone()).collect();
        files.sort();
        files.dedup();
        assert_eq!(files.len(), 6, "{os:?}");
        assert!(hs.iter().all(|h| !h.after_connect.is_empty() && !h.name.is_empty()));
    }
}

fn touch_server(d: &std::path::Path) {
    let name = if cfg!(windows) { "nyx-mcp.exe" } else { "nyx-mcp" };
    fs::create_dir_all(d).unwrap();
    fs::write(d.join(name), "x").unwrap();
}

fn ends_with(p: &std::path::Path, tail: &str) -> bool {
    p.ends_with(tail.replace('/', std::path::MAIN_SEPARATOR_STR) + if cfg!(windows) { ".exe" } else { "" })
}

#[test]
fn a_shipped_program_looks_only_at_the_override_and_its_own_folder() {
    let sb = Sandbox::new("find_shipped");
    let app = sb.0.join("install/bin/nyx");
    assert!(!is_cargo_build_tree(&app));
    // nothing next to it: the parents and a target/release around it are NOT searched
    touch_server(&sb.0.join("install/target/release"));
    touch_server(&sb.0.join("install"));
    assert_eq!(find_server_with(None, Some(&app)), None, "the caller then offers 'Locate nyx-mcp...'");
    // next to the program
    touch_server(&sb.0.join("install/bin"));
    let found = find_server_with(None, Some(&app)).unwrap();
    assert!(ends_with(&found, "install/bin/nyx-mcp"), "{}", found.display());
    // the override wins, and a wrong override finds nothing instead of falling back
    let elsewhere = sb.0.join("else");
    touch_server(&elsewhere);
    let name = if cfg!(windows) { "nyx-mcp.exe" } else { "nyx-mcp" };
    let found = find_server_with(Some(elsewhere.join(name).into_os_string()), Some(&app)).unwrap();
    assert!(ends_with(&found, "else/nyx-mcp"));
    assert_eq!(find_server_with(Some(sb.0.join("nope").into_os_string()), Some(&app)), None);
    // no `\?\` prefix leaks out (the hosts start the program by this path)
    assert!(!found.to_string_lossy().starts_with(r"\?\"));
    assert_eq!(find_server_with(None, None), None);
}

#[test]
fn a_developer_build_also_looks_through_its_parents_and_their_target_release() {
    let sb = Sandbox::new("find_dev");
    let app = sb.0.join("a/b/release/nyx");
    fs::create_dir_all(sb.0.join("a/b/release/deps")).unwrap();
    assert!(!is_cargo_build_tree(&app), "deps alone is not enough");
    fs::create_dir_all(sb.0.join("a/b/release/.fingerprint")).unwrap();
    assert!(is_cargo_build_tree(&app));
    assert_eq!(find_server_with(None, Some(&app)), None);
    // <root>/target/release, three folders up from the program's own folder
    touch_server(&sb.0.join("a/target/release"));
    assert!(ends_with(&find_server_with(None, Some(&app)).unwrap(), "a/target/release/nyx-mcp"));
    // a parent folder itself (it is looked at before its own target/release)
    touch_server(&sb.0);
    touch_server(&sb.0.join("a"));
    assert!(ends_with(&find_server_with(None, Some(&app)).unwrap(), "a/nyx-mcp"), "the nearer folder wins");
    // the program's own folder wins over everything
    touch_server(&sb.0.join("a/b/release"));
    assert!(ends_with(&find_server_with(None, Some(&app)).unwrap(), "b/release/nyx-mcp"));
    // five folders in total: the program's folder and four parents; a sixth is not looked at
    let far = sb.0.join("far/1/2/3/4/5/release/nyx");
    fs::create_dir_all(far.parent().unwrap().join("deps")).unwrap();
    fs::create_dir_all(far.parent().unwrap().join(".fingerprint")).unwrap();
    touch_server(&sb.0.join("far/1"));
    assert_eq!(find_server_with(None, Some(&far)), None, "far/1 is the sixth folder up from release/ and is not searched");
    touch_server(&sb.0.join("far/1/2"));
    assert!(find_server_with(None, Some(&far)).is_some(), "far/1/2 is the fifth");
}

#[test]
fn plain_path_drops_only_the_verbatim_prefix() {
    assert_eq!(plain_path(r"\\?\C:\x\nyx-mcp.exe".into()).to_string_lossy(), r"C:\x\nyx-mcp.exe");
    assert_eq!(plain_path("/opt/nyx-mcp".into()).to_string_lossy(), "/opt/nyx-mcp");
}

#[test]
fn the_snippet_and_the_claude_code_command_name_the_program() {
    let exe = std::path::Path::new("/opt/nyx-mcp");
    let v: Value = serde_json::from_str(&snippet(exe)).unwrap();
    assert_eq!(v, json!({ "mcpServers": { "nyx": { "command": "/opt/nyx-mcp", "args": [] } } }));
    assert_eq!(claude_code_command(exe), "claude mcp add nyx -- \"/opt/nyx-mcp\"");
}

#[test]
fn antigravity_is_found_by_its_config_folder_not_by_the_folder_gemini_cli_shares() {
    let sb = Sandbox::new("antigravity");
    let env = sb.env(Os::Windows);
    let exe = sb.exe();
    let h = host(&hosts(&env), "antigravity").config.clone();
    assert!(h.ends_with(".gemini/config/mcp_config.json") || h.ends_with(".gemini\\config\\mcp_config.json"), "{}", h.display());
    // Gemini CLI alone leaves `~/.gemini`: that is not Antigravity
    fs::create_dir_all(env.home.join(".gemini")).unwrap();
    assert_eq!(status_of(&env, "antigravity", &exe), Status::NotInstalled);
    fs::create_dir_all(env.home.join(".gemini/config")).unwrap();
    assert_eq!(status_of(&env, "antigravity", &exe), Status::NotConnected);
}

#[test]
fn antigravity_keeps_its_other_servers_and_repairs_an_entry_that_points_at_an_older_build() {
    let sb = Sandbox::new("antigravity_stale");
    let env = sb.env(Os::Windows);
    let exe = sb.exe();
    let cfg = host(&hosts(&env), "antigravity").config.clone();
    fs::create_dir_all(cfg.parent().unwrap()).unwrap();
    // the file as a user had it: nyx set up by hand for an older build, and another server
    fs::write(&cfg, r#"{"mcpServers": {"nyx": {"command": "C:\\old\\nyx-mcp.exe", "args": []}, "other": {"command": "x", "args": ["-y"]}}}"#).unwrap();
    assert!(matches!(status_of(&env, "antigravity", &exe), Status::Stale { .. }));
    let d = connect(&env, "antigravity", &exe).unwrap();
    assert!(d.changed && d.backup.is_some());
    assert_eq!(status_of(&env, "antigravity", &exe), Status::Connected);
    let v = read(&cfg);
    assert_eq!(v["mcpServers"]["nyx"], json!({ "command": exe.to_string_lossy(), "args": [] }), "written untyped, as Cursor and Windsurf are");
    assert_eq!(v["mcpServers"]["other"], json!({ "command": "x", "args": ["-y"] }));
    assert_eq!(v["mcpServers"].as_object().unwrap().keys().collect::<Vec<_>>(), ["nyx", "other"], "the order is kept");
    assert!(disconnect(&env, "antigravity").unwrap().changed);
    assert!(read(&cfg)["mcpServers"].get("other").is_some() && read(&cfg)["mcpServers"].get("nyx").is_none());
}

#[test]
fn the_instructions_for_an_ai_name_the_program_every_file_and_the_rules() {
    let sb = Sandbox::new("instructions");
    let env = sb.env(Os::Linux);
    let exe = sb.exe();
    fs::create_dir_all(env.home.join(".gemini/config")).unwrap();
    let t = agent_instructions(&env, &exe);
    assert!(t.contains(&exe.display().to_string()), "{t}");
    for h in hosts(&env) {
        assert!(t.contains(&h.config.display().to_string()) && t.contains(h.name), "{}: {t}", h.id);
        assert!(t.contains(&format!("under the key \"{}\"", h.servers_key())), "{}: {t}", h.id);
    }
    // the hosts whose entry has more than command and args say so; the plain ones do not repeat it
    assert!(t.contains("\"type\":\"stdio\""), "Claude Code and VS Code want a type: {t}");
    assert!(t.contains("Cursor") && t.contains("(not found on this computer)"), "{t}");
    let antigravity = t.lines().find(|l| l.starts_with("- Antigravity")).unwrap();
    assert!(!antigravity.contains("not found"), "its folder is there: {antigravity}");
    for rule in [".nyx-backup", "plain JSON", "restart", "script_run"] {
        assert!(t.contains(rule), "the text says '{rule}': {t}");
    }
}

fn probe(mode: &str, within: Duration) -> Result<ServerInfo, String> {
    let mut c = Command::new(fake_mcp());
    c.env("FAKE_MCP", mode);
    test_server_with(c, within)
}

#[test]
fn the_handshake_reports_name_version_protocol_and_tools_and_skips_notifications() {
    for mode in ["", "noise"] {
        let info = probe(mode, Duration::from_secs(10)).unwrap();
        assert_eq!(info, ServerInfo { name: "fake-nyx".into(), version: "9.9.9".into(), protocol: "2025-06-18".into(), tools: vec!["geometry_build".into(), "solve".into()] }, "mode '{mode}'");
    }
    assert!(test_server(&fake_mcp()).is_ok());
}

#[test]
fn a_hung_server_a_failing_server_and_a_missing_program_give_plain_errors() {
    let t = std::time::Instant::now();
    let e = probe("silent", Duration::from_millis(400)).unwrap_err();
    assert!(e.contains("no answer to initialize"), "{e}");
    assert!(t.elapsed() < Duration::from_secs(5), "the caller is not held up by a hung server");
    let e = probe("error", Duration::from_secs(10)).unwrap_err();
    assert!(e.contains("tools/list failed") && e.contains("boom"), "{e}");
    let e = test_server(std::path::Path::new("/definitely/not/here/nyx-mcp")).unwrap_err();
    assert!(e.contains("cannot start"), "{e}");
}
