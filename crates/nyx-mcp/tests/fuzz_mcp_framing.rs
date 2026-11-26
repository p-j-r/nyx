// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified `bolero` fuzz and regression tests for `nyx_mcp` bounded line framing, JSON-RPC request parsing, channel backlog, and session caps.

use bolero::TypeGenerator;
use nyx_mcp::{
    overload_reply, pump_reader, read_bounded_line, session_name, BoundedLine, Incoming, Server, CHANNEL_CAP,
    MAX_LINE_BYTES, MAX_SESSIONS, MAX_WAIT_S,
};
use serde_json::{json, Value};
use std::io::{BufReader, Cursor};
use std::path::PathBuf;
use std::sync::mpsc;
use std::time::{Duration, Instant};

fn test_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("nyx_mcp_fuzz_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

#[derive(Debug, TypeGenerator)]
struct FramingCase {
    buf_cap: u8,
    max_bytes: u16,
    stream: Vec<u8>,
}

#[derive(Debug, TypeGenerator)]
struct BacklogCase {
    chan_cap: u8,
    max_line: u16,
    stop_on_overload: bool,
    raw_lines: Vec<Vec<u8>>,
}

#[test]
fn fuzz_read_bounded_line_arbitrary_bytes_and_chunking() {
    bolero::check!().with_iterations(512).for_each(|bytes: &[u8]| {
        let cap = bytes.first().copied().unwrap_or(1).max(1) as usize;
        let max_bytes = usize::from(bytes.get(1).copied().unwrap_or(64)).max(1);
        let payload = if bytes.len() > 2 { &bytes[2..] } else { bytes };
        let mut reader = BufReader::with_capacity(cap, Cursor::new(payload));
        let mut count = 0usize;
        loop {
            count += 1;
            assert!(count <= payload.len() + 4, "bounded reader failed to make forward progress");
            match read_bounded_line(&mut reader, max_bytes).expect("in-memory Cursor read") {
                BoundedLine::Eof => break,
                BoundedLine::Line(line) => {
                    assert!(line.len() <= max_bytes);
                    assert!(!line.ends_with('\n') && !line.ends_with('\r'));
                }
                BoundedLine::Bad(_) => {}
                BoundedLine::TooLarge(n) => assert!(n > max_bytes),
            }
        }
    });

    bolero::check!()
        .with_type::<FramingCase>()
        .with_iterations(512)
        .for_each(|case: &FramingCase| {
            let cap = usize::from(case.buf_cap % 64).max(1);
            let max_bytes = usize::from(case.max_bytes).max(1);
            let mut reader = BufReader::with_capacity(cap, Cursor::new(&case.stream));
            loop {
                match read_bounded_line(&mut reader, max_bytes).expect("Cursor read") {
                    BoundedLine::Eof => break,
                    BoundedLine::Line(line) => {
                        assert!(line.len() <= max_bytes);
                        assert!(!line.ends_with('\n') && !line.ends_with('\r'));
                    }
                    BoundedLine::Bad(_) => {}
                    BoundedLine::TooLarge(n) => assert!(n > max_bytes),
                }
            }
        });
}

#[test]
fn fuzz_pump_reader_json_rpc_and_overload_invariants() {
    let dir = test_dir("pump");
    bolero::check!()
        .with_type::<BacklogCase>()
        .with_iterations(256)
        .for_each(|case: &BacklogCase| {
            let cap = usize::from(case.chan_cap % 32).max(1);
            let max_line = usize::from(case.max_line % 4096).max(8);
            let mut stream = Vec::new();
            for line in case.raw_lines.iter().take(64) {
                for &b in line.iter().take(512) {
                    stream.push(if b == b'\n' { b' ' } else { b });
                }
                stream.push(b'\n');
            }

            let (tx, rx) = mpsc::sync_channel(cap);
            let mut overloads = Vec::new();
            pump_reader(&mut Cursor::new(&stream), &tx, max_line, |reply| {
                overloads.push(reply.clone());
                !case.stop_on_overload
            });

            for ov in &overloads {
                assert_eq!(ov["jsonrpc"], json!("2.0"));
                assert_eq!(ov["error"]["code"], json!(-32001));
            }
            drop(tx);

            let (_stx, srx) = mpsc::channel();
            let mut server = Server::with_idle_ttl(
                dir.clone(),
                dir.clone(),
                dir.clone(),
                MAX_WAIT_S,
                srx,
                Duration::from_secs(60),
            );
            let mut delivered = 0usize;
            while let Ok(msg) = rx.try_recv() {
                delivered += 1;
                let safe_msg = match msg {
                    Incoming::Msg(mut v) => {
                        // Avoid spawning nyx-script-worker subprocesses while fuzzing framing.
                        if v.get("method").and_then(Value::as_str) == Some("tools/call") {
                            let tool = v
                                .get("params")
                                .and_then(|p| p.get("name"))
                                .and_then(Value::as_str)
                                .unwrap_or("");
                            if matches!(tool, "script_run" | "script_check") {
                                v["params"]["arguments"] = json!({"source": "1", "path": "x"});
                            }
                        }
                        Incoming::Msg(v)
                    }
                    other => other,
                };
                let _ = server.handle_incoming(safe_msg);
            }
            assert!(delivered <= cap, "sync_channel delivered {delivered} > cap {cap}");
        });
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn fuzz_session_name_backlog_and_session_cap_no_panic() {
    let dir = test_dir("sess");
    bolero::check!().with_iterations(512).for_each(|bytes: &[u8]| {
        let s = String::from_utf8_lossy(bytes);
        let args = json!({ "session": s.as_ref() });
        if let Ok(name) = session_name(&args) {
            assert!(!name.is_empty() && name.len() <= 80);
            assert!(name.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_'));
        }
    });

    // Server push_backlog cap and MAX_SESSIONS / TTL invariants.
    let (_tx, rx) = mpsc::channel();
    let mut server = Server::with_idle_ttl(
        dir.clone(),
        dir.clone(),
        dir.clone(),
        MAX_WAIT_S,
        rx,
        Duration::from_secs(10),
    );
    for i in 0..(CHANNEL_CAP + 32) {
        server.push_backlog(Incoming::Msg(json!({"jsonrpc":"2.0","id":i,"method":"ping"})));
        assert!(server.backlog.len() <= CHANNEL_CAP);
    }

    assert!(overload_reply(&Incoming::Msg(json!({"jsonrpc":"2.0","method":"ping"}))).is_none());
    assert!(overload_reply(&Incoming::Bad("bad".into())).is_some());
    assert!(overload_reply(&Incoming::TooLarge(MAX_LINE_BYTES + 1)).is_some());

    for hostile in ["", ".", "..", "../x", "a/b", "a\\b", "bad name", &"a".repeat(81)] {
        assert!(session_name(&json!({ "session": hostile })).is_err(), "must reject {hostile:?}");
    }
    assert_eq!(session_name(&json!({})).unwrap(), "default");

    for i in 1..MAX_SESSIONS {
        server.ensure_session(&format!("sess_{i}")).expect("under MAX_SESSIONS");
    }
    assert_eq!(server.sessions.len(), MAX_SESSIONS);
    let Err(err) = server.ensure_session("overflow") else {
        panic!("33rd session must fail");
    };
    assert_eq!(err.0, -32002);
    assert!(err.1.contains("SESSION_LIMIT_EXCEEDED"));

    server.reap_idle_sessions(Instant::now() + Duration::from_secs(20));
    assert_eq!(server.sessions.len(), 1);
    assert!(server.sessions.contains_key("default"));
    server.ensure_session("recovered").expect("slot freed after reap");

    let _ = std::fs::remove_dir_all(&dir);
}
