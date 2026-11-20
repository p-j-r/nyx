// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Named timing spans for the frame log (productivity plan P3). `NYX_FRAME_LOG=<file>.csv` writes one row per frame; with spans, the row also
//! says where the frame's time went (`poll=48.2;script=0.0;shell_ui=31.1`). A span is a name and the milliseconds a closure took. When the
//! variable is unset a span is just the closure: nothing is read, timed or stored.

use std::cell::RefCell;
use std::sync::OnceLock;
use std::time::Instant;

thread_local! {
    static SPANS: RefCell<Vec<(&'static str, f64)>> = const { RefCell::new(Vec::new()) };
}

fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("NYX_FRAME_LOG").is_ok_and(|p| p.to_lowercase().ends_with(".csv")))
}

/// Runs `f` and, when the frame log is on, adds `name` with the milliseconds it took to this frame's spans (the same name twice in a frame adds up).
pub fn span<T>(name: &'static str, f: impl FnOnce() -> T) -> T {
    if !on() {
        return f();
    }
    let t = Instant::now();
    let v = f();
    let ms = t.elapsed().as_secs_f64() * 1e3;
    SPANS.with(|s| {
        let mut s = s.borrow_mut();
        match s.iter_mut().find(|(n, _)| *n == name) {
            Some(e) => e.1 += ms,
            None => s.push((name, ms)),
        }
    });
    v
}

/// The spans of this frame as `name=ms;name=ms`, and the list is emptied for the next frame.
pub fn take() -> String {
    SPANS.with(|s| s.borrow_mut().drain(..).map(|(n, ms)| format!("{n}={ms:.2}")).collect::<Vec<_>>().join(";"))
}
