// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cards: short texts a host reads on demand through the `get_guide` command, embedded in the binary. They describe the commands, not a
//! physics: how to build a geometry, assign materials, solve and read a result.

pub const TOPICS: &[&str] = &["nyx-instructions", "geometry", "optical-modes", "cqed", "troubleshooting"];

pub fn card(topic: &str) -> Option<&'static str> {
    match topic {
        "nyx-instructions" => Some(include_str!("../cards/nyx-instructions.md")),
        "geometry" => Some(include_str!("../cards/geometry.md")),
        "optical-modes" => Some(include_str!("../cards/optical-modes.md")),
        "cqed" => Some(include_str!("../cards/cqed.md")),
        "troubleshooting" => Some(include_str!("../cards/troubleshooting.md")),
        _ => None,
    }
}

/// The text a server sends with `initialize`: just enough to find the guide (the `api` tool) and to know that the work is scripts.
pub const SERVER_INSTRUCTIONS: &str = "nyx runs finite-element photonics and phononics (waveguide modes, effective index, polarization, vibration modes, heat). You work by writing short Python scripts: call api first (it is the guide and the list of commands), then script_run; script_check checks a script without running it. Do a whole task in one script_run, not one call per step; keep wait_s under your tool timeout (at most 90) and continue a `running` answer with script_wait. Only what a script prints or returns comes back. Use nyx for numbers and state the assumptions you used.";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_topic_has_a_card_within_its_size_budget() {
        for t in TOPICS {
            let c = card(t).unwrap_or_else(|| panic!("no card for {t}"));
            // about 4 characters per token: 800 tokens at most
            assert!(c.len() < 3600, "{t} has {} characters", c.len());
        }
        assert!(card("nope").is_none());
    }
}
