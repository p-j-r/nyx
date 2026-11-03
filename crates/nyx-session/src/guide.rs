//! Cards: short texts a host reads on demand through `get_guide`, embedded in the binary.

pub const TOPICS: &[&str] = &["nyx-instructions", "optical-modes", "troubleshooting"];

pub fn card(topic: &str) -> Option<&'static str> {
    match topic {
        "nyx-instructions" => Some(include_str!("../cards/nyx-instructions.md")),
        "optical-modes" => Some(include_str!("../cards/optical-modes.md")),
        "troubleshooting" => Some(include_str!("../cards/troubleshooting.md")),
        _ => None,
    }
}

/// The text a server sends with `initialize`: just enough to find the guide.
pub const SERVER_INSTRUCTIONS: &str = "nyx runs finite-element photonics (waveguide modes, effective index, polarization). Call get_guide with topic \"nyx-instructions\" before other nyx tools. For a one-line question about a waveguide mode use recipe_run. Use nyx for numbers; state the assumptions listed in its results.";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_topic_has_a_card_within_its_size_budget() {
        for t in TOPICS {
            let c = card(t).unwrap_or_else(|| panic!("no card for {t}"));
            // about 4 characters per token: 700 tokens at most
            assert!(c.len() < 2800, "{t} has {} characters", c.len());
        }
        assert!(card("nope").is_none());
    }
}
