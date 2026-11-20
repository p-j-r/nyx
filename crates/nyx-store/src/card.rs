// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The concept card: a summary of an experiment generated from its spec, so nothing has to be typed to make it findable.
//!
//! What a project looks like is not the store's business: the engine makes the card (`nyx_engine::store_link::card_of`, ADR 39) and hands the
//! function to `WorkflowDir::create_with` / `open_with`. The store's own default, `plain`, reads only a `title`, so that a spec the store
//! knows nothing about still gets a card and can be found by its title.

use crate::model::{Card, GeoCard};
use serde_json::Value;

/// The card of a spec the store cannot read: its title, if it has one.
pub fn plain(spec: &Value) -> Card {
    let title = spec.get("title").and_then(Value::as_str).unwrap_or("untitled");
    Card { summary: title.to_string(), keywords: title.split_whitespace().map(|w| w.to_lowercase()).collect(), geometry: GeoCard::default(), ..Card::default() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn a_plain_card_carries_the_title_and_never_fails() {
        let c = plain(&json!({ "title": "Strip on sapphire", "anything": [1, 2] }));
        assert_eq!(c.summary, "Strip on sapphire");
        assert!(c.keywords.contains(&"sapphire".to_string()));
        assert_eq!(plain(&json!(3)).summary, "untitled");
    }
}
