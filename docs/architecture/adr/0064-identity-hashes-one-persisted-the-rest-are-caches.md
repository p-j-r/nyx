# 64. Identity: one persisted hash (BLAKE3 of the bytes); every other hash is a cache key

Date: 2026-11-13

## Status

Accepted (phase 9). Resolves ROUGH_SPOTS R08 ("three hash functions, none documented as the identity rule"). `nyx-session` and its FNV handles are gone (ADR 57), which leaves the two kinds below.

## Context

Hashes are used for two different jobs. One is *identity*: two things with the same hash are the same thing, and the hash is written to disk or compared across runs and machines (a frozen spec, a run's spec, an archive entry). The other is a *cache key*: a hash decides whether a computed value in memory is still valid and never leaves the process.

## Decision

* **Identity is `nyx_store::util::hash_bytes`**: the text `b3:` and the first 32 hex digits of the BLAKE3 hash of the bytes. It is the only hash that is persisted or compared across processes: the content address of a frozen spec (`specs/<hash>.nyx.json`), the spec of a version and of a run (`spec_hash`), the check of every entry of a `nyx` archive. The bytes are the text a project serialises to with sorted keys (`Project::to_json`, ADR 26), so the same project has the same hash in every build. A change of the file format (ADR 63) changes the hashes; that is intended, nothing migrates.
* **Cache keys may use any hash but are never stored or sent.** The engine keys its mesh and its results by `hash_bytes` of the geometry, mesh settings and sizing (`Engine::mesh_key`) and of the whole project (`Engine::stamp`), so an edit never has to remember to invalidate them. `geom::math::Hasher64` (64-bit, over the numeric value of each expression, the operation and the name) decides which features of an evaluated model changed (`geom::diff`, the arrangement cache); its numbers are meaningless outside the process and must never be written to a file or put in a contract. `femcore::adapt::topo::FxHasher` is a hash-map hasher, not an identity.
* **Run and workflow ids are not hashes**: they are short time-ordered random ids (`util::short_id`).
* A new persisted identity uses `hash_bytes`. A new cache key may use what is convenient, and says so in a comment.

## Consequences

* Nothing to unify: `geom` stays free of a hash dependency (R18, `tests/layering.rs`), the store keeps BLAKE3.
* Two projects that differ only in the order of keys or in white space are the same project after `to_json`; two that differ in one digit are different.

## Code

`crates/nyx-store/src/util.rs` (`hash_bytes`, `short_id`), `crates/nyx-engine/src/engine.rs` (`mesh_key`, `stamp`), `crates/geom/src/math.rs` (`Hasher64`), `crates/geom/src/diff.rs`.
