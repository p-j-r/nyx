# 100. nyx is dual licensed: AGPL-3.0 or a commercial licence

Date: 2026-11-20

## Status

Accepted on Paul's instruction. Built on the bookmark of the merge of the plotting stack with `main`.

## Context

nyx had no licence file. Paul's decision: nyx is licensed under the GNU Affero General Public License, version 3, and a commercial licence is available; people choose which one they want. Parts of the repository came with their own licence: `crates/geom` and the ported solver and mesher modules of `femcore` (BSD-2-Clause, `docs/ports`), the EB Garamond font (SIL OFL 1.1), and the third-party crates.

## Decision

* **The grant is `AGPL-3.0-only`** (version 3, not "or any later version"; Paul chose this), **or** a commercial licence from the copyright holder, at the licensee's choice. The SPDX expression is `AGPL-3.0-only OR LicenseRef-nyx-Commercial`.
* **Notices.** `LICENSE` is the verbatim text of the AGPL-3.0 (661 lines, from gnu.org). `NOTICE` says in one paragraph that the licensee chooses between the two licences, carries the copyright line (`Copyright (C) 2026 p-j-r`), the disclaimer of warranty, and lists what keeps its own licence. The README has a Licence section with the one-line summary. The About sheet of the window (Help, "About nyx and licences") shows the copyright, the choice, the absence of warranty and where the licence text is, which is what section 5(d) of the AGPL asks of an interactive interface.
* **Cargo metadata.** `[workspace.package] license` holds the expression and every crate uses it (`license.workspace = true`), except `geom` (`BSD-2-Clause`, unchanged) and `femcore` (`(AGPL-3.0-only OR LicenseRef-nyx-Commercial) AND BSD-2-Clause`, since it holds the ported modules).
* **Per-file headers.** Every Rust file of nyx's own code starts with `// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial` and `// Copyright (C) 2026 p-j-r` (344 files). Left without a header, because the BSD-2-Clause ports and nyx's own code are mixed there and the file-by-file split is not recorded: `crates/geom` (BSD-2-Clause) and the modules `linalg`, `adapt`, `meshing` and `mesh` of `femcore`; `NOTICE` covers them.
* **Packaging.** `tools/package_dist.sh` builds the release programs and copies them to `dist/<name>/` with `LICENSE`, `NOTICE` and `OFL-EBGaramond.txt`.
* **What keeps its own licence is not relicensed**: the BSD-2-Clause parts, the font (its licence text is shipped next to it and shown in the About sheet) and the third-party crates.

## Not in this decision

* The commercial licence itself (terms, contact) is not in the repository; the notices say it is available from the copyright holder.
* A policy for contributions from other people (dual licensing needs the right to relicense them) is not written.

## Consequences

* `cargo metadata` reports the licence of every crate; a tool that reads SPDX sees the choice.
* The AGPL applies to what is not listed above; section 13 (offering the source to users who interact with a modified version over a network) is the operator's duty for a modified nyx served remotely.

## Code

`LICENSE`, `NOTICE`, `README.md`, `Cargo.toml` and the `Cargo.toml` of the crates, `crates/app/src/sheets.rs` (`sheet_about`).
