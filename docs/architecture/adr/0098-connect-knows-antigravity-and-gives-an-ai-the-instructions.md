# 98. Connect knows Antigravity and gives an AI the instructions

Date: 2026-11-20

## Status

Accepted as the latency work L4 (Paul, 2026-11-20: "why cant we have both; that way someone can type into antigravity and it can find the right line; but also people can press connect"). Extends ADR 28 (the one copy of the Connect logic). The number was free after ADR 97.

## Context

Antigravity, where the benchmark session of ADR 97 ran, was not a host of `nyx-connect`: the Connect menu of the app listed only installed hosts, so on a machine with Antigravity and none of the five others it offered "Copy config JSON" and nothing else. The `nyx` entry in Antigravity's configuration had been written by hand and pointed at an older build (the one without the reader and the wait limit of ADR 97), which Connect could not have told the user.

On a machine that runs Antigravity the configuration is `~/.gemini/config/mcp_config.json`, an object `{"mcpServers": {"nyx": {"command", "args"}}}`: the shape Cursor and Windsurf use. Its folder `~/.gemini` is also the folder of Gemini CLI, so it says nothing about Antigravity. A configuration file in the workspace (`.agents/mcp_config.json`) is documented by the vendor and reported as ignored; it is not used.

The command-line program, started by a double click, opened a console that closed at once: with no argument it printed its help and ended.

## Decision

* **A host row.** `antigravity`, name `Antigravity`, file `~/.gemini/config/mcp_config.json` on every operating system, key `mcpServers`, the plain entry (`command`, `args: []`, no `type`). Found by its `config` folder (`~/.gemini/config`), not by `~/.gemini`. Connect, Reconnect (an entry that points at another build is `stale`), Disconnect, the backup and the JSON report work as for the other hosts. `connect_hosts` and `connect_report` go to 1.1 (a host and an enum value added).
* **Instructions an AI can follow.** `nyx_connect::agent_instructions(env, exe)`, shown by `nyx-connect instructions` and by the button "Copy instructions for an AI" of the Connect menu: a text for any AI app that can edit files. It names the program, the entry, the file and the key of every host (this machine's paths, hosts not found marked) and the rules `connect` itself follows: keep the other servers, back the file up once, keep it plain JSON, change only `command` of an existing entry, then say that the app must be restarted. Built from the host table, so it cannot disagree with the rows.
* **A double click is answered.** `nyx-connect` without an argument prints its commands and the state of every host, and on Windows, when it is alone in its console (`GetConsoleProcessList`), waits for Enter. It writes nothing.
* **No `--max-wait` per host.** The entry stays `"args": []`: the default of 90 s is under the 180 s of Antigravity. A host with a shorter limit gets `--max-wait` by hand (`docs/HOST_SETUP.md`).

## Consequences

* Three ways to the same entry: the button, the instructions pasted into the AI app, the snippet by hand. They write the same thing; the row is the only code that edits a file.
* Not tried on a real Antigravity: the file and the entry shape come from a machine that runs it, whether it reloads the file without a restart does not. The text after connecting says "restart or refresh".
* A fifth kind of host (a different shape, a file per workspace) is still one row and one test; nothing in the logic knows the names.
