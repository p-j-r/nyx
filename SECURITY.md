# Security Policy

## Supported Versions

Only the latest release on `main` receives security fixes.

| Version | Supported |
| ------- | --------- |
| `main` / latest `v0.x` release | Yes |
| Older releases | No |

## Reporting a Vulnerability

Please **do not** open a public GitHub issue for security vulnerabilities.

Instead, report vulnerabilities privately via GitHub Security Advisories:
https://github.com/p-j-r/nyx/security/advisories/new

Include:
- A description of the issue and affected crate/binary
- Minimal reproduction steps (input file, `.nys` script, or command)
- Expected vs. actual behavior and OS/version

We aim to acknowledge reports within 7 days and coordinate a fix before public disclosure.

## Security Boundaries of `nyx`

Reports are especially relevant for bypasses or crashes in:
- **Script & worker sandboxing** (`nyx-script`, `nyx-engine` `path_guard` — scripts must not read or write outside their workspace root)
- **Untrusted file & archive parsers** (`.gds`, `.fmsh`, `.nyx.json`, `results.h5`, `.nyxpack` / `nyx-setup` extraction)
- **Network fetch guards** (`material_fetch` SSRF and size limits)
- **MCP framing & local IPC** (`nyx-mcp`, `nyx-connect`)
