# Workflow command lists (productivity plan P0)

JSON-lines files of engine commands, replayed by `crates/nyx-engine/examples/workflow.rs` and checked by `crates/nyx-engine/tests/workflows.rs`.
One line is `{"cmd": "<command of Engine::catalog()>", "args": {...}}`; the optional keys `note`, `job_arg`, `after_stage` and `expect_error`
control the replay only (documented in the header of `workflow.rs`). The same files are the input of `nyx check` / `nyx run` (`nyx check` / `nyx run` take scripts, `*.nys`, not these files; the replay program is the one that reads them).

```
cargo build --release -p nyx-engine --examples
target/release/examples/peakrun target/release/examples/workflow tools/workflows/w1_soi_width.jsonl   # from the repository root; on Windows the programs end in .exe
```

| File | What it times |
|---|---|
| `w1_soi_width.jsonl` | 2D SOI wire: width changed five times, solve, read neff |
| `w2_soi_wavelength.jsonl` | the same wire, ten wavelengths by hand |
| `w3a_cancel_early.jsonl`, `w3b_cancel_factorisation.jsonl`, `w3c_cancel_eigensolver.jsonl` | 3D box cavity (a section swept along y): one cancel per process (so each cancelled solve is the first of its process), early, at the start of the factorisation, 3.5 s later in the iterations; time from cancel to `Cancelled` |
| `b1_blocking_2d.jsonl`, `b2_blocking_3d.jsonl` | the calls a window makes on its own thread: `project_open`, mesh, the `results.h5` write, `open_results` |
| `m1_cavity_large.jsonl` | the larger cavity (about 11 GB, 55 s): memory and results I/O at that size |

`b*` and `m1` write `runs/p0_*.h5` (git-ignored). Wall times drift by up to 1.8x between sessions on this machine: take the best of three and
compare against a control build.
