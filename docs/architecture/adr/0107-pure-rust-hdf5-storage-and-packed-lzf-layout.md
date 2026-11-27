# 107. Pure-Rust HDF5 storage (`hdf5-pure` + `zlib-rs`) and packed 5D Shuffle+LZF layout

Date: 2026-11-21

## Status

Accepted. Amends ADR 13 (pure-Rust build and HDF5 dependencies), ADR 26 (workspace crate layering), ADR 62 and 63 (`results.h5` writer and `HDF5_DIR`), and ADR 99 (results file layout version 4 and direct-path XDMF). Resolves `docs/RESULTS_FORMAT_STUDY.md` Spike B via the `storage_eval` benchmark and verification harness (13/13 evaluation gates passed).

## Context

`docs/RESULTS_FORMAT_STUDY.md` identified two coupled issues with how `nyx` stored and shipped HDF5 results (`crates/nyx-engine/src/h5io.rs`):

1. **Shipping and build fragility (`hdf5.dll` / `HDF5_DIR`).** Every binary (`nyx.exe`, `nyx-mcp.exe`, `nyx-worker.exe`) linked against `hdf5-metno` (`0.15.1` / `hdf5-metno-sys 0.12.5`), requiring `HDF5_DIR` at build time and `hdf5.dll` on `PATH` at run time (or a C compiler and CMake for a static `hdf5-src` + `libz-sys` build).
2. **Dataset fragmentation and memory scaling on multi-mode 3D solves.** Legacy layout `L1` (`nyx-results` v3/v4 per-solution groups) created $O(16 \times M)$ separate B-tree headers and `Deflate(3)` datasets per file (`/results/<physics>/solution_<k>/cell_fields/<comp>/{re,im}`). On 14-mode 3D quadratic (`Tetrahedron_10`) mechanical and electromagnetic solves (`69–600 MB`), whole-file reading buffered the entire file in memory, `Deflate(3)` achieved only $1.07\times$ to $2.45\times$ compression on `cell_fields`, and parameter sweeps had no in-place append mechanism.

In the `storage_eval` evaluation (`<parent folder>/storage_eval`), pure-Rust `hdf5-pure = "=0.48.0"` paired with workspace-unified pure-Rust SIMD `flate2` (`features = ["zlib-rs"]`, `zlib-rs 0.6.8`) passed all 13 evaluation gates (`EC-00` through `EC-12`), including `0 ULP` bitwise parity on every real `nyx2_benches` `.h5` benchmark (`sos_300nm_modes.h5`, `sos_300nm_rib_modes.h5`, `sos_300nm_mech_modes.h5`, `sos_300nm_strip_coupled_q0_a1.h5`, `swg_3d_bock2010_case1a.h5`), bidirectional `0 ULP` round-trip compatibility with Python `h5py 3.13.0`, bounded-RSS single-mode streaming (`~57.2 MB` working-set delta when slicing one `60 MB` mode from a `600 MB` 3D `Tet10` file), corruption/truncation fuzz resilience, and reproduction of the P2 SOS Brillouin photoelastic (`PE`) and moving-boundary (`MB`) overlap integrals (`B-SBS M1 G_B(v_g) = 37.30 W^-1 m^-1`, `F-SBS A1 G_B(v_g) = 6053.66 W^-1 m^-1`) to `0 ULP`.

## Decision

1. **Zero-C pure-Rust HDF5 dependency (`hdf5-pure = "=0.48.0"` + `zlib-rs`).**
   * `hdf5-metno` and `.cargo/config.toml` `HDF5_DIR` are removed from the workspace.
   * The workspace root `Cargo.toml` pins `hdf5-pure = { version = "=0.48.0", default-features = true }` (without `fast-deflate`, which would pull C `zlib-ng` and CMake) and unifies `flate2 = { version = "1.1", default-features = false, features = ["zlib-rs"] }` across all crates.
   * Compilation succeeds 100% offline with `rustc` alone (`zero-C`, `zero-CMake`, `zero-*-sys`), verified mechanically with poisoned `CC`, `CXX`, and `CMAKE` environment variables (`Gate EC-00`).
   * The `phdf5` cargo feature on `nyx-engine` is retained as a zero-C compatibility stub: `H5Mode::Parallel` returns an explicit unsupported error because `hdf5-pure` has no MPI-IO driver (`RESULTS_FORMAT_STUDY.md` U9 noted `phdf5` was never compiled or run).

2. **New leaf crate `crates/nyx-storage`.**
   * `nyx-storage` is added to the workspace as a leaf storage engine (`nyx-engine -> nyx-storage` in `crates/nyx-engine/tests/layering.rs`; `nyx-storage` depends on no other workspace crate).
   * Its five modules provide:
     * `lib.rs`: `FieldStorageLayout`, `ReadBudgetConfig`, `MeshData` (including optional curved-edge arrays `edge_node_pairs` `[N_e, 2]` and `edge_node_points` `[N_e, 3]`), `ModeHeader`, `PackedComponentMeta`, and `assert_f64_bits_eq`.
     * `reader.rs` (`NyxReader`): opens files via `hdf5_pure::File::open_streaming(path)` (positional `pread`, never buffering the whole file in RAM on open); enforces `ReadBudgetConfig` (default `max_dataset_bytes = 512 MiB`, `max_group_depth = 24`) using `checked_mul` across dataset dimensions before allocation (`Gate EC-11`); decodes both scalar and 1-element array string attributes (`as_str()` and `as_strings()`); supports both optical `eig_re`/`eig_im` and mechanical `value` (`eig_im = 0.0`) eigenvalue attributes; and transparently slices a single mode from either packed 5D datasets (`ds.read_f64_rows(mode_idx as u64, 1)`) or legacy per-solution `cell_fields/<comp>/{re,im}` datasets.
     * `writer.rs`: writes whole files (`write_nyx_packed`) with atomic `<stem>.h5.tmp` to `<stem>.h5` rename (`Gate EC-10`), and supports bounded-RSS incremental mode-by-mode streaming (`init_nyx_packed_streaming` + `append_solved_mode_in_place`) using `MaxExtent::Unlimited` on mode axis `0` and `File::open_rw` + `append_staged` + `session.commit()`.
     * `xdmf.rs`: `generate_direct_xdmf` and atomic `write_direct_xdmf`, emitting XDMF 3.0 XML referencing HDF5 dataset paths directly without `/viz` hard links.
     * `sweep.rs`: `SweepCheckpoint` for incremental in-place parameter-sweep checkpointing (`File::open_rw` + `append` / `append_staged` + `commit()`) and Windows NTFS atomic file replacement helpers.

3. **Three field storage layouts (`FieldStorageLayout`) with `L2ShuffleLzfPacked` as default.**
   * `L2ShuffleLzfPacked` (default): packs all $M$ modes of a physics into `/results/<physics>/cell_fields_packed` of shape `[M, C, N_c, N_v, R]` (`f64`, where `R = 1` when all cell-field components are real-only (`!has_im`) and `R = 2` (`[re, im]`) when any component is complex), chunked per mode as `[1, C, N_c, N_v, R]` with filter pipeline `Shuffle (id = 2) -> LZF (id = 32000)`, accompanied by `/results/<physics>/eigenvalues` (`f64` `[M, 2]`), UTF-8 JSON `u8` dataset `/results/<physics>/cell_fields_meta` encoding `Vec<PackedComponentMeta { name, unit, nloc, has_im }>`, and attribute `field_layout = "L2ShuffleLzfPacked"` (`FieldStorageLayout::as_str()`). Readable out of the box by vanilla Python `h5py` (`f["/results/ems/cell_fields_packed"][k]`) with zero external filter plugins, achieving $3.2\times$ smaller files ($16.3\times$ compression ratio on 3D `Tet10`) and $4.2\times$ faster writes than legacy `Deflate(3)`.
   * `L3ContiguousPacked`: same 5D `[M, C, N_c, N_v, R]` dataset stored contiguous and uncompressed (`field_layout = "L3ContiguousPacked"`) for scratch/hot-loop single-`pread` mode slicing ($1.98\text{--}3.20\ \text{GB/s}$).
   * `L1LegacyGzip3`: per-solution `/results/<physics>/solution_<k>/cell_fields/<comp>/{re,im}` datasets chunked with `Shuffle + Deflate(3)` (`field_layout = "L1LegacyGzip3"`) for structural backward compatibility.
   * In `L2ShuffleLzfPacked` and `L3ContiguousPacked`, `nyx-engine::h5io` writes `/results/<physics>/cell_fields_packed` and `/results/<physics>/cell_fields_meta` in place of the per-solution `solution_<k>/cell_fields` groups (replacing rather than duplicating them), while preserving per-solution `solution_<k>/fields` and `scalars` for XDMF and metadata queries. On read, `h5io` reconstructs `cell_fields` from `cell_fields_packed` (accepting `R = 1` for real-only and `R = 2` for complex) when present and falls back to per-solution `solution_<k>/cell_fields` on `L1LegacyGzip3` and legacy v3/v4 files.

## Consequences

* Building and shipping `nyx`, `nyx-mcp`, and `nyx-worker` require no C compiler, no CMake, no `HDF5_DIR`, and no `hdf5.dll`.
* All existing layout-3 and layout-4 `.h5` files remain readable bit-for-bit (`0 ULP`), while newly written results files carry the packed L2 `Shuffle + LZF` representation and can be sliced mode-by-mode with bounded memory (`ReadBudgetConfig`).
* `phdf5` parallel MPI-IO writing is retired in favour of single-file atomic writes and in-place `SweepCheckpoint` / streaming mode appends.

## Code

`crates/nyx-storage/Cargo.toml`, `crates/nyx-storage/src/{lib,reader,writer,xdmf,sweep}.rs`, `crates/nyx-storage/tests/gates.rs`, `crates/nyx-engine/src/{h5io,results_file,shared_queue}.rs`, `crates/nyx-engine/tests/{h5,layering}.rs`, `Cargo.toml`, `.cargo/config.toml`.

## Contract

`results_file` layout version 4 (additive packed `L2` datasets `cell_fields_packed`, `cell_fields_meta`, `eigenvalues`, and attribute `field_layout` in `docs/architecture/contracts/VERSIONS.md`); 5 mandatory verification gates (`EC-00`, `EC-01`, `EC-02`/`EC-06`, `EC-09`, `EC-11`/`EC-12`) in `crates/nyx-storage/tests/gates.rs`.
