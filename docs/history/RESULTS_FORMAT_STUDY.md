# Results format study: HDF5, Zarr and the rest, and how nyx should ship and use them

Status: study, not a decision. Written 2026-11-20 against `main` 219fe4e1. ADR 62 and 63 (HDF5 is central, paul's decision) stand until paul says otherwise; this document is the evidence for that call. Nothing in the code was changed.

## 1. Summary

1. **Two separate questions.** *How nyx is shipped* (today every exe imports `hdf5.dll` and the dist folder does not carry it) is not the same as *which format nyx writes*. Both a static libhdf5 and a pure-Rust HDF5 crate fix the first and keep the format. Zarr would only change the second.
2. **Recommendation: keep HDF5 as the format.** The use case that decides it is "one portable file you can send to someone". HDF5 is one file; Zarr is a folder of many files that needs packing (and the Rust zip adapter is read-only). Zarr's real advantages (chunk-parallel writes, JSON metadata, no C library) do not match what nyx does today: one writer per run, the file written once and renamed into place.
3. **Fix shipping with a static HDF5 build** (spike A, section 7). It needs cmake and a C compiler on the build machine only, and later GitHub Actions can do it. The pure-Rust HDF5 crate (spike B) is the interesting second route, but it is not a recommendation until it reads real nyx files bit for bit.
4. **The compression codec is not a reason to change format.** On three real files the data is mostly float noise: deflate-3 (what nyx writes now) gets 1.04x to 1.8x on the node fields and zstd-3 gets within 0 to 28% of that. The one big lever is a layout one: in the sms mechanical file, `cell_fields` are 108 of 123 MB raw and look like a five-fold duplicate of the node fields (section 5).
5. **AI handling is a layout question, not a format one.** An agent already reads results through `result_info`, `result_field` and `results_file_info`. What is missing is a small, stable, human- and machine-readable summary next to each `results.h5` (section 6).

## 2. What nyx does with results today (the use cases)

From `crates/nyx-engine/src/{h5io,results_file,worker,shared_queue}.rs` and ADR 62, 63:

| # | Use case | Where it shows |
|---|---|---|
| U1 | Write once at the end of a solve, atomically (a `.partial` file renamed when complete; the queue renames whole item folders) | `write_snapshot`, `shared_queue.rs` |
| U2 | Re-open a file: project, mesh, all fields; re-mesh and check it equals the stored mesh; bit-exact round trip | `open_results`, `tests/h5.rs` |
| U3 | Partial reads: scalars and field info without the bulk, a field sampled to `max_values` | `results_file_info`, `result_info`, `result_field` |
| U4 | Share one file (mail it, attach it, archive it) and have it open on someone else's machine without nyx | ADR 62 context; paul's stated reason for HDF5 |
| U5 | Visualise outside nyx (ParaView, VisIt) | `xdmf.rs` writes `results.xdmf` next to the file |
| U6 | Analyse in Python/MATLAB/Julia | h5py (3.13 is installed here), `h5read` in MATLAB |
| U7 | Run with no extra install | the complaint that started this study |
| U8 | An AI agent reads and reasons about results | MCP tools; see section 6 |
| U9 | Parallel writing (MPI-IO) | `phdf5` feature, **never compiled or run** (no MPI SDK), and every process solves the whole problem anyway: a footnote, not a driver |

Data shape: float64 node fields (`fields/<name>`), per-cell fields stored as re/im pairs (`cell_fields/<name>/re,im`), a mesh (points, cells, facets, labels, curved-edge data), the project as JSON text, many solutions per physics (a mode study has 14 or more). Typical files here are 20 to 90 MB (the biggest in `Documents` are 92 MB).

## 3. The shipping problem, separately

Facts checked:

* `dist/nyx-main-219fe4e1/{nyx,nyx-mcp,nyx-worker}.exe` all contain an import of `hdf5.dll`; the folder has no DLL. They run only when Anaconda's `Library\bin` is on PATH.
* `.cargo/config.toml` points `HDF5_DIR` at `<local HDF5 install>` (an environment value wins, ADR 63).
* `hdf5-metno` 0.15.1 has `default = []`, a feature `static = ["hdf5-sys/static"]` and a feature `zlib = ["hdf5-sys/zlib"]`. In `hdf5-metno-sys` 0.12.5, `static = ["dep:hdf5-src"]` builds HDF5 from source and `zlib = ["dep:libz-sys", "hdf5-src?/zlib"]` is **a separate switch**.
* **Consequence:** the static build must enable **both** `static` and `zlib` (`hdf5-metno = { features = ["static", "zlib"] }`, to be confirmed in the spike). A static HDF5 without zlib cannot read or write the deflate-compressed datasets of every existing `results.h5`, and `tests/h5.rs` would fail. The `threadsafe` feature of `-sys` exists too; the current code does not need it.
* Building HDF5 from source needs cmake and a C compiler on the build machine; neither is on PATH here. Users of the exe need nothing. First build is slow; the result is cached by cargo.

Options for shipping:

| Option | Format kept | Build needs | Run needs | Status |
|---|---|---|---|---|
| a. Today: dynamic libhdf5 | HDF5 | HDF5_DIR | `hdf5.dll` on PATH | what ships now |
| b. Ship the DLLs next to the exes | HDF5 | HDF5_DIR | the DLL files in dist | smallest change, still a C library to track |
| c. **Static libhdf5 (`static`+`zlib`)** | HDF5 | cmake, C compiler | nothing | recommended; spike A |
| d. Pure-Rust HDF5 crate | HDF5 | nothing | nothing | spike B; see 4.3 |
| e. A different format | changes | depends | depends | the rest of this document |

(b) is a fine stopgap while (c) is built.

## 4. Format comparison

### 4.1 Candidates

* **HDF5** (now). One file, 20+ years of tools, h5py/MATLAB/Julia/R readers, deflate built in to every build.
* **Zarr v3**. A hierarchy of JSON metadata plus chunk objects, normally a directory, optionally a zip. Rust: `zarrs` 0.23.14 (Oct 2026), by one maintainer at ANU, "fully up-to-date and conformant with the Zarr 3.1 specification", native `complex128`, sharding, codecs gzip/zstd/blosc/crc32c, filesystem and memory stores, read-only zip adapter (`zarrs_zip`: Read yes, Write no), HTTP/object stores. Its own docs link a list of "correctness issues with past versions". Python: zarr-python 3 (not installed here).
* **VTKHDF**. HDF5 files with a layout ParaView reads natively (unstructured grids, temporal data, compression). Replaces XDMF for visualisation only. The Nov 2025 Kitware status post lists **high-order elements as future work**; I did not verify whether quadratic cells (tri6, tet10) work today, so it cannot be assumed to replace XDMF.
* **Others, one line each.** netCDF-4 / Exodus II / CGNS are HDF5 under another schema and still need the C library. Parquet/Arrow are tabular and a poor fit for meshes. ADIOS2 is an HPC streaming library, not an archive format. SQLite is a single ACID file but no scientific ecosystem reads array blobs from it. A home-grown container (tar + zstd + raw arrays) is easy but nobody else opens it, which loses U4/U6.

### 4.2 HDF5 against Zarr, by use case

| Use case | HDF5 | Zarr v3 |
|---|---|---|
| U1 atomic write once | Write `.partial`, one rename. Works on Windows today. | A directory to rename (many files; on Windows an open handle or a scanner on any one file fails the rename) or a zip written then renamed. No Rust zip *writer* in `zarrs`; pack the folder yourself (stored zip is simple). |
| U2 bit-exact re-open | Tested. | Raw f64 chunks are exact; compression codecs here are lossless. Would need the same test again. |
| U3 partial reads | Chunked datasets; `result_field` already samples. | Equally good; per-chunk reads are native, metadata is plain JSON. |
| U4 share one file | **One file.** | A folder, or a zip you build. Zip store is read by zarr-python and `zarrs_zip` but not written by `zarrs`. |
| U5 ParaView | XDMF works today (quadratic cells included, per ADR 62/README, tested in Rust). VTKHDF would be the ParaView-native layout but high-order is unverified. | ParaView 6 reads Zarr only through a Python xarray plugin with CF-style grids; that does not describe an unstructured quadratic mesh. Nyx would still need to write XDMF or VTKHDF, which are HDF5 formats. |
| U6 Python/MATLAB/Julia | h5py, `h5read`, HDF5.jl, everywhere. | zarr-python, xarray, Zarr.jl. MATLAB support not checked. |
| U7 no install | Needs the static or pure-Rust route (section 3). | Pure Rust apart from C codecs (zstd and blosc are C under the hood; limit to gzip+zstd+crc32c, or pure-Rust gzip only, to avoid cmake). |
| U8 AI | Needs a reader; the JSON summary sidecar (section 6) carries the load. | JSON metadata is readable as text; a minor plus. |
| U9 parallel writers | MPI-IO driver, unbuilt. | Chunks are independent objects, so processes can write different chunks without MPI. Real advantage, but only if nyx ever has several writers of one array. |
| Compression | deflate always; zstd only with a filter plugin on the *reader* side (hdf5plugin for h5py), which undoes portability. | zstd and blosc are registered codecs in the spec. |
| Longevity, risk | HDF Group, huge user base. Corruption on a crash during a non-SWMR write (mitigated by the rename). | Spec 3.0 (2023) and 3.1 are young; one main Rust implementation. |

Result: Zarr wins U9 and, marginally, compression and metadata; HDF5 wins U1, U4, U5, U6 and has the tested code. The deciding case is U4.

### 4.3 Pure-Rust HDF5 (keeps the format, drops the library)

Search found three crates: `hdf5-pure` (0.47.0, commits through 11 Nov 2026, "no C dependencies, checked in CI against libhdf5", chunked storage, deflate/shuffle/fletcher32/LZF, compound types for complex, SWMR, partial reads), `rust-hdf5` and `hdf5-no-c` (claims an API drop-in compatible with the FFI crate). `hdf5-pure`'s own README warns that **other Rust HDF5 projects overlook format rules and "reads can silently return incorrect values"**. That is a reason to test it hard before trusting it, and a reason to read it with its own claim in mind. None of the three is from the HDF Group. Spike B decides.

## 5. Measurements on real files

Script: `docs/results_format_measure.py.txt` (kept as .txt because `*.py` is git-ignored; run with Python, h5py, zstandard; h5py + zlib level 3 + zstd level 3, chunk 32768 elements as nyx writes them). Compression here is of the float64 datasets of the solutions only.

| File | Physics | File size | `cell_fields` raw | deflate-3 | zstd-3 | node `fields` raw | deflate-3 | zstd-3 |
|---|---|---|---|---|---|---|---|---|
| `…guided-mech…/r0001_bd6` | sms eigenfrequency (14 modes) | 69.4 MB | 108.1 MB | 54.2 (2.00x) | 41.7 (2.59x) | 14.8 MB | 13.8 (1.07x) | 13.3 (1.11x) |
| `…optical-eigenmodes…/r0001_zag` | ems modes | 50.9 MB | 39.2 MB | 36.5 (1.07x) | 35.5 (1.10x) | 13.2 MB | 12.7 (1.04x) | 12.6 (1.04x) |
| `p0_3dlarge.h5` | ems cavity (3D) | 26.9 MB | 52.2 MB | 21.3 (2.45x) | 20.6 (2.53x) | 5.1 MB | 4.8 (1.08x) | 4.7 (1.08x) |

Findings:

* **The codec barely matters.** zstd beats deflate by 0 to 28 % on `cell_fields` and 3 to 4 % on node fields. It is not a reason to change format, and the part of that gain that HDF5 could take needs a plugin.
* **Byte-shuffle hurt** on the sms file in an earlier run (1.81x deflate, 1.18x with shuffle), which fits the data being gathered duplicates rather than smooth fields.
* **`cell_fields` dominate** (108 of 123 MB raw in the sms file) and in the sms file they are **about five-fold repeats**: 160 830 values per field but only 32 518 unique (the P2 node count). They look like node values gathered per element. In the ems files they are real distinct values (135 298 of 136 284 unique), and in the 3D cavity the imaginary parts are all zero (1 unique), which already compresses to almost nothing. **Not verified:** that the sms repeats are exactly the gathered node values; check that before removing anything.
* Consequence: for sms, **storing the node field once and gathering per element on read would cut the file by roughly four-fifths of `cell_fields`**, a bigger win than any codec or format change. That is a layout change (results_file layout 4), and it would touch ADR 62's bit-for-bit claim only in that cell fields would be rebuilt rather than read; the values would stay equal.
* Writing, in the same test: deflate-3 3.1 s and zstd-3 1.3 s for 123 MB single-threaded in Python; not a nyx timing.

Not measured: read speed of whole-file open or sampled `result_field` in either format. Needs spike B and a Zarr spike; I did not install `zarr` or fetch `zarrs` for this document.

## 6. An AI handling results: a recommended way of producing them

The format is not what limits an agent. It already has `result_info` (including `stale`, `age_s`, `validated`), `result_field` (sampled), `results_file_info` and `analyse`. Gaps:

1. **No portable summary next to the file.** Someone with only a folder cannot see what the run was without opening HDF5.
2. **Results of a sweep are separate files** with no index.
3. **The guide does not tell an agent how to read results.** (`guide.rs` has no results card.)

Recommendation (a layout and convention, independent of HDF5 vs Zarr):

* Every recorded run writes, next to `results.h5` and `results.xdmf`, a small **`summary.json`** (and optionally a `summary.md` rendering of it): project title and hash, physics and study, mesh size, per solution the label and scalars (frequencies, n_eff, Q, ...), `validated` / convergence state, solver seconds and peak memory, the nyx version and the results layout version, and the paths of the files. Generated from data the engine already holds; no new physics.
* A fixed **reading recipe** in the MCP guide: read `summary.json`, then `result_info`, then `result_field` with `max_values` for slices; open the full file only to re-run analyses.
* For sweeps, one **`index.json`** that lists run folders with their parameters and headline scalars (a table an agent or a person can filter).
* Keep fields out of JSON (ADR 63); the summary holds scalars only.

## 7. Spikes that settle the open points

Ask-first items: each fetches crates or installs tools, so none was run. Use a private `CARGO_TARGET_DIR` (not `../nyx2-target`, which `nyx2-latency` is building into).

**A. Static HDF5.** Install cmake and the MSVC build tools (or confirm the toolchain). In a scratch copy set `hdf5-metno = { features = ["static", "zlib"] }`, unset `HDF5_DIR`, build `nyx`, `nyx-mcp` and `nyx-worker`.
Pass: the `h5` tests pass; the dist exes start and write and re-open a results file with Anaconda removed from PATH; the file has deflate-compressed datasets and h5py reads it back identical to the dynamic build's output. Fail: any dataset unreadable, or the build needs something a CI runner lacks. Then GitHub Actions: a Windows runner with cmake, caching the build.

**B. Pure-Rust HDF5.** In a scratch crate, `hdf5-pure` only.
Pass: it reads three real files (sms, ems, 3D above) and every dataset and attribute equals h5py's reading bit for bit; a file it writes with chunked deflate, VarLenUnicode-equivalent text and u8 JSON datasets reads back in h5py and in `h5io::read_results`; read of one sampled field is as fast as today. Fail: any silent difference. Run the same files through the `hdf5-no-c` crate only if `hdf5-pure` fails on API grounds.

**C. (only if paul wants it) Zarr.** Write the same snapshot with `zarrs` (gzip+zstd+crc32c, sharding, stored-zip pack step), re-open, compare size and times with HDF5. It tests U1 and U4 on Windows with the shared queue, which is where I expect the cost.

**D. Layout 4.** Verify that sms `cell_fields` are gathered node values, then try storing node fields only and rebuild on read; run the benchmark test suite (the reductions read `cell_fields`).

## 8. What would change the recommendation

* If nyx becomes multi-writer or distributed, so that several processes write one result: Zarr (or ADIOS2) becomes attractive and U4 weakens.
* If sharing means "put it on a web server and let a viewer stream a slice": Zarr's HTTP/object-store reads fit better than a monolithic HDF5.
* If spike B passes and the pure-Rust crate stays maintained: no cmake, no DLL, same format. Then shipping becomes trivial and static libhdf5 can be dropped.

## 9. Implementation decision (2026-11-21, ADR 107)

* **Spike B passed (`storage_eval`, 13/13 evaluation gates).** In `<parent folder>/storage_eval`, pure-Rust `hdf5-pure = "=0.48.0"` paired with workspace-unified pure-Rust SIMD `flate2` (`features = ["zlib-rs"]`, `zlib-rs 0.6.8`) passed all 13 verification gates (`EC-00` through `EC-12`):
  * `0 ULP` bitwise read parity on every real `nyx2_benches` `.h5` benchmark (`sos_300nm_modes.h5`, `sos_300nm_rib_modes.h5`, `sos_300nm_mech_modes.h5`, `sos_300nm_strip_coupled_q0_a1.h5`, `swg_3d_bock2010_case1a.h5`).
  * Bidirectional `0 ULP` round-trip compatibility with Python `h5py 3.13.0` across Layout `L1` (`Deflate-3`), Layout `L2` (packed 5D `[M, C, N_c, N_v, 2]` chunked `[1, C, N_c, N_v, 2]` with `Shuffle (2) + LZF (32000)`), and Layout `L3` (contiguous packed 5D).
  * Bounded-RSS single-mode streaming via `File::open_streaming` and `read_f64_rows(k, 1)` (`~57.2 MB` working-set delta when slicing one `60 MB` mode from a `600 MB` 3D `Tetrahedron_10` file), `ReadBudgetConfig` (`512 MiB` default ceiling with `checked_mul` overflow protection), corruption/truncation fuzz resilience, and reproduction of the P2 SOS Brillouin photoelastic and moving-boundary overlap integrals (`B-SBS M1 G_B(v_g) = 37.30 W^-1 m^-1`, `F-SBS A1 G_B(v_g) = 6053.66 W^-1 m^-1`) to `0 ULP`.
* **Integrated into `nyx` via ADR 107.** `hdf5-metno` and `.cargo/config.toml` `HDF5_DIR` are replaced by pure-Rust `hdf5-pure = "=0.48.0"` and `flate2/zlib-rs` in the new leaf crate `crates/nyx-storage` (`reader`, `writer`, `xdmf`, `sweep`) and `crates/nyx-engine/src/h5io.rs`. Layout 4 files in `L2ShuffleLzfPacked` (default) replace per-solution `solution_<k>/cell_fields` with packed 5D `/results/<physics>/cell_fields_packed` (`[M, C, N_c, N_v, R]`, `R = 1` real-only or `R = 2` complex), UTF-8 JSON `u8` dataset `cell_fields_meta`, `eigenvalues`, and attribute `field_layout = "L2ShuffleLzfPacked"`, preserving `solution_<k>/fields` for XDMF and reading legacy v3 and v4 files bit-for-bit.
