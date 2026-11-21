# Level 3: components of nyx-storage

Crate `crates/nyx-storage`. Pure-Rust HDF5 storage engine (hdf5-pure + zlib-rs): packed 5D Shuffle+LZF and contiguous layouts, bounded-memory streaming reader/writer, direct-path XDMF, and in-place parameter-sweep appender.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_storage.mmd) · [PlantUML](../export/plantuml/structurizr-Components_storage.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-storage"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 86 ["nyx-storage"]
        style 86 fill:#ffffff,stroke:#2e6295,color:#2e6295

        87["<div style='font-weight: bold'>Storage types and read budget</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>FieldStorageLayout (L1<br />Deflate-3, L2 Shuffle+LZF<br />packed 5D, L3 Contiguous<br />packed 5D), ReadBudgetConfig<br />(512 MiB default ceiling),<br />MeshData, ModeHeader,<br />PackedComponentMeta and<br />bitwise f64 verification<br />helpers.</div>"]
        style 87 fill:#85bbf0,stroke:#5d82a8,color:#000000
        88["<div style='font-weight: bold'>Streaming HDF5 reader</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Bounded-memory streaming<br />reader over<br />hdf5_pure::File::open_streaming:<br />checked_mul read-budget<br />guards, O(1) mode-header<br />scan, mesh reader (with<br />optional curved edges), and<br />transparent single-mode<br />cell-field slicing across<br />L2/L3 packed 5D and legacy<br />L1/v3 files.</div>"]
        style 88 fill:#85bbf0,stroke:#5d82a8,color:#000000
        89["<div style='font-weight: bold'>Packed and streaming HDF5 writer</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Whole-file atomic writer<br />(write_nyx_packed) and<br />bounded-RSS incremental<br />mode-by-mode streaming<br />appender<br />(init_nyx_packed_streaming,<br />append_solved_mode_in_place<br />via File::open_rw and<br />append_staged).</div>"]
        style 89 fill:#85bbf0,stroke:#5d82a8,color:#000000
        90["<div style='font-weight: bold'>Direct-path XDMF generator</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Direct-path XDMF 3.0 XML<br />companion generator<br />(generate_direct_xdmf,<br />write_direct_xdmf)<br />referencing internal HDF5<br />dataset paths without /viz<br />hard links.</div>"]
        style 90 fill:#85bbf0,stroke:#5d82a8,color:#000000
        91["<div style='font-weight: bold'>Sweep checkpoint and atomic IO</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Incremental in-place<br />parameter-sweep checkpointing<br />(SweepCheckpoint over<br />unlimited chunked datasets<br />via File::open_rw) and<br />Windows NTFS atomic file<br />replacement helpers.</div>"]
        style 91 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      126[("<div style='font-weight: bold'>Results file</div><div style='font-size: 70%; margin-top: 0px'>[Container: File .h5 and .xdmf]</div><div style='font-size: 80%; margin-top:10px'>Project, mesh and every<br />solution of every solved<br />physics, compressed (packed<br />5D Shuffle+LZF or Deflate),<br />exact; with an XDMF<br />description next to it.</div>")]
      style 126 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    88-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    91-. "<div>uses reader</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    91-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    87-. "<div>uses reader</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    87-. "<div>uses sweep</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    87-. "<div>uses writer</div><div style='font-size: 70%'>[Rust paths]</div>" .->89
    87-. "<div>uses xdmf</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    89-. "<div>uses reader</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    89-. "<div>uses sweep</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    89-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    89-. "<div>uses xdmf</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    90-. "<div>uses sweep</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    90-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    88-. "<div>Streams and slices</div><div style='font-size: 70%'>[declared in spec]</div>" .->126
    89-. "<div>Writes and appends atomically</div><div style='font-size: 70%'>[declared in spec]</div>" .->126
    90-. "<div>Writes .xdmf index</div><div style='font-size: 70%'>[declared in spec]</div>" .->126

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Storage types and read budget](../4-code/storage/types.md) | 1 | 165 | 102, 62, 99 | FieldStorageLayout (L1 Deflate-3, L2 Shuffle+LZF packed 5D, L3 Contiguous packed 5D), ReadBudgetConfig (512 MiB default ceiling), MeshData, ModeHeader, PackedComponentMeta and bitwise f64 verification helpers. |
| [Streaming HDF5 reader](../4-code/storage/reader.md) | 1 | 770 | 102, 62, 99 | Bounded-memory streaming reader over hdf5_pure::File::open_streaming: checked_mul read-budget guards, O(1) mode-header scan, mesh reader (with optional curved edges), and transparent single-mode cell-field slicing across L2/L3 packed 5D and legacy L1/v3 files. |
| [Packed and streaming HDF5 writer](../4-code/storage/writer.md) | 1 | 484 | 102, 62, 99 | Whole-file atomic writer (write_nyx_packed) and bounded-RSS incremental mode-by-mode streaming appender (init_nyx_packed_streaming, append_solved_mode_in_place via File::open_rw and append_staged). |
| [Direct-path XDMF generator](../4-code/storage/xdmf.md) | 1 | 146 | 102, 62, 99 | Direct-path XDMF 3.0 XML companion generator (generate_direct_xdmf, write_direct_xdmf) referencing internal HDF5 dataset paths without /viz hard links. |
| [Sweep checkpoint and atomic IO](../4-code/storage/sweep.md) | 1 | 284 | 102, 87 | Incremental in-place parameter-sweep checkpointing (SweepCheckpoint over unlimited chunked datasets via File::open_rw) and Windows NTFS atomic file replacement helpers. |
