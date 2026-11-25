// SPDX-License-Identifier: BSD-2-Clause
// Copyright (C) 2026 p-j-r
//! Unified `bolero` fuzz and regression tests for `geom::gds` parser, flattened importer, path expansion, and `REC_BOX` / resource-cap limits.

use bolero::TypeGenerator;
use geom::gds::{
    encode_ibm_f64, expand_path_polygons, import_gds_bytes_with_limits, import_gds_library_with_limits,
    parse_gds_with_limits, GdsAref, GdsBoundary, GdsCell, GdsImportLimits, GdsImportOptions, GdsImportResult,
    GdsLibrary, GdsPath, GdsSref, GdsText,
};
use geom::math::V2;

const TEST_LIMITS: GdsImportLimits = GdsImportLimits {
    max_file_bytes: 64 * 1024,
    max_records: 2_000,
    max_depth: 12,
    max_instances: 256,
    max_polygons: 256,
    max_decompose_vertices: 256,
};

fn default_opts() -> GdsImportOptions {
    GdsImportOptions {
        layers: vec![("m1".into(), (1, 0)), ("m2".into(), (2, 0))],
        cell: None,
        split_islands: true,
        pin_layer: Some((10, 0)),
        simplify_tol: 0.0,
        arc_tol: 0.0,
    }
}

fn assert_result_finite(res: &GdsImportResult) {
    assert!(res.dbu_m.is_finite() && res.dbu_m > 0.0);
    for sh in &res.sheets {
        assert!(sh.area.is_finite(), "non-finite sheet area: {}", sh.area);
        assert!(sh.centroid[0].is_finite() && sh.centroid[1].is_finite());
        assert!(sh.bbox.iter().all(|v| v.is_finite()));
        for ring in sh.outers.iter().chain(&sh.holes) {
            for p in ring {
                assert!(p[0].is_finite() && p[1].is_finite(), "non-finite ring vertex: {p:?}");
            }
        }
        for (outer, holes) in &sh.islands {
            for p in outer.iter().chain(holes.iter().flatten()) {
                assert!(p[0].is_finite() && p[1].is_finite(), "non-finite island vertex: {p:?}");
            }
        }
    }
}

fn assert_rings_finite(rings: &[Vec<V2>]) {
    for ring in rings {
        for p in ring {
            assert!(p[0].is_finite() && p[1].is_finite(), "non-finite path polygon vertex: {p:?}");
        }
    }
}

fn push_record(out: &mut Vec<u8>, tag: u16, payload: &[u8]) {
    let mut len = 4usize + payload.len();
    let pad = len % 2;
    len += pad;
    let Ok(len_u16) = u16::try_from(len) else {
        return;
    };
    out.extend_from_slice(&len_u16.to_be_bytes());
    out.extend_from_slice(&tag.to_be_bytes());
    out.extend_from_slice(payload);
    if pad == 1 {
        out.push(0);
    }
}

fn push_i16_record(out: &mut Vec<u8>, tag: u16, vals: &[i16]) {
    let mut payload = Vec::with_capacity(vals.len() * 2);
    for &v in vals {
        payload.extend_from_slice(&v.to_be_bytes());
    }
    push_record(out, tag, &payload);
}

fn push_i32_record(out: &mut Vec<u8>, tag: u16, vals: &[i32]) {
    let mut payload = Vec::with_capacity(vals.len() * 4);
    for &v in vals {
        payload.extend_from_slice(&v.to_be_bytes());
    }
    push_record(out, tag, &payload);
}

fn push_str_record(out: &mut Vec<u8>, tag: u16, s: &str) {
    let mut b = s.as_bytes().to_vec();
    if b.len() % 2 == 1 {
        b.push(0);
    }
    push_record(out, tag, &b);
}

fn wrap_gds_stream(body: &[u8]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(128 + body.len());
    push_i16_record(&mut buf, 0x0002, &[600]); // HEADER
    push_i16_record(&mut buf, 0x0102, &[2026, 1, 1, 0, 0, 0, 2026, 1, 1, 0, 0, 0]); // BGNLIB
    push_str_record(&mut buf, 0x0206, "FUZZLIB"); // LIBNAME
    let mut units = [0u8; 16];
    units[..8].copy_from_slice(&encode_ibm_f64(1e-3));
    units[8..].copy_from_slice(&encode_ibm_f64(1e-9));
    push_record(&mut buf, 0x0305, &units); // UNITS
    buf.extend_from_slice(body);
    push_record(&mut buf, 0x0400, &[]); // ENDLIB
    buf
}

fn structured_cell_body(
    edges_tag: u16,
    sub_tag: u16,
    width: i32,
    pathtype: i16,
    bgnextn: i32,
    endextn: i32,
    xy: &[i32],
) -> Vec<u8> {
    let mut body = Vec::new();
    push_i16_record(&mut body, 0x0502, &[2026, 1, 1, 0, 0, 0, 2026, 1, 1, 0, 0, 0]); // BGNSTR
    push_str_record(&mut body, 0x0606, "TOP"); // STRNAME
    push_record(&mut body, edges_tag, &[]);
    push_i16_record(&mut body, 0x0D02, &[1]); // LAYER
    push_i16_record(&mut body, sub_tag, &[0]); // DATATYPE or BOXTYPE
    if edges_tag == 0x0900 {
        push_i16_record(&mut body, 0x2102, &[pathtype]); // PATHTYPE
        push_i32_record(&mut body, 0x0F03, &[width]); // WIDTH
        push_i32_record(&mut body, 0x3003, &[bgnextn]); // BGNEXTN
        push_i32_record(&mut body, 0x3103, &[endextn]); // ENDEXTN
    }
    push_i32_record(&mut body, 0x1003, xy); // XY
    push_record(&mut body, 0x1100, &[]); // ENDEL
    push_record(&mut body, 0x0700, &[]); // ENDSTR
    body
}

#[derive(Debug, TypeGenerator)]
struct FuzzGdsCase {
    pathtype: i16,
    width: i32,
    bgnextn: i32,
    endextn: i32,
    layer: i16,
    cols: i16,
    rows: i16,
    split_islands: bool,
    with_roi: bool,
    simplify_scale: u8,
    arc_scale: u8,
    use_box_record: bool,
    pts: [[i16; 2]; 5],
}

fn scaled_pt(p: [i16; 2]) -> [i32; 2] {
    [i32::from(p[0]) * 10, i32::from(p[1]) * 10]
}

#[test]
fn fuzz_gds_raw_bytes_no_panic() {
    let opts = default_opts();
    let mut seq_buf = Vec::new();
    bolero::check!().with_iterations(512).for_each(|bytes: &[u8]| {
        if let Ok(lib) = parse_gds_with_limits(bytes, &TEST_LIMITS) {
            let _ = import_gds_library_with_limits(&lib, &opts, &TEST_LIMITS, None);
        }
        if let Ok(res) = import_gds_bytes_with_limits(bytes, &opts, &TEST_LIMITS, Some([-1.0, -1.0, 1.0, 1.0])) {
            assert_result_finite(&res);
        }

        // Also feed bytes inside a minimal GDSII header+cell wrapper so deep record parsers run often.
        seq_buf.clear();
        let mut cell = Vec::new();
        push_i16_record(&mut cell, 0x0502, &[2026, 1, 1, 0, 0, 0, 2026, 1, 1, 0, 0, 0]);
        push_str_record(&mut cell, 0x0606, "TOP");
        cell.extend_from_slice(&bytes[..bytes.len().min(4096)]);
        push_record(&mut cell, 0x0700, &[]);
        seq_buf.extend_from_slice(&wrap_gds_stream(&cell));
        if let Ok(res) = import_gds_bytes_with_limits(&seq_buf, &opts, &TEST_LIMITS, None) {
            assert_result_finite(&res);
        }
    });
}

#[test]
fn fuzz_gds_structured_paths_rec_box_and_library_no_panic() {
    bolero::check!().with_type::<FuzzGdsCase>().with_iterations(512).for_each(|case: &FuzzGdsCase| {
        // 1. Direct path polygon expansion edge cases (pathtype 4, negative/zero width, coincident pts).
        let path_v2: Vec<V2> = case.pts.iter().map(|p| [f64::from(p[0]) * 1e-6, f64::from(p[1]) * 1e-6]).collect();
        for pt in [case.pathtype, 0, 1, 2, 4] {
            let w = if case.width == i32::MIN { -1e-6 } else { f64::from(case.width) * 1e-9 };
            let rings =
                expand_path_polygons(&path_v2, w, pt, f64::from(case.bgnextn) * 1e-9, f64::from(case.endextn) * 1e-9);
            assert_rings_finite(&rings);
        }

        // 2. Structured GdsLibrary covering boundary, path (pathtype 4 & negative width), SREF, AREF, TEXT.
        let mut ring: Vec<[i32; 2]> = case.pts[..4].iter().copied().map(scaled_pt).collect();
        ring.push(ring[0]);
        let leaf = GdsCell {
            name: "LEAF".into(),
            boundaries: vec![GdsBoundary { layer: case.layer.clamp(1, 2), datatype: 0, xy: ring }],
            paths: vec![GdsPath {
                layer: 1,
                datatype: 0,
                pathtype: if case.pathtype % 2 == 0 { 4 } else { case.pathtype },
                width: case.width.clamp(-10_000, 10_000),
                bgnextn: case.bgnextn.clamp(-2_000, 2_000),
                endextn: case.endextn.clamp(-2_000, 2_000),
                xy: case.pts.iter().copied().map(scaled_pt).collect(),
            }],
            srefs: Vec::new(),
            arefs: Vec::new(),
            texts: vec![GdsText {
                layer: 10,
                texttype: 0,
                string: "pin_0".into(),
                xy: scaled_pt(case.pts[0]),
                strans: 0,
                mag: 1.0,
                angle_deg: 0.0,
            }],
        };
        let top = GdsCell {
            name: "TOP".into(),
            boundaries: Vec::new(),
            paths: Vec::new(),
            srefs: vec![GdsSref { sname: "LEAF".into(), strans: 0, mag: 1.0, angle_deg: 0.0, xy: [0, 0] }],
            arefs: vec![GdsAref {
                sname: "LEAF".into(),
                strans: 0,
                mag: 1.0,
                angle_deg: 0.0,
                cols: case.cols.clamp(-4, 8),
                rows: case.rows.clamp(-4, 8),
                xy: [[0, 0], [100_000, 0], [0, 100_000]],
            }],
            texts: Vec::new(),
        };
        let lib =
            GdsLibrary { version: 600, name: "STRUCT".into(), user_unit: 1e-3, dbu_m: 1e-9, cells: vec![top, leaf] };
        let opts = GdsImportOptions {
            layers: vec![("m1".into(), (1, 0)), ("m2".into(), (2, 0))],
            cell: Some("TOP".into()),
            split_islands: case.split_islands,
            pin_layer: Some((10, 0)),
            simplify_tol: f64::from(case.simplify_scale) * 1e-8,
            arc_tol: f64::from(case.arc_scale) * 1e-8,
        };
        let roi = case.with_roi.then_some([-5e-4, -5e-4, 5e-4, 5e-4]);
        if let Ok(res) = import_gds_library_with_limits(&lib, &opts, &TEST_LIMITS, roi) {
            assert_result_finite(&res);
        }

        // 3. Binary REC_BOX (0x2D00) / REC_PATH (0x0900) round-trip through parse + import.
        let mut xy_flat = Vec::with_capacity(10);
        for p in case.pts[..4].iter().copied().map(scaled_pt) {
            xy_flat.push(p[0]);
            xy_flat.push(p[1]);
        }
        xy_flat.push(xy_flat[0]);
        xy_flat.push(xy_flat[1]);
        let body = if case.use_box_record {
            structured_cell_body(0x2D00, 0x2E02, 0, 0, 0, 0, &xy_flat)
        } else {
            structured_cell_body(
                0x0900,
                0x0E02,
                case.width.clamp(-5_000, 5_000),
                if case.pathtype % 2 == 0 { 4 } else { case.pathtype },
                case.bgnextn.clamp(-1_000, 1_000),
                case.endextn.clamp(-1_000, 1_000),
                &xy_flat,
            )
        };
        let stream = wrap_gds_stream(&body);
        if let Ok(res) = import_gds_bytes_with_limits(&stream, &default_opts(), &TEST_LIMITS, None) {
            assert_result_finite(&res);
        }
    });
}

#[test]
fn regression_rec_box_pathtype4_negative_width_and_limits() {
    // 1. REC_BOX imports as polygon with expected area.
    let box_body = structured_cell_body(0x2D00, 0x2E02, 0, 0, 0, 0, &[0, 0, 1000, 0, 1000, 2000, 0, 2000, 0, 0]);
    let box_bytes = wrap_gds_stream(&box_body);
    let box_res = import_gds_bytes_with_limits(&box_bytes, &default_opts(), &TEST_LIMITS, None).expect("REC_BOX valid");
    assert_eq!(box_res.sheets.len(), 1);
    assert!((box_res.sheets[0].area - 2e-12).abs() < 1e-18);

    // 2. PATH with pathtype == 4 + asymmetric extensions + negative WIDTH parses & imports cleanly without panic.
    let path_body = structured_cell_body(0x0900, 0x0E02, -400, 4, 150, 250, &[0, 0, 2000, 0]);
    let path_bytes = wrap_gds_stream(&path_body);
    let path_lib = parse_gds_with_limits(&path_bytes, &TEST_LIMITS).expect("pathtype 4 + negative width parses");
    assert_eq!(path_lib.cells[0].paths[0].pathtype, 4);
    assert_eq!(path_lib.cells[0].paths[0].width, -400);
    let path_res =
        import_gds_library_with_limits(&path_lib, &default_opts(), &TEST_LIMITS, None).expect("pathtype 4 imports");
    assert_result_finite(&path_res);
    assert!(expand_path_polygons(&[[0.0, 0.0], [1.0, 0.0]], -1.0, 4, 0.1, 0.2).is_empty());

    // 3. Keyhole/high-vertex polygon exceeds max_decompose_vertices and returns Err instead of OOM.
    let n_pts = 128i32;
    let mut poly = Vec::with_capacity((n_pts as usize) * 2 + 2);
    for i in 0..n_pts {
        let th = std::f64::consts::TAU * (i as f64) / (n_pts as f64);
        poly.push((10_000.0 * th.cos()).round() as i32);
        poly.push((10_000.0 * th.sin()).round() as i32);
    }
    poly.push(poly[0]);
    poly.push(poly[1]);
    let dense_body = structured_cell_body(0x0800, 0x0E02, 0, 0, 0, 0, &poly);
    let dense_bytes = wrap_gds_stream(&dense_body);
    let mut tight = TEST_LIMITS;
    tight.max_decompose_vertices = 32;
    let res = import_gds_bytes_with_limits(&dense_bytes, &default_opts(), &tight, None);
    assert!(res.is_err(), "expected vertex-budget error when max_decompose_vertices exceeded");

    // 4. Stream size and record limits reject before allocation.
    let mut tiny = TEST_LIMITS;
    tiny.max_file_bytes = 16;
    assert!(parse_gds_with_limits(&box_bytes, &tiny).is_err());

    let mut few_recs = TEST_LIMITS;
    few_recs.max_records = 3;
    assert!(parse_gds_with_limits(&box_bytes, &few_recs).is_err());
}
