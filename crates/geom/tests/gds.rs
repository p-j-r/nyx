// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for the pure-Rust binary GDSII parser, IBM float codec, hierarchy
//! flattener (`SREF`/`AREF` with reflection, magnification, and rotation), `PATH` expansion,
//! keyhole & overlapping boolean polygon union, hull/hole orientation, island ordering, and
//! `TEXT` pin matching.

use geom::gds::{
    decode_ibm_f64, encode_ibm_f64, import_gds_bytes, import_gds_bytes_with_limits, import_gds_library,
    import_gds_library_with_limits, parse_gds, parse_gds_with_limits, GdsCell, GdsImportLimits, GdsImportOptions,
    GdsLibrary, GdsSref, STRANS_REFLECT_X,
};
use geom::prim::polygon_area;

/// Helper builder for synthesizing spec-compliant binary GDSII streams in tests.
struct GdsWriter {
    buf: Vec<u8>,
}

impl GdsWriter {
    fn new() -> Self {
        Self { buf: Vec::new() }
    }

    fn rec(&mut self, rec_type: u8, data_type: u8, payload: &[u8]) {
        let len = (4 + payload.len()) as u16;
        self.buf.extend_from_slice(&len.to_be_bytes());
        self.buf.push(rec_type);
        self.buf.push(data_type);
        self.buf.extend_from_slice(payload);
    }

    fn rec_empty(&mut self, rec_type: u8) {
        self.rec(rec_type, 0x00, &[]);
    }

    fn rec_i16(&mut self, rec_type: u8, val: i16) {
        self.rec(rec_type, 0x02, &val.to_be_bytes());
    }

    fn rec_u16_bits(&mut self, rec_type: u8, val: u16) {
        self.rec(rec_type, 0x01, &val.to_be_bytes());
    }

    fn rec_i32(&mut self, rec_type: u8, val: i32) {
        self.rec(rec_type, 0x03, &val.to_be_bytes());
    }

    fn rec_real8(&mut self, rec_type: u8, val: f64) {
        self.rec(rec_type, 0x05, &encode_ibm_f64(val));
    }

    fn rec_ascii(&mut self, rec_type: u8, s: &str) {
        let mut p = s.as_bytes().to_vec();
        if !p.len().is_multiple_of(2) {
            p.push(0);
        }
        self.rec(rec_type, 0x06, &p);
    }

    fn rec_xy(&mut self, pts: &[[i32; 2]]) {
        if pts.is_empty() {
            self.rec(0x10, 0x03, &[]);
            return;
        }
        for chunk_pts in pts.chunks(8000) {
            let mut p = Vec::with_capacity(chunk_pts.len() * 8);
            for &[x, y] in chunk_pts {
                p.extend_from_slice(&x.to_be_bytes());
                p.extend_from_slice(&y.to_be_bytes());
            }
            self.rec(0x10, 0x03, &p);
        }
    }

    fn begin_lib(&mut self, name: &str, user_unit: f64, dbu_m: f64) {
        self.rec_i16(0x00, 600); // HEADER v6.0
        self.rec(0x01, 0x02, &[0u8; 24]); // BGNLIB
        self.rec_ascii(0x02, name); // LIBNAME
        let mut u = Vec::with_capacity(16);
        u.extend_from_slice(&encode_ibm_f64(user_unit));
        u.extend_from_slice(&encode_ibm_f64(dbu_m));
        self.rec(0x03, 0x05, &u); // UNITS
    }

    fn end_lib(&mut self) {
        self.rec_empty(0x04); // ENDLIB
    }

    fn begin_cell(&mut self, name: &str) {
        self.rec(0x05, 0x02, &[0u8; 24]); // BGNSTR
        self.rec_ascii(0x06, name); // STRNAME
    }

    fn end_cell(&mut self) {
        self.rec_empty(0x07); // ENDSTR
    }

    fn boundary(&mut self, layer: i16, datatype: i16, open_pts: &[[i32; 2]]) {
        self.rec_empty(0x08); // BOUNDARY
        self.rec_i16(0x0D, layer);
        self.rec_i16(0x0E, datatype);
        let mut closed = open_pts.to_vec();
        if closed.first() != closed.last() {
            closed.push(open_pts[0]);
        }
        self.rec_xy(&closed);
        self.rec_empty(0x11); // ENDEL
    }

    fn path(&mut self, layer: i16, datatype: i16, pathtype: i16, width: i32, pts: &[[i32; 2]]) {
        self.rec_empty(0x09); // PATH
        self.rec_i16(0x0D, layer);
        self.rec_i16(0x0E, datatype);
        self.rec_i16(0x21, pathtype);
        self.rec_i32(0x0F, width);
        self.rec_xy(pts);
        self.rec_empty(0x11); // ENDEL
    }

    fn sref(&mut self, sname: &str, xy: [i32; 2], reflect_x: bool, mag: Option<f64>, angle_deg: Option<f64>) {
        self.rec_empty(0x0A); // SREF
        self.rec_ascii(0x12, sname);
        if reflect_x || mag.is_some() || angle_deg.is_some() {
            self.rec_u16_bits(0x1A, if reflect_x { STRANS_REFLECT_X } else { 0 });
        }
        if let Some(m) = mag {
            self.rec_real8(0x1B, m);
        }
        if let Some(a) = angle_deg {
            self.rec_real8(0x1C, a);
        }
        self.rec_xy(&[xy]);
        self.rec_empty(0x11); // ENDEL
    }

    #[allow(clippy::too_many_arguments)]
    fn aref(
        &mut self,
        sname: &str,
        cols: i16,
        rows: i16,
        xy3: [[i32; 2]; 3],
        reflect_x: bool,
        mag: Option<f64>,
        angle_deg: Option<f64>,
    ) {
        self.rec_empty(0x0B); // AREF
        self.rec_ascii(0x12, sname);
        if reflect_x || mag.is_some() || angle_deg.is_some() {
            self.rec_u16_bits(0x1A, if reflect_x { STRANS_REFLECT_X } else { 0 });
        }
        if let Some(m) = mag {
            self.rec_real8(0x1B, m);
        }
        if let Some(a) = angle_deg {
            self.rec_real8(0x1C, a);
        }
        let mut cr = [0u8; 4];
        cr[0..2].copy_from_slice(&cols.to_be_bytes());
        cr[2..4].copy_from_slice(&rows.to_be_bytes());
        self.rec(0x13, 0x02, &cr); // COLROW
        self.rec_xy(&xy3);
        self.rec_empty(0x11); // ENDEL
    }

    fn text(&mut self, layer: i16, texttype: i16, xy: [i32; 2], label: &str) {
        self.rec_empty(0x0C); // TEXT
        self.rec_i16(0x0D, layer);
        self.rec_i16(0x16, texttype);
        self.rec_xy(&[xy]);
        self.rec_ascii(0x19, label);
        self.rec_empty(0x11); // ENDEL
    }
}

#[test]
fn ibm_float_codec_matches_known_vectors_and_roundtrips_bit_for_bit() {
    assert_eq!(decode_ibm_f64([0, 0, 0, 0, 0, 0, 0, 0]), 0.0);
    // 1.0 = (1/16) * 16^1 => exp = 65 = 0x41, mantissa = 0x10_0000_0000_0000
    assert_eq!(decode_ibm_f64([0x41, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]), 1.0);
    // -45.0 = -(45/256) * 16^2 => sign|exp = 0x80|0x42 = 0xC2, mantissa = 0x2D_0000_0000_0000
    assert_eq!(decode_ibm_f64([0xC2, 0x2D, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]), -45.0);
    // Real GDSII 1e-9 m (1 nm) DBU record from KLayout: 0x39 44 B8 2F A0 9B 5A 54
    let dbu_nm = decode_ibm_f64([0x39, 0x44, 0xB8, 0x2F, 0xA0, 0x9B, 0x5A, 0x54]);
    assert!((dbu_nm - 1.0e-9).abs() < 1.0e-24, "dbu_nm = {dbu_nm:e}");

    for &v in &[0.0, 1.0, -1.0, 0.5, -2.5, 90.0, -135.0, 1.0e-3, 1.0e-6, 1.0e-9, 2.5e-10, std::f64::consts::PI] {
        let enc = encode_ibm_f64(v);
        let dec = decode_ibm_f64(enc);
        assert_eq!(dec.to_bits(), v.to_bits(), "IBM float roundtrip mismatch for {v}");
    }
}

#[test]
fn overlapping_and_keyhole_polygons_extract_ccw_hulls_and_cw_holes_with_exact_si_coords() {
    // 1 nm DBU (1e-9 m)
    let mut w = GdsWriter::new();
    w.begin_lib("hulls_holes_lib", 1.0e-3, 1.0e-9);
    w.begin_cell("TOP");

    // Layer (1, 0): 4 overlapping rectangular frame bars forming a [-100..100] x [-100..100] um outer square
    // with an interior [-40..40] x [-30..30] um hole!
    w.boundary(1, 0, &[[-100_000, -100_000], [100_000, -100_000], [100_000, -30_000], [-100_000, -30_000]]);
    w.boundary(1, 0, &[[-100_000, 30_000], [100_000, 30_000], [100_000, 100_000], [-100_000, 100_000]]);
    w.boundary(1, 0, &[[-100_000, -100_000], [-40_000, -100_000], [-40_000, 100_000], [-100_000, 100_000]]);
    w.boundary(1, 0, &[[40_000, -100_000], [100_000, -100_000], [100_000, 100_000], [40_000, 100_000]]);

    // Layer (2, 0): KLayout-style single-BOUNDARY keyhole polygon:
    // Outer [0..80] x [0..60] um with CW keyhole cut to inner hole [20..60] x [15..45] um.
    w.boundary(
        2,
        0,
        &[
            [80_000, 60_000],
            [0, 60_000],
            [0, 15_000],      // cutline start A
            [20_000, 15_000], // cutline end B
            [20_000, 45_000], // CW hole walk
            [60_000, 45_000],
            [60_000, 15_000],
            [20_000, 15_000], // back to B
            [0, 15_000],      // reverse cutline B -> A
            [0, 0],
            [80_000, 0],
        ],
    );

    w.end_cell();
    w.end_lib();
    // Append 2048-byte zero tape padding after ENDLIB to verify trailing pad tolerance.
    w.buf.extend_from_slice(&[0u8; 128]);

    let opts = GdsImportOptions {
        layers: vec![("frame".into(), (1, 0)), ("keyhole".into(), (2, 0))],
        cell: None,
        split_islands: true,
        pin_layer: None,
        simplify_tol: 0.0,
        ..Default::default()
    };
    let res = import_gds_bytes(&w.buf, &opts).unwrap();
    assert_eq!(res.cell, "TOP");
    assert_eq!(res.sheets.len(), 2);

    // Check "frame": outer = 200x200 um = 40_000 um^2, hole = 80x60 um = 4_800 um^2 => net = 35_200 um^2
    let frame = &res.sheets[0];
    assert_eq!(frame.name, "frame");
    assert_eq!(frame.outers.len(), 1);
    assert_eq!(frame.holes.len(), 1);
    assert_eq!(frame.outers[0].len(), 4);
    assert_eq!(frame.holes[0].len(), 4);
    assert!(polygon_area(&frame.outers[0]) > 0.0, "exterior hull must be CCW");
    assert!(polygon_area(&frame.holes[0]) < 0.0, "interior hole must be CW");
    assert!((polygon_area(&frame.outers[0]) - 40_000e-12).abs() < 1e-22);
    assert!((polygon_area(&frame.holes[0]) + 4_800e-12).abs() < 1e-22);
    assert!((frame.area - 35_200e-12).abs() < 1e-22);
    assert_eq!(frame.bbox, [-100e-6, -100e-6, 100e-6, 100e-6]);

    // Check "keyhole": outer = 80x60 um = 4_800 um^2, hole = 40x30 um = 1_200 um^2 => net = 3_600 um^2
    let kh = &res.sheets[1];
    assert_eq!(kh.name, "keyhole");
    assert_eq!(kh.outers.len(), 1);
    assert_eq!(kh.holes.len(), 1);
    assert_eq!(kh.outers[0].len(), 4, "cutline seam vertex must be stripped cleanly: {:?}", kh.outers[0]);
    assert_eq!(kh.holes[0].len(), 4);
    assert!(polygon_area(&kh.outers[0]) > 0.0, "keyhole exterior hull must be CCW");
    assert!(polygon_area(&kh.holes[0]) < 0.0, "keyhole interior hole must be CW");
    assert!((kh.area - 3_600e-12).abs() < 1e-22);
    assert_eq!(kh.bbox, [0.0, 0.0, 80e-6, 60e-6]);
    assert_eq!(res.layer_bbox, [-100e-6, -100e-6, 100e-6, 100e-6]);
}

#[test]
fn sref_aref_hierarchy_reflection_mag_rotation_paths_and_pins() {
    let mut w = GdsWriter::new();
    w.begin_lib("hier_lib", 1.0e-3, 1.0e-9);

    // Unreferenced dummy macro with no requested layers (tests smart root-cell auto-detection!)
    w.begin_cell("UNUSED_LOGO");
    w.boundary(99, 0, &[[0, 0], [5_000, 0], [5_000, 5_000], [0, 5_000]]);
    w.end_cell();

    // Child cell "L_PAD": asymmetric L-shape in [0..20] x [0..10] um with a pin on (1, 20) at (5, 5) um.
    // Area = 20*10 - 10*5 = 150 um^2.
    w.begin_cell("L_PAD");
    w.boundary(1, 10, &[[0, 0], [20_000, 0], [20_000, 5_000], [10_000, 5_000], [10_000, 10_000], [0, 10_000]]);
    w.end_cell();

    // Child cell "UNIT_DOT": 6 um x 4 um rectangle [0..6] x [0..4] um (area = 24 um^2)
    w.begin_cell("UNIT_DOT");
    w.boundary(1, 10, &[[0, 0], [6_000, 0], [6_000, 4_000], [0, 4_000]]);
    w.end_cell();

    // Root cell "CHIP_TOP"
    w.begin_cell("CHIP_TOP");

    // 1. SREF of "L_PAD" with reflect_x = true, mag = 2.0, angle_deg = 90.0 at origin (100 um, 50 um).
    // Local point (x, y):
    //   reflect_x => (x, -y)
    //   mag = 2   => (2x, -2y)
    //   rot +90   => (2y, 2x)
    //   + (100, 50) => (100 + 2y, 50 + 2x) um!
    // So local [0..20] x [0..10] maps to world [100..120] x [50..90] um, area = 150 * 4 = 600 um^2!
    w.sref("L_PAD", [100_000, 50_000], true, Some(2.0), Some(90.0));
    w.text(1, 20, [105_000, 60_000], "Q1_pad");

    // 2. AREF of "UNIT_DOT" (cols = 3, rows = 2) rotated 90 deg CCW at (-100 um, 0 um):
    // Unrotated UNIT_DOT [0..6] x [0..4] rotated +90 deg => [-4..0] x [0..6] um.
    // Col step = (0, 20 um) => p1 = (-100, 60) um; Row step = (-20 um, 0) => p2 = (-140, 0) um.
    // Produces 6 disjoint 24 um^2 islands!
    w.aref("UNIT_DOT", 3, 2, [[-100_000, 0], [-100_000, 60_000], [-140_000, 0]], false, None, Some(90.0));

    // 3. PATH elements on layer (3, 0):
    // - Flush path (pathtype = 0): (0, 0) -> (50 um, 0), width = 10 um => area = 500 um^2, x in [0, 50], y in [-5, 5]
    w.path(3, 0, 0, 10_000, &[[0, 0], [50_000, 0]]);
    // - Half-width extend path (pathtype = 2): (0, 30 um) -> (40 um, 30 um), width = 10 um
    //   extends by 5 um on each end => x in [-5, 45] (len 50 um), y in [25, 35], area = 500 um^2
    w.path(3, 0, 2, 10_000, &[[0, 30_000], [40_000, 30_000]]);
    // - Round-cap path (pathtype = 1): (0, 60 um) -> (30 um, 60 um), width = 10 um
    //   area ~= 30*10 + pi * 5^2 = 378.54 um^2
    w.path(3, 0, 1, 10_000, &[[0, 60_000], [30_000, 60_000]]);

    w.end_cell();
    w.end_lib();

    let lib = parse_gds(&w.buf).unwrap();
    assert_eq!(lib.cells.len(), 4);

    let opts = GdsImportOptions {
        layers: vec![("sig".into(), (1, 10)), ("traces".into(), (3, 0))],
        cell: None, // auto-detects CHIP_TOP over UNUSED_LOGO!
        split_islands: true,
        pin_layer: Some((1, 20)),
        simplify_tol: 0.0,
        ..Default::default()
    };
    let res = import_gds_bytes(&w.buf, &opts).unwrap();
    assert_eq!(res.cell, "CHIP_TOP");

    // Layer "sig" has 1 large transformed L_PAD (600 um^2, named "Q1_pad" via TEXT pin!)
    // + 6 equal-area AREF dots (24 um^2 each, ordered by ascending centroid x then y!).
    let sig_names = &res.layers.iter().find(|(k, _)| k == "sig").unwrap().1;
    assert_eq!(sig_names.len(), 7);
    assert_eq!(sig_names[0], "Q1_pad");
    assert_eq!(sig_names[1..], ["sig_1", "sig_2", "sig_3", "sig_4", "sig_5", "sig_6"]);
    assert_eq!(res.pins, vec![("Q1_pad".to_string(), "Q1_pad".to_string())]);

    let q1 = &res.sheets[0];
    assert!((q1.area - 600e-12).abs() < 1e-22, "q1.area = {:e}", q1.area);
    assert!((q1.bbox[0] - 100e-6).abs() < 1e-15);
    assert!((q1.bbox[1] - 50e-6).abs() < 1e-15);
    assert!((q1.bbox[2] - 120e-6).abs() < 1e-15);
    assert!((q1.bbox[3] - 90e-6).abs() < 1e-15);

    // Check deterministic centroid ordering of the 6 equal-area AREF dots:
    // Row ir=1 is at x = -120 um (bbox [-124..-120], cx = -122 um), ic = 0, 1, 2 => cy = 3, 23, 43 um
    // Row ir=0 is at x = -100 um (bbox [-104..-100], cx = -102 um), ic = 0, 1, 2 => cy = 3, 23, 43 um
    let expected_centroids =
        [[-122e-6, 3e-6], [-122e-6, 23e-6], [-122e-6, 43e-6], [-102e-6, 3e-6], [-102e-6, 23e-6], [-102e-6, 43e-6]];
    for (k, want_c) in expected_centroids.iter().enumerate() {
        let dot = &res.sheets[1 + k];
        assert!((dot.area - 24e-12).abs() < 1e-22);
        assert!((dot.centroid[0] - want_c[0]).abs() < 1e-15, "dot {k} cx = {:?}", dot.centroid);
        assert!((dot.centroid[1] - want_c[1]).abs() < 1e-15, "dot {k} cy = {:?}", dot.centroid);
    }

    // Check PATH expansion on "traces":
    // Two 500 um^2 traces (tie-broken by centroid: flush cx=25, cy=0 vs extend cx=20, cy=30 => extend has smaller cx=20!)
    // followed by the ~378.54 um^2 round-cap path!
    let tr_sheets: Vec<_> = res.sheets.iter().filter(|s| s.layer_key == "traces").collect();
    assert_eq!(tr_sheets.len(), 3);
    assert!((tr_sheets[0].area - 500e-12).abs() < 1e-22);
    for (got, want) in tr_sheets[0].bbox.iter().zip([-5e-6, 25e-6, 45e-6, 35e-6]) {
        assert!((got - want).abs() < 1e-18, "pathtype 2 bbox: got {got:e} vs {want:e}");
    }
    for (got, want) in tr_sheets[1].bbox.iter().zip([0.0, -5e-6, 50e-6, 5e-6]) {
        assert!((got - want).abs() < 1e-18, "pathtype 0 bbox: got {got:e} vs {want:e}");
    }
    let expected_round = 300e-12 + std::f64::consts::PI * 25e-12;
    // 12-segment semicircle polygon area is within ~1.2% of exact disk
    assert!((tr_sheets[2].area - expected_round).abs() < 0.02 * expected_round);
}

#[test]
fn pin_inside_hole_does_not_match_outer_ring_and_simplify_tol_pins_shared_contacts() {
    let mut w = GdsWriter::new();
    w.begin_lib("pin_hole_lib", 1.0e-3, 1.0e-9);
    w.begin_cell("MAIN");

    // Outer ground ring [-100..100] x [-100..100] um with hole [-50..50] x [-50..50] um
    w.boundary(
        1,
        0,
        &[
            [100_000, 100_000],
            [-100_000, 100_000],
            [-100_000, -50_000],
            [-50_000, -50_000],
            [-50_000, 50_000],
            [50_000, 50_000],
            [50_000, -50_000],
            [-50_000, -50_000],
            [-100_000, -50_000],
            [-100_000, -100_000],
            [100_000, -100_000],
        ],
    );

    // Inner pad [-20..20] x [-20..20] um with a tiny 20 nm jitter at (-10 um, 20.01 um) (unshared -> simplified)
    // and a 30 nm shared contact deflection at (0, 20.03 um) shared with JJ on layer (2, 0)!
    w.boundary(
        1,
        10,
        &[
            [-20_000, -20_000],
            [20_000, -20_000],
            [20_000, 20_000],
            [0, 20_030],       // shared with JJ -> must survive simplify_tol = 100 nm!
            [-10_000, 20_010], // unshared 10 nm jitter -> must be simplified away!
            [-20_000, 20_000],
        ],
    );
    w.boundary(2, 0, &[[-2_000, 20_000], [2_000, 20_000], [0, 20_030], [0, 30_000]]);

    // Pin at (0, 0) lies INSIDE the ground hole and INSIDE the inner pad:
    // MUST match ONLY the inner pad, NEVER the outer ground ring!
    w.text(1, 20, [0, 0], "island_center");
    // Pin at (-80 um, 0) lies on the outer ground metal:
    w.text(1, 20, [-80_000, 0], "ground_plane");

    w.end_cell();
    w.end_lib();

    let opts = GdsImportOptions {
        layers: vec![("gnd".into(), (1, 0)), ("sig".into(), (1, 10)), ("jj".into(), (2, 0))],
        cell: Some("MAIN".into()),
        split_islands: true,
        pin_layer: Some((1, 20)),
        simplify_tol: 100e-9,
        ..Default::default()
    };
    let res = import_gds_bytes(&w.buf, &opts).unwrap();
    assert_eq!(res.sheets[0].name, "ground_plane");
    assert_eq!(res.sheets[1].name, "island_center");
    assert_eq!(res.sheets[2].name, "jj");

    let pad_outer = &res.sheets[1].outers[0];
    assert!(
        pad_outer.iter().any(|p| p[0].abs() < 1e-15 && (p[1] - 20.03e-6).abs() < 1e-15),
        "shared contact vertex (0, 20.03 um) must be pinned across layers: {pad_outer:?}"
    );
    assert!(
        !pad_outer.iter().any(|p| (p[0] + 10e-6).abs() < 1e-15 && (p[1] - 20.01e-6).abs() < 1e-15),
        "unshared 10 nm jitter must be simplified away: {pad_outer:?}"
    );
}

#[test]
fn diagnostics_reject_malformed_or_cyclic_streams_and_split_islands_false_merges_cleanly() {
    assert!(parse_gds(&[]).is_err());
    assert!(parse_gds(&[0, 6, 0, 2, 2, 88]).is_err()); // missing UNITS & ENDLIB

    // Cyclic SREF hierarchy A -> B -> A must return a clean GeomError instead of overflowing stack
    let mut cyc = GdsWriter::new();
    cyc.begin_lib("cyclic", 1.0e-3, 1.0e-9);
    cyc.begin_cell("A");
    cyc.sref("B", [0, 0], false, None, None);
    cyc.end_cell();
    cyc.begin_cell("B");
    cyc.sref("A", [0, 0], false, None, None);
    cyc.end_cell();
    cyc.end_lib();

    let opts = GdsImportOptions {
        layers: vec![("m1".into(), (1, 0))],
        cell: Some("A".into()),
        split_islands: true,
        pin_layer: None,
        simplify_tol: 0.0,
        ..Default::default()
    };
    let err = import_gds_bytes(&cyc.buf, &opts).unwrap_err().0;
    assert!(err.contains("64 levels") || err.contains("cyclic"), "unexpected cyclic error: {err}");

    // split_islands = false merges multiple disconnected islands under one sheet name
    let mut w = GdsWriter::new();
    w.begin_lib("merge_lib", 1.0e-3, 1.0e-9);
    w.begin_cell("TOP");
    w.boundary(1, 10, &[[0, 0], [10_000, 0], [10_000, 10_000], [0, 10_000]]);
    w.boundary(1, 10, &[[30_000, 0], [50_000, 0], [50_000, 10_000], [30_000, 10_000]]);
    w.end_cell();
    w.end_lib();

    let merged = import_gds_bytes(
        &w.buf,
        &GdsImportOptions {
            layers: vec![("sig".into(), (1, 10))],
            cell: None,
            split_islands: false,
            pin_layer: None,
            simplify_tol: 0.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(merged.sheets.len(), 1);
    assert_eq!(merged.sheets[0].name, "sig");
    assert_eq!(merged.sheets[0].outers.len(), 2);
    assert!((merged.sheets[0].area - 300e-12).abs() < 1e-22);
}

#[test]
fn aref_bomb_and_binary_dag_sref_bomb_are_rejected_without_oom() {
    // 1. Massive 30,000 x 30,000 AREF (900M instances) on requested layer must fail in O(1) before looping
    let mut w = GdsWriter::new();
    w.begin_lib("aref_bomb", 1.0e-3, 1.0e-9);
    w.begin_cell("DOT");
    w.boundary(1, 0, &[[0, 0], [100, 0], [100, 100], [0, 100]]);
    w.end_cell();
    w.begin_cell("TOP");
    w.aref("DOT", 30_000, 30_000, [[0, 0], [30_000_000, 0], [0, 30_000_000]], false, None, None);
    w.end_cell();
    w.end_lib();

    let opts = GdsImportOptions { layers: vec![("m1".into(), (1, 0))], ..Default::default() };
    let err = import_gds_bytes(&w.buf, &opts).unwrap_err().0;
    assert!(err.contains("max_instances"), "unexpected AREF bomb error: {err}");

    // 2. Binary DAG SREF bomb (depth = 20 < 64, 2^20 = 1,048,576 instances > 250,000)
    let mut dag = GdsWriter::new();
    dag.begin_lib("dag_bomb", 1.0e-3, 1.0e-9);
    dag.begin_cell("L0");
    dag.boundary(1, 0, &[[0, 0], [100, 0], [100, 100], [0, 100]]);
    dag.end_cell();
    for lvl in 1..=20 {
        dag.begin_cell(&format!("L{lvl}"));
        let prev = format!("L{}", lvl - 1);
        dag.sref(&prev, [0, 0], false, None, None);
        dag.sref(&prev, [200, 0], false, None, None);
        dag.end_cell();
    }
    dag.end_lib();

    let err_dag = import_gds_bytes(&dag.buf, &opts).unwrap_err().0;
    assert!(err_dag.contains("max_instances"), "unexpected binary DAG SREF error: {err_dag}");

    // 3. Explicit max_file_bytes, max_records, and max_polygons limits
    let small_file_lim = GdsImportLimits { max_file_bytes: 32, ..Default::default() };
    assert!(parse_gds_with_limits(&w.buf, &small_file_lim).unwrap_err().0.contains("max_file_bytes"));
    let small_rec_lim = GdsImportLimits { max_records: 5, ..Default::default() };
    assert!(parse_gds_with_limits(&w.buf, &small_rec_lim).unwrap_err().0.contains("max_records"));
    let small_poly_lim = GdsImportLimits { max_polygons: 3, ..Default::default() };
    let mut polys = GdsWriter::new();
    polys.begin_lib("many_polys", 1.0e-3, 1.0e-9);
    polys.begin_cell("TOP");
    for k in 0..5 {
        polys.boundary(1, 0, &[[k * 200, 0], [k * 200 + 100, 0], [k * 200 + 100, 100], [k * 200, 100]]);
    }
    polys.end_cell();
    polys.end_lib();
    assert!(import_gds_bytes_with_limits(&polys.buf, &opts, &small_poly_lim, None)
        .unwrap_err()
        .0
        .contains("max_polygons"));
}

#[test]
fn deep_sref_chain_above_max_depth_is_rejected() {
    let mut w = GdsWriter::new();
    w.begin_lib("deep_sref", 1.0e-3, 1.0e-9);
    w.begin_cell("LEAF");
    w.boundary(1, 0, &[[0, 0], [100, 0], [100, 100], [0, 100]]);
    w.end_cell();
    for lvl in 1..=66 {
        w.begin_cell(&format!("N{lvl}"));
        let child = if lvl == 1 { "LEAF".to_string() } else { format!("N{}", lvl - 1) };
        w.sref(&child, [0, 0], false, None, None);
        w.end_cell();
    }
    w.end_lib();

    let opts = GdsImportOptions { layers: vec![("m1".into(), (1, 0))], ..Default::default() };
    let err = import_gds_bytes(&w.buf, &opts).unwrap_err().0;
    assert!(err.contains("64 levels"), "unexpected deep SREF error: {err}");
}

#[test]
fn oversized_decompose_raw_polygon_is_rejected_before_quadratic_work() {
    let n = 17_000usize;
    let mut pts: Vec<[i32; 2]> = Vec::with_capacity(n + 1);
    for i in 0..n {
        let ang = std::f64::consts::TAU * (i as f64 / n as f64);
        let r = 50_000.0;
        pts.push([(r * ang.cos()).round() as i32, (r * ang.sin()).round() as i32]);
    }
    let mut w = GdsWriter::new();
    w.begin_lib("big_poly", 1.0e-3, 1.0e-9);
    w.begin_cell("TOP");
    w.boundary(1, 0, &pts);
    w.end_cell();
    w.end_lib();

    let opts = GdsImportOptions { layers: vec![("m1".into(), (1, 0))], ..Default::default() };
    let err = import_gds_bytes(&w.buf, &opts).unwrap_err().0;
    assert!(err.contains("max_decompose_vertices"), "unexpected decompose limit error: {err}");
}

#[test]
fn non_finite_or_non_positive_rec_units_and_xf_overflow_are_rejected() {
    // 1. Non-positive user_unit (0.0 and -1e-3) in binary REC_UNITS
    for bad_uu in [0.0, -1.0e-3] {
        let mut w = GdsWriter::new();
        w.begin_lib("bad_uu", bad_uu, 1.0e-9);
        w.begin_cell("TOP");
        w.boundary(1, 0, &[[0, 0], [100, 0], [100, 100], [0, 100]]);
        w.end_cell();
        w.end_lib();
        let err = parse_gds(&w.buf).unwrap_err().0;
        assert!(err.contains("user unit"), "expected user_unit error for {bad_uu}, got: {err}");
    }

    // 2. Non-finite user_unit / dbu_m via programmatic GdsLibrary
    let opts = GdsImportOptions { layers: vec![("m1".into(), (1, 0))], ..Default::default() };
    for bad_uu in [f64::NAN, f64::INFINITY, -0.0] {
        let lib = GdsLibrary {
            version: 600,
            name: "nan_lib".into(),
            user_unit: bad_uu,
            dbu_m: 1e-9,
            cells: vec![GdsCell { name: "TOP".into(), ..Default::default() }],
        };
        let err = import_gds_library(&lib, &opts).unwrap_err().0;
        assert!(err.contains("user_unit"), "expected library user_unit error for {bad_uu}, got: {err}");
    }

    // 3. Finite GdsXf overflow during composition (6 nested SREFs each with mag = 1e60 => 1e360 = inf)
    let mut mag_bomb = GdsWriter::new();
    mag_bomb.begin_lib("xf_overflow", 1.0e-3, 1.0e-9);
    mag_bomb.begin_cell("LEAF");
    mag_bomb.boundary(1, 0, &[[0, 0], [100, 0], [100, 100], [0, 100]]);
    mag_bomb.end_cell();
    for lvl in 1..=6 {
        mag_bomb.begin_cell(&format!("S{lvl}"));
        let child = if lvl == 1 { "LEAF".to_string() } else { format!("S{}", lvl - 1) };
        mag_bomb.sref(&child, [0, 0], false, Some(1.0e60), None);
        mag_bomb.end_cell();
    }
    mag_bomb.end_lib();
    let err_xf = import_gds_bytes(&mag_bomb.buf, &opts).unwrap_err().0;
    assert!(err_xf.contains("non-finite"), "expected GdsXf non-finite error, got: {err_xf}");

    // Also direct GdsLibrary with NaN mag
    let nan_sref_lib = GdsLibrary {
        version: 600,
        name: "nan_xf".into(),
        user_unit: 1e-3,
        dbu_m: 1e-9,
        cells: vec![
            GdsCell {
                name: "TOP".into(),
                srefs: vec![GdsSref { sname: "C".into(), strans: 0, mag: f64::NAN, angle_deg: 0.0, xy: [0, 0] }],
                ..Default::default()
            },
            GdsCell { name: "C".into(), ..Default::default() },
        ],
    };
    let lim = GdsImportLimits::default();
    let err_nan_sref = import_gds_library_with_limits(&nan_sref_lib, &opts, &lim, None).unwrap_err().0;
    assert!(err_nan_sref.contains("magnification") || err_nan_sref.contains("non-finite"), "{err_nan_sref}");
}

#[test]
fn roi_and_layer_lazy_aref_clipping_skips_unwanted_million_hole_arrays() {
    // 1. 1,000 x 1,000 = 1,000,000 AREF holes on unrequested layer (99, 0) + 1,000 x 1,000 AREF pads on (1, 0)
    //    spanning [0..10 mm] x [0..10 mm] (pitch 10 um, 4x4 um pads).
    let mut w = GdsWriter::new();
    w.begin_lib("flux_trap_lib", 1.0e-3, 1.0e-9);
    w.begin_cell("TRAP");
    w.boundary(99, 0, &[[0, 0], [2_000, 0], [2_000, 2_000], [0, 2_000]]);
    w.end_cell();
    w.begin_cell("PAD");
    w.boundary(1, 0, &[[0, 0], [4_000, 0], [4_000, 4_000], [0, 4_000]]);
    w.end_cell();
    w.begin_cell("TOP");
    w.boundary(1, 0, &[[-20_000, -20_000], [-10_000, -20_000], [-10_000, -10_000], [-20_000, -10_000]]);
    w.aref("TRAP", 1_000, 1_000, [[0, 0], [10_000_000, 0], [0, 10_000_000]], false, None, None);
    w.end_cell();
    w.end_lib();

    // Importing layer (1, 0) skips the 1,000,000-element TRAP AREF on unrequested layer (99, 0) before allocation!
    let opts = GdsImportOptions { layers: vec![("m1".into(), (1, 0))], ..Default::default() };
    let res_layer = import_gds_bytes(&w.buf, &opts).unwrap();
    assert_eq!(res_layer.sheets.len(), 1);
    assert!((res_layer.sheets[0].area - 100e-12).abs() < 1e-22);

    // 2. Now put the 1,000 x 1,000 AREF on requested layer (1, 0), rotated 45 deg around origin (wait: 0 deg first, then 45 deg):
    let mut w_roi = GdsWriter::new();
    w_roi.begin_lib("roi_aref_lib", 1.0e-3, 1.0e-9);
    w_roi.begin_cell("PAD");
    w_roi.boundary(1, 0, &[[0, 0], [4_000, 0], [4_000, 4_000], [0, 4_000]]);
    w_roi.end_cell();
    w_roi.begin_cell("TOP");
    // 1000 x 1000 array, pitch 10 um -> covers [0..9994 um] x [0..9994 um]
    w_roi.aref("PAD", 1_000, 1_000, [[0, 0], [10_000_000, 0], [0, 10_000_000]], false, None, None);
    w_roi.end_cell();
    w_roi.end_lib();

    // Without ROI, 1,000,000 instances exceeds max_instances (250,000)
    assert!(import_gds_bytes(&w_roi.buf, &opts).unwrap_err().0.contains("max_instances"));

    // With ROI [99e-6, 99e-6, 115e-6, 115e-6]: hits ic in {10, 11} and ir in {10, 11} => 4 pads (each 16 um^2)!
    let res_roi = import_gds_bytes_with_limits(
        &w_roi.buf,
        &opts,
        &GdsImportLimits::default(),
        Some([99.0e-6, 99.0e-6, 115.0e-6, 115.0e-6]),
    )
    .unwrap();
    assert_eq!(res_roi.sheets.len(), 4, "ROI must select exactly the 2x2 overlapping AREF pads");
    for sh in &res_roi.sheets {
        assert!((sh.area - 16e-12).abs() < 1e-22);
        assert!(sh.bbox[0] >= 99.0e-6 && sh.bbox[2] <= 115.0e-6);
        assert!(sh.bbox[1] >= 99.0e-6 && sh.bbox[3] <= 115.0e-6);
    }

    // 3. Rotated 45-degree 1000 x 1000 AREF clipped to a tight ROI near the first element (0, 0):
    let mut w_rot = GdsWriter::new();
    w_rot.begin_lib("roi_rot_aref", 1.0e-3, 1.0e-9);
    w_rot.begin_cell("PAD");
    w_rot.boundary(1, 0, &[[0, 0], [2_000, 0], [2_000, 2_000], [0, 2_000]]);
    w_rot.end_cell();
    w_rot.begin_cell("TOP");
    // Stepped at 10 um along (1, 1)/sqrt(2) and (-1, 1)/sqrt(2) (~7071 nm per step)
    w_rot.aref("PAD", 1_000, 1_000, [[0, 0], [7_071_000, 7_071_000], [-7_071_000, 7_071_000]], false, None, Some(45.0));
    w_rot.end_cell();
    w_rot.end_lib();
    let res_rot = import_gds_bytes_with_limits(
        &w_rot.buf,
        &opts,
        &GdsImportLimits::default(),
        Some([-3.0e-6, -1.0e-6, 3.0e-6, 4.0e-6]),
    )
    .unwrap();
    assert_eq!(res_rot.sheets.len(), 1, "rotated 45-deg AREF ROI should keep only the origin pad");
    assert!((res_rot.sheets[0].area - 4e-12).abs() < 1e-22);
}
