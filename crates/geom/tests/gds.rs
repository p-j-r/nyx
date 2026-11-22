// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for the pure-Rust binary GDSII parser, IBM float codec, hierarchy
//! flattener (`SREF`/`AREF` with reflection, magnification, and rotation), `PATH` expansion,
//! keyhole & overlapping boolean polygon union, hull/hole orientation, island ordering, and
//! `TEXT` pin matching.

use geom::gds::{decode_ibm_f64, encode_ibm_f64, import_gds_bytes, parse_gds, GdsImportOptions, STRANS_REFLECT_X};
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
        let mut p = Vec::with_capacity(pts.len() * 8);
        for &[x, y] in pts {
            p.extend_from_slice(&x.to_be_bytes());
            p.extend_from_slice(&y.to_be_bytes());
        }
        self.rec(0x10, 0x03, &p);
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
