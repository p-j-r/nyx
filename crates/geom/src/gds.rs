// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-Rust binary GDSII stream parser, hierarchy flattener, `PATH` expander, keyhole/boolean
//! polygon merger (`arrange` + `poly_batch`), island splitter, and `TEXT` pin matcher.

use crate::arrange::arrange;
use crate::complex::Complex;
use crate::error::{GeomError, GeomResult};
use crate::math::*;
use crate::model::{check_name, Model};
use crate::poly_batch::{clean_ring, shoelace_area, simplify_rings_with_pins};
use crate::prim::polygon_area;
use crate::scene::{build_scene, Cache};
use std::collections::{BTreeMap, HashMap, HashSet};

// GDSII record type codes.
const REC_HEADER: u8 = 0x00;
const REC_BGNLIB: u8 = 0x01;
const REC_LIBNAME: u8 = 0x02;
const REC_UNITS: u8 = 0x03;
const REC_ENDLIB: u8 = 0x04;
const REC_BGNSTR: u8 = 0x05;
const REC_STRNAME: u8 = 0x06;
const REC_ENDSTR: u8 = 0x07;
const REC_BOUNDARY: u8 = 0x08;
const REC_PATH: u8 = 0x09;
const REC_SREF: u8 = 0x0A;
const REC_AREF: u8 = 0x0B;
const REC_TEXT: u8 = 0x0C;
const REC_LAYER: u8 = 0x0D;
const REC_DATATYPE: u8 = 0x0E;
const REC_WIDTH: u8 = 0x0F;
const REC_XY: u8 = 0x10;
const REC_ENDEL: u8 = 0x11;
const REC_SNAME: u8 = 0x12;
const REC_COLROW: u8 = 0x13;
const REC_NODE: u8 = 0x15;
const REC_TEXTTYPE: u8 = 0x16;
const REC_STRING: u8 = 0x19;
const REC_STRANS: u8 = 0x1A;
const REC_MAG: u8 = 0x1B;
const REC_ANGLE: u8 = 0x1C;
const REC_PATHTYPE: u8 = 0x21;
const REC_BOX: u8 = 0x2D;
const REC_BOXTYPE: u8 = 0x2E;
const REC_BGNEXTN: u8 = 0x30;
const REC_ENDEXTN: u8 = 0x31;

// GDSII data type codes.
const DT_BIT_ARRAY: u8 = 0x01;
const DT_I16: u8 = 0x02;
const DT_I32: u8 = 0x03;
const DT_REAL4: u8 = 0x04;
const DT_REAL8: u8 = 0x05;
const DT_ASCII: u8 = 0x06;

/// STRANS bit 15: reflection across the local X-axis (`y -> -y`) before scaling and rotation.
pub const STRANS_REFLECT_X: u16 = 0x8000;

/// Decodes an 8-byte big-endian IBM System/360 hexadecimal floating-point value to IEEE-754 `f64`.
pub fn decode_ibm_f64(bytes: [u8; 8]) -> f64 {
    let raw = u64::from_be_bytes(bytes);
    let mantissa = raw & 0x00ff_ffff_ffff_ffff;
    if mantissa == 0 {
        return 0.0;
    }
    let sign = if (bytes[0] & 0x80) != 0 { -1.0 } else { 1.0 };
    let exp16 = (bytes[0] & 0x7f) as i32 - 64;
    sign * (mantissa as f64) * 2.0_f64.powi(4 * exp16 - 56)
}

/// Encodes an IEEE-754 `f64` into an 8-byte big-endian IBM hexadecimal floating-point value.
pub fn encode_ibm_f64(val: f64) -> [u8; 8] {
    if val == 0.0 || !val.is_finite() {
        return [0u8; 8];
    }
    let sign = if val.is_sign_negative() { 0x80u8 } else { 0x00u8 };
    let mut a = val.abs();
    let mut exp_adjust = 0i32;
    if a < f64::MIN_POSITIVE {
        a *= 2.0_f64.powi(54);
        exp_adjust = -54;
    }
    let bits = a.to_bits();
    let raw_exp2 = ((bits >> 52) & 0x7ff) as i32;
    if raw_exp2 == 0 {
        return [0u8; 8];
    }
    let exp2 = raw_exp2 - 1023 + exp_adjust;
    let sig53 = (bits & 0x000f_ffff_ffff_ffff) | (1u64 << 52);
    let e2 = exp2 + 1;
    let r = (-e2).rem_euclid(4);
    let exp16 = (e2 + r) / 4;
    let biased = exp16 + 64;
    if biased < 0 {
        return [0u8; 8];
    }
    if biased > 127 {
        return [sign | 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff];
    }
    let mantissa = sig53 << (3 - r as u32);
    let mut out = mantissa.to_be_bytes();
    out[0] = sign | (biased as u8);
    out
}

fn decode_ibm_f32(bytes: [u8; 4]) -> f64 {
    let raw = u32::from_be_bytes(bytes);
    let mantissa = raw & 0x00ff_ffff;
    if mantissa == 0 {
        return 0.0;
    }
    let sign = if (bytes[0] & 0x80) != 0 { -1.0 } else { 1.0 };
    let exp16 = (bytes[0] & 0x7f) as i32 - 64;
    sign * (mantissa as f64) * 2.0_f64.powi(4 * exp16 - 24)
}

/// Parsed GDSII `BOUNDARY` (or `BOX`) polygon in integer database units.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GdsBoundary {
    pub layer: i16,
    pub datatype: i16,
    pub xy: Vec<[i32; 2]>,
}

/// Parsed GDSII `PATH` element in integer database units.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GdsPath {
    pub layer: i16,
    pub datatype: i16,
    pub pathtype: i16,
    pub width: i32,
    pub bgnextn: i32,
    pub endextn: i32,
    pub xy: Vec<[i32; 2]>,
}

/// Parsed GDSII `SREF` structure instance.
#[derive(Clone, Debug, PartialEq)]
pub struct GdsSref {
    pub sname: String,
    pub strans: u16,
    pub mag: f64,
    pub angle_deg: f64,
    pub xy: [i32; 2],
}

/// Parsed GDSII `AREF` array structure instance.
#[derive(Clone, Debug, PartialEq)]
pub struct GdsAref {
    pub sname: String,
    pub strans: u16,
    pub mag: f64,
    pub angle_deg: f64,
    pub cols: i16,
    pub rows: i16,
    pub xy: [[i32; 2]; 3],
}

/// Parsed GDSII `TEXT` label element.
#[derive(Clone, Debug, PartialEq)]
pub struct GdsText {
    pub layer: i16,
    pub texttype: i16,
    pub string: String,
    pub xy: [i32; 2],
    pub strans: u16,
    pub mag: f64,
    pub angle_deg: f64,
}

/// Parsed GDSII structure (`BGNSTR` .. `ENDSTR`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GdsCell {
    pub name: String,
    pub boundaries: Vec<GdsBoundary>,
    pub paths: Vec<GdsPath>,
    pub srefs: Vec<GdsSref>,
    pub arefs: Vec<GdsAref>,
    pub texts: Vec<GdsText>,
}

/// Parsed GDSII binary library (`HEADER` .. `ENDLIB`).
#[derive(Clone, Debug, PartialEq)]
pub struct GdsLibrary {
    pub version: i16,
    pub name: String,
    /// Database unit in user units (e.g. `1e-3` for `1 nm` in `um`).
    pub user_unit: f64,
    /// Database unit in SI metres (e.g. `1e-9` for `1 nm`).
    pub dbu_m: f64,
    pub cells: Vec<GdsCell>,
}

/// Options controlling GDSII hierarchy flattening, boolean union, island splitting, and pin matching.
#[derive(Clone, Debug)]
pub struct GdsImportOptions {
    /// Ordered `(layer_key, (layer, datatype))` pairs to extract.
    pub layers: Vec<(String, (i16, i16))>,
    /// Explicit top cell name; `None` auto-detects the root cell not referenced by any `SREF`/`AREF`.
    pub cell: Option<String>,
    /// Split disconnected islands on each layer into separate sheets (`true` by default).
    pub split_islands: bool,
    /// Optional `(layer, texttype)` carrying `TEXT` pin labels used to name containing islands.
    pub pin_layer: Option<(i16, i16)>,
    /// Closed-ring Ramer-Douglas-Peucker simplification tolerance in metres (`0.0` disables).
    pub simplify_tol: f64,
    /// Circular-arc fitting tolerance in metres (`0.0` disables, keeping polygons).
    pub arc_tol: f64,
}

impl Default for GdsImportOptions {
    fn default() -> Self {
        Self { layers: Vec::new(), cell: None, split_islands: true, pin_layer: None, simplify_tol: 0.0, arc_tol: 0.0 }
    }
}

/// One connected planar island: `(outer_ccw, holes_cw)`.
pub type IslandRings = (Vec<V2>, Vec<Vec<V2>>);

/// One extracted conductive sheet (a single connected island when `split_islands` is `true`, or
/// the merged multi-island layer when `split_islands` is `false`).
#[derive(Clone, Debug, PartialEq)]
pub struct GdsSheet {
    /// Resolved sheet name (matched pin label, `<layer_key>`, or `<layer_key>_<idx>`).
    pub name: String,
    /// Requested layer key from `GdsImportOptions::layers`.
    pub layer_key: String,
    /// Exterior polygon rings oriented counter-clockwise (`polygon_area > 0`).
    pub outers: Vec<Vec<V2>>,
    /// Interior hole rings oriented clockwise (`polygon_area < 0`).
    pub holes: Vec<Vec<V2>>,
    /// Grouped `(outer_ccw, holes_cw)` per connected island belonging to this sheet.
    pub islands: Vec<IslandRings>,
    /// Net conductive area in `m^2` (`sum(outer areas) - sum(hole areas)`).
    pub area: f64,
    /// Area-weighted 2D centroid `[cx, cy]` in metres.
    pub centroid: V2,
    /// 2D bounding box `[xmin, ymin, xmax, ymax]` in metres.
    pub bbox: [f64; 4],
}

/// Flattened, boolean-unioned GDSII import result in SI metres.
#[derive(Clone, Debug, PartialEq)]
pub struct GdsImportResult {
    /// Top cell name resolved and flattened.
    pub cell: String,
    /// Database unit in SI metres.
    pub dbu_m: f64,
    /// Extracted sheets in deterministic order.
    pub sheets: Vec<GdsSheet>,
    /// Sheet names grouped per requested layer key `(layer_key, sheet_names)`.
    pub layers: Vec<(String, Vec<String>)>,
    /// 2D bounding box `[xmin, ymin, xmax, ymax]` across all imported `layers` in metres.
    pub layer_bbox: [f64; 4],
    /// 2D bounding box `[xmin, ymin, xmax, ymax]` across all geometry in the flattened top cell in metres.
    pub cell_bbox: [f64; 4],
    /// Matched `TEXT` pin labels mapped to their resolved sheet name `(pin_label, sheet_name)`.
    pub pins: Vec<(String, String)>,
}

#[derive(Clone, Debug)]
enum ElemKind {
    Boundary,
    Path,
    Sref,
    Aref,
    Text,
    BoxElem,
    Other,
}

#[derive(Clone, Debug)]
struct ElemBuilder {
    kind: ElemKind,
    layer: i16,
    datatype: i16,
    texttype: Option<i16>,
    pathtype: i16,
    width: i32,
    bgnextn: i32,
    endextn: i32,
    sname: Option<String>,
    string: Option<String>,
    strans: u16,
    mag: f64,
    angle_deg: f64,
    colrow: Option<(i16, i16)>,
    xy: Vec<[i32; 2]>,
}

impl ElemBuilder {
    fn new(kind: ElemKind) -> Self {
        Self {
            kind,
            layer: 0,
            datatype: 0,
            texttype: None,
            pathtype: 0,
            width: 0,
            bgnextn: 0,
            endextn: 0,
            sname: None,
            string: None,
            strans: 0,
            mag: 1.0,
            angle_deg: 0.0,
            colrow: None,
            xy: Vec::new(),
        }
    }
}

fn read_i16(payload: &[u8], rec_name: &str) -> GeomResult<i16> {
    if payload.len() < 2 {
        return Err(GeomError(format!("GDSII {rec_name} record truncated ({} bytes)", payload.len())));
    }
    Ok(i16::from_be_bytes([payload[0], payload[1]]))
}

fn read_u16(payload: &[u8], rec_name: &str) -> GeomResult<u16> {
    if payload.len() < 2 {
        return Err(GeomError(format!("GDSII {rec_name} record truncated ({} bytes)", payload.len())));
    }
    Ok(u16::from_be_bytes([payload[0], payload[1]]))
}

fn read_i32(payload: &[u8], rec_name: &str) -> GeomResult<i32> {
    if payload.len() < 4 {
        return Err(GeomError(format!("GDSII {rec_name} record truncated ({} bytes)", payload.len())));
    }
    Ok(i32::from_be_bytes([payload[0], payload[1], payload[2], payload[3]]))
}

fn read_real(payload: &[u8], data_type: u8, rec_name: &str) -> GeomResult<f64> {
    match data_type {
        DT_REAL8 => {
            if payload.len() < 8 {
                return Err(GeomError(format!("GDSII {rec_name} 8-byte real truncated ({} bytes)", payload.len())));
            }
            let mut b = [0u8; 8];
            b.copy_from_slice(&payload[..8]);
            Ok(decode_ibm_f64(b))
        }
        DT_REAL4 => {
            if payload.len() < 4 {
                return Err(GeomError(format!("GDSII {rec_name} 4-byte real truncated ({} bytes)", payload.len())));
            }
            let mut b = [0u8; 4];
            b.copy_from_slice(&payload[..4]);
            Ok(decode_ibm_f32(b))
        }
        other => Err(GeomError(format!("GDSII {rec_name} expected floating-point data type, got 0x{other:02x}"))),
    }
}

fn read_ascii(payload: &[u8]) -> GeomResult<String> {
    let end = payload.iter().rposition(|&b| b != 0).map_or(0, |i| i + 1);
    std::str::from_utf8(&payload[..end])
        .map(|s| s.to_string())
        .map_err(|e| GeomError(format!("GDSII string is not valid ASCII/UTF-8: {e}")))
}

/// Parses a binary GDSII stream into a structured `GdsLibrary`.
pub fn parse_gds(bytes: &[u8]) -> GeomResult<GdsLibrary> {
    if bytes.len() < 4 {
        return Err(GeomError("GDSII stream is empty or shorter than a 4-byte record header".into()));
    }
    let mut pos = 0usize;
    let mut version: Option<i16> = None;
    let mut lib_name = String::new();
    let mut units: Option<(f64, f64)> = None;
    let mut cells: Vec<GdsCell> = Vec::new();
    let mut cur_cell: Option<GdsCell> = None;
    let mut cur_elem: Option<ElemBuilder> = None;
    let mut saw_endlib = false;

    while pos < bytes.len() {
        if pos + 4 > bytes.len() {
            return Err(GeomError(format!("GDSII stream has trailing truncated header at byte offset {pos}")));
        }
        let rec_len = u16::from_be_bytes([bytes[pos], bytes[pos + 1]]) as usize;
        if rec_len == 0 {
            return Err(GeomError(format!("GDSII stream has invalid zero-length record at byte offset {pos}")));
        }
        if rec_len < 4 || pos + rec_len > bytes.len() {
            return Err(GeomError(format!(
                "GDSII record at byte offset {pos} declares length {rec_len}, exceeding stream size {}",
                bytes.len()
            )));
        }
        let rec_type = bytes[pos + 2];
        let data_type = bytes[pos + 3];
        let payload = &bytes[pos + 4..pos + rec_len];
        pos += rec_len;

        match rec_type {
            REC_HEADER => {
                if data_type != DT_I16 {
                    return Err(GeomError(format!("GDSII HEADER expected i16 data type, got 0x{data_type:02x}")));
                }
                version = Some(read_i16(payload, "HEADER")?);
            }
            REC_BGNLIB => {}
            REC_LIBNAME => {
                if data_type != DT_ASCII && !payload.is_empty() {
                    return Err(GeomError(format!("GDSII LIBNAME expected ASCII data type, got 0x{data_type:02x}")));
                }
                lib_name = read_ascii(payload)?;
            }
            REC_UNITS => {
                let (uu, dbu) = match data_type {
                    DT_REAL8 if payload.len() >= 16 => {
                        let mut b0 = [0u8; 8];
                        let mut b1 = [0u8; 8];
                        b0.copy_from_slice(&payload[..8]);
                        b1.copy_from_slice(&payload[8..16]);
                        (decode_ibm_f64(b0), decode_ibm_f64(b1))
                    }
                    DT_REAL4 if payload.len() >= 8 => {
                        let mut b0 = [0u8; 4];
                        let mut b1 = [0u8; 4];
                        b0.copy_from_slice(&payload[..4]);
                        b1.copy_from_slice(&payload[4..8]);
                        (decode_ibm_f32(b0), decode_ibm_f32(b1))
                    }
                    _ => {
                        return Err(GeomError(format!(
                            "GDSII UNITS record malformed (data_type=0x{data_type:02x}, payload={} bytes)",
                            payload.len()
                        )));
                    }
                };
                if !dbu.is_finite() || dbu <= 0.0 {
                    return Err(GeomError(format!("GDSII UNITS database unit in metres ({dbu}) must be positive")));
                }
                units = Some((uu, dbu));
            }
            REC_ENDLIB => {
                if cur_elem.is_some() || cur_cell.is_some() {
                    return Err(GeomError(
                        "GDSII ENDLIB encountered while a structure or element is still open".into(),
                    ));
                }
                saw_endlib = true;
                break;
            }
            REC_BGNSTR => {
                if cur_cell.is_some() {
                    return Err(GeomError("nested GDSII BGNSTR without matching ENDSTR".into()));
                }
                cur_cell = Some(GdsCell::default());
            }
            REC_STRNAME => {
                let c = cur_cell
                    .as_mut()
                    .ok_or_else(|| GeomError("GDSII STRNAME encountered outside BGNSTR..ENDSTR".into()))?;
                c.name = read_ascii(payload)?;
            }
            REC_ENDSTR => {
                if cur_elem.is_some() {
                    return Err(GeomError("GDSII ENDSTR encountered before ENDEL".into()));
                }
                let c = cur_cell.take().ok_or_else(|| GeomError("GDSII ENDSTR without matching BGNSTR".into()))?;
                if c.name.is_empty() {
                    return Err(GeomError("GDSII structure is missing STRNAME".into()));
                }
                cells.push(c);
            }
            REC_BOUNDARY => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII BOUNDARY in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::Boundary));
            }
            REC_PATH => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII PATH in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::Path));
            }
            REC_SREF => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII SREF in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::Sref));
            }
            REC_AREF => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII AREF in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::Aref));
            }
            REC_TEXT => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII TEXT in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::Text));
            }
            REC_BOX => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII BOX in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::BoxElem));
            }
            REC_NODE => {
                if cur_cell.is_none() || cur_elem.is_some() {
                    return Err(GeomError("GDSII NODE in invalid context".into()));
                }
                cur_elem = Some(ElemBuilder::new(ElemKind::Other));
            }
            REC_LAYER => {
                if let Some(el) = cur_elem.as_mut() {
                    el.layer = read_i16(payload, "LAYER")?;
                }
            }
            REC_DATATYPE | REC_BOXTYPE => {
                if let Some(el) = cur_elem.as_mut() {
                    el.datatype = read_i16(payload, "DATATYPE")?;
                }
            }
            REC_TEXTTYPE => {
                if let Some(el) = cur_elem.as_mut() {
                    el.texttype = Some(read_i16(payload, "TEXTTYPE")?);
                }
            }
            REC_PATHTYPE => {
                if let Some(el) = cur_elem.as_mut() {
                    el.pathtype = read_i16(payload, "PATHTYPE")?;
                }
            }
            REC_WIDTH => {
                if let Some(el) = cur_elem.as_mut() {
                    el.width = read_i32(payload, "WIDTH")?;
                }
            }
            REC_BGNEXTN => {
                if let Some(el) = cur_elem.as_mut() {
                    el.bgnextn = read_i32(payload, "BGNEXTN")?;
                }
            }
            REC_ENDEXTN => {
                if let Some(el) = cur_elem.as_mut() {
                    el.endextn = read_i32(payload, "ENDEXTN")?;
                }
            }
            REC_SNAME => {
                if let Some(el) = cur_elem.as_mut() {
                    el.sname = Some(read_ascii(payload)?);
                }
            }
            REC_STRING => {
                if let Some(el) = cur_elem.as_mut() {
                    el.string = Some(read_ascii(payload)?);
                }
            }
            REC_STRANS => {
                if data_type != DT_BIT_ARRAY && data_type != DT_I16 {
                    return Err(GeomError(format!("GDSII STRANS expected 2-byte bit array, got 0x{data_type:02x}")));
                }
                if let Some(el) = cur_elem.as_mut() {
                    el.strans = read_u16(payload, "STRANS")?;
                }
            }
            REC_MAG => {
                let m = read_real(payload, data_type, "MAG")?;
                if !m.is_finite() || m <= 0.0 {
                    return Err(GeomError(format!("GDSII MAG ({m}) must be finite and positive")));
                }
                if let Some(el) = cur_elem.as_mut() {
                    el.mag = m;
                }
            }
            REC_ANGLE => {
                let a = read_real(payload, data_type, "ANGLE")?;
                if !a.is_finite() {
                    return Err(GeomError(format!("GDSII ANGLE ({a}) must be finite")));
                }
                if let Some(el) = cur_elem.as_mut() {
                    el.angle_deg = a;
                }
            }
            REC_COLROW => {
                if payload.len() < 4 {
                    return Err(GeomError(format!("GDSII COLROW truncated ({} bytes)", payload.len())));
                }
                let cols = i16::from_be_bytes([payload[0], payload[1]]);
                let rows = i16::from_be_bytes([payload[2], payload[3]]);
                if cols <= 0 || rows <= 0 {
                    return Err(GeomError(format!("GDSII COLROW ({cols}, {rows}) must have positive counts")));
                }
                if let Some(el) = cur_elem.as_mut() {
                    el.colrow = Some((cols, rows));
                }
            }
            REC_XY => {
                if data_type != DT_I32 || !payload.len().is_multiple_of(8) {
                    return Err(GeomError(format!(
                        "GDSII XY expected 8-byte i32 pairs, got data_type=0x{data_type:02x}, len={}",
                        payload.len()
                    )));
                }
                if let Some(el) = cur_elem.as_mut() {
                    el.xy.reserve(payload.len() / 8);
                    for chunk in payload.as_chunks::<8>().0 {
                        let x = i32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                        let y = i32::from_be_bytes([chunk[4], chunk[5], chunk[6], chunk[7]]);
                        el.xy.push([x, y]);
                    }
                }
            }
            REC_ENDEL => {
                let el = cur_elem.take().ok_or_else(|| GeomError("GDSII ENDEL without open element".into()))?;
                let cell = cur_cell.as_mut().ok_or_else(|| GeomError("GDSII ENDEL outside BGNSTR".into()))?;
                match el.kind {
                    ElemKind::Boundary | ElemKind::BoxElem => {
                        if el.xy.len() < 3 {
                            return Err(GeomError(format!(
                                "GDSII BOUNDARY in cell '{}' has fewer than 3 XY points",
                                cell.name
                            )));
                        }
                        cell.boundaries.push(GdsBoundary { layer: el.layer, datatype: el.datatype, xy: el.xy });
                    }
                    ElemKind::Path => {
                        if el.xy.len() < 2 {
                            return Err(GeomError(format!(
                                "GDSII PATH in cell '{}' has fewer than 2 XY points",
                                cell.name
                            )));
                        }
                        cell.paths.push(GdsPath {
                            layer: el.layer,
                            datatype: el.datatype,
                            pathtype: el.pathtype,
                            width: el.width,
                            bgnextn: el.bgnextn,
                            endextn: el.endextn,
                            xy: el.xy,
                        });
                    }
                    ElemKind::Sref => {
                        let sname = el
                            .sname
                            .ok_or_else(|| GeomError(format!("GDSII SREF in cell '{}' lacks SNAME", cell.name)))?;
                        let &xy = el
                            .xy
                            .first()
                            .ok_or_else(|| GeomError(format!("GDSII SREF in cell '{}' lacks XY", cell.name)))?;
                        cell.srefs.push(GdsSref { sname, strans: el.strans, mag: el.mag, angle_deg: el.angle_deg, xy });
                    }
                    ElemKind::Aref => {
                        let sname = el
                            .sname
                            .ok_or_else(|| GeomError(format!("GDSII AREF in cell '{}' lacks SNAME", cell.name)))?;
                        let (cols, rows) = el
                            .colrow
                            .ok_or_else(|| GeomError(format!("GDSII AREF in cell '{}' lacks COLROW", cell.name)))?;
                        if el.xy.len() < 3 {
                            return Err(GeomError(format!(
                                "GDSII AREF in cell '{}' requires 3 XY points, got {}",
                                cell.name,
                                el.xy.len()
                            )));
                        }
                        cell.arefs.push(GdsAref {
                            sname,
                            strans: el.strans,
                            mag: el.mag,
                            angle_deg: el.angle_deg,
                            cols,
                            rows,
                            xy: [el.xy[0], el.xy[1], el.xy[2]],
                        });
                    }
                    ElemKind::Text => {
                        let &xy = el
                            .xy
                            .first()
                            .ok_or_else(|| GeomError(format!("GDSII TEXT in cell '{}' lacks XY", cell.name)))?;
                        let texttype = el.texttype.unwrap_or(el.datatype);
                        cell.texts.push(GdsText {
                            layer: el.layer,
                            texttype,
                            string: el.string.unwrap_or_default(),
                            xy,
                            strans: el.strans,
                            mag: el.mag,
                            angle_deg: el.angle_deg,
                        });
                    }
                    ElemKind::Other => {}
                }
            }
            _ => {}
        }
    }

    let version = version.ok_or_else(|| GeomError("GDSII stream is missing HEADER record".into()))?;
    let (user_unit, dbu_m) = units.ok_or_else(|| GeomError("GDSII stream is missing UNITS record".into()))?;
    if !saw_endlib {
        return Err(GeomError("GDSII stream ended before ENDLIB record".into()));
    }
    Ok(GdsLibrary { version, name: lib_name, user_unit, dbu_m, cells })
}

// ------------------------------------------------------------------------------------------------- hierarchy flattening

#[derive(Clone, Copy, Debug)]
struct GdsXf {
    m: [[f64; 2]; 2],
    t: V2,
}

fn sin_cos_deg(deg: f64) -> (f64, f64) {
    let d = deg.rem_euclid(360.0);
    if d < 1e-12 || (360.0 - d) < 1e-12 {
        (0.0, 1.0)
    } else if (d - 90.0).abs() < 1e-12 {
        (1.0, 0.0)
    } else if (d - 180.0).abs() < 1e-12 {
        (0.0, -1.0)
    } else if (d - 270.0).abs() < 1e-12 {
        (-1.0, 0.0)
    } else {
        deg.to_radians().sin_cos()
    }
}

impl GdsXf {
    const ID: Self = Self { m: [[1.0, 0.0], [0.0, 1.0]], t: [0.0, 0.0] };

    fn from_instance(strans: u16, mag: f64, angle_deg: f64, origin: V2) -> Self {
        let sy = if (strans & STRANS_REFLECT_X) != 0 { -1.0 } else { 1.0 };
        let (sin, cos) = sin_cos_deg(angle_deg);
        Self { m: [[mag * cos, -mag * sin * sy], [mag * sin, mag * cos * sy]], t: origin }
    }

    fn apply(&self, p: V2) -> V2 {
        [self.m[0][0] * p[0] + self.m[0][1] * p[1] + self.t[0], self.m[1][0] * p[0] + self.m[1][1] * p[1] + self.t[1]]
    }

    fn compose(&self, child: &Self) -> Self {
        let a = &self.m;
        let b = &child.m;
        Self {
            m: [
                [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
                [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
            ],
            t: self.apply(child.t),
        }
    }

    fn scale(&self) -> f64 {
        (self.m[0][0] * self.m[1][1] - self.m[0][1] * self.m[1][0]).abs().sqrt()
    }
}

#[derive(Clone, Debug)]
struct FlatPin {
    label: String,
    pos: V2,
}

struct FlatWorld {
    /// Raw polygons grouped by `(layer, datatype)` in SI metres.
    polys: HashMap<(i16, i16), Vec<Vec<V2>>>,
    /// Bounding box across all geometry on all layers in the flattened cell.
    all_bbox: Aabb,
    /// Flattened `TEXT` pins matching `pin_layer`.
    pins: Vec<FlatPin>,
}

impl Default for FlatWorld {
    fn default() -> Self {
        Self { polys: HashMap::new(), all_bbox: Aabb::EMPTY, pins: Vec::new() }
    }
}

/// Expands a GDSII `PATH` centerline into one or more closed counter-clockwise polygons in metres.
pub fn expand_path_polygons(centerline: &[V2], width: f64, pathtype: i16, bgnextn: f64, endextn: f64) -> Vec<Vec<V2>> {
    if width <= 0.0 || !width.is_finite() {
        return Vec::new();
    }
    let eps = 1e-13 * width.max(1e-12);
    let mut pts: Vec<V2> = Vec::with_capacity(centerline.len());
    for &p in centerline {
        if pts.last().is_none_or(|&q| dist2(p, q) > eps) {
            pts.push(p);
        }
    }
    if pts.len() < 2 {
        return Vec::new();
    }
    let hw = 0.5 * width;
    let m = pts.len() - 1;
    let mut tangents = Vec::with_capacity(m);
    let mut normals = Vec::with_capacity(m);
    let mut seg_lens = Vec::with_capacity(m);
    for i in 0..m {
        let d = sub2(pts[i + 1], pts[i]);
        let l = norm2(d);
        seg_lens.push(l);
        let t = [d[0] / l, d[1] / l];
        tangents.push(t);
        normals.push([-t[1], t[0]]);
    }

    // Build single mitered/capped contour first; fall back to segment-by-segment convex pieces if tight bends self-overlap.
    let mut right_rail: Vec<V2> = Vec::new();
    let mut left_rail: Vec<V2> = Vec::new();

    // Start point / extension
    let p_start = match pathtype {
        2 => sub2(pts[0], mul2(tangents[0], hw)),
        4 => sub2(pts[0], mul2(tangents[0], bgnextn)),
        _ => pts[0],
    };
    let p_end = match pathtype {
        2 => add2(pts[m], mul2(tangents[m - 1], hw)),
        4 => add2(pts[m], mul2(tangents[m - 1], endextn)),
        _ => pts[m],
    };

    right_rail.push(sub2(p_start, mul2(normals[0], hw)));
    left_rail.push(add2(p_start, mul2(normals[0], hw)));

    let mut clean_miter = true;
    for k in 1..m {
        let (t0, t1) = (tangents[k - 1], tangents[k]);
        let (n0, n1) = (normals[k - 1], normals[k]);
        let denom = 1.0 + dot2(t0, t1);
        let turn = cross2(t0, t1);
        if denom <= 0.2 {
            // Very sharp hairpin (> ~127 deg turn): bevel outer side, check inner projection
            clean_miter = false;
            break;
        }
        let miter = mul2(add2(n0, n1), hw / denom);
        let proj = hw * turn.abs() / denom;
        if proj > 0.95 * seg_lens[k - 1].min(seg_lens[k]) {
            clean_miter = false;
            break;
        }
        right_rail.push(sub2(pts[k], miter));
        left_rail.push(add2(pts[k], miter));
    }

    if clean_miter {
        let mut ring = Vec::new();
        ring.extend_from_slice(&right_rail);
        if pathtype == 1 {
            // Semicircular end cap around pts[m] from -n_{m-1} through +t_{m-1} to +n_{m-1}
            let n_steps = 12usize;
            let t_end = tangents[m - 1];
            let n_end = normals[m - 1];
            for s in 0..=n_steps {
                let phi = -std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (s as f64 / n_steps as f64);
                let (sn, cs) = phi.sin_cos();
                let dir = add2(mul2(t_end, cs), mul2(n_end, sn));
                let pt = add2(pts[m], mul2(dir, hw));
                if ring.last().is_none_or(|&q| dist2(pt, q) > eps) {
                    ring.push(pt);
                }
            }
        } else {
            ring.push(sub2(p_end, mul2(normals[m - 1], hw)));
            ring.push(add2(p_end, mul2(normals[m - 1], hw)));
        }
        for &lp in left_rail.iter().rev() {
            if ring.last().is_none_or(|&q| dist2(lp, q) > eps) {
                ring.push(lp);
            }
        }
        if pathtype == 1 {
            // Semicircular start cap around pts[0] from +n_0 through -t_0 to -n_0
            let n_steps = 12usize;
            let t0 = tangents[0];
            let n0 = normals[0];
            for s in 1..n_steps {
                let phi = std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (s as f64 / n_steps as f64);
                let (sn, cs) = phi.sin_cos();
                let dir = add2(mul2(t0, cs), mul2(n0, sn));
                let pt = add2(pts[0], mul2(dir, hw));
                if ring.last().is_none_or(|&q| dist2(pt, q) > eps) {
                    ring.push(pt);
                }
            }
        }
        let ring = clean_ring(&ring, eps);
        if ring.len() >= 3 && is_simple_ring(&ring, eps) && polygon_area(&ring) > 0.0 {
            return vec![ring];
        }
    }

    // Fallback for tight self-overlapping bends: emit one convex rectangle/capsule per segment plus round/square vertex joins
    // (boolean union in `arrange` merges them cleanly).
    let mut out = Vec::new();
    for i in 0..m {
        let a = if i == 0 { p_start } else { pts[i] };
        let b = if i + 1 == m { p_end } else { pts[i + 1] };
        if dist2(a, b) > eps {
            let n = normals[i];
            out.push(vec![sub2(a, mul2(n, hw)), sub2(b, mul2(n, hw)), add2(b, mul2(n, hw)), add2(a, mul2(n, hw))]);
        }
    }
    let disk = |c: V2| -> Vec<V2> {
        let n_steps = 16usize;
        (0..n_steps)
            .map(|s| {
                let a = TAU * (s as f64 / n_steps as f64);
                [c[0] + hw * a.cos(), c[1] + hw * a.sin()]
            })
            .collect()
    };
    for &p in &pts[1..m] {
        out.push(disk(p));
    }
    if pathtype == 1 {
        out.push(disk(pts[0]));
        out.push(disk(pts[m]));
    }
    out
}

fn seg_properly_crosses(a: V2, b: V2, c: V2, d: V2, tol: f64) -> bool {
    let ab = sub2(b, a);
    let cd = sub2(d, c);
    let d1 = cross2(ab, sub2(c, a));
    let d2 = cross2(ab, sub2(d, a));
    let d3 = cross2(cd, sub2(a, c));
    let d4 = cross2(cd, sub2(b, c));
    let scale = norm2(ab).max(norm2(cd)) * tol;
    if d1.abs() <= scale && point_on_open_seg(c, a, b, tol) {
        return true;
    }
    if d2.abs() <= scale && point_on_open_seg(d, a, b, tol) {
        return true;
    }
    if d3.abs() <= scale && point_on_open_seg(a, c, d, tol) {
        return true;
    }
    if d4.abs() <= scale && point_on_open_seg(b, c, d, tol) {
        return true;
    }
    (d1 > scale && d2 < -scale || d1 < -scale && d2 > scale) && (d3 > scale && d4 < -scale || d3 < -scale && d4 > scale)
}

fn point_seg_dist(p: V2, a: V2, b: V2) -> f64 {
    let ab = sub2(b, a);
    let l2 = dot2(ab, ab);
    let t = if l2 > 0.0 { (dot2(sub2(p, a), ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    dist2(p, add2(a, mul2(ab, t)))
}

fn point_on_open_seg(p: V2, a: V2, b: V2, tol: f64) -> bool {
    if dist2(p, a) <= tol || dist2(p, b) <= tol {
        return false;
    }
    point_seg_dist(p, a, b) <= tol
}

fn is_simple_ring(ring: &[V2], tol: f64) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        for j in i + 1..n {
            if dist2(ring[i], ring[j]) <= tol {
                return false;
            }
        }
    }
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        for j in i + 2..n {
            if i == 0 && j + 1 == n {
                continue;
            }
            let (c, d) = (ring[j], ring[(j + 1) % n]);
            if seg_properly_crosses(a, b, c, d, tol) {
                return false;
            }
        }
    }
    true
}

fn count_transitive(
    cells: &[GdsCell],
    cell_map: &HashMap<&str, usize>,
    idx: usize,
    req_layers: &HashSet<(i16, i16)>,
    memo: &mut HashMap<usize, (bool, usize)>,
    visiting: &mut HashSet<usize>,
) -> (bool, usize) {
    if let Some(&res) = memo.get(&idx) {
        return res;
    }
    if !visiting.insert(idx) {
        return (false, 0);
    }
    let c = &cells[idx];
    let mut has_req = c.boundaries.iter().any(|b| req_layers.contains(&(b.layer, b.datatype)))
        || c.paths.iter().any(|p| req_layers.contains(&(p.layer, p.datatype)));
    let mut total = c.boundaries.len() + c.paths.len() + c.texts.len();
    for s in &c.srefs {
        if let Some(&ci) = cell_map.get(s.sname.as_str()) {
            let (hr, t) = count_transitive(cells, cell_map, ci, req_layers, memo, visiting);
            has_req |= hr;
            total = total.saturating_add(t);
        }
    }
    for a in &c.arefs {
        if let Some(&ci) = cell_map.get(a.sname.as_str()) {
            let (hr, t) = count_transitive(cells, cell_map, ci, req_layers, memo, visiting);
            has_req |= hr;
            let mult = (a.cols as usize).saturating_mul(a.rows as usize);
            total = total.saturating_add(t.saturating_mul(mult));
        }
    }
    visiting.remove(&idx);
    memo.insert(idx, (has_req, total));
    (has_req, total)
}

fn resolve_top_cell(lib: &GdsLibrary, explicit: Option<&str>, req_layers: &HashSet<(i16, i16)>) -> GeomResult<usize> {
    if lib.cells.is_empty() {
        return Err(GeomError("GDSII library contains no structures (BGNSTR)".into()));
    }
    let cell_map: HashMap<&str, usize> = lib.cells.iter().enumerate().map(|(i, c)| (c.name.as_str(), i)).collect();
    if let Some(name) = explicit {
        return cell_map.get(name).copied().ok_or_else(|| {
            let names: Vec<&str> = lib.cells.iter().map(|c| c.name.as_str()).collect();
            GeomError(format!("GDSII library has no cell '{name}' (available cells: {})", names.join(", ")))
        });
    }
    let mut referenced: HashSet<&str> = HashSet::new();
    for c in &lib.cells {
        for s in &c.srefs {
            referenced.insert(s.sname.as_str());
        }
        for a in &c.arefs {
            referenced.insert(a.sname.as_str());
        }
    }
    let roots: Vec<usize> =
        (0..lib.cells.len()).filter(|&i| !referenced.contains(lib.cells[i].name.as_str())).collect();
    if roots.is_empty() {
        return Err(GeomError(
            "GDSII library has no unreferenced root cell (cyclic hierarchy); specify `cell` explicitly".into(),
        ));
    }
    if roots.len() == 1 {
        return Ok(roots[0]);
    }
    let mut memo = HashMap::new();
    let mut visiting = HashSet::new();
    Ok(*roots
        .iter()
        .max_by_key(|&&i| {
            let (has_req, total) = count_transitive(&lib.cells, &cell_map, i, req_layers, &mut memo, &mut visiting);
            (has_req, total, i)
        })
        .unwrap())
}

fn flatten_rec(
    lib: &GdsLibrary,
    cell_map: &HashMap<&str, usize>,
    cell_idx: usize,
    xf: GdsXf,
    depth: usize,
    pin_layer: Option<(i16, i16)>,
    out: &mut FlatWorld,
) -> GeomResult<()> {
    if depth > 64 {
        return Err(GeomError(format!(
            "GDSII hierarchy exceeds 64 levels at cell '{}' (cyclic SREF/AREF?)",
            lib.cells[cell_idx].name
        )));
    }
    let dbu_m = lib.dbu_m;
    let to_m = |p: [i32; 2]| -> V2 { [p[0] as f64 * dbu_m, p[1] as f64 * dbu_m] };
    let cell = &lib.cells[cell_idx];

    for b in &cell.boundaries {
        let mut pts: Vec<V2> = b.xy.iter().map(|&p| xf.apply(to_m(p))).collect();
        for &q in &pts {
            out.all_bbox.grow([q[0], q[1], 0.0]);
        }
        if pts.len() >= 3 {
            if polygon_area(&pts) < 0.0 {
                pts.reverse();
            }
            out.polys.entry((b.layer, b.datatype)).or_default().push(pts);
        }
    }

    for p in &cell.paths {
        if p.width == 0 {
            continue;
        }
        let scale = if p.width < 0 { 1.0 } else { xf.scale() };
        let w_m = (p.width.unsigned_abs() as f64) * dbu_m * scale;
        let bgn_m = (p.bgnextn as f64) * dbu_m * scale;
        let end_m = (p.endextn as f64) * dbu_m * scale;
        let center: Vec<V2> = p.xy.iter().map(|&q| xf.apply(to_m(q))).collect();
        for poly in expand_path_polygons(&center, w_m, p.pathtype, bgn_m, end_m) {
            for &q in &poly {
                out.all_bbox.grow([q[0], q[1], 0.0]);
            }
            out.polys.entry((p.layer, p.datatype)).or_default().push(poly);
        }
    }

    if let Some((pl, ptt)) = pin_layer {
        for t in &cell.texts {
            if t.layer == pl && t.texttype == ptt {
                let label = t.string.trim().to_string();
                if !label.is_empty() {
                    out.pins.push(FlatPin { label, pos: xf.apply(to_m(t.xy)) });
                }
            }
        }
    }

    for s in &cell.srefs {
        let &target = cell_map
            .get(s.sname.as_str())
            .ok_or_else(|| GeomError(format!("GDSII SREF in '{}' references unknown cell '{}'", cell.name, s.sname)))?;
        let child_xf = GdsXf::from_instance(s.strans, s.mag, s.angle_deg, to_m(s.xy));
        flatten_rec(lib, cell_map, target, xf.compose(&child_xf), depth + 1, pin_layer, out)?;
    }

    for a in &cell.arefs {
        let &target = cell_map
            .get(a.sname.as_str())
            .ok_or_else(|| GeomError(format!("GDSII AREF in '{}' references unknown cell '{}'", cell.name, a.sname)))?;
        let p0 = to_m(a.xy[0]);
        let p1 = to_m(a.xy[1]);
        let p2 = to_m(a.xy[2]);
        let col_step = mul2(sub2(p1, p0), 1.0 / (a.cols as f64));
        let row_step = mul2(sub2(p2, p0), 1.0 / (a.rows as f64));
        for ir in 0..a.rows {
            for ic in 0..a.cols {
                let origin = add2(p0, add2(mul2(col_step, ic as f64), mul2(row_step, ir as f64)));
                let child_xf = GdsXf::from_instance(a.strans, a.mag, a.angle_deg, origin);
                flatten_rec(lib, cell_map, target, xf.compose(&child_xf), depth + 1, pin_layer, out)?;
            }
        }
    }

    Ok(())
}

// ------------------------------------------------------------------------------------------------- keyhole & boolean polygon processing

fn strip_collinear(ring: &[V2], tol: f64) -> Vec<V2> {
    let mut pts = clean_ring(ring, tol);
    loop {
        let n = pts.len();
        if n <= 3 {
            break;
        }
        let mut removed = false;
        let mut next_pts = Vec::with_capacity(n);
        for i in 0..n {
            let prev = pts[(i + n - 1) % n];
            let cur = pts[i];
            let nxt = pts[(i + 1) % n];
            if !removed && pts.len() > 3 && dist2(prev, nxt) > tol && point_seg_dist(cur, prev, nxt) <= tol {
                removed = true;
            } else {
                next_pts.push(cur);
            }
        }
        pts = clean_ring(&next_pts, tol);
        if !removed {
            break;
        }
    }
    pts
}

fn canonicalize_ring(mut ring: Vec<V2>, want_ccw: bool) -> Vec<V2> {
    if ring.len() < 3 {
        return ring;
    }
    let area = polygon_area(&ring);
    if (want_ccw && area < 0.0) || (!want_ccw && area > 0.0) {
        ring.reverse();
    }
    if let Some(i0) =
        (0..ring.len()).min_by(|&a, &b| ring[a][0].total_cmp(&ring[b][0]).then(ring[a][1].total_cmp(&ring[b][1])))
    {
        ring.rotate_left(i0);
    }
    ring
}

fn point_in_or_on_ring(ring: &[V2], p: V2, tol: f64) -> bool {
    let n = ring.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        if point_seg_dist(p, a, b) <= tol {
            return true;
        }
        if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0] {
            inside = !inside;
        }
    }
    inside
}

fn point_strictly_in_ring(ring: &[V2], p: V2, tol: f64) -> bool {
    let n = ring.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        if point_seg_dist(p, a, b) <= tol {
            return false;
        }
        if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0] {
            inside = !inside;
        }
    }
    inside
}

/// Decomposes one raw GDSII `BOUNDARY` polygon (which may contain keyhole cutlines `u -> v` and
/// `v -> u` connecting an exterior hull to one or more interior CW holes) into simple `(outer_ccw, holes_cw)` pairs.
fn decompose_raw_polygon(raw: &[V2]) -> Vec<IslandRings> {
    let mut bb = Aabb::EMPTY;
    for &p in raw {
        bb.grow([p[0], p[1], 0.0]);
    }
    let diag = bb.diagonal();
    if !diag.is_finite() || diag <= 1e-30 {
        return Vec::new();
    }
    let snap_tol = 1e-12 * diag;
    let mut pts = clean_ring(raw, snap_tol);
    if pts.len() < 3 || shoelace_area(&pts) <= 1e-24 * diag * diag {
        return Vec::new();
    }
    if polygon_area(&pts) < 0.0 {
        pts.reverse();
    }
    if is_simple_ring(&pts, snap_tol) {
        let c = strip_collinear(&pts, snap_tol);
        if c.len() >= 3 && shoelace_area(&c) > 1e-24 * diag * diag {
            return vec![(canonicalize_ring(c, true), Vec::new())];
        }
        return Vec::new();
    }

    // Intern vertices within snap_tol.
    let mut verts: Vec<V2> = Vec::new();
    let intern = |p: V2, verts: &mut Vec<V2>| -> u32 {
        for (i, &q) in verts.iter().enumerate() {
            if dist2(p, q) <= snap_tol {
                return i as u32;
            }
        }
        verts.push(p);
        (verts.len() - 1) as u32
    };
    let vids: Vec<u32> = pts.iter().map(|&p| intern(p, &mut verts)).collect();
    let mut raw_segs: Vec<(u32, u32)> = Vec::new();
    for i in 0..vids.len() {
        let (u, v) = (vids[i], vids[(i + 1) % vids.len()]);
        if u != v {
            raw_segs.push((u, v));
        }
    }

    // Split any directed segment at vertices lying strictly on its interior (multi-hole scanline cutlines).
    let mut split_segs: Vec<(u32, u32)> = Vec::new();
    for (u, v) in raw_segs {
        let (a, b) = (verts[u as usize], verts[v as usize]);
        let ab = sub2(b, a);
        let l2 = dot2(ab, ab);
        let mut mids: Vec<(f64, u32)> = Vec::new();
        for (wi, &w) in verts.iter().enumerate() {
            let wi = wi as u32;
            if wi == u || wi == v {
                continue;
            }
            if point_seg_dist(w, a, b) <= snap_tol {
                let t = dot2(sub2(w, a), ab) / l2;
                if t > 0.0 && t < 1.0 {
                    mids.push((t, wi));
                }
            }
        }
        mids.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        mids.dedup_by_key(|x| x.1);
        let mut prev = u;
        for (_, m) in mids {
            if prev != m {
                split_segs.push((prev, m));
                prev = m;
            }
        }
        if prev != v {
            split_segs.push((prev, v));
        }
    }

    // Cancel opposite directed cutline edges (u, v) against (v, u).
    let mut counts: BTreeMap<(u32, u32), i32> = BTreeMap::new();
    for (u, v) in split_segs {
        if let Some(opp) = counts.get_mut(&(v, u)) {
            if *opp > 0 {
                *opp -= 1;
                if *opp == 0 {
                    counts.remove(&(v, u));
                }
                continue;
            }
        }
        *counts.entry((u, v)).or_insert(0) += 1;
    }

    #[derive(Clone)]
    struct Half {
        to: u32,
        ang: f64,
        used: bool,
    }
    let mut adj: Vec<Vec<Half>> = vec![Vec::new(); verts.len()];
    for ((u, v), cnt) in counts {
        let d = sub2(verts[v as usize], verts[u as usize]);
        let ang = d[1].atan2(d[0]);
        for _ in 0..cnt {
            adj[u as usize].push(Half { to: v, ang, used: false });
        }
    }
    for list in &mut adj {
        list.sort_by(|a, b| a.ang.total_cmp(&b.ang).then(a.to.cmp(&b.to)));
    }

    let mut raw_cycles: Vec<Vec<u32>> = Vec::new();
    for start_u in 0..verts.len() as u32 {
        while let Some(first_idx) = adj[start_u as usize].iter().position(|h| !h.used) {
            adj[start_u as usize][first_idx].used = true;
            let mut cur = adj[start_u as usize][first_idx].to;
            let mut prev = start_u;
            let mut stack: Vec<u32> = vec![start_u];
            let mut pos_in_stack: HashMap<u32, usize> = HashMap::from([(start_u, 0)]);

            loop {
                if let Some(&idx_seen) = pos_in_stack.get(&cur) {
                    let cyc: Vec<u32> = stack[idx_seen..].to_vec();
                    for &rem in &stack[idx_seen + 1..] {
                        pos_in_stack.remove(&rem);
                    }
                    stack.truncate(idx_seen + 1);
                    if cyc.len() >= 3 {
                        raw_cycles.push(cyc);
                    }
                    if stack.len() == 1 && cur == start_u {
                        break;
                    }
                    prev = if stack.len() >= 2 { stack[stack.len() - 2] } else { cur };
                } else {
                    pos_in_stack.insert(cur, stack.len());
                    stack.push(cur);
                }

                // Pick next unused outgoing edge turning most clockwise from backward vector (cur -> prev).
                let back = sub2(verts[prev as usize], verts[cur as usize]);
                let back_ang = back[1].atan2(back[0]);
                let next_edge = adj[cur as usize]
                    .iter()
                    .enumerate()
                    .filter(|(_, h)| !h.used)
                    .min_by(|(_, a), (_, b)| {
                        let da = {
                            let d = (back_ang - a.ang).rem_euclid(TAU);
                            if d < 1e-13 {
                                TAU
                            } else {
                                d
                            }
                        };
                        let db = {
                            let d = (back_ang - b.ang).rem_euclid(TAU);
                            if d < 1e-13 {
                                TAU
                            } else {
                                d
                            }
                        };
                        da.total_cmp(&db)
                    })
                    .map(|(i, _)| i);

                let Some(ni) = next_edge else { break };
                adj[cur as usize][ni].used = true;
                prev = cur;
                cur = adj[cur as usize][ni].to;
            }
        }
    }

    let area_eps = 1e-24 * diag * diag;
    let mut outers: Vec<Vec<V2>> = Vec::new();
    let mut holes: Vec<Vec<V2>> = Vec::new();
    for cyc in raw_cycles {
        let pts_c: Vec<V2> = cyc.into_iter().map(|v| verts[v as usize]).collect();
        let clean = strip_collinear(&pts_c, snap_tol);
        if clean.len() < 3 {
            continue;
        }
        let a = polygon_area(&clean);
        if a > area_eps {
            outers.push(canonicalize_ring(clean, true));
        } else if a < -area_eps {
            holes.push(canonicalize_ring(clean, false));
        }
    }

    if outers.is_empty() {
        return Vec::new();
    }
    let outer_areas: Vec<f64> = outers.iter().map(|o| shoelace_area(o)).collect();
    let mut grouped: Vec<IslandRings> = outers.into_iter().map(|o| (o, Vec::new())).collect();
    for h in holes {
        let sample = h[0];
        let mut cands: Vec<usize> =
            (0..grouped.len()).filter(|&oi| point_in_or_on_ring(&grouped[oi].0, sample, snap_tol)).collect();
        cands.sort_by(|&a, &b| outer_areas[a].total_cmp(&outer_areas[b]));
        if let Some(&best) = cands.first() {
            grouped[best].1.push(h);
        }
    }
    grouped
}

/// Pre-snaps vertices and T-junctions within `snap_r` across all rings of a layer so near-coincident
/// GDSII edges/corners from rotated `SREF`/`AREF` instances coincide cleanly inside `geom::arrange`.
fn snap_layer_pieces(pieces: Vec<IslandRings>, diag: f64) -> Vec<IslandRings> {
    let band = crate::arrange::GRAZE_K * 1e-9 * diag;
    let snap_r = 2.0 * band;
    let mut canon: Vec<V2> = Vec::new();
    let intern_pt = |p: V2, canon: &mut Vec<V2>| -> V2 {
        for &q in canon.iter() {
            if dist2(p, q) <= snap_r {
                return q;
            }
        }
        canon.push(p);
        p
    };

    let snap_ring_verts = |ring: &[V2], canon: &mut Vec<V2>| -> Vec<V2> {
        let snapped: Vec<V2> = ring.iter().map(|&p| intern_pt(p, canon)).collect();
        clean_ring(&snapped, snap_r)
    };

    let mut snapped_pieces: Vec<IslandRings> = Vec::new();
    for (outer, holes) in pieces {
        let o = snap_ring_verts(&outer, &mut canon);
        if o.len() < 3 {
            continue;
        }
        let o_perim: f64 = (0..o.len()).map(|i| dist2(o[i], o[(i + 1) % o.len()])).sum();
        if o_perim <= 0.0 || 2.0 * shoelace_area(&o) / o_perim <= 1.2 * band {
            continue;
        }
        let mut hs = Vec::new();
        for h in holes {
            let hc = snap_ring_verts(&h, &mut canon);
            if hc.len() < 3 {
                continue;
            }
            let h_perim: f64 = (0..hc.len()).map(|i| dist2(hc[i], hc[(i + 1) % hc.len()])).sum();
            if h_perim > 0.0 && 2.0 * shoelace_area(&hc) / h_perim > 1.2 * band {
                hs.push(hc);
            }
        }
        snapped_pieces.push((o, hs));
    }

    // Insert canonical vertices lying within `snap_r` of any open segment so T-junctions coincide.
    let split_ring_tjunctions = |ring: &[V2], canon: &[V2]| -> Vec<V2> {
        let mut out = Vec::with_capacity(ring.len());
        let n = ring.len();
        for i in 0..n {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            out.push(a);
            let ab = sub2(b, a);
            let l2 = dot2(ab, ab);
            if l2 <= snap_r * snap_r {
                continue;
            }
            let mut mids: Vec<(f64, V2)> = Vec::new();
            for &v in canon {
                if dist2(v, a) <= snap_r || dist2(v, b) <= snap_r {
                    continue;
                }
                if point_seg_dist(v, a, b) <= snap_r {
                    let t = dot2(sub2(v, a), ab) / l2;
                    if t > 0.0 && t < 1.0 {
                        mids.push((t, v));
                    }
                }
            }
            mids.sort_by(|x, y| x.0.total_cmp(&y.0));
            for (_, v) in mids {
                if out.last().is_none_or(|&prev| dist2(prev, v) > snap_r) && dist2(v, b) > snap_r {
                    out.push(v);
                }
            }
        }
        clean_ring(&out, snap_r)
    };

    snapped_pieces
        .into_iter()
        .filter_map(|(o, hs)| {
            let o2 = split_ring_tjunctions(&o, &canon);
            if o2.len() < 3 {
                return None;
            }
            let hs2: Vec<Vec<V2>> =
                hs.into_iter().map(|h| split_ring_tjunctions(&h, &canon)).filter(|h| h.len() >= 3).collect();
            Some((o2, hs2))
        })
        .collect()
}

fn loop_to_ring(cx: &Complex, li: usize, col_tol: f64, want_ccw: bool) -> Vec<V2> {
    let lp = &cx.loops[li];
    let mut pts = Vec::with_capacity(lp.edges.len());
    for &(ei, fwd) in &lp.edges {
        let ed = &cx.edges[ei as usize];
        let vid = if fwd { ed.v0 } else { ed.v1 };
        let p3 = cx.vertices[vid as usize].p;
        pts.push([p3[0], p3[1]]);
    }
    let clean = strip_collinear(&pts, col_tol);
    canonicalize_ring(clean, want_ccw)
}

/// Boolean-unions all decomposed polygons of a layer via `geom::arrange` and returns disconnected
/// `(outer_ccw, holes_cw)` islands.
fn union_layer_polygons(raw_polys: &[Vec<V2>]) -> GeomResult<Vec<IslandRings>> {
    let mut pieces: Vec<IslandRings> = Vec::new();
    let mut bb = Aabb::EMPTY;
    for raw in raw_polys {
        for (o, hs) in decompose_raw_polygon(raw) {
            for &p in &o {
                bb.grow([p[0], p[1], 0.0]);
            }
            pieces.push((o, hs));
        }
    }
    if pieces.is_empty() {
        return Ok(Vec::new());
    }
    let diag = bb.diagonal();
    if !diag.is_finite() || diag <= 1e-30 {
        return Ok(Vec::new());
    }

    let pieces = snap_layer_pieces(pieces, diag);
    if pieces.is_empty() {
        return Ok(Vec::new());
    }

    let mut m = Model::new_2d();
    let mut piece_ids = Vec::with_capacity(pieces.len());
    for (pi, (outer, holes)) in pieces.iter().enumerate() {
        let o_pts: Vec<(f64, f64)> = outer.iter().map(|p| (p[0], p[1])).collect();
        let out_id = m.polygon(&format!("p{pi}"), &o_pts);
        if holes.is_empty() {
            piece_ids.push(out_id);
        } else {
            let mut h_ids = Vec::with_capacity(holes.len());
            for (hi, h) in holes.iter().enumerate() {
                let h_pts: Vec<(f64, f64)> = h.iter().map(|p| (p[0], p[1])).collect();
                h_ids.push(m.polygon(&format!("p{pi}_h{hi}"), &h_pts));
            }
            piece_ids.push(m.difference(&format!("p{pi}_d"), out_id, &h_ids));
        }
    }
    let union_id = if piece_ids.len() == 1 { piece_ids[0] } else { m.union("u", &piece_ids) };
    m.overlay(&[("layer", union_id)]);

    let mut cache = Cache::new();
    let scene = build_scene(&m.list, &mut cache)?;
    let cx = arrange(&scene, &mut cache)?;
    if cx.regions.is_empty() {
        return Ok(Vec::new());
    }

    let col_tol = 1e-12 * diag;
    let reg = &cx.regions[0];
    let mut outer_loop_map: HashMap<u32, usize> = HashMap::new();
    let mut islands: Vec<IslandRings> = Vec::new();

    for &li in &reg.loops {
        if cx.loops[li as usize].area > 0.0 {
            let outer = loop_to_ring(&cx, li as usize, col_tol, true);
            if outer.len() >= 3 {
                outer_loop_map.insert(li, islands.len());
                islands.push((outer, Vec::new()));
            }
        }
    }
    for &li in &reg.loops {
        if cx.loops[li as usize].area < 0.0 {
            if let Some(parent_li) = cx.loops[li as usize].outer {
                if let Some(&isl_idx) = outer_loop_map.get(&parent_li) {
                    let hole = loop_to_ring(&cx, li as usize, col_tol, false);
                    if hole.len() >= 3 {
                        islands[isl_idx].1.push(hole);
                    }
                }
            }
        }
    }

    Ok(islands)
}

fn ring_area_and_centroid(ring: &[V2]) -> (f64, V2) {
    let n = ring.len();
    if n < 3 {
        return (0.0, [0.0, 0.0]);
    }
    // Shift to ring[0] for translation-invariant floating-point accuracy on distant coordinates.
    let o = ring[0];
    let mut a2 = 0.0;
    let mut cx6 = 0.0;
    let mut cy6 = 0.0;
    for i in 0..n {
        let p = sub2(ring[i], o);
        let q = sub2(ring[(i + 1) % n], o);
        let c = cross2(p, q);
        a2 += c;
        cx6 += (p[0] + q[0]) * c;
        cy6 += (p[1] + q[1]) * c;
    }
    let area = 0.5 * a2.abs();
    if a2.abs() <= 1e-300 {
        let avg = ring.iter().fold([0.0, 0.0], |acc, &p| add2(acc, p));
        return (0.0, mul2(avg, 1.0 / n as f64));
    }
    (area, [o[0] + cx6 / (3.0 * a2), o[1] + cy6 / (3.0 * a2)])
}

fn island_metrics(outer: &[V2], holes: &[Vec<V2>]) -> (f64, V2, [f64; 4]) {
    let (oa, oc) = ring_area_and_centroid(outer);
    let mut net_a = oa;
    let mut wx = oa * oc[0];
    let mut wy = oa * oc[1];
    for h in holes {
        let (ha, hc) = ring_area_and_centroid(h);
        net_a -= ha;
        wx -= ha * hc[0];
        wy -= ha * hc[1];
    }
    let centroid = if net_a > 1e-300 { [wx / net_a, wy / net_a] } else { oc };
    let mut xmin = f64::INFINITY;
    let mut ymin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for &p in outer {
        xmin = xmin.min(p[0]);
        ymin = ymin.min(p[1]);
        xmax = xmax.max(p[0]);
        ymax = ymax.max(p[1]);
    }
    (net_a.max(0.0), centroid, [xmin, ymin, xmax, ymax])
}

fn compare_area_centroid(a_area: f64, a_c: V2, b_area: f64, b_c: V2, extent: f64) -> std::cmp::Ordering {
    let area_scale = a_area.max(b_area).max(1e-30);
    if (a_area - b_area).abs() > 1e-12 * area_scale {
        return b_area.total_cmp(&a_area); // descending area
    }
    let pos_tol = 1e-12 * extent.max(1e-12);
    if (a_c[0] - b_c[0]).abs() > pos_tol {
        return a_c[0].total_cmp(&b_c[0]); // ascending x
    }
    a_c[1].total_cmp(&b_c[1]) // ascending y
}

/// Sanitizes a GDSII `TEXT` pin label into a valid `geom` feature name.
pub fn sanitize_pin_name(raw: &str) -> String {
    const FORBIDDEN: &str = ".|/,()[]&!:=\"'<>*";
    let mut out = String::with_capacity(raw.len());
    for c in raw.trim().chars() {
        if c.is_control() || c.is_whitespace() || FORBIDDEN.contains(c) {
            if !out.ends_with('_') {
                out.push('_');
            }
        } else {
            out.push(c);
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "pin".to_string()
    } else {
        trimmed.to_string()
    }
}

fn unique_sheet_name(base: &str, used: &mut HashSet<String>) -> String {
    if used.insert(base.to_string()) {
        return base.to_string();
    }
    for k in 1usize.. {
        let cand = format!("{base}_{k}");
        if used.insert(cand.clone()) {
            return cand;
        }
    }
    unreachable!()
}

/// Parses and imports a binary GDSII stream according to `opts`.
pub fn import_gds_bytes(bytes: &[u8], opts: &GdsImportOptions) -> GeomResult<GdsImportResult> {
    let lib = parse_gds(bytes)?;
    import_gds_library(&lib, opts)
}

/// Flattens `lib`, boolean-unions requested layers, splits/orders islands, applies optional
/// `simplify_tol`, and matches `TEXT` pins.
pub fn import_gds_library(lib: &GdsLibrary, opts: &GdsImportOptions) -> GeomResult<GdsImportResult> {
    if opts.layers.is_empty() {
        return Err(GeomError("import_gds requires at least one entry in `layers`".into()));
    }
    if opts.simplify_tol < 0.0 || !opts.simplify_tol.is_finite() {
        return Err(GeomError(format!("simplify_tol ({}) must be non-negative and finite", opts.simplify_tol)));
    }
    if opts.arc_tol < 0.0 || !opts.arc_tol.is_finite() {
        return Err(GeomError(format!("arc_tol ({}) must be non-negative and finite", opts.arc_tol)));
    }
    for (key, _) in &opts.layers {
        check_name("layer", key)?;
    }

    let req_set: HashSet<(i16, i16)> = opts.layers.iter().map(|(_, ld)| *ld).collect();
    let top_idx = resolve_top_cell(lib, opts.cell.as_deref(), &req_set)?;
    let top_name = lib.cells[top_idx].name.clone();
    let cell_map: HashMap<&str, usize> = lib.cells.iter().enumerate().map(|(i, c)| (c.name.as_str(), i)).collect();

    let mut world = FlatWorld::default();
    flatten_rec(lib, &cell_map, top_idx, GdsXf::ID, 0, opts.pin_layer, &mut world)?;

    // Deterministically sort flattened pins.
    world
        .pins
        .sort_by(|a, b| a.label.cmp(&b.label).then(a.pos[0].total_cmp(&b.pos[0])).then(a.pos[1].total_cmp(&b.pos[1])));
    world.pins.dedup_by(|a, b| a.label == b.label && dist2(a.pos, b.pos) <= 1e-15);

    // 1. Boolean-union each requested layer and sort islands deterministically.
    struct RawLayer {
        key: String,
        islands: Vec<(Vec<V2>, Vec<Vec<V2>>)>,
    }
    let mut raw_layers: Vec<RawLayer> = Vec::with_capacity(opts.layers.len());
    let mut layer_aabb = Aabb::EMPTY;

    for (key, ld) in &opts.layers {
        let mut islands = match world.polys.get(ld) {
            Some(polys) => union_layer_polygons(polys)?,
            None => Vec::new(),
        };
        let mut ext_bb = Aabb::EMPTY;
        for (o, _) in &islands {
            for &p in o {
                ext_bb.grow([p[0], p[1], 0.0]);
                layer_aabb.grow([p[0], p[1], 0.0]);
            }
        }
        let extent = ext_bb.diagonal();
        // Sort holes inside each island deterministically (descending hole area, then centroid x, y).
        for (_, hs) in &mut islands {
            hs.sort_by(|a, b| {
                let (aa, ac) = ring_area_and_centroid(a);
                let (ba, bc) = ring_area_and_centroid(b);
                compare_area_centroid(aa, ac, ba, bc, extent)
            });
        }
        // Sort islands deterministically (descending net area, then centroid x, y).
        islands.sort_by(|(ao, ahs), (bo, bhs)| {
            let (aa, ac, _) = island_metrics(ao, ahs);
            let (ba, bc, _) = island_metrics(bo, bhs);
            compare_area_centroid(aa, ac, ba, bc, extent)
        });
        raw_layers.push(RawLayer { key: key.clone(), islands });
    }

    if raw_layers.iter().all(|l| l.islands.is_empty()) {
        let specs: Vec<String> = opts.layers.iter().map(|(k, (l, d))| format!("{k}: ({l}, {d})")).collect();
        return Err(GeomError(format!(
            "none of the requested GDSII layers [{}] contain any polygon geometry in cell '{top_name}'",
            specs.join(", ")
        )));
    }

    // 2. Optional RDP simplification with cross-island / cross-layer shared-vertex pinning.
    // When `arc_tol > 0.0`, preserve raw circular-arc vertices so `Op::Polygons` fits exact arcs.
    if opts.simplify_tol > 0.0 && opts.arc_tol == 0.0 {
        let all_island_pts: Vec<Vec<V2>> = raw_layers
            .iter()
            .flat_map(|l| {
                l.islands
                    .iter()
                    .map(|(o, hs)| o.iter().copied().chain(hs.iter().flatten().copied()).collect::<Vec<V2>>())
            })
            .collect();
        let mut flat_idx = 0usize;
        for l in &mut raw_layers {
            for (o, hs) in &mut l.islands {
                let external_pins: Vec<V2> = all_island_pts
                    .iter()
                    .enumerate()
                    .filter(|&(idx, _)| idx != flat_idx)
                    .flat_map(|(_, pts)| pts.iter().copied())
                    .collect();
                let (mut sim_o, mut sim_hs) =
                    simplify_rings_with_pins(std::slice::from_ref(o), hs, &external_pins, opts.simplify_tol, 0.0)?;
                *o = canonicalize_ring(sim_o.remove(0), true);
                for h in &mut sim_hs {
                    let taken = std::mem::take(h);
                    *h = canonicalize_ring(taken, false);
                }
                *hs = sim_hs;
                flat_idx += 1;
            }
        }
    }

    // 3. Match TEXT pins and construct output sheets.
    let chip_diag = layer_aabb.diagonal().max(1e-9);
    let pin_tol = 1e-9 * chip_diag;
    let island_contains = |o: &[V2], hs: &[Vec<V2>], p: V2| -> bool {
        point_in_or_on_ring(o, p, pin_tol) && !hs.iter().any(|h| point_strictly_in_ring(h, p, pin_tol))
    };

    let mut used_names: HashSet<String> = HashSet::new();
    let mut sheets: Vec<GdsSheet> = Vec::new();
    let mut layers_out: Vec<(String, Vec<String>)> = Vec::with_capacity(raw_layers.len());
    let mut pins_map: BTreeMap<String, String> = BTreeMap::new();

    for rl in raw_layers {
        let n_isl = rl.islands.len();
        if n_isl == 0 {
            layers_out.push((rl.key, Vec::new()));
            continue;
        }
        let mut layer_sheet_names = Vec::new();

        if opts.split_islands {
            for (idx, (outer, holes)) in rl.islands.into_iter().enumerate() {
                let matched_labels: Vec<String> = world
                    .pins
                    .iter()
                    .filter(|pin| island_contains(&outer, &holes, pin.pos))
                    .map(|pin| pin.label.clone())
                    .collect();

                let base_name = if let Some(first_pin) = matched_labels.first() {
                    sanitize_pin_name(first_pin)
                } else if n_isl == 1 {
                    rl.key.clone()
                } else {
                    format!("{}_{idx}", rl.key)
                };
                let sheet_name = unique_sheet_name(&base_name, &mut used_names);
                for lbl in matched_labels {
                    pins_map.entry(lbl).or_insert_with(|| sheet_name.clone());
                }
                let (area, centroid, bbox) = island_metrics(&outer, &holes);
                layer_sheet_names.push(sheet_name.clone());
                sheets.push(GdsSheet {
                    name: sheet_name,
                    layer_key: rl.key.clone(),
                    outers: vec![outer.clone()],
                    holes: holes.clone(),
                    islands: vec![(outer, holes)],
                    area,
                    centroid,
                    bbox,
                });
            }
        } else {
            let sheet_name = unique_sheet_name(&rl.key, &mut used_names);
            let mut outers = Vec::with_capacity(n_isl);
            let mut holes = Vec::new();
            let mut total_area = 0.0;
            let mut wx = 0.0;
            let mut wy = 0.0;
            let mut s_bbox = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];

            for (o, hs) in &rl.islands {
                for pin in &world.pins {
                    if island_contains(o, hs, pin.pos) {
                        pins_map.entry(pin.label.clone()).or_insert_with(|| sheet_name.clone());
                    }
                }
                let (ia, ic, ib) = island_metrics(o, hs);
                total_area += ia;
                wx += ia * ic[0];
                wy += ia * ic[1];
                s_bbox[0] = s_bbox[0].min(ib[0]);
                s_bbox[1] = s_bbox[1].min(ib[1]);
                s_bbox[2] = s_bbox[2].max(ib[2]);
                s_bbox[3] = s_bbox[3].max(ib[3]);
                outers.push(o.clone());
                holes.extend(hs.iter().cloned());
            }
            let centroid = if total_area > 1e-300 { [wx / total_area, wy / total_area] } else { [0.0, 0.0] };
            layer_sheet_names.push(sheet_name.clone());
            sheets.push(GdsSheet {
                name: sheet_name,
                layer_key: rl.key.clone(),
                outers,
                holes,
                islands: rl.islands,
                area: total_area,
                centroid,
                bbox: s_bbox,
            });
        }

        layers_out.push((rl.key, layer_sheet_names));
    }

    let layer_bbox = [layer_aabb.lo[0], layer_aabb.lo[1], layer_aabb.hi[0], layer_aabb.hi[1]];
    let cell_bbox = if world.all_bbox.is_empty() {
        layer_bbox
    } else {
        [world.all_bbox.lo[0], world.all_bbox.lo[1], world.all_bbox.hi[0], world.all_bbox.hi[1]]
    };

    Ok(GdsImportResult {
        cell: top_name,
        dbu_m: lib.dbu_m,
        sheets,
        layers: layers_out,
        layer_bbox,
        cell_bbox,
        pins: pins_map.into_iter().collect(),
    })
}
