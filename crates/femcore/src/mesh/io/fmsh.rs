//! `.fmsh`: a compact, lossless, self-checking binary format for mixed-cell meshes and results.
//!
//! One file holds the mesh and, optionally, any number of result steps (time or frequency) with
//! named point or cell fields, so nothing else is needed next to it and the mesh is stored once.
//! Every array is lossless (floats round-trip bit for bit, including NaN payloads and `-0.0`) and
//! is stored with whichever coding is smallest. A CRC-32 catches truncation and corruption.
//!
//! Layout (little endian; `varint` = LEB128, `zigzag` = signed to unsigned):
//!
//! ```text
//! "FMSH" version:u8(=2) D:u8 | n_points n_cells n_facets n_steps (varint)
//! points   : per axis an f64 array (below)
//! cells    : block(run-length (cell type code, count)); block(connectivity: first node as zigzag
//!            delta against the previous cell's first node, other nodes as zigzag deltas against the
//!            cell's first node); block(run-length labels)
//! facets   : block(run-length node counts); block(connectivity as above); block(run-length labels)
//! steps    : per step value:f64, n_fields (varint), then per field: name (varint length + UTF-8),
//!            location (0 point, 1 cell), components:u8, f64 array of components * count values
//! crc32    : u32 over everything before it
//!
//! block    : tag:u8 (0 raw, 1 order-0 range coder, 2 order-1 range coder) raw_len (varint)
//!            [comp_len (varint)] bytes          -- the smallest of the three is kept
//! f64 array: flag:u8 then, whichever is smallest after block coding:
//!            0 byte planes (plane j holds byte j of every value; the exponent planes of smooth or
//!              structured data are nearly constant and collapse); 1 Gorilla-style XOR with the
//!              previous value (control byte 0xFF for "same", else lead_zero_bytes << 3 | (len-1),
//!              followed by the `len` significant bytes); 2 dictionary (count, distinct values in
//!              order of first appearance as byte planes, then one varint index per value)
//! ```
//!
//! Unstructured coordinates and solution values are essentially random in their low mantissa
//! bytes, so the gain there is bounded by the sign/exponent/top-mantissa bytes (about 15-25 %);
//! indices, labels, cell types and structured coordinates compress much better. Nothing is lossy.

use super::entropy::{self, Model};
use crate::meshers::cellmesh::{CellMesh, CellType};

const MAGIC: &[u8; 4] = b"FMSH";
const VERSION: u8 = 2;

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8 & 0x7F) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn zigzag(v: i64) -> u64 {
    ((v << 1) ^ (v >> 63)) as u64
}

fn unzigzag(u: u64) -> i64 {
    ((u >> 1) as i64) ^ -((u & 1) as i64)
}

/// CRC-32 (IEEE, reflected, polynomial 0xEDB88320), bitwise implementation.
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            crc = if crc & 1 != 0 { (crc >> 1) ^ 0xEDB8_8320 } else { crc >> 1 };
        }
    }
    !crc
}

struct Reader<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Reader<'a> {
    fn byte(&mut self) -> Result<u8, String> {
        let v = *self.b.get(self.i).ok_or("unexpected end of data")?;
        self.i += 1;
        Ok(v)
    }
    fn varint(&mut self) -> Result<u64, String> {
        let (mut v, mut shift) = (0u64, 0);
        loop {
            let b = self.byte()?;
            v |= ((b & 0x7F) as u64) << shift;
            if b & 0x80 == 0 {
                return Ok(v);
            }
            shift += 7;
            if shift > 63 {
                return Err("varint too long".into());
            }
        }
    }
    fn bytes(&mut self, n: usize) -> Result<&'a [u8], String> {
        let s = self.b.get(self.i..self.i + n).ok_or("unexpected end of data")?;
        self.i += n;
        Ok(s)
    }
}

/// Where a field lives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    Point,
    Cell,
}

/// A named array of `components` values per point or per cell (component-major per entity:
/// `data[entity * components + c]`).
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    pub name: String,
    pub location: Location,
    pub components: usize,
    pub data: Vec<f64>,
}

impl Field {
    pub fn scalar_on_points(name: &str, data: Vec<f64>) -> Self {
        Field { name: name.into(), location: Location::Point, components: 1, data }
    }
    pub fn scalar_on_cells(name: &str, data: Vec<f64>) -> Self {
        Field { name: name.into(), location: Location::Cell, components: 1, data }
    }
}

/// Fields at one time / frequency / iteration value.
#[derive(Clone, Debug, PartialEq)]
pub struct Step {
    pub value: f64,
    pub fields: Vec<Field>,
}

/// A mesh with optional result steps.
#[derive(Clone, Debug)]
pub struct Dataset<const D: usize> {
    pub mesh: CellMesh<D>,
    pub steps: Vec<Step>,
}

fn put_block(out: &mut Vec<u8>, data: &[u8]) {
    let o0 = entropy::compress(data, Model::Order0);
    let o1 = entropy::compress(data, Model::Order1);
    let (tag, best) = if o1.len() <= o0.len() { (2u8, &o1) } else { (1u8, &o0) };
    if best.len() + 2 < data.len() {
        out.push(tag);
        put_varint(out, data.len() as u64);
        put_varint(out, best.len() as u64);
        out.extend_from_slice(best);
    } else {
        out.push(0);
        put_varint(out, data.len() as u64);
        out.extend_from_slice(data);
    }
}

fn get_block(r: &mut Reader) -> Result<Vec<u8>, String> {
    let tag = r.byte()?;
    let raw_len = r.varint()? as usize;
    match tag {
        0 => Ok(r.bytes(raw_len)?.to_vec()),
        1 | 2 => {
            let cl = r.varint()? as usize;
            let body = r.bytes(cl)?;
            // the coder cannot produce more than 8 bytes per byte of input; bound the allocation
            if raw_len > cl.saturating_mul(1 << 14).saturating_add(64) {
                return Err("implausible block length".into());
            }
            Ok(entropy::decompress(body, raw_len, if tag == 1 { Model::Order0 } else { Model::Order1 }))
        }
        t => Err(format!("unknown block tag {t}")),
    }
}

fn planes(bits: &[u64]) -> Vec<u8> {
    let n = bits.len();
    let mut out = vec![0u8; n * 8];
    for (i, b) in bits.iter().enumerate() {
        for (j, byte) in b.to_le_bytes().into_iter().enumerate() {
            out[j * n + i] = byte;
        }
    }
    out
}

fn unplane(bytes: &[u8], n: usize) -> Result<Vec<u64>, String> {
    if bytes.len() != n * 8 {
        return Err("byte-plane length mismatch".into());
    }
    Ok((0..n).map(|i| u64::from_le_bytes(std::array::from_fn(|j| bytes[j * n + i]))).collect())
}

fn put_f64s(out: &mut Vec<u8>, vals: &[f64]) {
    let bits: Vec<u64> = vals.iter().map(|v| v.to_bits()).collect();
    // 0: byte planes
    let mut a = vec![0u8];
    put_block(&mut a, &planes(&bits));
    // 1: XOR with previous
    let mut g = Vec::new();
    let mut prev = 0u64;
    for &b in &bits {
        let x = b ^ prev;
        prev = b;
        if x == 0 {
            g.push(0xFF);
        } else {
            let lz = (x.leading_zeros() / 8) as u8;
            let tz = (x.trailing_zeros() / 8) as u8;
            let len = 8 - lz - tz;
            g.push((lz << 3) | (len - 1));
            g.extend_from_slice(&(x >> (8 * tz as u32)).to_le_bytes()[..len as usize]);
        }
    }
    let mut b = vec![1u8];
    put_block(&mut b, &g);
    // 2: dictionary (only worth building when values repeat)
    let mut best = if b.len() < a.len() { b } else { a };
    let mut index: std::collections::HashMap<u64, u64> = std::collections::HashMap::new();
    let mut table: Vec<u64> = Vec::new();
    let mut idx: Vec<u8> = Vec::new();
    for &bt in &bits {
        let k = *index.entry(bt).or_insert_with(|| {
            table.push(bt);
            table.len() as u64 - 1
        });
        put_varint(&mut idx, k);
        if table.len() * 2 > bits.len() + 16 {
            break;
        }
    }
    if table.len() * 2 <= bits.len() + 16 {
        let mut d = vec![2u8];
        put_varint(&mut d, table.len() as u64);
        put_block(&mut d, &planes(&table));
        put_block(&mut d, &idx);
        if d.len() < best.len() {
            best = d;
        }
    }
    out.extend_from_slice(&best);
}

fn get_f64s(r: &mut Reader, n: usize) -> Result<Vec<f64>, String> {
    match r.byte()? {
        0 => Ok(unplane(&get_block(r)?, n)?.into_iter().map(f64::from_bits).collect()),
        1 => {
            let stream = get_block(r)?;
            let mut s = Reader { b: &stream, i: 0 };
            let mut out = Vec::with_capacity(n);
            let mut prev = 0u64;
            for _ in 0..n {
                let c = s.byte()?;
                let x = if c == 0xFF {
                    0
                } else {
                    let (lz, len) = ((c >> 3) as u32, (c & 7) as usize + 1);
                    if lz as usize + len > 8 {
                        return Err("corrupt control byte".into());
                    }
                    let tz = 8 - lz - len as u32;
                    let mut buf = [0u8; 8];
                    buf[..len].copy_from_slice(s.bytes(len)?);
                    u64::from_le_bytes(buf) << (8 * tz)
                };
                prev ^= x;
                out.push(f64::from_bits(prev));
            }
            if s.i != stream.len() {
                return Err("trailing bytes in coordinate stream".into());
            }
            Ok(out)
        }
        2 => {
            let nt = r.varint()? as usize;
            if nt > n {
                return Err("dictionary larger than array".into());
            }
            let table = unplane(&get_block(r)?, nt)?;
            let stream = get_block(r)?;
            let mut s = Reader { b: &stream, i: 0 };
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                let k = s.varint()? as usize;
                out.push(f64::from_bits(*table.get(k).ok_or("dictionary index out of range")?));
            }
            Ok(out)
        }
        f => Err(format!("unknown f64 array flag {f}")),
    }
}

fn put_runs(out: &mut Vec<u8>, vals: impl Iterator<Item = i64>) {
    let mut runs: Vec<(i64, u64)> = Vec::new();
    for v in vals {
        match runs.last_mut() {
            Some((x, n)) if *x == v => *n += 1,
            _ => runs.push((v, 1)),
        }
    }
    let mut s = Vec::new();
    put_varint(&mut s, runs.len() as u64);
    for (v, n) in runs {
        put_varint(&mut s, zigzag(v));
        put_varint(&mut s, n);
    }
    put_block(out, &s);
}

fn get_runs(r: &mut Reader, n: usize) -> Result<Vec<i64>, String> {
    let stream = get_block(r)?;
    let mut s = Reader { b: &stream, i: 0 };
    let nr = s.varint()? as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..nr {
        let v = unzigzag(s.varint()?);
        let c = s.varint()? as usize;
        if out.len() + c > n {
            return Err("run-length overflow".into());
        }
        out.extend(std::iter::repeat_n(v, c));
    }
    if out.len() != n || s.i != stream.len() {
        return Err("run-length size mismatch".into());
    }
    Ok(out)
}

fn put_conn(out: &mut Vec<u8>, conn: &[u32], offsets: &[u32]) {
    let mut s = Vec::new();
    let mut prev_first = 0i64;
    for w in offsets.windows(2) {
        let vs = &conn[w[0] as usize..w[1] as usize];
        let first = vs[0] as i64;
        put_varint(&mut s, zigzag(first - prev_first));
        prev_first = first;
        for &v in &vs[1..] {
            put_varint(&mut s, zigzag(v as i64 - first));
        }
    }
    put_block(out, &s);
}

fn get_conn(r: &mut Reader, sizes: &[usize], n_points: usize) -> Result<Vec<u32>, String> {
    let stream = get_block(r)?;
    let mut s = Reader { b: &stream, i: 0 };
    let mut conn = Vec::with_capacity(sizes.iter().sum());
    let mut prev_first = 0i64;
    for &sz in sizes {
        let first = prev_first + unzigzag(s.varint()?);
        prev_first = first;
        conn.push(first as u32);
        for _ in 1..sz {
            conn.push((first + unzigzag(s.varint()?)) as u32);
        }
    }
    if s.i != stream.len() {
        return Err("trailing bytes in connectivity".into());
    }
    if conn.iter().any(|&v| v as usize >= n_points) {
        return Err("vertex id out of range".into());
    }
    Ok(conn)
}

fn offsets_of(sizes: &[usize]) -> Vec<u32> {
    std::iter::once(0)
        .chain(sizes.iter().scan(0u32, |s, &n| {
            *s += n as u32;
            Some(*s)
        }))
        .collect()
}

/// Encode a mesh without results.
pub fn encode<const D: usize>(m: &CellMesh<D>) -> Vec<u8> {
    encode_steps(m, &[])
}

/// Encode a mesh with result steps. Panics if a field's length does not match its location.
pub fn encode_steps<const D: usize>(m: &CellMesh<D>, steps: &[Step]) -> Vec<u8> {
    let mut o = Vec::new();
    o.extend_from_slice(MAGIC);
    o.push(VERSION);
    o.push(D as u8);
    put_varint(&mut o, m.n_points() as u64);
    put_varint(&mut o, m.n_cells() as u64);
    put_varint(&mut o, m.n_facets() as u64);
    put_varint(&mut o, steps.len() as u64);
    for a in 0..D {
        let col: Vec<f64> = m.points.iter().map(|p| p[a]).collect();
        put_f64s(&mut o, &col);
    }
    put_runs(&mut o, m.types.iter().map(|t| t.code() as i64));
    put_conn(&mut o, &m.conn, &m.offsets);
    put_runs(&mut o, m.label.iter().map(|&l| l as i64));
    put_runs(&mut o, (0..m.n_facets()).map(|f| m.facet(f).len() as i64));
    put_conn(&mut o, &m.facet_conn, &m.facet_offsets);
    put_runs(&mut o, m.facet_label.iter().map(|&l| l as i64));
    for s in steps {
        o.extend_from_slice(&s.value.to_bits().to_le_bytes());
        put_varint(&mut o, s.fields.len() as u64);
        for f in &s.fields {
            let count = if f.location == Location::Point { m.n_points() } else { m.n_cells() };
            assert_eq!(f.data.len(), count * f.components, "field {:?} has the wrong length", f.name);
            assert!((1..=255).contains(&f.components));
            put_varint(&mut o, f.name.len() as u64);
            o.extend_from_slice(f.name.as_bytes());
            o.push(if f.location == Location::Point { 0 } else { 1 });
            o.push(f.components as u8);
            put_f64s(&mut o, &f.data);
        }
    }
    let crc = crc32(&o);
    o.extend_from_slice(&crc.to_le_bytes());
    o
}

/// Decode a mesh (results, if any, are ignored), verifying the checksum.
pub fn decode<const D: usize>(bytes: &[u8]) -> Result<CellMesh<D>, String> {
    Ok(decode_dataset::<D>(bytes)?.mesh)
}

/// Decode a mesh and its result steps, verifying the checksum.
pub fn decode_dataset<const D: usize>(bytes: &[u8]) -> Result<Dataset<D>, String> {
    if bytes.len() < 4 + 2 + 4 {
        return Err("too short".into());
    }
    let (body, tail) = bytes.split_at(bytes.len() - 4);
    if crc32(body) != u32::from_le_bytes(tail.try_into().unwrap()) {
        return Err("checksum mismatch".into());
    }
    let mut r = Reader { b: body, i: 0 };
    if r.bytes(4)? != MAGIC {
        return Err("bad magic".into());
    }
    if r.byte()? != VERSION {
        return Err("unsupported version".into());
    }
    if r.byte()? as usize != D {
        return Err("dimension mismatch".into());
    }
    let (np, nc, nf, ns) = (r.varint()? as usize, r.varint()? as usize, r.varint()? as usize, r.varint()? as usize);
    let axes: Vec<Vec<f64>> = (0..D).map(|_| get_f64s(&mut r, np)).collect::<Result<_, _>>()?;
    let mut m = CellMesh::<D>::new();
    m.points = (0..np).map(|i| std::array::from_fn(|a| axes[a][i])).collect();
    m.types = get_runs(&mut r, nc)?
        .iter()
        .map(|&c| CellType::from_code(c as u8).filter(|t| t.dim() == D).ok_or("bad cell type".to_string()))
        .collect::<Result<_, _>>()?;
    let sizes: Vec<usize> = m.types.iter().map(|t| t.nv()).collect();
    m.conn = get_conn(&mut r, &sizes, np)?;
    m.offsets = offsets_of(&sizes);
    m.label = get_runs(&mut r, nc)?.into_iter().map(|v| v as i32).collect();
    let fsizes: Vec<usize> = get_runs(&mut r, nf)?.into_iter().map(|v| v as usize).collect();
    m.facet_conn = get_conn(&mut r, &fsizes, np)?;
    m.facet_offsets = offsets_of(&fsizes);
    m.facet_label = get_runs(&mut r, nf)?.into_iter().map(|v| v as i32).collect();
    let mut steps = Vec::new();
    for _ in 0..ns {
        let value = f64::from_bits(u64::from_le_bytes(r.bytes(8)?.try_into().unwrap()));
        let nfld = r.varint()? as usize;
        let mut fields = Vec::new();
        for _ in 0..nfld {
            let nl = r.varint()? as usize;
            let name = String::from_utf8(r.bytes(nl)?.to_vec()).map_err(|_| "field name is not UTF-8")?;
            let location = match r.byte()? {
                0 => Location::Point,
                1 => Location::Cell,
                _ => return Err("bad field location".into()),
            };
            let components = r.byte()? as usize;
            let count = if location == Location::Point { np } else { nc };
            let data = get_f64s(&mut r, count * components)?;
            fields.push(Field { name, location, components, data });
        }
        steps.push(Step { value, fields });
    }
    if r.i != body.len() {
        return Err("trailing bytes".into());
    }
    Ok(Dataset { mesh: m, steps })
}

/// Write a mesh to a file.
pub fn save<const D: usize>(path: &std::path::Path, m: &CellMesh<D>) -> std::io::Result<()> {
    std::fs::write(path, encode(m))
}

/// Write a mesh with result steps to a file.
pub fn save_steps<const D: usize>(path: &std::path::Path, m: &CellMesh<D>, steps: &[Step]) -> std::io::Result<()> {
    std::fs::write(path, encode_steps(m, steps))
}

/// Read a mesh from a file.
pub fn load<const D: usize>(path: &std::path::Path) -> Result<CellMesh<D>, String> {
    decode(&std::fs::read(path).map_err(|e| e.to_string())?)
}

/// Read a mesh and its result steps from a file.
pub fn load_dataset<const D: usize>(path: &std::path::Path) -> Result<Dataset<D>, String> {
    decode_dataset(&std::fs::read(path).map_err(|e| e.to_string())?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshers::fixtures::{box_tet, box_tri};

    fn same<const D: usize>(a: &CellMesh<D>, b: &CellMesh<D>) {
        assert_eq!(a.points.len(), b.points.len());
        for (p, q) in a.points.iter().zip(&b.points) {
            for k in 0..D {
                assert_eq!(p[k].to_bits(), q[k].to_bits());
            }
        }
        assert_eq!(a.types, b.types);
        assert_eq!(a.conn, b.conn);
        assert_eq!(a.offsets, b.offsets);
        assert_eq!(a.label, b.label);
        assert_eq!(a.facet_conn, b.facet_conn);
        assert_eq!(a.facet_offsets, b.facet_offsets);
        assert_eq!(a.facet_label, b.facet_label);
    }

    #[test]
    fn varint_and_zigzag_roundtrip() {
        for v in [0i64, 1, -1, 63, -64, 1 << 40, -(1 << 40), i64::MAX, i64::MIN] {
            assert_eq!(unzigzag(zigzag(v)), v);
            let mut o = Vec::new();
            put_varint(&mut o, zigzag(v));
            assert_eq!(unzigzag(Reader { b: &o, i: 0 }.varint().unwrap()), v);
        }
    }

    #[test]
    fn crc_known_value() {
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
    }

    #[test]
    fn roundtrip_2d_3d_and_mixed() {
        let a = CellMesh::from_simplex(&box_tri([[0.0, 1.0], [-2.0, 3.5]], 17, 9));
        same(&a, &decode::<2>(&encode(&a)).unwrap());
        let b = CellMesh::from_simplex(&box_tet([[0.0, 1.0], [0.0, 2.0], [0.0, 3.0]], [4, 3, 2]));
        same(&b, &decode::<3>(&encode(&b)).unwrap());
        let mut c = CellMesh::<3>::new();
        c.points = vec![
            [0.1, 0.2, 0.3],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
        ];
        c.push_cell(CellType::Hex, &[0, 1, 2, 3, 4, 5, 6, 7], -4);
        c.push_cell(CellType::Prism, &[1, 8, 2, 5, 9, 6], 12);
        c.push_facet(&[0, 1, 2, 3], 5);
        c.push_facet(&[1, 8, 2], 6);
        same(&c, &decode::<3>(&encode(&c)).unwrap());
    }

    #[test]
    fn special_floats_survive() {
        let mut m = CellMesh::<2>::new();
        m.points = vec![[0.0, -0.0], [f64::MIN_POSITIVE, 1e300], [f64::INFINITY, f64::NAN], [0.1, 0.1]];
        m.push_cell(CellType::Tri, &[0, 1, 3], 1);
        same(&m, &decode::<2>(&encode(&m)).unwrap());
    }

    #[test]
    fn detects_corruption_truncation_and_wrong_dimension() {
        let a = CellMesh::from_simplex(&box_tri([[0.0, 1.0], [0.0, 1.0]], 5, 5));
        let e = encode(&a);
        let mut bad = e.clone();
        bad[20] ^= 0x40;
        assert!(decode::<2>(&bad).is_err());
        assert!(decode::<2>(&e[..e.len() - 7]).is_err());
        assert!(decode::<3>(&e).is_err());
    }

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *seed >> 11
    }

    #[test]
    fn structured_meshes_compress_well() {
        let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 100, 100);
        let a = CellMesh::from_simplex(&m);
        let enc = encode(&a);
        // plain binary arrays: 8 bytes per coordinate + 4 per index + 4 per label
        let plain = a.n_points() * 16 + a.conn.len() * 4 + a.n_cells() * 4 + a.facet_conn.len() * 4 + a.n_facets() * 4;
        println!("fmsh {} B, plain binary {} B", enc.len(), plain);
        assert!(enc.len() * 8 < plain);
    }

    #[test]
    fn random_data_never_expands_noticeably() {
        let mut s = 11u64;
        let mut m = CellMesh::<3>::new();
        for _ in 0..4000 {
            m.points.push(std::array::from_fn(|_| {
                f64::from_bits(0x3FF0_0000_0000_0000 | (lcg(&mut s) & 0xF_FFFF_FFFF_FFFF))
            }));
        }
        let plain = m.points.len() * 24;
        assert!(encode(&m).len() < plain + plain / 100 + 64, "{} vs {plain}", encode(&m).len());
    }

    #[test]
    fn results_roundtrip_bit_exact() {
        let m = CellMesh::from_simplex(&box_tri([[0.0, 1.0], [0.0, 1.0]], 12, 10));
        let mut s = 5u64;
        let rnd =
            |s: &mut u64, n: usize| (0..n).map(|_| f64::from_bits(lcg(s) << 11 | (lcg(s) & 0x7FF))).collect::<Vec<_>>();
        let mut special = vec![0.0, -0.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE, 1.0 / 3.0];
        special.resize(m.n_points(), 2.5);
        let steps = vec![
            Step {
                value: 0.25,
                fields: vec![
                    Field {
                        name: "u".into(),
                        location: Location::Point,
                        components: 2,
                        data: rnd(&mut s, 2 * m.n_points()),
                    },
                    Field::scalar_on_cells("area", (0..m.n_cells()).map(|c| m.measure(c)).collect()),
                    Field::scalar_on_points("special", special),
                ],
            },
            Step { value: 0.5, fields: vec![Field::scalar_on_points("const", vec![1.0; m.n_points()])] },
        ];
        let enc = encode_steps(&m, &steps);
        let d = decode_dataset::<2>(&enc).unwrap();
        same(&m, &d.mesh);
        assert_eq!(d.steps.len(), 2);
        for (a, b) in steps.iter().zip(&d.steps) {
            assert_eq!(a.value.to_bits(), b.value.to_bits());
            assert_eq!(a.fields.len(), b.fields.len());
            for (f, g) in a.fields.iter().zip(&b.fields) {
                assert_eq!((&f.name, f.location, f.components), (&g.name, g.location, g.components));
                assert!(f.data.iter().zip(&g.data).all(|(x, y)| x.to_bits() == y.to_bits()), "{}", f.name);
            }
        }
        // a constant field costs next to nothing
        let with = encode_steps(&m, &steps[1..]).len();
        assert!(with < encode(&m).len() + 80, "{with}");
        // the plain decoder ignores results
        same(&m, &decode::<2>(&enc).unwrap());
    }

    #[test]
    fn smooth_fields_compress_but_stay_exact() {
        let m = CellMesh::from_simplex(&box_tri([[0.0, 1.0], [0.0, 1.0]], 40, 40));
        let u: Vec<f64> = m.points.iter().map(|p| (3.0 * p[0]).sin() * p[1]).collect();
        let enc = encode_steps(&m, &[Step { value: 0.0, fields: vec![Field::scalar_on_points("u", u.clone())] }]);
        let d = decode_dataset::<2>(&enc).unwrap();
        assert!(u.iter().zip(&d.steps[0].fields[0].data).all(|(a, b)| a.to_bits() == b.to_bits()));
        assert!(enc.len() < encode(&m).len() + u.len() * 8);
    }
}
