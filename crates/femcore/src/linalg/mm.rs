//! Matrix Market (`.mtx`) reading and writing: the export format of port item V1.
//!
//! Matrices are written in `coordinate` format (real or complex; `symmetric`/`hermitian`
//! storage is chosen automatically when the matrix has that structure, which halves the file),
//! vectors in dense `array` format. Numbers are printed with 17 significant digits so a
//! write/read round trip is bit exact. A reader for the same format is in the SLEPc harness
//! (`harness/c/common.h`), so both sides see identical input.
//!
//! Complex *symmetric* matrices (PML problems) are written as `symmetric`, which in Matrix
//! Market means `a_ij = a_ji` *without* conjugation — the reader in `common.h` honours that.

use crate::linalg::csr::Csr;
use crate::linalg::scalar::{C64, Scalar};
use std::fmt::Write as _;
use std::io::{self, BufRead, Write};

#[derive(Debug)]
pub enum MmError {
    Io(io::Error),
    Parse(String),
    /// The file holds complex numbers but a real matrix was requested.
    ComplexIntoReal,
}

impl std::fmt::Display for MmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MmError::Io(e) => write!(f, "I/O error: {e}"),
            MmError::Parse(m) => write!(f, "Matrix Market parse error: {m}"),
            MmError::ComplexIntoReal => write!(f, "file is complex but a real matrix was requested"),
        }
    }
}

impl std::error::Error for MmError {}

impl From<io::Error> for MmError {
    fn from(e: io::Error) -> Self {
        MmError::Io(e)
    }
}

fn parse_err<S: Into<String>>(s: S) -> MmError {
    MmError::Parse(s.into())
}

fn fmt_scalar<T: Scalar>(v: T) -> String {
    if T::IS_COMPLEX { format!("{:.17e} {:.17e}", v.re(), v.im()) } else { format!("{:.17e}", v.re()) }
}

/// Writes `a` in coordinate format, exploiting symmetry/Hermitian structure when exact.
pub fn write_matrix<T: Scalar, W: Write>(mut w: W, a: &Csr<T>) -> Result<(), MmError> {
    let field = if T::IS_COMPLEX { "complex" } else { "real" };
    let sym = if a.nrows() == a.ncols() && a.is_symmetric(0.0) {
        "symmetric"
    } else if T::IS_COMPLEX && a.nrows() == a.ncols() && a.is_hermitian(0.0) {
        "hermitian"
    } else {
        "general"
    };
    let lower_only = sym != "general";
    let mut entries = Vec::new();
    for i in 0..a.nrows() {
        let (cols, vals) = a.row(i);
        for (&c, &v) in cols.iter().zip(vals) {
            if !lower_only || (c as usize) <= i {
                entries.push((i, c as usize, v));
            }
        }
    }
    writeln!(w, "%%MatrixMarket matrix coordinate {field} {sym}")?;
    writeln!(w, "{} {} {}", a.nrows(), a.ncols(), entries.len())?;
    let mut buf = String::new();
    for (i, c, v) in entries {
        buf.clear();
        let _ = write!(buf, "{} {} {}", i + 1, c + 1, fmt_scalar(v));
        writeln!(w, "{buf}")?;
    }
    Ok(())
}

/// Writes column vectors (`n × ncols`, column-major) in dense array format.
pub fn write_vectors<T: Scalar, W: Write>(mut w: W, data: &[T], n: usize, ncols: usize) -> Result<(), MmError> {
    assert!(data.len() >= n * ncols);
    let field = if T::IS_COMPLEX { "complex" } else { "real" };
    writeln!(w, "%%MatrixMarket matrix array {field} general")?;
    writeln!(w, "{n} {ncols}")?;
    for v in &data[..n * ncols] {
        writeln!(w, "{}", fmt_scalar(*v))?;
    }
    Ok(())
}

struct Header {
    format: String,
    field: String,
    symmetry: String,
}

fn read_header<R: BufRead>(r: &mut R) -> Result<(Header, String), MmError> {
    let mut line = String::new();
    r.read_line(&mut line)?;
    let toks: Vec<String> = line.split_whitespace().map(|s| s.to_ascii_lowercase()).collect();
    if toks.len() < 5 || toks[0] != "%%matrixmarket" || toks[1] != "matrix" {
        return Err(parse_err("missing %%MatrixMarket matrix header"));
    }
    let header = Header { format: toks[2].clone(), field: toks[3].clone(), symmetry: toks[4].clone() };
    // Skip comments, return the size line.
    loop {
        line.clear();
        if r.read_line(&mut line)? == 0 {
            return Err(parse_err("unexpected end of file before the size line"));
        }
        if !line.trim_start().starts_with('%') && !line.trim().is_empty() {
            return Ok((header, line.clone()));
        }
    }
}

fn parse_value<T: Scalar>(field: &str, toks: &mut std::str::SplitWhitespace<'_>) -> Result<T, MmError> {
    let mut next = || -> Result<f64, MmError> {
        toks.next().ok_or_else(|| parse_err("missing value"))?.parse::<f64>().map_err(|e| parse_err(e.to_string()))
    };
    match field {
        "real" | "integer" => Ok(T::from_f64(next()?)),
        "complex" => {
            let (re, im) = (next()?, next()?);
            if !T::IS_COMPLEX && im != 0.0 {
                return Err(MmError::ComplexIntoReal);
            }
            Ok(T::from_parts(re, im))
        }
        "pattern" => Ok(T::ONE),
        other => Err(parse_err(format!("unsupported field '{other}'"))),
    }
}

/// Reads a coordinate-format matrix, expanding symmetric/skew/Hermitian storage.
pub fn read_matrix<T: Scalar, R: BufRead>(mut r: R) -> Result<Csr<T>, MmError> {
    let (h, size) = read_header(&mut r)?;
    if h.format != "coordinate" {
        return Err(parse_err("expected coordinate format"));
    }
    let dims: Vec<usize> =
        size.split_whitespace().map(|t| t.parse().map_err(|_| parse_err("bad size line"))).collect::<Result<_, _>>()?;
    if dims.len() != 3 {
        return Err(parse_err("size line must hold rows, cols, nnz"));
    }
    let (nrows, ncols, nnz) = (dims[0], dims[1], dims[2]);
    let mut trip: Vec<(usize, usize, T)> = Vec::with_capacity(if h.symmetry == "general" { nnz } else { 2 * nnz });
    let mut line = String::new();
    for _ in 0..nnz {
        line.clear();
        if r.read_line(&mut line)? == 0 {
            return Err(parse_err("fewer entries than announced"));
        }
        let mut toks = line.split_whitespace();
        let i: usize =
            toks.next().ok_or_else(|| parse_err("missing row"))?.parse().map_err(|_| parse_err("bad row"))?;
        let j: usize =
            toks.next().ok_or_else(|| parse_err("missing col"))?.parse().map_err(|_| parse_err("bad col"))?;
        if i == 0 || j == 0 || i > nrows || j > ncols {
            return Err(parse_err("index out of range"));
        }
        let v: T = parse_value(&h.field, &mut toks)?;
        trip.push((i - 1, j - 1, v));
        if i != j {
            match h.symmetry.as_str() {
                "symmetric" => trip.push((j - 1, i - 1, v)),
                "skew-symmetric" => trip.push((j - 1, i - 1, -v)),
                "hermitian" => trip.push((j - 1, i - 1, v.conj())),
                _ => {}
            }
        }
    }
    Csr::from_triplets(nrows, ncols, &trip).map_err(|e| parse_err(e.to_string()))
}

/// Reads dense `array` data; returns `(column-major values, rows, cols)`.
pub fn read_vectors<T: Scalar, R: BufRead>(mut r: R) -> Result<(Vec<T>, usize, usize), MmError> {
    let (h, size) = read_header(&mut r)?;
    if h.format != "array" {
        return Err(parse_err("expected array format"));
    }
    let dims: Vec<usize> =
        size.split_whitespace().map(|t| t.parse().map_err(|_| parse_err("bad size line"))).collect::<Result<_, _>>()?;
    if dims.len() != 2 {
        return Err(parse_err("size line must hold rows, cols"));
    }
    let (n, c) = (dims[0], dims[1]);
    let mut out = Vec::with_capacity(n * c);
    let mut line = String::new();
    while out.len() < n * c {
        line.clear();
        if r.read_line(&mut line)? == 0 {
            return Err(parse_err("fewer values than announced"));
        }
        if line.trim().is_empty() {
            continue;
        }
        out.push(parse_value(&h.field, &mut line.split_whitespace())?);
    }
    Ok((out, n, c))
}

/// Convenience: complex eigenvalue lists as `re im` text lines (used by the bundle's
/// `expected.txt` and by the harness).
pub fn format_c64_list(v: &[C64]) -> String {
    v.iter().map(|z| format!("{:.17e} {:.17e}\n", z.re, z.im)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr::{laplacian_2d, mass_1d};
    use crate::linalg::rng::random_vector;

    fn round_trip<T: Scalar>(a: &Csr<T>) -> Csr<T> {
        let mut buf = Vec::new();
        write_matrix(&mut buf, a).unwrap();
        read_matrix::<T, _>(buf.as_slice()).unwrap()
    }

    #[test]
    fn real_symmetric_round_trip_is_exact_and_halves_storage() {
        let a = laplacian_2d(9, 7);
        let mut buf = Vec::new();
        write_matrix(&mut buf, &a).unwrap();
        let text = String::from_utf8(buf.clone()).unwrap();
        assert!(text.starts_with("%%MatrixMarket matrix coordinate real symmetric"));
        let nnz_written: usize = text.lines().nth(1).unwrap().split_whitespace().nth(2).unwrap().parse().unwrap();
        assert_eq!(nnz_written, (a.nnz() + a.nrows()) / 2);
        assert_eq!(read_matrix::<f64, _>(buf.as_slice()).unwrap(), a);
        assert_eq!(round_trip(&mass_1d(30, 0.37)), mass_1d(30, 0.37));
    }

    #[test]
    fn complex_hermitian_symmetric_and_general() {
        let n = 20;
        let h = Csr::from_triplets(
            n,
            n,
            &(0..n)
                .flat_map(|i| {
                    let mut v = vec![(i, i, C64::new(2.0, 0.0))];
                    if i + 1 < n {
                        v.push((i, i + 1, C64::new(-1.0, 0.5)));
                        v.push((i + 1, i, C64::new(-1.0, -0.5)));
                    }
                    v
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let mut buf = Vec::new();
        write_matrix(&mut buf, &h).unwrap();
        assert!(String::from_utf8(buf).unwrap().starts_with("%%MatrixMarket matrix coordinate complex hermitian"));
        assert_eq!(round_trip(&h), h);
        // Complex symmetric (not Hermitian): PML-like.
        let s = h.map(|v| v * C64::new(0.0, 1.0));
        let t = Csr::axpby(C64::ONE, &laplacian_2d(4, 5).to_complex(), C64::new(0.0, 0.3), &Csr::identity(20));
        let mut buf = Vec::new();
        write_matrix(&mut buf, &t).unwrap();
        assert!(String::from_utf8(buf).unwrap().starts_with("%%MatrixMarket matrix coordinate complex symmetric"));
        assert_eq!(round_trip(&t), t);
        assert_eq!(round_trip(&s).nnz(), s.nnz());
        // General (non-symmetric).
        let g = Csr::from_triplets(3, 3, &[(0, 1, C64::new(1.0, 2.0)), (2, 0, C64::new(-3.0, 0.0))]).unwrap();
        assert_eq!(round_trip(&g), g);
    }

    #[test]
    fn vectors_round_trip_and_type_checks() {
        let v = random_vector::<C64>(17 * 3, 4);
        let mut buf = Vec::new();
        write_vectors(&mut buf, &v, 17, 3).unwrap();
        let (back, n, c) = read_vectors::<C64, _>(buf.as_slice()).unwrap();
        assert_eq!((n, c), (17, 3));
        assert_eq!(back, v);
        // A complex file cannot be read into a real matrix, but a real one widens to complex.
        let a = Csr::from_triplets(2, 2, &[(0, 0, C64::new(1.0, 1.0))]).unwrap();
        let mut buf = Vec::new();
        write_matrix(&mut buf, &a).unwrap();
        assert!(matches!(read_matrix::<f64, _>(buf.as_slice()), Err(MmError::ComplexIntoReal)));
        let r = laplacian_2d(3, 3);
        let mut buf = Vec::new();
        write_matrix(&mut buf, &r).unwrap();
        assert_eq!(read_matrix::<C64, _>(buf.as_slice()).unwrap(), r.to_complex());
        assert!(read_matrix::<f64, _>("garbage".as_bytes()).is_err());
    }
}
