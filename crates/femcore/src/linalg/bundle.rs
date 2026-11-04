//! Problem bundles: a directory that holds everything needed to reproduce a solve elsewhere.
//!
//! ```text
//! A.mtx        stiffness-like matrix (Matrix Market coordinate)
//! B.mtx        optional mass-like matrix
//! v0.mtx       optional start vector (Matrix Market array); lets two solvers start identically
//! meta.txt     `key = value` lines: scalar (real|complex), symmetry (hermitian|general), name,
//!              suggested sigma_re/sigma_im/nev, free-form notes
//! ```
//!
//! nyx exports its matrices with [`Bundle::write`] (convert its assembly type to [`Csr`] with
//! [`Csr::from_triplets`]); the SLEPc harness in `harness/` reads the same files.

use crate::linalg::csr::Csr;
use crate::linalg::mm::{self, MmError};
use crate::linalg::scalar::Scalar;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Bundle<T: Scalar> {
    pub a: Csr<T>,
    pub b: Option<Csr<T>>,
    pub v0: Option<Vec<T>>,
    pub meta: BTreeMap<String, String>,
}

impl<T: Scalar> Bundle<T> {
    pub fn new(a: Csr<T>, b: Option<Csr<T>>) -> Self {
        let mut meta = BTreeMap::new();
        meta.insert("scalar".into(), if T::IS_COMPLEX { "complex" } else { "real" }.into());
        meta.insert("n".into(), a.nrows().to_string());
        meta.insert(
            "symmetry".into(),
            if a.is_hermitian(0.0) && b.as_ref().is_none_or(|b| b.is_hermitian(0.0)) { "hermitian" } else { "general" }
                .into(),
        );
        Self { a, b, v0: None, meta }
    }

    pub fn set(&mut self, key: &str, value: impl ToString) {
        self.meta.insert(key.into(), value.to_string());
    }

    pub fn get_f64(&self, key: &str) -> Option<f64> {
        self.meta.get(key)?.parse().ok()
    }

    pub fn write(&self, dir: &Path) -> Result<(), MmError> {
        fs::create_dir_all(dir)?;
        mm::write_matrix(BufWriter::new(File::create(dir.join("A.mtx"))?), &self.a)?;
        if let Some(b) = &self.b {
            mm::write_matrix(BufWriter::new(File::create(dir.join("B.mtx"))?), b)?;
        }
        if let Some(v0) = &self.v0 {
            mm::write_vectors(BufWriter::new(File::create(dir.join("v0.mtx"))?), v0, v0.len(), 1)?;
        }
        let mut f = BufWriter::new(File::create(dir.join("meta.txt"))?);
        for (k, v) in &self.meta {
            writeln!(f, "{k} = {v}")?;
        }
        f.flush()?;
        Ok(())
    }

    pub fn read(dir: &Path) -> Result<Self, MmError> {
        let a = mm::read_matrix(BufReader::new(File::open(dir.join("A.mtx"))?))?;
        let b = match File::open(dir.join("B.mtx")) {
            Ok(f) => Some(mm::read_matrix(BufReader::new(f))?),
            Err(_) => None,
        };
        let v0 = match File::open(dir.join("v0.mtx")) {
            Ok(f) => Some(mm::read_vectors(BufReader::new(f))?.0),
            Err(_) => None,
        };
        Ok(Self { a, b, v0, meta: read_meta(dir)? })
    }
}

/// Reads `meta.txt` alone (to learn the scalar type before reading the matrices).
pub fn read_meta(dir: &Path) -> Result<BTreeMap<String, String>, MmError> {
    let text = fs::read_to_string(dir.join("meta.txt"))?;
    Ok(text
        .lines()
        .filter_map(|l| l.split_once('='))
        .map(|(k, v)| (k.trim().to_string(), v.trim().to_string()))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr::{laplacian_2d, mass_1d};
    use crate::linalg::rng::random_vector;
    use crate::linalg::scalar::C64;

    #[test]
    fn bundle_round_trip() {
        let dir = std::env::temp_dir().join(format!("femsolve-bundle-test-{}", std::process::id()));
        let a = Csr::axpby(C64::ONE, &laplacian_2d(6, 5).to_complex(), C64::new(0.0, 0.2), &Csr::identity(30));
        let mut bun = Bundle::new(a.clone(), Some(mass_1d(30, 1.0).to_complex()));
        bun.v0 = Some(random_vector::<C64>(30, 7));
        bun.set("sigma_re", 1.5);
        bun.write(&dir).unwrap();
        let back = Bundle::<C64>::read(&dir).unwrap();
        assert_eq!(back.a, a);
        assert_eq!(back.b, bun.b);
        assert_eq!(back.v0, bun.v0);
        assert_eq!(back.get_f64("sigma_re"), Some(1.5));
        assert_eq!(back.meta["scalar"], "complex");
        assert_eq!(back.meta["symmetry"], "general");
        let _ = fs::remove_dir_all(&dir);
    }
}
