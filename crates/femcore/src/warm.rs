// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Warm starts of the eigensolvers (productivity plan P5b, ADR 96): the wanted modes of one solve as the start block of the next solve of a
//! sweep that keeps its mesh. A shift-and-invert Krylov-Schur (or block Arnoldi) run that begins from the previous point's `nev` eigenvectors
//! needs far fewer iterations when the problem moved a little; the factorisation, which is a new matrix at every point, is not touched.
//!
//! **Opt-in, and off is unchanged.** Like `progress`, a `Warm` is installed on the solving thread (`install`); the two eigensolver entry points
//! ask `begin` for a `Ticket` and get `None` when nothing is installed, so a solve without a context takes exactly the path it always took
//! (same start block, same bits). A solve with a context records its modes (`Ticket::finish`) for the next one whether or not it was started
//! warm itself, which is how the first point of a sweep feeds the second.
//!
//! **Matching.** The eigensolves of a solve are numbered in the order they are made (a chain of physics, a retry that asks for more modes); the
//! i-th takes the i-th entry of the input if it is the same kind of solve, on the same number of unknowns, for the same number of modes. A
//! mismatch (another mesh, another count) is a cold start, never an error: a start block can only make a solve faster or slower, not wrong.
//!
//! **What is not guaranteed.** A start block made only of the old modes can miss a mode that enters the wanted window and is orthogonal to
//! them; rounding noise usually lets it grow (slowly), but it is not guaranteed. A warm result is therefore never `checked` (the engine marks it
//! `exploratory`, ADR 91 and 96) until it is solved again cold.
//!
//! The entries serialise to bytes (`to_bytes`, `from_bytes`), so a worker process or the queue worker can be given the vectors and hand back the
//! ones it found.

use crate::linalg::scalar::{Scalar, C64};
use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

/// The modes one eigensolve found: its kind, the number of unknowns, the number of modes asked for, and the vectors (columns).
#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub kind: String,
    pub n: usize,
    pub nev: usize,
    pub cols: Vec<Vec<C64>>,
}

/// What one eigensolve did with a warm context, for the run's record.
#[derive(Clone, Debug, PartialEq)]
pub struct Call {
    pub kind: String,
    pub n: usize,
    pub nev: usize,
    /// started from the previous solve's modes
    pub warm: bool,
    /// restarts (Krylov-Schur) or cycles (block Arnoldi) the run took
    pub iterations: usize,
}

/// The warm context of one solve: the entries of the previous solve (input), the ones this solve leaves (output) and what each call did.
pub struct Warm {
    input: Vec<Option<Entry>>,
    next: AtomicUsize,
    output: Mutex<Vec<Option<Entry>>>,
    calls: Mutex<Vec<Call>>,
}

impl Warm {
    /// A context whose solve starts from `input` (empty: a cold first point that only records).
    pub fn new(input: Vec<Option<Entry>>) -> Arc<Warm> {
        Arc::new(Warm { input, next: AtomicUsize::new(0), output: Mutex::new(vec![]), calls: Mutex::new(vec![]) })
    }

    /// The entries this solve left, in the order of its eigensolves (a solve that failed or was not made leaves a gap).
    pub fn output(&self) -> Vec<Option<Entry>> {
        self.output.lock().expect("warm output").clone()
    }

    /// What every eigensolve did, in order.
    pub fn calls(&self) -> Vec<Call> {
        self.calls.lock().expect("warm calls").clone()
    }

    /// True when at least one eigensolve of this solve started from the previous solve's modes.
    pub fn was_used(&self) -> bool {
        self.calls.lock().expect("warm calls").iter().any(|c| c.warm)
    }
}

thread_local! {
    static CURRENT: RefCell<Option<Arc<Warm>>> = const { RefCell::new(None) };
}

/// Makes `w` the warm context of the current thread until the guard is dropped.
pub fn install(w: &Arc<Warm>) -> InstallGuard {
    let prev = CURRENT.with(|c| c.replace(Some(w.clone())));
    InstallGuard { prev }
}

pub struct InstallGuard {
    prev: Option<Arc<Warm>>,
}

impl Drop for InstallGuard {
    fn drop(&mut self) {
        let prev = self.prev.take();
        CURRENT.with(|c| *c.borrow_mut() = prev);
    }
}

/// The context installed on this thread, if any.
pub fn current() -> Option<Arc<Warm>> {
    CURRENT.with(|c| c.borrow().clone())
}

/// One eigensolve's share of the context: what to start from and where to leave its modes.
pub struct Ticket {
    warm: Arc<Warm>,
    index: usize,
    kind: &'static str,
    n: usize,
    nev: usize,
    start: Option<Vec<Vec<C64>>>,
}

/// Takes the next place of the installed context for an eigensolve of `kind` on `n` unknowns asking for `nev` modes; `None` without a context.
pub fn begin(kind: &'static str, n: usize, nev: usize) -> Option<Ticket> {
    let warm = current()?;
    let index = warm.next.fetch_add(1, Ordering::SeqCst);
    let start = warm
        .input
        .get(index)
        .and_then(Option::as_ref)
        .filter(|e| e.kind == kind && e.n == n && e.nev == nev && !e.cols.is_empty() && e.cols.iter().all(|c| c.len() == n))
        .map(|e| e.cols.clone());
    Some(Ticket { warm, index, kind, n, nev, start })
}

impl Ticket {
    /// The columns of the previous solve's modes, one after the other (`k * n` numbers), or `None` for a cold start.
    pub fn start_columns<T: Scalar>(&self) -> Option<Vec<T>> {
        let cols = self.start.as_ref()?;
        Some(cols.iter().flat_map(|c| c.iter().map(|z| T::from_parts(z.re, z.im))).collect())
    }

    /// How many columns the start has (0 for a cold start).
    pub fn start_len(&self) -> usize {
        self.start.as_ref().map_or(0, Vec::len)
    }

    /// Records the modes this eigensolve found (the first `nev` of them) and how it went.
    pub fn finish<'a, T: Scalar + 'a>(self, vectors: impl Iterator<Item = &'a [T]>, iterations: usize) {
        let cols: Vec<Vec<C64>> = vectors.take(self.nev).map(|v| v.iter().map(|z| z.to_c64()).collect()).collect();
        let mut out = self.warm.output.lock().expect("warm output");
        if out.len() <= self.index {
            out.resize(self.index + 1, None);
        }
        out[self.index] = Some(Entry { kind: self.kind.to_string(), n: self.n, nev: self.nev, cols });
        self.warm.calls.lock().expect("warm calls").push(Call {
            kind: self.kind.to_string(),
            n: self.n,
            nev: self.nev,
            warm: self.start.is_some(),
            iterations,
        });
    }
}

const MAGIC: &[u8; 8] = b"NYXWARM1";

/// The entries as bytes (little endian), for a worker process: a count, then for each a presence byte and the kind, sizes and numbers.
pub fn to_bytes(entries: &[Option<Entry>]) -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(MAGIC);
    b.extend_from_slice(&(entries.len() as u32).to_le_bytes());
    for e in entries {
        match e {
            None => b.push(0),
            Some(e) => {
                b.push(1);
                b.extend_from_slice(&(e.kind.len() as u16).to_le_bytes());
                b.extend_from_slice(e.kind.as_bytes());
                b.extend_from_slice(&(e.n as u64).to_le_bytes());
                b.extend_from_slice(&(e.nev as u64).to_le_bytes());
                b.extend_from_slice(&(e.cols.len() as u32).to_le_bytes());
                for c in &e.cols {
                    for z in c {
                        b.extend_from_slice(&z.re.to_le_bytes());
                        b.extend_from_slice(&z.im.to_le_bytes());
                    }
                }
            }
        }
    }
    b
}

/// Reads what `to_bytes` wrote; an error says what is wrong (a file of another version, a truncated one).
pub fn from_bytes(bytes: &[u8]) -> Result<Vec<Option<Entry>>, String> {
    struct R<'a>(&'a [u8]);
    impl<'a> R<'a> {
        fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
            if self.0.len() < n {
                return Err("the warm-start data is truncated".into());
            }
            let (a, b) = self.0.split_at(n);
            self.0 = b;
            Ok(a)
        }
        fn u64(&mut self) -> Result<u64, String> {
            Ok(u64::from_le_bytes(self.take(8)?.try_into().map_err(|_| "bad warm-start data")?))
        }
        fn f64(&mut self) -> Result<f64, String> {
            Ok(f64::from_le_bytes(self.take(8)?.try_into().map_err(|_| "bad warm-start data")?))
        }
    }
    let mut r = R(bytes);
    if r.take(8)? != MAGIC {
        return Err("the warm-start data is not of this version".into());
    }
    let count = u32::from_le_bytes(r.take(4)?.try_into().map_err(|_| "bad warm-start data")?) as usize;
    let mut out = Vec::with_capacity(count.min(1024));
    for _ in 0..count {
        if r.take(1)?[0] == 0 {
            out.push(None);
            continue;
        }
        let klen = u16::from_le_bytes(r.take(2)?.try_into().map_err(|_| "bad warm-start data")?) as usize;
        let kind = String::from_utf8(r.take(klen)?.to_vec()).map_err(|_| "bad warm-start data")?;
        let (n, nev) = (r.u64()? as usize, r.u64()? as usize);
        let ncols = u32::from_le_bytes(r.take(4)?.try_into().map_err(|_| "bad warm-start data")?) as usize;
        if n.checked_mul(ncols).is_none_or(|t| t.checked_mul(16).is_none_or(|b| b > r.0.len())) {
            return Err("the warm-start data is truncated".into());
        }
        let mut cols = Vec::with_capacity(ncols);
        for _ in 0..ncols {
            let mut c = Vec::with_capacity(n);
            for _ in 0..n {
                let (re, im) = (r.f64()?, r.f64()?);
                c.push(C64::new(re, im));
            }
            cols.push(c);
        }
        out.push(Some(Entry { kind, n, nev, cols }));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(n: usize) -> Entry {
        Entry { kind: "k".into(), n, nev: 2, cols: (0..2).map(|j| (0..n).map(|i| C64::new(i as f64 + j as f64, -(i as f64))).collect()).collect() }
    }

    #[test]
    fn entries_round_trip_through_bytes_and_a_damaged_file_is_refused() {
        let v = vec![Some(entry(5)), None, Some(entry(3))];
        let b = to_bytes(&v);
        assert_eq!(from_bytes(&b).unwrap(), v);
        assert!(from_bytes(&b[..b.len() - 3]).is_err());
        assert!(from_bytes(b"nonsense").is_err());
        let mut other = b.clone();
        other[3] = b'X';
        assert!(from_bytes(&other).is_err());
    }

    #[test]
    fn a_ticket_starts_warm_only_from_the_same_kind_size_and_count() {
        let input = vec![Some(entry(4)), Some(entry(4)), Some(entry(4))];
        let w = Warm::new(input);
        let _g = install(&w);
        assert_eq!(begin("k", 4, 2).unwrap().start_len(), 2);
        assert_eq!(begin("other", 4, 2).unwrap().start_len(), 0, "another kind of solve");
        assert_eq!(begin("k", 5, 2).unwrap().start_len(), 0, "another number of unknowns");
    }

    #[test]
    fn without_a_context_there_is_no_ticket() {
        assert!(begin("k", 4, 2).is_none());
    }

    #[test]
    fn finish_leaves_the_modes_and_the_record_in_order() {
        let w = Warm::new(vec![]);
        let _g = install(&w);
        let a = begin("k", 3, 2).unwrap();
        let b = begin("k", 3, 1).unwrap();
        let v: Vec<f64> = vec![1.0, 2.0, 3.0];
        b.finish([v.as_slice()].into_iter(), 4);
        a.finish([v.as_slice(), v.as_slice(), v.as_slice()].into_iter(), 7);
        let out = w.output();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].as_ref().unwrap().cols.len(), 2, "the first nev modes");
        assert_eq!(out[1].as_ref().unwrap().cols.len(), 1);
        assert!(!w.was_used());
        assert_eq!(w.calls().len(), 2);
    }
}
