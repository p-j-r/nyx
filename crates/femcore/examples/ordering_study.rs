//! S2 experiment: does a nested-dissection ordering beat faer's default AMD for sparse Cholesky?
//! Graph-based nested dissection (BFS level-set separators, recursion down to 100 nodes) against AMD
//! on 2-D and 3-D grid Laplacians: factor non-zeros and numeric factorisation time.
//! `cargo run --release --example ordering_study`
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::perm::PermRef;
use faer::sparse::SparseColMatRef;
use faer::sparse::SymbolicSparseColMatRef;
use faer::sparse::linalg::cholesky::{CholeskySymbolicParams, SymmetricOrdering, factorize_symbolic_cholesky};
use faer::{Par, Side};
use femcore::linalg::Csr;
use femcore::linalg::csr::{laplacian_2d, laplacian_3d};
use std::collections::VecDeque;
use std::time::Instant;

fn bfs_levels(adj: &[Vec<u32>], nodes: &[usize], mark: &mut [u32], id: u32, start: usize) -> Vec<Vec<usize>> {
    let mut seen = std::collections::HashSet::new();
    let mut levels: Vec<Vec<usize>> = vec![vec![start]];
    seen.insert(start);
    let mut q = VecDeque::new();
    q.push_back((start, 0usize));
    while let Some((u, l)) = q.pop_front() {
        for &v in &adj[u] {
            let v = v as usize;
            if mark[v] == id && seen.insert(v) {
                if levels.len() <= l + 1 {
                    levels.push(Vec::new());
                }
                levels[l + 1].push(v);
                q.push_back((v, l + 1));
            }
        }
    }
    let _ = nodes;
    levels
}

/// Appends an elimination order of `nodes` (subdomains first, separator last).
fn nested_dissection(adj: &[Vec<u32>], nodes: Vec<usize>, mark: &mut [u32], next_id: &mut u32, out: &mut Vec<usize>) {
    if nodes.len() <= 100 {
        out.extend(nodes);
        return;
    }
    *next_id += 1;
    let id = *next_id;
    for &v in &nodes {
        mark[v] = id;
    }
    // pseudo-peripheral start: BFS twice
    let l0 = bfs_levels(adj, &nodes, mark, id, nodes[0]);
    let far = *l0.last().unwrap().first().unwrap();
    let levels = bfs_levels(adj, &nodes, mark, id, far);
    let reached: usize = levels.iter().map(|l| l.len()).sum();
    if reached < nodes.len() {
        // disconnected: split off the reached component
        let in_comp: std::collections::HashSet<usize> = levels.iter().flatten().copied().collect();
        let (a, b): (Vec<usize>, Vec<usize>) = nodes.iter().partition(|v| in_comp.contains(v));
        nested_dissection(adj, a, mark, next_id, out);
        nested_dissection(adj, b, mark, next_id, out);
        return;
    }
    // separator = the level where the cumulative count first reaches half
    let mut cum = 0;
    let mut mid = levels.len() / 2;
    for (i, l) in levels.iter().enumerate() {
        cum += l.len();
        if cum * 2 >= nodes.len() {
            mid = i;
            break;
        }
    }
    let mid = mid.clamp(1, levels.len().saturating_sub(2).max(1));
    let a: Vec<usize> = levels[..mid].iter().flatten().copied().collect();
    let b: Vec<usize> = levels[mid + 1..].iter().flatten().copied().collect();
    let sep = levels[mid].clone();
    if a.is_empty() || b.is_empty() {
        out.extend(nodes);
        return;
    }
    nested_dissection(adj, a, mark, next_id, out);
    nested_dissection(adj, b, mark, next_id, out);
    out.extend(sep);
}

fn study(name: &str, a: &Csr<f64>) {
    let n = a.nrows();
    // lower-triangle CSC pattern and values: column j holds the rows i ≥ j (= the upper triangle of row j by symmetry)
    let mut col_ptr = vec![0u32];
    let (mut row_idx, mut vals) = (Vec::new(), Vec::new());
    for j in 0..n {
        let (cs, vs) = a.row(j);
        for (&c, &v) in cs.iter().zip(vs) {
            if c as usize >= j {
                row_idx.push(c);
                vals.push(v);
            }
        }
        col_ptr.push(row_idx.len() as u32);
    }
    let adj: Vec<Vec<u32>> =
        (0..n).map(|i| a.row(i).0.iter().copied().filter(|&c| c as usize != i).collect()).collect();
    let mut order = Vec::with_capacity(n);
    let t = Instant::now();
    nested_dissection(&adj, (0..n).collect(), &mut vec![0u32; n], &mut 0, &mut order);
    let t_nd = t.elapsed().as_secs_f64();
    assert_eq!(order.len(), n);
    // PermRef: forward[new] = old
    let fwd: Vec<u32> = order.iter().map(|&o| o as u32).collect();
    let mut inv = vec![0u32; n];
    for (new, &old) in order.iter().enumerate() {
        inv[old] = new as u32;
    }
    let pattern = SymbolicSparseColMatRef::new_checked(n, n, &col_ptr, None, &row_idx);
    let run = |label: &str, ord: SymmetricOrdering<'_, u32>| {
        let t = Instant::now();
        let sym = factorize_symbolic_cholesky(pattern, Side::Lower, ord, CholeskySymbolicParams::default()).unwrap();
        let t_sym = t.elapsed().as_secs_f64();
        let mut lvals = vec![0.0f64; sym.len_val()];
        let mat = SparseColMatRef::new(pattern, &vals);
        let mut mem = MemBuffer::new(sym.factorize_numeric_llt_scratch::<f64>(Par::Seq, Default::default()));
        let t = Instant::now();
        sym.factorize_numeric_llt::<f64>(
            &mut lvals,
            mat,
            Side::Lower,
            Default::default(),
            Par::Seq,
            MemStack::new(&mut mem),
            Default::default(),
        )
        .unwrap();
        let t_num = t.elapsed().as_secs_f64();
        println!(
            "  {label:<10} factor nnz {:>10}   symbolic {t_sym:.3} s   numeric {t_num:.3} s (sequential)",
            sym.len_val()
        );
    };
    println!("{name} (n = {n}), nested dissection built in {t_nd:.3} s");
    run("AMD", SymmetricOrdering::Amd);
    let perm = PermRef::new_checked(&fwd, &inv, n);
    run("ND", SymmetricOrdering::Custom(perm));
}

fn main() {
    study("2-D Laplacian 200x200", &laplacian_2d(200, 200));
    study("3-D Laplacian 40^3", &laplacian_3d(40));
}
