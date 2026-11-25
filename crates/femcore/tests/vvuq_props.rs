// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Property-based (`proptest`) checks for 2D/3D `SimplexMesh` splits and refinements, `Csr` matrix algebra, and `formula::parse` AST round-trips.

use femcore::adapt::drive::{adapt2, AdaptOptions};
use femcore::adapt::drive3::adapt3;
use femcore::adapt::metric::Sym;
use femcore::adapt::ops3::AdaptOptions3;
use femcore::adapt::topo3::TopoOptions3;
use femcore::adapt::work::{simplex_from_adapt, simplex_to_adapt};
use femcore::fieldexpr::{Cmp, Fun, Fun2, Op, Tree};
use femcore::formula::{parse as parse_formula, Named, Names};
use femcore::linalg::csr::Csr;
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::fixtures::{brick_tet, rect_tri};
use femcore::meshing::refine3d::refine_to_size;
use num_complex::Complex64 as C64;
use proptest::prelude::*;
use std::collections::BTreeMap;

fn spacings_to_coords(spacings: &[f64]) -> Vec<f64> {
    let mut xs = Vec::with_capacity(spacings.len() + 1);
    let mut acc = 0.0;
    xs.push(acc);
    for &dx in spacings {
        acc += dx;
        xs.push(acc);
    }
    xs
}

fn arb_grid_axis(max_spans: usize) -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec((2u8..=8).prop_map(|k| f64::from(k) * 0.25), 1..=max_spans)
        .prop_map(|d| spacings_to_coords(&d))
}

fn push_2d_boundary_facets(cm: &mut CellMesh<2>, nx: usize, ny: usize) {
    let idx = |ix: usize, iy: usize| -> u32 { (iy * (nx + 1) + ix) as u32 };
    for ix in 0..nx {
        cm.push_facet(&[idx(ix, 0), idx(ix + 1, 0)], 1);
        cm.push_facet(&[idx(ix + 1, ny), idx(ix, ny)], 2);
    }
    for iy in 0..ny {
        cm.push_facet(&[idx(0, iy + 1), idx(0, iy)], 3);
        cm.push_facet(&[idx(nx, iy), idx(nx, iy + 1)], 4);
    }
}

fn push_3d_prism_boundary_facets(cm: &mut CellMesh<3>, nx: usize, ny: usize, nz: usize) {
    let idx = |ix: usize, iy: usize, iz: usize| -> u32 { ((iz * (ny + 1) + iy) * (nx + 1) + ix) as u32 };
    for iy in 0..ny {
        for ix in 0..nx {
            let (a, b, c, d) = (idx(ix, iy, 0), idx(ix + 1, iy, 0), idx(ix + 1, iy + 1, 0), idx(ix, iy + 1, 0));
            cm.push_facet(&[a, c, b], 1);
            cm.push_facet(&[a, d, c], 1);
            let (a, b, c, d) = (idx(ix, iy, nz), idx(ix + 1, iy, nz), idx(ix + 1, iy + 1, nz), idx(ix, iy + 1, nz));
            cm.push_facet(&[a, b, c], 2);
            cm.push_facet(&[a, c, d], 2);
        }
    }
    for iz in 0..nz {
        for ix in 0..nx {
            cm.push_facet(&[idx(ix, 0, iz), idx(ix + 1, 0, iz), idx(ix + 1, 0, iz + 1), idx(ix, 0, iz + 1)], 3);
            cm.push_facet(&[idx(ix, ny, iz), idx(ix, ny, iz + 1), idx(ix + 1, ny, iz + 1), idx(ix + 1, ny, iz)], 4);
        }
    }
    for iz in 0..nz {
        for iy in 0..ny {
            cm.push_facet(&[idx(0, iy, iz), idx(0, iy, iz + 1), idx(0, iy + 1, iz + 1), idx(0, iy + 1, iz)], 5);
            cm.push_facet(&[idx(nx, iy, iz), idx(nx, iy + 1, iz), idx(nx, iy + 1, iz + 1), idx(nx, iy, iz + 1)], 6);
        }
    }
}

fn push_3d_boundary_quad_facets(cm: &mut CellMesh<3>, nx: usize, ny: usize, nz: usize) {
    let idx = |ix: usize, iy: usize, iz: usize| -> u32 { ((iz * (ny + 1) + iy) * (nx + 1) + ix) as u32 };
    for iy in 0..ny {
        for ix in 0..nx {
            cm.push_facet(&[idx(ix, iy, 0), idx(ix, iy + 1, 0), idx(ix + 1, iy + 1, 0), idx(ix + 1, iy, 0)], 1);
            cm.push_facet(&[idx(ix, iy, nz), idx(ix + 1, iy, nz), idx(ix + 1, iy + 1, nz), idx(ix, iy + 1, nz)], 2);
        }
    }
    for iz in 0..nz {
        for ix in 0..nx {
            cm.push_facet(&[idx(ix, 0, iz), idx(ix + 1, 0, iz), idx(ix + 1, 0, iz + 1), idx(ix, 0, iz + 1)], 3);
            cm.push_facet(&[idx(ix, ny, iz), idx(ix, ny, iz + 1), idx(ix + 1, ny, iz + 1), idx(ix + 1, ny, iz)], 4);
        }
    }
    for iz in 0..nz {
        for iy in 0..ny {
            cm.push_facet(&[idx(0, iy, iz), idx(0, iy, iz + 1), idx(0, iy + 1, iz + 1), idx(0, iy + 1, iz)], 5);
            cm.push_facet(&[idx(nx, iy, iz), idx(nx, iy + 1, iz), idx(nx, iy + 1, iz + 1), idx(nx, iy, iz + 1)], 6);
        }
    }
}

fn triplets_to_dense(nrows: usize, ncols: usize, tris: &[(usize, usize, f64)]) -> Vec<f64> {
    let mut d = vec![0.0; nrows * ncols];
    for &(r, c, v) in tris {
        d[r * ncols + c] += v;
    }
    d
}

fn arb_triplets(nrows: usize, ncols: usize, len: usize) -> impl Strategy<Value = Vec<(usize, usize, f64)>> {
    prop::collection::vec(
        (0..nrows, 0..ncols, (-16i32..=16).prop_map(|k| f64::from(k) * 0.25)),
        0..=len,
    )
}

fn arb_c64_triplets(nrows: usize, ncols: usize, len: usize) -> impl Strategy<Value = Vec<(usize, usize, C64)>> {
    prop::collection::vec(
        (
            0..nrows,
            0..ncols,
            (-8i32..=8, -8i32..=8).prop_map(|(re, im)| C64::new(f64::from(re) * 0.25, f64::from(im) * 0.25)),
        ),
        0..=len,
    )
}

struct PropFormulaNames;

impl Names for PropFormulaNames {
    fn name(&mut self, name: &str) -> Result<Option<Named>, String> {
        Ok(match name {
            "x" => Some(Named::Scalar(Tree::Leaf(0))),
            "y" => Some(Named::Scalar(Tree::Leaf(1))),
            "z" => Some(Named::Scalar(Tree::Leaf(2))),
            "u" => Some(Named::Scalar(Tree::Leaf(3))),
            _ => None,
        })
    }

    fn material(&mut self, _name: &str) -> Result<Tree, String> {
        Ok(Tree::num(1.5))
    }

    fn region(&mut self, _name: &str) -> Result<Tree, String> {
        Ok(Tree::Leaf(4))
    }
}

#[derive(Clone, Debug)]
enum Ast {
    Num(f64),
    Leaf(&'static str),
    Neg(Box<Ast>),
    Bin(Op, Box<Ast>, Box<Ast>),
    Fun1(Fun, Box<Ast>),
    Fun2(Fun2, Box<Ast>, Box<Ast>),
    Cmp(Cmp, Box<Ast>, Box<Ast>),
    Where(Box<Ast>, Box<Ast>, Box<Ast>),
}

fn arb_ast() -> impl Strategy<Value = Ast> {
    let leaf = prop_oneof![
        (1i32..=16).prop_map(|k| Ast::Num(f64::from(k) * 0.25)),
        Just(Ast::Leaf("x")),
        Just(Ast::Leaf("y")),
        Just(Ast::Leaf("z")),
        Just(Ast::Leaf("u")),
    ];
    leaf.prop_recursive(3, 24, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|a| Ast::Neg(Box::new(a))),
            (
                prop_oneof![Just(Op::Add), Just(Op::Sub), Just(Op::Mul), Just(Op::And), Just(Op::Or)],
                inner.clone(),
                inner.clone(),
            )
                .prop_map(|(op, a, b)| Ast::Bin(op, Box::new(a), Box::new(b))),
            (inner.clone(), inner.clone()).prop_map(|(a, b)| Ast::Bin(
                Op::Div,
                Box::new(a),
                Box::new(Ast::Bin(Op::Add, Box::new(Ast::Num(1.0)), Box::new(Ast::Fun1(Fun::Abs2, Box::new(b))))),
            )),
            (
                prop_oneof![Just(Fun::Sin), Just(Fun::Cos), Just(Fun::Abs), Just(Fun::Abs2), Just(Fun::Tanh), Just(Fun::Step), Just(Fun::Not)],
                inner.clone(),
            )
                .prop_map(|(f, a)| Ast::Fun1(f, Box::new(a))),
            (prop_oneof![Just(Fun2::Min), Just(Fun2::Max)], inner.clone(), inner.clone())
                .prop_map(|(f, a, b)| Ast::Fun2(f, Box::new(a), Box::new(b))),
            (
                prop_oneof![Just(Cmp::Lt), Just(Cmp::Le), Just(Cmp::Gt), Just(Cmp::Ge), Just(Cmp::Eq), Just(Cmp::Ne)],
                inner.clone(),
                inner.clone(),
            )
                .prop_map(|(c, a, b)| Ast::Cmp(c, Box::new(a), Box::new(b))),
            (inner.clone(), inner.clone(), inner).prop_map(|(c, a, b)| Ast::Where(Box::new(c), Box::new(a), Box::new(b))),
        ]
    })
}

fn ast_to_str(ast: &Ast) -> String {
    match ast {
        Ast::Num(x) => format!("{x}"),
        Ast::Leaf(name) => (*name).to_string(),
        Ast::Neg(a) => format!("(-{})", ast_to_str(a)),
        Ast::Bin(op, a, b) => {
            let sym = match op {
                Op::Add => "+",
                Op::Sub => "-",
                Op::Mul => "*",
                Op::Div => "/",
                Op::And => "and",
                Op::Or => "or",
            };
            format!("({} {sym} {})", ast_to_str(a), ast_to_str(b))
        }
        Ast::Fun1(Fun::Not, a) => format!("(not {})", ast_to_str(a)),
        Ast::Fun1(f, a) => format!("{}({})", f.name(), ast_to_str(a)),
        Ast::Fun2(f, a, b) => format!("{}({}, {})", f.name(), ast_to_str(a), ast_to_str(b)),
        Ast::Cmp(c, a, b) => {
            let sym = match c {
                Cmp::Lt => "<",
                Cmp::Le => "<=",
                Cmp::Gt => ">",
                Cmp::Ge => ">=",
                Cmp::Eq => "==",
                Cmp::Ne => "!=",
            };
            format!("({} {sym} {})", ast_to_str(a), ast_to_str(b))
        }
        Ast::Where(c, a, b) => format!("where({}, {}, {})", ast_to_str(c), ast_to_str(a), ast_to_str(b)),
    }
}

fn tree_to_formula(t: &Tree) -> String {
    match t {
        Tree::Num(z) => {
            if z.im == 0.0 {
                format!("{:.17}", z.re)
            } else if z.re == 0.0 {
                format!("({:.17}j)", z.im)
            } else {
                format!("({:.17} + {:.17}j)", z.re, z.im)
            }
        }
        Tree::Leaf(i) => ["x", "y", "z", "u"][*i].to_string(),
        Tree::Neg(a) => format!("(-{})", tree_to_formula(a)),
        Tree::Bin(op, a, b) => {
            let sym = match op {
                Op::Add => "+",
                Op::Sub => "-",
                Op::Mul => "*",
                Op::Div => "/",
                Op::And => "and",
                Op::Or => "or",
            };
            format!("({} {sym} {})", tree_to_formula(a), tree_to_formula(b))
        }
        Tree::Call(Fun::Not, a) => format!("(not {})", tree_to_formula(a)),
        Tree::Call(f, a) => format!("{}({})", f.name(), tree_to_formula(a)),
        Tree::Call2(Fun2::Pow, a, b) => format!("({} ^ {})", tree_to_formula(a), tree_to_formula(b)),
        Tree::Call2(f, a, b) => format!("{}({}, {})", f.name(), tree_to_formula(a), tree_to_formula(b)),
        Tree::Pow(a, n) => format!("({} ^ {n})", tree_to_formula(a)),
        Tree::Cmp(c, a, b) => {
            let sym = match c {
                Cmp::Lt => "<",
                Cmp::Le => "<=",
                Cmp::Gt => ">",
                Cmp::Ge => ">=",
                Cmp::Eq => "==",
                Cmp::Ne => "!=",
            };
            format!("({} {sym} {})", tree_to_formula(a), tree_to_formula(b))
        }
        Tree::Where(c, a, b) => format!("where({}, {}, {})", tree_to_formula(c), tree_to_formula(a), tree_to_formula(b)),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn mesh2d_cell_splits_and_adapt2_preserve_check_and_area(
        xs in arb_grid_axis(3),
        ys in arb_grid_axis(3),
        quad_flags in prop::collection::vec(any::<bool>(), 9..=9),
        label_base in 1i32..=5,
        h_scale in 0.4_f64..1.2,
    ) {
        let nx = xs.len() - 1;
        let ny = ys.len() - 1;
        let mut cm = CellMesh::<2>::new();
        for &y in &ys {
            for &x in &xs {
                cm.points.push([x, y]);
            }
        }
        let idx = |ix: usize, iy: usize| -> u32 { (iy * (nx + 1) + ix) as u32 };
        let mut cell_idx = 0usize;
        for iy in 0..ny {
            for ix in 0..nx {
                let (v00, v10, v11, v01) = (idx(ix, iy), idx(ix + 1, iy), idx(ix + 1, iy + 1), idx(ix, iy + 1));
                let lab = label_base + ((ix ^ iy) as i32 % 2);
                if quad_flags[cell_idx % quad_flags.len()] {
                    cm.push_cell(CellType::Quad, &[v00, v10, v11, v01], lab);
                } else {
                    cm.push_cell(CellType::Tri, &[v00, v10, v11], lab);
                    cm.push_cell(CellType::Tri, &[v00, v11, v01], lab);
                }
                cell_idx += 1;
            }
        }
        push_2d_boundary_facets(&mut cm, nx, ny);
        cm.check().expect("2D mixed CellMesh valid");
        let area_before = cm.total_measure();

        let split = cm.split_to_simplices().expect("2D split succeeds");
        split.check().expect("split SimplexMesh<2>::check passes");
        let split_nyx = simplex_from_adapt(&split);
        split_nyx.check().expect("split NyxMesh::check passes");
        prop_assert!((CellMesh::from_simplex(&split).total_measure() - area_before).abs() < 1e-11);

        let rt = simplex_to_adapt::<2>(&split_nyx).expect("roundtrip 2D simplex mesh");
        rt.check().expect("roundtrip SimplexMesh<2>::check passes");

        // Refine/adapt with adapt2 and verify SimplexMesh::check + measure preservation.
        let tri_base = rect_tri(&xs, &ys);
        let base_area = CellMesh::from_simplex(&tri_base).total_measure();
        let inv_h2 = 1.0 / (h_scale * h_scale);
        let metric = move |_: &[f64; 2]| -> Sym<2> { [[inv_h2, 0.0], [0.0, inv_h2]] };
        let opt = AdaptOptions { max_iter: 2, max_cells: Some(200), ..AdaptOptions::default() };
        let adapted = adapt2(&tri_base, &metric, &opt).expect("adapt2 succeeds");
        adapted.mesh.check().expect("adapted SimplexMesh<2>::check passes");
        simplex_from_adapt(&adapted.mesh).check().expect("adapted NyxMesh::check passes");
        let adapted_area = CellMesh::from_simplex(&adapted.mesh).total_measure();
        prop_assert!((adapted_area - base_area).abs() < 1e-10 * (1.0 + base_area));
    }

    #[test]
    fn mesh3d_cell_splits_refine_and_adapt3_preserve_check_and_volume(
        xs in arb_grid_axis(2),
        ys in arb_grid_axis(2),
        zs in arb_grid_axis(2),
        mode in 0u8..3,
        h_scale in 0.7_f64..1.4,
    ) {
        let nx = xs.len() - 1;
        let ny = ys.len() - 1;
        let nz = zs.len() - 1;
        let mut cm = CellMesh::<3>::new();
        for &z in &zs {
            for &y in &ys {
                for &x in &xs {
                    cm.points.push([x, y, z]);
                }
            }
        }
        let idx = |ix: usize, iy: usize, iz: usize| -> u32 { ((iz * (ny + 1) + iy) * (nx + 1) + ix) as u32 };
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let v000 = idx(ix, iy, iz);
                    let v100 = idx(ix + 1, iy, iz);
                    let v110 = idx(ix + 1, iy + 1, iz);
                    let v010 = idx(ix, iy + 1, iz);
                    let v001 = idx(ix, iy, iz + 1);
                    let v101 = idx(ix + 1, iy, iz + 1);
                    let v111 = idx(ix + 1, iy + 1, iz + 1);
                    let v011 = idx(ix, iy + 1, iz + 1);
                    let lab = 1 + ((ix + iy + iz) as i32 % 2);
                    if mode == 0 {
                        cm.push_cell(CellType::Hex, &[v000, v100, v110, v010, v001, v101, v111, v011], lab);
                    } else {
                        cm.push_cell(CellType::Prism, &[v000, v100, v110, v001, v101, v111], lab);
                        cm.push_cell(CellType::Prism, &[v000, v110, v010, v001, v111, v011], lab);
                    }
                }
            }
        }
        if mode == 0 {
            push_3d_boundary_quad_facets(&mut cm, nx, ny, nz);
        } else {
            push_3d_prism_boundary_facets(&mut cm, nx, ny, nz);
        }
        cm.check().expect("3D Hex/Prism CellMesh valid");
        let vol_before = cm.total_measure();

        let split = cm.split_to_simplices().expect("3D split succeeds");
        split.check().expect("split SimplexMesh<3>::check passes");
        let split_nyx = simplex_from_adapt(&split);
        split_nyx.check().expect("split 3D NyxMesh::check passes");
        prop_assert!((CellMesh::from_simplex(&split).total_measure() - vol_before).abs() < 1e-10 * (1.0 + vol_before));

        // Check refine_to_size and adapt3 on brick_tet.
        let tet_base = brick_tet(&xs, &ys, &zs);
        let base_vol = CellMesh::from_simplex(&tet_base).total_measure();
        let tet_cm = CellMesh::<3>::from_simplex(&tet_base);
        let (refined_cm, _) = refine_to_size(&tet_cm, &|_| h_scale, &|mid, _, _, _| mid, 1.4, 4_000)
            .expect("refine_to_size succeeds");
        let refined = refined_cm.split_to_simplices().expect("refined split succeeds");
        refined.check().expect("refined SimplexMesh<3>::check passes");
        simplex_from_adapt(&refined).check().expect("refined 3D NyxMesh::check passes");
        prop_assert!((CellMesh::from_simplex(&refined).total_measure() - base_vol).abs() < 1e-10 * (1.0 + base_vol));

        let inv_h2 = 1.0 / (h_scale * h_scale);
        let metric = move |_: &[f64; 3]| -> Sym<3> {
            [[inv_h2, 0.0, 0.0], [0.0, inv_h2, 0.0], [0.0, 0.0, inv_h2]]
        };
        let opt = AdaptOptions3 { max_iter: 2, max_cells: Some(2_000), ..AdaptOptions3::default() };
        let adapted = adapt3(&tet_base, &metric, &TopoOptions3::default(), &opt).expect("adapt3 succeeds");
        adapted.mesh.check().expect("adapted SimplexMesh<3>::check passes");
        simplex_from_adapt(&adapted.mesh).check().expect("adapted 3D NyxMesh::check passes");
        let adapted_vol = CellMesh::from_simplex(&adapted.mesh).total_measure();
        prop_assert!((adapted_vol - base_vol).abs() < 1e-10 * (1.0 + base_vol));
    }
}

proptest! {
    #[test]
    fn csr_algebra_linearity_double_transpose_and_sum_all(
        (nrows, ncols) in (1usize..=12, 1usize..=12),
        alpha in (-8i32..=8).prop_map(|k| f64::from(k) * 0.25),
        beta in (-8i32..=8).prop_map(|k| f64::from(k) * 0.25),
        raw_a in arb_triplets(12, 12, 48),
        raw_b in arb_triplets(12, 12, 48),
        raw_c in arb_triplets(12, 12, 48),
        raw_x in prop::collection::vec((-8i32..=8).prop_map(|k| f64::from(k) * 0.25), 12..=12),
    ) {
        let map_t = |raw: &[(usize, usize, f64)]| -> Vec<(usize, usize, f64)> {
            raw.iter().map(|&(r, c, v)| (r % nrows, c % ncols, v)).collect()
        };
        let tris_a = map_t(&raw_a);
        let tris_b = map_t(&raw_b);
        let tris_c = map_t(&raw_c);
        let x = &raw_x[..ncols];

        let a = Csr::from_triplets(nrows, ncols, &tris_a).unwrap();
        let b = Csr::from_triplets(nrows, ncols, &tris_b).unwrap();
        let c = Csr::from_triplets(nrows, ncols, &tris_c).unwrap();
        prop_assert!(a.has_sorted_unique_rows());
        prop_assert_eq!(a.to_dense(), triplets_to_dense(nrows, ncols, &tris_a));

        // Double transpose exact identity.
        let att = a.transpose().transpose();
        prop_assert_eq!(att.nrows(), a.nrows());
        prop_assert_eq!(att.ncols(), a.ncols());
        prop_assert_eq!(att.indptr(), a.indptr());
        prop_assert_eq!(att.indices(), a.indices());
        prop_assert_eq!(att.values(), a.values());

        // Linearity: (alpha * A + beta * B) * x == alpha * (A * x) + beta * (B * x).
        let combo = Csr::axpby(alpha, &a, beta, &b);
        prop_assert!(combo.has_sorted_unique_rows());
        let mut y_combo = vec![0.0; nrows];
        let mut y_a = vec![0.0; nrows];
        let mut y_b = vec![0.0; nrows];
        combo.spmv(x, &mut y_combo);
        a.spmv(x, &mut y_a);
        b.spmv(x, &mut y_b);
        for i in 0..nrows {
            let rhs = alpha * y_a[i] + beta * y_b[i];
            prop_assert!((y_combo[i] - rhs).abs() <= 1e-11, "row {i}: {} != {rhs}", y_combo[i]);
        }

        // spmv_acc accumulates cleanly onto existing y.
        let mut y_acc = vec![1.5; nrows];
        a.spmv_acc(x, &mut y_acc);
        for i in 0..nrows {
            prop_assert!((y_acc[i] - (1.5 + y_a[i])).abs() <= 1e-12);
        }

        // Csr::sum_all matches pairwise addition.
        let sum3 = Csr::sum_all(vec![a.clone(), b.clone(), c.clone()]);
        let pair3 = Csr::axpby(1.0, &Csr::axpby(1.0, &a, 1.0, &b), 1.0, &c);
        prop_assert_eq!(sum3.to_dense(), pair3.to_dense());

        // Matmul agrees with sequential SpMV: (A * A^T) * v == A * (A^T * v).
        let at = a.transpose();
        let aat = a.matmul(&at);
        let mut v_mid = vec![0.0; ncols];
        let mut y_seq = vec![0.0; nrows];
        let mut y_mm = vec![0.0; nrows];
        let w: Vec<f64> = (0..nrows).map(|i| (i as f64) * 0.5 - 1.0).collect();
        at.spmv(&w, &mut v_mid);
        a.spmv(&v_mid, &mut y_seq);
        aat.spmv(&w, &mut y_mm);
        for i in 0..nrows {
            prop_assert!((y_mm[i] - y_seq[i]).abs() <= 1e-10);
        }
    }

    #[test]
    fn csr_complex_axpby_and_double_transpose(
        (nrows, ncols) in (1usize..=8, 1usize..=8),
        raw_a in arb_c64_triplets(8, 8, 24),
        raw_b in arb_c64_triplets(8, 8, 24),
        raw_x in prop::collection::vec((-4i32..=4, -4i32..=4).prop_map(|(r, i)| C64::new(f64::from(r) * 0.5, f64::from(i) * 0.5)), 8..=8),
    ) {
        let map_t = |raw: &[(usize, usize, C64)]| -> Vec<(usize, usize, C64)> {
            raw.iter().map(|&(r, c, v)| (r % nrows, c % ncols, v)).collect()
        };
        let a = Csr::from_triplets(nrows, ncols, &map_t(&raw_a)).unwrap();
        let b = Csr::from_triplets(nrows, ncols, &map_t(&raw_b)).unwrap();
        prop_assert_eq!(a.transpose().transpose().to_dense(), a.to_dense());

        let sum = Csr::axpby(C64::new(1.0, 0.0), &a, C64::new(1.0, 0.0), &b);
        let x = &raw_x[..ncols];
        let (mut ys, mut ya, mut yb) = (vec![C64::default(); nrows], vec![C64::default(); nrows], vec![C64::default(); nrows]);
        sum.spmv(x, &mut ys);
        a.spmv(x, &mut ya);
        b.spmv(x, &mut yb);
        for i in 0..nrows {
            let d = ys[i] - (ya[i] + yb[i]);
            prop_assert!(d.re.abs() <= 1e-11 && d.im.abs() <= 1e-11);
        }
    }

    #[test]
    fn formula_parse_ast_roundtrip_and_eval(
        ast in arb_ast(),
        coords in proptest::array::uniform4((-4.0_f64..4.0).prop_map(|x| (x * 4.0).round() * 0.25)),
    ) {
        let text = ast_to_str(&ast);
        let mut names = PropFormulaNames;
        let tree1 = parse_formula(&text, &mut names).expect("valid generated formula parses");
        tree1.check(&[false, false, false, false]).expect("real tree checks cleanly");

        let leaves = [
            C64::new(coords[0], 0.0),
            C64::new(coords[1], 0.0),
            C64::new(coords[2], 0.0),
            C64::new(coords[3], 0.0),
        ];
        let val1 = tree1.eval(&leaves);

        // Format parsed AST back to text and re-parse -> structural AST equality.
        let reprinted = tree_to_formula(&tree1);
        let tree2 = parse_formula(&reprinted, &mut names).expect("reprinted AST reparses");
        prop_assert_eq!(&tree2, &tree1, "AST round-trip mismatch for {}", reprinted);
        let val2 = tree2.eval(&leaves);
        prop_assert!((val1.re - val2.re).abs() <= 1e-12 && (val1.im - val2.im).abs() <= 1e-12);

        // Quick consistency check on BTreeMap unused import avoidance.
        let _ = BTreeMap::<String, f64>::new();
    }
}
