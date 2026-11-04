//! Storage study: `.fmsh` against plain binary arrays (and, when the tools are installed, against
//! gzip/xz/zstd run on the plain arrays; those are comparisons only, nothing links to them).
//! Every file is decoded again and compared bit for bit.
//!
//! Usage: `cargo run --release --example fmsh_size`

use femcore::adapt::drive::{AdaptOptions, adapt2};
use femcore::adapt::field::BackgroundMetric2;
use femcore::adapt::metric::*;
use femcore::meshing::pml_box::{PmlSpec, build_pml};
use femcore::mesh::cells::CellMesh;
use femcore::mesh::fixtures::{Field as _, Gaussian, box_tet, box_tri, rect_tri};
use femcore::mesh::io::fmsh::{Field, Step, decode_dataset, encode_steps};
use std::io::Write;
use std::process::{Command, Stdio};
use std::time::Instant;

fn plain<const D: usize>(m: &CellMesh<D>, steps: &[Step]) -> Vec<u8> {
    let mut o = Vec::new();
    for p in &m.points {
        for v in p {
            o.extend_from_slice(&v.to_le_bytes());
        }
    }
    for v in &m.conn {
        o.extend_from_slice(&v.to_le_bytes());
    }
    o.extend(m.types.iter().map(|t| t.code()));
    for v in &m.label {
        o.extend_from_slice(&v.to_le_bytes());
    }
    for v in &m.facet_conn {
        o.extend_from_slice(&v.to_le_bytes());
    }
    for v in &m.facet_label {
        o.extend_from_slice(&v.to_le_bytes());
    }
    for s in steps {
        for f in &s.fields {
            for v in &f.data {
                o.extend_from_slice(&v.to_le_bytes());
            }
        }
    }
    o
}

fn external(cmd: &str, args: &[&str], data: &[u8]) -> Option<usize> {
    let mut c =
        Command::new(cmd).args(args).stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::null()).spawn().ok()?;
    let mut si = c.stdin.take()?;
    let d = data.to_vec();
    let h = std::thread::spawn(move || {
        let _ = si.write_all(&d);
    });
    let out = c.wait_with_output().ok()?;
    let _ = h.join();
    out.status.success().then_some(out.stdout.len())
}

fn report<const D: usize>(name: &str, m: &CellMesh<D>, steps: &[Step]) {
    let raw = plain(m, steps);
    let t = Instant::now();
    let enc = encode_steps(m, steps);
    let te = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let d = decode_dataset::<D>(&enc).unwrap();
    let td = t.elapsed().as_secs_f64();
    assert!(d.mesh.points.iter().zip(&m.points).all(|(a, b)| (0..D).all(|k| a[k].to_bits() == b[k].to_bits())));
    assert_eq!((&d.mesh.conn, &d.mesh.label, &d.mesh.types), (&m.conn, &m.label, &m.types));
    for (a, b) in d.steps.iter().zip(steps) {
        for (f, g) in a.fields.iter().zip(&b.fields) {
            assert!(f.data.iter().zip(&g.data).all(|(x, y)| x.to_bits() == y.to_bits()));
        }
    }
    let ext = |c: &str, a: &[&str]| external(c, a, &raw).map(|n| n.to_string()).unwrap_or_else(|| "n/a".into());
    println!(
        "| {name} | {} | {} | {} | {:.2}x | {} | {} | {} | {:.2} / {:.2} |",
        m.n_cells(),
        raw.len(),
        enc.len(),
        raw.len() as f64 / enc.len() as f64,
        ext("gzip", &["-9", "-c"]),
        ext("xz", &["-9e", "-c"]),
        ext("zstd", &["-19", "-c"]),
        te,
        td
    );
}

fn main() {
    println!("| case | cells | plain B | fmsh B | ratio | gzip -9 B | xz -9e B | zstd -19 B | enc / dec s |");
    println!("|---|---|---|---|---|---|---|---|---|");

    // structured 2D mesh, decimal coordinates
    let s = CellMesh::from_simplex(&box_tri([[0.0, 1.0], [0.0, 1.0]], 150, 150));
    report("structured tri 150x150", &s, &[]);

    // adapted (unstructured) 2D mesh with an analytic field and a fake complex field
    let g = Gaussian::<2> { centre: [0.5, 0.5], width: [0.12, 0.03] };
    let lin = |n: usize| (0..=n).map(|i| i as f64 / n as f64).collect::<Vec<_>>();
    let bg = rect_tri(&lin(300), &lin(300));
    let hess: Vec<Sym<2>> = bg.points.iter().map(|p| g.hessian(p)).collect();
    let vol = bg.vertex_volumes();
    let opt = LpOptions {
        p: 2.0,
        target_complexity: complexity_for_elements(30000.0, 2),
        hmin: 1e-4,
        hmax: 0.5,
        ..Default::default()
    };
    let mut lp = lp_metrics(&hess, &vol, &opt);
    grade(&mut lp.metrics, &bg.points, &bg.edges(), 0.4, 50);
    let bgf = BackgroundMetric2::new(&bg, &lp.metrics);
    let start = box_tri([[0.0, 1.0], [0.0, 1.0]], 20, 20);
    let a = adapt2(&start, &bgf, &AdaptOptions::default()).unwrap();
    let am = CellMesh::from_simplex(&a.mesh);
    let u: Vec<f64> = am.points.iter().map(|p| g.value(p)).collect();
    let ev: Vec<f64> =
        am.points.iter().flat_map(|p| [g.value(p) * (9.0 * p[0]).cos(), g.value(p) * (9.0 * p[0]).sin()]).collect();
    let steps = vec![Step {
        value: 1.0,
        fields: vec![
            Field::scalar_on_points("u", u),
            Field { name: "E".into(), location: femcore::mesh::io::fmsh::Location::Point, components: 2, data: ev },
        ],
    }];
    report("adapted tri, mesh only", &am, &[]);
    report("adapted tri + 3 f64 point fields", &am, &steps);

    // 3D PML mesh (hexes, prisms, tets)
    let phys = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [16, 16, 16]);
    let pm = build_pml(&phys, &PmlSpec::<3>::uniform(0.5, 6)).unwrap();
    report("3D PML (hex/prism/tet), native", &pm.mesh, &[]);
    let split = CellMesh::from_simplex(&pm.mesh.split_to_simplices().unwrap());
    report("3D PML split to tets", &split, &[]);
}
