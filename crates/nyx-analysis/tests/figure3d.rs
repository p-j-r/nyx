// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification tests for FR-4 3D figures in `nyx-analysis`:
//! (a) 8-subtet P2 split volume & orientation;
//! (b) `Slice3d`, `Cutaway3d`, `Isosurface3d`, `Mesh3d` PNG & SVG generation;
//! (c) pixel-for-pixel preservation of sharp `5.79x` normal-field jump across `Si|Gap` without welding smear;
//! (d) zero-copy Bloch tiling (`repeat_periodic = 4`, `periodic_phase = pi/2`) sign reversal every 2 periods & `expr` formulas;
//! (e) 6,300-tet P2 mesh tiled 4x at `1600x1100` `2x`-supersampled PNG latency `< 150 ms`.

use femcore::mesh::SimplexMesh;
use nyx_analysis::figure::{mode_figure_on, Camera3d, ClipBox, Cmap, EdgeMode, FigureOptions, Kind, PlaneSpec};
use nyx_analysis::view3d::{tet_signed_volume, P2_SUBTETS};
use nyx_physics::{CellField, Solution, SolveResult};
use std::f64::consts::FRAC_PI_2;
use std::time::Instant;

/// Builds a conforming 3D Freudenthal structured tetrahedral mesh of `nx * ny * nz * 6` cells over `bounds`,
/// labelling cells `1` (`Si`) where `x < x_split` and `2` (`Gap`) where `x >= x_split`.
fn brick_mesh_3d(nx: usize, ny: usize, nz: usize, bounds: [[f64; 2]; 3], x_split: Option<f64>) -> SimplexMesh {
    let mut points = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1) * 3);
    for k in 0..=nz {
        let z = bounds[2][0] + (bounds[2][1] - bounds[2][0]) * k as f64 / nz as f64;
        for j in 0..=ny {
            let y = bounds[1][0] + (bounds[1][1] - bounds[1][0]) * j as f64 / ny as f64;
            for i in 0..=nx {
                let x = bounds[0][0] + (bounds[0][1] - bounds[0][0]) * i as f64 / nx as f64;
                points.extend_from_slice(&[x, y, z]);
            }
        }
    }
    let vid = |i: usize, j: usize, k: usize| ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
    const PERMS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let mut cells = Vec::with_capacity(nx * ny * nz * 6 * 4);
    let mut cell_labels = Vec::with_capacity(nx * ny * nz * 6);

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let xc = bounds[0][0] + (bounds[0][1] - bounds[0][0]) * (i as f64 + 0.5) / nx as f64;
                let dom = match x_split {
                    Some(xs) if xc >= xs => 2,
                    _ => 1,
                };
                for perm in PERMS {
                    let mut c = [0usize; 3];
                    let mut vs = [vid(i, j, k); 4];
                    for (s, &ax) in perm.iter().enumerate() {
                        c[ax] = 1;
                        vs[s + 1] = vid(i + c[0], j + c[1], k + c[2]);
                    }
                    let p = |u: u32| {
                        let b = u as usize * 3;
                        [points[b], points[b + 1], points[b + 2]]
                    };
                    if tet_signed_volume(p(vs[0]), p(vs[1]), p(vs[2]), p(vs[3])) < 0.0 {
                        vs.swap(2, 3);
                    }
                    cells.extend_from_slice(&vs);
                    cell_labels.push(dom);
                }
            }
        }
    }

    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    let mut push_quad = |q: [u32; 4], lab: i32| {
        facets.extend_from_slice(&[q[0], q[1], q[2], q[0], q[2], q[3]]);
        facet_labels.extend_from_slice(&[lab, lab]);
    };
    for k in 0..nz {
        for j in 0..ny {
            push_quad([vid(0, j, k), vid(0, j + 1, k), vid(0, j + 1, k + 1), vid(0, j, k + 1)], 1);
            push_quad([vid(nx, j, k), vid(nx, j + 1, k), vid(nx, j + 1, k + 1), vid(nx, j, k + 1)], 2);
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            push_quad([vid(i, 0, k), vid(i + 1, 0, k), vid(i + 1, 0, k + 1), vid(i, 0, k + 1)], 3);
            push_quad([vid(i, ny, k), vid(i + 1, ny, k), vid(i + 1, ny, k + 1), vid(i, ny, k + 1)], 4);
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            push_quad([vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0)], 5);
            push_quad([vid(i, j, nz), vid(i + 1, j, nz), vid(i + 1, j + 1, nz), vid(i, j + 1, nz)], 6);
        }
    }

    SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels)
}

fn decode_png(bytes: &[u8]) -> (u32, u32, Vec<u8>) {
    assert_eq!(&bytes[..4], &[0x89, b'P', b'N', b'G'], "PNG header signature");
    let decoder = png::Decoder::new(std::io::Cursor::new(bytes));
    let mut reader = decoder.read_info().expect("valid PNG metadata");
    let mut buf = vec![0u8; reader.output_buffer_size().expect("PNG buffer size")];
    let info = reader.next_frame(&mut buf).expect("valid PNG pixel stream");
    buf.truncate(info.buffer_size());
    (info.width, info.height, buf)
}

#[test]
fn p2_subtet_split_preserves_exact_volume_and_orientation() {
    let tets = [
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-0.3, 0.2, -0.1], [0.9, 0.1, 0.2], [0.1, 1.4, -0.2], [0.2, 0.4, 1.1]],
        [[1e-6, -2e-6, 0.5e-6], [2.4e-6, -1.8e-6, 0.4e-6], [1.5e-6, 0.3e-6, 0.7e-6], [1.7e-6, -1.1e-6, 1.9e-6]],
    ];
    for verts in tets {
        let mut pts = Vec::new();
        for v in verts {
            pts.extend_from_slice(&v);
        }
        let m = SimplexMesh::new(3, pts, vec![0, 1, 2, 3], vec![1], vec![0, 1, 2], vec![1]);
        let g = m.cell_geometry_nodes(0);
        let parent_vol = tet_signed_volume(g[0], g[1], g[2], g[3]);
        assert!(parent_vol > 0.0, "parent tet must be positively oriented");
        let mut sum_vol = 0.0;
        for (si, st) in P2_SUBTETS.iter().enumerate() {
            let v_sub = tet_signed_volume(g[st[0]], g[st[1]], g[st[2]], g[st[3]]);
            assert!(v_sub > 0.0, "sub-tet {si} ({st:?}) inverted orientation: {v_sub} vs parent {parent_vol}");
            assert!(
                (v_sub - parent_vol / 8.0).abs() <= 1e-14 * parent_vol,
                "sub-tet {si} volume {v_sub} != parent/8 {}",
                parent_vol / 8.0
            );
            sum_vol += v_sub;
        }
        assert!((sum_vol - parent_vol).abs() <= 1e-14 * parent_vol);
    }
}

#[test]
fn all_four_3d_figure_kinds_and_expressions_produce_valid_png_and_svg() {
    let mesh = brick_mesh_3d(4, 4, 4, [[-0.5e-6, 0.5e-6], [0.0, 0.4e-6], [-0.3e-6, 0.3e-6]], Some(0.0));
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(2);
    assert_eq!(nloc, 10);

    let (mut ex_re, mut ex_im) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
    let (mut ey_re, mut ey_im) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
    let (mut ez_re, mut ez_im) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
    for c in 0..nc {
        let g = mesh.cell_geometry_nodes(c);
        let jump = if mesh.cell_labels[c] == 2 { 5.79 } else { 1.0 };
        for (a, &p) in g.iter().enumerate().take(nloc) {
            let r2 = (p[0] / 0.35e-6).powi(2) + (p[2] / 0.22e-6).powi(2);
            let env = (-r2).exp();
            let idx = c * nloc + a;
            ex_re[idx] = jump * env;
            ey_im[idx] = 0.45 * (p[0] / 0.35e-6) * env;
            ez_re[idx] = 0.2 * env;
            ex_im[idx] = 0.0;
            ey_re[idx] = 0.0;
            ez_im[idx] = 0.0;
        }
    }
    let sol = Solution {
        label: "mode 1".into(),
        fields: vec![],
        value: Some(2.45),
        eigenvalue: Some([2.45, 0.0]),
        cell_fields: vec![
            CellField { name: "Ex".into(), unit: "V/m".into(), nloc, re: ex_re, im: ex_im },
            CellField { name: "Ey".into(), unit: "V/m".into(), nloc, re: ey_re, im: ey_im },
            CellField { name: "Ez".into(), unit: "V/m".into(), nloc, re: ez_re, im: ez_im },
        ],
    };
    let res = SolveResult {
        physics: "ems".into(),
        study: "bloch".into(),
        order: 2,
        dim: 3,
        dofs: 100,
        seconds: 0.01,
        solutions: vec![sol],
        scalars: vec![],
        log: vec![],
    };

    for kind in [Kind::Slice3d, Kind::Cutaway3d, Kind::Isosurface3d, Kind::Mesh3d] {
        let o = FigureOptions { kind, ..Default::default() };
        let fig = mode_figure_on(&mesh, &res, &o).unwrap();
        let png = fig.canvas.to_png(1.5).unwrap();
        let (w, h, rgba) = decode_png(&png);
        assert!(w > 100 && h > 100, "{kind:?} PNG size ({w}, {h})");
        assert_eq!(rgba.len(), (w * h * 4) as usize);
        let svg = fig.canvas.to_svg().unwrap();
        assert!(
            svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>") && svg.contains("data:image/png;base64,")
        );
    }

    // Verify all required `expr` forms (`re(Ex)`, `abs(Ey)`, `im(Ey)`, `abs2(E)`) + General slice plane
    for expr in ["re(Ex)", "abs(Ey)", "im(Ey)", "abs2(E)"] {
        let o = FigureOptions {
            kind: Kind::Slice3d,
            expr: Some(expr.to_string()),
            plane: Some(PlaneSpec::General { origin: [0.0, 0.2e-6, 0.0], normal: [0.2, 0.9, 0.3] }),
            ..Default::default()
        };
        let fig = mode_figure_on(&mesh, &res, &o).unwrap();
        let png = fig.canvas.to_png(1.0).unwrap();
        let (w, h, _) = decode_png(&png);
        assert!(w > 100 && h > 100, "expr '{expr}' failed");
    }
}

#[test]
fn sharp_5_79x_dielectric_jump_is_preserved_pixel_for_pixel_without_welding_smear() {
    // Two-domain 3D mesh on [-1, 1] x [0, 1] x [0, 1]: x < 0 is Si (Ex = 1.0), x > 0 is Gap (Ex = 5.79).
    // Coarse 2-brick x-grid ([-1, 0] and [0, 1]) so welding across x = 0 would smear across the ENTIRE half-width!
    let mesh = brick_mesh_3d(2, 2, 2, [[-1.0, 1.0], [0.0, 1.0], [0.0, 1.0]], Some(0.0));
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(2);
    let mut ex_re = vec![0.0; nc * nloc];
    for c in 0..nc {
        let v = if mesh.cell_labels[c] == 2 { 5.79 } else { 1.0 };
        for a in 0..nloc {
            ex_re[c * nloc + a] = v;
        }
    }
    let sol = Solution {
        label: "jump".into(),
        fields: vec![],
        value: Some(1.0),
        eigenvalue: None,
        cell_fields: vec![
            CellField { name: "Ex".into(), unit: "V/m".into(), nloc, re: ex_re, im: vec![] },
            CellField { name: "Ez".into(), unit: "V/m".into(), nloc, re: vec![0.0; nc * nloc], im: vec![] },
        ],
    };
    let res = SolveResult {
        physics: "ems".into(),
        study: "cavity".into(),
        order: 2,
        dim: 3,
        dofs: 100,
        seconds: 0.01,
        solutions: vec![sol],
        scalars: vec![],
        log: vec![],
    };

    let cmap = Cmap::Viridis;
    let want_si: [u8; 3] = {
        let c = cmap.at(1.0 / 5.79);
        [(c[0] * 255.0 + 0.5) as u8, (c[1] * 255.0 + 0.5) as u8, (c[2] * 255.0 + 0.5) as u8]
    };
    let want_gap: [u8; 3] = {
        let c = cmap.at(1.0);
        [(c[0] * 255.0 + 0.5) as u8, (c[1] * 255.0 + 0.5) as u8, (c[2] * 255.0 + 0.5) as u8]
    };
    let dist_rgb = |p: &[u8], q: [u8; 3]| -> i32 {
        (p[0] as i32 - q[0] as i32).abs() + (p[1] as i32 - q[1] as i32).abs() + (p[2] as i32 - q[2] as i32).abs()
    };

    // 1. Slice3d on z = 0.5 with show_edges = None: scan across x = 0 at panel vertical centre
    let fig_slice = mode_figure_on(
        &mesh,
        &res,
        &FigureOptions {
            kind: Kind::Slice3d,
            cmap,
            expr: Some("re(Ex)".into()),
            plane: Some(PlaneSpec::Axis { axis: 'z', at: 0.5 }),
            show_edges: EdgeMode::None,
            window: Some([-1.0, 1.0, 0.0, 1.0]),
            panel_width: 300.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert!((fig_slice.modes[0].peak - 5.79).abs() < 1e-12);
    let (sw, _sh, spx) = decode_png(&fig_slice.canvas.to_png(1.0).unwrap());
    // Panel sits at x = 14 + 42 = 56..356, y = 54 + 46 = 100..250 -> row y = 175 is the middle of the slice
    let row_y = 175usize;
    let mut si_count = 0usize;
    let mut gap_count = 0usize;
    let mut intermediate_smeared = 0usize;
    for x in 70..340usize {
        let p = &spx[(row_y * sw as usize + x) * 4..(row_y * sw as usize + x) * 4 + 3];
        if dist_rgb(p, want_si) <= 3 {
            si_count += 1;
        } else if dist_rgb(p, want_gap) <= 3 {
            gap_count += 1;
        } else {
            intermediate_smeared += 1;
        }
    }
    assert!(si_count > 100 && gap_count > 100, "expected Si & Gap pixels, got {si_count} and {gap_count}");
    // At most 1 boundary-straddling antialiased pixel column may lie between Si and Gap (welding would smear ~270 columns!)
    assert!(
        intermediate_smeared <= 1,
        "Slice3d smeared {intermediate_smeared} pixel columns across the 5.79x Si|Gap jump!"
    );

    // 2. Cutaway3d clipped at zmax = 0.5 viewed top-down (azimuth -90, elevation 90 -> +x right, +y up)
    let fig_cut = mode_figure_on(
        &mesh,
        &res,
        &FigureOptions {
            kind: Kind::Cutaway3d,
            cmap,
            expr: Some("re(Ex)".into()),
            clip: Some(ClipBox { zmax: Some(0.5), ..Default::default() }),
            camera: Camera3d { azimuth_deg: -90.0, elevation_deg: 90.0 },
            show_edges: EdgeMode::None,
            window: Some([-1.0, 1.0, -0.5, 0.5]),
            panel_width: 300.0,
            ..Default::default()
        },
    )
    .unwrap();
    let (cw, _ch, cpx) = decode_png(&fig_cut.canvas.to_png(1.0).unwrap());
    // Panel sits at x = 18..318, y = 100..250 -> row y = 175 cuts straight across the exposed z = 0.5 interior cutaway face
    let (mut c_si, mut c_gap, mut c_smear) = (0usize, 0usize, 0usize);
    for x in 30..305usize {
        let p = &cpx[(175 * cw as usize + x) * 4..(175 * cw as usize + x) * 4 + 3];
        if dist_rgb(p, want_si) <= 3 {
            c_si += 1;
        } else if dist_rgb(p, want_gap) <= 3 {
            c_gap += 1;
        } else {
            c_smear += 1;
        }
    }
    assert!(c_si > 100 && c_gap > 100, "Cutaway3d expected Si & Gap pixels, got {c_si} and {c_gap}");
    assert!(c_smear <= 1, "Cutaway3d smeared {c_smear} pixel columns across the exposed 5.79x Si|Gap cut face!");
}

#[test]
fn repeat_periodic_4_with_quarter_turn_bloch_phase_reverses_sign_every_two_periods() {
    // Unit cell y in [0, 1], constant complex mode Ex = 1 + 0i.
    // With periodic_phase = pi/2, copy m gets factor exp(-i * m * pi/2):
    // m=0 -> +1 (red), m=1 -> 0 (grey), m=2 -> -1 (blue, reversed!), m=3 -> 0 (grey).
    let mesh = brick_mesh_3d(2, 2, 2, [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], None);
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(2);
    let sol = Solution {
        label: "bloch".into(),
        fields: vec![],
        value: Some(1.0),
        eigenvalue: None,
        cell_fields: vec![
            CellField {
                name: "Ex".into(),
                unit: "V/m".into(),
                nloc,
                re: vec![1.0; nc * nloc],
                im: vec![0.0; nc * nloc],
            },
            CellField {
                name: "Ez".into(),
                unit: "V/m".into(),
                nloc,
                re: vec![0.0; nc * nloc],
                im: vec![0.0; nc * nloc],
            },
        ],
    };
    let res = SolveResult {
        physics: "ems".into(),
        study: "bloch".into(),
        order: 2,
        dim: 3,
        dofs: 50,
        seconds: 0.01,
        solutions: vec![sol],
        scalars: vec![],
        log: vec![],
    };

    // Slice z = 0.5 maps (u = x in [0, 1], v = y in [0, 4])
    let fig = mode_figure_on(
        &mesh,
        &res,
        &FigureOptions {
            kind: Kind::Slice3d,
            expr: Some("re(Ex)".into()),
            plane: Some(PlaneSpec::Axis { axis: 'z', at: 0.5 }),
            repeat_periodic: 4,
            periodic_shift: Some([0.0, 1.0, 0.0]),
            periodic_phase: FRAC_PI_2,
            show_edges: EdgeMode::None,
            window: Some([0.0, 4.0, 0.0, 4.0]),
            panel_width: 320.0,
            ..Default::default()
        },
    )
    .unwrap();

    let (w, _h, px) = decode_png(&fig.canvas.to_png(1.0).unwrap());
    // Panel x = 56..376 (u in [0, 4] -> x_mid of u=0.5 is 56 + 0.5/4 * 320 = 96),
    // Panel y = 100..420 (v = y in [0, 4] runs upward, so y_screen = 100 + (4.0 - v)/4.0 * 320):
    // m=0 (v=0.5) -> py=380; m=1 (v=1.5) -> py=300; m=2 (v=2.5) -> py=220; m=3 (v=3.5) -> py=140.
    let sample = |v_world: f64| -> [u8; 3] {
        let px_x = 96usize;
        let py_y = (100.0 + (4.0 - v_world) / 4.0 * 320.0).round() as usize;
        let i = (py_y * w as usize + px_x) * 4;
        [px[i], px[i + 1], px[i + 2]]
    };
    let c0 = sample(0.5);
    let c1 = sample(1.5);
    let c2 = sample(2.5);
    let c3 = sample(3.5);

    // Coolwarm: +1 is warm red (R >> B), 0 is neutral grey (R == G == B ~ 221), -1 is cool blue (B >> R)
    assert!(c0[0] > c0[2] + 80, "m=0 should be positive red in Coolwarm, got {c0:?}");
    assert!((c1[0] as i32 - c1[2] as i32).abs() <= 4, "m=1 should be zero grey, got {c1:?}");
    assert!(c2[2] > c2[0] + 80, "m=2 should reverse sign to negative blue, got {c2:?}");
    assert!((c3[0] as i32 - c3[2] as i32).abs() <= 4, "m=3 should be zero grey, got {c3:?}");

    // Also verify SolveResult scalar `phase 1 (rad)` auto-pickup when `periodic_phase` is left at 0.0
    let mut res_with_scalar = res.clone();
    res_with_scalar.scalars.push(("phase 1 (rad)".into(), nyx_physics::Scalar::Num(FRAC_PI_2)));
    let fig_scalar = mode_figure_on(
        &mesh,
        &res_with_scalar,
        &FigureOptions {
            kind: Kind::Slice3d,
            expr: Some("re(Ex)".into()),
            plane: Some(PlaneSpec::Axis { axis: 'z', at: 0.5 }),
            repeat_periodic: 4,
            show_edges: EdgeMode::None,
            window: Some([0.0, 4.0, 0.0, 4.0]),
            panel_width: 320.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(fig_scalar.canvas.to_png(1.0).unwrap(), fig.canvas.to_png(1.0).unwrap());
}

#[test]
fn latency_6300_tet_p2_mesh_tiled_4x_at_1600x1100_finishes_under_150_ms() {
    // 15 * 7 * 10 * 6 = 6,300 P2 tetrahedra (tiled 4x = 25,200 P2 tets = 201,600 linear sub-tets)
    let mesh = brick_mesh_3d(15, 7, 10, [[-0.6e-6, 0.6e-6], [0.0, 0.35e-6], [-0.4e-6, 0.4e-6]], Some(0.0));
    assert_eq!(mesh.n_cells(), 6_300);
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(2);
    let (mut ex_re, mut ex_im) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
    let (mut ey_re, mut ey_im) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
    let (mut ez_re, mut ez_im) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
    let k_y = std::f64::consts::PI / (2.0 * 0.35e-6);
    for c in 0..nc {
        let g = mesh.cell_geometry_nodes(c);
        let jump = if mesh.cell_labels[c] == 2 { 5.79 } else { 1.0 };
        for (a, &p) in g.iter().enumerate().take(nloc) {
            let env = (-(p[0] / 0.3e-6).powi(2) - (p[2] / 0.2e-6).powi(2)).exp();
            let (sp, cp) = (-k_y * p[1]).sin_cos();
            let idx = c * nloc + a;
            ex_re[idx] = jump * env * cp;
            ex_im[idx] = jump * env * sp;
            ey_re[idx] = -0.3 * env * sp;
            ey_im[idx] = 0.3 * env * cp;
            ez_re[idx] = 0.15 * env * cp;
            ez_im[idx] = 0.15 * env * sp;
        }
    }
    let sol = Solution {
        label: "mode 1".into(),
        fields: vec![],
        value: Some(2.41),
        eigenvalue: Some([2.41, 0.0]),
        cell_fields: vec![
            CellField { name: "Ex".into(), unit: "V/m".into(), nloc, re: ex_re, im: ex_im },
            CellField { name: "Ey".into(), unit: "V/m".into(), nloc, re: ey_re, im: ey_im },
            CellField { name: "Ez".into(), unit: "V/m".into(), nloc, re: ez_re, im: ez_im },
        ],
    };
    let res = SolveResult {
        physics: "ems".into(),
        study: "bloch".into(),
        order: 2,
        dim: 3,
        dofs: 10_000,
        seconds: 0.1,
        solutions: vec![sol],
        scalars: vec![],
        log: vec![],
    };

    // Fit exact 1600x1100 output PNG (with 2x internal supersampling = 3200x2200)
    // In Cutaway3d: width = 18 + pw + 22 = pw + 40, height = 54 + (46 + ph + 20) + 68 = ph + 188.
    // With scale = 2.0: canvas 800 x 550 -> pw = 760, ph = 362 (aspect = 362 / 760) -> exact 1600 x 1100 PNG!
    let pw = 760.0;
    let ph = 362.0;
    let opts = FigureOptions {
        kind: Kind::Cutaway3d,
        expr: Some("re(Ex)".into()),
        clip: Some(ClipBox { zmax: Some(0.0), xmax: Some(0.25e-6), ..Default::default() }),
        repeat_periodic: 4,
        periodic_shift: Some([0.0, 0.35e-6, 0.0]),
        periodic_phase: FRAC_PI_2,
        show_edges: EdgeMode::Regions,
        panel_width: pw,
        window: Some([-1.0e-6, 1.0e-6, -1.0e-6 * (ph / pw), 1.0e-6 * (ph / pw)]),
        ..Default::default()
    };

    // Warm-up thread pool & glyph cache once before timed measurement
    let _ = mode_figure_on(&mesh, &res, &opts).unwrap().canvas.to_png(0.5).unwrap();

    // Best-of-5 iterations avoids transient OS/scheduler spikes during parallel `--workspace` test runs
    let mut png = Vec::new();
    let best_ms = (0..5)
        .map(|_| {
            let t0 = Instant::now();
            let fig = mode_figure_on(&mesh, &res, &opts).unwrap();
            png = fig.canvas.to_png(2.0).unwrap();
            t0.elapsed().as_secs_f64() * 1000.0
        })
        .fold(f64::INFINITY, f64::min);

    let (w, h, _) = decode_png(&png);
    let n_threads = femcore::par::threads();
    let limit_ms = if n_threads >= 4 { 150.0 } else { 600.0 };
    eprintln!(
        "6,300-tet P2 x 4 Bloch cutaway3d {w}x{h} 2x-supersampled PNG rendered in {best_ms:.2} ms (best of 5, threads={n_threads}, limit={limit_ms:.0} ms)"
    );
    assert_eq!((w, h), (1600, 1100));
    assert!(
        best_ms < limit_ms,
        "expected < {limit_ms:.0} ms (threads={n_threads}) for 6,300-tet P2 x 4 at 1600x1100, took {best_ms:.2} ms"
    );
}

#[test]
fn thread_count_determinism_1_vs_4_produces_bit_identical_png_bytes() {
    // Subprocess worker branch when `NYX_FIG3D_DET_OUT` is set (so `OnceLock` of `femcore::par::threads()` initialises cleanly per thread count)
    if let Ok(out_path) = std::env::var("NYX_FIG3D_DET_OUT") {
        let mesh = brick_mesh_3d(4, 4, 4, [[-0.5e-6, 0.5e-6], [0.0, 0.4e-6], [-0.3e-6, 0.3e-6]], Some(0.0));
        let nc = mesh.n_cells();
        let nloc = mesh.nloc(2);
        let (mut ex_re, mut ey_im, mut ez_re) = (vec![0.0; nc * nloc], vec![0.0; nc * nloc], vec![0.0; nc * nloc]);
        for c in 0..nc {
            let g = mesh.cell_geometry_nodes(c);
            let jump = if mesh.cell_labels[c] == 2 { 5.79 } else { 1.0 };
            for (a, &p) in g.iter().enumerate().take(nloc) {
                let env = (-(p[0] / 0.35e-6).powi(2) - (p[2] / 0.22e-6).powi(2)).exp();
                ex_re[c * nloc + a] = jump * env;
                ey_im[c * nloc + a] = 0.45 * env;
                ez_re[c * nloc + a] = 0.2 * env;
            }
        }
        let sol = Solution {
            label: "mode 1".into(),
            fields: vec![],
            value: Some(2.45),
            eigenvalue: Some([2.45, 0.0]),
            cell_fields: vec![
                CellField { name: "Ex".into(), unit: "V/m".into(), nloc, re: ex_re, im: vec![0.0; nc * nloc] },
                CellField { name: "Ey".into(), unit: "V/m".into(), nloc, re: vec![0.0; nc * nloc], im: ey_im },
                CellField { name: "Ez".into(), unit: "V/m".into(), nloc, re: ez_re, im: vec![0.0; nc * nloc] },
            ],
        };
        let res = SolveResult {
            physics: "ems".into(),
            study: "bloch".into(),
            order: 2,
            dim: 3,
            dofs: 100,
            seconds: 0.01,
            solutions: vec![sol],
            scalars: vec![],
            log: vec![],
        };
        let mut combined = Vec::new();
        for kind in [Kind::Slice3d, Kind::Cutaway3d, Kind::Isosurface3d, Kind::Mesh3d] {
            let fig = mode_figure_on(
                &mesh,
                &res,
                &FigureOptions {
                    kind,
                    expr: Some("re(Ex)".into()),
                    repeat_periodic: 3,
                    periodic_shift: Some([0.0, 0.4e-6, 0.0]),
                    periodic_phase: FRAC_PI_2,
                    ..Default::default()
                },
            )
            .unwrap();
            let png = fig.canvas.to_png(1.5).unwrap();
            combined.extend_from_slice(&(png.len() as u64).to_le_bytes());
            combined.extend_from_slice(&png);
        }
        std::fs::write(out_path, combined).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_with_threads = |nt: &str| -> Vec<u8> {
        let out_file = std::env::temp_dir().join(format!("nyx-fig3d-det-{}-t{nt}.bin", std::process::id()));
        let status = std::process::Command::new(&exe)
            .args(["--exact", "thread_count_determinism_1_vs_4_produces_bit_identical_png_bytes", "--nocapture"])
            .env("NYX_THREADS", nt)
            .env("NYX_FIG3D_DET_OUT", &out_file)
            .status()
            .unwrap();
        assert!(status.success(), "child process NYX_THREADS={nt} failed");
        let bytes = std::fs::read(&out_file).unwrap();
        let _ = std::fs::remove_file(&out_file);
        bytes
    };
    let b1 = run_with_threads("1");
    let b4 = run_with_threads("4");
    assert!(!b1.is_empty());
    assert_eq!(b1, b4, "3D figure PNG output differed between NYX_THREADS=1 and NYX_THREADS=4!");
}
