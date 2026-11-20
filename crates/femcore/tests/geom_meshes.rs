// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The sweep of a section along y. (The stress concentration of a plate with a hole, which needed the solid solver, moved to `nyx-sms`.)
use femcore::meshing::mesher::{extrude_y, Layer};
use std::collections::HashMap;

fn assert_conforming(m: &femcore::mesh::SimplexMesh) {
    let mut faces: HashMap<[u32; 3], usize> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            f.sort_unstable();
            *faces.entry([f[0], f[1], f[2]]).or_insert(0) += 1;
        }
    }
    assert!(faces.values().all(|&n| n <= 2));
    let mut boundary: Vec<[u32; 3]> = faces.iter().filter(|(_, &n)| n == 1).map(|(f, _)| *f).collect();
    let mut given: Vec<[u32; 3]> = (0..m.n_facets())
        .map(|f| {
            let mut v = m.facet(f).to_vec();
            v.sort_unstable();
            [v[0], v[1], v[2]]
        })
        .collect();
    boundary.sort();
    given.sort();
    assert_eq!(boundary, given);
}

#[test]
fn a_section_swept_along_y_is_conforming() {
    let p = femcore::meshing::geomdrive::plate(1.0, 0.4, 0.1).unwrap();
    let layers = vec![Layer { thickness: 0.2, n: 2, domain_offset: 0, fill: None }, Layer { thickness: 0.1, n: 1, domain_offset: 10, fill: None }];
    let m3 = extrude_y(&p.mesh, &layers, 0.0).unwrap();
    assert_eq!(m3.n_cells(), p.mesh.n_cells() * 3 * 3);
    let vol: f64 = femcore::fem::geometry(&m3).vol.iter().sum();
    assert!((vol - 1.0 * 0.4 * 0.3).abs() < 1e-12);
    assert_conforming(&m3);
    assert_eq!(m3.domain_labels(), vec![1, 11]);
    // the walls keep the labels of the section, the ends are one more and two more than the largest
    let nb = p.mesh.facet_labels.iter().copied().max().unwrap();
    let mut want = p.mesh.boundary_labels();
    want.extend([nb + 1, nb + 2]);
    assert_eq!(m3.boundary_labels(), want);
}
