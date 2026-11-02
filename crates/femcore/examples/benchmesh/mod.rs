//! Meshes of boxes for the timing studies, made from geom models (`geomdrive::plate` / `block`) with the element size of the finer of
//! the given counts. The boundary numbers are the ones the studies have always used: `rect` 1 x = 0, 2 y = 0, 3 y = ly, 4 x = lx (the
//! second axis is z of the x-z section); `block` 1 x = 0, 2 y = 0, 3 z = 0, 4 z = lz, 5 y = ly, 6 x = lx.
#![allow(dead_code)]
use femcore::mesh::SimplexMesh;

fn relabel(mesh: &mut SimplexMesh, map: &[(&[i32], i32)]) {
    let old = mesh.facet_labels.clone();
    for (f, l) in mesh.facet_labels.iter_mut().enumerate() {
        *l = map.iter().find(|(from, _)| from.contains(&old[f])).map_or(old[f], |(_, to)| *to);
    }
}

pub fn rect(lx: f64, ly: f64, nx: usize, ny: usize) -> SimplexMesh {
    let p = femcore::geomdrive::plate(lx, ly, (lx / nx as f64).min(ly / ny as f64)).expect("a plate");
    let mut m = p.mesh;
    relabel(&mut m, &[(&p.xmin, 1), (&p.zmin, 2), (&p.zmax, 3), (&p.xmax, 4)]);
    m
}

pub fn block(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> SimplexMesh {
    let b = femcore::geomdrive::block(lx, ly, lz, (lx / nx as f64).min(lz / nz as f64), ny).expect("a block");
    let mut m = b.mesh;
    relabel(&mut m, &[(&b.xmin, 1), (&b.ymin, 2), (&b.zmin, 3), (&b.zmax, 4), (&b.ymax, 5), (&b.xmax, 6)]);
    m
}
