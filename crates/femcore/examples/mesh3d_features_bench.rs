//! Timing of the 3D features of ADR 50 to 53 (no control run exists: the parent commit refuses them): wavelength sizing with the simple
//! mesher, adaptive meshing of a sphere, structured absorbing layers, the Maxwell indicators and one adaptation step of the loop on a
//! cavity. Prints `name seconds cells`.
use femcore::adapt::maxwell::{maxwell_indicators, CurlField, MaxwellOptions};
use femcore::linalg::csr::C;
use femcore::meshing::physics_size::PhysicsSizing;
use femcore::meshing::{mesh_for_sizing, mesh_model, MeshPlan, PmlMode, Sizing};
use geom::model::Model;
use std::time::Instant;

fn time<T>(name: &str, f: impl FnOnce() -> (T, usize)) {
    let t = Instant::now();
    let (_, cells) = f();
    println!("{name:34} {:9.3} s {cells:9} cells", t.elapsed().as_secs_f64());
}

struct Zero;
impl CurlField for Zero {
    fn sample(&self, _cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
        ([C::new(lam[1], 0.0), C::new(lam[2], 0.0), C::new(lam[3], 0.0)], [C::new(0.0, 0.0); 3])
    }
}

fn main() {
    let mut ball = Model::new_3d();
    let air = ball.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let b = ball.sphere("ball", (2.0, 2.0, 2.0), 0.9);
    ball.overlay(&[("air", air), ("ball", b)]);
    let waves = Sizing { order: 2, headroom: 1.0, wavelengths: [("ball".to_string(), 1.8), ("air".to_string(), 5.0)].into_iter().collect() };
    time("wavelength sizing, simple, 0.5", || {
        let (g, _) = mesh_for_sizing(&ball.list, &MeshPlan::simple(0.5), &waves, &PhysicsSizing::default()).unwrap();
        (0, g.mesh.n_cells())
    });
    time("adaptive sphere 0.6 x 0.5", || {
        let g = mesh_model(&ball.list, &MeshPlan::adaptive(0.6, 0.5)).unwrap();
        (0, g.mesh.n_cells())
    });
    time("adaptive + waves 0.8", || {
        let (g, _) = mesh_for_sizing(&ball.list, &MeshPlan::adaptive(0.8, 1.0), &Sizing { order: 2, headroom: 1.0, wavelengths: [("ball".to_string(), 1.8), ("air".to_string(), 5.0)].into_iter().collect() }, &PhysicsSizing::default()).unwrap();
        (0, g.mesh.n_cells())
    });
    let mut w = Model::new_3d();
    let air = w.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let sphere = w.sphere("ball", (0.0, 0.0, 0.0), 0.7);
    let win = w.box_("window", (-1.5, -1.5, -1.5), (1.5, 1.5, 1.5));
    w.overlay(&[("air", air), ("ball", sphere)]);
    w.pml("pml", win, &[], 0.5, 2.0, 2.0);
    time("structured PML, simple 0.1", || {
        let g = mesh_model(&w.list, &MeshPlan { pml: PmlMode::Structured, ..MeshPlan::simple(0.1) }).unwrap();
        (0, g.mesh.n_cells())
    });
    time("inline PML, simple 0.1", || {
        let g = mesh_model(&w.list, &MeshPlan::simple(0.1)).unwrap();
        (0, g.mesh.n_cells())
    });
    let g = femcore::meshing::geomdrive3::sphere_cavity(0.05, 0.004).unwrap();
    let n = g.mesh.n_cells();
    let ones = vec![1.0; n];
    time("Maxwell indicators, sphere 0.004", || {
        let opt = MaxwellOptions { k0_sq: 1.0, eps_r: &ones, mu_r: &ones, skip: &|_| false, natural: &|_| false };
        let ind = maxwell_indicators(&g.mesh, &Zero, &opt).unwrap();
        (ind.total, n)
    });
}
