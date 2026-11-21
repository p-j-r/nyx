//! AF-6 tests for conformal finite metal `thickness` (`[z_s, z_s + t_m]`) and substrate `etch_depth`
//! trenches (`[z_s - t_e, z_s]`) on `Op::Sheet` in `geom`.

use geom::complex::{EdgeGeom, FaceGeom};
use geom::eval::Tier;
use geom::expr::Expr;
use geom::model::Model;
use geom::select::{EntityKind, Selector};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-11 * a.abs().max(b.abs()).max(1.0)
}

#[test]
fn zero_thickness_and_etch_preserves_baseline_identically() {
    let build = |use_thick_call: bool| {
        let mut m = Model::new_3d();
        let sub = m.box_("Si", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
        let vac = m.box_("VacTop", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));
        let p0 = m.rect("p0", (-60e-6, -40e-6), (-10e-6, 40e-6));
        let p1 = m.rect("p1", (10e-6, -40e-6), (60e-6, 40e-6));
        if use_thick_call {
            m.sheet_thick("c0", p0, 0.0, 0.0, 0.0);
            m.sheet_thick("c1", p1, 0.0, 0.0, 0.0);
        } else {
            m.sheet("c0", p0, 0.0);
            m.sheet("c1", p1, 0.0);
        }
        m.overlay(&[("Si", sub), ("VacTop", vac)]);
        m.eval(Tier::Complex).unwrap()
    };

    let ev_base = build(false);
    let ev_zero = build(true);

    assert_eq!(ev_base.list, ev_zero.list);
    assert_eq!(ev_base.topology_signature(), ev_zero.topology_signature());
    assert_eq!(ev_base.layer_stack().unwrap().levels(), ev_zero.layer_stack().unwrap().levels());
    assert_eq!(ev_base.faces().len(), ev_zero.faces().len());
    assert_eq!(ev_base.edges3().len(), ev_zero.edges3().len());

    let lm = ev_base.label_map().unwrap();
    assert!(lm.sheet_thickness_m.is_empty());
    assert!(lm.sheet_etch_depth_m.is_empty());
}

#[test]
fn conformal_thick_metal_and_trench_topology_volumes_and_selectors() {
    let mut m = Model::new_3d();
    let die_w = 400e-6;
    let die_h = 400e-6;
    let h_sub = 350e-6;
    let h_vac = 350e-6;
    let t_m = 150e-9;
    let t_e = 80e-9;

    let sub = m.box_("Si", (-200e-6, -200e-6, -h_sub), (200e-6, 200e-6, 0.0));
    let vac = m.box_("VacTop", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, h_vac));

    // c0: 50 um x 80 um = 4.0e-9 m^2, perimeter = 260 um
    let p0 = m.rect("p0", (-60e-6, -40e-6), (-10e-6, 40e-6));
    // c1: 50 um x 80 um = 4.0e-9 m^2, perimeter = 260 um
    let p1 = m.rect("p1", (10e-6, -40e-6), (60e-6, 40e-6));
    // Add an upper zero-thickness bridge at z = 25 um overlapping the top half of c0 (y in [0, 40e-6])
    // so c0 is split into 2 overlay cells and we verify internal seam exclusion on thick sheets!
    let p_bridge = m.rect("p_bridge", (-60e-6, 0.0), (-10e-6, 40e-6));

    m.sheet_thick("c0", p0, 0.0, t_m, t_e);
    m.sheet_thick("c1", p1, 0.0, t_m, t_e);
    m.sheet("bridge", p_bridge, 25e-6);
    m.overlay(&[("Si", sub), ("VacTop", vac)]);

    let ev = m.eval(Tier::Complex).unwrap();
    let cx = ev.complex();
    cx.check().unwrap();

    let st = ev.layer_stack().unwrap();
    let levels = st.levels();
    assert_eq!(levels.len(), 6);
    assert!(close(levels[0], -h_sub));
    assert!(close(levels[1], -t_e));
    assert!(close(levels[2], 0.0));
    assert!(close(levels[3], t_m));
    assert!(close(levels[4], 25e-6));
    assert!(close(levels[5], h_vac));

    // Check explicit is_conductor_cavity tagging in [0, t_m] and VacTop inheritance in trench [-t_e, 0]
    let si_id = 0u32;
    let vac_id = 1u32;
    let trench_layer = st.layers.iter().find(|l| (l.z0 + t_e).abs() < 1e-15 && l.z1.abs() < 1e-15).unwrap();
    let metal_layer = st.layers.iter().find(|l| l.z0.abs() < 1e-15 && (l.z1 - t_m).abs() < 1e-15).unwrap();

    let mut masked_cells = 0usize;
    let mut unmasked_cells = 0usize;
    for c in 0..st.overlay_region_count() {
        if metal_layer.is_conductor_cavity[c] {
            masked_cells += 1;
            assert_eq!(metal_layer.region_of[c], None, "conductor cavity cell must be None in metal layer");
            assert_eq!(trench_layer.region_of[c], Some(si_id), "masked cell must retain Si in trench layer");
            assert!(!trench_layer.is_conductor_cavity[c]);
        } else {
            unmasked_cells += 1;
            assert_eq!(metal_layer.region_of[c], Some(vac_id));
            assert_eq!(
                trench_layer.region_of[c],
                Some(vac_id),
                "unmasked gap cell in [-t_e, 0] must inherit VacTop dielectric"
            );
        }
    }
    assert!(masked_cells >= 3 && unmasked_cells >= 1);

    // Check exact 3D region volumes:
    let a_die = die_w * die_h;
    let a_c0 = 50e-6 * 80e-6;
    let a_c1 = 50e-6 * 80e-6;
    let a_masked = a_c0 + a_c1;
    let a_unmasked = a_die - a_masked;
    let perim_c0 = 2.0 * (50e-6 + 80e-6);
    let perim_c1 = 2.0 * (50e-6 + 80e-6);

    let want_v_si = a_die * h_sub - a_unmasked * t_e;
    let want_v_vac = a_die * h_vac + a_unmasked * t_e - a_masked * t_m;
    assert!(close(ev.volume(si_id), want_v_si), "V_Si = {:e} vs {:e}", ev.volume(si_id), want_v_si);
    assert!(close(ev.volume(vac_id), want_v_vac), "V_Vac = {:e} vs {:e}", ev.volume(vac_id), want_v_vac);

    // Check total surface area of c0, c1 (bottom cap + top cap + vertical sidewalls) and Si|VacTop (trench floor + trench sidewalls)
    let area_of = |name: &str| -> f64 { ev.faces().iter().filter(|f| f.name == name).map(|f| f.area).sum() };
    let want_a_c0 = 2.0 * a_c0 + perim_c0 * t_m;
    let want_a_c1 = 2.0 * a_c1 + perim_c1 * t_m;
    let want_a_iface = a_unmasked + (perim_c0 + perim_c1) * t_e;
    assert!(close(area_of("c0"), want_a_c0), "area(c0) = {:e} vs {:e}", area_of("c0"), want_a_c0);
    assert!(close(area_of("c1"), want_a_c1), "area(c1) = {:e} vs {:e}", area_of("c1"), want_a_c1);
    assert!(
        close(area_of("Si|VacTop"), want_a_iface),
        "area(Si|VacTop) = {:e} vs {:e}",
        area_of("Si|VacTop"),
        want_a_iface
    );

    // Check both caps and lateral sidewalls exist for c0
    let c0_caps = ev.faces().iter().filter(|f| f.name == "c0" && matches!(f.geom, FaceGeom::Cap { .. })).count();
    let c0_lats = ev.faces().iter().filter(|f| f.name == "c0" && matches!(f.geom, FaceGeom::Lateral { .. })).count();
    assert!(c0_caps >= 4, "c0 split by upper_aux should have bottom and top cap fragments");
    assert!(c0_lats >= 4, "c0 should have vertical lateral sidewalls across [0, t_m]");

    // 1. Select::Outer must exclude interior conductor cavities (c0 and c1)!
    let outer_set = ev.resolve(&Selector::Outer).unwrap().set;
    assert_eq!(outer_set.kind, EntityKind::Faces);
    for &fid in &outer_set.ids {
        let f = &ev.faces()[fid as usize];
        assert!(
            f.name != "c0" && f.name != "c1" && !f.is_conductor_cavity(),
            "Selector::Outer must not include conductor cavity face '{}'",
            f.name
        );
    }
    let outer_area: f64 = outer_set.ids.iter().map(|&fid| ev.faces()[fid as usize].area).sum();
    let want_outer_area = 2.0 * a_die + 2.0 * (die_w + die_h) * (h_sub + h_vac);
    assert!(close(outer_area, want_outer_area), "outer_area = {:e} vs {:e}", outer_area, want_outer_area);

    // 2. "<sheet>:edge" must select the upper (z = t_m) and lower (z = 0) corner rim horizontal edges uniquely,
    // excluding internal fragment seams and vertical corners!
    let c0_edges = ev.resolve(&Selector::parse("c0:edge").unwrap()).unwrap().set;
    assert_eq!(c0_edges.kind, EntityKind::Edges);
    let mut len_z0 = 0.0;
    let mut len_ztm = 0.0;
    for &eid in &c0_edges.ids {
        let e = &ev.edges3()[eid as usize];
        match e.geom {
            EdgeGeom::Horizontal { z, .. } => {
                if z.abs() < 1e-12 {
                    len_z0 += e.length;
                } else if (z - t_m).abs() < 1e-12 {
                    len_ztm += e.length;
                } else {
                    panic!("unexpected z = {z:e} on c0:edge");
                }
            }
            _ => panic!("c0:edge must only select horizontal rim edges"),
        }
    }
    assert!(close(len_z0, perim_c0), "lower rim length = {len_z0:e} vs {perim_c0:e}");
    assert!(close(len_ztm, perim_c0), "upper rim length = {len_ztm:e} vs {perim_c0:e}");

    // 3. LabelMap and Complex carry sheet_thickness_m and sheet_etch_depth_m per sheet boundary label
    assert_eq!(cx.sheet_thickness("c0"), Some(t_m));
    assert_eq!(cx.sheet_etch_depth("c0"), Some(t_e));
    let lm = ev.label_map().unwrap();
    let c0_labs = lm.boundary("c0");
    let c1_labs = lm.boundary("c1");
    assert!(!c0_labs.is_empty() && !c1_labs.is_empty());
    for &id in c0_labs.iter().chain(c1_labs.iter()) {
        assert_eq!(lm.sheet_thickness_m.get(&id), Some(&t_m));
        assert_eq!(lm.sheet_etch_depth_m.get(&id), Some(&t_e));
    }
}

#[test]
fn nanometer_1nm_film_succeeds_and_invalid_conflicts_are_rejected() {
    // 1. 1 nm explicit film on a 10 mm die succeeds without triggering GRAZE_K * tol
    let mut m = Model::new_3d();
    let sub = m.box_("Si", (-5e-3, -5e-3, -500e-6), (5e-3, 5e-3, 0.0));
    let vac = m.box_("VacTop", (-5e-3, -5e-3, 0.0), (5e-3, 5e-3, 500e-6));
    let p = m.rect("p", (-100e-6, -100e-6), (100e-6, 100e-6));
    m.sheet_thick("c_1nm", p, 0.0, 1.0e-9, 1.0e-9);
    m.overlay(&[("Si", sub), ("VacTop", vac)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert_eq!(ev.complex().sheet_thickness("c_1nm"), Some(1.0e-9));

    // 2. Negative thickness or etch_depth is rejected
    let mut m_neg = Model::new_3d();
    let s1 = m_neg.box_("Si", (-100e-6, -100e-6, -50e-6), (100e-6, 100e-6, 0.0));
    let v1 = m_neg.box_("VacTop", (-100e-6, -100e-6, 0.0), (100e-6, 100e-6, 50e-6));
    let p1 = m_neg.rect("p", (-20e-6, -20e-6), (20e-6, 20e-6));
    m_neg.sheet_thick("c_bad", p1, 0.0, -10e-9, 0.0);
    m_neg.overlay(&[("Si", s1), ("VacTop", v1)]);
    assert!(m_neg.eval(Tier::Complex).is_err());

    // 3. Same-plane conflicting etch_depth is rejected
    let mut m_conf_etch = Model::new_3d();
    let s2 = m_conf_etch.box_("Si", (-100e-6, -100e-6, -50e-6), (100e-6, 100e-6, 0.0));
    let v2 = m_conf_etch.box_("VacTop", (-100e-6, -100e-6, 0.0), (100e-6, 100e-6, 50e-6));
    let pa = m_conf_etch.rect("pa", (-50e-6, -20e-6), (-10e-6, 20e-6));
    let pb = m_conf_etch.rect("pb", (10e-6, -20e-6), (50e-6, 20e-6));
    m_conf_etch.sheet_thick("ca", pa, 0.0, 150e-9, 80e-9);
    m_conf_etch.sheet_thick("cb", pb, 0.0, 150e-9, 50e-9);
    m_conf_etch.overlay(&[("Si", s2), ("VacTop", v2)]);
    assert!(m_conf_etch.eval(Tier::Complex).is_err());

    // 4. Overlapping sheets on the same plane with differing thickness is rejected
    let mut m_conf_thick = Model::new_3d();
    let s3 = m_conf_thick.box_("Si", (-100e-6, -100e-6, -50e-6), (100e-6, 100e-6, 0.0));
    let v3 = m_conf_thick.box_("VacTop", (-100e-6, -100e-6, 0.0), (100e-6, 100e-6, 50e-6));
    let pa3 = m_conf_thick.rect("pa", (-30e-6, -20e-6), (10e-6, 20e-6));
    let pb3 = m_conf_thick.rect("pb", (-10e-6, -20e-6), (30e-6, 20e-6));
    m_conf_thick.sheet_thick("ca", pa3, 0.0, 150e-9, 0.0);
    m_conf_thick.sheet_thick("cb", pb3, 0.0, 200e-9, 0.0);
    m_conf_thick.overlay(&[("Si", s3), ("VacTop", v3)]);
    assert!(m_conf_thick.eval(Tier::Complex).is_err());

    // 5. Sheet touching outer z_max or z_min is rejected
    let mut m_outer = Model::new_3d();
    let s4 = m_outer.box_("Si", (-100e-6, -100e-6, -50e-6), (100e-6, 100e-6, 0.0));
    let v4 = m_outer.box_("VacTop", (-100e-6, -100e-6, 0.0), (100e-6, 100e-6, 50e-6));
    let p4 = m_outer.rect("p", (-20e-6, -20e-6), (20e-6, 20e-6));
    m_outer.sheet_thick("c_top", p4, Expr::Num(50e-6), 100e-9, 0.0);
    m_outer.overlay(&[("Si", s4), ("VacTop", v4)]);
    assert!(m_outer.eval(Tier::Complex).is_err());
}
