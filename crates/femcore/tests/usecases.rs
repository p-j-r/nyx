//! The photonic use cases: modes of a wire on 300 nm SOI and Bragg gratings in it, each with the polarization analysis
//! (Stokes parameters, Poincare sphere) that every mode solver now returns.
use femcore::usecases::bragg::*;
use femcore::usecases::soi::*;

fn modes(r: &femcore::problem::SolveResult) -> Vec<(f64, femcore::analysis::ModeRow)> {
    let a = r.analysis.as_ref().expect("the mode solver returns a polarization analysis");
    let mut v: Vec<_> = r.solutions.iter().zip(&a.rows).map(|(s, row)| (s.neff.unwrap().0, row.clone())).collect();
    v.sort_by(|x, y| y.0.total_cmp(&x.0));
    v
}

#[test]
fn a_300_nm_soi_wire_has_quasi_te_and_quasi_tm_modes_between_the_slab_limits() {
    let wire = SoiWire::default();
    let r = wire.solve().unwrap();
    let m = modes(&r);
    assert_eq!(m.len(), 2);
    let (te, tm) = (&m[0], &m[1]);
    // the 500 nm x 300 nm wire in oxide: neff(TE0) = 2.67, neff(TM0) = 2.33, below the indices of the infinite 300 nm slab (3.049, 2.616)
    let (slab_te, slab_tm) = (wire.slab_reference(SlabPol::Te).unwrap(), wire.slab_reference(SlabPol::Tm).unwrap());
    assert!(te.0 > 2.5 && te.0 < slab_te, "TE0 {} against slab {slab_te}", te.0);
    assert!(tm.0 > 2.2 && tm.0 < slab_tm, "TM0 {} against slab {slab_tm}", tm.0);
    // polarization: the first is x polarized (the x pole of the sphere), the second y polarized, both nearly fully polarized, no handedness
    let (rte, rtm) = (te.1.report.as_ref().unwrap(), tm.1.report.as_ref().unwrap());
    assert!(rte.te_fraction > 0.97 && rte.dop > 0.9 && rte.sphere[0] > 0.9, "{rte:?}");
    assert!(rtm.te_fraction < 0.06 && rtm.dop > 0.9 && rtm.sphere[0] < -0.9, "{rtm:?}");
    assert!(rte.sphere[2].abs() < 1e-6 && rtm.sphere[2].abs() < 1e-6, "lossless modes of a symmetric wire are not circular");
    assert_eq!(rte.class.label(), "quasi-TE (x)");
    assert_eq!(rtm.class.label(), "quasi-TM (y)");
    // the longitudinal field is a sizeable part of the energy at this index contrast
    assert!(rte.long_fraction > 0.05 && rte.long_fraction < 0.4);
    // the result carries the summary and the table
    assert!(r.scalars.iter().any(|(n, _)| n == "polarization 1"));
    assert!(r.log.iter().any(|l| l.contains("quasi-TE (x)")));
    // exports
    let a = r.analysis.as_ref().unwrap();
    assert_eq!(a.csv().lines().count(), 3);
    assert!(a.sphere_svg(300.0).starts_with("<svg"));
}

#[test]
fn the_te_mode_of_a_wide_wire_tends_to_the_te_mode_of_the_slab() {
    // the distance to the slab value falls as the wire widens (the field leaks less and less to the sides); the shift sits just above the
    // slab index so that the solver returns the fundamental mode and not one of the many lateral orders in between
    let err = |w: f64| -> f64 {
        let slab = SoiWire::default().slab_reference(SlabPol::Te).unwrap();
        let wire = SoiWire { width: w, n_modes: 1, n_target: slab + 0.02, ..Default::default() };
        let m = modes(&wire.solve().unwrap());
        (slab - m[0].0).abs()
    };
    let (narrow, wide) = (err(1.0e-6), err(2.0e-6));
    assert!(wide < 0.5 * narrow, "{narrow} then {wide}");
    assert!(wide < 0.06, "TE0 of a 2 um wire within 0.06 of the slab: {wide}");
}

#[test]
fn the_te_mode_of_the_wire_is_converged_under_mesh_refinement() {
    let neff = |h: f64| modes(&SoiWire { h_core: h, h_near: 2.0 * h, ..Default::default() }.solve().unwrap())[0].0;
    let (a, b) = (neff(50e-9), neff(30e-9));
    assert!((a - b).abs() < 1.5e-3, "{a} against {b}");
}

#[test]
fn a_sloped_wall_on_one_side_makes_the_modes_hybrid_and_tilts_their_polarization() {
    // a 300 nm square core in oxide has a nearly degenerate TE/TM pair; a mirror-symmetric change keeps the modes x and y polarized,
    // a sloped left wall mixes them into two orthogonal linear states tilted by about 45 degrees
    let sym = SoiWire { width: 300e-9, ..Default::default() }.solve().unwrap();
    // the pair is degenerate (TE and TM of a square core in uniform oxide), so the eigensolver returns some mixture of the two and the
    // individual points on the sphere mean nothing; the analysis finds the group and resolves it into the x and y polarized states
    let a = sym.analysis.as_ref().unwrap();
    assert_eq!(a.groups.len(), 1, "{}", a.table());
    assert!(a.groups[0].invariant.as_ref().unwrap().dop < 0.2, "the group Stokes vector of a degenerate pair is unpolarized");
    let resolved = a.groups[0].resolved.as_ref().expect("a pair is resolved");
    let (r0, r1) = (resolved.reports[0].as_ref().unwrap(), resolved.reports[1].as_ref().unwrap());
    assert!(r0.sphere[0].abs() > 0.8 && r1.sphere[0].abs() > 0.8 && r0.sphere[0] * r1.sphere[0] < 0.0, "{:?} {:?}", r0.sphere, r1.sphere);
    let hyb = SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), ..Default::default() }.solve().unwrap();
    let m = modes(&hyb);
    let (a, b) = (m[0].1.report.as_ref().unwrap(), m[1].1.report.as_ref().unwrap());
    assert!(a.sphere[1].abs() > 0.3 && b.sphere[1].abs() > 0.3, "tilted: {:?} {:?}", a.sphere, b.sphere);
    // the two states are (nearly) antipodal on the sphere
    assert!(a.sphere[1] * b.sphere[1] < 0.0 && a.sphere[0] * b.sphere[0] < 0.0, "{:?} {:?}", a.sphere, b.sphere);
    let angle = femcore::polarization::sphere_distance(a.sphere, b.sphere);
    assert!(angle > 0.85 * std::f64::consts::PI, "angle on the sphere {angle}");
    // mixed, so each has a TE share between those of the pure modes
    assert!(a.te_fraction > 0.1 && a.te_fraction < 0.9 && b.te_fraction > 0.1 && b.te_fraction < 0.9);
}

/// Quick Bragg cells: first-order elements on a coarse mesh (some seconds each); the frequencies are only good to a few percent, the
/// comparisons between gratings are what is tested.
fn quick(kind: GratingKind) -> BraggGrating {
    BraggGrating { kind, order: 1, ..BraggGrating::coarse() }
}

#[test]
#[ignore = "slow (about 1.5 minutes): three 3D Floquet cells; run with --ignored"]
fn a_uniform_wire_has_no_gap_a_symmetric_grating_has_and_a_one_sided_one_has_less() {
    let solve = |kind| {
        let g = quick(kind);
        let r = g.solve().unwrap();
        let b = g.band_edges_of(&r).unwrap();
        (r, b)
    };
    let (ru, bu) = solve(GratingKind::Uniform);
    let (rs, bs) = solve(GratingKind::Symmetric);
    let (ra, ba) = solve(GratingKind::Asymmetric);
    // the unperturbed wire has forward and backward waves at the same frequency: a degenerate pair, one group in the analysis
    assert!(bu.gap_fraction < 1e-5, "{bu:?}");
    assert!(ru.analysis.as_ref().unwrap().groups.iter().any(|g| g.members.len() == 2));
    // the grating splits them; teeth on one side only give a smaller, but not tiny, split (about half as much perturbation)
    assert!(bs.gap_fraction > 0.002 && ba.gap_fraction > 0.001, "{bs:?} {ba:?}");
    let ratio = ba.gap_fraction / bs.gap_fraction;
    assert!(ratio > 0.3 && ratio < 0.8, "one-sided / two-sided gap = {ratio}");
    // the Bragg wavelength is near the design wavelength, the stop band is nanometres wide, and the coupling strength from the gap is
    // below the coupled-mode value that assumes long sections
    assert!((bs.bragg_wavelength / 1550e-9 - 1.0).abs() < 0.05, "{}", bs.bragg_wavelength);
    assert!(bs.kappa > 0.0 && bs.kappa < bs.kappa_estimate, "kappa {} against bound {}", bs.kappa, bs.kappa_estimate);
    assert!(bs.stopband > 1e-9 && bs.stopband < 40e-9, "{}", bs.stopband);
    // the band-edge modes are quasi-TE standing waves, with the Poincare points of a quasi-TE mode; the next two are quasi-TM
    for r in [&rs, &ra] {
        let a = r.analysis.as_ref().unwrap();
        for row in a.rows.iter().take(2) {
            let m = row.report.as_ref().unwrap();
            assert!(m.te_fraction > 0.9 && m.sphere[0] > 0.8, "{m:?}");
        }
        assert!(a.rows[2].report.as_ref().unwrap().sphere[0] < -0.8);
    }
}
