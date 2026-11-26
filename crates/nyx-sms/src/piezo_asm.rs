// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Coupled piezoelectric `(u, phi)` matrix assembly and constraint reduction across both
//! [`SimplexMesh`] and [`MixedDiscretisation`] (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`), with
//! connected-component potential gauge pinning, `Plane::Stress` electromechanical condensation,
//! static `InfiniteElement` (Rational & Exponential+Robin), curvilinear `PmlJacobian`, Floquet
//! periodic pairs, Kelvin-Voigt viscosity + `LowReflecting` damping, and thermal-piezo loads.

use crate::floquet::mixed_nodal_relations;
use crate::lrb::{damping_matrix, damping_matrix_mixed};
use crate::piezo_mat::{cell_constitutive_blocks, resolve_piezo_mats_mixed, resolve_piezo_mats_simplex, PiezoCellMat};
use crate::pml::{
    b_matrix_c, outer_clamped_dofs_mixed, outer_clamped_dofs_simplex, simplex_barycentric, ResolvedPmlDomains,
};
use crate::problem::{BcKind, ConductorRole, Plane, Problem, Study, ThermalSource};
use crate::solid::{
    add_boundary_load_mixed, assemble_exponential_ie_robin, b_matrix, collect_exponential_ie_faces,
    eval_raw_face_jacobian, mixed_cell_tensor, mixed_check_boundary_labels, mixed_untouched_dofs, raw_face_dofs,
    thermal_strain, ResolvedIeDomains, Thermal,
};
use femcore::fem::{
    add_boundary_load, assemble_mass, boundary_dofs, check_boundary_labels, geometry, quad_rule, shape,
    shape_gradients, untouched_dofs, RefElement,
};
use femcore::fem_mixed::{
    assemble_mass_mixed, eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, lagrange_shape_into,
    quad_rule_mixed, IsoparametricCell, MixedLagrangeDofMap,
};
use femcore::linalg::csr::{c, CCoo, CCsr, Coo, RCsr, C};
use femcore::linalg::par;
use femcore::material::M3;
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::mesh::SimplexMesh;
use femcore::periodic::{congruence_fast, congruence_real, match_vertices, nodal_relations, phase_factor, Reduction};
use femcore::props::TensorKind;
use femcore::section::{plane_tensor, wafer_axes, wafer_point};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// One electrical electrode group (fixed potential, ground, floating boundary, or conductor domain)
/// tracked for Stage 3 free-charge / capacitance reporting and potential visualization.
#[derive(Clone, Debug)]
pub struct ElectrodeInfo {
    pub label: String,
    pub bc_index: usize,
    pub dofs: Vec<usize>,
    pub fixed_voltage: Option<f64>,
    pub prescribed_charge: Option<f64>,
    pub is_floating: bool,
    pub root_dof: Option<usize>,
}

/// Driven multi-port terminal electrode (`BcKind::Terminal` or `ConductorRole::Terminal`).
#[derive(Clone, Debug)]
pub struct PiezoTerminal {
    pub name: String,
    pub z0_ohm: f64,
    pub dofs: Vec<usize>,
}

/// Complete reduced constraints, block sizes, and electrode metadata for a coupled solve.
#[derive(Clone, Debug)]
pub struct PiezoConstraints {
    pub d: usize,
    pub n1: usize,
    pub n_mech: usize,
    pub n_tot: usize,
    pub red_all: Reduction,
    pub red_u: Reduction,
    pub g_all: Vec<C>,
    pub n_u_red: usize,
    pub n_phi_red: usize,
    pub fixed_phi_vals: Vec<(usize, f64)>,
    pub gauge_pins: Vec<usize>,
    pub electrodes: Vec<ElectrodeInfo>,
    pub conductor_domains: Vec<i32>,
    pub terminals: Vec<PiezoTerminal>,
    pub terminal_lifts: Vec<Vec<C>>,
}

/// Assembled static coupled piezoelectric problem (real arithmetic).
pub struct PiezoStaticAssembly {
    pub constraints: PiezoConstraints,
    pub mats: Vec<PiezoCellMat>,
    pub by_dom: BTreeMap<i32, PiezoCellMat>,
    pub k_full: RCsr,
    pub k_red: RCsr,
    pub f_full: Vec<f64>,
    pub f_mech_elec_load_only: Vec<f64>,
    pub f_th_only: Vec<f64>,
    pub rhs_red: Vec<f64>,
    pub thermal: Option<Thermal>,
    pub(crate) ie_domains: ResolvedIeDomains,
    pub diagnostics: Vec<String>,
}

/// Assembled eigenfrequency coupled piezoelectric problem (real lossless or complex/damped).
pub struct PiezoEigenAssembly {
    pub constraints: PiezoConstraints,
    pub mats: Vec<PiezoCellMat>,
    pub by_dom: BTreeMap<i32, PiezoCellMat>,
    pub k_full_real: Option<RCsr>,
    pub k_red_real: Option<RCsr>,
    pub m_uu_red_real: Option<RCsr>,
    pub k_full_c: Option<CCsr>,
    pub k_red_c: Option<CCsr>,
    pub m_uu_red_c: Option<CCsr>,
    pub c_uu_red_c: Option<CCsr>,
    pub k_pp_red_real: Option<RCsr>,
    pub k_pp_red_c: Option<CCsr>,
    pub k_pu_red_real: Option<RCsr>,
    pub k_pu_red_c: Option<CCsr>,
    pub m1_phys_all: Option<RCsr>,
    pub m1_phys_pml: Option<RCsr>,
    pub(crate) resolved_pml: ResolvedPmlDomains,
    pub is_complex: bool,
    pub is_hermitian: bool,
    pub has_damping: bool,
    pub has_pml: bool,
    pub has_interior: bool,
    pub n_modes: usize,
    pub f_search: f64,
    pub n_pairs: usize,
    pub diagnostics: Vec<String>,
}

fn check_labels_simplex(mesh: &SimplexMesh, spec: &Problem) -> Result<(), String> {
    check_boundary_labels(
        mesh,
        spec.bcs.iter().filter(|b| !b.kind.is_domain_condition()).flat_map(|b| b.labels.iter().copied()),
    )?;
    let doms = mesh.domain_labels();
    let bad: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| b.kind.is_domain_condition())
        .flat_map(|b| b.labels.iter().copied())
        .filter(|l| !doms.contains(l))
        .collect();
    if !bad.is_empty() {
        return Err(format!("domain number(s) {bad:?} do not exist in this mesh (it has {doms:?})"));
    }
    Ok(())
}

fn collect_conductor_domain_dofs_simplex(mesh: &SimplexMesh, order: usize, labels: &[i32]) -> Vec<usize> {
    let mut set = BTreeSet::new();
    let mut scratch = Vec::new();
    for ci in 0..mesh.n_cells() {
        if labels.contains(&mesh.cell_labels[ci]) {
            mesh.cell_dofs(ci, order, &mut scratch);
            set.extend(scratch.iter().copied());
        }
    }
    set.into_iter().collect()
}

fn collect_conductor_domain_dofs_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    labels: &[i32],
) -> Vec<usize> {
    let mut set = BTreeSet::new();
    for ci in 0..mixed.n_cells() {
        if labels.contains(&mixed.cell_labels[ci]) {
            set.extend(dof_map.cell_dofs(ci).iter().copied());
        }
    }
    set.into_iter().collect()
}

#[allow(clippy::too_many_arguments)]
fn build_piezo_constraints(
    dim: usize,
    n1: usize,
    spec: &Problem,
    cell_labels: &[i32],
    get_cell_dofs: &dyn Fn(usize, &mut Vec<usize>),
    get_boundary_dofs: &dyn Fn(&[i32]) -> Vec<usize>,
    get_conductor_dofs: &dyn Fn(&[i32]) -> Vec<usize>,
    mech_untouched: Vec<usize>,
    outer_clamped_nodes: Vec<usize>,
    ie_robin_nodes: Vec<usize>,
    periodic_rels: Vec<(usize, usize, C)>,
) -> Result<PiezoConstraints, String> {
    let n_mech = dim * n1;
    let n_tot = (dim + 1) * n1;
    let conductor_domains = spec.conductor_domains();

    let mut relations: Vec<(usize, usize, C)> = Vec::new();
    let mut fixed_vals: Vec<(usize, C)> = Vec::new();
    let mut mech_rels: Vec<(usize, usize, C)> = Vec::new();
    let mut mech_fixed: Vec<(usize, C)> = Vec::new();

    for dof in mech_untouched {
        fixed_vals.push((dof, c(0.0, 0.0)));
        mech_fixed.push((dof, c(0.0, 0.0)));
    }
    for &node in &outer_clamped_nodes {
        for comp in 0..dim {
            let gd = comp * n1 + node;
            fixed_vals.push((gd, c(0.0, 0.0)));
            mech_fixed.push((gd, c(0.0, 0.0)));
        }
        fixed_vals.push((n_mech + node, c(0.0, 0.0)));
    }

    let axes = wafer_axes(dim);
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, values } = &bc.kind {
            let bd = get_boundary_dofs(&bc.labels);
            for (comp, &w) in axes.iter().enumerate() {
                if fixed[w] {
                    let val = c(values[w], 0.0);
                    for &dof in &bd {
                        let gd = comp * n1 + dof;
                        fixed_vals.push((gd, val));
                        mech_fixed.push((gd, val));
                    }
                }
            }
        }
    }

    // Separate Floquet periodic relations between mechanical and electrical blocks
    for &(dest, src, fac) in &periodic_rels {
        relations.push((dest, src, fac));
        if dest < n_mech && src < n_mech {
            mech_rels.push((dest, src, fac));
        }
    }

    // Electrical connectivity graph on scalar nodes `0..n1` for gauge BFS
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n1];
    let mut elec_active = vec![false; n1];
    let mut elec_anchored = vec![false; n1];
    let mut cell_scratch = Vec::new();
    for (ci, &dom) in cell_labels.iter().enumerate() {
        if conductor_domains.contains(&dom) {
            continue;
        }
        get_cell_dofs(ci, &mut cell_scratch);
        if let Some(&anchor) = cell_scratch.first() {
            elec_active[anchor] = true;
            for &other in &cell_scratch[1..] {
                elec_active[other] = true;
                adj[anchor].push(other);
                adj[other].push(anchor);
            }
        }
    }
    let has_dielectric_cell = elec_active.clone();
    for &node in &outer_clamped_nodes {
        elec_anchored[node] = true;
    }
    for &node in &ie_robin_nodes {
        elec_anchored[node] = true;
    }
    for &(dest, src, fac) in &periodic_rels {
        if dest >= n_mech && src >= n_mech {
            let ud = dest - n_mech;
            let us = src - n_mech;
            adj[ud].push(us);
            adj[us].push(ud);
            if (fac - c(1.0, 0.0)).norm() > 1e-9 {
                elec_anchored[ud] = true;
                elec_anchored[us] = true;
            }
        }
    }

    let mut electrodes = Vec::new();
    let mut floating_roots_check: Vec<(usize, usize, f64, String)> = Vec::new();
    let mut fixed_phi_map: BTreeMap<usize, f64> = BTreeMap::new();
    let mut terminals: Vec<PiezoTerminal> = Vec::new();
    let mut non_term_dirichlet: BTreeSet<usize> = BTreeSet::new();

    for (bc_index, bc) in spec.bcs.iter().enumerate() {
        match &bc.kind {
            BcKind::Ground => {
                let bd = get_boundary_dofs(&bc.labels);
                if !bd.is_empty() && bd.iter().all(|&dof| !has_dielectric_cell[dof]) {
                    return Err(format!(
                        "electrical boundary condition 'ground' on boundaries {:?} has no adjacent active dielectric cells",
                        bc.labels
                    ));
                }
                for &dof in &bd {
                    elec_anchored[dof] = true;
                    non_term_dirichlet.insert(dof);
                    fixed_vals.push((n_mech + dof, c(0.0, 0.0)));
                    fixed_phi_map.entry(dof).or_insert(0.0);
                }
                electrodes.push(ElectrodeInfo {
                    label: format!("ground_{bc_index}"),
                    bc_index,
                    dofs: bd,
                    fixed_voltage: Some(0.0),
                    prescribed_charge: None,
                    is_floating: false,
                    root_dof: None,
                });
            }
            BcKind::Potential(v) => {
                let bd = get_boundary_dofs(&bc.labels);
                if !bd.is_empty() && bd.iter().all(|&dof| !has_dielectric_cell[dof]) {
                    return Err(format!(
                        "electrical boundary condition 'potential' on boundaries {:?} has no adjacent active dielectric cells",
                        bc.labels
                    ));
                }
                for &dof in &bd {
                    elec_anchored[dof] = true;
                    non_term_dirichlet.insert(dof);
                    fixed_vals.push((n_mech + dof, c(*v, 0.0)));
                    fixed_phi_map.entry(dof).or_insert(*v);
                }
                electrodes.push(ElectrodeInfo {
                    label: format!("potential_{bc_index}"),
                    bc_index,
                    dofs: bd,
                    fixed_voltage: Some(*v),
                    prescribed_charge: None,
                    is_floating: false,
                    root_dof: None,
                });
            }
            BcKind::SurfaceCharge(_) => {
                let bd = get_boundary_dofs(&bc.labels);
                if !bd.is_empty() && bd.iter().all(|&dof| !has_dielectric_cell[dof]) {
                    return Err(format!(
                        "electrical boundary condition 'surface_charge' on boundaries {:?} has no adjacent active dielectric cells",
                        bc.labels
                    ));
                }
            }
            BcKind::FloatingPotential { charge } => {
                let bd = get_boundary_dofs(&bc.labels);
                if !bd.is_empty() && bd.iter().all(|&dof| !has_dielectric_cell[dof]) {
                    return Err(format!(
                        "electrical boundary condition 'floating_potential' on boundaries {:?} has no adjacent active dielectric cells",
                        bc.labels
                    ));
                }
                if let Some(&root) = bd.first() {
                    elec_active[root] = true;
                    for &dof in &bd[1..] {
                        elec_active[dof] = true;
                        adj[root].push(dof);
                        adj[dof].push(root);
                        relations.push((n_mech + dof, n_mech + root, c(1.0, 0.0)));
                    }
                    floating_roots_check.push((bc_index, root, *charge, format!("floating_potential_{bc_index}")));
                    electrodes.push(ElectrodeInfo {
                        label: format!("floating_potential_{bc_index}"),
                        bc_index,
                        dofs: bd,
                        fixed_voltage: None,
                        prescribed_charge: Some(*charge),
                        is_floating: true,
                        root_dof: Some(root),
                    });
                }
            }
            BcKind::Terminal { name, z0_ohm } => {
                let bd = get_boundary_dofs(&bc.labels);
                if bd.is_empty() {
                    return Err(format!(
                        "electrical terminal '{name}' on boundaries {:?} matched zero nodes",
                        bc.labels
                    ));
                }
                if bd.iter().all(|&dof| !has_dielectric_cell[dof]) {
                    return Err(format!(
                        "electrical terminal '{name}' on boundaries {:?} has no adjacent active dielectric cells",
                        bc.labels
                    ));
                }
                for &dof in &bd {
                    elec_anchored[dof] = true;
                    fixed_phi_map.entry(dof).or_insert(0.0);
                }
                terminals.push(PiezoTerminal { name: name.clone(), z0_ohm: *z0_ohm, dofs: bd });
            }
            BcKind::Conductor { role, value, charge, name, z0_ohm } => {
                let cd = get_conductor_dofs(&bc.labels);
                match role {
                    ConductorRole::Ground => {
                        for &dof in &cd {
                            elec_anchored[dof] = true;
                            non_term_dirichlet.insert(dof);
                            fixed_vals.push((n_mech + dof, c(0.0, 0.0)));
                            fixed_phi_map.entry(dof).or_insert(0.0);
                        }
                        electrodes.push(ElectrodeInfo {
                            label: format!("conductor_ground_{bc_index}"),
                            bc_index,
                            dofs: cd,
                            fixed_voltage: Some(0.0),
                            prescribed_charge: None,
                            is_floating: false,
                            root_dof: None,
                        });
                    }
                    ConductorRole::Fixed => {
                        for &dof in &cd {
                            elec_anchored[dof] = true;
                            non_term_dirichlet.insert(dof);
                            fixed_vals.push((n_mech + dof, c(*value, 0.0)));
                            fixed_phi_map.entry(dof).or_insert(*value);
                        }
                        electrodes.push(ElectrodeInfo {
                            label: format!("conductor_fixed_{bc_index}"),
                            bc_index,
                            dofs: cd,
                            fixed_voltage: Some(*value),
                            prescribed_charge: None,
                            is_floating: false,
                            root_dof: None,
                        });
                    }
                    ConductorRole::Terminal => {
                        if cd.is_empty() {
                            return Err(format!(
                                "conductor terminal '{name}' on domain(s) {:?} matched zero nodes",
                                bc.labels
                            ));
                        }
                        for &dof in &cd {
                            elec_anchored[dof] = true;
                            fixed_phi_map.entry(dof).or_insert(0.0);
                        }
                        terminals.push(PiezoTerminal { name: name.clone(), z0_ohm: *z0_ohm, dofs: cd });
                    }
                    ConductorRole::Floating => {
                        if let Some(&root) = cd.first() {
                            elec_active[root] = true;
                            for &dof in &cd[1..] {
                                elec_active[dof] = true;
                                adj[root].push(dof);
                                adj[dof].push(root);
                                relations.push((n_mech + dof, n_mech + root, c(1.0, 0.0)));
                            }
                            floating_roots_check.push((
                                bc_index,
                                root,
                                *charge,
                                format!("conductor_floating_{bc_index}"),
                            ));
                            electrodes.push(ElectrodeInfo {
                                label: format!("conductor_floating_{bc_index}"),
                                bc_index,
                                dofs: cd,
                                fixed_voltage: None,
                                prescribed_charge: Some(*charge),
                                is_floating: true,
                                root_dof: Some(root),
                            });
                        }
                    }
                }
            }
            BcKind::Fixed { .. } | BcKind::Load(_) | BcKind::LowReflecting => {}
        }
    }

    // Terminals cannot overlap fixed/grounded potentials or each other (would short-circuit port voltage excitation).
    let mut term_owner: BTreeMap<usize, &str> = BTreeMap::new();
    for term in &terminals {
        for &dof in &term.dofs {
            if non_term_dirichlet.contains(&dof) {
                return Err(format!("terminal '{}' overlaps a fixed/grounded potential", term.name));
            }
            if let Some(prev) = term_owner.insert(dof, term.name.as_str()) {
                if prev != term.name {
                    return Err(format!("terminal '{}' overlaps terminal '{prev}'", term.name));
                }
            }
        }
    }

    // Floating electrodes cannot overlap explicit fixed/grounded potentials or terminals (checked before gauge pinning)
    let explicit_fixed_phi = fixed_phi_map.clone();
    for el in &electrodes {
        if el.is_floating && el.dofs.iter().any(|dof| explicit_fixed_phi.contains_key(dof)) {
            return Err(format!("floating electrode '{}' overlaps a fixed/grounded potential or terminal", el.label));
        }
    }

    // Connected-component electrical gauge BFS: pin 1 DOF per floating electrical component
    let mut visited = vec![false; n1];
    let mut queue = VecDeque::new();
    for i in 0..n1 {
        if elec_anchored[i] {
            visited[i] = true;
            queue.push_back(i);
        }
    }
    while let Some(u) = queue.pop_front() {
        for &v in &adj[u] {
            if !visited[v] {
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }
    let mut gauge_pins = Vec::new();
    for seed in 0..n1 {
        if !elec_active[seed] || visited[seed] {
            continue;
        }
        gauge_pins.push(seed);
        fixed_vals.push((n_mech + seed, c(0.0, 0.0)));
        fixed_phi_map.entry(seed).or_insert(0.0);
        visited[seed] = true;
        queue.push_back(seed);
        while let Some(u) = queue.pop_front() {
            for &v in &adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }
    }

    if matches!(spec.study, Study::Static) && !gauge_pins.is_empty() {
        return Err(format!(
            "static coupled piezoelectric solve has {} active electrical component(s) with no Dirichlet potential reference ('ground', 'potential', or grounded/fixed 'conductor'); add a potential reference (or use an eigenfrequency study)",
            gauge_pins.len()
        ));
    }

    // Anchor any completely unused electrical DOF (e.g., isolated vertices)
    for i in 0..n1 {
        if !elec_active[i] && !elec_anchored[i] {
            fixed_vals.push((n_mech + i, c(0.0, 0.0)));
            fixed_phi_map.entry(i).or_insert(0.0);
        }
    }

    let mut all_fixed_vals = fixed_vals.clone();
    for term in &terminals {
        for &dof in &term.dofs {
            all_fixed_vals.push((n_mech + dof, c(0.0, 0.0)));
        }
    }
    let (red_all, g_all) = Reduction::try_build_with_values(n_tot, &relations, &all_fixed_vals)?;
    let (red_u, _) = Reduction::try_build_with_values(n_mech, &mech_rels, &mech_fixed)?;

    let mut terminal_lifts = Vec::with_capacity(terminals.len());
    if !terminals.is_empty() {
        if spec.periodic.is_empty() {
            for term in &terminals {
                let mut g = vec![c(0.0, 0.0); n_tot];
                for &dof in &term.dofs {
                    g[n_mech + dof] = c(1.0, 0.0);
                }
                terminal_lifts.push(g);
            }
        } else {
            for (q_idx, _) in terminals.iter().enumerate() {
                let mut port_fixed = fixed_vals.clone();
                for (ti, term) in terminals.iter().enumerate() {
                    let val = if ti == q_idx { c(1.0, 0.0) } else { c(0.0, 0.0) };
                    for &dof in &term.dofs {
                        port_fixed.push((n_mech + dof, val));
                    }
                }
                let (_, g_port) = Reduction::try_build_with_values(n_tot, &relations, &port_fixed)?;
                terminal_lifts.push(g_port);
            }
        }
    }

    for (_, root, charge, name) in &floating_roots_check {
        if red_all.index_of[n_mech + root] == usize::MAX
            && (explicit_fixed_phi.contains_key(root) || charge.abs() > 0.0)
        {
            return Err(format!(
                "floating electrode '{name}' (charge {charge} C) is shorted to a fixed/grounded potential"
            ));
        }
    }

    let n_u_red = red_u.n_reduced();
    let n_phi_red = red_all.n_reduced().saturating_sub(n_u_red);
    debug_assert_eq!(red_all.roots[..n_u_red], red_u.roots[..], "reduction must preserve [u_red; phi_red] partition");
    debug_assert!(red_all.roots[n_u_red..].iter().all(|&r| r >= n_mech));

    let fixed_phi_vals = fixed_phi_map.into_iter().collect();

    Ok(PiezoConstraints {
        d: dim,
        n1,
        n_mech,
        n_tot,
        red_all,
        red_u,
        g_all,
        n_u_red,
        n_phi_red,
        fixed_phi_vals,
        gauge_pins,
        electrodes,
        conductor_domains,
        terminals,
        terminal_lifts,
    })
}

/// Evaluates local stiffness/coupling/permittivity blocks at one real quadrature point.
#[inline]
#[allow(clippy::too_many_arguments)]
fn accumulate_piezo_qp_real(
    mat: &PiezoCellMat,
    d: usize,
    plane: Plane,
    nloc: usize,
    grad: &[f64],
    w_eff: f64,
    b: &mut [f64],
    ke_uu: &mut [f64],
    ke_up: &mut [f64],
    ke_pp: &mut [f64],
) {
    let nd = d * nloc;
    if mat.mech_active {
        b_matrix(grad, nloc, d, b);
    }
    if d == 2 {
        let blk = cell_constitutive_blocks(mat, 2, plane);
        if mat.mech_active {
            for r in 0..3 {
                for j in 0..nd {
                    let db_rj: f64 = (0..3).map(|t| blk.d_uu[r * 3 + t] * b[t * nd + j]).sum();
                    let wb = w_eff * db_rj;
                    for i in 0..=j {
                        ke_uu[i * nd + j] += b[r * nd + i] * wb;
                    }
                }
            }
            if mat.elec_active && blk.e_pu.iter().any(|&v| v != 0.0) {
                for i in 0..nd {
                    for p in 0..2 {
                        let be_ip: f64 = (0..3).map(|s| b[s * nd + i] * blk.e_pu[p * 3 + s]).sum();
                        let wb = w_eff * be_ip;
                        for a in 0..nloc {
                            ke_up[i * nloc + a] += wb * grad[a * 2 + p];
                        }
                    }
                }
            }
        }
        if mat.elec_active {
            for a in 0..nloc {
                for b_idx in a..nloc {
                    let mut acc = 0.0;
                    for p in 0..2 {
                        for q in 0..2 {
                            acc += grad[a * 2 + p] * blk.eps_ps_re[p * 2 + q] * grad[b_idx * 2 + q];
                        }
                    }
                    ke_pp[a * nloc + b_idx] -= w_eff * acc;
                }
            }
        }
    } else {
        let ns = 6;
        if mat.mech_active {
            for r in 0..ns {
                for j in 0..nd {
                    let db_rj: f64 = (0..ns).map(|t| mat.c6[r][t] * b[t * nd + j]).sum();
                    let wb = w_eff * db_rj;
                    for i in 0..=j {
                        ke_uu[i * nd + j] += b[r * nd + i] * wb;
                    }
                }
            }
            if mat.elec_active && mat.e.iter().flatten().any(|&v| v != 0.0) {
                for i in 0..nd {
                    for p in 0..3 {
                        let be_ip: f64 = (0..ns).map(|s| b[s * nd + i] * mat.e[p][s]).sum();
                        let wb = w_eff * be_ip;
                        for a in 0..nloc {
                            ke_up[i * nloc + a] += wb * grad[a * 3 + p];
                        }
                    }
                }
            }
        }
        if mat.elec_active {
            for a in 0..nloc {
                for b_idx in a..nloc {
                    let mut acc = 0.0;
                    for p in 0..3 {
                        for q in 0..3 {
                            acc += grad[a * 3 + p] * mat.eps_re[p][q] * grad[b_idx * 3 + q];
                        }
                    }
                    ke_pp[a * nloc + b_idx] -= w_eff * acc;
                }
            }
        }
    }
}

/// Evaluates local complex stiffness/coupling/permittivity blocks at one complex quadrature point
/// (supporting curvilinear PML stretching, adiabatic stiffness loss, and complex `eps = eps' - i eps''`).
#[inline]
#[allow(clippy::too_many_arguments)]
fn accumulate_piezo_qp_complex(
    mat: &PiezoCellMat,
    d: usize,
    plane: Plane,
    nloc: usize,
    gc: &[C],
    w_k: C,
    w_coup: C,
    bc: &mut [C],
    ke_uu: &mut [C],
    ke_up: &mut [C],
    ke_pp: &mut [C],
) {
    let nd = d * nloc;
    if mat.mech_active {
        b_matrix_c(gc, nloc, d, bc);
    }
    if d == 2 {
        let blk = cell_constitutive_blocks(mat, 2, plane);
        if mat.mech_active {
            for r in 0..3 {
                for j in 0..nd {
                    let db_rj: C = (0..3).fold(c(0.0, 0.0), |acc, t| acc + bc[t * nd + j] * blk.d_uu[r * 3 + t]);
                    let wb = w_k * db_rj;
                    for i in 0..=j {
                        ke_uu[i * nd + j] += bc[r * nd + i] * wb;
                    }
                }
            }
            if mat.elec_active && blk.e_pu.iter().any(|&v| v != 0.0) {
                for i in 0..nd {
                    for p in 0..2 {
                        let be_ip: C = (0..3).fold(c(0.0, 0.0), |acc, s| acc + bc[s * nd + i] * blk.e_pu[p * 3 + s]);
                        let wb = w_coup * be_ip;
                        for a in 0..nloc {
                            ke_up[i * nloc + a] += wb * gc[a * 2 + p];
                        }
                    }
                }
            }
        }
        if mat.elec_active {
            for a in 0..nloc {
                for b_idx in a..nloc {
                    let mut acc = c(0.0, 0.0);
                    for p in 0..2 {
                        for q in 0..2 {
                            let eps_c = c(blk.eps_ps_re[p * 2 + q], -blk.eps_ps_im[p * 2 + q]);
                            acc += gc[a * 2 + p] * eps_c * gc[b_idx * 2 + q];
                        }
                    }
                    ke_pp[a * nloc + b_idx] -= w_coup * acc;
                }
            }
        }
    } else {
        let ns = 6;
        if mat.mech_active {
            for r in 0..ns {
                for j in 0..nd {
                    let db_rj: C = (0..ns).fold(c(0.0, 0.0), |acc, t| acc + bc[t * nd + j] * mat.c6[r][t]);
                    let wb = w_k * db_rj;
                    for i in 0..=j {
                        ke_uu[i * nd + j] += bc[r * nd + i] * wb;
                    }
                }
            }
            if mat.elec_active && mat.e.iter().flatten().any(|&v| v != 0.0) {
                for i in 0..nd {
                    for p in 0..3 {
                        let be_ip: C = (0..ns).fold(c(0.0, 0.0), |acc, s| acc + bc[s * nd + i] * mat.e[p][s]);
                        let wb = w_coup * be_ip;
                        for a in 0..nloc {
                            ke_up[i * nloc + a] += wb * gc[a * 3 + p];
                        }
                    }
                }
            }
        }
        if mat.elec_active {
            for a in 0..nloc {
                for b_idx in a..nloc {
                    let mut acc = c(0.0, 0.0);
                    for p in 0..3 {
                        for q in 0..3 {
                            let eps_c = c(mat.eps_re[p][q], -mat.eps_im[p][q]);
                            acc += gc[a * 3 + p] * eps_c * gc[b_idx * 3 + q];
                        }
                    }
                    ke_pp[a * nloc + b_idx] -= w_coup * acc;
                }
            }
        }
    }
}

/// Evaluates local mechanical and dielectric thermal load vectors `(fe_u, fe_phi)` at one quadrature point.
#[inline]
#[allow(clippy::too_many_arguments)]
fn accumulate_piezo_thermal_qp(
    mat: &PiezoCellMat,
    d: usize,
    plane: Plane,
    nloc: usize,
    grad: &[f64],
    alpha: &M3,
    dt: f64,
    w_vol: f64,
    b: &mut [f64],
    fe_u: &mut [f64],
    fe_phi: &mut [f64],
) {
    let nd = d * nloc;
    let eth = thermal_strain(alpha, dt);
    b_matrix(grad, nloc, d, b);
    if d == 2 {
        match plane {
            Plane::Stress => {
                let blk = cell_constitutive_blocks(mat, 2, Plane::Stress);
                let eth_in = [eth[0], eth[1], eth[3]];
                for i in 0..nd {
                    let mut fu = 0.0;
                    for s in 0..3 {
                        let sig_s: f64 = (0..3).map(|t| blk.d_uu[s * 3 + t] * eth_in[t]).sum();
                        fu += b[s * nd + i] * sig_s;
                    }
                    fe_u[i] += w_vol * fu;
                }
                // `f_phi^{th} = -int B_phi^T e_o (C_oo)^-1 (C^E S_th)_o dOmega` (only out-of-plane `yy` eliminated in 2D)
                if mat.elec_active && mat.c6[2][2].abs() > 1e-30 && (mat.e[0][2] != 0.0 || mat.e[1][2] != 0.0) {
                    let c_eth_o: f64 = (0..4).map(|k| mat.c6[2][k] * eth[k]).sum();
                    let s_o_th = c_eth_o / mat.c6[2][2];
                    for p in 0..2 {
                        let wd = -w_vol * mat.e[p][2] * s_o_th;
                        for a in 0..nloc {
                            fe_phi[a] += wd * grad[a * 2 + p];
                        }
                    }
                }
            }
            Plane::Strain => {
                for i in 0..nd {
                    let mut fu = 0.0;
                    for (s_in, &row_idx) in [0usize, 1, 3].iter().enumerate() {
                        let sig_s: f64 = (0..4).map(|k| mat.c6[row_idx][k] * eth[k]).sum();
                        fu += b[s_in * nd + i] * sig_s;
                    }
                    fe_u[i] += w_vol * fu;
                }
            }
        }
    } else {
        for i in 0..nd {
            let mut fu = 0.0;
            for s in 0..6 {
                let sig_s: f64 = (0..6).map(|t| mat.c6[s][t] * eth[t]).sum();
                fu += b[s * nd + i] * sig_s;
            }
            fe_u[i] += w_vol * fu;
        }
    }
}

#[inline]
fn push_piezo_cell_blocks_real(
    coo: &mut Coo,
    d: usize,
    n1: usize,
    cell_nodes: &[usize],
    mat: &PiezoCellMat,
    ke_uu: &[f64],
    ke_up: &[f64],
    ke_pp: &[f64],
    u_dofs: &mut Vec<usize>,
    p_dofs: &mut Vec<usize>,
) {
    let nloc = cell_nodes.len();
    let nd = d * nloc;
    u_dofs.clear();
    u_dofs.reserve(nd);
    for comp in 0..d {
        for &gd in cell_nodes {
            u_dofs.push(comp * n1 + gd);
        }
    }
    p_dofs.clear();
    p_dofs.extend(cell_nodes.iter().map(|&gd| d * n1 + gd));

    if mat.mech_active {
        for i in 0..nd {
            let ri = u_dofs[i];
            for j in i..nd {
                let v = ke_uu[i * nd + j];
                if v != 0.0 {
                    let cj = u_dofs[j];
                    coo.push(ri, cj, v);
                    if i != j {
                        coo.push(cj, ri, v);
                    }
                }
            }
        }
        if mat.elec_active && ke_up.iter().any(|&v| v != 0.0) {
            for (r_i, &ri) in u_dofs.iter().enumerate() {
                for (c_j, &cj) in p_dofs.iter().enumerate() {
                    let v = ke_up[r_i * nloc + c_j];
                    if v != 0.0 {
                        coo.push(ri, cj, v);
                        coo.push(cj, ri, v);
                    }
                }
            }
        }
    }
    if mat.elec_active {
        for a in 0..nloc {
            let ra = p_dofs[a];
            for b_idx in a..nloc {
                let v = ke_pp[a * nloc + b_idx];
                if v != 0.0 {
                    let cb = p_dofs[b_idx];
                    coo.push(ra, cb, v);
                    if a != b_idx {
                        coo.push(cb, ra, v);
                    }
                }
            }
        }
    }
}

#[inline]
fn push_piezo_cell_blocks_complex(
    coo: &mut CCoo,
    d: usize,
    n1: usize,
    cell_nodes: &[usize],
    mat: &PiezoCellMat,
    ke_uu: &[C],
    ke_up: &[C],
    ke_pp: &[C],
    u_dofs: &mut Vec<usize>,
    p_dofs: &mut Vec<usize>,
) {
    let nloc = cell_nodes.len();
    let nd = d * nloc;
    u_dofs.clear();
    u_dofs.reserve(nd);
    for comp in 0..d {
        for &gd in cell_nodes {
            u_dofs.push(comp * n1 + gd);
        }
    }
    p_dofs.clear();
    p_dofs.extend(cell_nodes.iter().map(|&gd| d * n1 + gd));

    if mat.mech_active {
        for i in 0..nd {
            let ri = u_dofs[i];
            for j in i..nd {
                let v = ke_uu[i * nd + j];
                if v != c(0.0, 0.0) {
                    let cj = u_dofs[j];
                    coo.push(ri, cj, v);
                    if i != j {
                        coo.push(cj, ri, v);
                    }
                }
            }
        }
        if mat.elec_active && ke_up.iter().any(|&v| v != c(0.0, 0.0)) {
            for (r_i, &ri) in u_dofs.iter().enumerate() {
                for (c_j, &cj) in p_dofs.iter().enumerate() {
                    let v = ke_up[r_i * nloc + c_j];
                    if v != c(0.0, 0.0) {
                        coo.push(ri, cj, v);
                        coo.push(cj, ri, v);
                    }
                }
            }
        }
    }
    if mat.elec_active {
        for a in 0..nloc {
            let ra = p_dofs[a];
            for b_idx in a..nloc {
                let v = ke_pp[a * nloc + b_idx];
                if v != c(0.0, 0.0) {
                    let cb = p_dofs[b_idx];
                    coo.push(ra, cb, v);
                    if a != b_idx {
                        coo.push(cb, ra, v);
                    }
                }
            }
        }
    }
}

/// Deterministic chunked parallel real coupled matrix assembly on [`SimplexMesh`].
pub fn assemble_piezo_real_simplex(
    mesh: &SimplexMesh,
    mats: &[PiezoCellMat],
    order: usize,
    plane: Plane,
    thickness: f64,
) -> RCsr {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let n = (d + 1) * n1;
    let nc = mesh.n_cells();
    let thk = if d == 2 { thickness } else { 1.0 };
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };

    let chunk_size = 512usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let parallel = par::worth_parallel(nc);

    let blocks: Vec<Coo> = par::chunk_map(nchunks, parallel, |cid| {
        let mut coo = Coo::new(n, n);
        let mut scratch = Vec::new();
        let mut u_dofs = Vec::new();
        let mut p_dofs = Vec::new();
        let mut grad = vec![0.0; nloc * d];
        let mut b = vec![0.0; ns * nd];
        let mut ke_uu = vec![0.0; nd * nd];
        let mut ke_up = vec![0.0; nd * nloc];
        let mut ke_pp = vec![0.0; nloc * nloc];

        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            let mat = &mats[ci];
            if !mat.mech_active && !mat.elec_active {
                continue;
            }
            ke_uu.fill(0.0);
            ke_up.fill(0.0);
            ke_pp.fill(0.0);
            let vol_thk = geo.vol[ci] * thk;
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, ci, q, &mut grad);
                accumulate_piezo_qp_real(
                    mat,
                    d,
                    plane,
                    nloc,
                    &grad,
                    re.w[q] * vol_thk,
                    &mut b,
                    &mut ke_uu,
                    &mut ke_up,
                    &mut ke_pp,
                );
            }
            mesh.cell_dofs(ci, order, &mut scratch);
            push_piezo_cell_blocks_real(
                &mut coo,
                d,
                n1,
                &scratch[..nloc],
                mat,
                &ke_uu,
                &ke_up,
                &ke_pp,
                &mut u_dofs,
                &mut p_dofs,
            );
        }
        coo
    });

    let mut all = Coo::new(n, n);
    for blk in blocks {
        all.extend(&blk);
    }
    all.to_csr()
}

/// Deterministic chunked parallel real coupled matrix assembly on [`MixedDiscretisation`]
/// (supporting straight/curved cells and static `InfiniteElement` domains with optional outer exponential Robin).
pub(crate) fn assemble_piezo_real_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    mats: &[PiezoCellMat],
    plane: Plane,
    thickness: f64,
    ie_domains: &ResolvedIeDomains,
    inactive: &[i32],
    conductor_domains: &[i32],
) -> RCsr {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = (d + 1) * n1;
    let nc = mixed.n_cells();
    let thk = if d == 2 { thickness } else { 1.0 };
    let ns = if d == 2 { 3 } else { 6 };

    let chunk_size = 256usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let parallel = par::worth_parallel(nc);

    let blocks: Vec<Coo> = par::chunk_map(nchunks, parallel, |cid| {
        let mut coo = Coo::new(n, n);
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut grad = Vec::new();
        let mut b = Vec::new();
        let mut ke_uu = Vec::new();
        let mut ke_up = Vec::new();
        let mut ke_pp = Vec::new();
        let mut u_dofs = Vec::new();
        let mut p_dofs = Vec::new();

        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            let mat = &mats[ci];
            if !mat.mech_active && !mat.elec_active {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let ct = mixed.cell_type(ci);
            let iso = IsoparametricCell::from_mesh(mixed, ci);
            let in_ie = ie_domains.contains(dom);
            let deg = if in_ie {
                10
            } else if iso.curved {
                2 * order + 4
            } else {
                2 * order + 2
            };
            let rule = quad_rule_mixed(ct, deg);
            let cdofs = dof_map.cell_dofs(ci);
            let nloc = cdofs.len();
            let nd = d * nloc;

            grad.resize(nloc * d, 0.0);
            b.resize(ns * nd, 0.0);
            ke_uu.fill(0.0);
            ke_uu.resize(nd * nd, 0.0);
            ke_up.fill(0.0);
            ke_up.resize(nd * nloc, 0.0);
            ke_pp.fill(0.0);
            ke_pp.resize(nloc * nloc, 0.0);

            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                let w_eff = if in_ie {
                    let ie_j = ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, d);
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        for i in 0..d {
                            let mut acc = 0.0;
                            for k in 0..d {
                                acc += ie_j.j_inv[k][i] * gp[k];
                            }
                            grad[a * d + i] = acc;
                        }
                    }
                    wq * jg.det.abs() * ie_j.det * thk
                } else {
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        grad[a * d..a * d + d].copy_from_slice(&gp[..d]);
                    }
                    wq * jg.det.abs() * thk
                };
                accumulate_piezo_qp_real(mat, d, plane, nloc, &grad, w_eff, &mut b, &mut ke_uu, &mut ke_up, &mut ke_pp);
            }
            push_piezo_cell_blocks_real(&mut coo, d, n1, cdofs, mat, &ke_uu, &ke_up, &ke_pp, &mut u_dofs, &mut p_dofs);
        }
        coo
    });

    let mut all = Coo::new(n, n);
    for blk in blocks {
        all.extend(&blk);
    }
    let k_vol = all.to_csr();
    if !ie_domains.has_exponential_robin() {
        return k_vol;
    }

    // Mechanical spring Robin on outer faces of Exponential IE
    let c6: Vec<_> = mats.iter().map(|m| m.c6).collect();
    let k_mech_robin = assemble_exponential_ie_robin(mixed, dof_map, &c6, ie_domains, inactive, thk);

    // Electrical 1/R_tilde decay Robin on outer faces of Exponential IE (`-(n . eps . n / R_tilde) phi`)
    let mut k_elec_coo = Coo::new(n, n);
    let faces = collect_exponential_ie_faces(mixed, ie_domains, &[], conductor_domains);
    if !faces.is_empty() {
        let deg = 2 * order + 4;
        let mut fdofs = Vec::new();
        for (ci, ft, fv) in faces {
            let eps_m3 = mats[ci].eps_re;
            raw_face_dofs(mixed, dof_map, ft, &fv, &mut fdofs);
            let nloc = fdofs.len();
            let rule = facet_rule_mixed(ft, deg);
            let mut blk = vec![0.0; nloc * nloc];
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let (x_m, ds, n_hat) = eval_raw_face_jacobian(mixed, ft, &fv, xi);
                if ds <= 0.0 {
                    continue;
                }
                let x_w = wafer_point(&x_m[..d], d);
                let ie_j = ie_domains.eval_jacobian(mixed.cell_labels[ci], x_w, d);
                let mut n_tilde = [0.0; 3];
                for i in 0..d {
                    let s: f64 = (0..d).map(|k| ie_j.j_inv[k][i] * n_hat[k]).sum();
                    n_tilde[i] = ie_j.det.abs() * s;
                }
                let area_scale = (0..d).map(|i| n_tilde[i] * n_tilde[i]).sum::<f64>().sqrt();
                if area_scale <= 0.0 {
                    continue;
                }
                let nn: [f64; 3] = std::array::from_fn(|i| if i < d { n_tilde[i] / area_scale } else { 0.0 });
                let mut eps_n = 0.0;
                for i in 0..d {
                    for j in 0..d {
                        eps_n += nn[i] * eps_m3[i][j] * nn[j];
                    }
                }
                let coeff = wq * ds * area_scale * thk * eps_n / ie_j.r_tilde.max(1e-30);
                let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
                for a in 0..nloc {
                    let wa = coeff * nv[a];
                    for b in a..nloc {
                        blk[a * nloc + b] -= wa * nv[b];
                    }
                }
            }
            for a in 0..nloc {
                let ra = d * n1 + fdofs[a];
                for b in a..nloc {
                    let v = blk[a * nloc + b];
                    if v != 0.0 {
                        let cb = d * n1 + fdofs[b];
                        k_elec_coo.push(ra, cb, v);
                        if a != b {
                            k_elec_coo.push(cb, ra, v);
                        }
                    }
                }
            }
        }
    }

    // Pad mechanical Robin `(d*n1)` into `(d+1)*n1`
    let mut k_mr_coo = Coo::new(n, n);
    for r in 0..k_mech_robin.nrows() {
        let (cols, vals) = k_mech_robin.row(r);
        for (&col, &val) in cols.iter().zip(vals) {
            k_mr_coo.push(r, col as usize, val);
        }
    }
    let k_mr_padded = k_mr_coo.to_csr();
    let k_er = k_elec_coo.to_csr();
    RCsr::axpby(1.0, &RCsr::axpby(1.0, &k_vol, 1.0, &k_mr_padded), 1.0, &k_er)
}

fn thermal_load_piezo_simplex(
    mesh: &SimplexMesh,
    mats: &[PiezoCellMat],
    alpha: &[M3],
    t_field: &[f64],
    t_ref: f64,
    order: usize,
    plane: Plane,
    thickness: f64,
) -> Vec<f64> {
    let d = mesh.dim;
    let re = RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let n1 = mesh.n_dofs(order);
    let n = (d + 1) * n1;
    let nc = mesh.n_cells();
    let thk = if d == 2 { thickness } else { 1.0 };

    let chunk_size = 512usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let parts = par::chunk_map(nchunks, par::worth_parallel(nc), |cid| {
        let mut f = vec![0.0; n];
        let mut scratch = Vec::new();
        let mut grad = vec![0.0; nloc * d];
        let mut b = vec![0.0; ns * nd];
        let mut fe_u = vec![0.0; nd];
        let mut fe_p = vec![0.0; nloc];

        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            let mat = &mats[ci];
            if !mat.mech_active {
                continue;
            }
            mesh.cell_dofs(ci, order, &mut scratch);
            fe_u.fill(0.0);
            fe_p.fill(0.0);
            let vol_thk = geo.vol[ci] * thk;
            for q in 0..re.w.len() {
                let tq: f64 = (0..nloc).map(|a| re.n[q][a] * t_field[scratch[a]]).sum();
                shape_gradients(&re, &geo, ci, q, &mut grad);
                accumulate_piezo_thermal_qp(
                    mat,
                    d,
                    plane,
                    nloc,
                    &grad,
                    &alpha[ci],
                    tq - t_ref,
                    re.w[q] * vol_thk,
                    &mut b,
                    &mut fe_u,
                    &mut fe_p,
                );
            }
            for comp in 0..d {
                for a in 0..nloc {
                    f[comp * n1 + scratch[a]] += fe_u[comp * nloc + a];
                }
            }
            if mat.elec_active {
                for a in 0..nloc {
                    f[d * n1 + scratch[a]] += fe_p[a];
                }
            }
        }
        f
    });

    let mut out = vec![0.0; n];
    for p in parts {
        for i in 0..n {
            out[i] += p[i];
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn thermal_load_piezo_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    mats: &[PiezoCellMat],
    alpha: &[M3],
    t_field: &[f64],
    t_ref: f64,
    plane: Plane,
    thickness: f64,
    ie_domains: &ResolvedIeDomains,
) -> Vec<f64> {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = (d + 1) * n1;
    let nc = mixed.n_cells();
    let thk = if d == 2 { thickness } else { 1.0 };
    let ns = if d == 2 { 3 } else { 6 };

    let chunk_size = 256usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let parts = par::chunk_map(nchunks, par::worth_parallel(nc), |cid| {
        let mut f = vec![0.0; n];
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut grad = Vec::new();
        let mut b = Vec::new();
        let mut fe_u = Vec::new();
        let mut fe_p = Vec::new();

        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            let mat = &mats[ci];
            if !mat.mech_active {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let ct = mixed.cell_type(ci);
            let iso = IsoparametricCell::from_mesh(mixed, ci);
            let in_ie = ie_domains.contains(dom);
            let deg = if in_ie {
                10
            } else if iso.curved {
                2 * order + 4
            } else {
                2 * order + 2
            };
            let rule = quad_rule_mixed(ct, deg);
            let cdofs = dof_map.cell_dofs(ci);
            let nloc = cdofs.len();
            let nd = d * nloc;

            grad.resize(nloc * d, 0.0);
            b.resize(ns * nd, 0.0);
            fe_u.fill(0.0);
            fe_u.resize(nd, 0.0);
            fe_p.fill(0.0);
            fe_p.resize(nloc, 0.0);

            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                let tq: f64 = cdofs.iter().enumerate().map(|(a, &gd)| vals[a] * t_field[gd]).sum();
                let w_vol = if in_ie {
                    let ie_j = ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, d);
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        for i in 0..d {
                            let mut acc = 0.0;
                            for k in 0..d {
                                acc += ie_j.j_inv[k][i] * gp[k];
                            }
                            grad[a * d + i] = acc;
                        }
                    }
                    wq * jg.det.abs() * ie_j.det * thk
                } else {
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        grad[a * d..a * d + d].copy_from_slice(&gp[..d]);
                    }
                    wq * jg.det.abs() * thk
                };
                accumulate_piezo_thermal_qp(
                    mat,
                    d,
                    plane,
                    nloc,
                    &grad,
                    &alpha[ci],
                    tq - t_ref,
                    w_vol,
                    &mut b,
                    &mut fe_u,
                    &mut fe_p,
                );
            }
            for comp in 0..d {
                for a in 0..nloc {
                    f[comp * n1 + cdofs[a]] += fe_u[comp * nloc + a];
                }
            }
            if mat.elec_active {
                for a in 0..nloc {
                    f[d * n1 + cdofs[a]] += fe_p[a];
                }
            }
        }
        f
    });

    let mut out = vec![0.0; n];
    for p in parts {
        for i in 0..n {
            out[i] += p[i];
        }
    }
    out
}

fn check_static_mech_support(spec: &Problem, has_exp_robin: bool) -> Result<(), String> {
    let has_mech_fixed = spec.bcs.iter().any(|bc| match bc.kind {
        BcKind::Fixed { fixed, .. } => fixed.iter().any(|&f| f),
        _ => false,
    });
    if !has_mech_fixed && !has_exp_robin && spec.infinite_elements.is_empty() {
        return Err("a static solid problem needs a fixed boundary (otherwise it floats freely)".into());
    }
    Ok(())
}

/// Assembles Stage 1 static coupled piezoelectric problem on [`SimplexMesh`].
pub fn assemble_piezo_static_simplex(mesh: &SimplexMesh, spec: &Problem) -> Result<PiezoStaticAssembly, String> {
    crate::problem::check_simplex_straight(mesh)?;
    if !matches!(spec.study, Study::Static) {
        return Err("assemble_piezo_static requires Study::Static".into());
    }
    check_labels_simplex(mesh, spec)?;
    check_static_mech_support(spec, false)?;

    let d = mesh.dim;
    let order = spec.order;
    let n1 = mesh.n_dofs(order);
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let (by_dom, mats) = resolve_piezo_mats_simplex(mesh, spec, false)?;
    let mech_untouched = untouched_dofs(mesh, order, d, &spec.inactive);
    let get_cell = |ci: usize, out: &mut Vec<usize>| {
        mesh.cell_dofs(ci, order, out);
    };
    let get_bd = |labels: &[i32]| boundary_dofs(mesh, labels, order);
    let get_cd = |labels: &[i32]| collect_conductor_domain_dofs_simplex(mesh, order, labels);
    let constraints = build_piezo_constraints(
        d,
        n1,
        spec,
        &mesh.cell_labels,
        &get_cell,
        &get_bd,
        &get_cd,
        mech_untouched,
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )?;

    let k_full = assemble_piezo_real_simplex(mesh, &mats, order, spec.plane, spec.thickness);
    let k_red = congruence_real(&constraints.red_all, &k_full)?;

    let mut f_mech_elec = vec![0.0; constraints.n_tot];
    let axes = wafer_axes(d);
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::Load(t) => {
                for (comp, &w) in axes.iter().enumerate() {
                    if t[w] != 0.0 {
                        let tc = t[w] * thk;
                        add_boundary_load(mesh, &bc.labels, &|_| tc, order, &mut f_mech_elec, d + 1, comp);
                    }
                }
            }
            BcKind::SurfaceCharge(sigma_s) if *sigma_s != 0.0 => {
                let qc = -sigma_s * thk;
                add_boundary_load(mesh, &bc.labels, &|_| qc, order, &mut f_mech_elec, d + 1, d);
            }
            _ => {}
        }
    }
    for el in &constraints.electrodes {
        if let (true, Some(q), Some(root)) = (el.is_floating, el.prescribed_charge, el.root_dof) {
            if q != 0.0 {
                f_mech_elec[constraints.n_mech + root] -= q;
            }
        }
    }

    let mut thermal = None;
    let mut f_th = vec![0.0; constraints.n_tot];
    if let Some(th) = &spec.thermal {
        let alpha: Vec<M3> = spec
            .materials
            .cell_tensor_skipping(mesh, TensorKind::Expansion, Some(0.0), &spec.inactive)?
            .iter()
            .map(|a| plane_tensor(a, d))
            .collect();
        let t = match &th.source {
            ThermalSource::Uniform(tv) => vec![*tv; n1],
            ThermalSource::Imported(t) if t.len() == n1 => t.clone(),
            ThermalSource::Imported(_) => {
                return Err("no heat-transfer result on this mesh and element order; solve Heat transfer first".into())
            }
        };
        f_th = thermal_load_piezo_simplex(mesh, &mats, &alpha, &t, th.t_ref, order, spec.plane, spec.thickness);
        thermal = Some(Thermal { alpha, t, tref: th.t_ref });
    }

    let f_full: Vec<f64> = f_mech_elec.iter().zip(&f_th).map(|(&a, &b)| a + b).collect();
    let g_real: Vec<f64> = constraints.g_all.iter().map(|z| z.re).collect();
    let rhs_red = constraints.red_all.lift_rhs(&k_full, &f_full, &g_real);

    Ok(PiezoStaticAssembly {
        constraints,
        mats,
        by_dom,
        k_full,
        k_red,
        f_full,
        f_mech_elec_load_only: f_mech_elec,
        f_th_only: f_th,
        rhs_red,
        thermal,
        ie_domains: ResolvedIeDomains::default(),
        diagnostics: Vec::new(),
    })
}

/// Assembles Stage 1 static coupled piezoelectric problem natively on [`MixedDiscretisation`].
pub fn assemble_piezo_static_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
) -> Result<PiezoStaticAssembly, String> {
    if !matches!(spec.study, Study::Static) {
        return Err("assemble_piezo_static_mixed requires Study::Static".into());
    }
    mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let ie_domains = ResolvedIeDomains::from_mixed(mixed, spec);
    check_static_mech_support(spec, ie_domains.has_exponential_robin())?;

    let (by_dom, mats) = resolve_piezo_mats_mixed(mixed, spec, false)?;
    let mech_untouched = mixed_untouched_dofs(mixed, dof_map, d, &spec.inactive);
    let ie_clamped = ie_domains.outer_clamped_dofs(mixed, dof_map, spec);
    let conductor_domains = spec.conductor_domains();
    let ie_robin_nodes: Vec<usize> = if ie_domains.has_exponential_robin() {
        let mut s = BTreeSet::new();
        let mut fdofs = Vec::new();
        for (_ci, ft, fv) in collect_exponential_ie_faces(mixed, &ie_domains, &[], &conductor_domains) {
            raw_face_dofs(mixed, dof_map, ft, &fv, &mut fdofs);
            s.extend(fdofs.iter().copied());
        }
        s.into_iter().collect()
    } else {
        Vec::new()
    };

    let get_cell = |ci: usize, out: &mut Vec<usize>| {
        out.clear();
        out.extend_from_slice(dof_map.cell_dofs(ci));
    };
    let get_bd = |labels: &[i32]| dof_map.boundary_dofs(mixed, labels);
    let get_cd = |labels: &[i32]| collect_conductor_domain_dofs_mixed(mixed, dof_map, labels);
    let constraints = build_piezo_constraints(
        d,
        n1,
        spec,
        &mixed.cell_labels,
        &get_cell,
        &get_bd,
        &get_cd,
        mech_untouched,
        ie_clamped,
        ie_robin_nodes,
        Vec::new(),
    )?;

    let k_full = assemble_piezo_real_mixed(
        mixed,
        dof_map,
        &mats,
        spec.plane,
        spec.thickness,
        &ie_domains,
        &spec.inactive,
        &constraints.conductor_domains,
    );
    let k_red = congruence_real(&constraints.red_all, &k_full)?;

    let mut f_mech_elec = vec![0.0; constraints.n_tot];
    let axes = wafer_axes(d);
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::Load(t) => {
                for (comp, &w) in axes.iter().enumerate() {
                    if t[w] != 0.0 {
                        let tc = t[w] * thk;
                        add_boundary_load_mixed(mixed, dof_map, &bc.labels, tc, &mut f_mech_elec, comp);
                    }
                }
            }
            BcKind::SurfaceCharge(sigma_s) if *sigma_s != 0.0 => {
                let qc = -sigma_s * thk;
                add_boundary_load_mixed(
                    mixed,
                    dof_map,
                    &bc.labels,
                    qc,
                    &mut f_mech_elec[constraints.n_mech..constraints.n_tot],
                    0,
                );
            }
            _ => {}
        }
    }
    for el in &constraints.electrodes {
        if let (true, Some(q), Some(root)) = (el.is_floating, el.prescribed_charge, el.root_dof) {
            if q != 0.0 {
                f_mech_elec[constraints.n_mech + root] -= q;
            }
        }
    }

    let mut thermal = None;
    let mut f_th = vec![0.0; constraints.n_tot];
    if let Some(th) = &spec.thermal {
        let all_alpha = mixed_cell_tensor(&spec.materials, mixed, TensorKind::Expansion, Some(0.0))?;
        let alpha: Vec<M3> = (0..mixed.n_cells())
            .map(|ci| {
                if spec.inactive.contains(&mixed.cell_labels[ci]) {
                    [[0.0; 3]; 3]
                } else {
                    plane_tensor(&all_alpha[ci], d)
                }
            })
            .collect();
        let t = match &th.source {
            ThermalSource::Uniform(tv) => vec![*tv; n1],
            ThermalSource::Imported(t) if t.len() == n1 => t.clone(),
            ThermalSource::Imported(_) => {
                return Err("no heat-transfer result on this mesh and element order; solve Heat transfer first".into())
            }
        };
        f_th = thermal_load_piezo_mixed(
            mixed,
            dof_map,
            &mats,
            &alpha,
            &t,
            th.t_ref,
            spec.plane,
            spec.thickness,
            &ie_domains,
        );
        thermal = Some(Thermal { alpha, t, tref: th.t_ref });
    }

    let f_full: Vec<f64> = f_mech_elec.iter().zip(&f_th).map(|(&a, &b)| a + b).collect();
    let g_real: Vec<f64> = constraints.g_all.iter().map(|z| z.re).collect();
    let rhs_red = constraints.red_all.lift_rhs(&k_full, &f_full, &g_real);
    let diagnostics = ie_domains.diagnostics.clone();

    Ok(PiezoStaticAssembly {
        constraints,
        mats,
        by_dom,
        k_full,
        k_red,
        f_full,
        f_mech_elec_load_only: f_mech_elec,
        f_th_only: f_th,
        rhs_red,
        thermal,
        ie_domains,
        diagnostics,
    })
}

fn assemble_piezo_complex_simplex(
    mesh: &SimplexMesh,
    mats: &[PiezoCellMat],
    order: usize,
    plane: Plane,
    thickness: f64,
    resolved_pml: &ResolvedPmlDomains,
) -> CCsr {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let n = (d + 1) * n1;
    let nc = mesh.n_cells();
    let thk = if d == 2 { thickness } else { 1.0 };
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let re_pml = RefElement::new(d, order, 2 * order + 2);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };

    let chunk_size = 512usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let blocks: Vec<CCoo> = par::chunk_map(nchunks, par::worth_parallel(nc), |cid| {
        let mut coo = CCoo::new(n, n);
        let mut scratch = Vec::new();
        let mut u_dofs = Vec::new();
        let mut p_dofs = Vec::new();
        let mut grad = vec![0.0; nloc * d];
        let mut gc = vec![c(0.0, 0.0); nloc * d];
        let mut bc = vec![c(0.0, 0.0); ns * nd];
        let mut ke_uu = vec![c(0.0, 0.0); nd * nd];
        let mut ke_up = vec![c(0.0, 0.0); nd * nloc];
        let mut ke_pp = vec![c(0.0, 0.0); nloc * nloc];

        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            let mat = &mats[ci];
            if !mat.mech_active && !mat.elec_active {
                continue;
            }
            let dom = mesh.cell_labels[ci];
            ke_uu.fill(c(0.0, 0.0));
            ke_up.fill(c(0.0, 0.0));
            ke_pp.fill(c(0.0, 0.0));
            let vol_thk = geo.vol[ci] * thk;
            if resolved_pml.contains(dom) {
                let cv = mesh.cell(ci);
                for q in 0..re_pml.w.len() {
                    shape_gradients(&re_pml, &geo, ci, q, &mut grad);
                    let lam = simplex_barycentric(&re_pml, q, d, order);
                    let mut x_phys = [0.0; 3];
                    for corner in 0..=d {
                        let pt = mesh.point(cv[corner] as usize);
                        for k in 0..d {
                            x_phys[k] += lam[corner] * pt[k];
                        }
                    }
                    let x_w = wafer_point(&x_phys[..d], d);
                    let (jac, loss) = resolved_pml.eval_jacobian(dom, x_w, d, false);
                    for a in 0..nloc {
                        for i in 0..d {
                            let mut acc = c(0.0, 0.0);
                            for k in 0..d {
                                acc += jac.j_inv[k][i] * grad[a * d + k];
                            }
                            gc[a * d + i] = acc;
                        }
                    }
                    let w_coup = jac.det * (re_pml.w[q] * vol_thk);
                    let w_k = w_coup * loss;
                    accumulate_piezo_qp_complex(
                        mat, d, plane, nloc, &gc, w_k, w_coup, &mut bc, &mut ke_uu, &mut ke_up, &mut ke_pp,
                    );
                }
            } else {
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, ci, q, &mut grad);
                    for a in 0..nloc * d {
                        gc[a] = c(grad[a], 0.0);
                    }
                    let w = c(re.w[q] * vol_thk, 0.0);
                    accumulate_piezo_qp_complex(
                        mat, d, plane, nloc, &gc, w, w, &mut bc, &mut ke_uu, &mut ke_up, &mut ke_pp,
                    );
                }
            }
            mesh.cell_dofs(ci, order, &mut scratch);
            push_piezo_cell_blocks_complex(
                &mut coo,
                d,
                n1,
                &scratch[..nloc],
                mat,
                &ke_uu,
                &ke_up,
                &ke_pp,
                &mut u_dofs,
                &mut p_dofs,
            );
        }
        coo
    });

    let mut all = CCoo::new(n, n);
    for blk in blocks {
        all.extend(&blk);
    }
    all.to_csr()
}

fn assemble_piezo_complex_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    mats: &[PiezoCellMat],
    plane: Plane,
    thickness: f64,
    resolved_pml: &ResolvedPmlDomains,
) -> CCsr {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = (d + 1) * n1;
    let nc = mixed.n_cells();
    let thk = if d == 2 { thickness } else { 1.0 };
    let ns = if d == 2 { 3 } else { 6 };

    let chunk_size = 256usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let blocks: Vec<CCoo> = par::chunk_map(nchunks, par::worth_parallel(nc), |cid| {
        let mut coo = CCoo::new(n, n);
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut gc = Vec::new();
        let mut bc = Vec::new();
        let mut ke_uu = Vec::new();
        let mut ke_up = Vec::new();
        let mut ke_pp = Vec::new();
        let mut u_dofs = Vec::new();
        let mut p_dofs = Vec::new();

        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            let mat = &mats[ci];
            if !mat.mech_active && !mat.elec_active {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let ct = mixed.cell_type(ci);
            let iso = IsoparametricCell::from_mesh(mixed, ci);
            let in_pml = resolved_pml.contains(dom);
            let deg = if in_pml || iso.curved { 2 * order + 4 } else { 2 * order + 2 };
            let rule = quad_rule_mixed(ct, deg);
            let cdofs = dof_map.cell_dofs(ci);
            let nloc = cdofs.len();
            let nd = d * nloc;

            gc.resize(nloc * d, c(0.0, 0.0));
            bc.resize(ns * nd, c(0.0, 0.0));
            ke_uu.fill(c(0.0, 0.0));
            ke_uu.resize(nd * nd, c(0.0, 0.0));
            ke_up.fill(c(0.0, 0.0));
            ke_up.resize(nd * nloc, c(0.0, 0.0));
            ke_pp.fill(c(0.0, 0.0));
            ke_pp.resize(nloc * nloc, c(0.0, 0.0));

            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                let (w_k, w_coup) = if in_pml {
                    let x_w = wafer_point(&jg.x[..d], d);
                    let (jac, loss) = resolved_pml.eval_jacobian(dom, x_w, d, false);
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        for i in 0..d {
                            let mut acc = c(0.0, 0.0);
                            for k in 0..d {
                                acc += jac.j_inv[k][i] * gp[k];
                            }
                            gc[a * d + i] = acc;
                        }
                    }
                    let wc = jac.det * (wq * jg.det.abs() * thk);
                    (wc * loss, wc)
                } else {
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        for i in 0..d {
                            gc[a * d + i] = c(gp[i], 0.0);
                        }
                    }
                    let w = c(wq * jg.det.abs() * thk, 0.0);
                    (w, w)
                };
                accumulate_piezo_qp_complex(
                    mat, d, plane, nloc, &gc, w_k, w_coup, &mut bc, &mut ke_uu, &mut ke_up, &mut ke_pp,
                );
            }
            push_piezo_cell_blocks_complex(
                &mut coo,
                d,
                n1,
                cdofs,
                mat,
                &ke_uu,
                &ke_up,
                &ke_pp,
                &mut u_dofs,
                &mut p_dofs,
            );
        }
        coo
    });

    let mut all = CCoo::new(n, n);
    for blk in blocks {
        all.extend(&blk);
    }
    all.to_csr()
}

fn assemble_pml_mass_simplex(
    mesh: &SimplexMesh,
    rho: &[f64],
    order: usize,
    thickness: f64,
    resolved_pml: &ResolvedPmlDomains,
) -> CCsr {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let n = d * n1;
    let thk = if d == 2 { thickness } else { 1.0 };
    let rm = RefElement::new(d, order, 2 * order);
    let re_pml = RefElement::new(d, order, 2 * order + 2);
    let geo = geometry(mesh);
    let nloc = rm.nloc;

    let mut mloc = vec![0.0; nloc * nloc];
    for q in 0..rm.w.len() {
        for a in 0..nloc {
            for b in 0..nloc {
                mloc[a * nloc + b] += rm.w[q] * rm.n[q][a] * rm.n[q][b];
            }
        }
    }

    let nc = mesh.n_cells();
    let chunk_size = 512usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let blocks = par::chunk_map(nchunks, par::worth_parallel(nc), |cid| {
        let mut mc = CCoo::new(n, n);
        let mut scratch = Vec::new();
        let mut me = vec![c(0.0, 0.0); nloc * nloc];
        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            if rho[ci] <= 0.0 {
                continue;
            }
            let dom = mesh.cell_labels[ci];
            let vol_thk = geo.vol[ci] * thk;
            mesh.cell_dofs(ci, order, &mut scratch);
            if resolved_pml.contains(dom) {
                me.fill(c(0.0, 0.0));
                let cv = mesh.cell(ci);
                for q in 0..re_pml.w.len() {
                    let lam = simplex_barycentric(&re_pml, q, d, order);
                    let mut x_phys = [0.0; 3];
                    for corner in 0..=d {
                        let pt = mesh.point(cv[corner] as usize);
                        for k in 0..d {
                            x_phys[k] += lam[corner] * pt[k];
                        }
                    }
                    let x_w = wafer_point(&x_phys[..d], d);
                    let (jac, _) = resolved_pml.eval_jacobian(dom, x_w, d, false);
                    let w_m = jac.det * (re_pml.w[q] * vol_thk * rho[ci]);
                    for a in 0..nloc {
                        for b in 0..nloc {
                            me[a * nloc + b] += w_m * (re_pml.n[q][a] * re_pml.n[q][b]);
                        }
                    }
                }
            } else {
                for a in 0..nloc * nloc {
                    me[a] = c(rho[ci] * vol_thk * mloc[a], 0.0);
                }
            }
            for comp in 0..d {
                let dofs: Vec<usize> = scratch[..nloc].iter().map(|&x| x + comp * n1).collect();
                mc.push_block(&dofs, &dofs, &me);
            }
        }
        mc
    });
    let mut all = CCoo::new(n, n);
    for b in blocks {
        all.extend(&b);
    }
    all.to_csr()
}

fn assemble_pml_mass_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    rho: &[f64],
    thickness: f64,
    resolved_pml: &ResolvedPmlDomains,
) -> CCsr {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = d * n1;
    let thk = if d == 2 { thickness } else { 1.0 };
    let nc = mixed.n_cells();
    let chunk_size = 256usize;
    let nchunks = nc.div_ceil(chunk_size).max(1);
    let blocks = par::chunk_map(nchunks, par::worth_parallel(nc), |cid| {
        let mut mc = CCoo::new(n, n);
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut me = Vec::new();
        let c0 = cid * chunk_size;
        let c1 = ((cid + 1) * chunk_size).min(nc);
        for ci in c0..c1 {
            if rho[ci] <= 0.0 {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let ct = mixed.cell_type(ci);
            let iso = IsoparametricCell::from_mesh(mixed, ci);
            let in_pml = resolved_pml.contains(dom);
            let deg = if in_pml || iso.curved { 2 * order + 4 } else { 2 * order + 2 };
            let rule = quad_rule_mixed(ct, deg);
            let cdofs = dof_map.cell_dofs(ci);
            let nloc = cdofs.len();
            me.fill(c(0.0, 0.0));
            me.resize(nloc * nloc, c(0.0, 0.0));
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                let w_m = if in_pml {
                    let x_w = wafer_point(&jg.x[..d], d);
                    let (jac, _) = resolved_pml.eval_jacobian(dom, x_w, d, false);
                    jac.det * (wq * jg.det.abs() * thk * rho[ci])
                } else {
                    c(wq * jg.det.abs() * thk * rho[ci], 0.0)
                };
                for a in 0..nloc {
                    for b in 0..nloc {
                        me[a * nloc + b] += w_m * (vals[a] * vals[b]);
                    }
                }
            }
            for comp in 0..d {
                let dofs: Vec<usize> = cdofs.iter().map(|&x| x + comp * n1).collect();
                mc.push_block(&dofs, &dofs, &me);
            }
        }
        mc
    });
    let mut all = CCoo::new(n, n);
    for b in blocks {
        all.extend(&b);
    }
    all.to_csr()
}

/// Extracts a sub-block CSR `a[rows, cols]` where `rows`, `cols` are ascending local indices.
pub(crate) fn sub_csr<T: femcore::linalg::scalar::Scalar>(
    a: &femcore::linalg::csr::Csr<T>,
    r_start: usize,
    r_end: usize,
    c_start: usize,
    c_end: usize,
) -> femcore::linalg::csr::Csr<T> {
    let nr = r_end - r_start;
    let nc = c_end - c_start;
    let mut indptr = Vec::with_capacity(nr + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut values = Vec::new();
    for r in r_start..r_end {
        let (cols, vals) = a.row(r);
        for (&col_u32, &val) in cols.iter().zip(vals) {
            let col = col_u32 as usize;
            if col >= c_start && col < c_end {
                indices.push(col - c_start);
                values.push(val);
            }
        }
        indptr.push(indices.len());
    }
    femcore::linalg::csr::Csr::from_raw(nr, nc, indptr, indices, values)
}

/// Assembles Stage 1 eigenfrequency coupled piezoelectric problem on [`SimplexMesh`].
pub fn assemble_piezo_eigen_simplex(mesh: &SimplexMesh, spec: &Problem) -> Result<PiezoEigenAssembly, String> {
    crate::problem::check_simplex_straight(mesh)?;
    let Study::Eigenfrequency { n_modes, f_search } = spec.study else {
        return Err("assemble_piezo_eigen requires Study::Eigenfrequency".into());
    };
    check_labels_simplex(mesh, spec)?;

    let d = mesh.dim;
    let order = spec.order;
    let n1 = mesh.n_dofs(order);
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let (by_dom, mats) = resolve_piezo_mats_simplex(mesh, spec, true)?;
    let rho: Vec<f64> = mats.iter().map(|m| m.rho).collect();
    let resolved_pml = ResolvedPmlDomains::from_simplex(mesh, spec, false, f_search)?;
    let any_pml = !resolved_pml.domains.is_empty();
    let lrb_labels: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| matches!(b.kind, BcKind::LowReflecting))
        .flat_map(|b| b.labels.iter().copied())
        .collect();
    let has_visc = mats.iter().any(|m| m.mech_active && m.eta.iter().flatten().any(|&v| v != 0.0));
    let has_damping = has_visc || !lrb_labels.is_empty();
    let has_dielectric_loss = mats.iter().any(|m| m.elec_active && m.eps_im.iter().flatten().any(|&v| v.abs() > 0.0));
    let has_periodic = !spec.periodic.is_empty();
    let is_hermitian = !any_pml && !has_dielectric_loss && !has_damping;
    let is_complex = any_pml || has_dielectric_loss || has_damping || has_periodic;

    let mech_untouched = untouched_dofs(mesh, order, d, &spec.inactive);
    let outer_clamped =
        if any_pml { outer_clamped_dofs_simplex(mesh, order, spec, &resolved_pml, false, 0.0) } else { Vec::new() };

    let mut periodic_rels = Vec::new();
    for p in &spec.periodic {
        let mt = match_vertices(mesh, p)?;
        periodic_rels.extend(nodal_relations(mesh, order, d + 1, &mt, p.phase)?);
    }

    let get_cell = |ci: usize, out: &mut Vec<usize>| {
        mesh.cell_dofs(ci, order, out);
    };
    let get_bd = |labels: &[i32]| boundary_dofs(mesh, labels, order);
    let get_cd = |labels: &[i32]| collect_conductor_domain_dofs_simplex(mesh, order, labels);
    let constraints = build_piezo_constraints(
        d,
        n1,
        spec,
        &mesh.cell_labels,
        &get_cell,
        &get_bd,
        &get_cd,
        mech_untouched,
        outer_clamped,
        Vec::new(),
        periodic_rels,
    )?;

    if n_modes < 1 || n_modes + 1 >= constraints.n_u_red {
        return Err("number of modes must be between 1 and the number of free mechanical dofs".into());
    }

    let diagnostics = resolved_pml.diagnostics.clone();
    let n_u_red = constraints.n_u_red;
    let n_tot_red = constraints.red_all.n_reduced();

    let (
        k_full_real,
        k_red_real,
        m_uu_red_real,
        k_full_c,
        k_red_c,
        m_uu_red_c,
        c_uu_red_c,
        k_pp_red_real,
        k_pp_red_c,
        k_pu_red_real,
        k_pu_red_c,
    ) = if !is_complex {
        let kf = assemble_piezo_real_simplex(mesh, &mats, order, spec.plane, spec.thickness);
        let kr = congruence_real(&constraints.red_all, &kf)?;
        let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let m1 = assemble_mass(mesh, order, &rho_t);
        let blocks: Vec<&RCsr> = (0..d).map(|_| &m1).collect();
        let m_uu = RCsr::block_diag(&blocks);
        let mr = congruence_real(&constraints.red_u, &m_uu)?;
        let k_pp_r = sub_csr(&kr, n_u_red, n_tot_red, n_u_red, n_tot_red);
        let k_pu_r = sub_csr(&kr, n_u_red, n_tot_red, 0, n_u_red);
        (Some(kf), Some(kr), Some(mr), None, None, None, None, Some(k_pp_r), None, Some(k_pu_r), None)
    } else {
        let kf_c = assemble_piezo_complex_simplex(mesh, &mats, order, spec.plane, spec.thickness, &resolved_pml);
        let kr_c = congruence_fast(&constraints.red_all, &kf_c);
        let m_uu_c = if any_pml {
            assemble_pml_mass_simplex(mesh, &rho, order, spec.thickness, &resolved_pml)
        } else {
            let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
            let m1 = assemble_mass(mesh, order, &rho_t);
            let blocks: Vec<&RCsr> = (0..d).map(|_| &m1).collect();
            RCsr::block_diag(&blocks).to_complex()
        };
        let mr_c = congruence_fast(&constraints.red_u, &m_uu_c);

        let c_uu_red_c = if has_damping {
            let mut c_uu_r = RCsr::zeros(constraints.n_mech, constraints.n_mech);
            if has_visc {
                let eta_wafer: Vec<_> = mats.iter().map(|m| m.eta).collect();
                let cv = crate::solid::assemble_solid(mesh, &eta_wafer, order, spec.plane, spec.thickness);
                c_uu_r = RCsr::axpby(1.0, &c_uu_r, 1.0, &cv);
            }
            if !lrb_labels.is_empty() {
                let c6_arr: Vec<_> = mats.iter().map(|m| m.c6).collect();
                let cl = damping_matrix(mesh, &c6_arr, &rho, &lrb_labels, order, d);
                let cl_sc = if d == 2 { thk } else { 1.0 };
                c_uu_r = RCsr::axpby(1.0, &c_uu_r, cl_sc, &cl);
            }
            Some(congruence_fast(&constraints.red_u, &c_uu_r.to_complex()))
        } else {
            None
        };

        let k_pp_r = sub_csr(&kr_c, n_u_red, n_tot_red, n_u_red, n_tot_red);
        let k_pu_r = sub_csr(&kr_c, n_u_red, n_tot_red, 0, n_u_red);
        (None, None, None, Some(kf_c), Some(kr_c), Some(mr_c), c_uu_red_c, None, Some(k_pp_r), None, Some(k_pu_r))
    };

    let pml_cells: Vec<bool> = (0..mesh.n_cells()).map(|ci| resolved_pml.contains(mesh.cell_labels[ci])).collect();
    let has_interior = pml_cells.iter().zip(&rho).any(|(&p, &r)| !p && r > 0.0);
    let (m1_phys_all, m1_phys_pml) = if any_pml {
        let rho_all_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let rho_pml_t: Vec<f64> = rho.iter().zip(&pml_cells).map(|(&r, &p)| if p { r * thk } else { 0.0 }).collect();
        (Some(assemble_mass(mesh, order, &rho_all_t)), Some(assemble_mass(mesh, order, &rho_pml_t)))
    } else {
        (None, None)
    };

    Ok(PiezoEigenAssembly {
        constraints,
        mats,
        by_dom,
        k_full_real,
        k_red_real,
        m_uu_red_real,
        k_full_c,
        k_red_c,
        m_uu_red_c,
        c_uu_red_c,
        k_pp_red_real,
        k_pp_red_c,
        k_pu_red_real,
        k_pu_red_c,
        m1_phys_all,
        m1_phys_pml,
        resolved_pml,
        is_complex,
        is_hermitian,
        has_damping,
        has_pml: any_pml,
        has_interior,
        n_modes,
        f_search,
        n_pairs: spec.periodic.len(),
        diagnostics,
    })
}

/// Assembles Stage 1 eigenfrequency coupled piezoelectric problem natively on [`MixedDiscretisation`].
pub fn assemble_piezo_eigen_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
) -> Result<PiezoEigenAssembly, String> {
    let Study::Eigenfrequency { n_modes, f_search } = spec.study else {
        return Err("assemble_piezo_eigen_mixed requires Study::Eigenfrequency".into());
    };
    mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let (by_dom, mats) = resolve_piezo_mats_mixed(mixed, spec, true)?;
    let rho: Vec<f64> = mats.iter().map(|m| m.rho).collect();
    let resolved_pml = ResolvedPmlDomains::from_mixed(mixed, spec, false, f_search)?;
    let any_pml = !resolved_pml.domains.is_empty();
    let lrb_labels: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| matches!(b.kind, BcKind::LowReflecting))
        .flat_map(|b| b.labels.iter().copied())
        .collect();
    let has_visc = mats.iter().any(|m| m.mech_active && m.eta.iter().flatten().any(|&v| v != 0.0));
    let has_damping = has_visc || !lrb_labels.is_empty();
    let has_dielectric_loss = mats.iter().any(|m| m.elec_active && m.eps_im.iter().flatten().any(|&v| v.abs() > 0.0));
    let has_periodic = !spec.periodic.is_empty();
    let is_hermitian = !any_pml && !has_dielectric_loss && !has_damping;
    let is_complex = any_pml || has_dielectric_loss || has_damping || has_periodic;

    let mech_untouched = mixed_untouched_dofs(mixed, dof_map, d, &spec.inactive);
    let outer_clamped =
        if any_pml { outer_clamped_dofs_mixed(mixed, dof_map, spec, &resolved_pml, false, 0.0) } else { Vec::new() };

    let mut periodic_rels = Vec::new();
    for p in &spec.periodic {
        periodic_rels.extend(mixed_nodal_relations(mixed, dof_map, d + 1, p)?);
    }

    let get_cell = |ci: usize, out: &mut Vec<usize>| {
        out.clear();
        out.extend_from_slice(dof_map.cell_dofs(ci));
    };
    let get_bd = |labels: &[i32]| dof_map.boundary_dofs(mixed, labels);
    let get_cd = |labels: &[i32]| collect_conductor_domain_dofs_mixed(mixed, dof_map, labels);
    let constraints = build_piezo_constraints(
        d,
        n1,
        spec,
        &mixed.cell_labels,
        &get_cell,
        &get_bd,
        &get_cd,
        mech_untouched,
        outer_clamped,
        Vec::new(),
        periodic_rels,
    )?;

    if n_modes < 1 || n_modes + 1 >= constraints.n_u_red {
        return Err("number of modes must be between 1 and the number of free mechanical dofs".into());
    }

    let diagnostics = resolved_pml.diagnostics.clone();
    let n_u_red = constraints.n_u_red;
    let n_tot_red = constraints.red_all.n_reduced();

    let (
        k_full_real,
        k_red_real,
        m_uu_red_real,
        k_full_c,
        k_red_c,
        m_uu_red_c,
        c_uu_red_c,
        k_pp_red_real,
        k_pp_red_c,
        k_pu_red_real,
        k_pu_red_c,
    ) = if !is_complex {
        let kf = assemble_piezo_real_mixed(
            mixed,
            dof_map,
            &mats,
            spec.plane,
            spec.thickness,
            &ResolvedIeDomains::default(),
            &spec.inactive,
            &constraints.conductor_domains,
        );
        let kr = congruence_real(&constraints.red_all, &kf)?;
        let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let m1 = assemble_mass_mixed(mixed, dof_map, &rho_t);
        let blocks: Vec<&RCsr> = (0..d).map(|_| &m1).collect();
        let m_uu = RCsr::block_diag(&blocks);
        let mr = congruence_real(&constraints.red_u, &m_uu)?;
        let k_pp_r = sub_csr(&kr, n_u_red, n_tot_red, n_u_red, n_tot_red);
        let k_pu_r = sub_csr(&kr, n_u_red, n_tot_red, 0, n_u_red);
        (Some(kf), Some(kr), Some(mr), None, None, None, None, Some(k_pp_r), None, Some(k_pu_r), None)
    } else {
        let kf_c = assemble_piezo_complex_mixed(mixed, dof_map, &mats, spec.plane, spec.thickness, &resolved_pml);
        let kr_c = congruence_fast(&constraints.red_all, &kf_c);
        let m_uu_c = if any_pml {
            assemble_pml_mass_mixed(mixed, dof_map, &rho, spec.thickness, &resolved_pml)
        } else {
            let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
            let m1 = assemble_mass_mixed(mixed, dof_map, &rho_t);
            let blocks: Vec<&RCsr> = (0..d).map(|_| &m1).collect();
            RCsr::block_diag(&blocks).to_complex()
        };
        let mr_c = congruence_fast(&constraints.red_u, &m_uu_c);

        let c_uu_red_c = if has_damping {
            let mut c_uu_r = RCsr::zeros(constraints.n_mech, constraints.n_mech);
            if has_visc {
                let eta_wafer: Vec<_> = mats.iter().map(|m| m.eta).collect();
                let cv = crate::solid::assemble_solid_mixed(mixed, dof_map, &eta_wafer, spec.plane, spec.thickness);
                c_uu_r = RCsr::axpby(1.0, &c_uu_r, 1.0, &cv);
            }
            if !lrb_labels.is_empty() {
                let c6_arr: Vec<_> = mats.iter().map(|m| m.c6).collect();
                let cl = damping_matrix_mixed(mixed, dof_map, &c6_arr, &rho, &lrb_labels, d, thk);
                c_uu_r = RCsr::axpby(1.0, &c_uu_r, 1.0, &cl);
            }
            Some(congruence_fast(&constraints.red_u, &c_uu_r.to_complex()))
        } else {
            None
        };

        let k_pp_r = sub_csr(&kr_c, n_u_red, n_tot_red, n_u_red, n_tot_red);
        let k_pu_r = sub_csr(&kr_c, n_u_red, n_tot_red, 0, n_u_red);
        (None, None, None, Some(kf_c), Some(kr_c), Some(mr_c), c_uu_red_c, None, Some(k_pp_r), None, Some(k_pu_r))
    };

    let pml_cells: Vec<bool> = (0..mixed.n_cells()).map(|ci| resolved_pml.contains(mixed.cell_labels[ci])).collect();
    let has_interior = pml_cells.iter().zip(&rho).any(|(&p, &r)| !p && r > 0.0);
    let (m1_phys_all, m1_phys_pml) = if any_pml {
        let rho_all_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let rho_pml_t: Vec<f64> = rho.iter().zip(&pml_cells).map(|(&r, &p)| if p { r * thk } else { 0.0 }).collect();
        (Some(assemble_mass_mixed(mixed, dof_map, &rho_all_t)), Some(assemble_mass_mixed(mixed, dof_map, &rho_pml_t)))
    } else {
        (None, None)
    };

    Ok(PiezoEigenAssembly {
        constraints,
        mats,
        by_dom,
        k_full_real,
        k_red_real,
        m_uu_red_real,
        k_full_c,
        k_red_c,
        m_uu_red_c,
        c_uu_red_c,
        k_pp_red_real,
        k_pp_red_c,
        k_pu_red_real,
        k_pu_red_c,
        m1_phys_all,
        m1_phys_pml,
        resolved_pml,
        is_complex,
        is_hermitian,
        has_damping,
        has_pml: any_pml,
        has_interior,
        n_modes,
        f_search,
        n_pairs: spec.periodic.len(),
        diagnostics,
    })
}

#[allow(dead_code)]
fn _touch_eval(_: impl Fn()) {
    let _ = eval_facet_jacobian;
    let _ = phase_factor;
    let _ = quad_rule;
    let _ = shape;
}

/// Unreduced full complex operators for Stage 1 of the coupled piezoelectric driven frequency solver.
pub(crate) struct PiezoDrivenFull {
    pub constraints: PiezoConstraints,
    pub k_full: CCsr,
    pub m_uu_full: CCsr,
    pub c_uu_full: Option<CCsr>,
    pub resolved_pml: ResolvedPmlDomains,
    pub m1_phys_all: Option<RCsr>,
    pub m1_phys_pml: Option<RCsr>,
    pub is_complex: bool,
    pub has_pml: bool,
    pub has_damping: bool,
    pub has_interior: bool,
    pub diagnostics: Vec<String>,
}

/// Assembles full coupled operators (`K`, `M_uu`, `C_uu`) and terminal constraints for `Study::Driven` on [`SimplexMesh`].
pub(crate) fn assemble_piezo_driven_simplex(mesh: &SimplexMesh, spec: &Problem) -> Result<PiezoDrivenFull, String> {
    crate::problem::check_simplex_straight(mesh)?;
    let f_ref = spec.study.reference_freq();
    check_labels_simplex(mesh, spec)?;

    let d = mesh.dim;
    let order = spec.order;
    let n1 = mesh.n_dofs(order);
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let (_by_dom, mats) = resolve_piezo_mats_simplex(mesh, spec, true)?;
    let rho: Vec<f64> = mats.iter().map(|m| m.rho).collect();
    let resolved_pml = ResolvedPmlDomains::from_simplex(mesh, spec, false, f_ref)?;
    let any_pml = !resolved_pml.domains.is_empty();
    let lrb_labels: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| matches!(b.kind, BcKind::LowReflecting))
        .flat_map(|b| b.labels.iter().copied())
        .collect();
    let has_visc = mats.iter().any(|m| m.mech_active && m.eta.iter().flatten().any(|&v| v != 0.0));
    let has_damping = has_visc || !lrb_labels.is_empty();
    let has_dielectric_loss = mats.iter().any(|m| m.elec_active && m.eps_im.iter().flatten().any(|&v| v.abs() > 0.0));
    let has_periodic = !spec.periodic.is_empty();
    let is_complex = any_pml || has_dielectric_loss || has_damping || has_periodic;

    let mech_untouched = untouched_dofs(mesh, order, d, &spec.inactive);
    let outer_clamped =
        if any_pml { outer_clamped_dofs_simplex(mesh, order, spec, &resolved_pml, false, 0.0) } else { Vec::new() };

    let mut periodic_rels = Vec::new();
    for p in &spec.periodic {
        let mt = match_vertices(mesh, p)?;
        periodic_rels.extend(nodal_relations(mesh, order, d + 1, &mt, p.phase)?);
    }

    let get_cell = |ci: usize, out: &mut Vec<usize>| {
        mesh.cell_dofs(ci, order, out);
    };
    let get_bd = |labels: &[i32]| boundary_dofs(mesh, labels, order);
    let get_cd = |labels: &[i32]| collect_conductor_domain_dofs_simplex(mesh, order, labels);
    let constraints = build_piezo_constraints(
        d,
        n1,
        spec,
        &mesh.cell_labels,
        &get_cell,
        &get_bd,
        &get_cd,
        mech_untouched,
        outer_clamped,
        Vec::new(),
        periodic_rels,
    )?;

    if constraints.n_u_red == 0 && constraints.n_phi_red == 0 {
        return Err("driven study has no free degrees of freedom after constraints".into());
    }

    let k_full = if any_pml || has_dielectric_loss {
        assemble_piezo_complex_simplex(mesh, &mats, order, spec.plane, spec.thickness, &resolved_pml)
    } else {
        assemble_piezo_real_simplex(mesh, &mats, order, spec.plane, spec.thickness).to_complex()
    };

    let m_uu_full = if any_pml {
        assemble_pml_mass_simplex(mesh, &rho, order, spec.thickness, &resolved_pml)
    } else {
        let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let m1 = assemble_mass(mesh, order, &rho_t);
        let blocks: Vec<&RCsr> = (0..d).map(|_| &m1).collect();
        RCsr::block_diag(&blocks).to_complex()
    };

    let c_uu_full = if has_damping {
        let mut c_uu_r = RCsr::zeros(constraints.n_mech, constraints.n_mech);
        if has_visc {
            let eta_wafer: Vec<_> = mats.iter().map(|m| m.eta).collect();
            let cv = crate::solid::assemble_solid(mesh, &eta_wafer, order, spec.plane, spec.thickness);
            c_uu_r = RCsr::axpby(1.0, &c_uu_r, 1.0, &cv);
        }
        if !lrb_labels.is_empty() {
            let c6_arr: Vec<_> = mats.iter().map(|m| m.c6).collect();
            let cl = damping_matrix(mesh, &c6_arr, &rho, &lrb_labels, order, d);
            let cl_sc = if d == 2 { thk } else { 1.0 };
            c_uu_r = RCsr::axpby(1.0, &c_uu_r, cl_sc, &cl);
        }
        Some(c_uu_r.to_complex())
    } else {
        None
    };

    let pml_cells: Vec<bool> = (0..mesh.n_cells()).map(|ci| resolved_pml.contains(mesh.cell_labels[ci])).collect();
    let has_interior = pml_cells.iter().zip(&rho).any(|(&p, &r)| !p && r > 0.0);
    let (m1_phys_all, m1_phys_pml) = if any_pml {
        let rho_all_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let rho_pml_t: Vec<f64> = rho.iter().zip(&pml_cells).map(|(&r, &p)| if p { r * thk } else { 0.0 }).collect();
        (Some(assemble_mass(mesh, order, &rho_all_t)), Some(assemble_mass(mesh, order, &rho_pml_t)))
    } else {
        (None, None)
    };

    let diagnostics = resolved_pml.diagnostics.clone();
    Ok(PiezoDrivenFull {
        constraints,
        k_full,
        m_uu_full,
        c_uu_full,
        resolved_pml,
        m1_phys_all,
        m1_phys_pml,
        is_complex,
        has_pml: any_pml,
        has_damping,
        has_interior,
        diagnostics,
    })
}

/// Assembles full coupled operators (`K`, `M_uu`, `C_uu`) and terminal constraints for `Study::Driven` natively on [`MixedDiscretisation`].
pub(crate) fn assemble_piezo_driven_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
) -> Result<PiezoDrivenFull, String> {
    let f_ref = spec.study.reference_freq();
    mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let (_by_dom, mats) = resolve_piezo_mats_mixed(mixed, spec, true)?;
    let rho: Vec<f64> = mats.iter().map(|m| m.rho).collect();
    let resolved_pml = ResolvedPmlDomains::from_mixed(mixed, spec, false, f_ref)?;
    let any_pml = !resolved_pml.domains.is_empty();
    let lrb_labels: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| matches!(b.kind, BcKind::LowReflecting))
        .flat_map(|b| b.labels.iter().copied())
        .collect();
    let has_visc = mats.iter().any(|m| m.mech_active && m.eta.iter().flatten().any(|&v| v != 0.0));
    let has_damping = has_visc || !lrb_labels.is_empty();
    let has_dielectric_loss = mats.iter().any(|m| m.elec_active && m.eps_im.iter().flatten().any(|&v| v.abs() > 0.0));
    let has_periodic = !spec.periodic.is_empty();
    let is_complex = any_pml || has_dielectric_loss || has_damping || has_periodic;

    let mech_untouched = mixed_untouched_dofs(mixed, dof_map, d, &spec.inactive);
    let outer_clamped =
        if any_pml { outer_clamped_dofs_mixed(mixed, dof_map, spec, &resolved_pml, false, 0.0) } else { Vec::new() };

    let mut periodic_rels = Vec::new();
    for p in &spec.periodic {
        periodic_rels.extend(mixed_nodal_relations(mixed, dof_map, d + 1, p)?);
    }

    let get_cell = |ci: usize, out: &mut Vec<usize>| {
        out.clear();
        out.extend_from_slice(dof_map.cell_dofs(ci));
    };
    let get_bd = |labels: &[i32]| dof_map.boundary_dofs(mixed, labels);
    let get_cd = |labels: &[i32]| collect_conductor_domain_dofs_mixed(mixed, dof_map, labels);
    let constraints = build_piezo_constraints(
        d,
        n1,
        spec,
        &mixed.cell_labels,
        &get_cell,
        &get_bd,
        &get_cd,
        mech_untouched,
        outer_clamped,
        Vec::new(),
        periodic_rels,
    )?;

    if constraints.n_u_red == 0 && constraints.n_phi_red == 0 {
        return Err("driven study has no free degrees of freedom after constraints".into());
    }

    let k_full = if any_pml || has_dielectric_loss {
        assemble_piezo_complex_mixed(mixed, dof_map, &mats, spec.plane, spec.thickness, &resolved_pml)
    } else {
        assemble_piezo_real_mixed(
            mixed,
            dof_map,
            &mats,
            spec.plane,
            spec.thickness,
            &ResolvedIeDomains::default(),
            &spec.inactive,
            &constraints.conductor_domains,
        )
        .to_complex()
    };

    let m_uu_full = if any_pml {
        assemble_pml_mass_mixed(mixed, dof_map, &rho, spec.thickness, &resolved_pml)
    } else {
        let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let m1 = assemble_mass_mixed(mixed, dof_map, &rho_t);
        let blocks: Vec<&RCsr> = (0..d).map(|_| &m1).collect();
        RCsr::block_diag(&blocks).to_complex()
    };

    let c_uu_full = if has_damping {
        let mut c_uu_r = RCsr::zeros(constraints.n_mech, constraints.n_mech);
        if has_visc {
            let eta_wafer: Vec<_> = mats.iter().map(|m| m.eta).collect();
            let cv = crate::solid::assemble_solid_mixed(mixed, dof_map, &eta_wafer, spec.plane, spec.thickness);
            c_uu_r = RCsr::axpby(1.0, &c_uu_r, 1.0, &cv);
        }
        if !lrb_labels.is_empty() {
            let c6_arr: Vec<_> = mats.iter().map(|m| m.c6).collect();
            let cl = damping_matrix_mixed(mixed, dof_map, &c6_arr, &rho, &lrb_labels, d, thk);
            c_uu_r = RCsr::axpby(1.0, &c_uu_r, 1.0, &cl);
        }
        Some(c_uu_r.to_complex())
    } else {
        None
    };

    let pml_cells: Vec<bool> = (0..mixed.n_cells()).map(|ci| resolved_pml.contains(mixed.cell_labels[ci])).collect();
    let has_interior = pml_cells.iter().zip(&rho).any(|(&p, &r)| !p && r > 0.0);
    let (m1_phys_all, m1_phys_pml) = if any_pml {
        let rho_all_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
        let rho_pml_t: Vec<f64> = rho.iter().zip(&pml_cells).map(|(&r, &p)| if p { r * thk } else { 0.0 }).collect();
        (Some(assemble_mass_mixed(mixed, dof_map, &rho_all_t)), Some(assemble_mass_mixed(mixed, dof_map, &rho_pml_t)))
    } else {
        (None, None)
    };

    let diagnostics = resolved_pml.diagnostics.clone();
    Ok(PiezoDrivenFull {
        constraints,
        k_full,
        m_uu_full,
        c_uu_full,
        resolved_pml,
        m1_phys_all,
        m1_phys_pml,
        is_complex,
        has_pml: any_pml,
        has_damping,
        has_interior,
        diagnostics,
    })
}
