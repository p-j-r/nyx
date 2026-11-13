//! The visualisation layer of a results file and its XDMF description (ADR 62), recovered from the writer of the old `femcore::h5io` and
//! moved to the new result types. `/viz/<physics>` holds the nodes (dofs) and cells in the node order XDMF reads; an XDMF file (XML, written
//! with `quick-xml`) next to the `.h5` names those datasets, so ParaView and VisIt open the results without copying the data.

use crate::h5io::{link_name, Io};
use femcore::mesh::SimplexMesh;
use nyx_physics::SolveResult;
use std::path::Path;

/// Name usable in a plain HDF5/XDMF path: letters, digits and underscores.
pub fn sanitize(name: &str) -> String {
    let s: String = name.chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
    if s.chars().next().is_none_or(|c| c.is_ascii_digit()) {
        format!("f_{s}")
    } else {
        s
    }
}

/// XDMF topology type and the local node order of the cells for a mesh of the given order. Quadratic cells use the VTK/XDMF numbering
/// (corners, then edge midpoints), so the nodes are the dofs of the P2 field and the attributes are nodal.
fn topology(dim: usize, order: usize) -> (&'static str, Vec<usize>) {
    match (dim, order) {
        (2, 1) => ("Triangle", vec![0, 1, 2]),
        (3, 1) => ("Tetrahedron", vec![0, 1, 2, 3]),
        // ours: v0 v1 v2 e01 e02 e12  ->  v0 v1 v2 e01 e12 e02
        (2, _) => ("Triangle_6", vec![0, 1, 2, 3, 5, 4]),
        // ours: v0..v3 e01 e02 e03 e12 e13 e23  ->  e01 e12 e02 e03 e13 e23
        (_, _) => ("Tetrahedron_10", vec![0, 1, 2, 3, 4, 7, 5, 6, 8, 9]),
    }
}

/// `/viz/<physics>`: the nodes and cells in the node order XDMF reads, and for each solution the fields it points at: a scalar field is a
/// hard link to the dataset under `/results`, a vector field gets a three-component copy (`v3_<name>`), since XDMF vectors have three.
pub(crate) fn write_viz(io: &Io, code: &str, mesh: &SimplexMesh, r: &SolveResult) -> Result<(), String> {
    let order = if r.order == 2 { 2 } else { 1 };
    let viz = match io.file.group("viz") {
        Ok(g) => g,
        Err(_) => io.group(&io.file, "viz")?,
    };
    let g = io.group(&viz, &link_name(code))?;
    let (topo, perm) = topology(mesh.dim, order);
    let mut conn: Vec<u32> = Vec::with_capacity(mesh.n_cells() * perm.len());
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        conn.extend(perm.iter().map(|&i| dofs[i] as u32));
    }
    io.put(&g, "points", &mesh.dof_coordinates(order), mesh.dim)?;
    io.put(&g, "cells", &conn, perm.len())?;
    io.attr_text(&g, "topology", topo)?;
    for (k, s) in r.solutions.iter().enumerate() {
        let sg = io.group(&g, &format!("sol_{k}"))?;
        for (i, f) in s.fields.iter().enumerate() {
            if f.components.is_empty() {
                let src = format!("/results/{}/solution_{k}/fields/{i:02}_{}", link_name(code), link_name(&f.name));
                io.link(&sg, &src, &sanitize(&f.name))?;
            } else {
                let m = f.components.len();
                let n = f.values.len() / m.max(1);
                let mut v3 = vec![0.0; n * 3];
                for (j, c) in f.components.iter().enumerate() {
                    // wafer axes in the order of the coordinates of the mesh: a 2D section is (x, z), a 3D mesh (x, y, z)
                    let slot = match (c.as_str(), mesh.dim) {
                        ("x", _) => 0,
                        ("y", _) => 1,
                        ("z", 2) => 1,
                        _ => 2,
                    };
                    for p in 0..n {
                        v3[p * 3 + slot] = f.values[p * m + j];
                    }
                }
                io.put(&sg, &format!("v3_{}", sanitize(&f.name)), &v3, 3)?;
            }
        }
    }
    Ok(())
}

/// Writes the XDMF description of the results in `h5_path` next to it (same name, extension `.xdmf`). Every physics is a collection of
/// grids, one per solution (time step, mode).
pub fn write_xdmf(
    h5_path: &Path,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<std::path::PathBuf, String> {
    use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
    use quick_xml::Writer;
    let h5name = h5_path.file_name().and_then(|n| n.to_str()).ok_or("bad file name")?.to_string();
    let xe = |x: std::io::Error| x.to_string();
    let mut w = Writer::new_with_indent(Vec::new(), b' ', 2);
    w.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None))).map_err(xe)?;
    w.write_event(Event::DocType(BytesText::from_escaped("Xdmf SYSTEM \"Xdmf.dtd\" []"))).map_err(xe)?;
    let mut x = BytesStart::new("Xdmf");
    x.push_attribute(("Version", "3.0"));
    w.write_event(Event::Start(x)).map_err(xe)?;
    w.write_event(Event::Start(BytesStart::new("Domain"))).map_err(xe)?;
    let item = |w: &mut Writer<Vec<u8>>, dims: &str, number: &str, prec: &str, path: &str| -> Result<(), String> {
        let mut d = BytesStart::new("DataItem");
        d.push_attribute(("Dimensions", dims));
        d.push_attribute(("NumberType", number));
        d.push_attribute(("Precision", prec));
        d.push_attribute(("Format", "HDF"));
        w.write_event(Event::Start(d)).map_err(|x| x.to_string())?;
        w.write_event(Event::Text(BytesText::new(&format!("{h5name}:{path}")))).map_err(|x| x.to_string())?;
        w.write_event(Event::End(BytesEnd::new("DataItem"))).map_err(|x| x.to_string())
    };
    for (code, result) in results {
        let order = if result.order == 2 { 2 } else { 1 };
        let (topo, perm) = topology(mesh.dim, order);
        let np = mesh.n_dofs(order);
        let nc = mesh.n_cells();
        let base = format!("/viz/{}", link_name(code));
        // time steps when a transient study orders its solutions by a value that rises; modes are numbered
        let times_ok = result.solutions.len() > 1
            && result.study.contains("transient")
            && result.solutions.iter().all(|s| s.value.is_some())
            && result.solutions.windows(2).all(|p| p[1].value >= p[0].value);
        let collection = result.solutions.len() > 1;
        if collection {
            let mut c = BytesStart::new("Grid");
            c.push_attribute(("Name", code.as_str()));
            c.push_attribute(("GridType", "Collection"));
            c.push_attribute(("CollectionType", "Temporal"));
            w.write_event(Event::Start(c)).map_err(xe)?;
        }
        for (k, s) in result.solutions.iter().enumerate() {
            let mut g = BytesStart::new("Grid");
            g.push_attribute(("Name", format!("{code}_solution_{k}").as_str()));
            g.push_attribute(("GridType", "Uniform"));
            w.write_event(Event::Start(g)).map_err(xe)?;
            if collection {
                let mut t = BytesStart::new("Time");
                let tv = if times_ok { s.value.unwrap_or(k as f64) } else { k as f64 };
                t.push_attribute(("Value", format!("{tv:e}").as_str()));
                w.write_event(Event::Empty(t)).map_err(xe)?;
            }
            let mut t = BytesStart::new("Topology");
            t.push_attribute(("TopologyType", topo));
            t.push_attribute(("NumberOfElements", nc.to_string().as_str()));
            w.write_event(Event::Start(t)).map_err(xe)?;
            item(&mut w, &format!("{nc} {}", perm.len()), "UInt", "4", &format!("{base}/cells"))?;
            w.write_event(Event::End(BytesEnd::new("Topology"))).map_err(xe)?;
            let mut ge = BytesStart::new("Geometry");
            ge.push_attribute(("GeometryType", if mesh.dim == 2 { "XY" } else { "XYZ" }));
            w.write_event(Event::Start(ge)).map_err(xe)?;
            item(&mut w, &format!("{np} {}", mesh.dim), "Float", "8", &format!("{base}/points"))?;
            w.write_event(Event::End(BytesEnd::new("Geometry"))).map_err(xe)?;
            for f in &s.fields {
                let vector = !f.components.is_empty();
                let mut a = BytesStart::new("Attribute");
                a.push_attribute(("Name", f.name.as_str()));
                a.push_attribute(("AttributeType", if vector { "Vector" } else { "Scalar" }));
                a.push_attribute(("Center", "Node"));
                w.write_event(Event::Start(a)).map_err(xe)?;
                if vector {
                    item(
                        &mut w,
                        &format!("{np} 3"),
                        "Float",
                        "8",
                        &format!("{base}/sol_{k}/v3_{}", sanitize(&f.name)),
                    )?;
                } else {
                    item(&mut w, &np.to_string(), "Float", "8", &format!("{base}/sol_{k}/{}", sanitize(&f.name)))?;
                }
                w.write_event(Event::End(BytesEnd::new("Attribute"))).map_err(xe)?;
            }
            w.write_event(Event::End(BytesEnd::new("Grid"))).map_err(xe)?;
        }
        if collection {
            w.write_event(Event::End(BytesEnd::new("Grid"))).map_err(xe)?;
        }
    }
    w.write_event(Event::End(BytesEnd::new("Domain"))).map_err(xe)?;
    w.write_event(Event::End(BytesEnd::new("Xdmf"))).map_err(xe)?;
    let xml = String::from_utf8(w.into_inner()).map_err(|x| x.to_string())?;
    let out = h5_path.with_extension("xdmf");
    std::fs::write(&out, xml).map_err(|x| format!("{}: {x}", out.display()))?;
    Ok(out)
}
