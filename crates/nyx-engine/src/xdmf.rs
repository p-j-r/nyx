// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The XDMF description of a results file (ADR 62, ADR 99): a small XML file next to the `.h5` that is the light index of the heavy data in
//! it. It names the datasets of the file as they are stored, per physics `points` and `cells` and per solution `fields/<name>`, with their
//! shapes, so a reader finds a field and reads the part of it that it needs without opening the rest. Nothing is copied into the HDF5 file for
//! it: the XDMF refers to the nodes (`points`), the node of every cell node (`cells`, in the cell node order of nyx, which the `Information`
//! element `CellNodeOrder` of the topology states) and the stored fields. A reader that assumes the VTK order of quadratic cells would read
//! the mid-edge nodes of a quadratic cell in the wrong order; the file is made for nyx and its tools.

use crate::h5io::{cell_order_text, link_name};
use femcore::mesh::SimplexMesh;
use nyx_physics::SolveResult;
use std::path::Path;

/// Written at the top of every XDMF file: who the file is for.
const NOTE: &str = "Index of a nyx results file, made for nyx and its tools only. Do not open it with ParaView, VisIt or another XDMF viewer: the cells of a quadratic mesh are listed in the cell node order of nyx (see CellNodeOrder), not in the order those programs assume, so they would draw them wrongly.";

/// The XDMF topology type of the cells of a mesh of the given order (the node order within a cell is that of `cell_order_text`).
fn topology(dim: usize, order: usize) -> &'static str {
    match (dim, order) {
        (2, 1) => "Triangle",
        (3, 1) => "Tetrahedron",
        (2, _) => "Triangle_6",
        (_, _) => "Tetrahedron_10",
    }
}

/// Writes the XDMF description of the results in `h5_path` next to it (same name, extension `.xdmf`). Every physics is a collection of
/// grids, one per solution (time step, mode).
pub fn write_xdmf(
    h5_path: &Path,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<std::path::PathBuf, String> {
    if mesh.n_cells() == 0 || mesh.points.is_empty() {
        return Err("cannot write XDMF index for a results file with no saved mesh/fields".into());
    }
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
    w.write_event(Event::Comment(BytesText::new(NOTE))).map_err(xe)?;
    w.write_event(Event::Start(BytesStart::new("Domain"))).map_err(xe)?;
    {
        let mut i = BytesStart::new("Information");
        i.push_attribute(("Name", "Note"));
        i.push_attribute(("Value", NOTE));
        w.write_event(Event::Empty(i)).map_err(xe)?;
        if mesh.has_curved() {
            let mut c = BytesStart::new("Information");
            c.push_attribute(("Name", "CurvedEdgesIgnored"));
            c.push_attribute((
                "Value",
                "mesh contains curved mid-edge geometry nodes that XDMF linear/quadratic cell topology indexes without curved surface projection",
            ));
            w.write_event(Event::Empty(c)).map_err(xe)?;
        }
    }
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
    let info = |w: &mut Writer<Vec<u8>>, name: &str, value: &str| -> Result<(), String> {
        let mut i = BytesStart::new("Information");
        i.push_attribute(("Name", name));
        i.push_attribute(("Value", value));
        w.write_event(Event::Empty(i)).map_err(|x| x.to_string())
    };
    for (code, result) in results {
        let order = if result.order == 2 { 2 } else { 1 };
        let topo = topology(mesh.dim, order);
        let nloc = mesh.nloc(order);
        let np = mesh.n_dofs(order);
        let nc = mesh.n_cells();
        let base = format!("/results/{}", link_name(code));
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
            info(&mut w, "CellNodeOrder", cell_order_text(mesh.dim, order))?;
            item(&mut w, &format!("{nc} {nloc}"), "UInt", "4", &format!("{base}/cells"))?;
            w.write_event(Event::End(BytesEnd::new("Topology"))).map_err(xe)?;
            let mut ge = BytesStart::new("Geometry");
            ge.push_attribute(("GeometryType", if mesh.dim == 2 { "XY" } else { "XYZ" }));
            w.write_event(Event::Start(ge)).map_err(xe)?;
            item(&mut w, &format!("{np} {}", mesh.dim), "Float", "8", &format!("{base}/points"))?;
            w.write_event(Event::End(BytesEnd::new("Geometry"))).map_err(xe)?;
            for (i, f) in s.fields.iter().enumerate() {
                let m = f.components.len();
                let mut a = BytesStart::new("Attribute");
                a.push_attribute(("Name", f.name.as_str()));
                a.push_attribute((
                    "AttributeType",
                    match m {
                        0 => "Scalar",
                        3 => "Vector",
                        _ => "Matrix",
                    },
                ));
                a.push_attribute(("Center", "Node"));
                w.write_event(Event::Start(a)).map_err(xe)?;
                if m > 0 {
                    // the components in the order of the stored values, which are interleaved by node
                    info(&mut w, "Components", &f.components.join(" "))?;
                }
                let dims = if m == 0 { np.to_string() } else { format!("{np} {m}") };
                let src = format!("{base}/solution_{k}/fields/{i:02}_{}", link_name(&f.name));
                item(&mut w, &dims, "Float", "8", &src)?;
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
