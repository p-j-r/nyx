// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Direct-path ParaView XDMF 3.0 companion XML generator eliminating HDF5 `/viz` hard links.

use crate::sweep::atomic_write_bytes;
use crate::{MeshData, ModeHeader};
use anyhow::Result;
use std::path::Path;

/// Sanitizes a `nyx` field dataset name like `"00_|E|"` or `"03_Re Ex"` into its legacy `/viz` link name (`"_E_"`, `"Re_Ex"`).
pub fn viz_link_name(field_dataset_name: &str) -> String {
    let raw = match field_dataset_name.find('_') {
        Some(idx) if field_dataset_name[..idx].chars().all(|c| c.is_ascii_digit()) => &field_dataset_name[idx + 1..],
        _ => field_dataset_name,
    };
    raw.chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect()
}

fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            _ => out.push(ch),
        }
    }
    out
}

/// Generates a ParaView-ready XDMF 3.0 companion XML referencing `<h5_filename>` datasets directly,
/// eliminating the need for `/viz` hard links (`SCHEMA_AND_XDMF_SPEC.md`).
///
/// Each entry in `field_datasets` is `(dataset_key, display_name, num_components)` where `num_components` is `1` (Scalar) or `3` (Vector).
pub fn generate_direct_xdmf(
    h5_filename: &str,
    mesh: &MeshData,
    physics: &str,
    headers: &[ModeHeader],
    field_datasets: &[(String, String, usize)],
) -> String {
    let n_cells = mesh.cells_shape[0];
    let n_verts = mesh.cells_shape[1];
    let n_pts = mesh.points_shape[0];
    let dim = mesh.points_shape[1];
    let geom_type = if dim == 2 { "XY" } else { "XYZ" };
    let h5_esc = xml_escape(h5_filename);
    let topo_esc = xml_escape(&mesh.cell_topology);
    let phys_esc = xml_escape(physics);

    let mut xml = String::with_capacity(8192);
    xml.push_str("<?xml version=\"1.0\"?>\n");
    xml.push_str("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    xml.push_str("<Xdmf Version=\"3.0\">\n");
    xml.push_str("  <Domain>\n");

    // 1. Base mesh grid
    xml.push_str("    <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    xml.push_str(&format!("      <Topology TopologyType=\"{}\" NumberOfElements=\"{}\">\n", topo_esc, n_cells));
    xml.push_str(&format!(
        "        <DataItem Dimensions=\"{} {}\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">{}:/mesh/cells</DataItem>\n",
        n_cells, n_verts, h5_esc
    ));
    xml.push_str("      </Topology>\n");
    xml.push_str(&format!("      <Geometry GeometryType=\"{}\">\n", geom_type));
    xml.push_str(&format!(
        "        <DataItem Dimensions=\"{} {}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{}:/mesh/points</DataItem>\n",
        n_pts, dim, h5_esc
    ));
    xml.push_str("      </Geometry>\n");
    xml.push_str("      <Attribute Name=\"cell_labels\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    xml.push_str(&format!(
        "        <DataItem Dimensions=\"{}\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">{}:/mesh/cell_labels</DataItem>\n",
        n_cells, h5_esc
    ));
    xml.push_str("      </Attribute>\n");
    xml.push_str("    </Grid>\n");

    // 2. Per-mode solution grids
    if !headers.is_empty() {
        xml.push_str(&format!(
            "    <Grid Name=\"{}\" GridType=\"Collection\" CollectionType=\"Temporal\">\n",
            phys_esc
        ));
        for h in headers {
            xml.push_str(&format!("      <Grid Name=\"mode_{}\" GridType=\"Uniform\">\n", h.mode_index));
            xml.push_str(&format!("        <Time Value=\"{:.16e}\"/>\n", h.eig_re));
            xml.push_str(&format!(
                "        <Topology TopologyType=\"{}\" NumberOfElements=\"{}\">\n",
                topo_esc, n_cells
            ));
            xml.push_str(&format!(
                "          <DataItem Dimensions=\"{} {}\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">{}:/mesh/cells</DataItem>\n",
                n_cells, n_verts, h5_esc
            ));
            xml.push_str("        </Topology>\n");
            xml.push_str(&format!("        <Geometry GeometryType=\"{}\">\n", geom_type));
            xml.push_str(&format!(
                "          <DataItem Dimensions=\"{} {}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{}:/mesh/points</DataItem>\n",
                n_pts, dim, h5_esc
            ));
            xml.push_str("        </Geometry>\n");

            for (ds_key, display_name, n_comp) in field_datasets {
                let (attr_type, dims) =
                    if *n_comp == 3 { ("Vector", format!("{} 3", n_pts)) } else { ("Scalar", format!("{}", n_pts)) };
                xml.push_str(&format!(
                    "        <Attribute Name=\"{}\" AttributeType=\"{}\" Center=\"Node\">\n",
                    xml_escape(display_name),
                    attr_type
                ));
                xml.push_str(&format!(
                    "          <DataItem Dimensions=\"{}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{}:/results/{}/solution_{}/fields/{}</DataItem>\n",
                    dims,
                    h5_esc,
                    phys_esc,
                    h.mode_index,
                    xml_escape(ds_key)
                ));
                xml.push_str("        </Attribute>\n");
            }
            xml.push_str("      </Grid>\n");
        }
        xml.push_str("    </Grid>\n");
    }

    xml.push_str("  </Domain>\n");
    xml.push_str("</Xdmf>\n");
    xml
}

/// Atomically writes a direct-path ParaView `.xdmf` companion file to `xdmf_path`.
pub fn write_direct_xdmf(
    xdmf_path: &Path,
    h5_filename: &str,
    mesh: &MeshData,
    physics: &str,
    headers: &[ModeHeader],
    field_datasets: &[(String, String, usize)],
) -> Result<()> {
    let xml = generate_direct_xdmf(h5_filename, mesh, physics, headers, field_datasets);
    atomic_write_bytes(xdmf_path, xml.as_bytes())
}
