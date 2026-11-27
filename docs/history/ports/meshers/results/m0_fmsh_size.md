# `.fmsh` storage study

`cargo run --release --example fmsh_size`. "plain" is the same arrays as raw binary (f64 coordinates, u32 indices and
labels, one byte per cell type, f64 fields). The gzip/xz columns run the system tools on the plain arrays for
comparison only (zstd is not installed here); nothing links to them. Every file is decoded again and compared bit for bit.

| case | cells | plain B | fmsh B | ratio | gzip -9 B | xz -9e B | zstd -19 B | enc / dec s |
|---|---|---|---|---|---|---|---|---|
| structured tri 150x150 | 45000 | 1137016 | 37710 | 30.15x | 216502 | 69336 | n/a | 0.09 / 0.01 |
| adapted tri, mesh only | 34709 | 868161 | 410684 | 2.11x | 472245 | 380580 | n/a | 0.13 / 0.03 |
| adapted tri + 3 f64 point fields | 34709 | 1284945 | 769885 | 1.67x | 856549 | 754372 | n/a | 0.31 / 0.07 |
| 3D PML (hex/prism/tet), native | 51648 | 1955640 | 141869 | 13.78x | 354024 | 136136 | n/a | 0.12 / 0.02 |
| 3D PML split to tets | 131712 | 3351288 | 268148 | 12.50x | 582728 | 263612 | n/a | 0.24 / 0.05 |

Reading the table:

* Structured and PML meshes (few distinct coordinates, regular connectivity) shrink 12-30x; xz needs 1.0-1.8x more space on them.
* An adapted unstructured mesh is the hard case: 2.1x. About 219 KB of the 411 KB are coordinates (17k vertices x 16 B raw; only
  the sign/exponent/top-mantissa bytes compress) and 191 KB are connectivity (5.5 B per triangle). Vertex ids of an adapted
  mesh have almost no spatial locality (a 128-entry recent-vertex cache hits 12 % of references), so index coding gains little.
  xz is 7 % smaller here. A topological connectivity coder (edgebreaker style, about 2 bits per triangle plus the id permutation)
  would cut the 191 KB to roughly 40 KB, but it is a much larger piece of code; not done.
* Solution fields are full-entropy doubles: 1.67x including the mesh, i.e. the fields themselves save only the top-byte
  redundancy. Lossy quantisation would cut them much further, but it is not offered since the requirement is no loss.
* Storing the native hex/prism/tet PML mesh takes half the space of the split-to-tets mesh (142 KB vs 268 KB).
