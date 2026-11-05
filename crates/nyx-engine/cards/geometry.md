Geometry is a feature list: shapes, then operations on them, then one overlay that names the regions.

- A 2D model is the x-z cross-section (x lateral, z up). add_rect takes the two corners {min:[x,z], max:[x,z]}; add_circle {center, r}; add_polygon {pts:[[x,z],...]}. list_operations gives every operation with its fields, defaults and a one-line description.
- Features have names (give "name" or accept rect1, rect2). Later features refer to earlier ones by name: add_difference {a:"outer", b:["hole"]}, add_union {ids:[...]}, add_intersect {ids:[...]}.
- add_overlay is the root. regions are in painter order, each {name, input:<feature name>}; where regions overlap the later one wins. So list the big background first (a buried oxide, a cladding) and the core last. Regions named in the overlay are what materials are assigned to.
- Parameters: set_param {name, value} defines a number or an expression other fields can use (width, 2*thickness). update_feature changes fields of a feature; remove_feature, rename_feature, suppress_feature.
- Boundaries get names from the regions they separate ("BOX|Si", the interface), or from a feature edge (box.bottom). geometry_names lists them. A selector is also accepted anywhere a boundary is wanted: "outer & plane(x = 0)".
- Absorbing layers: add_pml around a window (thickness, element layers, grading); the pieces are named region.pml_* and region* in a material assignment covers the region and its pieces.
- Mesh size hints: add_size {sel, h} asks for element size h at the selected boundaries or regions. The global size is set_mesh h.
- A 3D model has dim 3 (project_new {dim:3}); a 2D section can also be swept along y with set_mesh {sweep:[{thickness, n}]}, the two ends are then the boundaries ymin and ymax.
- Check before meshing: check_project, geometry_names. Errors name the feature and the reason.
