# Addendum 1 to GEOM_HANDOFF.md: support for parameter sweeps

Read after `GEOM_HANDOFF.md`. Nothing in the original brief changes; this adds five small requirements so that the integration side can reuse meshes
and solver setup across parameter sweeps (it will morph a mesh instead of remeshing when the topology is unchanged). Where this conflicts with the brief, this wins.
Fit them into the existing steps: items 1, 2, 4 and 5 belong with step B, item 3 with the association contract in section 5 (build it in step A for 2D).

## Requirements

1. **`topology_signature() -> u64`** on the evaluated model. A hash of the combinatorial structure only: regions, faces, edges, vertices, loops, incidence,
   surface and curve *kinds*, and region kinds. It must **not** depend on parameter values. Two evaluations of the same feature list with different parameters
   have equal signatures exactly when their topology is the same (a feature disappearing, a face splitting, or an intersection appearing changes it).

2. **Stable ids within a signature.** If two evaluations have equal signatures, they have identical sets of `RegionId`, `FaceId`, `EdgeId`, `VertexId`,
   in identical order, and each id means the same entity. (The brief already requires determinism and stability under *unrelated* edits; this extends it to
   *parameter* changes.) Test: sweep a width over a range with unchanged topology and assert id sets and per-entity names are identical at every point.

3. **Parametric positions for boundary nodes.** The mesher side will store, for each node on the boundary, `(entity id, parametric fraction)`: one fraction in
   `[0,1]` along the parameter range of an edge, two in `[0,1]` for a face (normalised to the face's parameter bounding box). Provide on every evaluation:
   - `fractions_of(entity, point) -> Fractions` (the fractions of the closest point of that entity to `point`),
   - `position_at(entity, fractions) -> [f64;3]` (the exact point).
   Round trip within one evaluation must return the point to `1e-12 * L`. Used to move boundary nodes onto the new exact geometry when a parameter changes.

4. **`diff(old: &Evaluated, new: &Evaluated) -> ChangeReport`** with
   `ChangeReport { changed_features: Vec<FeatureId>, signature_equal: bool, kind: None | ZOnly | Morph | Remesh }`.
   `ZOnly` means the 2D overlay of a layer stack is unchanged and only z-extents differ; `Morph` means equal signature with moved geometry; `Remesh` means
   the signature differs. Computed from the content hashes the cache already keeps, so it should cost almost nothing.

5. **Continuous parameterisation.** Edge and face parameters (start point of a circle, seam of a cylinder or torus, direction of a line) must vary
   continuously with the feature parameters: no jump of start point or seam when a radius or width changes slightly, and ranges `[t0, t1]` change continuously.
   Document the convention per curve and surface kind. Test: for a small parameter step, `fractions_of(entity, position_at_old(entity, f))` on the new evaluation moves by a
   distance of the order of the step, never by a seam jump.

## Tests to add
- Signature equality/inequality cases (width sweep with no change; a width large enough to remove a feature; adding an unrelated feature).
- Id stability across a sweep (item 2) and name stability.
- Fraction round trips and continuity (items 3 and 5), including circles (seam), cylinders (theta seam), and the periodic end faces of a Bragg cell.
- `diff` kinds on: a z-only layer thickness change, a width change, a feature removal, no change.

## Not asked of you
Mesh morphing, caching of meshes or solver data, sweep orchestration, UI: all on the integration side (see `SWEEP_CACHING_PROPOSAL.md` if you are curious; you do not need it).
