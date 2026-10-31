# References used so far

Status of each citation: **verified** = I know the result well enough to have re-derived or tested
it here; **recalled** = cited from memory, not re-read, treat as unconfirmed.

## Metric algebra (`adapt/metric.rs`)

* Edge length with geometrically varying size, `l = (la - lb) / ln(la / lb)`: **verified**. It is the
  integral of `|e| / h(t)` with `h(t) = ha^(1-t) hb^t`; the test `edge_length_matches_integral`
  checks it numerically. Attributed in the literature to Frey and George, *Mesh Generation* (recalled).
* Intersection by simultaneous reduction: **verified** by the tests (result dominates both inputs,
  determinant equals `det(M1) * prod max(1, lambda_k)`, is symmetric and idempotent). Standard in
  the anisotropic-adaptation literature (Alauzet; Frey and George, recalled).
* Log-Euclidean mean: Arsigny, Fillard, Pennec, Ayache, "Log-Euclidean metrics for fast and simple
  calculus on diffusion tensors", MRM 2006 (recalled).
* Gradation with a linear size-growth law (`h' = h + g * distance` per eigendirection): my
  parameterisation. Hypothesis, not checked against the paper: it matches Alauzet, "Size gradation
  control of anisotropic meshes", FEAD 2010, for an isotropic field with an appropriate `g`.
* L^p optimal metric `M = c det(|H|)^(-1/(2p+D)) |H|` with complexity normalisation
  `C = integral sqrt(det M)`: Loseille and Alauzet, "Continuous mesh framework" parts I and II,
  SIAM J. Numer. Anal. 2011 (recalled; the formula is checked here for the exponent and the
  complexity normalisation by `lp_exponent_matches_formula` and `lp_hits_target`, which test the
  algebra, not the optimality claim).
* Metric volume of the ideal simplex (`sqrt(3)/4`, `sqrt(2)/12`), mean-ratio quality formulas
  `4 sqrt(3) A / sum l^2` and `12 (3V)^(2/3) / sum l^2`: **verified** (equal 1 on the regular
  simplex, `regular_simplices_have_quality_one`).
* Fraction of edges with metric length in `[1/sqrt(2), sqrt(2)]` as a conformity figure: common
  convention in the Mmg/Omega_h documentation (recalled, definition only).

## PML

* Complex coordinate stretching, polynomial profile, design reflection
  `R = exp(-2 sigma_max t / (m + 1))`: Berenger 1994 and the standard textbook derivation
  (recalled). `validate/helmholtz.rs` uses `sigma_max = -(m+1) ln(R0) / (2 t)`.

## Still to read (needed by later milestones, from the handoff)

arXiv 2201.02806 (ParMmg in PETSc); arXiv 1308.2480 (PRAgMaTIc threads); CDT3D
(10.1007/s00366-024-01994-0); Omega_h SCOREC report 2016-25; Goalie docs; Loseille and Alauzet
parts I and II; "Very High Order Anisotropic Metric-Based Mesh Adaptation in 3D"; Dobrzynski and
Frey 2008; Becker and Rannacher DWR; Beck, Hiptmair, Hoppe, Wohlmuth; Demkowicz; Labelle and
Shewchuk 2007. None has been read in this session.
