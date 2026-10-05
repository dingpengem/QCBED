# One stored slice per simulation-cell period

This audit supersedes the earlier adapter policy requiring two slabs. The supported
strict simulation cell now accepts any finite, strictly increasing fractional
boundary list starting at 0 and ending at 1, with at least two values. `(0, 1)`
explicitly selects one slice of thickness `|C|` Å; `(0, .5, 1)` selects two, and
`(0, .25, .5, .75, 1)` selects four. Existing multi-slice inputs are unchanged.

One slice is a slicing approximation: all atoms in the period contribute to one
projected potential before propagation over that period. It is not asserted to
produce the same wave as a finer partition. Slicing and mesh convergence remain
specimen-dependent.

## Source and path audit

| Path | One-slice behavior / retained invariant |
|---|---|
| `MultisliceSimulationCell` | Unchanged exact integer-C-shear proof of transverse lattice translations; P=A−mC, Q=B−nC. Unsupported geometries still fail before slicing. |
| Boundary parser / `build_depth_resolved_sample` | Minimum changed from three endpoints to two. Finite, increasing endpoints exactly spanning [0,1] remain mandatory. |
| Atomic assignment | Recast physical depth z+m*x+n*y; subtract floor(depth), assign half-open intervals using right-sided search. Each atom appears exactly once; boundary z=1 belongs to next periodic copy at z=0. Site order, x/y, occupancy, B, group ownership and recorded C translations are retained. |
| `build_zone_axis_simulation_sample` | Builds one `SliceGeometry(|C|)` and one atom collection. Elastic/absorptive sums and phase-grating routines receive that collection without special cases, scaling or normalization changes. |
| Stored count / repeat handling | One stored type, an explicit sequence of zero-based ID 0 repeated as required; no fixed capacities. |
| `depth_geometry.requested_depth_to_sample_geometry` | Repeated sequential addition while total < requested depth; final whole-slice overshoot and floating-point addition order unchanged. Effective depth is the accumulated sum, not an interpolated depth or a substituted ceil formula. |
| `private predecessor implementation`, `private predecessor implementation` | Actual source MSEQ loops cycle `MOD(I,RNSL)`, replace zero by RNSL, and add DELTAZ×APPLIEDDILATION. RNSL=1 has no zero divisor or missing neighboring slice. No private predecessor implementation-only declaration is used as authority. |
| `simulation`, `multislice`, `tilt_series` | Collections require positive stored count, not two; propagation indexes grating and propagator by sequence ID. complex128, Fourier signs/order, per-slice propagation and final reciprocal sampling are unchanged. |
| Depth refinement | `apply_crystallographic_parameters` rebuilds from the same source partition at changed depth; `refinement_trial._prepare_trial` uses the same cyclic helper. Existing simplex/bounds initialization supports the resulting positive depth and count. |
| Cache/reuse | Prepared snapshots include stored geometries, gratings, parameters and explicit sequence. Modern `evaluate_grouped` rebuilds the physical snapshot per vector; it does not reuse a cached two-slice partition. No new cache or optimization introduced. |
| Resume | Canonical identity hashes derived model/state/layout/run settings, ownership and groups; grouped identity also includes experimental data, bounds and execution settings. One-vs-two-slice partitions yield different identities. Unchanged one-slice Save/reload and checkpoint Resume match; changed partition is rejected. No checkpoint version change needed. |

The exact transverse-periodicity restriction remains in force even for `(0,1)`.
This change does not introduce a fallback that projects unsupported skew geometry
into a fictitious periodic slice.

## What Mesh X/Y sample

`simulation_sample.py` takes `P,Q` from the proven recast simulation cell and
constructs the legacy transverse metric `a=|P|`, `b=|Q|`, `gamma=angle(P,Q)`.
`simulation.py` builds its reciprocal dual and samples H/K on that metric.
`projected_potential_direct` evaluates fractional coordinates
`x=column/MeshX`, `y=row/MeshY`; complex arrays have shape `(MeshY, MeshX)`.
Thus the real-space sampling point is `r = x P + y Q`.

The propagation axis is `z_beam=C/|C|` along the selected [u,v,w]. The frame places
P along x_beam and y_beam=z_beam×x_beam. Q can be oblique in that transverse plane:
Mesh Y indexes the Q lattice coordinate, not necessarily pure Cartesian y_beam.
Step lengths are `|P|/MeshX` and `|Q|/MeshY`; only for gamma=90° are these Cartesian
dx/dy. These are projected-potential/wavefield samples, not CBED detector pixels,
and are not generally the original [001] cell mesh.

The non-[001] synthetic regression uses original cell diag(3,4,5) Å and zone [100].
The proven transverse vectors are P=(0,−4,0), Q=(0,0,−5), C=(3,0,0) Å. Mesh 8×6
therefore samples a 4×5 Å transverse cell with array shape (6,8), steps 0.5 and
5/6 Å, and a one-slice propagation period of 3 Å. The original [001] transverse
cell would instead be 3×4 Å.

## Validation scope

`tests/test_one_slice_simulation.py` covers one/two/four slices, periodic atom
ownership, finite complex projected potentials and gratings with/without
absorption, complex-wave comparisons after every occurrence against an explicit
pre-sliced model, non-[001] rectangular sampling, unchanged rejection of
incompatible geometry, depth refinement/rebuild, Save/reload and checkpoint
compatibility. An optional gfortran test extracts and compiles the actual ONEIT
sequence block into a temporary standalone harness; legacy source is read-only.

The Python pre-sliced comparison is not a claim of full historical Fortran-wave
agreement. Available saved Fortran scalar/listing tests are run separately; no
saved one-slice historical complex-wave fixture is present. Measured test and
build results are recorded in [macOS validation](../packaging/macos/VALIDATION.md).
