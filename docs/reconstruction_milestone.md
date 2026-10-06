# QCBED Python reconstruction milestone — release closure

> Historical milestone/audit record. For current controls and the A/B/C/D workflow,
> use the [User Manual](user_manual.md); dated counts and replay inputs below retain
> their original provenance.

Accepted 2026-10-04. This milestone preserves the legacy numerical conventions
while providing a canonical Python/YAML and Apple Silicon desktop workflow.

- Live Patch A begins with raw CBED, without a historical prepared manifest.
- An optional aperture supplies threshold geometry and the existing PSF/OTF
  correction and noise-characterization path. Without calibration, explicit
  manual disk geometry and labelled Poisson uncertainty remain available.
- Interactive centrosymmetric feature selection, arc-based lattice detection,
  live A10 corners and ordered beam extraction feed Pattern/Sigma/Weight views.
  Manual bad-pixel exclusion changes weights only.
- Sequential A/B/C validation connects an immutable validated preparation to
  conventional-cell crystallography, explicit hkl mapping, atomic model/site
  groups and the supported modern requested-depth/shared-B refinement workflow.
- Conventional and gradient backends, Run/Resume, cooperative Cancel and Results
  remain supported. Prepared descriptors are generated replay/cache artifacts,
  not required primary inputs.
- The frozen native arm64 `QCBED.app` retains the Python/CLI development workflow.

## Accepted evidence

128 focused tests, 102 GUI tests and 2880 full-suite tests passed. Diff whitespace
check was clean. Real raw FCC Al preparation and frozen automated Qt acceptance
reached canonical run initialization with **zero optimizer evaluations**. The
intentional pre-optimizer stop produces a failed test checkpoint, not a completed
refinement. No refinement campaign was repeated for acceptance or release closure.

The live seed is (81,9), (-9,82), with explicit 38 px radii for this case and
central/symmetry-center origin (127,126). It is not forced to match reviewed historical geometry.
See [the workflow](gui_abc.md) for exact corners and handoff details.

## Known limitations

- No matched real aperture-image end-to-end acceptance yet; component and
  orchestration tests are distinct from that outstanding real-data check.
- Poisson fallback is approximate when aperture calibration is unavailable.
- Real A7 global lattice refinement can remain unavailable when independent arc
  support is sparse; missing confirmation remains nonfatal subject to ROI coverage.
- Final acceptance used automated Qt interactions, not a human Finder/file-picker
  session.
- The app is ad-hoc signed, not Developer ID signed or notarized.

## Build provenance and artifacts

[Recorded provenance](../packaging/macos/RELEASE_PROVENANCE.json) captures the
actual build base commit, dirty-tree flag, packaged-source fingerprint, UTC build
time, interpreter/dependency versions, architecture, artifact checksums, signature
status and accepted test counts. The accepted build predates release closure:
its base commit must not be misrepresented as a clean release source revision.
The packaged-source fingerprint was recomputed and matches the accepted build.
Closure changes documentation only; no rebuild or numerical tests were needed.

Final local artifacts (generated and intentionally not committed):

- `dist/QCBED.app`
- `dist/QCBED-arm64.zip`

The editable `examples/fcc_al.yml` replay example references locally generated
`tmp/` preparation products; those products are not bundled as portable sample
data. New users should follow raw-CBED preparation rather than assume that replay
cache exists in a fresh checkout.

Local output artifacts and the editable example's output-path customization are
separate from the accepted milestone source. They are preserved, not deleted or
silently included in a release commit. No history rewriting or remote push is
part of closure.
