# Complex structure-factor refinement: pre-implementation audit

Status: **historical pre-implementation audit**. The user subsequently selected
per-slice projected elastic coefficient E_g, identified by original hkl and physical
slice interval. See the [implemented workflow](projected_elastic_coefficients.md).
The measurements below document the original mathematical contract, not release
acceptance of the later feature.

## Existing mathematical contract

**Projected elastic coefficient $E_g^{(s)}$**: the slice-resolved transverse
Fourier coefficient of the elastic potential in the zone-axis-oriented
multislice simulation cell. Reflections are displayed using their original
crystallographic hkl identities.

For stored slice s and internal reciprocal index G = (H, K), the elastic kernel
computes the complex, source-scaled projected coefficient

```text
E_s(G) = sum_j [f_source(Z_j, |g|/2) / area_s]
                 * occupancy_j * exp(-B_j * |g|^2/4)
                 * exp(+i * 6.283185307 * (H*x_j + K*y_j)).
```

Modern fixed occupancy comes from the physical site group; vacancy is zero and
dilation is unity. B is in Å². `f_source` is the existing interpolated electron
atomic factor with its potential-factor scaling, not an unscaled electron
scattering length. There is no elastic voltage argument, slice-thickness factor,
or 3D cell-volume normalization in this sum. Voltage enters propagation and the
interaction constant; absorptive atomic factors additionally depend on voltage.

The current units contract is **source-scaled projected coefficient**. V·Å is the
existing dimensional interpretation (atomic potential factor interpreted as V·Å³,
divided by area in Å²), explicitly documented as an inference. Relabeling these
numbers as conventional electron structure factors in Å would be incorrect.
See [structure factors](structure_factors.md) and [atomic scattering](scattering.md).

The separate absorptive kernel computes

```text
A_s(G) = -sum_j area_s * f_abs(Z_j, |g|/2, voltage, B_j, a_j, b_j, c_j)
                   * occupancy_j * exp(-B_j * |g|^2/4)
                   * exp(+i * 6.283185307 * (H*x_j + K*y_j)).
C_s(G) = E_s(G) + i*A_s(G).
```

The absorptive area multiplication and sign above are the existing implemented
contract; they must not be replaced by the elastic area division. Cutoff and
Nyquist masks are applied to C. An unnormalized negative-exponent Fourier sum
constructs the projected potential V_s, and the grating is `exp(-i*sigma*V_s)`.
Existing propagation interleaves each slice grating and its propagation kernel.

The compatibility refinement parameters are real, dimensionless, region/reflection
multipliers:

```text
E_s,effective(G) = m_elastic(G) * E_s,IAM(G)
A_s,effective(G) = m_absorptive(G) * A_s,IAM(G)
```

A region multiplier applies to every stored slice in that region. It is neither
an absolute F_g nor an independent real/imaginary component of C. Reflection
matching uses exact internal H/K first, then configured symmetry, with the first
matching record winning. Duplicate h/k records keep distinct configuration IDs
but are not independent controls of the same grid coefficient.

## Trace through the current implementation

| Layer | Implementation and finding |
|---|---|
| Atomic factors | [scattering.py](../src/qcbed/scattering.py): electron branch, existing interpolation and source scaling |
| Elastic IAM | [structure_factors.py](../src/qcbed/structure_factors.py): complex128 per-slice sum, positive atomic phase, occupancy and applied B |
| Absorptive IAM | [absorption.py](../src/qcbed/absorption.py): separate complex sum, voltage/B/absorption coefficients, existing sign and area convention |
| Grating and propagation | [phase_grating.py](../src/qcbed/phase_grating.py), [simulation.py](../src/qcbed/simulation.py): E+iA, masks, Fourier reconstruction, grating, then sequential slice propagation |
| Multiplier application | [refinement_trial.py](../src/qcbed/refinement_trial.py): `_reflection_scales`, separate real multipliers before grating construction |
| Original reflection identity | [index_transform.py](../src/qcbed/crystallography/index_transform.py): exact integer `HKL = hkl @ T.T`, checked inverse; [beam_indexing.py](../src/qcbed/crystallography/beam_indexing.py) retains original hkl and stable beam ID |
| Configuration/vector | [config_compiler.py](../src/qcbed/config_compiler.py), [refinement.py](../src/qcbed/refinement.py): `(region_id, reflection_id, elastic/absorptive)` becomes ordered region/record/component coordinates |
| Modern bridge | [crystallography_config.py](../src/qcbed/crystallography_config.py), [site_groups.py](../src/qcbed/crystallography/site_groups.py): no physical reflection alias; strict prepared sample has unity multipliers |
| YAML | [configuration.py](../src/qcbed/configuration.py): pre-sliced records store IDs, h/k and two multipliers; [modern_yaml.py](../src/qcbed/modern_yaml.py) has no physical complex-F_g schema |
| Modern optimizer/Resume | [crystallography_refinement.py](../src/qcbed/crystallography_refinement.py), [grouped_run.py](../src/qcbed/grouped_run.py): physical aliases, explicit bounds, original ownership and configuration/data identity; new F_g metadata must participate in this identity |
| Results | [canonical_run.py](../src/qcbed/canonical_run.py) records compatibility multiplier labels; [grouped_results.py](../src/qcbed/gui/grouped_results.py) reconstructs modern physical-vector histories, currently without complex F_g metadata |

## Conjugacy, symmetry and zero coefficients

Real atomic amplitudes, occupancies and B imply separately
`E_s(-G) = conj(E_s(G))` and `A_s(-G) = conj(A_s(G))`. However,
`C_s(-G) = conj(E_s(G)) + i*conj(A_s(G))`, generally **not** `conj(C_s(G))`.
Elastic real/imaginary refinement should explicitly link Friedel partners; it
must not conjugate the entire combined elastic-plus-absorptive coefficient.
Legacy independently selected multipliers need not preserve these relationships,
so automatic conversion cannot silently repair arbitrary compatibility inputs.

The strict modern bridge does not infer a space-group reflection-refinement
constraint system. The old NONE/INVE/RECT/SQUA/HEXA matching options are not a
validated original-cell space-group ownership model. Equal IAM values alone do
not establish permanent refinement links. Further crystallographic constraints
must be explicit if imposed; automatic full symmetry grouping is not implemented.

Elastic E(000) is real. Changing it alone introduces a uniform slice phase and
does not identify an intensity-fit parameter. Imaginary E(000) violates the
elastic conjugacy constraint. The absorptive zero coefficient has a different
role. The proposed elastic GUI should therefore keep (000) read-only.

A real multiplier only moves E along a line through its complex IAM value. It
cannot in general implement independent complex refinement, or change only Re(E)
while retaining a nonzero Im(E). A complex ratio could do so for nonzero IAM E,
but is not the existing real-multiplier contract. At zero/near-zero IAM E, division
is invalid/unstable. Direct replacement of the chosen elastic coefficient before
the existing grating calculation avoids division, once its slice ownership is
specified. Magnitude-only mode has no IAM phase when IAM E is exactly zero; that
case needs an explicit unavailable state or a separately chosen phase, not an
implicit phase fabricated from a multiplier.

Bounds must be physical component bounds. For a real IAM value and a real
multiplier interval, transformed endpoints must be sorted when the IAM value is
negative. A complex IAM value maps a real multiplier interval to a line segment,
not an independent rectangular real/imaginary bound pair. No universal bounds
are established by the existing implementation.

## FCC Al teaching-model calculation

These are newly calculated **IAM projected elastic coefficients**, not refinement
results. Inputs match the teaching model: [Al_FCC.cif](../reference_data/Al_FCC.cif),
a = 4.03893 Å, zone [001], g1 = (200), g2 = (020), B = 0.86 Å², occupancy 1,
200 kV, 16×16 mesh, cutoff 2.5 Å⁻¹, and one stored slice [0,1]. Absorption is
disabled, as in [fcc_al.yml](../examples/fcc_al.yml). Transverse area is
16.312955544899996 Å². The separate historical preflight extxyz uses a = 4.049 Å;
its coefficients differ and must not be substituted for the teaching CIF values.

| Original hkl | Internal HK | Re(E_IAM) | Im(E_IAM) |
|---|---|---:|---:|
| (0,0,0) | (0,0) | 69.13611558002403 | 0 |
| (-2,0,0) | (-2,0) | 19.591754648483377 | +1.7592060366859781e-9 |
| (2,0,0) | (2,0) | 19.591754648483377 | -1.7592060366859781e-9 |
| (0,-2,0) | (0,-2) | 19.591754648483377 | +1.7592060366859781e-9 |
| (0,2,0) | (0,2) | 19.591754648483377 | -1.7592060366859781e-9 |
| (-2,-2,0) | (-2,-2) | 11.731316291670417 | +2.1067844926464262e-9 |
| (2,2,0) | (2,2) | 11.731316291670417 | -2.1067844926464262e-9 |
| (-2,2,0) | (-2,2) | 11.73131629158207 | 0 |
| (2,-2,0) | (2,-2) | 11.73131629158207 | 0 |

All nine mappings round-trip through the exact inverse. Maximum imaginary
magnitude is 2.107e-9; it is below 1e-8 in the current coefficient convention and
consistent with the deliberately retained printed phase constant. Real-only is
a suitable example default, retaining rather than discarding these small imaginary
values. This is not a global real-only assumption for other structures/origins.

Reconstructing the grating with these same full IAM coefficient arrays and the
production grating function gives bitwise equality with the prepared baseline:
maximum absolute grating difference **0.0**. No optimizer evaluations were used.
This establishes the audited source of the numbers; it does not validate a new
refinement implementation, GUI display or YAML conversion.

## Definition required before implementation

The existing simulation stores E_s(G) for **each stored slice**. A single number
per hkl could instead mean a sum over a repeat cell or a conventional 3D unit-cell
structure factor. These choices have different normalization and ownership.
For multiple slices, specifying only the total does not determine the individual
slice gratings; propagation depends on their distribution along the beam.

The least-assumptive extension of the existing kernel is a **per-slice projected
elastic coefficient**, bound by original hkl plus a stable physical slice interval,
with linked Friedel partners. If the desired public parameter is instead one
whole-cell F_g per hkl, a distribution rule must first be agreed, including its
zero-IAM and changed-slicing behavior. Neither choice has been implemented or
silently selected in this audit.

Once that definition is settled, implementation still needs the physical target
schema and stale-IAM tracking; three refinement modes; original-hkl and conjugacy
ownership; signed/zero-safe optimizer coordinates; GUI, compact YAML, Resume and
Results integration; and end-to-end regression and frozen-app acceptance.
AMOEBA and PARAMIN currently reject zero physical starts. Complex imaginary
components commonly start at zero, so adding a catalogue entry alone is not
sufficient: an explicit modern coordinate adapter is needed without changing the
legacy packing floor or the actual physical start.

## Validation performed for this audit

The following existing focused tests passed: **152 passed in 0.21 s**.

```sh
python -m pytest -q tests/test_structure_factors.py tests/test_absorption.py \
  tests/test_refinement_execution.py::test_reference_absorption_and_nonunit_elastic_absorptive_multipliers \
  tests/test_refinement_execution.py::test_multiple_stored_slices_and_reflection_corrections \
  tests/test_configuration.py::test_duplicate_beams_and_distinct_duplicate_reflection_targets
```

They cover independent hand-calculated complex sums, nonzero imaginary parts,
occupancy/B effects, Friedel relationships, saved numerical listings, existing
multiplier behavior, and duplicate reflection identity. They are not tests of an
implemented modern complex-F_g feature. No scientific kernels were modified.
The full suite and app build/signing are deferred until the feature exists and its
focused acceptance checks are stable. No expensive refinement campaign was run.
