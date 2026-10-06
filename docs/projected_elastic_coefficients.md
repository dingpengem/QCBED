# Projected elastic coefficient E_g refinement

**Projected elastic coefficient $E_g^{(s)}$**: the slice-resolved transverse
Fourier coefficient of the elastic potential in the zone-axis-oriented
multislice simulation cell. Reflections are displayed using their original
crystallographic hkl identities.

QCBED currently refines the slice-resolved projected elastic Fourier coefficients
used directly by the multislice calculation. These are **not conventional 3D
unit-cell structure factors**. The physical identity is the original crystallographic
hkl **and** the fractional physical slice interval, for example `(200), [0, 0.5]`.
Internal array indices and detector beam order do not define a saved target.

## Physical convention

IAM E_g is calculated by the same complex128 atomic-sum function used in sample
preparation: source-scaled electron atomic factors divided by transverse cell area,
occupancy, exp(-B|g|²/4), and positive atomic phase using the retained printed
constant 6.283185307. The numerical units are the **source-scaled projected
coefficient**; V·Å is the documented dimensional interpretation, not an independently
established conversion to conventional electron structure factors in Å.

The grating uses **E_g + i A_g**. E_g and A_g are separately complex elastic and
absorptive coefficients. This editor changes E_g only. Absorption a/b/c remains a
separate site-group refinement family and is recalculated through its existing
kernel. Elastic IAM E_g has no direct voltage argument; voltage still affects
propagation, the interaction constant and absorptive scattering.

For the elastic coefficient, E(-g) = conj(E(g)). Each pair has one owner: the hkl
whose first nonzero component is positive. The GUI shows its linked negative
partner. A_g retains its own independently calculated conjugacy; the combined
E_g + i A_g generally does **not** obey the elastic conjugacy equation. No new
space-group equivalence is inferred from equal IAM values. Elastic (000) remains
read-only because a uniform elastic phase does not supply an intensity-fit target.
Cutoff/Nyquist-excluded pairs are unavailable.

## C → Refinement Targets → Projected Elastic Coefficients

Validate Simulation Setup, then use the dedicated table. It shows the owner and
linked reflection, slice interval, complex IAM value, selected mode, starting value,
physical hard bounds and status. On smaller windows, scroll the table horizontally
to reach Mode, starting values, bounds and Reset to IAM. Select a row to inspect
IAM and starting Re, Im, magnitude and phase (radians) in the details panel.
The panel provides individually labeled real/imaginary lower and upper bound
fields in Complex mode, synchronized with the table.

- **Real only:** optimize Re(E_g), retaining the committed imaginary component.
  This is the default when |Im(IAM E_g)| ≤ 1e-8 in the current coefficient convention.
- **Complex (Re + Im):** optimize two independent real-valued vector components.
  Enter a starting value such as `19.5 +0i`; enter separate bounds as `real, imaginary`.
- **Magnitude only:** optimize |E_g|, retaining the IAM phase. Its lower bound must
  be nonnegative. A zero refined magnitude is valid when IAM is nonzero and retains
  the IAM phase for subsequent proposals; Results reports the zero value’s phase
  as unavailable. Exactly zero IAM E_g has no phase, so this mode is unavailable there.

Magnitude + phase refinement is deferred; phase is not an independent optimizer
parameter in this milestone.

No universal physical bounds are supplied. Bounds must be finite, ordered and
contain the starting components. A positive `step_size` in the coefficient's units
is persisted in the configuration (default 0.01); it is a numerical search setting,
not a physical bound or a convergence guarantee.

**Apply Projected Elastic Coefficients** commits the table; validation/Save also
applies pending coefficient edits. **Reset to IAM** restores the current atomic
reference. Explicit starting overrides retain a model-reference fingerprint and
become stale when structure, lattice, occupancy, B, absorption, voltage or slicing
changes. Their values are not silently replaced; Reset or explicitly recommit the
starting value before Run. Targets with no starting override track IAM automatically.
Changed slice intervals require explicit removal/rebinding; **Clear saved coefficient
targets** removes coefficient records without changing the atomic model.

The initial GUI catalogue covers measured reflections and their partners. YAML
may also specify another original in-zone hkl within the mesh/cutoff. The public
parameter remains tied to its physical hkl and slice, regardless of beam ordering.

## Backend and optimizer integration

An absolute complex E_g is substituted at the selected reciprocal position and its
conjugate partner before the existing grating function. Absorption, masks, Fourier
signs, normalization, and multislice propagation are unchanged. IAM-valued settings
use the exact baseline path. Direct replacement handles zero/near-zero IAM values
without dividing by IAM E_g.

The old real-multiplier compatibility path remains unchanged. Where multiplication
is equivalent, `E_refined = m * E_IAM`; a real scalar cannot represent general
independent complex components. Legacy records do not carry original-cell hkl and
physical slice ownership, so old YAML stays on its compatibility path rather than
being silently reinterpreted as modern E_g targets.

Components append deterministically after the existing physical vector, ordered by
original hkl, slice interval and real/imaginary component. Arbitrary subsets can mix
with depth, B, absorption, voltage (without absorption) and scan targets. Physical
component values, bounds and histories remain signed and include exact zero.

For vectors containing E_g targets, ANNEAL, AMOEBA and PARAMIN use an explicit affine
adapter `x = stage_initial + (u - 1) * bound_width`. This avoids multiplicative search
singularities at zero without changing physical initial values or the legacy
engines. AMOEBA/PARAMIN use the stored per-component steps as physical displacements;
PARAMIN requires steps smaller than bound widths. SIMPARA uses physical components
and its existing absolute uncertainty widths. These new-vector coordinate semantics
do not change older runs without E_g targets. Legacy D1/D2 simplex formulas should
not be applied to the new affine coordinates.

## Compact YAML and Resume

Records live under `qcbed.refinement.projected_elastic`. For example, in a one-slice
model (bounds here are illustrative, **not universal recommendations**):

```yaml
projected_elastic:
  - hkl: [2, 0, 0]
    interval: [0.0, 1.0]
    mode: complex
    refine: true
    lower: [0.0, -5.0]
    upper: [40.0, 5.0]
    step_size: 0.01
```

Omitting `initial` uses current IAM. Explicit overrides serialize `initial: [real,
imaginary]` and an editor-generated `reference` fingerprint. The YAML contains no
second project file, coefficient grids or opaque reciprocal array identity. Verbose
canonical YAML stores the same records under `crystallography.projected_elastic`.
Original hkl/slice ownership, mode, physical starting values and bounds participate
in checkpoint identity. Incompatible targets, source models and reordered experimental
data are rejected for Resume rather than reassociated by column number.

## Results

The **Projected Elastic Coefficients** Results tab derives IAM, initial, displayed,
latest, best and retained values from the same recorded physical-vector journal.
It shows original hkl, slice interval, mode, Re(E_g), Im(E_g), magnitude, phase in
radians, bounds and the fixed-imaginary constraint for Real-only. Unavailable states
remain unavailable after initialization/cancellation; viewing never reevaluates an
objective. Selecting a row selects its component history and the matching original-hkl
beam diagnostics when that beam was measured. The existing Pattern Matching tab
contains the experimental/simulated comparison.

## Future scientific extension

**Conventional 3D structure-factor refinement** is a separate future feature. It
requires an explicit mapping from a unit-cell F_hkl to the slice-resolved projected
coefficients, including normalization, slice distribution, zero coefficients and
changed-slicing behavior. This milestone does not claim that mapping or automatic
space-group constraint inference.

## Optional user equality groups

Coefficient targets are **Independent** by default. Select two or more owner rows
and choose **Share One Refinement Parameter**, then name the **Shared refinement parameter group**. **Make Independent**
removes the whole selected group and restores independent coordinates. This is a
user constraint, not an automatically inferred space-group symmetry rule.

| Mode | Shared optimizer coordinates | What remains member-specific |
|---|---|---|
| Real only | One real component | Each fixed imaginary component |
| Complex | One real and one imaginary component | IAM reference and beam diagnostics |
| Magnitude only | One magnitude | Each fixed IAM phase |

Friedel partners are still mandatory conjugate partners, separate from optional
user groups. A group has at least two owners with the same physical slice interval,
mode and Refine activation. Different current starts require an explicit **Use first
selected** choice; no average is invented. Alternatively, edit starts to a common
value before sharing. Bounds are the intersection of member domains; existing
member bounds remain stored. Empty intersections or a start outside the intersection
are rejected with an explanation. No universal group bounds are introduced.

Group membership is stored directly on each physical record using `equality_group`:

```yaml
projected_elastic:
  - hkl: [2, 0, 0]
    interval: [0.0, 1.0]
    mode: real
    refine: true
    equality_group: group_1
    lower: [19.0, null]
    upper: [20.0, null]
  - hkl: [0, 2, 0]
    interval: [0.0, 1.0]
    mode: real
    refine: true
    equality_group: group_1
    lower: [19.0, null]
    upper: [20.0, null]
```

These bounds are illustrative. Omitting the group field preserves independent
behavior and existing YAML compatibility. Target generation resolves original hkl
and slice, never row positions. Changed modes/slices, missing members or incompatible
starts make a group incomplete rather than reassigning it. Explicit starts retain
the existing stale-IAM rules after upstream atomic-model changes. Groups survive
beam reordering; incompatible checkpoint identities are rejected.

The optimizer layout contains one shared coordinate set per group. Its search
step uses the first member in canonical original-hkl order; changing member order
in YAML does not change that choice. Results labels
include the group and members, while individual coefficient rows retain each IAM
reference, fixed component/phase, initial/latest/best/retained states and shared bounds.
Member rows select the shared history and their own separate beam diagnostics.

See the [25-evaluation audit and stopping semantics](25_iteration_audit.md) for
iterations versus evaluations, exact optional iteration caps, and Resume behavior.


### Parameter sharing in the editor and Results

The **Parameter sharing** column shows **Independent** or **Shared: Group 1**.
The live summary reports active physical targets, shared parameter groups and
optimizer parameters. Draft estimates are labeled separately from the validated compiled layout. Mandatory Friedel partners
are excluded from these counts. Sharing is enforced through common optimizer
coordinates throughout refinement, not by averaging final results.

The confirmation previews the dimensionality change: depth plus two Real-only
coefficients changes from three coordinates to two. Two Complex coefficients
share two coordinates (one real and one imaginary), not one scalar.
**Shared real component** leaves each imaginary component fixed;
**Shared complex coefficient** shares both components; **Shared magnitude**
retains each member's IAM phase. **Make Independent** restores separate coordinates.

Results labels these as **Shared refinement parameter group**, while keeping
per-reflection IAM values and beam diagnostics separate. Existing YAML
`equality_group` keys and checkpoint identities remain compatible; stored historical
labels are translated for display without modifying checkpoint files.


### Resolved refinement parameter layout

After **Validate Refinement Targets**, the **Resolved Refinement Parameter Layout**
shows each shared group, its original hkl/slice members, mode, optimizer coordinates,
starting values and intersected hard bounds. Independent targets and mandatory
Friedel relationships are listed separately. The summary and D dashboard distinguish
**Active physical targets** from **Optimizer parameters** using the compiled layout.
Depth plus two shared Real-only coefficients therefore has three physical targets
and two optimizer coordinates. Multiple groups remain separate; Complex groups
contribute two coordinates each. Draft edits require revalidation.

Coefficient-only selections are supported, including category-filtered configurations:
the model preview may contain zero ordinary coordinates before elastic coordinates
are appended. Sharing remains enforced in the optimizer layout and writeback, not
by copying independently optimized results afterward. YAML and Resume use the
existing physical identities and group representation.


### Sharing through the validation buttons

**Validate Refinement Targets** commits the visible coefficient checkboxes, modes,
starting values and bounds before validating the complete layout. No separate
**Apply Projected Elastic Coefficients** step is required. **Share One Refinement
Parameter** and **Make Independent** immediately commit the sharing change, clear
refinement/C validation and disable Run/Resume readiness. Their feedback reports
the actual compiled vector length before and after the change and asks for
revalidation. For depth plus two Real-only coefficients this is **3 → 2 → 3**.

After validation, the optimizer count is the compiled starting-vector length. Before validation, compatible draft selections have a labeled estimate. The earlier hypothetical “if independent” count is no longer presented as the active optimizer dimension.
