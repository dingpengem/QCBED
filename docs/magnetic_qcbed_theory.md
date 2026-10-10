# Magnetic QCBED: theory, conventions, and inverse refinement

**Status:** Research-stage derivation, 2026-10-10. This document specifies a future
QCBED research module; it does **not** claim that induced magnetization is already
quantitatively recoverable from experiment. Conventions and neglected terms are
part of the model specification, not implementation details.

## 1. Scope and observables

The objective is to simulate and (eventually) infer magnetic order parameters
from quantitative convergent-beam electron diffraction (CBED) by combining
strong dynamical **electrostatic** scattering with much weaker elastic
**magnetic** scattering. NiO provides a Type-II antiferromagnetic test case;
there is no NiO-specific assumption in the generic propagator. Magnetic elastic
CBED differs from EMCD (an inelastic dichroism measurement).

The measured data are the detector intensities indexed by scan/tilt, beam, and
pixel, with experimental standard deviation `sigma` and **binary** final weight
`w=0/1`, consistent with the existing QCBED data-preparation contract.

The proposed primary quantities are (1) local staggered sublattice moment;
(2) magnetic order-axis direction; (3) canting angle; (4) change between
applied-field conditions; and, later, (5) constrained magnetic Fourier
components. Field-induced *net* magnetization must be distinguished from a
reorientation of a pre-existing antiferromagnetic order parameter.

## 2. Electron phase convention, SI units, and paraxial limit

Let an electron have charge `-e`, `e>0`, incident kinetic energy `e U`,
free momentum `p0`, relativistic velocity `v`, and `k0=p0/hbar=2 pi/lambda`
(**radians/m**, not cycles/m). Its full wave has carrier `exp(+i k0 z)` and
envelope `Psi(r_perp,z)`. A scalar relativistic small-angle treatment using
minimal coupling yields, to the indicated paraxial order [1],

\[
 i\partial_z\Psi=\left[\frac{(i\nabla_\perp+e\mathbf A_\perp/\hbar)^2}{2k_0}
 +\frac{e A_z}{\hbar}-\frac{e\Phi}{\hbar v}\right]\Psi.
 \tag{1}
\]

Here `Phi` is electric scalar potential in V, `A` is magnetic vector potential
in T m (Wb/m), and `B=curl A` is magnetic induction in T. Equation (1)
neglects electron spin coupling and higher-order relativistic/paraxial terms;
see [1,2] for conditions and more complete Pauli-multislice treatments.
This choice of carrier, electron charge and gauge fixes the **signs**.

At `A_perp=0`, the free-space operator is `-nabla_perp^2/(2k0)`; a transverse
plane wave `exp(i k_perp.r_perp)` acquires

\[
 P(\Delta z,\mathbf k_\perp)=\exp[-i\Delta z|\mathbf k_\perp|^2/(2k_0)].
 \tag{2}
\]

An idealized thin transmission slice with electrostatic potential and
longitudinal vector potential has

\[
 T_s(\mathbf r_\perp)=\exp\left[
 +i\frac{e}{\hbar v}\int_s\Phi\,dz
 -i\frac{e}{\hbar}\int_s A_z\,dz
 \right].\tag{3}
\]

Thus `sigma_e=e/(hbar*v)` [rad/(V m)] and `sigma_A=e/hbar` [rad/Wb].
The projected longitudinal magnetic field integral `int A_z dz` has **Wb**
units, not T, V, or V angstrom.

**Compatibility caution:** Legacy QCBED builds a grating from its own
source-scaled projected `E_g + i A_g`, with a negative electrostatic grating
phase by the repository's conventions. Do **not** silently identify the
legacy projected potential with `Phi` in (3) without establishing its sign,
normalization and units. The research adapter multiplies an already-built
existing QCBED phase grating by the magnetic factor only; it preserves its
nonmagnetic propagator, FFT ordering, absorption and other validated physics.

## 3. Magnetostatics: from magnetization to vector potential

Define magnetization `M(r)` in A/m and write, for a **three-dimensionally
periodic** magnetic cell,

\[
 \mathbf M(\mathbf r)=\sum_\mathbf q\mathbf M_\mathbf q
 e^{i\mathbf q\cdot\mathbf r},\qquad \mathbf q=2\pi\mathbf g.
 \tag{4}
\]

Here `g` is in cycles/m (`1/angstrom` when converted from crystallographic
reciprocal vectors), while `q` is in radians/m. Magnetostatics with no free
current gives

\[
 \nabla\times\mathbf B=\mu_0\nabla\times\mathbf M,
 \quad\mathbf B=\nabla\times\mathbf A,
 \quad\nabla\cdot\mathbf A=0.\tag{5}
\]

Using `nabla^2 A=-mu0 curl M`, for each `q != 0`:

\[
 \boxed{\mathbf A_\mathbf q=
 \frac{i\mu_0}{q^2}\big(\mathbf q\times\mathbf M_\mathbf q\big)}
 \qquad(q\ne0),\tag{6}
\]

\[
 \boxed{\mathbf B_\mathbf q=\mu_0\left[
 \mathbf M_\mathbf q-\frac{\mathbf q(\mathbf q\cdot\mathbf M_\mathbf q)}{q^2}
 \right]}.\tag{7}
\]

Equations (6–7) enforce `q dot A_q=0`, `q dot B_q=0`, Hermitian symmetry
`A_-q=A_q*`, and magnetic reversal `A[-M]=-A[M]`. They express the
transverse projection of magnetization. Longitudinal magnetization Fourier
components parallel to `q` do not contribute to `B_q` in this periodic
magnetostatic treatment; consequently not every magnetic Fourier component
is measurable independently.

### Why `q=0` is not an optional detail

The formula above **does not define** the `q=0` vector potential. Its finite
limit is generally not obtained by dividing by a tiny `q`. For uniform
magnetization, a specimen's demagnetizing field depends on geometry,
boundaries, external applied field, and magnetic domains. A uniform magnetic
induction may be represented, for example, by a nonperiodic vector
potential `A = (B cross r)/2`, which cannot be represented by a periodic
`A_z` grating alone. Even if an AFM has no net `M` at zero field, spin
canting introduces a nonzero spatial average. **Therefore the initial
periodic prototype cannot refine total induced net magnetization** without
a separate magnetostatic boundary/beam-optics solution. The prototype
returns the mean magnetization explicitly but excludes it from its computed
`A` and `B` grids, rather than silently setting a physical uniform field to 0.

For a finite foil, even nonzero periodic harmonics and their boundary fields
need a 3D field solution or benchmarked approximation. Periodic-z FFTs are
suitable for an *infinite periodic bulk reference*, not the vacuum/foil
boundary problem.

## 4. Magnetic unit cell, structure, and Fourier form factors

A magnetic propagation vector can enlarge the magnetic unit cell beyond the
chemical unit cell. Source magnetic sites are represented by position `r_j`,
local moment vector `mu_j` in A m² (`mu_B`), site occupancy and a normalized
magnetic density / form factor `f_j(q)`:

\[
 \mathbf M_\mathbf q=\frac{1}{\Omega}\sum_j
 \boldsymbol\mu_j f_j(\mathbf q)e^{-i\mathbf q\cdot\mathbf r_j}.
 \tag{8}
\]

`Omega` is the magnetic cell volume in m³. A normalization with
`f_j(0)=1` assigns the specified integrated magnetic moment to each site.
The present demonstration uses an **isotropic Gaussian** form factor
`f(q)=exp(-sigma_r² q²/2)` strictly as a geometry/units fixture; it must be
replaced by element-resolved, validated magnetic form factors, DFT spin
density or experimentally justified density before quantitative NiO claims.
NiO orbital contributions, noncollinearity, thermal magnetization effects and
local oxidation/bonding can matter. Published parametrized magnetic potentials
and fields [3] offer a physically motivated future source.

For ideal Type-II NiO, alternation between consecutive ferromagnetic
`(111)` Ni planes can be represented on the Ni fcc sublattice as
`eta_j=(-1)^(x_j+y_j+z_j)` for positions in units of the cubic lattice
constant; these sums are integers for the fcc Ni sites. A 2x2x2
conventional cubic supercell is a simple commensurate **illustrative** fixture,
not a unique minimal magnetic cell. Different T-domains correspond to
symmetry-related ordering wavevectors and different S-domains to order-axis
choices; each must be modeled in the beam-aligned coordinates.

## 5. Two-sublattice canting and the inverse question

For a simple symmetric two-sublattice order, let `n` and `h` be perpendicular
unit directions and `theta` a cant angle:

\[
 \mathbf\mu_A=\mu_0(\cos\theta\,\mathbf n+
 \sin\theta\,\mathbf h),\quad
 \mathbf\mu_B=\mu_0(-\cos\theta\,\mathbf n+
 \sin\theta\,\mathbf h).\tag{9}
\]

Then, **per Ni site** (not for the total magnetic cell),

\[
 m_{\rm AFM}=\mu_0\cos\theta,\qquad
 m_{\rm ind}=\mu_0\sin\theta.\tag{10}
\]

The vector `h` is the canting direction; its relationship with applied `B`
depends on anisotropy, magnetic susceptibility and domain state. The very
simple Eq. (9) must not be mistaken for a field-calibrated spin-flop model.

The staggered periodic component in (9) scales with `mu0 cos theta`, while
its homogeneous component scales with `mu0 sin theta`. If the simulator
omits `q=0`, the magnetic data generally cannot separately determine `mu0`
and `theta` from the staggered amplitude alone. At `theta=0`,
`d(mu0 cos theta)/dtheta=0`: the leading sensitivity to small canting is
quadratic in theta (and can be even less favorable after imaging effects).
The forward model or data must provide additional independently sensitive
measurements to recover `m_ind`, rather than merely report an optimizer
best-fit. This degeneracy is a central falsification test for the project.

## 6. Magnetic multislice model, first practical release

A symmetric split operator for the longitudinal scalar approximation is

\[
 \Psi_{s+1}\approx P(\Delta z_s/2)
 T_s(\Phi,A_z)P(\Delta z_s/2)\Psi_s .\tag{11}
\]

Different splitting order from the existing QCBED engine is **not** by
itself evidence of exact numerical compatibility. The first direct
integration into the mature engine uses the existing stored gratings:

\[
 T_s^{\rm combined} = T_s^{\rm QCBED,electric}
 \exp[-i e\mathcal A_{z,s}/\hbar],\quad
 \mathcal A_{z,s}=\int_s A_z\,dz,\tag{12}
\]

and then runs the **unmodified** QCBED sequence simulator, slice sequence,
propagators and beam extraction. The research prototype approximates the
projection by midpoint `A_z(z_mid) Δz`. Convergence with finer z slicing
is mandatory before physics acceptance. QCBED internal HK and original
crystallographic `hkl` remain distinct; magnetic satellites require extending
the reciprocal basis to the magnetic supercell and cannot always be mapped
to integral hkl in the chemical-cell basis.

**Limitations:** Retaining `A_z` alone is fixed to a specified Coulomb gauge
and is **not** generally gauge invariant. It is justified only when the
transverse vector-potential terms are demonstrably small for the beam/field
geometry being modeled. The complete transverse kinetic operator in Eq. (1),
spin/Pauli coupling (including spinor transport), foil boundaries,
macroscopic applied field, Lorentz deflection, magnetostriction, atomic
magnetic form factors, phonon/magnon diffuse scattering, detector PSF and
noise, and domain mixing are **not** contained in this scalar prototype.
Structural absorption may still be represented by the existing grating.

## 7. QCBED magnetic inverse problem

For settings `b` (field/temperature/probe orientation) and active detector
samples `i`, introduce structural nuisance parameters `s`, instrument/
scale parameters `n_b`, and magnetic parameters `m_b`:

\[
 I^{\rm pred}_{bi}=S_b\,I_{bi}(\mathbf s,\mathbf m_b,\mathbf n_b)
 +b_{bi}.\tag{13}
\]

The current QCBED final weight is a **binary inclusion mask**; experimental
variance (not a fractional weight) supplies the denominator:

\[
 \chi^2=\sum_{b,i} w_{bi}
 \frac{(I^{\rm obs}_{bi}-I^{\rm pred}_{bi})^2}{\sigma_{bi}^2}
 +R(\mathbf s,\mathbf m,\mathbf n).\tag{14}
\]

`R` may implement explicit priors or physical constraints. It must be
reported separately from data misfit. The implemented conventional/gradient
objective modes have their own existing normalization/stencil semantics;
the magnetic research module must not bypass them when joining canonical
refinement. For raw low-count magnetic spots, a Poisson likelihood may be
preferable to the Gaussian approximation: choose the statistical model
before interpreting errors. Correlated camera noise, normalization drift,
and diffuse background also matter.

Use a *joint* field-series model with common chemical structure, geometry,
thickness (where valid), occupancies and Debye–Waller factors, while allowing
field-dependent magnetic moments/domains and separately calibrated
field-dependent beam shifts/illumination. Uniform field from the objective
lens must not be confused with a specimen-intrinsic magnetic change.

### Local sensitivity and identifiability

Let `p` be the fitted **physical** parameters and `J_ij=dI_i/dp_j`. After
selecting active pixels and using experimental sigma, form the whitened
Jacobian `Jw_ij=J_ij/sigma_i` and Fisher/Gauss–Newton matrix

\[
 F=J_w^T J_w.\tag{15}
\]

The scaled singular values and null vectors of `Jw`, not simply whether the
optimizer stops, establish local identifiability. To test whether magnetic
parameters `m` are identifiable **after nuisance adjustments**, profile out
nuisance parameters `n`:

\[
 F_{m|n}=F_{mm}-F_{mn}F_{nn}^{+}F_{nm}.\tag{16}
\]

`+` denotes a pseudoinverse; a deficient profile or a small eigenvalue
indicates degeneracy. Choose meaningful scales for parameters (e.g. `mu_B`
for moment, radians for angle, angstroms for thickness) because raw singular
values depend on units. Gaussian covariance `F^{-1}` is meaningful only for
full-rank regular local problems with a correctly calibrated noise model.
**Pseudoinverse diagonal entries of a singular Fisher matrix are not valid
finite confidence intervals along null directions.** Ultimately perform
profile likelihood or synthetic-recovery tests with nuisance variation.

Needed diagnostics: separate sensitivity per disk, Fisher condition number,
magnetic-to-nuisance correlations, field-reversal differences, convergence
with slice thickness/grid/beam geometry, and out-of-sample synthetic recovery.

## 8. Lorentz-CBED experimental design and controls

A low-field Lorentz setting permits nearly field-free measurements of AFM
order, not automatic measurement of magnetism. Its reduced numerical aperture
and optical reconfiguration may affect CBED probe size, intensity and camera
geometry. Calibrate actual field *at the specimen*, probe semiangle,
illumination, position stability and sample thickness in each mode.

1. Start with field-free magnetic CBED feasibility at favorable magnetic
   superlattice reflections and known zone axes; compare with structural CBED.
2. Measure the same region with carefully calibrated external field or
   temperature and matched electron optics, alongside nonmagnetic reference
   regions. Record domain geometry and whether magnetic period changes.
3. Use field reversal, sample reversal/domain variations and measurements
   above the ordering temperature (where feasible) as discriminating controls.
4. Separate systematic disk displacement/distortion from internal magnetic
   intensity modulation; model objective-lens and holder fields along the path.
5. For quantitative parameters, measure enough distinct zone axes, foil
   thicknesses or field conditions to break degeneracies, and report limits
   of detection even when no induced moment can be uniquely recovered.

The approach builds on direct NiO magnetic diffraction [4] and the magnetic
point-group CBED calculations of Zaiets *et al.* [1], but **quantitative
induced-moment refinement is an additional, currently unvalidated hypothesis**.

## 9. Acceptance gates for physics development

**R0 (implemented research preview)**: periodic bulk Fourier `M -> A,B`
with correct SI units and algebraic invariants; optional scalar `A_z`
split-step; existing QCBED-grating compositional adapter; NiO toy magnetic
cell; finite-difference Fisher diagnostics. No field-induced net moment claim.

**R1 (required before a published magnetic CBED simulation)**: match
published magnetic Bragg intensities/symmetries to a benchmark with calibrated
magnetic form factors, correct magnetic reciprocal cell and converged mesh,
slice/thickness/zone axes; investigate neglected `A_perp` and spin terms.

**R2 (required before field-dependent refinement)**: finite foil and vacuum,
external plus demagnetizing field (`q=0`), diffraction optics/beam deflection,
and a field-calibrated magnetic state model; joint field-series objective.

**R3 (required before scientific parameter claims)**: Fisher/profile
identifiability and synthetic recovery under realistic noise, systematic
uncertainty and domain mixtures; ideally independent sample/magnetometry
constraints; real data and held-out measurements.

Do not expose an unvalidated field-induced moment in the production GUI or
serialize research coefficients as interchangeable with elastic `E_g^(s)`.

## References

1. O. Zaiets *et al.*, “Determination of magnetic symmetries by convergent
   beam electron diffraction,” *Physical Review Research* **8**, 013094
   (2026). <https://doi.org/10.1103/hf8x-79ry>.
2. A. Edström, A. Lubk, and J. Rusz, “Magnetic effects in the paraxial
   regime of elastic electron scattering,” *Physical Review B* **94**,
   174414 (2016). <https://doi.org/10.1103/PhysRevB.94.174414>.
3. K. Lyon and J. Rusz, “Parameterization of magnetic vector potentials
   and fields for efficient multislice calculations of elastic electron
   scattering,” *Acta Crystallographica A* **77**, 509–518 (2021).
   <https://doi.org/10.1107/S2053273321008792>.
4. J. C. Loudon, “Antiferromagnetism in NiO observed by transmission
   electron diffraction,” *Physical Review Letters* **109**, 267204
   (2012). <https://doi.org/10.1103/PhysRevLett.109.267204>.
5. “Simulations of magnetic Bragg scattering in transmission electron
   microscopy,” *Ultramicroscopy* (2023),
   <https://www.sciencedirect.com/science/article/pii/S0304399123000153>.