# Shared EXX MLWF controls

HSE06 and PBEh40 use the same localization controls in `&functional`:

```fortran
  exx_mlwf_interval=5
  exx_mlwf_maxiter=100
  exx_mlwf_tolerance=1d-7
  exx_mlwf_radius=0d0
  exx_local_fft='auto'
```

The unchanged defaults are interval 10, maximum iterations 200, tolerance
1e-6 and radius 0. Interval counts changed-source exchange refreshes, not ionic
or electronic time steps. Tolerance is always in atomic units (bohr squared).
The three old names `hse_mlwf_interval`, `hse_mlwf_maxiter`, and
`hse_mlwf_tolerance` remain compatible inputs. Supplying both names with different
values is an error; matching values are accepted. Logs use the canonical names.

## Meaning of the support radius

`exx_mlwf_radius` uses the length unit selected by `unit_system` (bohr in `a.u.`,
angstrom in `A_eV_fs`); it is logged in bohr. Zero keeps full support. A positive
radius applies the existing LCFO periodic sphere-mask geometry to each
occupation-weighted Wannier source `Q=Psi sqrt(f/2) U`, about its periodic center
in the Born–von Karman supercell (fragment supercell for DC). Both source
appearances in the Fock action use this mask. Targets are not masked and sources
are not renormalized. Exactly zero pair densities are skipped, and compact
support accelerates pair formation and output accumulation through the existing
Wannier backend. With `exx_local_fft='auto'` (default), compact sources use a padded local-box
FFT whenever its volume is smaller than the global grid. `'off'` retains the
full-grid FFT for comparison. The local kernel is sampled from the inverse FFT
of the original global multiplier, including G=0; it is not a new local periodic
Coulomb model. Padding at least 2*m-1 along every occupied box axis makes the
restricted convolution equal to the global discrete action, up to roundoff.
The source mask remains the only additional approximation in this comparison.

Periodic boxes are unwrapped across the largest empty gap on each axis. Broad
sources fall back to the global FFT. The local path currently runs serially per
pair and caches one box shape; the global path keeps its OpenMP/batched FFTs.
`EXX_FFT` logs local/global pair counts and actual/full pair-grid points. Those
counts exclude kernel/plan setup, source gathering and other SCF work. Global
source/target arrays and an O(G) periodic kernel are still stored. Dense target
states can still overlap every source: this does not establish overall linear
scaling, nor does it enable the missing PBEh optical/DC-MD adapters.

Centers use circular density moments. As in the existing LCFO path, a source
with normalized moment below 0.1 on any axis stays uncut because its center is
ambiguous. Zero-norm sources are also protected. `EXX_MLWF` output reports the
radius, protected count, total discarded norm fraction, and maximum discarded
fraction of any source at localization updates. A large radius reproduces full
support. Compare several radii against radius 0 for energy, density and spectrum;
a small norm loss alone is not an exchange-energy error bound.

This is an explicit, gauge-dependent source approximation. It preserves the
Hermiticity of the exchange action and reuses ACE, but it does not establish a
variational energy functional or consistent ionic forces. Positive radius is
currently limited to static DFT (including static DC), with MD, relaxation,
real-time propagation, restart and legacy snapshot export rejected. Any force
values printed by a static finite-radius run are not validated for use in MD or
optimization. Unconverged localization is reported during iteration and does not silently
disable the requested mask. If a source is actually truncated and the last
spread minimization was not converged, the final SCF result is rejected, even
when the density criterion passed. Increase localization iterations or use full
support; density convergence alone cannot certify a gauge-dependent cutoff.
Finite-radius initialization uses physical k indices rather than k-rank layout
for random seeds. This removes an avoidable MPI dependence of the initial state,
but does not guarantee a unique localization minimum or radius convergence.

The separate `pbeh_coulomb_radius` truncates the Coulomb **interaction kernel**,
not orbitals. Changing it is a different approximation. The existing LCFO-RT
`hse_lcfo_wf_radius` remains a separate control in bohr; it is not an alias for
`exx_mlwf_radius` and is unaffected by this SCF extension.

## Adaptive support on the spatial mesh

For Gamma-point native mesh calculations, use:

```fortran
&functional
  exx_mlwf_norm_fraction=0.999d0
  exx_mlwf_radius=0d0
  exx_mlwf_interval=5
  exx_mlwf_maxiter=100
  exx_mlwf_tolerance=1d-7
  exx_local_fft='auto'
/
```

`exx_mlwf_norm_fraction` defaults to 0 (disabled); 0.999 retains at least
99.9% of each occupation-weighted source's squared norm. Each periodic sphere
has its own radius, determined from local grid rows and small reductions.
Sources with ambiguous circular centers remain uncut. There is no source
renormalization and no cutoff on arbitrary target wavefunctions. A value of
1 is the full-support reference through the same spatial backend. A positive
fixed `exx_mlwf_radius` and positive norm fraction cannot be combined.

The masked source appears in both factors of the exchange action. Compact
convolution uses the same discrete global kernel (including G=0), but stores
only its distributed grid rows and gathers compact displacement/source/target
tiles. Target work is batched in four columns; no full wavefunction grid is
gathered. Broad sources fall back to pencil FFTs. `exx_local_fft='off'` is the
full-grid convolution reference for the **same masked sources**.
`EXX_ADAPTIVE` reports retained fraction, maximum radius and norm loss, and
local/global pair counts with local FFT point counts for orbital group zero.
Only exactly zero pair densities are skipped. Dense ACE/gauge algebra and
source-target pair work remain; this is not a claim of linear overall scaling.

In SCF, full support is used until the selected convergence residual is below
`sqrt(threshold)`. Masking then requires converged MLWF localization. This
allows occupied-only calculations without adding empty states just to measure
a gap. Mode changes reset the density mixing history and invalidate that
iteration's convergence result. Failed gauge transport or localization restores
full support; a result without established adaptive localization is rejected.
In DC, all source geometry and exchange convolutions use the buffered fragment
cell and its communicators. SCF switching readiness is synchronized across
fragments so global mixing history is reset consistently.

Adaptive support is also admitted for fixed-ion, DC-initialized native mesh
RT with the existing Taylor4/ACE predictor-corrector. RT starts localization
immediately, without a gap or electronic-temperature criterion, and recomputes
support on exchange refreshes. This is distinct from propagating in an LCFO
basis. Moving ions, ionic relaxation, retained-LCFO RT, restart and snapshot
export remain rejected for adaptive support. The fixed-radius legacy option
keeps its previous static-only contract. As with existing finite support,
energy/force variational consistency and long-time optical accuracy are not
established by norm retention alone.
