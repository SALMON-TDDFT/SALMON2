# Wannier90 post-localization covariance validation

## Problem

The production DMN supplies a block representation of the affine generators.
Within a spectral basin block, a generator may mix internal channels densely.
Wannier90 can therefore preserve the supplied symmetry without mapping each
individual Wannier centre to exactly one other centre.  The current immediate
post-Wannier centre-orbit gate incorrectly assumes a monomial representation
and rejects the converged Si64 result.

## Design

Replace the immediate centre-orbit gate with a direct generator covariance
measurement.  For every streamed affine generator, measure

`U^H D_band(g) U - D_wann(g)`

where `U` is the returned Wannier90 transform, `D_band` is the retained-basis
generator action, and `D_wann` is the same action expressed in the spectral
trial gauge written to the DMN.  Reject unless the maximum defect over all
generators is within the accepted symmetry tolerance.

The validation must process one generator at a time.  It may allocate small
dense work matrices on the coordinator, but must not retain all generator
matrices or all character sectors.  The maximum defect, checked workspace
receipt, and deterministic provenance fingerprint are retained for downstream
checkpoint binding.

The later post-character periodic-position localization and its individual
centre-orbit validation remain authoritative.  That stage intentionally fixes
the internal block gauge and is therefore the correct place to require
monomial centre behavior.

## Testing

- A focused MPI fixture accepts a dense internal block representation when the
  transform satisfies covariance.
- The fixture rejects a unitary transform that violates covariance.
- The route contract requires the covariance gate before downstream W90
  anchoring and forbids the premature centre-orbit gate.
- Existing MPI 1/2/4/8 W90 and construction fixtures remain green.
- Si64 must pass the new immediate covariance gate and proceed to the later
  canonical centre validation.
