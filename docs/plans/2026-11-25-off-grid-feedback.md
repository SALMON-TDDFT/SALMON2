# Off-grid molecular feedback

Approved scope: extend the existing dipole feedback first, using reciprocal trilinear weights at macro sample sites. Fortran builds the eight-site map, samples the molecular drive, and deposits current with inverse physical cell volumes. Existing Yee edge sampling/deposition then composes with this map. The own-vacuum field uses exactly the same map. Keep rollback, step numbering and checkpoint contracts.

Implementation sequence: failing affine/work identity tests; Fortran map and sampling; multi-site feedback and checkpoint integration; kernel build and oscillator MPI validation; Japanese/English guides, AI catalog and checklist. No new GS, physical reproduction CI or implicit extrapolation. Extended Iwasa spatial/magnetic response is not reduced to a dipole without an explicit approximation; its reciprocal full-current coupling remains separate work.

Completed: Fortran mapping, multi-site feedback/checkpoint integration, vacuum/LD oscillator tests, real C2H2 process-exit restart, old/new single-site comparison, final wheel92tests(90pass/2skip). See docs/off-grid-dipole-feedback-validation-2026-11-25.md. Full spatial/magnetic reciprocal response remains separate work.
