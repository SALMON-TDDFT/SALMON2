# Stored trajectory optical diagnostics (2026-11-21)
Approved optical queue: duration/window sensitivity and post-hoc k contributions.
Extend existing analysis.py and convergence_kernel.f90, no new Python component.
Native explicit-frequency transform (none/hann/causal cubic), common requested
frequencies for all prefixes; no interpolation/zero padding/resampling or solver.
Finite uniform samples, positive/increasing frequencies including optional DC,
Nyquist masking. Large N*F*M costs explicit; no implicit full-grid quadratic work.
Sensitivity returns spectra, counts/durations/resolution and differences versus
longest prefix with first requested window; no automatic convergence certificate.
k reduction consumes existing weighted DecomposedCurrent, explicit charge/matter
sign and AU, dt/time0; native sums band/spin entries into k labels with no extra
weights, exposes total and optional comparison to simultaneous total charge
current. Optional native Fourier k spectra. Missing decomposition is not inferred,
partial bands labeled partial; never declare a sampling convergence result.
Tests: native Fourier vs independent analytic/geometric sum, prefix/window effects,
zero/weak drive and masks; weighted cancellation/sign/total consistency, invalid
NaN/Inf/time/index rejection and immutable independent outputs. Source-to-wheel
and synthetic usage example, paired guides/AI JSON; no GS or physical CI.

Implementation is connected to existing analysis.py and convergence_kernel.f90.
Two consolidated regression cases cover analytic transforms, actual prefix/window
differences, cancellation/sign and invalid inputs. Source suite665 tests passed,
one optional skip. Fresh sdist/wheel and stored-file usage validation follow.
No solver/source inventory modification or new GS data.

Final verification: source and fresh installed wheel each completed665 tests
(one optional skip, zero failures). Final native comparison normalization also
runs in Fortran. Saved-NPZ CLI and direct API tested from installed wheel; a
subprocess trap confirms no process launch, and MPI remains uninitialized.
Weighted synthetic matter-current decomposition reconstructs physical charge
current bitwise; source NPZ SHA256 unchanged. Existing optical_linearity accepts
new optical spectra even when optional peak metadata is unavailable.
Paired guide/AI JSON synchronized. Canonical record/artifacts:
calc2026/salmonpy-trajectory-diagnostics-20261121/record.json.
No GS, remote CI or new real-material optical-accuracy claim.
