# Non-LD workflow validation

Implementation: 0.1.0a2, code commit 045f34faef532ceb106ea00274319a70d27eca8d.
CI: https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43357199430

Implemented: isolated FDTD batch input/execution/observation reading, box and
ellipsoid material-grid preview with boundary clearance, Acos2 pulse composition,
SALMON text table reading, explicitly normalized spectra, and the common A-only
external-backend protocol/orchestration. Built-in TDDFT/SBE adapters remain absent.

Local integration uses a previously source/object-audited GNU serial SALMON binary:
SHA256 36ddc1b5d6527c89f0718be0e81cec0ae9fb2fe86fc5632f0b0551e961086c79.
`tools/test_fdtd_entrypoint.py` compares two physical inputs with independently
handwritten SALMON namelists, run through the same executable at one OpenMP thread:
- Dielectric ellipsoid, 32x32x32 grid, 200 steps.
- Dielectric box, 33x32x34 grid, 200 steps (odd/even coordinate conventions).
Both observation tables and shape.cube files match byte-for-byte. Python rasterized
material IDs match the SALMON cube values exactly. An LD impulse case also ran
successfully (200 samples), separate from these dielectric regressions.

Fixed findings:
- Serialize al_em + dl_em: num_rgrid_em + dl_em hits a legacy global-input check.
- Reject dielectric + impulse; use plane-wave Acos2 for constant dielectrics.
- Keep material clear of PML interface coefficient planes as well as layer cells.
- Revalidate/copy directly supplied LD configs into immutable tuples.
- Map common TDDFT Acos2 to FDTD CEP with +0.25 cycles; independent derivative test.
- Bound accumulated backend-time roundoff without accepting an incorrect timestep.

The wheel suite tests input/ownership, scalar rejection, independence, LD physical
response, numerical convergence, shape semantics, Fourier phase/normalization,
Fortran exponent parsing, pulse convention, FDTD serialization/process failure,
and A-only backend contracts. It does not run a full SALMON FDTD binary on CI;
that integration is the separate local check above. No external GS data collected.

Final CI result: all five jobs succeeded. Windows AMD64, Linux x86_64, macOS
Intel and macOS arm64 each passed CPython 3.10–3.14 installed-wheel tests:
27 tests per wheel, 540 checks across 20 wheels. Source-distribution rebuilding
and wheel-record tests also succeeded. Artifacts are attached to the linked run.
