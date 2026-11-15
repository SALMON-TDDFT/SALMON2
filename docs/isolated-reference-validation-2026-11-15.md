# Isolated midpoint reference preparation

This is an F3 prerequisite, not an implementation or acceptance test of a new
E-driven TDDFT worker. The reference tool itself changes no native solver code
or public API; subsequent worker validation is recorded separately.

## Reproduce locally

With NumPy available, an explicitly supplied compatible SALMON executable, and a
new output directory outside this repository:

```sh
python tools/test_isolated_reference.py /path/to/salmon /path/to/new-reference
```

The tool uses bundled C2H2 exercise inputs and their C/H pseudopotentials, retains
source files unchanged, and records executable, pseudopotential and generated
input SHA256 hashes. It uses one process and one thread, with a 180-second timeout
per invocation. It requires SCF convergence and two identical 20-step midpoint
runs from separate copies of the same GS. It checks finite output, sample count,
time grid, nonzero E/dipole response and forward-difference E consistency.
No download, build, MPI launch or CI workflow is performed. This is not a hard
memory limit or a physical convergence study.

## Local measurement

Apple Silicon, existing explicit runtime executable SHA256:
`29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed`.

- C2H2 PZ, 24 x 24 x 32 grid, spacing 0.35 Angstrom.
- SCF converged at iteration 149, residual 7.4018414e-10.
- Midpoint: 20 steps, dt=0.001 fs; original exercise pulse retained.
- Two independent batch runs: dipole/field and energy text arrays exactly equal.
- E_tot versus forward difference of output A: maximum absolute discrepancy
  1.553978308879711e-10 V/Angstrom (output conversion/rounding included).
- Combined GS and two RT runs: about 3.8 seconds, one thread. RSS not measured.

The short, coarse calculation establishes reproducible wiring data only. It does
not establish converged optical response or independent interactive worker states.

## Important interface finding

`initialization_rt.f90` prepares E(n)=-(A(n+1)-A(n))/dt. The isolated midpoint
path uses `calcVbox(..., itt, ...)`; AETRS additionally uses itt+1 and requires a
separate input contract. The new worker must preserve these sample choices.

In this reference, the file's E_ext columns are zero, while E_tot contains the E
used in propagation. `calcVbox` resets `system%vec_E_ext` without assigning it for
the ordinary isolated branch; `write.f90` writes that field. The periodic branch
assigns it separately. Consequently the E_ext header is insufficient to identify
the consumed field. The regression uses E_tot and checks it against A. This audit
records the behavior; it does not silently change native output or dynamics.

Remaining: versioned E-worker capability/protocol, isolated state loader,
explicit dipole observables, interleaved state isolation and replay tests, then
comparison with this reference. Maxwell current-density normalization remains a
separate explicit contract.
