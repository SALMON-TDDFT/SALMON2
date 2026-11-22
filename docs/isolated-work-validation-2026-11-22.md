# Isolated TDDFT internal energy / dipole work

Measured-electron-number work committed/pushed as416b1a83 to independent SALMONpy.
This follow-up is uncommitted. Reuse the explicitly authorized C2H2 GS; no new SCF,
no manual heavy CI. Results: canonical calc2026/salmonpy-c2h2-feedback-20261122.

## Observation contract

E3 returns six dipole values, electron number, field-free internal energy and a
0/1 freshness flag. Preserve the existing calc_eigen_energy condition (step1,
capacity, multiples of out_rt_energy_step); initial energy is also fresh. Do not
force energy evaluations or introduce a reduction. For isolated fixed-ion PZ,
update_vlocal removes the external Vbox before calc_eigen_energy and total-energy
assembly, so energy%E_tot excludes external dipole interaction. Public energy is
None when not fresh. Revision3 evidence/checkpoints reject older forms.

Fortran accumulates held-endpoint work ΣE[n]·Δp and midpoint work
Σ(E[n]+E[n-1])/2·Δp separately. This worker actually consumes held endpoint E;
midpoint is a distinct field-work convention. Compare energy only at fresh times.
Fresh data are explicitly indexed; stale raw energies cannot be treated as a
physical time series. Finite residuals are diagnostic; no stop threshold added.
No automatic subtraction of the independent zero-field energy drift.

## Actual measurements

PZ C2H2,24×24×32,0.35Å; outer MPI4,one field rank per context,64 field steps,
dt=.003au. Field time=electronic time+dt. Last fresh energy at electronic step60
(time.18au), rather than the final electronic step63. Old E/p/j/Ne histories are
bit-identical. E3/batch reference dipole difference3.8553e-12e·bohr; internal energy
at common output times differs at most2.13163e-14Hartree.

| Quantity (Hartree) | default amplitude.001 | explicit signal control.1 |
|---|---:|---:|
| held-endpoint work |6.98380840377e-10|6.98427541762e-6|
| midpoint work |9.85414518807e-10|9.85425956741e-6|
| internal-energy change |5.76303307298e-8|7.04692036635e-6|
| held-work residual |5.69319498894e-8|6.26449487250e-8|
| midpoint residual |5.66449162110e-8|-2.80733920106e-6|
| zero-field drift (not subtracted) |5.74708671763e-8|5.74708671763e-8|

The default signal is below background drift and does not certify energy accuracy.
The signal control's held-work relative residual is0.00896943; midpoint residual
relative to midpoint work is0.284886. Thus full coupled energy conservation is
NOT established. The next task is the molecular/field time-placement agreement,
then full Maxwell/LD/boundary energy balance, while preserving baseline evidence.
Do not hide the mismatch by subtracting drift or changing the propagator.

## Implementation verification

Analytic held-work/freshness and reversed-sign controls pass. Scalar IEEE NaN/Inf
and overflow rejection leave output buffers unchanged. Snapshot/replay includes
energy and freshness. Actual memory/disk replay, interleaved independent workers
and freshness/None semantics pass. Actual source continuous versus wheel injected
rejection at interval30/retry: all retained energy/work/field/dipole/norm arrays
exactly equal. This proves isolation/retry, not physical conservation.

Source/wheel686 tests each, one optional skip each. Source/object prohibition
audits and Git-only runtime preparation/build pass. New-source GS-free NLCC XC
and core-force derivatives MPI1/2×OMP1/2 pass before source identity advancement;
no new whole-application force receipt. Paired guides/AI JSON synchronized.
