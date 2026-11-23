# Ohmic field-energy observation — 2026-11-22

The Fortran observer reports nonnegative `ohmic_work` separately from LD polarization work, external-current work, six-face Poynting flux and constitutive boundary exchange. Both residuals and the summary include it. Propagation and checkpoint data are unchanged; no communication is added.

## Coefficient identity and step numbers

At each observed electric edge, let c1 be its actual electric multiplier, cJ its negative current coefficient and C2_q the actual implicit LD coefficients. For the accepted update from field step n to n+1:

    c4 = -4 pi dt / cJ
    epsilon = c4 (1+c1)/2
    sigma = -(1-c1)/cJ - sum_q C2_q
    E_mid = (E_n + E_(n+1))/2
    W_Ohmic = volume dt sum_edges sigma E_mid²

The LD implicit term must be removed; treating the whole attenuation coefficient as Ohmic loss double-counts polarization work. Only subtraction roundoff within64 double-precision spacings at1 times max(abs(total_loss),abs(LD_loss)) is set to zero. Material coefficients and propagation are not changed by that diagnostic treatment. Finite-value/negative-conductivity failures leave result buffers unchanged. All IEEE inquiries are scalar.

The current field observation is float64[19], with constitutive exchange at18 and Ohmic work at19 (one-based). The summary consumes(19,n) and returns12 components, with separately summed constitutive exchange and Ohmic work at11/12. The standalone conductive kernel consumes12 coefficient pointers and returns its17 base values plus Ohmic work. Explicitly rebuild field and lightweight libraries together; no old-format fallback is added. Zero-conductivity configurations skip the additional evaluation and preserve previous values, adding only ohmic_work=0.

## Evidence

An independent manufactured uniform-field step with directional conductivity and nonzero implicit LD coefficients verifies positive analytic Joule heat, LD work separation, zero-conductivity control, spatial partition additivity and NaN/Inf/negative-coefficient rejection with unchanged output. The existing summary check now verifies separate constitutive and Ohmic terms.

Field fixtures:40³ grid, spacing2 au, dt0.003 au,64 steps; sigma0.2 au. Conductive dielectric epsilon2.5/mu1.4; a two-body conductor+LD fixture adds background epsilon4, plasma frequency20 and damping2; a unit epsilon/mu conductor tests the dispatch condition independently of dielectric weighting.

| Profile | Max absolute modified residual (Hartree) | Total Ohmic work (Hartree) |
|---|---:|---:|
| Conductive dielectric |3.676269358376949e-14|0.0016279559076623148|
| Conductor + LD |1.865174681370263e-14|0.0008482573727115383|
| Unit background conductor |4.462988138775881e-14|0.004436447909295241|

- Mixed fixture: observation/no observation and previous runtime final per-rank E/H, LD, PML and current checkpoint bytes match exactly.
- Save at step32, terminate process and directly resume in a separate MPI job to64: all histories/budgets and final state bytes match continuous execution exactly.
- Field MPI1/2 and x/y/z partitions agree at rtol2e-12, atol2e-13.
- Zero-conductivity legacy, reciprocal tensor background and subpixel controls preserve all previous observations and state bytes; Ohmic work is exactly zero.
- Source and installed wheel suites:695 tests each, one optional skip each. Final wheel continuous/restarted histories and state bytes match source exactly.
- Current-source preparation and source/application-object prohibition audits pass. Wheel implementation, paired guides, AI JSON, inventory and changed solver files match the checkout byte-for-byte. Python -I loads the separately installed wheel; dependencies use the existing validation environment. This is not offline dependency-installation evidence.

These are bounded CPU discrete-balance checks, not physical convergence, arbitrary-material coverage or cross-platform validation. LD polarization work includes reversible storage and is not solely irreversible heat. No fresh GS/SCF or optional remote physical CI was run.

Canonical records: calc2026/salmonpy-ohmic-field-energy-20261122/summary.json and verified-artifacts/. Individual field cases remain in the existing field-energy campaign identified by the summary. Local OneDrive storage is confirmed, cloud synchronization was not checked.
