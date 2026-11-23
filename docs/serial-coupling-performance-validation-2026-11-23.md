# Serial field/SBE performance mode — 2026-11-23

`execution_mode="performance"` now propagates through `maxwell.couple` to the
field and every SBE state. Four field components support the same option for
standalone stepping. Setup contracts, input/CFL/capacity and explicit output /
snapshot/restore checks remain. Per-step finite scans and clock comparisons
are omitted. Update expressions and ordering are unchanged; no new Python
numerical algorithm was introduced. Failed performance updates have no automatic
rollback guarantee. TDDFT/HSE and MPI-session performance mode remain pending.

Nonzero multi-component drive: all four fields with SBE agree exactly with
checked mode for A/current; 1D has two independent points. Direct field snapshot
restoration and continuation agree. A nonfinite trial propagates in performance
mode as intended. Unsupported custom opt-in is reported; extension-backed
factories need not expose Python introspection metadata.

The final sdist-built wheel was installed separately. Its field coupling example
ran transverse, longitudinal and 3-point 1D propagation for 4000 steps each in
both modes: all time/A/current arrays were exactly equal and current was nonzero.
Source and installed wheel each ran 700 tests: 699 passed, one optional skip.
Wheel code, sample, paired field guides and AI JSON match source bytes. Build
retains mandatory source/object eigensolver checks. No GS data or remote CI.
This is numerical/ownership validation, not physical convergence or speedup.

Records: calc2026/salmonpy-serial-coupling-performance-20261123.
Local OneDrive storage does not confirm cloud synchronization.
Wheel SHA256: `e7e693a9ae696ec2a1efb37be26a6147bc262609f9b620b40c1026fffc58a93f`.
