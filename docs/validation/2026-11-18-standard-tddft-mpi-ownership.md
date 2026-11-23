# Standard embedded TDDFT initialization MPI ownership

## Change and boundary

The standard communication module now uses the same initialization-child owner
tracking verified for fast-PR HSE. Scope covers embedded TDDFT initialization and
its initial-current evaluation. Close discards a pending numeric restore candidate
and protected payload, frees request caches and XC, then retires only this state's
children. The host parent and unscoped children keep their existing ownership.
No propagation arithmetic/update order or checkpoint encoding was changed.
One native lifetime/process remains enforced; no standard multi-context or
heterogeneous reentrancy, initialization/MPI/allocation-failure recovery claim.

Optional raw ABI: sp_tddft_context_mpi_protocol()==1 and sp_tddft_resources(int[4])
(live tracked children, registry slots, pending candidate, allocated state).
Diagnostics do not measure total MPI or process memory. Vacant registry slots
remain until detach. Original source headers/provenance and fixed binaries retained;
standard tree inventory and the dependent experimental step-config base pin updated.
Shared LD/SBE kernels remain byte-equal to native/.

## Local verification

macOS26 Apple Silicon, Python3.14, Open MPI, GNU15, Libxc4.3.4; OMP/BLAS1;
Mac64GiB/18 cores. Fresh external source/build and source/object prohibition audits.
The MPI C wrapper's GNU selection (OMPI_CC) is recorded with the build command;
initial Clang/OpenMP discovery failure is preserved separately. No Actions run.

- Old runtime fails the new probe's explicit missing-capability check.
- Five XC inputs: PZ, Libxc-PZ, Libxc-PBE, TB-mBJ automatic and fixed-c; each with
  two orbital or two real-space ranks.20 MPI jobs compare old/new,10 pairs.
  Every rank's initial plus32 propagated charge currents match exactly.
- Before init: zero tracked resources/payload/state. Orbital split creates17
  tracked children; real-space split creates15. Staging a direct checkpoint adds
  one pending candidate. Close returns live children, pending and state to0,
  leaves the host collective usable, and a second init still returns status12.
  The PZ fixture explicitly changes only XC from the synthetic PBE fixture;
  inputs/hashes and transformation are retained, not a physical GS certificate.
- Four ordinary protected-step regressions (PBE, Taylor HSE, PT-CN HSE and PBE/HSE
  combined),32 steps, two field ranks, match preceding numerical baselines.
- Normal full-k mixed FDTD/two-LD/SBE/PBE/shared-HSE128-step run with6 ranks matches
  preceding A/E/j. Save at63 and direct fresh resume executes65 new steps with no
  replay. Handled last-HSE candidate staging failure preserves all field/LD/SBE/
  electronic originals and output streams; retry matches continuous propagation.
- Source and isolated installed-wheel suites:603 tests each, one external-condition
  skip each. Source distribution→wheel/native checksum test passes. Installed
  TB-mBJ-auto real-space MPI resource probe passes. Installed normal mixed sample
  also completes continuous/save/failed-restore/retry/resume against its reference.
- Source preparation and experimental explicit-step-config preparation pass both
  from checkout and wheel. The first suite detected a stale dependent base pin;
  its hash was updated without relaxing source authentication. Read-only review
  found no blocker; paired guides/AI JSON and active checklist synchronized.

Evidence: developer outputs salmonpy-standard-mpi-20261118. Raw native checkpoints
are not canonical portable archives. Larger systems, other native environments,
unscoped/common MPI/I/O/timers and different-input state ownership remain open.
