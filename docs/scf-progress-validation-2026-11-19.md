# Standard SCF progress ownership verification

The new s_scf_iteration_progress is owned by s_mixing. Previous density/local
potential, electron normalization, no-mixing convergence counter and convergence
flag no longer disappear as SCF stack locals. Exact final iteration/residual and
ion-update values are retained too. Fresh-entry initialization still resets only
this context; pre-existing mixing rates and histories are untouched by that reset.
There is no new public resume API, C signature, checkpoint format or continuation
entry. Band/DC/molecular and fast-HSE continuation are outside this slice.

Verification:

- The native ownership probe first failed against the old module's missing type.
  With the new full build it preserves nontrivial bounds and independent copies,
  and repeated progress reset leaves the existing rate/density history unchanged.
- Actual before/after convergence routines are extracted, compiled and compared
  bitwise for all five criteria, spin/electron-count branches, zero/tiny/finite
  perturbations and MPI1/2 x OMP1/2. They perform no orbital solve or SCF call.
- The entire iteration executable body matches after reversing ownership names
  and excluding newly appended final scalar captures. This is an implementation
  equivalence check, not real-material SCF interrupted/continuous validation.
- Complete standard runtime compiled/linked with GNU/OpenMPI/Libxc on this Mac.
  Source/all application-object prohibition audits pass. No IEEE inquiries were
  added/changed. Fujitsu and other platforms remain untested.
- Synthetic PZ orbital2 and Libxc-PBE spatial2 on MPI4 preserve independent RT
  states and direct no-replay resume. Current hashes exactly match the previous
  runtime. Installed wheel's native MPI runner and SCF probe also pass.
-628 source and628 installed-package test methods pass, one skip each. No Python
  method added. Current optional RT-step profile's structures and manifest/base
  pins updated; source preparation/audits and offline sdist→wheel succeed.

One concurrent worker timeout test failed during startup; its isolated7 tests
and a full628 rerun passed. No timeout policy changed. Earlier missing GNU MPI C
backend configuration, stale optional-profile pin and incorrect installed-test
discovery are retained as diagnostic evidence, not hidden as physics failures.

The new previous-field arrays remain allocated with a retained GS until reset/
close (two local scalar fields). No resource budget is turned into a hard stop.
Full SCF continuation/complete checkpoint remain the next unfinished work.

Canonical records (no collected GS arrays):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-scf-progress-20261119

Local writes verified, cloud synchronization unconfirmed. No CI dispatched.
