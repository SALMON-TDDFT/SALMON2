# SOC exchange helper integration

Goal: promote the existing validated native SOC r/rprime worker candidate onto c2b85dff, preserving all newer KY, periodic indexing and PZHEEV guards. User explicitly authorized selecting this newer MPI implementation and SOC automatic pair selection.

Architecture: all ranks parse/broadcast normal namelist input before splitting. Existing explicit nproc_k determines SCF owner ranks for supported static SOC HSE06 full-q k-only runs; extra launch ranks service exchange only, without local SCF/ACE orbitals. Default nproc_k=0 retains ordinary decomposition. Bounded two-axis worker kernel comes from work/soc-native-workers-20261021 (final scoped revision); owner SCF/ACE keeps its communicator. No active binaries/jobs changed. No old source tree is built.

Tasks:
1. Port scoped worker/service modules and replace environment tuning with explicit configuration arguments. Preserve final numerical and low-idle-CPU fixes; MPI/HSE guards, line limits and explicit imports.
2. Add post-input communicator split, owner/helper entry lifecycle, namelist worker divisor/groups (safe existing defaults16/1), and no-helper serial/HSE-off paths. Validate unsupported inputs collectively before service entry.
3. Integrate native exchange adapter using latest candidate SOC symmetry routines without replacing newer unrelated fixes. Make eligible SOC auto choose full pair; semilocal/unsupported auto remains ordinary and explicit unsupported modes reject. Preserve explicit reference and old overrides.
4. RED/GREEN selector and communicator/worker tests; source ban gate; fresh serial/MPI HSE builds and all linked objects ban gate. Native same-input 4owner vs4owner+helpers SCF equivalence, failure collective guards; no GS payload collection. Review scientific numerical scope honestly; older Fugaku compiler/layout discrepancy remains recorded, not claimed resolved.
5. Update docs, code review, final verification and commit reviewed code. No automatic production deployment.
