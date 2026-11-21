# Status-bearing owner transport

The managed, unanimously selected coupled-response path carries group status with A/packet delivery and J return. Standalone/custom-dispatch paths retain the existing protocol.

A: subgroup preflight gathers local failures; leaders exchange one status per source/destination together with only the field source's destination-owned values. Numeric A uses typed Alltoallv (one float status header per source), packets use object alltoall. Each subgroup receives a combined status/value message. This preserves whole-job pre-callback rejection. Only on failure, broadcast the explanatory message once from the globally selected reporter after every rank has learned its rank code from the transport headers. Recipient packet validation remains a distinct pre-callback vote.

J: group callback/reduction votes already establish local outcome. Leaders gather each group's status header and currents together. A short leader outcome notification remains necessary so remote groups know whether to commit or roll back. Field subgroups receive status and their spatially owned J in the same Scatterv; material subgroups receive the outcome. No J/current or clock is committed on failure. The context owns one duplicated failure-only communicator, independent of caller lifetime, and frees it at close. Normal steps never use it. Error descriptions use one reporter-rooted Bcast only on failure.

Keep the pre-field clock vote, pre-H sample vote, native numerical/physical reductions, and protected capture/rollback/cleanup votes. Successful separate hierarchical outer votes become A3→1, packets without H4→2, with H5→3. Header/data transports and final acknowledgements still communicate; measure actual Python MPI operations as well, do not equate this count with all MPI traffic or speedup. Status Alltoallv has O(groups) headers per leader (O(groups²) total), like replicated leader status distribution; multi-node scaling is unproven.

Validate failed count first, then compare preceding wheel and current source: normal fields/current, direct restart, protected late/deep failure and retry, single-rank invalid fields/payloads/receipts, packet H ordering, custom backend and custom public-response fallback. Preserve all old numerical/order/checkpoint contracts. No new GS, native solver rebuild or manual CI. Update Japanese/English guides, AI JSON and checklist.

## Verified completion

- Normal MPI4 A transport now has one separate hierarchical vote and four Python leader collectives, versus three/five in the preceding installed wheel. Packet no-H counts are two/five versus four/six; H counts are three/six versus five/seven. These are scoped call counts, not full native MPI or speedup measures.
- Normal reporter-domain Bcasts: zero. Each supported A/J failure issues exactly one explanatory Bcast from the rank selected in headers. Native/SBE protected traces cover seven late/deep failures with reporters1,2,0. The necessary pre-H failure decision retains its old hierarchical route before data transport.
- Dedicated error communicator is owned by ParallelContext. MPI3 test frees the caller's parent before use, uses a middle-position field group, successfully returns a backend error and leaves caller MPI alive after context close.
- MPI4 native FDTD+two field-local LD materials+two SBE points: all state bytes, aliases and output files preserved after supported failures; protected retry and same-runtime direct checkpoint resume match continuous execution with zero replay. MPI6 protected LD/SBE/custom, standalone packets and independent response overrides retain their contracts.
- Seventeen A/E/J comparisons are bitwise equal to preceding controls. Final source and installed-wheel suites run668 tests each, with one optional skip and no failures. Installed-wheel MPI tests repeat protected failure/retry, resume, H/no-H counts, reporter Bcasts, parent lifetime and custom callbacks.
- Japanese/English guides, AI JSON and checklist synchronized. No new GS/TDDFT/HSE material-physics run, no full-runtime rebuild and no manual heavy CI dispatch. Evidence and wheel/sdist: canonical `calc2026/salmonpy-status-transport-20261121/record.json`. A diagnostic output-path collision and initial failing contracts are retained/documented there.
- Separately repaired the reported CI failure: Python3.12 throws TypeError for a frozen slotted dataclass property assignment where Python3.14 throws AttributeError. The test now accepts either rejection and checks the unchanged value. Test-only commit8d5d90f5 was pushed; lightweight run44852912505 passed, platform matrix skipped. Communication changes remain uncommitted.

Manual transport-count probe (ordinary sample arguments):

```sh
mpiexec -n 4 python tools/test_outer_agreement_mpi.py --expected-votes 1 --expected-leader-ops 4 --report NEW_REPORT.json BUILD RUNTIME BUNDLE NEW_OUTPUT --layers SBE_ONLY_LAYOUT --field-grid 2 1 1 --tddft-ranks 1 --steps 96
```

The report includes successful/failed steps, scoped leader operations and reporter-domain Bcast roots. Use `--driver SCRIPT` for an explicit alternative contract driver and `--failures-only` for rejection-only probes. Broken MPI/process loss and multi-node scaling remain outside this evidence.
