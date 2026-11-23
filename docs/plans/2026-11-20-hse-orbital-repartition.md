# Fast-PR HSE orbital rank-count migration

Opt-in `NativeGSBackend(..., hse_orbital=True)` makes the full-support Gamma
orbital GS algorithm explicit even for one rank. It records the existing Gamma
exchange selector in input text and requires HSE-orbital protocol2. Default k
and existing same-layout profiles keep their arithmetic. Source and target use
the same explicit profile and input except the generated parallel block.

Orbital-repartition protocol2 adds dense HSE. Fortran maps ACE factor columns,
source/action/work columns and spatial source/previous columns by global orbital
index. Replicated gauge and occupation state stays intact. Preserve mixing and
SCF/phase state, bind communicator IDs, stage before publication; no SCF/ACE
reconstruction. Packed/distributed ACE, tiled gauges, sparse localized support,
RT factor histories and changed-runtime migration are not this dense profile.

Verify Gamma1↔2 and PBE-stage cuts, exact cut orbitals/zero bootstrap and
continued source-initialized references at existing1e-9/1e-10; unequal orbital
partitions if capacity permits. Same-layout HSE, previous k-HSE/ordinary parity,
protected independent owners, source/installed suites, sample and wheel. Update
JP/EN/AI, native header and active checklist. No new modules or routine CI.
