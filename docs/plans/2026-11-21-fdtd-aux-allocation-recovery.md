# Material and boundary coefficient vector recovery

Bounded R3 continuation; PML bulk split-field coefficient arrays already have
1118 protection. Keep small automatic PML workspaces outside blanket OOM claims.

- Three material and six E/H boundary vectors use STAT and normalized1120.
  Skip remaining local allocations after failure and agree before use/communication.
- Coalesce the material gate with LD admission,retaining1119 and allocation/
  arithmetic order. One new field-group startup agreement with LD,two without
  LD;no new step/application-wide collective. Budget excess still warns.
- Existing testing-only selectors28/29 fail material slot2/boundary slot2 on
  one field rank after partial allocation. Production contains no injection hook.
  Existing1118/1119 regressions remain part of the consolidated fixture.

Final record: `salmonpy-fdtd-aux-allocation-recovery-v2-20261121/summary.json` in
canonical calc2026. GNU Mac MPI1/2 scalar/non-diagonal-tensor LD,source and wheel:
partial failure cleanup,retry,zero-replay direct continuation and prior-runtime
byte comparison passed.12 new controls match prior fields/observations/shape/
sidecars/checkpoints;both suites679 run,one optional skip,no failures. Native
source/object audits and Git-only preparation/build passed. No GS collection or
optional CI dispatch. Previous uncoalesced runs remain preserved as evidence.

Not covered:arbitrary tensor/subpixel/observation/scratch allocations,small
implicit automatic allocations,OS kill/physical OOM,broken MPI or process failure.
Next R3 work:bounded supported output-IO errors and selected auxiliary buffers.
