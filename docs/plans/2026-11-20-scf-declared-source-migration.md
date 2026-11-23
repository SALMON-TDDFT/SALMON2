# Declared SCF source-version transitions, without conversion

Keep strict and same-source modes. Add compatible_runtime='declared': accept
only an explicit packaged source-contract pair plus equal codec/type inventories,
native stream/platform proof and declared numerical-state semantics. Matching
byte layout alone is not blanket semantic compatibility. The catalog is a
reviewed allow-list with transition ID and evidence; never infer arbitrary
source-version compatibility. Codec fingerprints include complete numerical type
source files and both generated IO routines. Generator --check checks inventory.

Native state revision1 is added without changing stream1, serialized fields or
solver arithmetic. Baseline is the previous audited same-source rebuild runtime;
the current source differs by this metadata export. Validate ordinary/HSE, full/
reduced-k mixed rank changes, zero SCF replay and continuous agreement; strict/
same-source refuse the changed source, unknown pairs/revision/codec mismatch
refuse before creation. Register only validated exact pairs, include IDs in
state/save provenance. Same-source v1 checkpoints remain valid. No conversion
or cross-platform claims; no routine CI. JP/EN/AI/catalog/wheel/sample parity.

## Accepted local checks (2026-11-20)

- Only source difference from the previous audited prepared tree is the four-line
  state-revision export in src/gs/scf_disk_native.inc. Codec fingerprints match;
  real binary/source digests differ. Exact ordinary/HSE pairs are registered.
- Source tests637:636 pass,1 skip. Installed-wheel tests637:636 pass,1 skip.
- Source ordinary same-layout/4→2 and HSE same-layout/2→8/reduced8→4/pre-PBE2→8
  pass exact cut,bootstrap0 and unchanged continuation rtol1e-9/atol1e-10.
- Source/stream/platform/codec/revision/missing proof and default/same-source
  refusal controls pass before owner creation. Corrupt/late-read controls preserve
  originals. Previous normal ordinary and full/reduced HSE arrays are bit exact.
- Installed ordinary/HSE and public reduced-k8→4 example pass. Historical
  scf-rebuild/1 archives are checked using independent copies and matching input
  profiles; no source archive metadata is edited.
- Evidence: calc2026/salmonpy-scf-declared-record-20261120. Campaign-specific GS
  storage was authorized. No CI dispatched; cloud synchronization not confirmed.
- This validates a metadata-only source transition. Arithmetic changes, changed
  numerical meaning, cross-platform migration and format conversion remain open.
