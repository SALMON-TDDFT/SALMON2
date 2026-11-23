# Model-owned parser workspace and defaults

Approved continuation of model-input ownership. Separate standard native parser
work areas/default records without enabling heterogeneous live propagation yet.

- An optional s_input_parser passed to read_input owns scratch units, coordinate
  block flags, namelist status, parsed values and a fresh unit-dependent default
  record. Remove parser handle/status module globals; retain a deliberate legacy
  workspace for unchanged no-argument callers. Legacy named atomic temporary files
  remain for MD, which reads them later; owned native parsing uses unnamed scratch.
- Extract existing default assignment block into an explicit defaults factory,
  preserving order and conversions. Reset owned coordinate views before each read
  so default snapshots contain no preceding atoms. Store defaults and final values
  separately in each parser; unit strings remain per-model.
- Native rt_model owns its parser; read_input(path,model%parser), then explicit
  initialization receives parser%values. Existing C/Python admission stays one model.
- Test actual Fortran A→B→A, omitted functional/k/pulse/mixing defaults, atom count
  and Cartesian/reduced/external/absent sources, units, simultaneous scratch owners,
  diagnostics and released unit handles. Run MPI1/2 and compare separate-process
  parsing, plus five-XC old/new native/MPI/direct/mixed regression.
- Audit new source before build and all application objects before use; update
  manifests and stage base pin; paired manuals/AI/checklist, source/wheel tests.
  Preserve historical files and legacy behavior. No Actions, no new GS collection.

## Completed verification

- Actual parser MPI1/2 tests pass for A→B→A, workspace reuse, omitted values,
  differing units/atom sources, simultaneous scratch owners, closed handles and
  preservation of preexisting named temp files. Fresh legacy logs match exactly.
- New/old native parity:5 XCs ×3 MPI configurations,15 cases; both points' returned
  current, RT numeric columns and energy columns match exactly; direct restart passes.
- Two-PBE/two-TB-mBJ normal mixture:10 independent vs6 shared ranks, continuous
 128 and63→128 protected no-replay restore/failure-retry exactly agree.
- Raw internal-MPI admission/retry passes. Source/installed617 tests each, one skip;
  sdist→wheel, inventories/native hashes and installed parser/mixed restore pass.
- Source/all application-object prohibition checks retained and passed. Selected
  records and runtime saved under the designated calc2026; GS data/native checkpoint
  payloads excluded. No Actions, no Fujitsu compilation or cloud-sync claim.

Next work is heterogeneous model registration/activation, remaining common derived
resources and Python handles. Sequential parser ownership is complete within the
selected scalar-TDDFT scope; heterogeneous live propagation remains disabled.
