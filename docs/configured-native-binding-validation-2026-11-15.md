# Configured LD/TDDFT binding validation (2026-11-15)

M5f dedicated provider bridge. bind_configured_ld_endpoints explicitly reuses
LDMaterialProvider/LDResponseBackend and E_new=-(A_new-A_old)/dt. It is not the
native semi-implicit E/H route. bind_configured_tddft borrows an existing configured
TDDFTInitialState, verifies managed runtime/executable identity and reuses
TDDFTResponseBackend. Binding/prepare starts no worker and makes no snapshot copy.
The caller retains ownership and closes the model after sessions.

Provider identities include numerical parameters or the initial manifest/input,
executable hash and worker settings. Immutable material records flow through
response_identity and existing MPI identity checks/output. No change to native
updates, current sign, worker recovery or MPI group capabilities.

## Tests

- Three focused tests: LD no early state allocation and exact direct E parity;
  managed TDDFT runtime mismatch rejection and borrowed initial identity; wrong
  model types rejected.
- Source suite: 381 tests, one filesystem-condition skip.
- Fresh installed wheel suite: 381 tests, two skips (filesystem and optional plot).
- Built from sdist to wheel with max two jobs; installed imports isolated from source.
- Installed native LD: three MPI ranks, two independent points driven with opposite
  A waveforms; all eight steps exactly match independent direct LD E calculations.
- Installed real PZ Si TDDFT: fresh GS from packaged starter, relocated/pinned
  GS/1 selected by ID, three MPI ranks, two independent material workers. Both
  opposite-drive point histories match independent direct TDDFT runs exactly for
  20 steps. Positive-drive history also matches packaged-starter reference path
  exactly (maximum absolute current difference 0).
- Saved material/initial/point/run metadata read back correctly.
- Review caught rank-local waveform-length drift in the example. The sample now
  agrees dt/shape/content digest before constructing the session; a bounded real
  three-rank differing-length regression rejects before stepping on every rank.

## Resources and limitations

macOS Apple Silicon, Python 3.14, local MPI/mpi4py. Maximum three MPI ranks and two
simultaneous serial TDDFT child workers, each one thread; root direct reference is
run after MPI material sessions close. No peak memory measurement/hard limit.
Local MPI socket access needed sandbox escalation. Each MPI launch bounded to
45 seconds. No full solver build, GitHub Actions, commit or push.

Existing managed executable SHA256:
`29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed`.
Si input is 16³, k mesh 2³, 20 short steps. This is connection/update-order evidence,
not physical convergence, HSE/other-XC coverage or cross-platform accuracy.
One rank per material group remains required; this does not expose collective
native TDDFT or multi-instance native reentrancy. Restart identity remains open.

Local arrays/report: `/tmp/salmonpy-configured-native-mpi-20261115`.
Fresh starter preparation/report: `/tmp/salmonpy-configured-native-inputs-20261115`.
Generated restart payloads removed after verification; regenerate a new GS before
reusing its saved bundle. Reports, current arrays and manifests retained.
User reproduction: examples/mpi_configured_native.py with explicit --mode as in
CONFIGURED_NATIVE.md/en. Waveform failure worker: tools/test_configured_waveform_mpi.py.
