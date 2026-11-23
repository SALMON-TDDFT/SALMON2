# Ordinary XC on the fast source profile — 2026-11-19

Local Apple Silicon/GNU15.2/OpenMPI/Libxc/OpenBLAS CPU validation. Existing synthetic
8³,2-band,8-k fixtures are copied only to temporary test directories. No new
SCF/physical GS dataset is computed or collected. This is implementation parity,
not real-material convergence or platform portability certification.

## Shared implementation

The fast prepared tree now consumes the standard embedded-TDDFT/endpoint-worker/
MPI-bridge files directly. SALMONPY_FAST_CORE selects a distinct internal module
name and the fast retained-GS seed provider. Fast numerical owners already contain
stencil and XC scratch, so no legacy global-workspace move is performed there.
The standard build keeps its original ownership path. Fast initialization accepts
memory seeds by name and fills the imposed initial A for the embedded endpoint
protocol; ordinary HSE initialization is unchanged when that worker is inactive.

HSE-disabled host builds use the existing embedded schema and declare
source_profile=fast-pr, ordinary_rt_abi=1, and import-only gs_abi/1. NativeGSBackend
keeps unknown imported energy/residual/convergence as None and still refuses SCF/
forces/eigenvalues/export. The existing GS→RT adapter selects the ordinary registry
for this profile; fast-HSE memory creation retains its separate declaration/path.
GS interleaving uses the ordinary park/resume hooks. No extra Python module,
sample file or unittest method is added.

## Tests and failures retained

The old fast runtime lacks the ordinary model/memory-GS ABI (red evidence).
The first file-based heterogeneous test failed because the fast reader captured
its parser snapshot before API dt/capacity/MPI overrides. It now captures after
worker_configure, matching the standard reader. No numerical update is changed.

The existing memory-GS tool adds explicit control-build and interleave-fixture
options and a k-only split. With a different control binary, comparisons use
rtol3e-13/atol3e-15 and record actual max differences; default same-build controls
retain strict equality. Exact direct-resume and sibling/fresh-state assertions
are unchanged. The final controls compare16 three-component current steps with
an independently pinned old standard binary and resume8→16 without replay:

| XC | Internal MPI2 split | Max current difference (au) |
|---|---|---:|
| PZ | k | 0 |
| TB-mBJ | k | 0 |
| Libxc-PZ | orbital | 2.710505431213761e-20 |
| Libxc-PBE | spatial | 5.421010862427522e-20 |

External k tables, missing-asset rejection, input metadata and different-XC GS
interleaving also pass. The heterogeneous-registry tool runs MPI6, TB-mBJ and
Libxc-PBE with different dt/capacity/orbital-spatial partitions. Independent and
restart trajectories agree exactly, including56 steps, slot reuse, model-close
lifetime and rejected registration without disturbing accepted owners.

Fast-HSE MPI6 registry controls retain exact currents/energy/exchange history/
direct restart against the prior fixed runtime. HSE memory initialization and
standard MPI memory initialization are also checked. All new source and complete
application-object trees pass the mandatory eigensolver audits, including external
fast host adapters. Ordinary HSE-disabled host linking has no FFTW dependency.

A first FDTD/LD/SBE/ordinary-TDDFT run reached128 steps, but its independent control
was incorrectly given an MPI executable inside a host-MPI parent. Inherited MPI
runtime environment caused shared-memory initialization failure. That local test
child was terminated; no production job or fixed binary was changed. A separate
serial control build exposed a pre-existing missing communication_dummy GS-cleanup
interface. Minimal serial no-resource helpers now allow that storage to compile;
there is still no serial host-MPI attachment. Failures/logs are retained.

## Remaining scope

Ordinary RT is connected in the HSE-disabled profile. Computed GS lifecycle/export,
combined ordinary/HSE registry ownership within one GS session, and final removal
of the old maintained source tree remain. Different material MPI groups remain
available for normal mixed calculations. No implicit compile/download, changed
physical convergence criteria, GPU claim or new runtime Release is introduced.

Diagnostic/distribution records are local under
calc2026/salmonpy-ordinary-fast-rt-20261119-final (the base directory preserves
pre-admission records). GS, wavefunction and direct-checkpoint
payloads are excluded from that collected record set; cloud synchronization and
other machines/compilers are not verified. GitHub Actions was not triggered.


## Final distribution and coupling checks

- Source and installed wheel suites:628 methods,627 pass and1 external-runtime skip.
- sdist→wheel and installed fast-source preparation reproduce all261 src hashes.
  Installed retained_gs example computes two independent Libxc-PBE currents after
  GS close while retaining a different TB-mBJ GS owner. Installed spatial-MPI2
  Libxc-PBE memory/file/direct-resume control passes.
- HSE-disabled serial build now excludes MPI-only embedded_gs and supplies only
  no-resource cleanup helpers for legacy GS storage. It links without MPI/FFTW;
  all application objects pass the mandatory audit. No SCF was executed.
- FDTD/LD/SBE/PZ128 compares to independent serial TDDFT with maximum current
  difference2.574980159653073e-19au (SBE exact). Direct save63→resume65 reproduces
  continuous128 A, E and charge current exactly, with no step replay.
- Different-cell/XC/dt/capacity/MPI-partition model ownership also passes the
  MPI6 registry test. The changed-cell fixture intentionally reuses synthetic
  wavefunctions; it is not a physical ground state of that cell.
- Guides/AI catalog and existing sample updated; no extra Python module/sample/
  unittest method. Computed GS/export and final tree merge remain the next work.


The final ordinary-only frontend rejects HSE06 before native GS initialization;
ordinary RT creation also rejects XC names outside the implemented endpoint-worker
contract while leaving GS arrays available. The existing unittest method and
MPI memory-GS tool check that unsupported hybrid admission preserves GS resources
and accepted current, then verify future/restart trajectories. This is refusal of
an unavailable compiled engine, not a new convergence/resource cutoff.
