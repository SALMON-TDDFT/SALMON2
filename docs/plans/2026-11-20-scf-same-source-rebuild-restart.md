# Explicit same-source rebuilt-runtime SCF restart

Keep strict binary identity by default. An explicit compatible_runtime=True may
accept another runtime SHA only when the checkpoint records an equal audited
source/adapter fingerprint, native stream representation probe, build schema and
OS/architecture. The native probe hashes exact storage sizes/bytes for integer,
integer8, real8, complex8, true/false and character plus HSE/FFTW layout flags.
No numerical state or files are touched by the probe. Same-source means all
recorded solver/adapter files, not only a filename or user-declared version.

Keep input/assets/material checks and MPI admission, checksums, protected staging
and zero-replay semantics. Record source/current runtime provenance, never edit
old manifests. Older manifests/runtimes retain exact-binary restore only.

Verify two real separately built binaries from one audited source (different
SHA, same format contract), ordinary/HSE, same-layout and mixed changed ranks,
exact cut,bootstrap0,continuous comparison; rejected default/source/ABI mismatch
before mutation. Unit tests, installed-wheel/sample and JP/EN/AI/list. No CI.
Only same-platform rebuilds are certified here; source-version migration and
portable archive encoding remain separate work.

Accepted: real HSE/ordinary A/B binary hashes differ and format/source contracts
match; probes work without MPI/GS initialization or file output. Ordinary/HSE
same-layout,ordinary4→2,HSE2→8,reduced8→4,pre-PBE2→8: exact cut,bootstrap0,
continuation within unchanged1e-9/1e-10. Default/missing/source/stream/platform
refusals happen before owner creation; corrupt/late-read controls preserve source
files. Normal HSE full/reduced and ordinary results exact vs previous runtimes.
Source/installed637 tests each636pass/1skip; wheel HSE reduced8→4, ordinary
same-layout and public rebuilt-runtime8→4 sample pass. Source/wheel/docs/header/
audited-native bytes match. Source/wheel ordinary+HSE GS export explicitly
authorized; evidence in calc2026/salmonpy-scf-rebuild-record-20261120. No CI.

One harness attempt incorrectly reconstructed periodic halos for same-layout
load; it now checks the saved rank's exact work values. No solver arithmetic or
tolerances changed; failed diagnostic output is retained. Initial/final wheel
iterations are retained; the final sample fixes an obsolete same-build label
and emits actual runtime provenance. Source-version migration remains separate.
