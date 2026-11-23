# Native detached restore candidates — 2026-11-17

## Implementation and scope

Rebuilt standard E/H and TDDFT libraries and the fast-PR HSE host-MPI adapter
expose checkpoint_protocol=1. The existing checkpoint function accepts mode2 to
read a detached numerical candidate, mode3 to publish prepared memory, and mode4
to discard. Modes0/1 and archive encodings remain compatible. Standard adapters
agree modes and pending readiness within their solver group before branching.
HSE candidates bind their originating context; another context cannot consume
them. The Python internal helper also rejects concurrent candidates for one owner.

HSE numeric exchange/predictor histories are read and validated once into detached
storage. Publication never reopens the checkpoint. Staging requires its immutable
kernel already initialized and does not refresh exchange into the live context.
MPI communicators, Libxc handles, output units and FFT resources are not cloned.

This is an experimental backend-owned numerical candidate protocol, not general
native reentrancy. Hosts pause evolution, prepare Python metadata, stage through
the common checkpoint entry, agree success, explicitly publish native candidates,
then install prepared Python clocks/A/current. The common holder switch alone
does not perform native publication. Old aliases cannot continue an old native
trajectory after successful publication. The ordinary mixed-native example still
uses atomic=False; it has not silently acquired this contract.

## Verification

`tools/test_native_restore_candidates_mpi.py` uses four host ranks and real
FDTD→A→standard Libxc-PBE / fast-PR HSE→j→FDTD feedback. Initial feedback is zero,
as in the existing coordinator. At step8 every owner saves its state. At step10
a later HSE shard is truncated. Candidate staging fails collectively, all
candidates are discarded, and every original native shard saved before/after
the attempted restore is byte-identical. Originals continue to step24.

Retry from step8 succeeds without replay; the following16 field/current steps
match the continuous trajectory exactly. Tested Taylor layouts (field, standard
TDDFT, fast-PR HSE rank counts): (1,1,2), (1,2,1), (2,1,1). These also reject
mixed publish/discard modes within a two-rank group without hanging. PT-CN is
tested with (1,1,2), including the installed-wheel path. Taylor traces from all
four ranks match the previous unchanged native binaries exactly (maximum error0),
as do installed-wheel PT-CN traces. The legacy two-rank HSE immediate-load
save/resume test also passes (split7→12, dt=.001).

Do not confuse failed-restore byte preservation with post-success archive byte
identity: native restoration deliberately marks orbital overlap caches dirty.
Final electronic archives differ in the unused/cache validity flags (2 or4 bytes
per electronic rank in this fixture); sampled A/E/j agree exactly and the field
archive agrees in full. Full before/after and final streams are retained for
inspection. These short fixtures verify ownership/update order/restart behavior,
not converged optical properties or cross-platform physical accuracy.

Source and installed-wheel Python suites each run500 tests (one optional skip).
New tests cover stage/discard, failure cleanup, no checkpoint I/O during publish,
old-library rejection and preventing a second candidate from discarding the first.
Source-distribution-to-wheel build, installation, paired guide synchronization
and source inventory verification are included in the local campaign.

Artifacts/logs/build receipts:
`/Users/otobetoshihito/.codex/outputs/salmonpy-native-candidates-20261117`.
Final libraries are standard-build-v4 and hse-build-v4. New output libraries
reuse unchanged core objects read-only with recorded compile/link commands and
source/adapter hashes; existing fixed binaries are never modified. Builds use2
workers and tests at most4 MPI ranks, OpenMP/BLAS1. No remote Actions run. RSS is
not a hard enforced resource cap.

## Remaining limits

Publication still includes allocating assignments. Process death, OOM, failed
communicators, Fortran STOP and a caller evolving states during restoration are
outside handled-exception recovery. There is one pending candidate per native
lifetime, same binary/input/MPI topology only. Stable public native ABI, automatic
whole-session integration into the ordinary mixed-native coordinator, topology
redistribution and unrestricted simultaneous native contexts remain separate work.
