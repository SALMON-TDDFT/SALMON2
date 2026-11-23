# Group5: transient output units and restored conversion views

Reproduce stale conversion factors after a model unit A→B→A switch, and file
unit1/77/100/101 collisions in synthetic DOS/PDOS/band/split+unsplit AVS calls.
Keep conversion literals/formulas and output formats unchanged. Use NEWUNIT for
short-lived files; AVS closes its own unit rather than literal1. Retain its fp
argument for existing callers while acquiring a temporary unit internally.

Move derived conversion declarations/initialization unchanged into a small logical
module alongside salmon_global in the existing file. inputoutput re-exports the
same legacy names; model_inputs activation rebuilds the conversion view from its
owned selectors. Validate selectors before changing any model globals. No new
serialized fields or Python numerical work. Keep checkpoint codec bytes intact.
This is sequential activation, not a threaded/reentrant I/O engine. Existing
Python work-directory scopes and owner input names remain the output routing.

Test old failures and new synthetic MPI2 fixtures, compare output bytes with the
unchanged audited previous runtime, then ordinary/full/reduced HSE in-RAM normal
and protected-restore parity. Test parser/omitted defaults and input ownership;
source/wheel tests, docs and packaged fixture parity. Do not collect GS datasets
or dispatch CI. Broader process-global I/O/CWD/threading/MPI remain separate.

## Accepted checks

Old audited ordinary runtime fails units(code12) and DOS/PDOS/band/AVS ownership
(code21) synthetic controls. New ordinary/HSE and installed ordinary MPI2 pass,
with previous-output bytes identical and host sentinels unchanged. Parser-owned
A→B→A/defaults/tempfile/legacy logs pass; input record ownership passes.
Ordinary MPI2 and full/reduced-HSE MPI4 (three groups) and wheel ordinary/reduced
HSE: normal/protected restore/retry arrays bit exact to previous runtimes, restore
SCF work0, no GS datasets collected. Conversion routine body is byte-identical
and numerical stream writer files unchanged. Source/wheel637 tests each636pass,
1skip. Source/object audits and package/guide sync pass. Evidence:
calc2026/salmonpy-gs-io-record-20261120. No CI; cloud sync not established.
