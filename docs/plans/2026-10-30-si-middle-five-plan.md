# Si middle five-intensity Implementation Plan

> Execute inline using executing-plans; human authorized option1 and submission.

**Goal:** Submit one validated middle allocation for ten Si pump/probe cases.
**Architecture:** Reuse validated balanced binary and existing physical templates.
One GS group, strict shared dependency gate, then10 disjoint concurrent RT groups.
**Tech Stack:** Python3.6, bash, Fujitsu MPI vcoordfile, authenticated nativeFirefoxOOD.

1. Read actual templates/gates and current PJM; preserve old queued job until ready.
2. Add independent tools/si-middle-five-20261030 runner/prepare/verify/submit scripts.
   Test layout coverage/disjoint coordinates/30-31 pairs and reject overlaps,
   missing case/invalid K/sourceinput changes. Fix K8-fixed SVD count with exact K.
3. NativeOOD prepare new remote independent root, copy only static dependencies;
   validated binary SHA, source provenance, 10 physics variants, GS/shared gates,
   syntax/line lengths/static SHA/quota, all-rank gated runtime. No new executable.
4. Run local+remote fixture tests and inspect all job/launch paths and vcoordfiles.
5. Immediately verify old52968000 QUE/unstarted/no outputs then cancel once,
   record receipt and history; submit new job once with lock and raw receipt.
   No ambiguous retry. Independent PJM name/nodes/MPI/24h confirmation.
6. Save nonGS preparation receipt originals canonicalcalc2026 with all sizeSHA,
   update existing96 monitor and keep all old evidence. Report new jobID.
