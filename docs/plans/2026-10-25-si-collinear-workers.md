# Si collinear exchange workers implementation plan

> Execute inline with executing-plans and test-driven-development. Human approved existing worker extension on 2026-10-25.

**Goal:** Distribute symmetry-reduced full-q Si exchange over MPI r/rprime workers and OpenMP in both GS and learned ACE RT.

**Architecture:** Preserve 80 k-owner ranks, their SCF communicator and ACE histories. Generalize existing bounded worker tiles to scalar spin component without altering the default two-component SOC route. Source expansion retains every little-group projector and normalization. Helpers participate only in exact exchange; commands and metadata remain world-collective outside OMP. Reject unsupported geometry/backend/mode before service entry.

**Tech stack:** Fortran, MPI FUNNELED, OpenMP, serial FFTW plans, existing CMake forbidden-eigensolver guards.

## Tasks
1. Add scalar worker test comparing actual full-q action with direct distributed reference, including nontrivial little-group averaging, zero-owning helpers, OMP1/4, rprime groups1/2/4 and source NaN/Inf collective rejection. Observe failure before code changes.
2. Generalize worker/FFT component counts and metadata, retaining old SOC defaults. Connect occupied scalar source expansion to service adapter before original orbit route. Extend explicit eligibility to static collinear full-support symmetric HSE06 GS and Taylor4 RT (off/learned only); keep all other guards.
3. Build with source prohibition audit then all application objects audit. Run scalar/SOC worker regressions plus same-owner GS and learned-RT reference/helper comparison, finite/end/all child exit, current/energy/SVD and predicted/exact diagnostics.
4. Review changes and transfer independent source candidate through native authenticated Firefox OOD only. Fujitsu build and full linked object audit. Submit one minimal feature verification job (not production retries). Require all actual streams/rank exits/finite/accuracy/OMP/memory/SVD and final PJM pass.
5. After verification, submit new independent Nr16/Nr20 Nk8 GS+Pump-only RT once each, 160/320 MPI, OMP12 (40/80 nodes), 80 k owners, dt .03 au, 1.55eV 1e12W/cm2 pulse21.34538166fs plus1fs, original yyn group, no separate startup. Preserve all existing jobs. Register in same96 monitor, collect nonGS originals and common-window current/energy/HHG differences to OneDrive-direct calc2026. No auto retry/threshold relaxation; no GS collection; cloud sync unverified.
