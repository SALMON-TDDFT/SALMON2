# Standard TDDFT same-model native handles

## Delivered scope

Experimental C/C++ handle/model protocol1 for one immutable full-k periodic,
unpolarized, fixed-ion semilocal model (PZ/TB-mBJ/Libxc-PZ/PBE), middlepoint and
transverse. Common dt/capacity/decomposition/initial A. Each pointer-owned context
has independent wavefunctions/density/A/current/clock/XC, protected buffers and
output metadata, pending direct restore and serializer bound. Registry growth
copies pointers, never owning contexts; positive int64 generation tokens reject
stale or inconsistent selections before native computation. New handle boundary
returns charge j; legacy matter-J ABI and native checkpoint-v2 encoding stay.

Input read once. New states load original GS with unique construction-only labels;
no hidden numerical template/per-state GS directory copy. RT/current and energy
streams are owned, optional global-name analysis outputs refuse this new profile.
Legacy output scope stays. Existing stream files are not overwritten on handle
construction. Output label storage matches existing256-character sysname rather
than adding an80-character restriction. Model outer timers stop with last close;
idle recreation is supported. Context initialization MPI children retire at close.

Native C only in this task. Existing Python create_state and ordinary standard
MPI plans retain one native lifetime/point per group. Python shared-model factory,
normal multipoint integration, different dt/capacity or heterogeneous inputs,
other native initializers in one process, concurrent entry and portable archives
remain separate work. One attachment lifetime with explicit host cleanup.

## Verification

Local macOS26 Apple Silicon, Python3.14, Open MPI/GNU15/Libxc4.3.4, OMP/BLAS1;
64GiB/18 cores. Explicit fresh source snapshots/build outputs, mandatory source/
object prohibition audits and shared LD/SBE byte-equality checks. No Actions.

- Old runtime fails explicit missing-handle-capability test.
- Final numerical campaign:5 XC cases (PZ, Libxc-PZ/PBE, TB-mBJ auto/fixed) ×
  MPI1/orbital2/rgrid2 =15 cases. Three solver groups compare interleaved A/B
  with independent preceding-runtime singleton controls.32 steps, three nonzero
  Cartesian components including rotating xy. Initial and propagated j match
  exactly on all ranks; finite/nonzero responses required. Both original RT
  energy streams match numeric values in the old independent streams exactly.
- Six simultaneous states force registry growth. Extra siblings close without
  changing originals. Both owners hold protected snapshots simultaneously;
  releasing B permits its pending candidate while A protection remains. A trial
  step/rollback and output restoration leave B's candidate intact. Both pending
  candidates coexist; canceling A does not unblock B until its own cancellation.
- Save at8/9, fresh direct load and remaining24/23 steps match original trajectories
  with zero accepted-step replay. A's protection survives sibling creation/load/
  propagation/close. Closed/reused token rejects12, mixed valid rank tokens reject11,
  label capacity disagreement rejects13 collectively. B's pending candidate
  survives A close. All children retire at final close; idle model fresh creation
  starts original GS and matches its first step. Host MPI remains usable.
-96-character first output label passes the same full MPI2 probe against old
  controls. C and C++ header syntax checks pass.
- Legacy protected mixed regressions (PBE, Taylor HSE, PT-CN HSE, combined),32
  steps/two field ranks, match previous baselines. Ordinary PBE/two-LD/SBE/shared
  HSE/FDTD128-step direct63→128 and handled late candidate failure/retry pass;
  only65 new calls, no replay. These legacy probes used the same-kernel predecessor
  before the final label-storage-only widening; final15-case controls use the
  released runtime.
- Source/isolated installed suites603 tests each, one skip each. sdist→wheel,
  native checksums and installed MPI2/rgrid TB-mBJ-auto shared probe pass.
  Source and experimental-step preparation from installed final wheel both pass.
  Paired guides, AI catalog and checklist synchronized; read-only review found
  no remaining blocker. No physical-convergence or cross-platform certificate.

Initial missing import/link failure, timer-lifetime correction, raw unselected
checkpoint guard and collective label-capacity review are retained in developer
logs along with final successful evidence. Complete source inventories, original
upstream licensing/provenance and fixed binaries are preserved. Evidence outside
Git: developer outputs salmonpy-standard-handles-20261118.
