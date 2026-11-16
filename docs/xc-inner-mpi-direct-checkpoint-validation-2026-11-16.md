# Non-PZ inner-MPI direct checkpoint — 2026-11-16

Code:ea22fe15, installed macOS arm64 cp314 wheel from xc-direct-final-dist. No solver or Python numerical changes were needed for this campaign. Host library SHA256:98f77a50cbfe66f18df7a96c20078a0216b5bb975823e9ba6a2a97ce7fb25627; embedded TDDFT source SHA256:f41d4bad46ae752542478f3c33e32a3a20ad3a06c044cbd430b4d7ece18c4a7f. Libxc4.3.4 as documented in xc-direct-checkpoint-validation-2026-11-16.md.

## Actual execution

Apple Silicon/macOS, Python3.14, OMP/BLAS1. Four outer MPI ranks: field rank0 holding two local LD materials, SBE rank1, TDDFT ranks2/3 in a separate communicator. Actual native log reports nproc_k=1,nproc_ob=2,nproc_rgrid=1,1,1, rather than inferring parallelism from a command argument.

For each matching-XC converged Si SCF fixture (16³ grid,2³k,8 atoms), execute continuous128, save63, then a NEW job/output loading internal arrays and continuing through128. dt.003, same40³ field and mixed layer geometry. Sample options: `--tddft-ranks 2 --rank-budget 4 --layers LAYOUT_JSON`; save adds `--direct-save-at 63`, resume adds `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. Use the same host library, serial reference executable and GS/RT bundle across the three jobs. See COUPLED guides for the full command entry.

| XC | Restored step | Actual scheduled native continuation steps | Max complete-history difference A/c,E,j | Max TDDFT2 vs independent serial j difference |
|---|---|---|---|---|
| TB-mBJ automatic cval | 63 | 65 | 0,0,0 | 6.911788849595091e-19 |
| TB-mBJ fixed cval1.2 | 63 | 65 | 0,0,0 | 7.436949276892757e-19 |
| Libxc-PBE | 63 | 65 | 0,0,0 | 6.776263578034403e-19 |

The comparison explicitly uses np.testing.assert_array_equal against the SAME two-rank continuous trajectory. The serial reference is a separate validation, with small decomposition roundoff; it is not used to reconstruct the resumed trajectory. SBE serial replay is exact. Existing native/electronic clock and previous-current feedback checks pass.

All checkpoints contain distinct `electronic-2.bin` and `electronic-3.bin`, plus field state, per-rank orchestration/SBE state and saved observable frames. Every file hash was rechecked against complete.json and ranks=4. Payload sizes82,191,003/82,191,084/82,191,025 bytes respectively (excluding manifest). These fixture sizes are not aggregate RSS measurements or memory caps. No past propagation is performed in the resumed native loop; its range is restored_step..steps and result.json reports63/65. Static setup and independent serial reference validation still occur.

Persistent evidence beneath `~/.codex/outputs/salmonpy-state-inputs-20261116`: `xc-mpi2-{tbmbj-auto,tbmbj-fixed,libxc-pbe}-{continuous,save,resume}`, their `.log` files and `xc-mpi2-direct-comparison.json`. These are developer evidence paths, not user prerequisites. Focused direct-checkpoint/manual/source-manifest tests22 pass, guide synchronization and whitespace checks pass. Sdist-to-wheel build in `xc-inner-mpi-dist` passes; wheel Japanese/English guides and AI JSON match source bytes. No GitHub Actions used.

## Scope and remaining work

Same binary, initial data, environment and topology. Does not establish changed-topology restart, HSE resources, other inner decompositions, generic custom-backend persistence, atomic cross-group rollback, production scaling or optical convergence. Field2+SBE1+TDDFT2 requires five ranks and remains untested within the four-core budget. No Linux/Windows/Intel-Mac native validation is claimed. Historical native text output file positions are not resumed. Molecular work remains paused.

The active checklist's old delivery/input-validation/layer-setup/checkpoint bullets were reconciled against completed commits. Broader GS/TDDFT ownership, recoverable native errors, memory estimates/caps, public generic checkpoint contracts and environment/topology coverage remain open.
