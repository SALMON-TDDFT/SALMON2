# Molecular feedback public entry implementation plan

Goal: move the verified one-molecule MPI coupling schedule out of the development
fixture into salmonpy.dipole_source.MolecularFeedback. The caller owns MPI, both
fields, the isolated state, mapping and output paths. No implicit preparation,
SCF, spawning, download or field lifetime ownership. Existing isolated state
creation may explicitly start its normal worker.

Architecture: construct collectively with comm, local field (or None), a state
only on total_root, DipoleMap, explicit total_root/reference_root and explicit
policy=environment_scattered_only. Initially one molecule/point, midpoint only.
start() performs the declared zero-current hold. step() preserves the existing
predict/rollback/current/validate order and returns an immutable root-only frame.
Optional existing field energy observers are read before backup release. No
trajectory/control calculation or implicit file output. Existing field-group
backups and quantum replay are reused; no extra world agreement on successful
steps. Failure-only recovery votes report whether restoration succeeded.

1. Add a failing orchestration test to the existing isolated tests. Test start,
   rejection, preserved clocks/dipole, retry and caller-owned lifetimes. Synthetic
   responses are protocol checks, not physical evidence.
2. Implement in the existing module and public inventory. Numerical subtraction,
   average and current remain in the existing Fortran kernel. No new solver or
   native context selector.
3. Connect the actual MPI fixture to the public entry. Remove its duplicated
   coupling transaction loop; retain separate validation/replay/background work.
4. Reuse the authorized C2H2 GS, audited worker and field library. Compare the
   source continuous and installed-wheel rejection/retry trajectories and energy
   observations with the previous six-rank baseline. Source/wheel suites.
5. Add a complete caller-owned MPI example, paired human guides and AI JSON.
   Record limitations: one molecule, lab-axis dipole, own-vacuum subtraction
   removes vacuum radiation reaction, electronic restore replays history,
   no arbitrary OS/MPI failure recovery or physical convergence claim.
