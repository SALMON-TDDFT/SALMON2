# Mixed LD/SBE/TDDFT validation — 2026-11-16

Follows the LD/SBE test. A new standard PZ runtime was explicitly built from the bundled pinned sources (MPI OFF, Libxc OFF, GNU Fortran/OpenBLAS,2 build workers), source/object audited and packaged with dependency licenses/hashes. Existing binaries were not patched. The installed wheel ran the mixed sample against the relocated managed runtime.

The fixed Si starter input generated a real SCF GS:8 atoms,16³ grid,2×2×2 k sampling,54 iterations,energy -31.74653905752293 Hartree. The verified TDDFT loader snapshot was copied to persistent storage for MPI use; the preparation script explicitly uses an internal snapshot attribute for this developer fixture, not a new public export API. This small GS is a reproducibility starting point, not converged optical response.

Four outer MPI ranks: field0,LD1,SBE2,TDDFT3. LD and SBE parameters are synthetic; TDDFT uses the Si GS and a managed serial native child worker, one OpenMP/BLAS thread. dt=.02,100 steps,3 material points. Current feeds the shared Multiscale1D field at every step. A/c and charge current j match the independent serial field/raw-material schedule exactly at every point and step. Maximum |j| is [1.367256868541919e-5,2.1255286504850842e-5,1.8292756772865497e-5] atomic units. No zero-response parity claim.

Invocation (substitute your own prepared compatible paths):
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT --runtime RUNTIME --bundle GS_BUNDLE --steps 100
```
Output retains point IDs10/20/30 and material keys, A/j histories and per-rank communication timings. Timings include collective wait; they are not isolated network latency. Aggregate RSS was not measured or capped, so no measured4GiB guarantee is claimed. Parent TDDFT rank waits for its serial worker; this is not internally parallel native TDDFT.

Limits: one field rank, one point per material rank; endpoint LD uses finite-difference E, not native E/H semi-implicit LD. Production field-local many-LD topology, MPI inside TDDFT groups, whole-session restart, arbitrary material convergence and cross-platform validation remain open. No GitHub Actions were run. The preceding changes were pushed as a0ab5f4d; this follow-up evidence is separate.

Persistent evidence root: `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116` (provenance only, not a required user path).

- `mixed-runtime-build/build.json`: `7d8260d2b9eacfec61f18f5153007917790787e20bd9284e9b77bceb4d6f29d5`
- `mixed-runtime-package-v2/runtime/runtime.json`: `52efbafeeef18befdc1b4a5805b5a09ed2d617fa05770b4400e73b4cf6918e7b`
- `mixed-si-preparation.json`: `827e52c4f512f91a768f353f18377b24c743a689efcee12cf166f65039d2b114`
- `mixed-three-si/result.json`: `1285dd0d1ec5f7be2b71be854f5d7d54718d5e861f719d150a2d90e03e531801`
- `mixed-three-si/response.npz`: `64673168378a40c1c18c76c7973cf86a0a4f870e62ac79769e70e0660e9ab7ba`
- `prepare_mixed_si.py`: `73639b3e5cba7163c4fd7636b47436be07aa15faa8f98594a9221b34875f55d1`
