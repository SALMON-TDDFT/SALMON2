# Distributed LCFO exchange construction implementation plan

**Goal:** Remove repeated global coefficient rotations and exchange/ACE construction from native LCFO RT, preserving the physical operator and fragment x orbital MPI.

**Architecture:** Store C, transported frames, W and ACE factors by core rows. Exchange only the rows required by each fragment's basis; reverse that exchange to sum fragment operator contributions back to their row owners. Keep the fallback exchange as fragment blocks, including both midpoint endpoints. Assemble orbital-space metrics from local contributions, solve small problems on the root, and use ScaLAPACK for >=128 occupied states when built with it. Orbital-space rotations/metrics remain dense; this is distributed computation, not linear-scaling sparse orbital algebra. Initial MLWF seeding and the existing diagnostic dump may gather coefficients to one root once.

**Tech Stack:** Fortran, MPI Alltoallv, LAPACK/ScaLAPACK, existing SALMON BLAS and native RT.

The user approved this design and execution. Reuse the clean dedicated dc-hse-mlwf-ace checkout. Preserve the frozen b1773c78 diamond binary and all queued comparisons. Run no numerical regression alongside those comparisons; preparation and compilation may proceed.

1. Add a standalone MPI regression `testsuites/unit_lcfo_rt/test_distributed_build.py` and `distributed_build_probe.f90`: unequal row counts, noncontiguous halo requests, adjoint accumulation, polar transport, ACE action versus dense reference, singular/zero rejection, and a >=128-state distributed decomposition. Compile against proposed modules and confirm the missing-interface failure before implementation.
2. Add `src/xc/lcfo_dist_rows.f90`: static request plan, coefficient-row exchange, adjoint contribution sum, and root-only initialization gather. Test the actual sum of overlapping Hermitian fragment blocks, not merely a mirrored pack/unpack.
3. Add `src/xc/lcfo_dist_dense.f90`: local overlap and ACE metrics, collective decisions, root LAPACK and ScaLAPACK SVD/eigensolver. Preserve strict metric rejection and zero exchange. Register both modules in `src/xc/CMakeLists.txt`.
4. Convert `lcfo_rt_wannier.f90` to local coefficient/frame rows. Keep initial file format and deterministic seeding. Use distributed transport and halo rows for reconstruction. Update standalone localization fixtures to the local-row API.
5. Convert `hse_lcfo_rt.f90`: local coefficient collection, collective cache decision, fragment-block Hx, W by adjoint accumulation, local ACE factors, and distributed fallback application. Initialize matching halo plans for every orbital group. Preserve impulse/laser, predictor rollback and midpoint averaging. Add phase timings and local storage diagnostics.
6. After diamond completes, run the MPI probe with/without ScaLAPACK, existing transport/sphere/action tests, and native MPI2/MPI4 regressions. Cover ACE-invalid fallback and unequal orbital columns. Fix failures from evidence.
7. Build the ScaLAPACK-enabled executable separately. Compare frozen versus new Si64/Si128 RT on the same GS, R9/ACE4, MPI8/16 OMP1 BLAS1, one numerical job at a time. Check initial gauge, current, density, energy, timings and sampled memory.
8. Request review, document measured improvements and remaining dense orbital-space costs, commit and push the existing branch. Preserve all diamond results in their original study.
