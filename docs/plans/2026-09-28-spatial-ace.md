# Spatial ACE implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement and certify the spatially distributed ACE algebra as the first stage of the approved native mesh parallelization.

**Architecture:** Store only local grid rows of factors and targets. Reduce the small metric before factorization and overlaps before action, via an optional reduction callback preserving serial callers. MPI peers share orbital/k dimensions, dv and call order; grid row counts may differ or be zero. Collective validation prevents a bad local value from stranding peers. No LCFO projection. Native RT admission remains guarded until distributed MLWF refresh and exchange generation are connected.

**Tech Stack:** Fortran complex BLAS/LAPACK, MPI executable regression, Python unittest.

## Task 1: Distributed algebra oracle
Create testsuites/651_hybrid_exchange/ace/test_spatial.py and spatial_driver.f90. Build with mpifort and OpenBLAS. Use a deterministic complex negative operator, independent dense action, uneven 1/2/4 rank row partitions including a zero-row rank, and multiple k/target columns. Test source interpolation, arbitrary targets, midpoint average, zero exchange and invalid local input. Observe missing callback API fail first.

## Task 2: Implementation
Extend src/xc/hse_ace.f90 build/apply with optional sum_grid complex matrix callback. Sum local metrics and overlaps; handle zero rows and globally zero exchange; collectively reject nonfinite data. Preserve serial interface and average representation. Run the new oracle and existing ACE regression.

## Task 3: Verification
Build HSE ON/OFF, independent read-only review, document exact memory/reduction behavior and remaining RT gate. Commit this usable algebra stage without claiming distributed native RT is enabled. Native integration requires the subsequent distributed MLWF/exchange stage rather than a root-gather fallback.
