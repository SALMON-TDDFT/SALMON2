# Standalone Si512 owner exchange modules

This directory preserves the adopted standalone Si512 implementation. It is not
connected to the SALMON runtime or included in its CMake source list. Integration
into SALMON's native exchange/restart interfaces remains separate work.

The source snapshot includes ordered return point bundles, source/return MPI
datatype reuse, cached return descriptors and task tables, nonparticipant task
skipping, and double-buffered return communication/computation overlap.
The source snapshot main module SHA256 is
`bde5b2a0f907dc5d20ea1f452539b30dd299db9aec9a05cce63c8e0d2c652b7a`.

It excludes the unassessed preparation-OMP, source-plan-cache, and input-pipeline
candidates submitted as Fugaku jobs 53989239, 53989489, and 53990460. No executable,
GS/restart data, job-management scripts, Python runners or benchmark output is
included here. Normal SALMON build and runtime behavior are unchanged by this
standalone snapshot.

The snapshot removes trailing whitespace from the adopted source; its code and
behavior are unchanged.

## Validation

Fugaku job 53986991 (2026-11-12, 1024 MPI, one rank/node, 48 OMP) compared eight
sampled tiles against the previously adopted return-overlap version. Both
returned status zero; all compared values were bitwise identical, with 4,198,400
evaluated pairs and equal source/return payloads. Eight-tile wall time was
3.742255 versus 3.253934 seconds. This single ordered comparison includes
communication variation and does not establish full EXX or RT speedup.

Local boundschecked MPI comparisons covered multiple return/prefetch windows,
empty orbital ownership, one rank, prefetch disabled, large blocks and mixed OMP
team sizes. All compared values matched bitwise. The standalone module has not
been validated as a full Si512 RT trajectory or a native SOC integration.

## Compile independently

From an out-of-source directory, with `SOURCE` set to this directory and `FFTW_INCLUDE`
set to the FFTW include directory:

```sh
mpifort -O2 -fopenmp -I"$FFTW_INCLUDE" -c \
  "$SOURCE/wannier_pair_owner.f90" \
  "$SOURCE/wannier_patch_owner.f90" \
  "$SOURCE/wannier_neighbor_owner.f90" \
  "$SOURCE/orbital_owner_layout.f90" \
  "$SOURCE/spatial_patch_adapter.f90"
```

Consumers need MPI and FFTW3 at link time. Preserve MPI-thread-funneled execution
and complete pending requests before reusing or releasing communication buffers.
