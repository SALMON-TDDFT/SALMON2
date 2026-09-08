# Wannier-first DG flux diagnostic

## Goal

Check the static spectrum of the DG Hamiltonian after the basis has been
converted to Wannier functions.  The purpose is to avoid tuning the flux strength
by hand and instead make the boundary derivative well defined in the chosen
Wannier basis.

## Principle

The no-flux DC/LCFO eigenstates are transformed to global Wannier functions.
Each fragment owns the Wannier functions whose centers belong to that fragment.
Before RT propagation, SALMON reconstructs those Wannier functions on the
fragment buffer box and exports their values and outward normal derivatives on
the six core faces.

The DG surface Hamiltonian is then assembled from these Wannier traces:

```text
u_I(face),  du_I/dn_I(face)
```

This keeps the derivative definition tied to the actual real-space Wannier
function, rather than to a core-truncated LCFO fragment basis.

## Implementation Steps

1. Add a per-fragment trace file:
   `buffer_periodic_wannier_trace.bin`

   The file stores:
   - fragment id
   - core size, buffer size, box size
   - grid spacing and volume element
   - number of kept Wannier functions
   - for each of the six faces, the face values and outward normal derivatives

2. Write the trace file from the global Wannier BPW path.

   The writer reconstructs

   ```text
   w_phi(r, iw) = sum_io phi_lcfo_buffer(r, io) * wcoef(io, iw)
   ```

   and differentiates on the buffer box using the same SALMON stencil.

3. Add a small diagnostic projector.

   The projector reads the Wannier BPW files and the new trace files, assembles
   the static DG surface terms in Wannier space, and reports the HOMO-LUMO gap
   for:
   - local H only
   - surface self only
   - surface cross only
   - local H plus full DG surface

4. Run the 2x2x2 diamond case with `OMP_NUM_THREADS=2` and `mpirun -np 8`.

   The first acceptance criterion is that the Wannier-DG Hamiltonian gap remains
   physically close to the no-flux LCFO gap without manual scaling of flux terms.

## Non-goals

- No lambda-cut tuning to match the gap.
- No production RT wiring until the static Wannier-DG spectrum is credible.
- No all-to-all position matrix or global density matrix in the RT path.
