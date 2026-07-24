# DG-DC Ground-State Production Adapter Design

## Scope

Complete Task 4 of the local-periodic DG-DC plan by connecting the adaptive
continuation driver to the production DC ground-state route. The route remains
default-off and must not invoke LCFO, WPW, checkpoint publication, or RT.

## Architecture

Add a GS-native adapter between the common nodal state and the existing DC
ground-state objects. The adapter applies

```text
H(lambda) = H(volume + nonlocal) + lambda * H(physical SIPG faces)
```

so only the complete Hermitian physical interface operator is continued.
Auxiliary local-box periodic wrap is never treated as a physical face.

The adapter owns no independent density-functional implementation. It reuses
the existing DC density and potential route:

1. Reconstruct each fragment's spin density from core-owned nodal orbitals and
   occupations.
2. Assemble the total-system density through `calc_rho_total_dcdft`.
3. Apply `dg_dc_gs_density_mix_rate` as a simple transactional linear density
   mix. Do not continue the pre-handoff `s_mixing` history.
4. Reuse `update_density_and_potential` for the existing DC Hartree and LDA-XC
   calculation, producing `dc%Vh_tot`, `dc%Vxc_tot`, and `dc%vloc_tot`.
5. Reuse `calc_vlocal_fragment_dcdft` to distribute the updated total local
   potential to fragment storage, then map it to the nodal operator layout.

The continuation begins immediately after Task 3 candidate materialization
and before any output consumer. It runs collectively on `dc%icomm_tot`.

## State and ownership

The adapter records the configured candidate count and the valid local
candidate count separately. Padded columns used only to make MPI array shapes
uniform are excluded from occupations, density reconstruction, normalization,
and projector diagnostics.

Rollback restores nodal orbitals, density, potential, continuation lambda, and
adapter iteration state. Mixing reset means replacing the linear-mixing input
density with the last accepted density; no hidden DC mixing history survives
the handoff or rollback.

## Diagnostics and failure behavior

Every production callback returns diagnostics measured from the same
Hamiltonian and density epoch:

- complete-H orbital residual;
- occupied-subspace/projector change;
- orthogonality defect;
- SIPG Hermiticity and physical-face balance defects;
- integrated electron count;
- eigensolver convergence and iteration count.

Fingerprint changes, nonfinite data, incomplete ownership, MPI disagreement,
or failure of the unmixed `lambda=1` fixed-point check fail collectively.
Failure does not publish a ground state and does not fall through to LCFO,
WPW, checkpoint, or RT.

## Verification

Extend the two-rank continuation fixture with a production-shaped adapter
fixture that proves:

- only the SIPG face action is lambda-scaled;
- the DC density/Hartree/LDA-XC/vlocal call order is preserved;
- padded candidates contribute neither charge nor subspace diagnostics;
- reset and rollback restore density and potential together;
- `main_dft` invokes continuation only for the default-off DG-DC route.

Run the focused continuation, SIPG, complete-H, nonlocal, density, and normal
DC route regressions. Complete specification and code-quality reviews and
resolve all Critical/Important findings before committing Task 4.
