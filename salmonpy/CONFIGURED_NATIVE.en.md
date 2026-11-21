# Dedicated configured LD and TDDFT bindings

Use separate factories to retain each calculation implementation's numerical and runtime contract:

```python
from salmonpy.materials import (
    bind_configured_ld_endpoints, bind_configured_tddft, response_identity,
)

# Choose exactly one for the corresponding ConfiguredMaterial:
ld_binding = bind_configured_ld_endpoints(ld_material, key='ld-cache')
tddft_binding = bind_configured_tddft(tddft_initial, key='gs-cache', runtime=runtime_directory)
```

LD explicitly converts endpoints to `E_new = -(A_new-A_old)/dt` using the existing
LDResponseBackend. This is not the solver semi-implicit E/H coupling and does not
replace the distributed FDTD-local LD path. Initial field/history are zero.

TDDFT requires a ConfiguredMaterial holding TDDFTInitialState (from the configured
GS bridge or tddft_from_catalog). The explicit managed runtime must support serial
TDDFT and its executable checksum must match the snapshot. The runtime is checked
at binding construction, prepare and by the existing calculation implementation at worker creation.
Binding and prepare borrow the initial snapshot: no extra copy, SCF, process start,
network or compilation. Actual state creation uses existing independent child
workers and preserves solver updates, A/c input and charge current j=-J_SALMON.
The prepared model must stay open until all sessions close. These borrowed providers
do not close it; close the ConfiguredMaterial explicitly afterward.

Both bindings contain immutable material/optional initial-state records. Their
provider identities also include numerical data: LD provider parameters, or TDDFT
snapshot manifest/input, executable hash, max_steps, threads and timeout. Absolute
runtime paths and object addresses are not identity. Point/run metadata still has
the ordinary material/initial digests; calculation implementation settings remain in the provider
identity and should be recorded separately if needed for full provenance.

```python
bindings = {10: tddft_binding, 90: tddft_binding}
point_map = {10: 'positive', 90: 'negative'}
record = response_identity(bindings, point_map, run_id='experiment-1')
session = ResponseSession(context, bindings, dt=0.02,
                          point_record=record, point_map=point_map)
```

Follow [MPI_IDENTITY.en.md](MPI_IDENTITY.en.md) for collective lifetime and output.
Each material group has one MPI rank. This does not enable group-parallel solver
TDDFT or in-process multiple TDDFT states. The factory does not certify resource
budgets, physics convergence or restart portability. LD direct E and SBE/custom
endpoint-A APIs remain unchanged.

Example (3 ranks; waveform NPZ has scalar positive dt and finite A with shape
(steps+1,3), A[0]=0, atomic units):

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_configured_native \
  inputs/materials.json metal-A waveform.npz /new/ld-output --mode ld-endpoints
mpiexec -n 3 python -m salmonpy.examples.mpi_configured_native \
  inputs/materials.json Si-ground waveform.npz /new/td-output \
  --mode tddft --runtime /local/runtime
```

TDDFT uses an initial_state_id, LD a material_id. Use an explicitly prepared local
catalog/GS, the compatible runtime and solver LD wheel. The sample drives two
points with opposite waveforms, saves IDs/current, and compares against independent
direct states after MPI sessions close. One thread per TDDFT worker; at most two
simultaneous workers. No hard memory limit. The Si starter remains only a short
connection regression, not a converged optical calculation.
