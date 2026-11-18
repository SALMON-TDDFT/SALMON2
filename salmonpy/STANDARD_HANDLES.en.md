# Standard TDDFT native same-model handles

[日本語](STANDARD_HANDLES.md) · Topic: `standard_tddft_handles` · Revision: `2026-11-18.1`

An explicitly built standard host-MPI runtime now exposes
`sp_tddft_handle_protocol()==1` and `sp_tddft_model_protocol()==1`. The matching
header is [native_adapters/standard_tddft_handles.h](native_adapters/standard_tddft_handles.h).
This is the experimental **native C/C++ interface**. The existing Python
`NativeTDDFTBackend.create_state` and normal MPI plans still use one point/group
and one native lifetime/process; a Python shared-model factory is subsequent work.

One process/model uses an immutable periodic, non-symmetric, unpolarized, fixed-ion
PZ/TB-mBJ/Libxc-PZ/PBE input with middlepoint transverse propagation. **dt, capacity,
decomposition and initial A are common.** This does not allow different material
inputs, symmetry groups, per-state dt/capacity, concurrent threads or other native
solver initializers in the same process. Use a fresh directory with shared,
verified GS/PP/input files and keep them unchanged throughout the attachment.
Only RT/current and energy streams are owned by this profile. Optional analysis
outputs using global names are refused, rather than silently dropped or overwritten.
The legacy single-state route retains its documented output/profile scope.

| Entry (`sp_tddft_` prefix) | Meaning |
|---|---|
| `handle_init(input,dt,capacity,A0,decomposition,&h,j)` | Read model once and create its first state; mutually exclusive with legacy init. |
| `model_new(&h,j)` | Fresh independent state from original GS, including after all earlier states close; no per-state GS directory copy. |
| `handle_step(h,Aold,Anew,j)` | Advance the selected owner using the unchanged kernels. Aold must equal its accepted endpoint. |
| `handle_info(h,&step,A,j)` | Read accepted endpoint and current. Time is step×model dt. |
| `handle_checkpoint(h,path,mode)` | 0 save,1 fresh direct load,2 stage,3 commit,4 cancel. Files are rank-local and host-authenticated. |
| `handle_transaction(h,mode)` | 0 capture,1 numeric rollback/close own streams,3 reopen,2 release. Host restores owner stream bytes between1 and3. |
| `handle_label(h,buffer,capacity)` | Output prefix; prepend the common base_directory and append `_rt.data` / `_rt_energy.data`. Buffer includes trailing null. |
| `handle_close(h)` | Release this owner, pending/protected data and initialization MPI children; siblings and host MPI stay live. |

Handle calls return **charge-current density j** in atomic units, minus SALMON
matter current. Text files retain the existing SALMON conventions; read their
column definitions. A is Cartesian A/c in atomic units; dt is atomic time. The
**legacy** `sp_tddft_init/step` ABI still returns matter J, and Python performs its
existing sign conversion. Do not negate a handle result again. x/y/z and rotating
waveforms are permitted by this full-k profile.

Host initializes MPI and calls `sp_attach(MPI_Comm_c2f(comm))`. All group ranks call
in the same order on the main thread. Tokens are opaque positive int64 values,
not user material IDs or pointers. Reuse increments generation; stale or mixed
rank tokens reject before numerical selection. Label buffer errors agree across
ranks. Close every handle before `sp_detach`, then let the host finalize MPI.
One attachment lifetime/process is the supported contract. Importing SALMONpy
neither builds nor downloads this runtime.

Wavefunctions, density, current/A/clock, XC payload, protected snapshots/output
metadata, staged restore and serializer bounds belong to each context. A parked
pointer selector does not copy any owning object. Construction temporarily selects
a unique output name; model settings are not switched at each step. Outer timers
are model-owned and stop with the last state. Idle recreation restarts them.
Initialization child communicators retire only after that owner's requests/XC
close; registry memory is bounded by peak live states, not total creations.

Direct restart requires the same binary, immutable inputs, dt/capacity and MPI
partition. Load into a fresh state; history resumes without accepted-step replay.
Numeric restore does not merge previous text logs. Host owns identity/file hashes
and per-owner output journaling. These files are not portable canonical archives.
Known numerical and transactional rejection retains its existing guarantees;
initialization/allocation/I/O/MPI/process failure does not gain a general retry
contract from token identity.

Developers can run the installed validation tool from
`salmonpy/_runtime_source/tools/test_standard_tddft_handles_mpi.py`, or the source
checkout's `tools/` entry:

```sh
mpiexec -n 6 python tools/test_standard_tddft_handles_mpi.py \
  NEW_BUILD OLD_BUILD GS NEW_OUTPUT --group-size 2 --split orbital
```

Three MPI groups compare two interleaved handle states to two independent legacy
controls. The test also grows the registry, saves at8/9 and resumes fresh states,
interleaves protected/staged candidates, closes siblings, checks stale tokens and
recreates after an idle model. This verifies ownership/interfaces, not material
convergence. See [coupling](COUPLED.en.md) for the existing Python route and
[HSE sharing](SHARED_HSE.en.md) for its separate Python shared-model API.
