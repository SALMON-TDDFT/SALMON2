# solver HSE handles

[日本語](NATIVE_HANDLES.md)

An explicitly built MPI managed by the calling program HSE runtime exports
`sp_fast_hse_handle_protocol()==1` and the declarations in
`native_adapters/fast_hse_handles.h`. This is a bounded C/C++ entry for multiple
independent electronic states of **one immutable model** in one attached MPI
group and working directory. Python legacy `create_state` retains its one-solver-lifetime limit;
[open_model](SHARED_HSE.en.md) now shares the model through independent calculation states. The new entry does not change the ordinary mixed
sample's rank layout or automatically use multiple contexts in one process.

The calling program initializes MPI, calls `sp_attach(MPI_Comm_c2f(comm))`, and keeps that
communicator alive. Every operation except the protocol query is collective,
serial and on the main thread, in the same order on all ranks. Calls never
initialize/finalize MPI, change the working directory or compile a runtime.

| Operation | Meaning |
|---|---|
| `handle_init(path,dt,capacity,symmetry,&h,j)` | First and only process-model initialization; returns its first state. Use this or the legacy singleton initializer, exclusively. |
| `model_new(dt,capacity,&h,j)` | With model protocol1, create after zero live states in the retained initialized model; no hidden template. |
| `handle_new(model,dt,capacity,&h,j)` | Fresh independent state from the original GS of a **live** model token; does not copy its current excitation or propagation. |
| `handle_step(h,Aold,Anew,j,&energy)` | Existing HSE endpoint update on the selected state. |
| `handle_info(h,&step,A,j)`, `handle_history(h,values)` | Accepted endpoint and eight rank-local history diagnostics. |
| `handle_transaction(h,mode)` | Modes0 capture,1 numeric rollback/close logs,2 reopen/flush/release,3 reopen. Failed reconnect retains ownership metadata. calling program restores saved file bytes before reopening. |
| `handle_checkpoint(h,path,mode)` | Modes0 save,1 fresh direct load,2 stage,3 commit,4 cancel; rank-local paths may differ. Existing identity/partition/authentication requirements apply. |
| `handle_label(h,buffer,capacity)` | Null-terminated output basename for calling program journaling; append `_rt.data` / `_rt_energy.data`. |
| `handle_close(h)` | Release this numerical/FFT/XC owner; does not detach MPI managed by the calling program or invalidate other states. |

All names above have the `sp_fast_hse_` prefix. Tokens are positive opaque int64
values, never addresses or user material IDs. A slot generation prevents an old
closed token from selecting a new state reusing the slot. The registry grows as
needed; allocation errors are returned collectively. Status0 means success,
11 inconsistent rank arguments,12 invalid/stale token or attachment,13 too-small
label buffer,16 registry allocation failure; other statuses keep existing HSE
meanings. Invalid/mixed initial dt/capacity is rejected before model lifetime
consumption, so corrected initialization can be retried.

`dt` is atomic time, A is Cartesian **A/c in atomic units**, j is **charge-current
density** in atomic units (minus SALMON matter current), and energy is total
energy in Hartree. Aold must exactly equal the accepted endpoint: query info
after initialization or restore. x/y/z and rotating/circular input are supported;
retained symmetry must support the actual field. The HSE implementation is the only
HSE engine. Same-model symmetry contexts require shared_symmetry_protocol1; different symmetry groups still require further ownership work.

Use [examples/native_fast_hse_handles.c](examples/native_fast_hse_handles.c) with
the same MPI/toolchain as the explicitly built library. Compile explicitly:

```sh
mpicc -std=c11 -I/path/to/salmonpy/native_adapters \
  /path/to/salmonpy/examples/native_fast_hse_handles.c \
  /path/to/libsalmonpy_fast_hse.so -lm -o native_handles
mpiexec -n 2 ./native_handles rt.inp
```

On macOS select the `.dylib` instead. Run in a **new** directory containing a
zero-initial-field full-k RT/GS fixture prepared for that MPI size. The sample
uses two dt values and different polarizations, closes the first state, creates
another from the survivor, and proves the old token stays invalid. It is a
composition demonstration, not a physically converged material calculation.

Shared model input, symmetry, initialization, I/O/timer and MPI ownership are
still bounded. With context_mpi_protocol1, initialization child communicators retire at state
close; unscoped children stay until `sp_detach`. This is not a promise of unlimited
create/close cycling. Close all
handles before detaching, then let the calling program finalize MPI. Heterogeneous models,
arbitrary communicator changes and concurrent/threaded entry are not enabled.
Stable token identity does not establish portable checkpoint files: inactive
legacy workspace is not canonical across independently initialized states.
Numerical restart verification compares accepted observables, exchange histories
and future propagation from restored state; it does not treat whole-file byte
identity across independent states as a physical requirement.
