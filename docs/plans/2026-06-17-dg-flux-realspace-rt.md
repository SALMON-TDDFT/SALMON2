# DG Flux Real-Space RT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Initialize RT from Flux-inclusive DC-LCFO occupied eigenstates and propagate them as real-space wavefunctions with a Flux-inclusive Hamiltonian.

**Architecture:** Add an explicit RT mode that reconstructs occupied `spsi%zwf` from DC-LCFO-Flux files, then uses the existing real-space RT Taylor path. Refactor the existing GS real-space Flux correction into a reusable module and add a complex `zwf` version so RT Hamiltonian applications include the same DG Flux boundary operator.

**Tech Stack:** Fortran 2008, SALMON RT/GS modules, MPI, OpenMP, CMake, EigenExa-enabled SALMON build.

---

### Task 1: Add The Mode Flag With Guard Rails

**Files:**
- Modify: `src/io/salmon_global.f90:151-165`
- Modify: `src/io/inputoutput.f90:311-323`
- Modify: `src/io/inputoutput.f90:689-790`
- Modify: `src/io/inputoutput.f90:1139-1187`
- Modify: `src/io/inputoutput.f90:1360-1370`
- Modify: `src/io/inputoutput.f90:2283-2290`
- Modify: `src/io/inputoutput.f90:2798-3180`

**Step 1: Add the new global variable**

Add near the existing DG RT propagation variables:

```fortran
character(1) :: yn_dg_flux_realspace_rt
```

**Step 2: Add it to the `&propagation` namelist**

Append it after `yn_dg_fragment_rt`:

```fortran
& yn_dg_fragment_rt, &
& yn_dg_flux_realspace_rt, &
& time_integrator_dg_fragment, &
```

**Step 3: Add default, read, broadcast, log, and argument checks**

Use default:

```fortran
yn_dg_flux_realspace_rt = 'n'
```

Broadcast:

```fortran
call comm_bcast(yn_dg_flux_realspace_rt,nproc_group_global)
```

Log:

```fortran
write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_dg_flux_realspace_rt', yn_dg_flux_realspace_rt
```

Validate with:

```fortran
call yn_argument_check(yn_dg_flux_realspace_rt)
```

**Step 4: Add unsupported-case stops**

Add a new validation block after the existing DG-Fragment RT checks:

```fortran
if(yn_dg_flux_realspace_rt == 'y') then
  if(theory /= 'tddft_pulse' .and. theory /= 'tddft_response' .and. theory /= 'single_scale_maxwell_tddft') &
  & stop "DG-Flux real-space RT: theory must be tddft_pulse, tddft_response, or single_scale_maxwell_tddft"
  if(yn_conventional_from_dcdft /= 'y') &
  & stop "DG-Flux real-space RT: set yn_conventional_from_dcdft='y' to reconstruct DC-LCFO wavefunctions"
  if(yn_dg_fragment_rt == 'y') &
  & stop "DG-Flux real-space RT: disable yn_dg_fragment_rt; this mode uses real-space RT propagation"
  if(yn_spinorbit == 'y') &
  & stop "DG-Flux real-space RT: spin-orbit mode is not implemented"
  if(product(num_kgrid) /= 1) &
  & stop "DG-Flux real-space RT: # of k-points must be 1"
  if(nproc_k /= 1) &
  & stop "DG-Flux real-space RT: nproc_k must be 1"
end if
```

**Step 5: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds. No runtime behavior changes yet because the flag defaults to `n`.

**Step 6: Commit**

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90
git commit -m "feat: add DG flux real-space RT mode flag"
```

---

### Task 2: Move The Real-Space Flux Operator To A Shared Module

**Files:**
- Create: `src/common/dc_lcfo_flux_operator.f90`
- Modify: `src/gs/conjugate_gradient.f90:1030-1311`
- Modify: the CMake source list that owns `src/common` modules

**Step 1: Create the shared module**

Move the helper types and routines used by `apply_dc_lcfo_flux_hpsi_rwf` out of `src/gs/conjugate_gradient.f90` into:

```fortran
module dc_lcfo_flux_operator
  implicit none
  private
  public :: apply_dc_lcfo_flux_hpsi_rwf
contains
  ! moved helper routines and apply_dc_lcfo_flux_hpsi_rwf
end module dc_lcfo_flux_operator
```

Keep the routine body unchanged in this task. The goal is a pure refactor.

**Step 2: Import the shared module from GS CG**

At the top of `src/gs/conjugate_gradient.f90`, add:

```fortran
use dc_lcfo_flux_operator, only: apply_dc_lcfo_flux_hpsi_rwf
```

Remove the moved helper definitions from the bottom of the file.

**Step 3: Update CMake ordering**

Add `common/dc_lcfo_flux_operator.f90` before `gs/conjugate_gradient.f90` and before any future RT file that uses it.

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 5: Run a GS smoke test using the existing Flux DC path**

Use the already working diamond Flux input in a temporary copy:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs > /tmp/dg_flux_refactor_gs.log 2>&1
```

Expected:

- The run reaches the DC-LCFO-Flux path.
- No new compile/link/runtime error from moving the Flux operator.

**Step 6: Commit**

```bash
git add src/common/dc_lcfo_flux_operator.f90 src/gs/conjugate_gradient.f90 src/common/CMakeLists.txt
git commit -m "refactor: share DC LCFO flux operator"
```

---

### Task 3: Add A Complex `zwf` Flux Operator

**Files:**
- Modify: `src/common/dc_lcfo_flux_operator.f90`

**Step 1: Add a public complex entry point**

Add:

```fortran
public :: apply_dc_lcfo_flux_hpsi_zwf
```

**Step 2: Implement `apply_dc_lcfo_flux_hpsi_zwf`**

Duplicate the `rwf` routine structure, but use complex buffers and `zwf`:

```fortran
complex(8), allocatable :: buf_send(:,:,:,:,:), buf_recv(:,:,:,:,:)
```

Replace every `psi%rwf(...,io,1,1)` read with:

```fortran
psi%zwf(ixg,iyg,izg,ispin,io,1,1)
```

Replace every `hpsi%rwf(...) = hpsi%rwf(...) + flux_coef * buf_recv(...)` update with:

```fortran
hpsi%zwf(il(1),il(2),il(3),ispin,io,1,1) = &
& hpsi%zwf(il(1),il(2),il(3),ispin,io,1,1) + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
```

Keep the same geometry, tag, and cache semantics as the real routine.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/common/dc_lcfo_flux_operator.f90
git commit -m "feat: add complex DC LCFO flux operator"
```

---

### Task 4: Route RT Hamiltonian Applications Through The Flux Operator

**Files:**
- Modify: `src/common/hamiltonian.f90:27-809`
- Modify: `src/rt/taylor.f90:21-99`
- Modify: `src/rt/time_evolution_step.f90:19-208`

**Step 1: Add optional DC data to `hpsi` or add a wrapper**

Preferred minimal interface:

```fortran
SUBROUTINE hpsi(tpsi,htpsi,info,mg,V_local,system,stencil,srg,ppg,ttpsi,dc_flux)
  type(s_dcdft), optional, intent(in) :: dc_flux
```

After the normal `zwf` Hamiltonian work finishes, add:

```fortran
if (present(dc_flux)) then
  if (yn_dg_flux_realspace_rt == 'y') then
    call apply_dc_lcfo_flux_hpsi_zwf(mg,system,info,stencil,dc_flux,tpsi,htpsi)
  end if
end if
```

Do the analogous real branch only if needed for testing.

**Step 2: Thread DC data through Taylor**

Extend `taylor`:

```fortran
subroutine taylor(..., rt, dc_flux)
  type(s_dcdft), optional, intent(in) :: dc_flux
```

Pass `dc_flux` to every `hpsi` call.

**Step 3: Thread DC data from `time_evolution_step`**

Add an optional `dc_flux` argument to `time_evolution_step` and pass it to `taylor`.

**Step 4: Update all call sites**

Search:

```bash
rg -n "call time_evolution_step|call taylor|call hpsi\\("
```

Update only the RT path that needs `dc_flux`. All existing callers should compile without changes because the new arguments are optional.

**Step 5: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 6: Commit**

```bash
git add src/common/hamiltonian.f90 src/rt/taylor.f90 src/rt/time_evolution_step.f90
git commit -m "feat: support Flux operator in real-space RT hpsi"
```

---

### Task 5: Reconstruct Occupied Flux Eigenstates Into `spsi%zwf`

**Files:**
- Modify: `src/gs/dc/lcfo.f90:766-950`
- Modify: `src/rt/initialization_rt.f90:287-303`

**Step 1: Split reconstruction into a mode-aware helper**

Add a new public routine near `init_conventional_from_dcdft`:

```fortran
subroutine init_flux_realspace_rt_from_dcdft(lg,mg,system,info,spsi)
```

It should share the file reading and basis reconstruction logic with `init_conventional_from_dcdft`, but force the Flux-compatible basis data path and write occupied orbitals into `spsi%zwf`.

**Step 2: Restrict reconstruction to occupied orbitals**

Loop over the local RT orbital range:

```fortran
do io = info%io_s, info%io_e
  if (io > system%no) cycle
  ...
end do
```

Use the DC coefficient column corresponding to the same occupied state.

**Step 3: Preserve normalization**

After reconstruction, compute local norms and reduce them. Stop if a norm is non-positive. Print the min/max reconstructed norm on the root rank for the new mode.

**Step 4: Hook initialization**

In `src/rt/initialization_rt.f90`, change the restart branch:

```fortran
if (yn_dg_flux_realspace_rt == 'y') then
  call init_flux_realspace_rt_from_dcdft(lg,mg,system,info,spsi_in)
else if (yn_dg_fragment_rt == 'y' .and. yn_dg_fragment_from_dcdft == 'y') then
  ...
```

**Step 5: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 6: Commit**

```bash
git add src/gs/dc/lcfo.f90 src/rt/initialization_rt.f90
git commit -m "feat: initialize DG flux RT from occupied real-space states"
```

---

### Task 6: Pass DC Fragment Context Into Real-Space RT

**Files:**
- Modify: `src/rt/main_tddft.f90:116-186`
- Modify: `src/rt/main_tddft.f90:127-131`
- Modify: `src/rt/main_tddft.f90:323-454`

**Step 1: Initialize DC fragment metadata for the new mode**

Reuse the existing DC-LCFO seed setup enough to populate `dc` with fragment geometry, communicator, and `id_frag`.

Do not build coefficient-space DG matrices when `yn_dg_flux_realspace_rt == 'y'`.

**Step 2: Pass `dc` to real-space propagation**

Change the ordinary call:

```fortran
call time_evolution_step(..., singlescale)
```

to:

```fortran
if (yn_dg_flux_realspace_rt == 'y') then
  call time_evolution_step(..., singlescale, dc_flux=dc)
else
  call time_evolution_step(..., singlescale)
end if
```

Use the final argument name chosen in Task 4.

**Step 3: Keep coefficient-space DG-RT disabled**

Ensure `yn_dg_flux_realspace_rt='y'` does not call:

- `init_dg_fragment_rt`
- `calculate_hamiltonian_matrix`
- coefficient Taylor/RK propagators
- DG coefficient density reconstruction

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 5: Commit**

```bash
git add src/rt/main_tddft.f90
git commit -m "feat: drive real-space RT with DG flux context"
```

---

### Task 7: Add A Minimal Diamond Flux Real-Space RT Input

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_flux_realspace_smoke`

**Step 1: Copy the existing RT smoke input**

Start from the existing short RT input used for the diamond Flux case.

**Step 2: Change mode settings**

Set:

```text
yn_conventional_from_dcdft = 'y'
yn_dg_fragment_rt = 'n'
yn_dg_flux_realspace_rt = 'y'
```

Keep the same impulse field, fragment layout, and DC data directory assumptions.

**Step 3: Use a tiny first run**

Set `nt` to a small value such as `2` or `5`.

**Step 4: Commit**

```bash
git add samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_flux_realspace_smoke
git commit -m "test: add DG flux real-space RT smoke input"
```

---

### Task 8: Validate The New Path

**Files:**
- Test only; no source edits expected unless failures identify a bug.

**Step 1: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 2: Run DC from scratch if needed**

Run:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs > /tmp/dg_flux_realspace_gs.log 2>&1
```

Expected:

- DC converges.
- `data_dcdft/fragments/*/wavefunctions.bin` exists.
- DC-LCFO-Flux diagonalization/export is reported.

**Step 3: Run the RT smoke**

Run:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_flux_realspace_smoke > /tmp/dg_flux_realspace_rt_nt5.log 2>&1
```

Expected:

- Initialization prints the new DG-Flux real-space RT mode.
- Occupied real-space reconstruction norm min/max is close to 1.
- The run reaches at least `itt=2`.
- No coefficient-space DG Hamiltonian build messages appear.

**Step 4: Inspect the current**

Extract `J_tot,z`, `J_para,z`, and `J_dia,z` from the log or current output.

Expected qualitative behavior:

- Initial `J_tot,z` is set by the kick.
- `J_para,z` begins cancelling `J_dia,z`.
- `J_tot,z` moves toward zero, not toward a persistent `J_dia` offset.

**Step 5: Commit any fixes**

If a bug is found, fix it in the smallest source set and commit with a message describing the failing validation.

---

### Task 9: Run The Longer Impulse Check

**Files:**
- Test only, unless a physics mismatch identifies a bug.

**Step 1: Create a temporary long input**

Use the smoke input and increase `nt` enough to see the first few oscillations.

**Step 2: Run RT**

Run:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < /tmp/inputfile_rt_flux_realspace_long > /tmp/dg_flux_realspace_rt_long.log 2>&1
```

Expected:

- `J_tot,z` oscillates around zero.
- The negative peak should not grow much larger than the positive peak.
- No large unphysical `J_tot,x` or `J_tot,y` appears for a z kick.

**Step 3: Record the result**

Add a short note to the final response with:

- `J_tot,z` first value.
- approximate max/min over the sampled window.
- whether the response is centered near zero.

Do not claim full long-time drift quality yet; that is explicitly a later diagnostic.
