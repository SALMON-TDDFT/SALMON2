# Tau Mixing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add opt-in `tau` mixing for tau-using XC functionals so meta-GGA SCF can stabilize the `tau -> vtau -> Hamiltonian` feedback loop without affecting non-tau functionals.

**Architecture:** Compute raw `tau` before XC, store it in SCF-side mixing state, and mix it alongside `rho`. For `simple` and `simple_potential`, use independent linear `tau` mixing with `tau_mixrate`. For `broyden` and `pulay`, append a weighted `tau` block to the `rho` mixing vector and recompute `vtau` fresh from the mixed `tau`.

**Tech Stack:** Fortran 2008, SALMON GS SCF path, source grep tests, local GCC build, Wisteria MPI runtime validation.

---

### Task 1: Add failing source tests for the public interface

**Files:**
- Modify: `tests/source/CMakeLists.txt`
- Create: `tests/source/check_tau_mixing_input_keys.sh`
- Create: `tests/source/check_tau_mixing_broyden_pack.sh`
- Create: `tests/source/check_tau_mixing_xc_override.sh`

**Step 1: Write the failing tests**

- `check_tau_mixing_input_keys.sh` should assert:
  - `inputoutput.f90` accepts `yn_tau_mixing`, `tau_mixrate`, `tau_metric_weight`
  - `variables.log` output includes the three keys
- `check_tau_mixing_broyden_pack.sh` should assert:
  - `s_mixing` owns `tau` history/state
  - `wrapper_broyden` and `pulay` pack/unpack a `tau` block
  - `tau_metric_weight` is applied through scaling
- `check_tau_mixing_xc_override.sh` should assert:
  - `exchange_correlation` accepts a `tau_override`
  - the old `mix_xc_operator_payload`-style `vtau` mixing path is absent

**Step 2: Run tests to verify they fail**

Run:
```bash
sh tests/source/check_tau_mixing_input_keys.sh
sh tests/source/check_tau_mixing_broyden_pack.sh
sh tests/source/check_tau_mixing_xc_override.sh
```

Expected: FAIL because none of the new `tau` mixing hooks exist yet.

### Task 2: Add input keys and logging

**Files:**
- Modify: `src/io/inputoutput.f90`

**Step 1: Write minimal implementation**

- Add defaults:
  - `yn_tau_mixing = 'n'`
  - `tau_mixrate`
  - `tau_metric_weight`
- Parse, validate, broadcast, and emit them to `variables.log`
- Guard usage so non-`tau` XC does not allocate or run extra work

**Step 2: Run targeted tests**

Run:
```bash
sh tests/source/check_tau_mixing_input_keys.sh
```

Expected: PASS

### Task 3: Extract reusable raw `tau` calculator

**Files:**
- Modify: `src/xc/salmon_xc.f90`
- Optionally create: `src/xc/xc_aux_fields.f90`

**Step 1: Write failing source expectation**

- Extend `check_tau_mixing_xc_override.sh` so it requires a public/helper `calc_tau` path callable from SCF

**Step 2: Verify the test fails**

Run:
```bash
sh tests/source/check_tau_mixing_xc_override.sh
```

Expected: FAIL because `calc_tau` is still internal to `exchange_correlation`

**Step 3: Write minimal implementation**

- Move `calc_tau` out of the internal procedure scope
- Keep the existing math and MPI behavior
- Preserve the `j` computation path needed by tau-using XC code

**Step 4: Re-run the test**

Run:
```bash
sh tests/source/check_tau_mixing_xc_override.sh
```

Expected: still FAIL because `tau_override` is not yet wired into XC

### Task 4: Add `tau` state to SCF mixing

**Files:**
- Modify: `src/common/structures.f90`
- Modify: `src/gs/mixing.f90`
- Modify: `src/io/checkpoint_restart.f90`

**Step 1: Write minimal data model**

- Extend `s_mixing` with `tau_in/tau_out` state/history
- Allocate only when `yn_tau_mixing='y'` and the active XC uses `tau`
- Add checkpoint/restart read/write for `tau` history

**Step 2: Run targeted tests**

Run:
```bash
sh tests/source/check_tau_mixing_broyden_pack.sh
```

Expected: still FAIL because packing logic is not implemented yet.

### Task 5: Implement `simple` and `simple_potential` tau mixing

**Files:**
- Modify: `src/gs/scf_iteration.f90`
- Modify: `src/gs/mixing.f90`

**Step 1: Write minimal implementation**

- In the GS loop, compute raw `tau` before XC
- For `simple` and `simple_potential`, linearly mix `tau` with `tau_mixrate`
- Pass mixed `tau` into XC

**Step 2: Run targeted tests**

Run:
```bash
sh tests/source/check_tau_mixing_xc_override.sh
```

Expected: partial progress; Broyden/Pulay test still FAILS.

### Task 6: Implement Broyden/Pulay composite vector mixing

**Files:**
- Modify: `src/gs/mixing.f90`
- Modify: `src/common/broyden.f90` only if wrapper-only packing proves insufficient

**Step 1: Write minimal implementation**

- Pack `rho` and `sqrt(tau_metric_weight) * tau` into one vector
- Reuse existing Broyden solver on the packed vector
- Unpack mixed `rho` and `tau`
- Mirror the same strategy for Pulay

**Step 2: Run targeted tests**

Run:
```bash
sh tests/source/check_tau_mixing_broyden_pack.sh
```

Expected: PASS

### Task 7: Make XC consume mixed `tau` and remove `vtau` mixing

**Files:**
- Modify: `src/xc/salmon_xc.f90`
- Modify: `src/gs/mixing.f90`
- Modify: `tests/source/check_builtin_r2scan_vtau_damping.sh`
- Modify: `tests/source/check_tau_mixing_xc_override.sh`

**Step 1: Write minimal implementation**

- Add optional `tau_override` to `exchange_correlation`
- If present, skip raw `tau` recomputation inside XC
- Recompute `vtau` fresh from mixed `tau`
- Remove remaining SCF-side `vtau` mixing helper

**Step 2: Run targeted tests**

Run:
```bash
sh tests/source/check_tau_mixing_xc_override.sh
sh tests/source/check_builtin_r2scan_vtau_damping.sh
```

Expected: PASS

### Task 8: Local verification

**Files:**
- Modify: none unless tests fail

**Step 1: Run source tests**

Run:
```bash
ctest --test-dir build-gcc15-env -R 'check_(tau_mixing_input_keys|tau_mixing_broyden_pack|tau_mixing_xc_override|builtin_r2scan_vtau_damping|builtin_r2scan_loop_perf_structure)' --output-on-failure
```

Expected: PASS

**Step 2: Rebuild**

Run:
```bash
cmake --build build-gcc15-env -j4 --target salmon
```

Expected: PASS

### Task 9: Wisteria runtime validation

**Files:**
- Modify: test input copies under Wisteria campaign directories only if needed

**Step 1: Push and rebuild on Wisteria**

Run:
```bash
git push origin feature/stress-tensor
```

Then on the persistent Wisteria session:
```bash
git pull --ff-only origin feature/stress-tensor
cd /work/go33/o33000/src/stress/SALMON2/build
make -j4
make install
```

Expected: PASS

**Step 2: Re-run primitive exact case**

Run on Wisteria in a fresh copied directory derived from:
`/work/go33/o33000/salmon/virial/pSi/test/campaign_20260406/primitive_r2scan_vtau_exact_a555_n24_k8`

Expected:
- no regression in startup
- compare `iter ~ 40, 60, 100+`
- determine whether plateau onset is delayed or removed

**Step 3: Run conventional sanity**

Use a conventional `r2scan` Si GS case with the same binary and confirm no obvious SCF regression.
