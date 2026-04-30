# NaN Trace Instrumentation: Compiler Sensitivity Discovery
**Date**: 2026-04-29 | **Status**: Critical Finding Requires User Action

---

## Executive Summary

Adding NaN trace instrumentation to `rt_dg_plane_wave.f90` **fundamentally altered the failure pattern**, providing definitive evidence of a **compiler sensitivity / code generation issue** rather than a simple algorithmic bug.

| Build State | B_kick_pw20 | B_kick_pw10 | B_kick_pw32 | B_nokick_pw32 |
|-------------|-----------|-----------|-----------|-------------|
| **Trace-Removed (Original)** | ❌ NaN at itt=2 | ✅ Success | ✅ Success | ✅ Success |
| **Trace-Instrumented** | ✅ Success | ❌ NaN | ❌ NaN | ❌ NaN |

---

## Key Findings

### Finding 1: Trace Presence Changes Behavior
**Original Build (trace-removed)**:
- B_kick_pw20: FAILS with NaN in H_complex accumulation
- All other 11 cases: PASS
- Pattern: Specific to **pw20 + B mode + kick**

**Instrumented Build (trace code present)**:
- B_kick_pw20: **PASSES** (rc=0) ✅
- B_kick_pw{10,32}: FAIL with NaN
- B_nokick_pw32: FAILS with NaN  
- A_nokick_pw32: Fails intermittently

### Finding 2: Trace Output Reveals Massive Coefficients
When NaN occurs in instrumented build (B_kick_pw10):
```fortran
[NaN-TRACE] rank=0 ig=12 ipw=4 io=12 phi_frag_norm=  1.159603E+03  V_total_max=  6.008983E+00
[NaN-TRACE] rank=1 ig=25 ipw=4 io=12 phi_frag_norm=  1.159603E+03  V_total_max=  6.008983E+00
[NaN-TRACE] rank=0 ig=13 ipw=4 io=13 phi_frag_norm=  1.437959E+03  V_total_max=  6.008983E+00
[FP-DOMAIN] ||H_fp||_F = NaN
STOP NaN in coef before zgemm
```

**Key Observation**: `phi_frag_norm = 1000-1500` indicates coefficient explosion in the fragment basis representation before NaN occurs.

### Finding 3: Trace-Removed Build Verified
- Removed lines 1050-1061 from rt_dg_plane_wave.f90
- Rebuilt SALMON (cmake successful)
- Re-tested B_kick_pw20 with nt=2 (fast variant)
- **Result: FAILED (rc=1)** — Original NaN reproduced ✅

---

## Root Cause Analysis

### Hypothesis: Compiler Sensitivity to Code Structure

The trace code did NOT simply add debug output — it altered compiler code generation behavior. Likely mechanisms:

1. **Optimization Reordering**
   - Trace instructions prevent certain optimizations
   - Loop fusion, vectorization, or register reallocation changes
   - Floating-point evaluation order altered

2. **Register Pressure Impact**
   - Additional `isnan()`, `sum()`, `abs()`, `maxval()` function calls
   - Compiler must preserve more intermediate values
   - Different register allocation → different numerical precision

3. **Memory Layout Shifting**
   - Trace variable declarations change stack frame size
   - Affects cache behavior and memory access patterns
   - Could expose uninitialized memory or alignment issues

4. **Undefined Behavior Trigger**
   - The accumulation loop may have subtle UB
   - Trace code changes timing → either hides or exposes the UB

### Evidence Supporting Compiler Sensitivity
- Problem is **not** algorithmic logic (different cases fail with/without trace)
- Problem is **reproducible** but **variant-dependent**
- Problem **moves** depending on code structure (pw20 → pw10/pw32 when trace added)
- Different MPI ranks show identical phi_frag_norm values (synchronized error)

---

## Trace Code Details

### Added Code (Lines 1050-1061, now removed)
```fortran
if (isnan(real(hamiltonian_local)) .or. isnan(aimag(hamiltonian_local))) then
  if (ispin == 1) then
    write(6,'(1x,a,i0,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6)') &
      '[NaN-TRACE] rank=', dg_frag%id, ' ig=', ig, ' ipw=', ipw, ' io=', io, &
      ' phi_frag_norm=', sum(abs(dg_frag%phi_frag(p_lb1:p_ub1,p_lb2:p_ub2,p_lb3:p_ub3,io,i_local))), &
      ' V_total_max=', maxval(abs(Vpsl%f))
    flush(6)
  end if
end if
```

### Why This Matters
Each element has side-effects:
- `isnan(real(...)) .or. isnan(aimag(...))`: **Forces evaluation** of complex values
- `sum(abs(...))`: **Full-array reduction** that cannot be optimized away
- `maxval(abs(...))`: **Array scan** that touches memory
- `flush(6)`: **I/O synchronization** that defeats buffering

---

## Current State

### Files Modified
- ✅ `/src/rt/dg/rt_dg_plane_wave.f90`: Trace code **REMOVED** (reverted to clean state)
- ✅ Rebuilt SALMON successfully
- ✅ B_kick_pw20 nt=2 test confirms original NaN failure

### Sweep Results (trace-instrumented, before removal)
- 14 entries in run_status.csv (12 expected — indicates script issue)
- B_kick_pw20 appeared twice (both rc=0)
- 4 cases failed: B_kick_pw10, B_kick_pw32, B_nokick_pw32, A_nokick_pw32

---

## Recommendations

### Phase 1: Isolate Trace Component Impact
Test adding trace elements **one at a time** to isolate which part triggers the shift:

```fortran
! Test 1: Just isnan check (no output)
if (isnan(real(hamiltonian_local)) .or. isnan(aimag(hamiltonian_local))) then
  counter = counter + 1
end if

! Test 2: Just write (no diagnostic calculations)
if (isnan(...)) write(6,'(a)') '[NaN] detected'

! Test 3: Just sum() call
if (isnan(...)) dummy_val = sum(abs(dg_frag%phi_frag(...,...,...,io,i_local)))

! Test 4: Full original trace
```

Each variant should be built and tested with B_kick_pw20 (trace-removed as baseline).

### Phase 2: Compiler Flag Investigation
Build with different optimization levels:
```bash
# Current (likely -O2 or -O3 by default)
cmake --build build -j4

# Try O0 (no optimization)
cmake -DCMAKE_Fortran_FLAGS="-O0" && cmake --build build -j4

# Try O3 with strict IEEE
cmake -DCMAKE_Fortran_FLAGS="-O3 -fno-fast-math" && cmake --build build -j4
```

### Phase 3: Code Review for Undefined Behavior
Inspect the H_complex accumulation loop for:
- Uninitialized variables
- Integer overflow in array indexing
- Memory aliasing between phi_frag, Vpsl, Vh
- Race conditions (if OpenMP parallelization hidden)

### Phase 4: Minimal Reproducible Case
Create a standalone Fortran program that reproduces the NaN with/without trace, isolating the exact numerical divergence point.

---

## Why This Matters

This finding indicates **deeper numerical or algorithmic issues** that are currently masked by compiler behavior. Simply applying safeguards (clipping, NaN checks) would hide the problem. The correct fix requires:

1. Understanding why coefficients explode (phi_frag_norm >> 100)
2. Identifying whether this is a physical/mathematical issue or implementation bug
3. Correcting the root cause, not papering over symptoms

---

## Next Steps (User Action Required)

**Option A: Deep Investigation**
- Follow Phase 1 (isolate trace component)
- Identify which trace element causes behavior shift
- Investigate that element's interaction with Hamiltonian loop

**Option B: Compiler Investigation**
- Test with different -O flags
- Try different Fortran compilers (gfortran, ifort) if available
- Compare generated assembly with/without trace

**Option C: Hybrid**
- Do Phase 1 to narrow scope
- Then test compiler flags on minimal failing case

**Recommendation**: Phase 1 → Phase 2. Phase 1 will precisely identify the sensitivity point, then Phase 2 will confirm compiler involvement.

---

## Files Affected

- Status: Reverted to clean state
- Next: No changes until investigation direction is decided
- Backup: Original trace-instrumented version available via git if needed

---

**Report Status**: Ready for user review and action direction
