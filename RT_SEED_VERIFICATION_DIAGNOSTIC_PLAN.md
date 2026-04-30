# RT Seed Verification Diagnostic Code Plan

## Objective

Verify two critical points affecting initial eigenstate quality:

1. **Point 1**: Do basis_functions_buffered.bin, wavefunctions_buffered.bin, eigenvalues_buffered.bin have **consistent state ordering**?
2. **Point 2**: Do H,S matrices reconstructed at RT startup match DC-LCFO definitions?

---

## Point 1: State Metadata Consistency Check

### Scope
Compare metadata across three buffered seed files after loading:
- basis_functions_buffered.bin (via read_buffered_basis_data)
- wavefunctions_buffered.bin (via read_wavefunctions_data)
- eigenvalues_buffered.bin (via read_buffered_esp_data)

### Implementation Location
File: `rt_dg_fragment.f90`
After: subroutine `read_buffered_basis_data` returns (line ~1900)
Before: subroutine `rt_dg_fragment_set_esp` called

### Verification Items
```
1. n_frag: consistency across basis/coef/esp files
2. nspin: consistency across basis/coef/esp files
3. nstate_frag: consistency across basis/coef/esp files
4. nstate_tot: consistency across basis/coef/esp files
5. n_mat: consistency across basis/coef/esp files
6. n_basis(:,:): consistency across basis/coef/esp files
7. index_basis(:,:,:): checksum comparison across basis/coef/esp
```

### Output Format
```
[SEED-CONSISTENCY-CHECK] Point 1 verification:
  n_frag:      basis=8  coef=8  esp=8  [OK]
  nspin:       basis=1  coef=1  esp=1  [OK]
  nstate_frag: basis=30 coef=30 esp=30 [OK]
  nstate_tot:  basis=256 coef=256 esp=256 [OK]
  n_mat:       basis=[64,0] coef=[64,0] esp=[64,0] [OK]
  n_basis checksum: basis=ABCD1234 coef=ABCD1234 [OK]
  index_basis checksum: basis=12345678 coef=12345678 [OK]
[SEED-CONSISTENCY-OK] All metadata aligned
```

**OR if mismatch:**
```
[SEED-CONSISTENCY-MISMATCH]
  nstate_tot: basis=256 but coef=254 (MISMATCH!)
  → State set does not correspond
[SEED-CONSISTENCY-FAIL]
```

### Code Stub
```fortran
subroutine check_seed_metadata_consistency(dg_frag, bdir_frag, enable_check)
  use communication, only: comm_is_root
  implicit none
  type(s_dg_fragment_rt), intent(in) :: dg_frag
  character(*), intent(in) :: bdir_frag
  logical, intent(in) :: enable_check
  
  if (.not. enable_check) return
  if (.not. dg_frag%use_buffered_coef_files) return
  if (dg_frag%id /= 0) return  ! Root only
  
  ! Read metadata from each file
  ! Compare checksums
  ! Output mismatch warnings
  
  if (any_mismatch) then
    write(*,'(1x,a)') "[SEED-CONSISTENCY-FAIL] State ordering inconsistency detected!"
    write(*,'(1x,a)') "  This explains large frozen-H residuals at RT startup."
  else
    write(*,'(1x,a)') "[SEED-CONSISTENCY-OK]"
  end if
end subroutine
```

---

## Point 2: H,S Eigenstate Problem Check at itt=0

### Scope
At first RT iteration (itt=0), compute and analyze:
1. **C† S C** (overlap normalization)
2. **C† H C** (Hamiltonian expectation)
3. **||H c_i - esp_i S c_i||** (generalized eigenvalue residual - MOST DIRECT)

### Expected Behavior
If coef are true eigenstates of generalized H,S problem:
- Diag(C† S C) ≈ 1.0 (normalized in S-metric)
- Offdiag_rms(C† S C) ≈ 0 (orthogonal in S-metric)
- Diag(C† H C) ≈ eigenvalues (correct energy expectations)
- **||H c_i - esp_i S c_i|| ≈ 0 (residuals tiny)**

### Implementation Location
File: `rt_dg_iteration.f90` or `rt_dg_fragment.f90`
Time: After H,S matrix assembly at itt=0, before first RK step

### Diagnostic Logic
```
For each state i=1..nstate_keep:
  1. Compute overlap metric: s_ii = <psi_i|S|psi_i> (should be 1.0)
  2. Compute cross-state: s_ij, i≠j (should be 0)
  3. Compute H-expec: h_ii = <psi_i|H|psi_i> (should match esp_i)
  4. Compute Rayleigh: R_i = h_ii / s_ii (should match esp_i / 1.0)
  5. Compute residual: r_i = ||H c_i - esp_i S c_i||  [MOST DIRECT]
  
  If r_i >> 1e-6 * |esp_i|:     c_i NOT eigenvector of (H,S)
  If offdiag_s >> 1e-3:         not orthogonal in S-metric
  If (h_ii - esp_i) >> 1e-3:    eigenvalue stale
```

### Output Format (DIAGNOSTIC ORDER)
```
[HS-EIGENSTATE] itt=0 spin=1 state=1-24:
  diag_CpSC_min  = 0.921  (should be 1.0)
  diag_CpSC_max  = 1.087
  offdiag_rms    = 0.052  (should be <1e-3)
  diag_CpHC_1    =-1.234  (sample state)
  esp_1          =-0.312  (mismatch!)
  rayleigh_1     =-1.339  (H_ii / S_ii)
  residual_rms   = 0.876  (||H c_i - esp_i S c_i||_rms)  ← MOST IMPORTANT
  residual_max   = 2.145  (should be << 1e-6)
[HS-EIGENSTATE-WARN] State 1: residual=2.145 >> 1e-6 → NOT eigenstate
[HS-EIGENSTATE-WARN] State 5: offdiag_S > 0.1 → not orthogonal to state 3
```

### Code Stub
```fortran
subroutine check_hs_eigenstate_at_startup(dg_frag, H_mat, S_mat, esp, itt, enable_check)
  use communication, only: comm_is_root, comm_summation
  implicit none
  type(s_dg_fragment_rt), intent(in) :: dg_frag
  real(8), intent(in) :: H_mat(:,:,:), S_mat(:,:,:), esp(:,:)
  integer, intent(in) :: itt
  logical, intent(in) :: enable_check
  
  real(8) :: s_ii, h_ii, rayleigh_i, residual_i
  real(8) :: diag_CpSC_min, diag_CpSC_max, offdiag_rms
  real(8) :: residual_rms, residual_max
  integer :: istate, io
  
  if (.not. enable_check) return
  if (itt /= 0) return  ! itt=0 only
  
  ! 1. Compute C† S C for normalization check
  ! 2. Compute C† H C for H expectation check
  ! 3. Compute ||H c_i - esp_i S c_i|| for eigenstate residual (MOST IMPORTANT)
  ! 4. Compute Rayleigh quotient R_i = H_ii / S_ii
  
  ! Output in order: diag_CpSC, offdiag_rms, diag_CpHC, rayleigh, residual_rms/max
  
  write(*,'(1x,a)') "[HS-EIGENSTATE] itt=0 spin=1 state=1-24:"
  write(*,'(1x,a,f10.4)') "  diag_CpSC_min  = ", diag_CpSC_min
  write(*,'(1x,a,f10.4)') "  diag_CpSC_max  = ", diag_CpSC_max
  write(*,'(1x,a,f10.4)') "  offdiag_rms    = ", offdiag_rms
  ! ... output sample states (first 3-5)
  write(*,'(1x,a,f10.4)') "  residual_rms   = ", residual_rms
  write(*,'(1x,a,f10.4)') "  residual_max   = ", residual_max
  
  if (residual_max > 1.0d-6) then
    write(*,'(1x,a)') "[HS-EIGENSTATE-FAIL] Large residuals → c_i NOT eigenstates of (H,S)"
    write(*,'(1x,a)') "  → Indicates DC→RT seed H,S mismatch or startup transform effect"
  else
    write(*,'(1x,a)') "[HS-EIGENSTATE-OK]"
  end if
end subroutine
```

---

## Activation

Add environment variables:
- `SALMON_DG_SEED_METADATA_CONSISTENCY_CHECK=1` → Enable Point 1 (metadata check)
- `SALMON_DG_HS_EIGENSTATE_CHECK=1` → Enable Point 2 (eigenstate residual check)

---

## Expected Outcomes

### Scenario A: Point 1 FAILS (metadata mismatch)
```
[SEED-CONSISTENCY-FAIL] nstate_tot: basis=256 coef=254
```
**Diagnosis**: DC→RT file metadata mismatch
**Root Cause**: DC-LCFO output error or RT seed reading bug
**Action**: Fix DC output or RT seed reading logic

### Scenario B: Point 2 FAILS (large residuals)
```
[HS-EIGENSTATE] residual_rms = 0.876, residual_max = 2.145 >> 1e-6
```
**Diagnosis**: c_i NOT eigenvectors of frozen H,S problem
**Root Cause**: Startup transform (Lowdin, stationary projection) altered coef without re-diagonalizing H,S
**Action**: Add startup H,S re-diagonalization or weaken transform

### Scenario C: Both PASS
**Diagnosis**: Problem lies in later stages
**Possibilities**:
- RK4 propagation numerical issue
- Frozen H matrix has nonphysical features
- esp not updated correctly in startup projection
**Action**: Enable other diagnostics (Rayleigh quotient, C†SC time evolution check)

---

## Implementation Summary: Key Points

### Point 1: Metadata 7-Item Checklist
```
✓ n_frag match
✓ nspin match
✓ nstate_frag match
✓ nstate_tot match
✓ n_mat match
✓ n_basis checksum match
✓ index_basis checksum match
```

### Point 2: Generalized Eigenstate Diagnostic (5-metric order)
1. **diag(C† S C)** → normalization in S-metric
2. **offdiag_rms(C† S C)** → orthogonality in S-metric
3. **diag(C† H C)** → H expectation values
4. **Rayleigh = H_ii / S_ii** → should match esp_i
5. **||H c_i - esp_i S c_i||** ← **MOST DIRECT INDICATOR** (main figure of merit)

If residual_max >> 1e-6: c_i NOT eigenstates → FAIL diagnosis

---

## Files to Modify

1. **rt_dg_fragment.f90**: 
   - Insert `check_seed_metadata_consistency` after metadata setup (line ~1900)
   - Add env var parsing for `SALMON_DG_SEED_METADATA_CONSISTENCY_CHECK`

2. **rt_dg_iteration.f90**: 
   - Insert `check_hs_eigenstate_at_startup` at itt==0 (line ~150-160)
   - Add env var parsing for `SALMON_DG_HS_EIGENSTATE_CHECK`

---

## Timeline

1. Code integration: ~30 min
2. Rebuild: ~5-10 min
3. Test run: ~10 min
4. Analysis: ~10 min
