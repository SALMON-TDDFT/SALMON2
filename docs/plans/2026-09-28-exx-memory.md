# Native EXX memory reduction implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reduce native mesh EXX working/storage memory without a new physical approximation.
**Architecture:** Remove native midpoint-factor concatenation/output/cached-action copies, then represent source ACE by exact-nonzero W and its metric inverse. Retain dense fallback and original mesh propagation.
**Tech Stack:** Fortran/MPI/BLAS/LAPACK/FFTW/ScaLAPACK; Python validation.

### Task1: Exact buffer reduction
Files: src/xc/hse_native.f90; testsuites/653_functional/test_exx_memory.py.
Write a regression that detects the redundant native midpoint/cached-action/output allocation in a temporary instrumented executable and compares native results to the frozen executable. RED must show redundant storage. Implement endpoint-sum midpoint action and direct accumulation to mesh; remove RT-only cached_action; transfer local into cached_source after energy calculation. Run source/native/fallback and DC regression. Expected numerical equivalence, eliminated buffers. Commit.

### Task2: Sparse source ACE
Files: src/xc/hse_ace.f90, src/xc/exx_orbitals.f90, src/xc/hse_native.f90; testsuites/651_hybrid_exchange/ace/exchange_driver.f90.
Write RED packed=true build/apply tests comparing dense factors on arbitrary complex targets and nonorthogonal support vectors, zero action, rejected state and empty/uneven MPI layouts. Extend ACE state with sparse column offsets/rows/values, metric inverse and dimensions. Factorization validation unchanged; source native builds select packing. Apply via streamed target columns and two overlap/action reductions across orbital/spatial groups; no mesh-by-all-global-orbitals temporary. Clear both representations on failure/rebuild. Release masked sources after fixed-ion source-ACE refresh; transport previous remains. Expected27 exchange configurations and native source/fallback regressions pass. Commit.

### Task3: Memory measurements and review
Freeze old executable/hash before edits. Compare old/new16-step results at8/64/128H2 using same validated seeds and existing old measurements; measure only new source runs, sequentially. Archive exact output/input/rank RSS and note theoretical allocated bytes vs process peak. Run independent fresh review using requesting-code-review skill, fix important findings RED/GREEN. Build ScaLAPACK, run relevant unit/native/DC regression. Document achieved RSS/time and remaining dense storage, commit; no push requested.

## Review focus
Midpoint state lifetime and collective call order; native/full/LCFO routes remain distinct; sparse W complex conjugation and dv; source-only release safe for cached refresh and full-action DC; empty owners; fallback invalidates all representations; exact zeros only; metric storage scaling disclosed.
