# Default Hermitian Exchange Implementation Plan

**Goal:** Prefer eligible crystal FFT reuse, otherwise default to exact Hermitian exchange pairs, including finite-temperature DC orbital MPI.

**Architecture:** Add a bounded Gamma orbital-pair action using raw target orbitals and separate occupation/2 weights. One convolution updates both directions with the appropriate source weight; no division by occupations, so empty states are safe. Preserve reference opt-out and existing wave/SR route. Eligible dense-k crystal reuse takes precedence in automatic selection. Existing submitted binaries/jobs stay unchanged.

**Tech Stack:** Fortran, MPI orbital/spatial communicators, existing distributed FFT, Python regression harness, CMake/CTest.

1. Add failing tests for fractional/zero occupations, unequal/empty orbital ownership, complex orbitals, pair counts and collective invalid-input rejection. Use the actual production module with existing FFT/communication dependencies.
2. Implement full-grid Gamma pair action, scalar IEEE validation, bounded four-column FFT scratch and reverse orbital reductions. Compare to the full ordered-pair reference on identical grids/kernel.
3. Connect only canonical Gamma source/target layouts with explicit occupations; retain other target actions. Default pair switch ON without forcing SR/wave-only settings. Resolve eligible crystal reuse before scalar pairs.
4. Test selection policies and production DC SCF on/off with orbital MPI; inspect per-rank exits and memory. Run prohibition source check before builds and all application-object audit after builds. Do not build historical sources.
5. Record numerical differences, FFT counts, memory scope and limitations under canonical calc2026. No remote production changes or submission without further applicable authorization.
