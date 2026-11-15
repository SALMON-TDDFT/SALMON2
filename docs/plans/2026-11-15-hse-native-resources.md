# Native HSE resource routing

User authorization: connect independent HSE resources to native TDDFT, preserve fast-HSE results, and test noninterference.

Preserve the existing optimized fast-PR dispatch. Do not replace it with the standalone generic kernel primitive, and do not switch process-global active pointers. Give the explicit RT owner its own kernel and symmetry resources. Bind that owner to the native system as a borrowed reference, and pass it through the real Hamiltonian, XC/energy, Taylor and PT-CN call sites. Default calls remain supported. The first supported independent native route is serial, full-grid, non-spin-polarized, full-support HSE with common immutable process configuration; spatial/Wannier/GPU and different global configurations must be rejected or remain on the legacy path. No unrestricted multi-instance TDDFT API is implied.

1. Add a failing real refresh/action/interleaving probe against the old audited build.
2. Move fast kernel ownership into explicit states; preserve call order, pair/orbit dispatch and ACE/history algorithms. Implement explicit destruction without copying FFT owners.
3. Route native TDDFT call sites via a system-owned reference. Test initialization and legacy sentinels, interleaved sources, refresh, trial histories, close/reopen, and unsupported paths.
4. Prepare new external source/build with updated inventory hashes; audit source and all objects. Compare the full native TDDFT current/energy/prediction traces with the immutable fast reference (32 cases, multiple k meshes, dt, propagators).
5. Run source tests and packaging validation, update English/Japanese/AI documentation with the verified scope and remaining global-configuration limits.

## Completed local verification

See ../hse-native-owners-validation-2026-11-15.md. Explicit-owner and default-owner TDDFT each pass 32 reference cases; native owner actions are isolated versus separate-process references and a shared-kernel mutation is detected. Full TDDFT context interleaving remains outside this verified stage. No remote CI was started, per the user’s Actions-budget instruction.
