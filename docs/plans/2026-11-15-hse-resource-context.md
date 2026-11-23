# Independent HSE kernel resource context

**Goal:** Establish independently owned configuration, FFT kernel and accepted ACE/source data, with explicit refresh and release, before routing the full native adapter.

**Architecture:** An opt-in small Fortran adapter links to the pinned fast-PR kernels. It never changes the existing prepared source, executable, default exchange path or process-global active context. Configuration is fixed at construction; FFT plans are destroyed explicitly, never copied. This is a serial full-support rectangular orthogonal-cell exchange primitive, not a TDDFT solver. Wannier/spatial/symmetry and factor prediction routing remain a later step.

**Ownership audit:** exx_native has shared symmetry_map, kernel, spatial and wannier. exx_k_kernel owns FFTW plans, thread plans and optional accelerator; exx_wannier owns main/worker/local FFT plans. Assignment of these types is unsafe. RT-only data can be copied but resource owners cannot. The new context uses explicit open/close; no copying, serialization or concurrent calls are supported. Keep resources separate from rt_state_clear's value assignment.

**Alternatives:** Moving all shared resources/configuration at once greatly expands regression scope. Swapping a global active context violates the approved ownership direction. First validate a standalone ownership primitive using the same kernel/ACE routines, then connect the native refresh dispatcher explicitly.

1. Add a native probe for two different shapes, spacings, screening/mixing parameters and sources. Compare independent kernel+ACE reference calls; interleave refresh/apply; close in both orders; reject invalid refresh without changing accepted output.
2. Implement an opt-in adapter with fixed configuration, private resource storage, candidate ACE on refresh, copied accepted sources and explicit cleanup. Document errors and the absence of full allocation-failure recovery.
3. Extend the existing audited object-link probe runner to compile this adapter and probe; preserve the existing explicit-routing runner behavior. Audit new sources and all application objects before linking.
4. Run probes at 1/2 OMP, source tests and sdist-to-wheel packaging/installed checks. Include adapter/probe in distribution. Document measured scope in Japanese/English and synchronize AI JSON.

No GS payload collection, upstream changes, general MPI/GPU/Wannier or full multi-material TDDFT claims. New outputs go outside this repository.

## Completed stage verification

- Two contexts differ in grid shape (2x2x2 / 3x2x2), spacings, k count (1 / 2), screening and mixing fraction. Reviewed probe passes at 1/2 OMP and from installed wheel at 2 OMP.
- Occupied-subspace rotations change a fixed target's action while the other context remains bitwise unchanged. Arbitrary-target ACE results match separately constructed reference contexts within 1e-12. Zero/nonfinite/wrong-shape refresh rejects with accepted result unchanged. Duplicate k construction failure permits repeated close and successful reopen.
- A stale-ACE mutation in an external source copy fails at STOP30, proving refresh coverage. Failed binaries/logs are retained only as failure evidence.
- Existing explicit-routing probe passes at 2 OMP. Source suite: 294 tests OK. Installed final wheel: 294 tests OK, 1 skipped. sdist-to-wheel build succeeds; packaged adapter, probe and runner bytes match source. Guides/AI JSON are synchronized.
- Source and all application-object audits precede probe linking. Existing full solver executable still matches recorded SHA256 d2b09f234785118034465e76e4f9d4bd61c7a5a4633c82931370320adac69aa7. No default full-solver path was modified; no new full-solver physics/speed claim.
- Review findings on invariant refresh and partial construction tests were corrected and re-reviewed. Native contexts remain explicitly noncopyable by contract, not enforced by Fortran intrinsic assignment.

Next: connect the independently owned resources and fixed configuration to native exx_refresh/exx_add_action, including optimized dispatch and prediction history; separately extract Wannier/spatial/symmetry ownership. This stage is not full multi-material TDDFT.
