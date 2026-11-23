# Native SCF disk restart implementation

Goal: persist same-topology SCF numerical state in Fortran and resume in a fresh host process without replaying accepted iterations. Keep converged-GS reusable archives separate.

1. Add a developer regression that requires disk ABI and fails on current runtime.
2. Generate checked rank-local stream IO for all numerical snapshot fields; exclude pointers, FFT/backend owners and MPI communicator IDs. Version and authenticate through Python manifest; preserve input/assets and actual decomposition.
3. Add native zero-iteration bootstrap and staged disk read; collect MPI status before publication. Rebind live resources and exchange communicators.
4. Bind explicit save_scf/load_scf entry points without implicit builds. Validate missing/corrupt/incompatible input before native initialization where possible.
5. Compare independent continuous N and fresh-process K-save/load+(N-K), including ordinary XC and fast HSE. Record only requested GS data/results outside checkout.
6. Update paired guides, AI task JSON and checklist; verify sdist/wheel and no CI dispatch.
