# Spatial HSE DC SCF plan

Continue the approved shared exchange migration. Preserve the unweighted MLWF frame Phi=Psi U for temporal transport, but form exchange factors Q=Psi sqrt(f/2) U so fractional and zero occupations give the correct density matrix. Invalidate native caches for occupation-only changes, weight exchange energy by f/2, and reduce DC core exchange across spatial ranks. Admit Gamma full-support HSE DC DFT on y/z fragment pencils; retain other legacy modes.

1. Extend distributed exchange oracle with fractional/zero occupations and an occupation-only update; verify against serial Wannier at multiple omega and ranks.
2. Implement optional occupation in hse_spatial, cache/energy changes in hse_native, DC dispatch and fragment-grid input checks.
3. Add finite-temperature DC serial-fragment vs distributed-fragment energy, core exchange and electron-count comparisons, including reconstruction from resulting LCFO output.
4. Build ON/OFF, run HSE regressions, independent review, document verified scope and remaining migrations.

## Validation and findings

- Distributed exchange oracle passes 1/2/4 ranks, omega 0/.11/.3, full/fractional/changed/all-zero occupations. Unweighted temporal transport remains consistent.
- H4 DC at 10000 K, six states/fragment: SCF density residual <1e-10 on 1/2/4 spatial ranks per fragment (373/1269/101 iterations). Total printed energy -78.972108 eV agrees; core EXX max difference 1.852e-11 Ha; charge error <7e-14.
- Initial SCF energy parity guards the hybrid Gaussian seed correction. Fragment eigen outputs confirm nontrivial partial occupations.
- DC-LCFO files from all three layouts reconstruct and propagate mesh RT, with energy-history agreement <1e-7 Ha.
- HSE suite six tests plus legacy multi-k thermal-charge test pass; existing PBEh Ehrenfest 16 tests and thermal solver unit test pass. ON/OFF builds and diff checks pass.
- Independent review found no blockers in occupation algebra, cache, core reduction, thermal solver or Gaussian seed changes.
- Debugging found spatial Gaussian centers previously depended on local grid origin; spatial hybrids now use the serial seed. HSE positive-temperature DC shares PBEh's charge-converged Fermi solver and capacity checks.
- 30000 K stress runs did not consistently converge by 500 iterations in legacy or spatial routes; changing mixing or orbital iterations did not establish robust convergence. The 10000 K verified case does not certify general high-temperature convergence. No performance claim.
- Legacy multi-k, finite-support and projected exchange routes are retained. HSE MD force validation remains separate.
