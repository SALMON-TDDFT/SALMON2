# LD regression oracle

legacy_ld.f90 freezes the arithmetic of upstream 3322d056 fdtd_eh (eh_store_old, eh_update_ld, coefficient loop). Only captured inputs and state names were made explicit. Do not regenerate the oracle from modified production code. Full application comparisons provide an independent integration oracle.

SBE legacy procedures are frozen from the same pinned upstream revision; the
regeneration tool no longer reads this repository's Git history. Run
`python tools/test_sbe_native_regression.py` after building the package.

Compile `native/lorentz_drude.f90`, `tests/native/legacy_ld.f90` and
`tests/native/test_ld.f90` in a temporary build directory for the LD oracle test.
`test_pzheevd.f90` is an optional ScaLAPACK check, not a wheel dependency.

The coordinate-initialization regression uses the exact dependency-free
`set_coo_em` routine from the maintained source. Run explicitly with
`python tools/test_fdtd_coordinates.py NEW_EXTERNAL_OUTPUT --compiler gfortran`.
It audits sources/application objects, builds only a tiny probe, and checks
periodic/nonperiodic odd/even coordinates and all unused padding after NaN/Inf
poisoning. No SALMON job, MPI, GS data or network is involved. GNU-tested;
Fujitsu compilation remains unverified. Existing native-contract CI includes
this small synthetic check in its existing job, not a new physical campaign.
