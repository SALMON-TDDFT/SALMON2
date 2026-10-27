# Spin-polarized Libxc NLCC support

Frozen unpolarized core density is added equally to the two local spin
channels. Density gradients include the same half-core contribution in
spin-polarized calculations, or the full core in scalar calculations.
The gradient preparation is enabled for both Libxc and native HSE builds.
A missing core-density argument with NLCC enabled is rejected.

Validation on macOS (2026-10-27):
- Release MPI build with USE_LIBXC=ON, USE_HSE=ON, Libxc 5.2.3.
- Source and all application-object PZHEEV prohibition audits passed.
- Polarized PBE exchange/correlation with a fixed unpolarized core:
  finite-difference density and gradient derivatives agreed with Libxc
  potentials to 3.061e-10 absolute error in a synthetic check.
  This check used the installed Libxc C API and does not by itself validate
  SALMON's complete noncollinear discretization.
- Ni primitive Nr18, SOC, NLCC, full 2x2x2 k mesh: two SCF steps ended
  normally with finite energies after the previous unsupported-path stop.
- Production PBE SCF did not converge within 3000 iterations at mixing
  0.03, or a further 6000 iterations at mixing 0.005. No converged PBE
  magnetic moment or PBE real-time trajectory is claimed.
- Fujitsu compilation and Libxc-only numerical execution were not tested.

Existing validated production executables and running MGGA23 trajectories
are separate from this build. Ground-state data and logs are not committed.
