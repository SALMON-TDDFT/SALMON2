# HSE spinor stationarity diagnostic

Run `python3 testsuites/unit_exx_stationarity/run.py --checker /absolute/path/to/tools/check_forbidden_eigensolvers.py`.
The test covers exact eigenstates, fractional occupations, equal-occupation gauge freedom, residual outside the trial span, projected non-Hermiticity, orthogonality defects and NaN/Infinity rejection.

For a native SOC HSE GS input with full real-space/orbital storage and k-only MPI, set `SALMON_HSE_STATIONARITY=1`. After normal initialization the executable refreshes physical exchange, applies H, prints three scalar metrics, and returns **before SCF**. This is a diagnostic, not a converged-state calculation. It does not export orbitals. Existing ordinary executions are unchanged when the environment variable is absent.

For each local k block, S=Psi^H Psi*dV, B=Psi^H H(Psi)*dV, Q=H(Psi)-Psi B. The output reports max|S-I|, max|B-B^H| in Hartree and the maximum over k of

`sqrt((2 sum_i f_i^2 ||Q_i||^2*dV + sum_ij (f_i-f_j)^2 |B_ij|^2) / sum_i f_i^2)`.

The commutator interpretation assumes an orthonormal trial span and a Hermitian operator; the separately reported defects must be checked first. Projected Hermiticity is not a proof for vectors outside the tested span. No acceptance threshold is imposed and no state is certified by this diagnostic. Initial mixed/restart-density input must be identified in the provenance; use an ordinary reset-wfn seed to reconstruct its density for a propagation-start check.
