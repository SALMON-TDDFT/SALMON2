# Owned GS convergence metric validation

## Change and boundary

Added `convergence_metric` to the internal `s_gs_run_config` in the fast-HSE
explicit-routing source overlay. `gs_input_config` captures it once; `gs_run`
checks the existing five names before activation and uses the selected value for
pre-SCF semantic restrictions. The SCF routine uses an owned local metric for
residual calculation, normalization selection and labels. Actual numerical sums,
normalization formulas and collective reduction order are unchanged. Legacy
callers without owned GS input still capture the old global selection.

This does not isolate mixing, normalization's electron-count inputs, grid/k/XC
initialization, communicator or timer resources. It is not evidence of general
different-material or concurrent GS contexts. Molecular work remains paused.

## Tests and debugging evidence

- Extended native probe initially failed to compile against the old context type:
  `/tmp/salmonpy-gs-metric-red-20261116/build.log`.
- Captured settings are used after poisoning the shared convergence name and the
  four earlier controls. Unknown metrics and potential-metric pre-SCF settings
  are refused before activation. Mutating caller configuration does not change
  retained A while B runs/closes. Closing twice remains safe.
- All five metric modes were compared with the original fast-HSE full executable:
  rho_dne, norm_rho, norm_rho_dng (four pre-SCF cases each); norm_pot,
  norm_pot_dng (one non-pre-SCF case each). All 14 cases match per-iteration
  convergence labels/residuals, scalar outputs and all generated restart bytes.
  Large thresholds make these control-flow regressions, not physical SCF accuracy.
- Independent review found no required Fortran changes; its request for direct
  per-iteration residual comparisons was implemented to avoid relying solely on
  early-converged final states.
- First build accidentally omitted Libxc; comparison differed by ~3e-14 eV.
  Inspection identified the configuration mismatch. A subsequent generic Libxc
  lookup attempted its build fallback and was stopped. Final build pins the same
  existing `/tmp/salmonpy-libxc-install` and requires package discovery; no fallback
  is used. All final comparisons pass without relaxing equality.

Final source/build: `/tmp/salmonpy-gs-metric-source-20261116`,
`/tmp/salmonpy-gs-metric-build-20261116.json`. GNU/OpenBLAS/Libxc, serial CPU HSE,
-O2; no free-line-length override. Source and application object prohibition
checks pass. At most two build jobs. No Actions, publication or new XC claim.

Cases: `/tmp/salmonpy-gs-metric-pinned-{metric}-20261116/result.json`.
Each stores input/probe/build/reference provenance, checkpoint hashes and the
reference residual trace. Original baseline commit c9e8ce431f43eefd4a88998134adb4796f8ec837.
Source Python suite: 435 tests, one skip, no failures
(`/tmp/salmonpy-gs-metric-python-20261116.log`).

Follow-up regression: eight native HSE TDDFT interleaving cases (2/3 k points,
existing propagators/factor modes, distinct dt/capacity/OpenMP choices) match
original separate execution: `/tmp/salmonpy-gs-metric-tddft-20261116/result.json`.
This is otherwise common-input synthetic TDDFT, not arbitrary-material isolation.

sdist-to-wheel build and installed tests pass (435 tests, two skips). Packaged
three changed native sources and both inventory records match checkout bytes.
Evidence: `/tmp/salmonpy-gs-metric-dist-20261116`,
`/tmp/salmonpy-gs-metric-installed-tests-20261116.log`.
Dependencies were reused offline. No fresh MPI/DC, HSE-OFF or other-platform
execution was done for this increment; earlier validation is not re-labelled as
new evidence. JP/EN guides and AI JSON are synchronized.
