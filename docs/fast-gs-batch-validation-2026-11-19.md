# Fast-source GS batch scalar reports — 2026-11-19

The maintained fast-source GS entry now writes the existing SALMONPY_GS_V1 status
and SALMONPY_FORCES_V1 force reports for fixed-ion non-DC DFT. Reporting reads
native convergence/residual/threshold/energy/forces after their existing calculations;
SCF arithmetic and force operation order are unchanged. Imported GS still returns
before computation/output. Fresh serial/no-HSE and MPI/fast-HSE builds passed
workspace source and all-application-object prohibition audits. Fixed reference
executables were reused after recorded SHA256 checks; no PZHEEV was executed.

DFT.run gains per-call save_restart=False, default True preserving existing output.
Energy and force results remain readable without restart arrays. restart_saved
records whether all three native restart files actually exist. Unconverged SCF
warns and retains results. Unexpected extra output from an older runtime warns,
records actual output and does not invalidate a completed calculation. Legacy
manifests default to saved=True. Warm-start/band/TDDFT/archive consumers report
that restart arrays are needed before attempting unavailable data transfer.

The existing GS lifecycle tool adds --batch instead of adding another test module.
Real public serial PZ (unconverged) and Libxc-PBE (large-threshold converged) runs
return energy difference0 against the prior fixed fast CLI. Force comparisons
respect the old report's decimal printing interval plus3e-13au, rather than claiming
exact recovery from rounded text. Old fast CLI omits the Python reports; the
reader rejection was checked. An initial old-public run also exposed the older
-1 default mixing admission, preserved in its log; reference controls explicitly
select the equivalent0 no-mixing-loop setting.

A first numerical comparison used modern Python CODATA eV/Hartree and differed
by1.66966344e-7Ha despite byte-identical native energy text. Legacy gs_info.data
uses inputoutput.f90's fixed27.21138505 conversion. The test now explicitly uses
that convention; no physical constants or calculated values were changed and
no tolerance was widened. Cross-interface constant convention audit remains separate.

HSE CLI reports with matching2-rank inputs parse correctly and old gs_info/eigen
text remains byte-identical. An earlier attempt used one rank for a2-rank input;
its missing-output log is retained. MPI6 HSE different-material retained GS→RT
and direct8→16 current/energy resume remain exact. Existing learned full-k Taylor/
reduced-k PT-CN current/energy/history/fresh/direct-resume controls match the
previous correctly linked fixed build exactly. Public batch HSE selection still
needs explicit fast-PR provenance integration; use NativeGSBackend for HSE.

Source and installed8-atom Si16^3/2^3-k samples ran4 SCF iterations without restart
output. They are intentionally unconverged execution examples, not physical-force
or optical-accuracy validation. No GS wavefunction/occupation datasets were
collected. The source/installed Python suites each run628 tests (627pass,1skip);
new output cases reuse existing test methods. The sdist→wheel→new environment
route, installed batch/sample execution, JP/EN/AI/source byte parity, manual
synchronization and whitespace checks pass. CI/commit/push were not performed.
Local OneDrive storage does not certify cloud synchronization.

Remaining: explicit batch-HSE provenance connection, final maintained solver-tree
consolidation, one-session ordinary/HSE registry unification and larger physical/
platform validation. Default saved-restart execution is covered by existing mocked
contracts; it was not run here because the user prohibits collecting GS datasets.
