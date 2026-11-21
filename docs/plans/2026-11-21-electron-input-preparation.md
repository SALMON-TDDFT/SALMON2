# Preparation-only DFT and TDDFT inputs

Extend the existing preparation module; no new implementation module or run-capable prepared object. DFT prepares the existing immutable model's native input and PP/symmetry bytes. Share the restart-output text helper between preparation and run; retain the existing DFT.input_text signature so custom renderer overrides remain compatible. Existing warm starts, HSE fast-PR provenance selection and numerical execution remain unchanged; initial-guess staging is outside the new fresh-DFT helper.

TDDFT prepares an explicitly supplied native RT input byte-for-byte and explicitly listed local files. Do not infer/reset dt, pulse, gauge, parallel layout or restart paths. Do not scan/copy/load GS automatically or require a binary/GS presence on the management node. This is not a replacement for TDDFTInitialState or a universal physics validator. Use streamed file copies/hashes for explicit assets.

Both entrypoints return a new Path and a prepared-not-run manifest; no subprocess, MPI initialization, compilation, worker, submission or automatic continuation. Preserve overwrite protection/staged publication and input bytes.

TDD: input-only contracts under MPI/scheduler environment and forbidden subprocess/worker callbacks; exact source/native-run inputs and explicit TDDFT assets; side effects absent. Numerical comparison uses unchanged audited executable, DFT save_restart=False (no new GS dataset) and TDDFT referencing existing GS without copying it or writing restart output. Verify no new GS records and compare native result bytes. Repeat sdist/wheel tests and typical examples. Update JP/EN/API/AI guides and checklist. No manual heavy CI.

## Verified completion

- Prior FDTD preparation committed/pushed as549aa120; its automatic lightweight CI44856367263 passed. No manual heavy CI dispatch.
- Shared prepared-file staging streams selected source files and hashes, preserves the previous FDTD input/manifest bytes and rejects overwrite. Existing DFT.input_text signature and custom overrides remain intact; run and prepare share only the unchanged restart-output text transform.
- Management-node tests forbid subprocess.run/Popen and GS snapshot loader calls under MPI/scheduler environment. Both prepared Paths stop at prepared-not-run. The installed prepare_electrons dft/tddft CLI behaves the same and creates no stdout.log or GS directory.
- Source and sdist-built installed-wheel suites each run671 tests, one optional skip and zero failures. Two new preparation contract methods, no new implementation modules; one typical DFT/TDDFT preparation example.
- Unchanged audited executable SHA256 e6f75bd582a034bcaf9246f996477ca8f3184fa55d52ef795435822319427b61 reused after hash verification; no new full-runtime build or repeated startup source/object audit.
- PZ DFT: short eight-iteration SCF with save_restart=False, including the preceding installed wheel control, has identical gs.inp, gs_status.dat, gs_forces.dat, gs_eigen.data, gs_k.data and gs_info.data across all paths. It is intentionally unconverged (warning and retained last iterate), not a material/convergence certificate. No new GS state retained.
- PZ TDDFT: supplied impulse input dt.003/nt8 references an existing GS in place, without snapshots or copies. Nonzero matter current and energy output files are byte-identical across raw/native and prepared source/wheel routes. All existing GS file hashes remain unchanged; no wfn.bin, occupation.bin or rho_inout.bin collected in the new result directory. Worker propagation and sign interfaces were not modified.
- Paired preparation/API guides, DFT link, AI JSON, examples and checklist updated. Canonical evidence/distributions: calc2026/salmonpy-electron-preparation-20261121/record.json. OneDrive cloud synchronization unconfirmed.

Fresh DFT input staging and explicit native RT-file staging are complete. General initial-guess staging, native runtime/GS compatibility validation on a real site and multi-node deployment remain separate work. HSE still requires the fast-PR execution profile. Changes remain uncommitted.
