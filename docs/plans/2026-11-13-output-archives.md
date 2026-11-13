# Selectable postprocessing outputs

Approved direction: retain native outputs with an index, units and time conventions;
large spatial data are opt-in. Same-host process separation is already tested,
remote-node transport and CPU affinity are not implemented by this task.

Implementation:
1. Add a versioned archive writer/reader for named NumPy arrays and JSON metadata,
   and preserve native files with checksums and table headers. Never overwrite a run.
2. Add TDDFTState.finish(directory): close the worker cleanly to flush native files,
   save accepted A/J histories, original input and initial-file hashes, executable
   identity and generated native outputs; omit unchanged initial-state assets.
   Finish closes that state, preserving factory reuse. Failed/replayed attempts
   are not mixed into the final accepted output directory.
3. Add optional output_directory to common crystal.respond; TDDFT uses finish;
   other A backends save the common response. Keep no-output behavior unchanged.
4. Select native density/current/ELF/electrostatic/decomposed-current output through
   existing immutable &analysis settings and strides; no default volume dumps.
   Explicitly mark unsupported wavefunction checkpoints/projections and distinguish
   recorded arrays from derived observables. Do not invent units or time locations.
5. Test archive roundtrip/no-overwrite, close/finalization behavior, live TDDFT raw
   output preservation and native comparison with spatial flags enabled; document
   availability/limitations and reproducible postprocessing examples.
