# SALMONpy development

Read AI_USAGE.md before using or extending these interfaces. First use: install a
matching wheel and run `python -m salmonpy.quickstart`; no separate solver is needed.
LD/SBE use the bundled calculation library. DFT/TDDFT/FDTD/bands require the
corresponding additional executable or shared library. Current supported
conditions belong to the topic guides, not dated implementation milestones.
Crystal TDDFT and SBE are driven by A/c; molecular length-gauge TDDFT is E-driven.

Start with `python -m salmonpy.manual` for the first calculation or `python -m salmonpy.manual QUERY` to find tasks, APIs, examples and constraints. Consult the registered SALMON-inputs repository for published input examples, checking version and required files.

`manual.json` is exclusively the AI-facing task catalog. Human-facing documentation is separate; do not present JSON search results as the user manual.

For beginner workflows, apply manual.json beginner_intake before choosing calculation scale: ask only missing memory/CPU budgets and intent, then state measured/estimated/unknown capability limits.

This is the independent SALMONpy repository; do not create a PR or merge
into SALMON itself. Use README.en.md / README.md for English / Japanese onboarding.

This is a standalone SALMONpy layout. Bundled calculation routines live in native/. Complete SALMON sources are bundled in runtime/sources, with only fast-PR HSE supported. Full SALMON is an explicitly built runtime; see runtime/README.md in the source checkout.

For MPI-owned-by-caller coupling use [COUPLED.en.md](COUPLED.en.md).
Independent same-model states and explicit different-input model registration
are supported under the conditions in SHARED_TDDFT.md, TDDFT_REGISTRY.md and
SHARED_HSE.md. Sequential main-thread operation is not concurrent threaded execution.
Do not infer unrestricted initialization, restart or physical-model support.

New numerical work belongs in Fortran/C. Python owns inputs, scheduling, dispatch and records. Timestep diagnostics use measured electron-number conservation, default relative target1e-5; target exceedance is a reported diagnostic, not an automatic stop. Stop only the numerically broken trial.
