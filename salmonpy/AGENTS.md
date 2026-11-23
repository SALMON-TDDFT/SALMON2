# SALMONpy development

Read AI_USAGE.md before using or extending these interfaces. For the first installed-wheel test, run `python -m salmonpy.quickstart`; no solver build or additional input is needed. Crystalline TDDFT and SBE interfaces are A-driven only. The LD Python preview is available; FDTD has a batch executable adapter; the built-in velocity-gauge SBE adapter and text loader are available; a process-isolated PZ/TB-mBJ/Libxc-PZ/Libxc-PBE middlepoint TDDFT adapter is available with a separately audited worker executable. Read WORKFLOWS.md for supported scope.

Start with `python -m salmonpy.manual` for the first calculation or `python -m salmonpy.manual QUERY` to find tasks, APIs, examples and constraints. Consult the registered SALMON-inputs repository for published input examples, checking version and required files.

`manual.json` is exclusively the AI-facing task catalog. Human-facing documentation is separate; do not present JSON search results as the user manual.

For beginner workflows, apply manual.json beginner_intake before choosing calculation scale: ask only missing memory/CPU budgets and intent, then state measured/estimated/unknown capability limits.

This project is intended for a separate repository; do not create a PR or merge
into SALMON itself. Use README.en.md / README.md for English / Japanese onboarding.

This is a standalone SALMONpy layout. Bundled calculation routines live in native/. Complete SALMON sources are bundled in runtime/sources, with only fast-PR HSE supported. Full SALMON is an explicitly built runtime; see runtime/README.md in the source checkout.

Experimental MPI managed by the calling program TDDFT and FDTD updates of E and H coupling now have a separate entry: [COUPLED.en.md](COUPLED.en.md). Legacy coupling has one point/group and one solver initialization per process. SHARED_HSE.md describes explicit same-model full-k Python independent calculation states; explicit points_per_group now connects same-model full/reduced-k shared HSE to normal sessions; do not infer heterogeneous reentrancy or automatic rank assignment.

New numerical work belongs in Fortran/C. Python owns inputs, scheduling, dispatch and records. Timestep diagnostics use measured electron-number conservation, default relative target1e-5; target exceedance is a reported diagnostic, not an automatic stop. Stop only the numerically broken trial.
