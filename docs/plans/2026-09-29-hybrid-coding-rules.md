# Hybrid coding-rules compliance

Approved scope: repair the reviewed coding-rule issues, then centralize functional classification. Preserve equations, data format and ongoing calculation binaries. No unrelated legacy refactoring or exchange-state redesign.

1. Add explicit-name global-hybrid/hybrid predicates with tests; migrate the repeated gates and narrow imports.
2. Add explicit implicit-none to affected hybrid procedures, wrap changed long lines and use yn_argument_check.
3. Register small numbered PBE0 DC-GS and native RT tests, with explicit producer-consumer fixtures.
4. Build separately; run unit and numbered MPI regression tests, check a non-HSE build, document untested optional configurations.
5. Prepare corresponding SALMON-DOCS input documentation. Record locations and verification.

Completed: classification/style fixes, two numbered tests, MPI and non-HSE serial builds, 3 unit and 20 integration tests, 6 CTest stages. Serial numerical verification passed with python3 after legacy shebang could not find python. HSE/non-MPI configuration checked. SALMON-DOCS manual update committed separately as be6aa8e; added RST section parses without warnings. Full Sphinx site build and accelerator configurations untested. See docs/hybrid-coding-rules-ja.md.
