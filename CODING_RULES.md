# SALMON2 Coding Rules

This file defines repository-wide coding and verification rules for SALMON2.
These rules apply to new and modified code. Existing legacy code may not fully conform; do not refactor unrelated code solely to bring it into compliance.

## Project Structure

- SALMON is a first-principles simulation code for optical response and electron dynamics based on TDDFT.
- Source code is organized mainly under `src/` by feature area, with tests under `testsuites/`.
- The SALMON user and developer manuals are maintained in the [SALMON-DOCS repository](https://github.com/SALMON-TDDFT/SALMON-DOCS).

## General Development Rules

- Read the relevant existing implementation before changing code.
- Follow this document for matters it specifies. Otherwise, prefer the style, module boundaries, data structures, and helper routines already used near the target code.
- Keep changes narrowly scoped. Avoid unrelated refactoring, formatting churn, and generated-file changes.
- Preserve serial/CPU and all existing optional parallel, accelerator, and third-party-library build paths, including MPI, OpenMP, OpenACC/CUDA, ScaLAPACK, EigenExa, FFTW, and LibXC, unless the change explicitly targets them.
- Do not modify MPI, OpenMP, OpenACC/CUDA, vectorization, or other optimization-specific code as part of an unrelated change. Modify such code only when the change explicitly targets the corresponding parallelization or optimization mechanism.

## Fortran Style

- Use `implicit none` in every program unit.
- Place subroutines and functions in modules whenever practical so callers have explicit interfaces.
- Prefer existing derived types from `src/common/structures.f90` over ad hoc shared state.
- In new and modified code, avoid mutable global or module state except for input variables and other deliberately shared framework state. Pass routine inputs and outputs explicitly through arguments and declare their `intent` whenever practical.
- Keep module variables `private` unless they are intentionally exposed as part of a well-defined module interface. For new modules, prefer default `private` visibility and declare exported procedures, types, constants, and exceptional shared variables explicitly as `public`.
- Put `use` statements in the narrowest practical scope, preferably inside the subroutine or function that needs the dependency. Use module-scope imports only when required by module-level declarations, interfaces, or genuinely shared implementation, and use `only` where practical.
- Follow existing SALMON naming conventions: use `s_*` for shared core derived types, while module-private helper types follow the local convention; use `yn_*` for yes/no character flags and `nproc_*`, `icomm_*`, `id_*`, `isize_*` for parallel decomposition data.
- Use existing precision conventions such as `real(8)` and `complex(8)` unless a surrounding module already uses a different abstraction.
- When adding preprocessor-dependent code, guard it with the existing feature macros from `config.h`.
- Keep new or modified free-form Fortran source lines at 132 characters or fewer, using portable continuation lines where needed.

## Input Variables

- Define user input variables in `src/io/salmon_global.f90` and read them through `src/io/inputoutput.f90`.
- SALMON input uses Fortran namelists. Add a new variable to an appropriate existing namelist, or define a new namelist when no existing group is suitable.
- In `read_input_common`, provide a default value, read the variable with its namelist, and broadcast it to the relevant MPI process group. When adding a new namelist, also add the corresponding `inml_*` status handling.
- For dimensional input variables, support `unit_system` consistently by using the existing `*_to_au` and `*_from_au` conversion conventions for defaults and values read from namelists.
- Add case-insensitive character input variables to the existing `string_lowercase` normalization in `read_input_common`.
- Add each new input variable to the `variables.log` output in `dump_input_common` so its defaulted or read value can be inspected. When adding a new namelist, also output its read status.
- Add value-range, validity, and cross-variable compatibility checks to `check_bad_input` as appropriate.
- For `yn_*` yes/no input variables, use the existing `yn_argument_check` procedure. Use `yyynnn_argument_check` for multi-character sequences of yes/no values where applicable.
- When adding or changing user-facing input variables, update the corresponding manual content in the SALMON-DOCS repository.

## C and CUDA Style

- Keep C source compatible with the repository's C99 build setting and follow the conventions used in the surrounding C or CUDA code.
- Preserve Fortran interoperability through the existing `iso_c_binding`, `bind(C)`, and wrapper interfaces when changing code across language boundaries.
- Keep compiler-, accelerator-, and feature-specific code behind the existing configuration and preprocessor guards.

## Parallelization

- Use the communication abstraction in `src/parallel/communication.f90` when an existing wrapper covers the operation.
- Respect the decomposition stored in `s_parallel_info`, including local ranges such as `im_s:im_e`, `ik_s:ik_e`, and `io_s:io_e`.
- Root-only output or filesystem operations should use existing root checks such as `comm_is_root`.
- When adding MPI or file I/O logic, verify behavior for serial, MPI, and relevant restart/checkpoint layouts.

## I/O and Restart

- Use existing filesystem helpers and file-handle utilities where available.
- Keep checkpoint and restart readers/writers symmetric in file names, local index ranges, and data layout.
- When changing restart formats, preserve backward-compatible behavior unless an intentional format change is documented in the PR.
- Add validation or compatibility checks when a restart file layout depends on parallel decomposition.

## Tests

- Add or update `testsuites/` cases for user-visible behavior, restart/checkpoint changes, and parallel-layout changes.
- Follow the existing test directory pattern: numbered directory, `CMakeLists.txt`, `inputfile`, `preparation`, and `verification`.
- Register new tests in `testsuites/CMakeLists.txt` for the appropriate CPU/GPU sections.
- Keep tests small enough for routine CI while still checking the behavior that changed.
- For restart tests, make the producer and consumer relationship explicit in `preparation`.

## Build and Verification

- Use an out-of-source build directory; do not generate CMake build files in the source tree.
- At minimum, verify that the changed code builds with a standard Fortran compiler.
- For MPI-sensitive changes, also verify an MPI-enabled build and a representative MPI run.
- For optional-library paths, test the affected configuration when feasible and state untested combinations clearly.
- Before submitting a PR, report the build options and tests that were run.
