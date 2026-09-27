# Ordinary CMake for Fugaku

Goal: configure the HSE/MLWF/ACE branch on a Fugaku login host with cmake -S . -B build, then cmake --build build, without handwritten compiler/library flags.

Design: before project(), auto-select the canonical Fugaku toolchain only on Linux when both mpifrtpx and mpifccpx are available and no explicit toolchain/compiler (including CC/FC) was supplied. SALMON_PLATFORM=auto|fugaku|generic provides deterministic override. Explicit conflicting selections fail clearly. The canonical toolchain defaults to Release, MPI and vendor ScaLAPACK; explicit USE_* options remain authoritative. Keep the old fujitsu-a64fx-ea name as a compatibility include. HSE dependencies continue installed-library compile/link detection or automatic source builds using the target toolchain. Do not run target programs during configuration.

Alternatives: requiring a toolchain command does not meet ordinary CMake intent; guessing by hostname is brittle. Compiler availability plus conservative override rules gives reproducible selection without overriding user choices.

Validation: add script-mode CMake tests with fake compiler paths (selection only; never claim Fujitsu compilation). Cover auto Linux, non-Linux, missing compiler, explicit compiler/toolchain/environment, generic/fugaku overrides, invalid/conflicting mode, compatibility toolchain defaults and user OFF settings, and propagation into dependency arguments. Observe failure before implementation. Fresh normal local MPI/HSE configure/build and HSE-off configure/build validate existing paths. Review the diff, record the exact Fugaku command and remaining host validation, and push.

Scope: no 3D numerical jobs or remote login in this task. Dependencies need network on first fallback build or compatible installed libraries. This task does not certify the Fujitsu compiler without access to it.

Ledger: missing selection module failed the first test as expected. Selection/default tests passed including fake Linux compiler availability, explicit overrides, missing/invalid/conflicting modes, preserved Debug, repeated include, MPI-off, old alias, and dependency toolchain propagation. Fresh local MPI/HSE and HSE-off builds completed with installed dependencies and no manual dependency paths. Read-only review found no code blocker; it identified the existing case-insensitive CMakeFiles ignore rule hiding the new module, so the module was explicitly added to Git. Suggested additional mode/default tests were added. Fujitsu execution remains unverified.
