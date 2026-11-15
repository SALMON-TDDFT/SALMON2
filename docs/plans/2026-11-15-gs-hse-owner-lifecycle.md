# GS HSE owner lifecycle

Continue the approved explicit-resource design. Main GS creates a local target HSE owner, binds it after system/grid initialization only when the existing native-owner eligibility predicate accepts the configuration, and unbinds/clears after all output and checkpoint work. Reuse the existing predicate rather than broadening support. Keep unsupported configurations, especially DC/distributed/thermal/localized exchange, on the legacy owner. Do not claim independent whole-GS contexts or reentrant globals.

Validate branch selection explicitly in the serial four-case GS regression; compare the original immutable fast-HSE output and transitions. Test eligible/ineligible binding and independent owners in the native probe, keep DC/MPI comparisons on the legacy path, verify TDDFT contexts and package consistency. New source/build directories, local only, maximum four total execution threads.

Completed locally: explicit-owner branch verified, native eligibility/ownership probe, four serial GS and eight DC/MPI comparisons, eight TDDFT context cases, source and fresh-wheel tests. See `docs/gs-hse-owner-validation-2026-11-15.md` for scope and evidence.
