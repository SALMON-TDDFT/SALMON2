# Native k-exchange automatic selection implementation plan

**Goal:** Omit optimization switches in normal inputs; retain one explicit comparison selector.

**Architecture:** Add exx_k_algorithm=auto/reference/pair/orbits in &functional. A stateless resolver derives the three existing switches after broadcast/normalization. Auto preserves conservative spinor/reference and existing eligible scalar-pair policy; spatial orbit reuse and spinor pairs remain manual. Blank internal defaults distinguish omitted legacy switches; explicit old switches retain compatibility, conflicting explicit new/old settings fail. Actual native selection and its reason are logged. Physical controls, field symmetry and all remote binaries/jobs are unchanged.

**Tech stack:** Fortran, CMake/CTest, GNU MPI/OpenMP/HSE, Python GS/RT fixture.

1. Add a failing compiled selector test covering defaults, each manual mode, old inputs, invalid values and conflicts.
2. Implement private stateless selection module; register unconditionally for input parsing without USE_HSE. Add default/read/broadcast/lowercase/log/check handling and contextual manual-mode guards. Derive flags once; record actual selected route at native initialization.
3. Register numbered testsuite445 for selector cases. Extend GS/RT fixture to omit switches for auto and support explicit comparison modes. Validate auto/reference/pair from identical saved GS, both polarizations and symmetry subgroups, plus invalid/conflicting inputs.
4. Update local and SALMON-DOCS manuals and example inputs. Build out of source, run standard tests and existing regressions, run GS and 20-step RT validation. Independently review and address findings. Record untested optional paths; commit results. No Fugaku UI/network/mutation is required.
