# SCF exchange-control owner routing

Extend the approved system-bound HSE owner design. Add internal owner-aware queries for adaptive readiness, frozen exchange and pending support changes, plus a setter with optional fields. Route GS SCF reads/writes and localization validation through the system. Preserve old public aliases for unconverted callers and resolve unbound systems to the legacy owner. Do not duplicate flags in s_dft_system: their lifetime belongs to the exchange resources.

Preserve collective support-event reduction, reset timing and temporary freeze/restore ordering. Test independent bound owners and unbound compatibility, then repeat serial GS and two-/four-rank DC comparisons with the immutable original fast-HSE executable. Keep physical convergence, full GS resource construction and arbitrary concurrent calls outside scope. Local tests only; maximum four total execution threads.

Completed: native owner tests, serial GS four cases, DC/MPI eight cases, independent TDDFT eight cases, Python and fresh-wheel checks. Details and limits: `docs/scf-owner-validation-2026-11-15.md`.
