# Si wheel-to-current introductory workflow

Goal: supply a reproducible beginner route from an installed binary wheel to a
real SCF-converged Si ground state and A-driven charge current, with explicit
runtime/pseudopotential provenance and numerical acceptance criteria.

Approved scope follows the user's request. No import-time build/download, public
release invention, synthetic GS substitution, or upstream source-tree restore.

1. Correct AI task current signs and guard against stale electronic-current labels.
2. Bundle one pinned upstream Si pseudopotential and small fixed PZ input record,
   preserving Apache license/provenance; describe this as a numerical starter,
   not a spatial/k/time-converged optical result.
3. Supply an explicit preparation/run example: validate the exact successful
   Actions runtime artifact, install via existing runtime APIs, copy inputs into
   a new directory, run DFT then TDDFT sequentially, and compare expected data.
   Retain logs and separate setup, calculation, validation and measurement outcomes.
4. Add a fresh-venv binary-wheel-only test driver, clear English/Japanese guides,
   AI task metadata, package data rules and failure tests.
5. Download the tested CI runtime, build the changed package's distributable
   wheel explicitly, then install in a fresh venv. Execute without repository
   imports or local solver build; measure resources separately. Record exact
   versions, hashes, inputs, resulting currents and limitations.

The initial reference comes from the released branch's tested prebuilt runtime.
Numerical tolerances will be stated, never silently widened on failures. Platform
availability and actual local verification are reported separately. The current
Actions assets expire; future stable release hosting remains a separate action.
