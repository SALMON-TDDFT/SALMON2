# Fast-PR HSE generation-checked handles

Goal: expose real independent same-model native contexts through a versioned C
entry, retaining the tested numerical kernels and old singleton ABI.

Architecture: first handle initializes the fixed process model using the existing
explicit attached-MPI route. New handles start from that same GS with individual
dt/capacity. A dynamically growing slot registry owns context pointers, accepted
A and failure flags; positive int64 tokens combine slot index and generation.
Closed/stale/mixed-rank tokens are collectively rejected before native kernels.
Reusing a slot never revalidates an old token. Close never detaches host MPI.

Use a synchronous main-thread dispatch adapter around the existing legacy C
operations, saving/restoring only the selected pointer/accepted A/flags; never
copy whole contexts/FFT/MPI/XC handles. No callbacks or concurrent entry. Legacy
and handle initialization are exclusive for a process lifetime. The handle route
selects a context only during a call, leaving the legacy singleton inaccessible
between calls. Each new context owns a unique output basename; expose a bounded
label query for host journaling. Numeric checkpoint encoding remains unchanged.

This is one shared immutable model/communicator/working directory, not arbitrary
heterogeneous/threaded execution. Public Python NativeFastHSEBackend retains its
one-lifetime limit until shared leases and its host directory lifecycle are wired.
Symmetry multi-context initialization and common input/I/O/timers remain future.

1. Test missing symbols against the preceding build, then add header/protocol,
   registry, checked selection and init/new/step/info/history/tx/checkpoint/close.
2. Test old/new stale generations, root closed while child continues, slot reuse,
   invalid and mixed-rank IDs, concurrent transactions and pending candidates,
   direct restore and ordinary trajectory/log equivalence to isolated controls.
   Use internal MPI1/2, independent material groups and fast Taylor/learned PT-CN.
3. Regression old singleton/mixed and deep recovery routes. Build only to new
   external directories; source/installed suites and sdist→wheel/hash check.
4. Update JP/EN/AI, checklist and persistent evidence. No Actions.
