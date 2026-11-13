# TDDFT context extraction checkpoint

The distributed response service now connects the existing process-isolated
TDDFT worker. This is deliberately not a claim of embedded multi-rank TDDFT.
A material group currently has one Python MPI rank and several resident serial
workers (one per point), called in sequence. Provider and response implementations
remain independently replaceable through TDDFTPreparedState plus a verified
managed serial runtime. No native numerical routine changed in this checkpoint.

Inspection of the prepared runtime identifies these next extraction boundaries:

| Current source / owner | Proposed explicit ownership | Equivalence requirement |
|---|---|---|
| src/rt/main_tddft.f90 local system/rt/energy/md and iteration counters | One electronic-state context per macro point | Preserve endpoint A assignment and odd/even input/output orbital swapping |
| spsi_in/spsi_out/tpsi, rho/rho_s, Vh and Vh_stock1/2, Vxc/V_local | Independent mutable arrays in each point context | Alternating point calls must equal separate uninterrupted runs |
| info, lg/mg, stencil, srg/srg_scalar, reciprocal grid and Poisson objects | Initially owned by the context; share only proven immutable data later | Subcommunicators must descend from the supplied solver group |
| xc_func, pp/ppg/ppn, ewald and nonlocal-k updates | Explicit lifetime tied to point/initial-data ownership | No release of another point's XC handles or mutable k-dependent projector state |
| salmon_global and input/configuration modules | Immutable validated per-material/run config passed explicitly | Switching material/XC/input cannot change another point's parameters |
| ofl, fixed unit 030, base_directory and worker stdin/stdout protocol | Context-owned output resources; callable step API separate from stream frontend | No process-wide working-directory switch or fixed-unit collision |
| worker_initial_A module save and shared timer counters | Per-context A/history and timing | Repeated creation and interleaving do not inherit prior excitation/timing state |

The present worker's time evolution already passes many physical arrays explicitly
into time_evolution_step, but its orchestration and transitive calls still use
module state. Do not simply lift the visible local arrays into a derived type
and declare reentrancy achieved. Audit transitive SAVE/module variables, I/O units,
MPI defaults and finalizers while extracting create/advance/close entry points.

Next native milestone: two contexts in one process at one MPI rank, with identical
initial state but different A, alternating steps and independently closing one.
Then run the same contexts with a supplied multi-rank solver communicator. The
existing process-isolated service is a reference for both comparisons, not a
replacement for implementing internal MPI/state ownership.
