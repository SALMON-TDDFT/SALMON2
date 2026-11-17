# Transactional native TDDFT step recovery

Approved continuation of native solver failure handling. Standard native periodic,
unpolarized, fixed-ion middlepoint PZ/TB-mBJ/Libxc-PZ/PBE only initially. Before
an accepted step mutates, retain the existing numerical checkpoint payload in
memory plus XC payload/cval. Add per-RT enabled/error fields (default disabled
for normal batch callers). Finite checks at wavefunction, density, Hartree, XC
and endpoint stages agree across the attached communicator before any rank
returns. Returned deep errors restore the original payload and PP phases without
propagation replay; Python raises NumericalError with a retryable state. Reject
unsupported protected profiles before mutation, never silently weaken recovery.

Explicit deep configuration errors use status instead of STOP on protected
paths. Maintain default batch behavior and valid arithmetic order. All checks
precede irreversible RT outputs. Existing external checkpoint encoding stays
unchanged. Memory overhead is a per-step numerical copy; expose uncertainty in
advice, retain warning-only budgets. OOM, MPI/process death and uncontrolled
third-party STOP are outside the handled contract. Whole FDTD-session step
rollback is separate; this guarantee is the individual TDDFT group/state.

Test actual overflow and compile-gated one-rank faults at multiple deep stages,
serialized state equality and output preservation, successful direct retry vs
old/new continuous results, independent groups. Production build has no injected
fault entry. Test PZ/TB-mBJ/Libxc and internal MPI2, ordinary mixed regression,
source/installed suites and sdist→wheel. Jobs<=18 cores, advisory48GiB. No CI.
