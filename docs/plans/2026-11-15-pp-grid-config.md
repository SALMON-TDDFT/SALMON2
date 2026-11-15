# Owned scalar pseudopotential grid setup

Continue the approved shared-input extraction. Keep the existing arithmetic/order and native TDDFT common-material guard.

- Capture periodicity, Poisson/FFT selection, scalar/SOC and +U selection and diagnostic setting in a per-grid value config at construction. Explicit config bypasses legacy input. Later updates retain this config.
- Obtain atom count/species from the system, species count from loaded PP arrays. Cache the fixed mapping; reject incompatible topology/config updates before changing grid arrays. Geometry updates remain supported by the existing sequence.
- Remove global atom/species dependencies from cache, interpolation and projector communication-table helpers. NLCC uses loaded PP periodicity and system atom mapping.
- The explicit independent path is scalar, without +U; existing SOC/+U legacy setup remains supported but has unextracted shared state and is not claimed independent. MPI communicator ownership and global timers remain outside scope.
- Verify a failing poisoned-global cache probe first, full scalar isolated grid A/B/A and moved-geometry updates, NLCC mapping, native periodic HSE driver/state regression, Python and installed wheel checks. All runs local and resource-bounded.

Loaded radial numerical tables and grid geometry are immutable for a grid lifetime. Guards compare configuration/layout, not full table contents; no per-step full-table copy/hash is required.
