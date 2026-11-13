# ASE structure to reusable ground state

Approved continuation: make the initial Python workflow complete by accepting
ASE Atoms, explicit pseudopotentials/electron count/grid/k mesh/functional and
returning a converged GS usable for bands and TDDFT.

Architecture: optional ASE adapter constructs immutable, ASE-independent DFT model;
model owns pseudopotential bytes, geometry and Hamiltonian settings. SCF stopping
settings and execution controls belong to each run. Native SCF emits a versioned
status with its existing convergence flag, iteration count, residual and energy;
no numerical update changes. Successful result manifests hash native files and
permit relocation/reopening. Separate initial-state conversions reuse existing
snapshot loaders and never mutate saved GS.

Scope: fixed-ion periodic unpolarized PZ/Libxc-PZ/PBE, scalar KY or NC UPF input,
explicit electron count, CPU one-rank launch. No Calculator, force/stress, geometry
optimization, pseudopotential downloading or automatic physical convergence study.

Tests: immutable model/pseudopotential snapshot; unsupported ASE charge/spin/PBC;
explicit native nonconvergence; archive tamper rejection; generated input units;
real Si GS versus direct native, band and independent TDDFT consumers; wheel import
without ASE plus optional ASE tests. Update AI JSON and separate human topic guide.
