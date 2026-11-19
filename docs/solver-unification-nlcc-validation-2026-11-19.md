# NLCC migration prerequisite for solver-source unification — 2026-11-19

The effective fast tree lacked two parts of the standard NLCC correction: Libxc
energy density was weighted by valence density while its potential used total
density, and calc_force lacked the sampled-core derivative. Existing fast total/
half-core gradients, spin handling, model-owned scratch and HSE projector-force
logic were retained rather than replacing whole numerical modules with standard
versions. Before the fix, prescribed-density Libxc-PZ/PBE variation failed by
about0.122au; after the density-weight fix the largest measured derivative error
is3.06e-10 (PZ/Libxc-PZ/PBE, core on/off, OMP1/2).

calc_nlcc_force now reads model-owned species and PP boundary configuration. V1
retains its historical cutoff/interpolation; V2 matches the density's full table
support, uniform lookup and lattice-plane periodic-image bound, including skew
cells and cores extending across multiple cells. All five force callers supply
their current Vxc. Status/force reductions use all solver replicas for rejection
and grid-only groups for the force, avoiding duplicate orbital contributions.
Scalar IEEE inquiries inspect every used table/potential/force element. V2 unused
padding is not inspected; an added negative-padding test first failed and now passes.
Fast HSE MD/optimization scope restrictions and external Calculator gates remain.

Broad V2 force tests initially crossed piecewise-linear radial knots at1e-6au
and large common-energy differences amplified cancellation. The final finite-
difference check uses pointwise density differences before integration and1e-7au
displacement; the original1e-7 force-error limit is unchanged. V1 still checks
1.5x spin averaging; broad V2 also uses exactly representable0.5x scaling to
separate interpolation tests from accumulation roundoff. Maximum final V2 force
variation error is5.36e-8au; V1 is below8e-10au in these controls. These are
prescribed-density implementation checks, not converged material-force evidence.

The original MPI1e-13 threshold caught a1.8205e-13 broad-V2 summation difference.
Only the new V2 core-force integral now uses compensated summation. The same
threshold passes with2.5018e-17 V2 error and6.2450e-17 V1 error across MPI2/4 and
OMP1/2. One bad spatial/orbital replica rejects all ranks; NaN/Inf and overflow
checks remain. Failed pre-compensation/padding probes and logs are retained.

Fresh final HSE-off serial and fast-HSE MPI builds passed workspace source and
all-application-object prohibition checks. No PZHEEV was called. Existing fixed
binaries were not replaced. No-NLCC PZ/PBE batch energies differ by0 from the
prior build, force reports respect native printing precision. No-NLCC HSE GS
metadata and all five array hashes match the prior build exactly for both materials
and local ranks. Six-rank GS→RT current/energy/direct-resume controls, independent
states, learned full-k Taylor/reduced-k PT-CN/history regressions also pass exactly.
No GS restart or numerical-array dataset is collected.

Fast MPI/worker adapter fallback sources now come from the prepared tree rather
than reading runtime/sources/standard at link time. The current worker's fallback
selection was inspected, not separately exercised with a new serial HSE worker.
Preparation still references shared standard files through verified manifests;
legacy source deletion and default-profile migration are deliberately not claimed.
The three standard GS internals overlap fast retained-GS functionality; conventional
standard-HSE modules are not migration targets. Byte differences alone do not prove
physical parity, which is why this missing NLCC fix was addressed before deletion.

Source/installed628 tests each pass627 with1skip; native cases extend existing
files, no new unittest methods. sdist→wheel→new environment, installed HSE MPI6
restart, JP/EN/AI/build-tool byte parity and installed preparation736-file parity
pass. No Actions, commit or push in this slice. This does not grant ASE NLCC-force
capability records, certify material convergence or validate Fujitsu compilation.
Local OneDrive storage is not proof of cloud synchronization.

Remaining: rehome shared source files into one current tree, migrate default
preparation/packaging provenance and verify remaining archive/NLCC runtime capability
routes, then retire duplicate maintained trees while retaining original history,
licenses and fixed comparison evidence.
