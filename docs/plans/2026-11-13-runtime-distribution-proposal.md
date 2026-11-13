# SALMONpy runtime distribution proposal

Status: design proposal, not a new installer or runtime-selection API.

Imports must remain fast, offline and free of compilation or environment changes.
A missing native component is reported when that component is requested. No silent
source-build fallback at import or first calculation. Wheel installation avoids
building on the user's machine when a compatible binary wheel is available:
https://packaging.python.org/en/latest/specifications/binary-distribution-format/

Current implementation: the LD/SBE wheel contains a native library and runtime
dependencies; Mac uses Accelerate for that small kernel. Full FDTD/TDDFT use an
explicit SalmonExecutable(path, sha256). No full electronic-structure runtime is
currently bundled. The existing source helper builds LD/SBE only.

Proposed standard distribution: keep the lightweight Python/LD/SBE distribution;
provide an optional OS/architecture-specific full SALMON runtime package with
validated Libxc and OpenBLAS BLAS/LAPACK and required compiler runtimes. Pin library
versions/ABI (including LP64), record licenses, build options, artifact hashes and
functional support. Choose the supported compatibility floor in CI, not at import.
For GNU Fortran full SALMON on Mac, use OpenBLAS: the prior full-GS test exposed
an Accelerate complex-return ABI failure. Do not borrow NumPy's private BLAS or
arbitrarily locate a system library. An installed Libxc does not automatically
make a functional supported by the worker.

User flow:
1. Install the desired prebuilt runtime in the environment (standard path).
2. Import SALMONpy with no download/build side effects.
3. Select a runtime when preparing a calculation; check requested functional,
   features, binary hash and ABI before starting a costly calculation.
4. If unsupported/missing, show the exact missing capability and explicit setup
   instructions. Never silently change functional, pseudopotential or backend.
5. HPC users may explicitly select their audited system/MPI/GPU build. Preparation
   is an explicit command/job, optionally cached by full source/config/compiler/
   library/architecture identity; switching BLAS is a different build, not an
   in-process knob. Preserve PZHEEVD-only requirements for ScaLAPACK builds.

Model inputs (XC, pseudopotential and TB-mBJ cval policy) are immutable in each
initial-state snapshot. Runtime capability metadata is distinct from GS provenance;
certifying GS/RT consistency still needs a dedicated provenance schema.
The old f90 Libxc interface currently requires the tested compatible Libxc 5.2.3
build; migrating SALMON to a modern interface is separate work, not an import-time
installation responsibility. Libxc is unnecessary for native TB-mBJ.
