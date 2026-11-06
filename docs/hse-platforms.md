# HSE platform status

Use the ordinary `configure.py`/CMake platform selection. CPU HSE is optional;
the generic OpenACC route and the experimental exchange-cuFFT backend have
separate build gates. Do not infer GPU support from a CPU build passing.

The submitted integration was built locally with GNU Fortran15 on Apple Silicon,
serial and MPI, with FFTW, Libxc and OpenBLAS. A non-HSE MPI build was also checked
during upstream integration. Detailed scope is in [validation](PR_VALIDATION.md).

The complete submitted tree has not been numerically certified on Fujitsu/A64FX,
NVHPC/OpenACC/cuFFT, or every optional library/functional combination. Earlier
Fugaku experiments used different immutable source snapshots and do not certify
this merged tree. No remote production executable was changed by PR preparation.

All platforms retain the source/object eigensolver guards and scalar IEEE policy
in [CODING_RULES.md](../CODING_RULES.md). The current testsuite invocation and
library requirements are in [build instructions](hse-build.md).
