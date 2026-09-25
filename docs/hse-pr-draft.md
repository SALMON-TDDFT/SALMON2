# Add native CPU HSE with Taylor4 + ACE and configurable short-range screening

SALMON gains optional self-consistent HSE ground-state and real-time calculations
using distributed screened Fock exchange and ACE. HSE RT defaults to Taylor4
with predictor/corrector; standard inputs need no propagator setting. The new
hse_omega input follows the selected inverse-length units and is shared by
native Fock exchange, the Libxc semilocal remainder and restart validation.

The implementation includes k-point MPI distribution, symmetry reconstruction,
OpenMP packing and BLAS actions, bounded automatic FFT layout selection, and
full-k linear Acos2 pulses with optional impulse probes. Unsupported configurations
are rejected; pulse RT restart remains intentionally unavailable. Developer CN
and full-action Taylor paths remain available only by explicit selection.

Ordinary CPU CMake builds enable HSE and discover compatible libraries or build
hash-pinned Libxc, FFTW and BLAS/LAPACK automatically. Existing vendor toolchains
retain their compiler and numerical-library settings.

Validation on Apple Silicon/GNU Fortran 15: HSE-enabled MPI and HSE-disabled serial
builds; eight native numerical tests; thirteen bounded input/restart cases; short
laser+probe and two negative restart checks; seven CTest checks for a converged
Si Nk=4^3 GS and 64-step impulse response with MPI4 × OpenMP2. Automatically built dependencies also passed the Nk=4^3
GS/impulse tests with MPI4 × OpenMP1. An eigenvector regression guards the
selected LAPACK; a GNU 15/AArch64 workaround is confined to the Netlib fallback.
See docs/merge-preparation.md for reproduction, exact scope and remaining CI/manual
work. This candidate is for review, not a claim of cross-platform certification.
