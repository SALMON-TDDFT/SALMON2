# CheFSI finite checks without array conversion temporaries

Run `python3 developer_tests/652_dc_lcfo/finite/run.py` (GNU Fortran by default).
FC and FFLAGS override the compiler and flags. The default flags deliberately
put array temporaries on the stack to detect this regression.

The test extracts the actual production Hamiltonian finite-check expression,
generic interface and helpers. A 400 x 400 x 1 x 26 complex array is checked
under an 8 MiB stack limit. The original REAL/AIMAG array-argument implementation
fails this test; the complex scalar checks pass. It also checks rank 2/3 helpers,
NaN in the real component, infinity in the imaginary component, empty arrays,
and a noncontiguous rank-4 section. No MPI or SALMON binary is needed.

This is a host-side compiler/unit regression, not a complete DC-LCFO run.
Fugaku full-module compilation used mpifrtpx tcsds-1.2.43, -Kopenmp -Nfjomplib
-Cpp -SCALAPACK -SSL2BLAMP -Kfast -Kocl -Ncheck_std=03s -Nalloc_assign.
The fixed GS binary and active/queued jobs were not modified.
Compute-node execution, complete link and LCFO output integrity remain unverified.
