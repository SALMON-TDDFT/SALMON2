# Si box and G-cutoff diagnostics

Requested conditions: Nr=12^3, Nk=Nq=4^3, TB-mBJ screening reference held fixed.
Six unique cases: (box, cutoff Ha)=(3,3),(4,3),(6,3),(12,3),(6,1.5),(6,6).
Use a common TB-mBJ GS (four occupied plus four empty states), four occupied
RT orbitals, dt=0.05 au, nt=16384, midpoint predictor/corrector and individual
zero-field controls. Run on Fugaku using the unchanged validated response-fix
executable, 64 MPI x 12 OMP (16 nodes). Verify its recorded SHA only at startup.

The small-case stage makes the common GS and TB-mBJ-only controls, then runs
five smaller cases. A following PJM step runs box=12 using the same GS, with a
three-step response pilot before its two full RT runs. Wall limits are 8 h
for the small-case stage and 24 h for box12; these are limits, not measured runtimes. Log rank peak memory
where /usr/bin/time is available. Keep failed outputs. Do not collect GS data.

Code findings: box=12 removes the box-constant orbital projection and box
averages; finite grid, G cutoff, Gamma closure and q sampling remain. G terms
above the cutoff are omitted, not restored as a bare tail. Pair-owned W is
not the only workspace: k-owned sigma/pk/action matrices and FFT transpose
buffers remain dense. Box12 W is 1.424 GiB globally; with 64 ranks, W itself
is 22.78 MiB/rank and each local complex k matrix is 45.56 MiB. These are
component sizes, not a complete memory bound. Set an adequate OpenMP stack
and measure the pilot rather than assuming total memory equals W/ranks.

Preflight: all vertices for cutoff 1.5, 3 and 6 Ha on this mesh satisfy the
strict 12-point-grid Nyquist check (maximum abs index 2,3,4 respectively).
Compare identical 19.8155 fs cubic windows, Re/Im spectra from 1 to 6 eV,
first stable strong peaks, heights, adjacent structures, common direct-gap
offsets, and 75-percent-window stability. Treat Re epsilon(1 eV) as a finite
frequency response. Small LDA/TB W differences alone do not establish which
spatial range dominates, and a box trend is not automatically a unique
short-range diagnosis. No operator optimization or functional change is
part of this scan.
