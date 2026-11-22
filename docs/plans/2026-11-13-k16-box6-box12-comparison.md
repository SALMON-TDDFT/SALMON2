# Si Nk=Nq=16^3 box comparison

User-authorized extension: Nr=12^3, Nk=Nq=16^3, cutoff=3 Ha, boxes 6 and 12.
Hold TB-mBJ, eta=0.001 Ha, timestep=0.05 au, impulse=1e-4 au,
predictor/corrector, orthonormalization interval=1 and the 16384-step target
fixed. GS has eight states; real-space RT propagates four occupied states.

Reuse the existing converged Fugaku q16 TB-mBJ GS as the read-only seed.
Each case refines that seed with identical 1024-rank SCF settings to 1e-11,
then computes W with the unchanged validated Fortran executable. RT reads
its own generating GS so TB-mBJ c/Vxc match the W header exactly. Require
identical reference metadata in both W headers, excluding the different
box-count field, before accepting the comparison. No GS data are collected.

Both cases use 256 nodes, 1024 MPI, 12 OpenMP threads/rank, with 256 MiB
OpenMP stacks. Box12 W is 91.125 GiB globally, 91.125 MiB/rank. Each local
k-owned complex matrix or FFT buffer is 182.25 MiB/rank; these are individual
components, not a full memory bound. W generation and 16-step RT pilots
precede full RT. Record per-rank RSS when /usr/bin/time is available.

Submit four sequential PJM subjobs, each with a 24-hour limit:
1. Both W preparations, reference-metadata equality, both RT pilots.
2. Box6 kicked and zero-field RT.
3. Box12 kicked RT.
4. Box12 zero-field RT.
Each RT stops gracefully after 20 hours if necessary and writes one
restart checkpoint; periodic full-W checkpoints are disabled. A box12
active-W plus reference checkpoint is approximately 273.375 GiB, so do not
write it every few steps. Partial shutdown is labelled explicitly and never
counts as a complete 16384-step spectrum. W/pilot failure blocks later RT
through explicit completion guards. Preserve existing jobs and binaries.

Compare identical time windows with individual zero-field subtraction,
SALMON cubic damping, 1–6 eV Re/Im spectra, stable first peaks and adjacent
structures, peak heights and common direct-gap offsets. Box12 eliminates
box projection, but not finite-grid/G-cutoff/q-sampling approximations.

At preparation, the user reports Fugaku switching systems and requests a
roughly ten-minute pause. Wait before further connection/submission attempts;
report submission only after actual PJM acceptance and state verification.
