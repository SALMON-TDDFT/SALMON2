# Retain transposed rVV10 spectra implementation plan

> Use superpowers:executing-plans inline, with one final independent review.

Goal: reduce FFTW redistribution cost without changing the discrete functional.
User continuation authorizes the previously identified transpose/data-layout optimization.

Design: keep FFTW forward spectra in Z pencils (memory order z,x,y), apply the
kernel with the corresponding global wavevectors, then inverse in Z,Y,X order.
This reduces each pair from eight to four redistributions. Preserve the default
X-to-X transform API and FFTE reference. Fine-grained blocked packing is a later
alternative; replacing communication libraries is unnecessary for this stage.
No new input switch or changed default. Coordinate-ordered communicators and
existing fallback constraints remain required.

1. Extend distributed_probe.f90 with optional transposed-spectrum calls. Compare
   each Fourier coefficient with the assembled X-layout FFTE reference; test
   normalized inverse and plan reuse for 2/3/4/5 channels and all existing layouts.
   Observe RED before adding the optional API.
2. In fftw_pencils.f90 add optional spectral_z. Forward skips steps3/4; inverse
   executes axes3,2,1 with steps3/4 between axes. Plans and buffers are unchanged.
3. In rvv10_distributed.f90 use Z spectra for FFTW; map flattened z,x,y coordinates
   to physical x,y,z wavevectors. Preserve native real-space output and FFTE path.
   Run the 54 potential cases and native25 tests plus LCFO y/z fixtures.
4. Extend fftw_benchmark.f90 to time old X spectra and new Z spectra in the same
   process with cached plans. Repeat grids32x24x16 and64x48x32 three times, serially,
   recording setup, phases and complete functional results. Retain FFTE default.
5. Independent review, document evidence and limits, commit locally without push.
