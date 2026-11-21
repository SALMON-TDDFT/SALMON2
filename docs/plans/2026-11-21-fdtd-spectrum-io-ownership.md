# E/H spectral output file-unit ownership

### Spectral output file-unit ownership

E/H finalization owns local newunit handles for linear-response trajectories and
spectra, frequency-resolved E/H/J planes, and ASE/ART spectra. The two simultaneous
windowed/unwindowed ASE/ART files use two local handles, not base-unit arithmetic.
The periodic linear-response response file is explicitly closed after its last
record; its earlier implicit close at process termination is no longer needed.
No module-wide file selector is added. Fourier transforms, normalizations, field
updates, loop/file ordering, names, record formats and MPI operations are retained.

The existing native I/O tool has --fdtd-spectrum-io, using actual small LD+FDTD
lifecycles with no electronic GS. GNU MPI1/2 controls cover isolated and periodic
impulse response, isolated non-unit-background raw response, Fourier planes, and
separate ASE and ART cases with three frequencies and both window variants. Caller-owned
units600/601/602 remain connected to their original files. Named spectral outputs
are closed before return; expected numeric outputs are nonempty and finite.
Unoccupied, occupied and previous-runtime outputs are compared byte for byte.
These are output/host-contract regressions in atomic units, not convergence or
material-property certification. Native E/H's existing supported-mode limits are
not lifted: the batch host used here enables spectral modes, not the embedded
step interface. No extra runtime collective, automatic build or startup run is
introduced. TTM file-unit ownership and general IO-error recovery remain separate.
