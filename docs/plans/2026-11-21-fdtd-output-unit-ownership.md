# Ordinary native E/H output unit ownership

A host-owned unit600 survives rejected shape reads. Ordinary point headers,
point appends, final observation information, plane snapshots/integrals and the
material-ID CUBE writer now open a local newunit and close only that owned unit.
Filenames, record formats, ordering and numerical/communication operations are
unchanged. The generated-shape writer no longer closes the caller's requested
unit for an unimplemented output format. Existing helper signatures are retained.

GNU MPI1/2 tests keep a host-owned unit600 open throughout initialization,
32 steps, direct checkpoint restore and finalization. Field arrays, point files,
checkpoint bytes and generated shapes match prior controls; enabled plane output
also matches. This is file-unit ownership isolation, not output-error recovery:
write failures may still use existing Fortran error handling. Legacy restart,
linear-response/spectral/ASE/ART and TTM output units remain separate work.

Scope: no checkpoint schema or public API changes. Test flags in the existing E/H recovery tool: --host-unit-library DIR --host-unit-lifetime; --observation-planes enables legacy spatial/plane-integral outputs. Keep physical calculations out of required CI.
