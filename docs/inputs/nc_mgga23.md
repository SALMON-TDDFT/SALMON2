# nc_mgga23

Use `xc='nc_mgga23'` in `&functional`; gamma=0.8 is fixed and recorded in
variables.log. This is the corrected 2023 noncollinear BR-CS gKS functional.

Requires complex noncollinear/SOC representation, conventional fixed-ion CPU
calculation and no symmetry reduction. Nonzero NLCC requires a V2 KY file
with consistent auxiliary-Pauli core density, kinetic density and spin current;
legacy density-only NLCC is rejected. Supported RT path:
`propagator='middlepoint',yn_predictor_corrector='y',yn_fix_func='n'`.
Full local, kinetic-density and current derivatives are included. Unsupported
combinations fail instead of falling back to LSDA.

See `samples/nc_mgga23/README.md` for paired input fragments, new diagnostics,
comparison command and validation limitations. The SALMON-DOCS external manual
repository was not modified in this task; this local input reference records
the new selector alongside the implementation.
