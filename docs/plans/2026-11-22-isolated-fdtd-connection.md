# Isolated TDDFT to FDTD implementation plan

Goal: continue the approved fixed-ion single-molecule E-driven feedback design.
Architecture: use the existing isolated endpoint-E worker without changing its
integration order. Trial the two fields, subtract own vacuum E in Fortran, advance
one electronic interval, map dipole difference to charge current in Fortran,
rollback fields and apply accepted current. Validate the actual selected endpoint
before committing. Cancel the electronic interval using supported replay restore.
This is not direct electronic restore and is not a claim of midpoint quantum
energy conservation. Keep one explicit field hold and report electronic time
as field time minus dt. Require initially zero selected E, matching worker E0.

1. Add a failing protocol test for propose/validate/commit/cancel, map sign,
   independence, and rejection followed by retry. Implement a private adapter in
   the existing dipole_source module; retain state ownership in isolated_tddft.
2. Extend the existing molecular MPI fixture with explicit isolated binary,
   validation receipt and input bundle arguments. No implicit download or SCF.
3. Exercise native FDTD in vacuum and with LD using a labelled dummy worker;
   compare accepted E drive to independently replayed dipoles. Inject an endpoint
   rejection; rollback both fields/electrons and compare retry to continuous.
4. Run isolated/dipole tests and full suite with the audited lightweight kernel.
   Update paired guides, AI catalog and checklist, distinguishing dummy tests
   from real TDDFT. No new GS data may be collected for this phase.
5. Real isolated TDDFT/FDTD remains pending a validated E worker and existing GS.
   Quantum work/norm and total electromagnetic energy need separate evidence.

## Follow-up: measured electron number

Expose the actual occupied-density integral from Fortran with each isolated
endpoint, including initial state. Use the existing calc_dip integration and
Hvol; do not infer electron number from input or dipole. Move the unpublished
isolated wire/evidence/checkpoint contract to revision2 with seven endpoint
values (six dipole components plus measured electron number); reject old frames.
Include this scalar in replay snapshots and compare it after rollback/retry.
Numerical work stays Fortran. Reuse the authorized C2H2 GS; no new SCF campaign.
Validate finite/NaN rejection, state independence and disk/memory replay, then
actual C2H2 in the existing field fixture. Norm loss1e-5 is diagnostic, not a stop.

## Follow-up: electronic internal-energy / work observations

Preserve the existing eigen-energy cadence and propagation order. E3 carries
six dipole values, electron number, field-free electronic energy and a 0/1 fresh
flag. Only steps0,1,capacity and multiples of out_rt_energy_step are fresh; expose
None otherwise. Do not force an energy recomputation or add a reduction.
Update unpublished evidence/checkpoint revisions without old-format fallbacks.
Fortran accumulates endpoint-held E·delta-p (the actually consumed molecular
field) and midpoint E·delta-p (field-work convention) separately. Compare both
against internal energy change only at fresh endpoints. Report residuals without
claiming conservation or imposing a new energy threshold. Use analytic sign/time
controls, finite rejection, snapshot/replay and the authorized existing C2H2 GS.
Keep quantum work distinct from full Maxwell energy and vacuum reaction physics.

## Follow-up: align electronic drive with field midpoint work

Use Fortran to average the selected left/right fields after own-source vacuum
subtraction. Feed that interval-constant midpoint E to the unchanged isolated
worker. Keep endpoint-held driving as an explicit numerical comparison option,
not an automatic fallback. Default the development coupled fixture to midpoint;
standalone IsolatedTDDFTState.step(E) and all solver update orders stay unchanged.
Commit the accepted right field only after validating the actual midpoint drive;
cancellation preserves both quantum state and previous selected field. Report
work of the actual applied drive separately from field endpoint/midpoint work,
with raw zero-field drift retained. Verify independent prescribed-drive replay,
injected reject/restore/retry, and quantum/field work alignment with the existing
C2H2 GS. Do not claim full-volume Maxwell energy/flux conservation yet.
