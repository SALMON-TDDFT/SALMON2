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
