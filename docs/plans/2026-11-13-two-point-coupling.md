# Two-point native multiscale connection validation

Use the unchanged SALMON 1D Weyl-gauge FDTD as the electromagnetic component and
connect two independent SALMONpy A-driven TDDFT workers. Compare the full closed
loop with ordinary two-rank native multi_scale_maxwell_tddft for PZ and TB-mBJ.
The validation-only process bridge calls the native weyl_init/weyl_calc and
calc_Ac_ext_t routines from audited objects; it does not reimplement Maxwell.
Use identity crystal rotation, two adjacent material macro cells, vacuum on both
sides, native absorbing boundaries, and the same small Si GS in both cells.
Keep the current sign and time staggering from main_ms: FDTD uses previous-step
charge current, then TDDFT takes endpoint A and returns electronic matter current;
store its negative as charge current for the next step. Native E/H output is
computed from A with the existing centered time/spatial differences.

Compare every grid point's A/E/H and old charge-current output and both material
points' new electronic currents. Include a no-feedback control and replay each
recorded A history into an independent TDDFT state. This is validation only;
there is no new public E/H FDTD step API, rotation/general geometry or production
multiscale orchestration in this change.
