# Ni LSDA spatial-grid convergence before k convergence

User authorizes starting this while the old HSE10fs comparison continues.
Vary Nr=8^3,10^3,12^3 at fixed full k mesh2^3 (mirror6 representatives).
Same primitive cell, SOC+NLCC KY potential,300K,24spinor states,16electrons,
PauliM[1,-1,0],E001,1.55eV,1e11W/cm2,CEP0, Acos2 tw1=27.4fs,
RT30fs,dt=.0005fs,Taylor8 middlepoint. Compare pulse and no-pulse trajectories,
magnetization and absorbed energy per atom; establish Nr before varying k.
The selected Nr/k will be shared among functionals; user requests LSDA-only
convergence assessment. No assumption of convergence after only3grids.

Reuse unchanged audited v5 executable SHA256
f50f1bbcd918430fa57b60d84acbfce59cb0f7acbffa03fd9ea84263a99b86b7.
No new SALMON build. Nr8 uses existing converged mirror LSDA GS.
Nr10/12 use periodic complex Fourier interpolation of Nr8 solely as initial
guess; no physical result accepted without freshSCF residual<1e-6/endSALMON.
Fourier embedding preserves complex modes including the original negative
Nyquist mode; all spinor Gram matrix changes checked <1e-11 (actual<8e-16).
Occupation/k ordering is unchanged. All GS files remain /private/tmp.

Driver /private/tmp/ni-lsda-nr-30fs-20261024/run.py, session14999.
Status status.json, stdout driver.log. First SCF10,12 then RT8,10,12;
pulse/no-pulse run sequentially to bound added CPU demand alongside HSE.
Do not duplicate or overwrite runs. Reference script code copy in outer
workspace work/ni-lsda-nr-30fs-20261024. Canonical final RT/analysis store:
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/
ni-lsda-nr-30fs-20261024. Collect RT/analysis only, no GS data/logs.

Initial progress: Nr10 SCF81 iterations residual9.38465e-7/endSALMON;
Nr12 SCF started. Final RT comparisons pending.
