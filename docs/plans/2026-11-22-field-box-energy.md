# Closed-box field/material energy observation

Goal: verify the authorized C2H2/LD/Maxwell case without changing solver updates.
Use a source-free box strictly inside PML; six signed face fluxes account for
all scattering directions. No inference of PML heat or missing source work.

Pure shared Fortran kernel reads old/new raw Yee fields, old/new LD currents and
external current. Compute leapfrog modified energy (E_n²+H_minus·H_plus)/8pi,
physical centered energy, LD work and source work, and independent six-face
Poynting transport with the matching discrete time/space staggering. Report both
energy residuals. Initially scope the reader to scalar epsilon_inf=mu=1,sigma=0,
no tensor/subpixel/PML/source in the chosen box; unsupported diagnostics must not
prevent the underlying calculation. LD work includes pole storage plus heat.

Use existing accepted-attempt backup as old state; no new physical caches or file
outputs. Observe only accepted steps (including initial hold), piggyback scalar
observations on existing accepted-E packets. One field-group packed reduction,
no new application-wide collective. Rollback/retry discards attempted observations.

TDD: manufactured source-free/source/LD cases, reversed-source sign, partition
additivity, finite/NaN rejection and unmodified outputs. Then actual vacuum/LD
C2H2 plus retry, unchanged prior trajectories, source/wheel suites and manuals.
Joint matter balance explicitly includes own-vacuum work; separate pointwise
quantum error from larger electromagnetic flux/storage scales. Reuse existing GS.
