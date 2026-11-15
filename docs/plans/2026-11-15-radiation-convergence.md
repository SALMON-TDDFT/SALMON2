# One-way radiation refinement plan

Continue R4 verification before choosing self-field feedback. No solver/API
changes. Use a compact sin^4 charge dipole at a fixed source and observe its
transverse electric field at fixed separation. Compare to full retarded point
-dipole field (near, induction and radiation terms), not a far-field-only formula.
Refine dt independently at fixed dx, then dx at fixed dt; use bounded two-rank
native runs, one field rank, within four cores. Keep pulse, physical positions
and amplitude fixed. Record geometry, finite-window errors and source work.

The continuum dipole radiation energy is an analytic benchmark; finite-time
source work is a separate measured quantity, not closed-surface flux. Do not
claim exact energy balance, converged physical molecular spectrum or self-field
renormalization from these tests. Add reproducible tool and paired documentation,
record failures and observed rates honestly. No Actions or push.
