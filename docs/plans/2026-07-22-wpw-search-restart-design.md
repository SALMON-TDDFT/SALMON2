# WPW Search-history Restart Comparison Design

## Decision

Add a default-on `yn_dg_wpw_search_history` control.  With the default `y`,
the current LOBPCG recurrence and `Z=[Q,R,P]` are unchanged.  With explicit
`n`, clear `P` after every reduced Rayleigh-Ritz update, so the next iteration
uses `Z=[Q,R,0]`, equivalent to a complete two-block restart.

The comparison uses the same normalized occupied-W/P basis, deterministic
extra states, metric cutoff, tolerance, iteration cap, and optional
preconditioner setting.  The first physical comparison uses
`yn_dg_wpw_preconditioner='n'` in both cases.  Existing structural, fixed-H,
checkpoint, and RT gates remain unchanged.

## Interpretation

- Restart improves convergence: the accumulated `P` recurrence is defective.
- Restart behaves similarly: investigate retained Rayleigh-Ritz update or
  operator/basis conditioning.
- Restart worsens convergence: retain history, then test periodic restart.

No production default changes in this task.
